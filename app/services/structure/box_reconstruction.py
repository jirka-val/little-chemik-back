"""
Rekonstrukce periodického boxu z rozložení vody, když nahrané PDB nemá
CRYST1 (export z VMD/PyMOL/Chimery, převod z .gro/.dcd bez boxu ...).

Vstup jsou jen kyslíky vod. Solvatovaný box je vodou vyplněný rovnoměrně až
k hranám, takže tvar buňky jde z oblaku bodů poznat:

1. Šířková funkce w(u) = max(P.u) - min(P.u) na husté polosféře směrů.
   Periodická buňka má protilehlé stěny rovnoběžné, lokální minima w(u)
   jsou proto normály párů stěn. Zpřesní se lokálním hledáním a proložením
   rovin body okrajových vrstev. Falešná minima (rohy, šum) odfiltruje
   "plochost" - kolik bodů leží v tenké vrstvě u obou stěn.
2. Počet párů stěn určí tvar:
     3 = rovnoběžnostěn (krychle, kvádr, kosý kvádr - i komolý oktaedr
         nebo dodekaedr zapsaný jako triklinická buňka),
     6 = rhombický dodekaedr (Wigner-Seitzova buňka FCC),
     7 = komolý oktaedr (Wigner-Seitzova buňka BCC).
   Výběr podle nejmenšího objemu by nefungoval - obecný kosý kvádr obsahuje
   všechny ostatní tvary jako speciální případ, vyhrál by vždycky.
3. Natočení: u TO/RD se ideální těleso přiloží na naměřené normály
   (Kabsch s hledáním přiřazení stěn), u rovnoběžnostěnu se mřížka
   spočítá z normál a šířek a pak se tvar rozpozná z redukované báze.
4. Velikost: voda na hraně není až u hrany buňky, box je o mezeru větší
   než oblak. Rozměr se dorovná testem periodicity - nejbližší vzdálenosti
   vod přes hranu (proti periodickým obrazům) musí odpovídat vodě uvnitř.
   Porovnává se dolní kvartil (citlivější na překryvy než medián).
5. Výsledek: tvar, parametry CRYST1 a transformace souřadnic do standardní
   orientace PDB (a podél x, b v rovině xy): X_std = (X - center) @ R^T + std_center.

Jen numpy (scipy v produkčním image není). Ověřeno na boxech z builderu
(krychle, komolý oktaedr; hrany do 0.2 A) a syntetických boxech všech tvarů
v náhodné orientaci.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

MIN_WATERS = 300
COS_TO = 1.0 / 3.0          # |cos 109.4712|
QUANTILE = 25.0             # kvantil nejbližších vzdáleností: hrana vs. uvnitř
CLASH_DISTANCE = 2.2        # A, O-O blíž = překryv
GAP_DISTANCE = 3.6          # A, nejbližší O-O dál = díra
# o kolik smí být na hraně víc překryvů/děr než uvnitř vody
QUALITY_MARGIN = 0.05
# kolik nejslabších kandidátů na stěnu se smí zahodit (falešná minima na hranách)
MAX_DROPPED_FACES = 2
# průměrná odchylka normál od ideálního TO/RD, 1 - cos(8 deg)
MAX_ORIENTATION_ERROR = 0.01

_COMBOS = np.array([m for m in itertools.product((-1, 0, 1), repeat=3) if any(m)], dtype=float)

IDEAL_PARAMS = {
    "cubic": {"alpha": 90.0, "beta": 90.0, "gamma": 90.0},
    "rectangular": {"alpha": 90.0, "beta": 90.0, "gamma": 90.0},
    # AMBER konvence: a=b=c, všechny úhly 109.4712
    "truncated octahedron": {"alpha": 109.4712, "beta": 109.4712, "gamma": 109.4712},
    # GROMACS konvence (xy čtverec)
    "rhombic dodecahedron": {"alpha": 60.0, "beta": 60.0, "gamma": 90.0},
}


@dataclass
class BoxReconstruction:
    ok: bool
    reason: str = ""
    shape: str = ""                  # cubic / rectangular / truncated octahedron / rhombic dodecahedron / triclinic
    n_face_pairs: int = 0
    params: Dict[str, float] = field(default_factory=dict)
    volume: float = 0.0
    rotation: Optional[np.ndarray] = None
    center: Optional[np.ndarray] = None
    std_center: Optional[np.ndarray] = None
    clash_fraction: float = 0.0
    gap_fraction: float = 0.0
    interior_clash_fraction: float = 0.0
    interior_gap_fraction: float = 0.0

    def transform(self, coords: np.ndarray) -> np.ndarray:
        """Souřadnice do standardní orientace rekonstruovaného boxu."""
        return (np.asarray(coords, dtype=float) - self.center) @ self.rotation.T + self.std_center

    def summary(self) -> Dict[str, object]:
        return {
            "shape": self.shape,
            "face_pairs": self.n_face_pairs,
            "rotation": np.round(self.rotation, 6).tolist() if self.rotation is not None else None,
            "translation": (np.round(self.std_center - self.center @ self.rotation.T, 3).tolist()
                            if self.rotation is not None else None),
            "boundary_clash_fraction": round(self.clash_fraction, 3),
            "boundary_gap_fraction": round(self.gap_fraction, 3),
            "interior_clash_fraction": round(self.interior_clash_fraction, 3),
            "interior_gap_fraction": round(self.interior_gap_fraction, 3),
        }


# ---------------------------------------------------------------------------
# Mřížky
# ---------------------------------------------------------------------------

def standard_vectors(p: Dict[str, float]) -> np.ndarray:
    """Mřížkové vektory z parametrů CRYST1 (a podél x, b v rovině xy)."""
    a, b, c = p["a"], p["b"], p["c"]
    al, be, ga = (math.radians(p[k]) for k in ("alpha", "beta", "gamma"))
    bx, by = b * math.cos(ga), b * math.sin(ga)
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return np.array([[a, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz]])


def cell_params(B: np.ndarray) -> Dict[str, float]:
    def angle(x, y):
        return math.degrees(math.acos(float(np.clip(np.dot(x, y) / np.linalg.norm(x) / np.linalg.norm(y), -1, 1))))
    a, b, c = (float(np.linalg.norm(v)) for v in B)
    return {"a": a, "b": b, "c": c, "alpha": angle(B[1], B[2]), "beta": angle(B[0], B[2]), "gamma": angle(B[0], B[1])}


def reduce_basis(B: np.ndarray) -> np.ndarray:
    """3 nejkratší lineárně nezávislé mřížkové vektory (řádky)."""
    combos = np.array([m for m in itertools.product(range(-2, 3), repeat=3) if any(m)])
    vecs = combos @ B
    basis: List[np.ndarray] = []
    for v in vecs[np.argsort(np.linalg.norm(vecs, axis=1), kind="stable")]:
        if not basis:
            basis.append(v)
        elif len(basis) == 1:
            if np.linalg.norm(np.cross(basis[0], v)) > 1e-6 * np.dot(v, v):
                basis.append(v)
        elif abs(np.linalg.det(np.array(basis + [v]))) > 1e-6 * np.linalg.norm(v) ** 3:
            basis.append(v)
            break
    return np.array(basis)


def classify_lattice(B: np.ndarray, len_tol: float = 0.03, cos_tol: float = 0.04) -> Tuple[str, np.ndarray]:
    """Tvar z redukované báze + báze v konvenci pro CRYST1."""
    R = reduce_basis(B)
    L = np.linalg.norm(R, axis=1)
    equal = np.ptp(L) < len_tol * L.mean()
    cs = np.array([np.dot(R[1], R[2]) / L[1] / L[2], np.dot(R[0], R[2]) / L[0] / L[2],
                   np.dot(R[0], R[1]) / L[0] / L[1]])
    if np.all(np.abs(cs) < cos_tol):
        return ("cubic" if equal else "rectangular"), R
    if equal and np.all(np.abs(np.abs(cs) - COS_TO) < cos_tol):
        for signs in itertools.product((1, -1), repeat=3):
            S = R * np.array(signs)[:, None]
            if all(np.dot(S[i], S[j]) < 0 for i, j in ((0, 1), (0, 2), (1, 2))):
                return "truncated octahedron", S
        return "truncated octahedron", R
    if equal and np.all((np.abs(np.abs(cs) - 0.5) < cos_tol) | (np.abs(cs) < cos_tol)):
        for perm in itertools.permutations(range(3)):
            for signs in itertools.product((1, -1), repeat=3):
                S = R[list(perm)] * np.array(signs)[:, None]
                p = cell_params(S)
                if abs(p["alpha"] - 60) < 2 and abs(p["beta"] - 60) < 2 and abs(p["gamma"] - 90) < 2:
                    return "rhombic dodecahedron", S
    return "triclinic", R


def _kabsch(Q: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Rotace Rm s Q @ Rm ~ M (řádkové vektory)."""
    u, _, vt = np.linalg.svd(Q.T @ M)
    d = np.sign(np.linalg.det(u @ vt))
    return u @ np.diag([1.0, 1.0, d]) @ vt


def _ideal_params(shape: str, conv: np.ndarray) -> Dict[str, float]:
    L = np.linalg.norm(conv, axis=1)
    if shape == "rectangular":
        return {"a": L[0], "b": L[1], "c": L[2], **IDEAL_PARAMS[shape]}
    if shape in IDEAL_PARAMS:
        return {"a": L.mean(), "b": L.mean(), "c": L.mean(), **IDEAL_PARAMS[shape]}
    return cell_params(conv)


# ---------------------------------------------------------------------------
# Stěny
# ---------------------------------------------------------------------------

def _hemisphere(m: int) -> np.ndarray:
    i = np.arange(m) + 0.5
    z = i / m
    phi = i * math.pi * (3.0 - math.sqrt(5.0))
    r = np.sqrt(1.0 - z * z)
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)


def _width(P: np.ndarray, u: np.ndarray) -> float:
    p = P @ u
    return float(p.max() - p.min())


def _widths(P: np.ndarray, U: np.ndarray, chunk: int = 256) -> np.ndarray:
    out = np.empty(len(U))
    for s in range(0, len(U), chunk):
        proj = P @ U[s:s + chunk].T
        out[s:s + chunk] = proj.max(axis=0) - proj.min(axis=0)
    return out


def _refine_min(P: np.ndarray, u: np.ndarray, step_deg: float = 2.0, tol_deg: float = 0.005) -> np.ndarray:
    best, wb = u / np.linalg.norm(u), _width(P, u)
    step = math.radians(step_deg)
    dirs = [(math.cos(t), math.sin(t)) for t in np.linspace(0, 2 * math.pi, 8, endpoint=False)]
    while step > math.radians(tol_deg):
        a = np.array([1.0, 0, 0]) if abs(best[0]) < 0.9 else np.array([0, 1.0, 0])
        t1 = np.cross(best, a)
        t1 /= np.linalg.norm(t1)
        t2 = np.cross(best, t1)
        for c, s in dirs:
            cand = best + step * (c * t1 + s * t2)
            cand /= np.linalg.norm(cand)
            w = _width(P, cand)
            if w < wb - 1e-9:
                best, wb = cand, w
                break
        else:
            step *= 0.5
    return best


def _refine_plane(P: np.ndarray, u: np.ndarray, layer: float = 2.5, iters: int = 3) -> np.ndarray:
    """Normála ze dvou rovnoběžných rovin proložených body okrajových vrstev."""
    for _ in range(iters):
        p = P @ u
        top, bot = P[p > p.max() - layer], P[p < p.min() + layer]
        if len(top) < 10 or len(bot) < 10:
            return u
        n = np.linalg.svd(np.vstack([top - top.mean(0), bot - bot.mean(0)]), full_matrices=False)[2][-1]
        u = n if np.dot(n, u) > 0 else -n
    return u


def _face_score(P: np.ndarray, u: np.ndarray, layer: float = 1.5) -> float:
    """Jak moc je směr opravdu stěna: počet bodů v okrajových vrstvách krát
    "plošnost" vrstvy (poměr 2. a 1. singulární hodnoty). Hrana mezi stěnami
    dá taky lokální minimum šířky, ale málo bodů téměř na přímce."""
    p = P @ u
    score = 0.0
    for X in (P[p > p.max() - layer], P[p < p.min() + layer]):
        if len(X) < 3:
            return 0.0
        s = np.linalg.svd(X - X.mean(axis=0), compute_uv=False)
        score += len(X) * (s[1] / s[0] if s[0] > 0 else 0.0)
    return score


def find_face_pairs(P: np.ndarray, m: int = 3000, k_nb: int = 12, merge_deg: float = 4.0,
                    layer: float = 1.5, flat_rel: float = 0.12) -> Tuple[List[np.ndarray], List[float]]:
    """Normály párů rovnoběžných stěn oblaku bodů a jeho šířky v jejich směru,
    seřazené od nejvýraznější stěny (_face_score)."""
    U = _hemisphere(m)
    W = _widths(P, U)
    G = np.abs(U @ U.T)                     # antipodální směry jsou stejné
    np.fill_diagonal(G, -1)
    nb = np.argpartition(-G, k_nb, axis=1)[:, :k_nb]
    refined: List[np.ndarray] = []
    for i in np.flatnonzero(W <= W[nb].min(axis=1)):
        u = _refine_min(P, U[i])
        if all(abs(np.dot(u, r)) < math.cos(math.radians(merge_deg)) for r in refined):
            refined.append(u)
    if not refined:
        return [], []
    flat = []
    for u in refined:
        p = P @ u
        flat.append(int(np.sum(p > p.max() - layer) + np.sum(p < p.min() + layer)))
    keep = [u for u, f in zip(refined, flat) if f >= flat_rel * max(flat)]
    normals = [_refine_plane(P, u) for u in keep]
    normals.sort(key=lambda u: -_face_score(P, u))
    return normals, [_width(P, u) for u in normals]


# ---------------------------------------------------------------------------
# Test periodicity
# ---------------------------------------------------------------------------

def _stat(x: np.ndarray) -> float:
    return float(np.percentile(x, QUANTILE))


def _face_nn(P: np.ndarray, n: np.ndarray, lattice: np.ndarray, layer: float = 1.5, reach: float = 6.0) -> np.ndarray:
    """Nejbližší vzdálenosti vod v okrajové vrstvě stěny (+n): skutečné vody
    + periodické obrazy přes všechny sousední buňky."""
    p = P @ n
    hi = p.max()
    L = P[p > hi - layer]
    imgs = (P[None, :, :] + (_COMBOS @ lattice)[:, None, :]).reshape(-1, 3)
    pi = imgs @ n
    nbrs = np.vstack([P[p > hi - reach], imgs[(pi > hi - reach) & (pi < hi + reach)]])
    out = np.empty(len(L))
    for s in range(0, len(L), 256):
        d2 = ((L[s:s + 256, None, :] - nbrs[None, :, :]) ** 2).sum(-1)
        d2[d2 < 1e-12] = np.inf
        out[s:s + 256] = np.sqrt(d2.min(axis=1))
    return out


def _boundary_nn(P: np.ndarray, normals: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    return np.concatenate([_face_nn(P, s * n, lattice) for n in normals for s in (1, -1)])


def _residuals(P: np.ndarray, normals: np.ndarray, lattice: np.ndarray, target: float) -> np.ndarray:
    return np.array([_stat(np.concatenate([_face_nn(P, n, lattice), _face_nn(P, -n, lattice)]))
                     for n in normals]) - target


def _interior_nn(P: np.ndarray, rng: np.random.Generator, sample: int = 600) -> np.ndarray:
    idx = rng.choice(len(P), min(sample, len(P)), replace=False)
    out = np.empty(len(idx))
    for s in range(0, len(idx), 100):
        part = idx[s:s + 100]
        d2 = ((P[part, None, :] - P[None, :, :]) ** 2).sum(-1)
        d2[np.arange(len(part)), part] = np.inf
        out[s:s + 100] = np.sqrt(d2.min(axis=1))
    return out


# ---------------------------------------------------------------------------
# Rovnoběžnostěn (3 páry stěn): obecná mřížka
# ---------------------------------------------------------------------------

def _parallelepiped(normals: List[np.ndarray], W) -> np.ndarray:
    """Mřížka s a_i . n_j = W_i delta_ij (řádky)."""
    return np.linalg.inv(np.array(normals)).T * np.asarray(W)[:, None]


def _fit_parallelepiped(P, normals, target, w0, iters: int = 15, tol: float = 0.02) -> np.ndarray:
    W = np.asarray(w0, dtype=float) + 0.5
    Nm = np.array(normals)
    for _ in range(iters):
        r = _residuals(P, Nm, _parallelepiped(normals, W), target)
        W = W - r
        if np.max(np.abs(r)) < tol:
            break
    return _parallelepiped(normals, W)


# ---------------------------------------------------------------------------
# Komolý oktaedr / rhombický dodekaedr (7 / 6 párů stěn)
# ---------------------------------------------------------------------------

def _ws_face_vectors(params: Dict[str, float], n_pairs: int) -> np.ndarray:
    """Mřížkové vektory kolmé na stěny Wigner-Seitzovy buňky (jeden z páru)."""
    combos = np.array([m for m in itertools.product(range(-2, 3), repeat=3) if any(m)])
    vecs = combos @ standard_vectors(params)
    vecs = vecs[np.argsort(np.linalg.norm(vecs, axis=1), kind="stable")][:2 * n_pairs]
    out: List[np.ndarray] = []
    for v in vecs:
        if all(np.linalg.norm(v + o) > 1e-6 for o in out):
            out.append(v)
    return np.array(out)


def _orient_ws(normals: List[np.ndarray], widths: List[float], shape: str):
    """Natočení a velikost ideálního TO/RD podle naměřených normál stěn."""
    unit = {"a": 1.0, "b": 1.0, "c": 1.0, **IDEAL_PARAMS[shape]}
    V = _ws_face_vectors(unit, len(normals))
    Q = V / np.linalg.norm(V, axis=1)[:, None]
    Nm = np.array(normals)
    m1, m2 = Nm[0], Nm[1]
    c12 = abs(np.dot(m1, m2))
    best_err, best = np.inf, None
    for i, j in itertools.permutations(range(len(Q)), 2):
        if abs(abs(np.dot(Q[i], Q[j])) - c12) > 0.08:
            continue
        for s1, s2 in itertools.product((1, -1), repeat=2):
            q1, q2 = s1 * Q[i], s2 * Q[j]
            Rm = _kabsch(np.array([q1, q2, np.cross(q1, q2)]), np.array([m1, m2, np.cross(m1, m2)]))
            QR = Q @ Rm
            match = np.argmax(np.abs(Nm @ QR.T), axis=1)
            if len(set(match.tolist())) != len(Nm):
                continue
            err = float(np.sum(1 - np.abs(np.einsum("ij,ij->i", Nm, QR[match]))))
            if err < best_err:
                best_err, best = err, (Rm, match)
    if best is None or best_err / len(Nm) > MAX_ORIENTATION_ERROR:
        return None
    Rm, match = best
    signs = np.sign(np.einsum("ij,ij->i", Nm, (Q @ Rm)[match]))
    Rm = _kabsch(Q[match] * signs[:, None], Nm)           # přes všechny stěny
    size = float(np.mean(np.asarray(widths) / np.linalg.norm(V[match], axis=1)))
    return {**unit, "a": size, "b": size, "c": size}, Rm


# ---------------------------------------------------------------------------
# Dorovnání velikosti ideálního tělesa
# ---------------------------------------------------------------------------

def _fit_size(P, normals, widths, shape, params, Rm, target, iters: int = 12, tol: float = 0.01):
    """Pevné natočení, velikost (L, u kvádru a/b/c) z testu periodicity."""
    p = dict(params)
    Nm = np.array(normals)
    for _ in range(iters):
        lat = standard_vectors(p) @ Rm
        r = _residuals(P, Nm, lat, target)
        if shape == "rectangular":
            axes = np.argmax(np.abs(Nm @ (lat / np.linalg.norm(lat, axis=1)[:, None]).T), axis=1)
            for ax, key in enumerate(("a", "b", "c")):
                if (axes == ax).any():
                    p[key] -= float(np.mean(r[axes == ax]))
            if np.max(np.abs(r)) < tol:
                break
        else:
            scale = 1.0 - float(np.mean(r / np.asarray(widths)))
            for key in ("a", "b", "c"):
                p[key] *= scale
            if abs(float(np.mean(r))) < tol:
                break
    return p


# ---------------------------------------------------------------------------
# Veřejné API
# ---------------------------------------------------------------------------

def reconstruct_box(water_oxygens: np.ndarray, seed: int = 0) -> BoxReconstruction:
    """Box z kyslíků vod. ok=False (s důvodem), když tvar nejde spolehlivě poznat."""
    P = np.asarray(water_oxygens, dtype=float)
    if len(P) < MIN_WATERS:
        return BoxReconstruction(False, f"too few water molecules ({len(P)}) to recognize a box")
    normals, wcloud = find_face_pairs(P)
    k = len(normals)
    target = _stat(_interior_nn(P, np.random.default_rng(seed)))
    # Kandidáti jsou seřazení od nejvýraznější stěny; nejvýš MAX_DROPPED_FACES
    # nejslabších (typicky hrana mezi stěnami) se smí zahodit.
    first_failure = None
    for n_pairs in (7, 6, 3):
        if 0 <= k - n_pairs <= MAX_DROPPED_FACES:
            fit = _fit_faces(P, normals[:n_pairs], wcloud[:n_pairs], target, seed)
            if fit.ok:
                return fit
            first_failure = first_failure or fit
    return first_failure or BoxReconstruction(
        False, f"the water does not form a periodic cell ({k} pairs of parallel faces)", n_face_pairs=k)


def _fit_faces(P: np.ndarray, normals: List[np.ndarray], wcloud: List[float], target: float,
               seed: int) -> BoxReconstruction:
    k = len(normals)
    Nm = np.array(normals)
    proj = P @ Nm.T
    center = np.linalg.lstsq(Nm, (proj.min(axis=0) + proj.max(axis=0)) / 2, rcond=None)[0]

    if k in (6, 7):
        shape = "rhombic dodecahedron" if k == 6 else "truncated octahedron"
        oriented = _orient_ws(normals, wcloud, shape)
        if oriented is None:
            return BoxReconstruction(False, f"the faces do not match a {shape}", n_face_pairs=k)
        params, Rm = oriented
    else:
        lattice = _fit_parallelepiped(P, normals, target, wcloud)
        shape, conv = classify_lattice(lattice)
        if np.linalg.det(conv) < 0:
            conv = -conv
        params = _ideal_params(shape, conv)
        Rm = _kabsch(standard_vectors(params), conv)
        if shape == "triclinic":
            params = cell_params(conv)
    if shape != "triclinic":
        params = _fit_size(P, normals, wcloud, shape, params, Rm, target)

    S = standard_vectors(params)
    lattice = S @ Rm
    boundary = _boundary_nn(P, Nm, lattice)
    interior = _interior_nn(P, np.random.default_rng(seed + 1))
    fit = BoxReconstruction(
        ok=True, shape=shape, n_face_pairs=k, params={key: float(v) for key, v in params.items()},
        volume=float(abs(np.linalg.det(S))), rotation=Rm, center=center, std_center=S.sum(axis=0) / 2,
        clash_fraction=float(np.mean(boundary < CLASH_DISTANCE)),
        gap_fraction=float(np.mean(boundary > GAP_DISTANCE)),
        interior_clash_fraction=float(np.mean(interior < CLASH_DISTANCE)),
        interior_gap_fraction=float(np.mean(interior > GAP_DISTANCE)),
    )

    def worse(boundary_fraction: float, interior_fraction: float) -> bool:
        return boundary_fraction > interior_fraction + max(QUALITY_MARGIN, 0.25 * interior_fraction)

    if worse(fit.clash_fraction, fit.interior_clash_fraction) or worse(fit.gap_fraction, fit.interior_gap_fraction):
        fit.ok = False
        fit.reason = (f"a {shape} was found, but the water does not continue across its faces "
                      f"({fit.clash_fraction * 100:.0f} % overlaps, {fit.gap_fraction * 100:.0f} % gaps at the boundary)")
    return fit
