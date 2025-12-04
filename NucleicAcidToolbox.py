from math import fabs, sqrt, asin, pi, isclose
import numpy as np

########################################################################################################
#
# general functions of Nucleic Acid Toolbox
#
########################################################################################################
ref_coords = {
    "A": {
        "C1'": [0.000, 0.000, 0.000],
        "C5": [0.614, 3.572, 0.000],
        "C6": [1.559, 4.612, 0.000],
        "N1": [2.866, 4.271, 0.000],
        "C2": [3.189, 2.973, 0.000],
        "N3": [2.395, 1.905, 0.000],
        "C4": [1.104, 2.277, 0.000]
    },
    "C": {
        "C1'": [0.000, 0.000, 0.000],
        "N1": [-0.000, 1.470, 0.000],
        "C2": [1.232, 2.128, 0.000],
        "N3": [1.259, 3.481, 0.000],
        "C4": [0.116, 4.170, 0.000],
        "C5": [-1.153, 3.525, 0.000],
        "C6": [-1.165, 2.186, 0.000]
    },
    "G": {
        "C1'": [0.000, 0.000, 0.000],
        "C5": [0.608, 3.565, 0.000],
        "C6": [1.521, 4.650, 0.000],
        "N1": [2.840, 4.211, 0.000],
        "C2": [3.232, 2.894, 0.000],
        "N3": [2.392, 1.873, 0.000],
        "C4": [1.104, 2.277, 0.000],
    },
    "U": {
        "C1'": [0.000, 0.000, 0.000],
        "N1": [-0.000, 1.470, 0.000],
        "C2": [1.218, 2.121, 0.000],
        "N3": [1.142, 3.491, 0.000],
        "C4": [-0.005, 4.259, 0.000],
        "C5": [-1.224, 3.512, 0.000],
        "C6": [-1.182, 2.175, 0.000]
    }
}

classification_combinations = {
    "AA": ["cWW", "tWW", "tWH"],
    "AC": ["cWW", "tWW"],
    "AG": ["cWW", "cWH", "tWH"],
    "AU": ["cWW", "tWW", "cWH"],
    "CA": ["cWW", "tWW", "tWH"],
    "CC": ["cWW", "tWW", "cWH", "tWH"],
    "CG": ["cWW", "tWW", "cWH", "tWH"],
    "CU": ["cWW", "tWW", "cWH"],
    "GA": ["cWW", "cWH"],
    "GC": ["cWW", "tWW"],
    "GG": ["tWW", "cWH", "tWH"],
    "GU": ["cWW", "tWW", "tWH"],
    "UA": ["cWW", "tWW", "cWH", "tWH"],
    "UC": ["cWW", "tWW"],
    "UG": ["cWW", "tWW", "cWH", "tWH"],
    "UU": ["cWW", "tWW", "cWH", "tWH"]
}

# average values taken from excel sheet
avg_vector = {
    "cWW": {
        "AU": np.array([8.57362, 5.92499, 0.0, 0.89515, -1.29077, 0.0]),
        "UA": np.array([8.55411, 5.95312, 0.0, 0.89515, -1.29077, 0.0]),
        "GC": np.array([8.53586, 6.32895, 0.0, 0.91060, -1.27992, 0.0]),
        "CG": np.array([8.77779, 5.98990, 0.0, 0.91060, -1.27992, 0.0]),
        "GU": np.array([6.98013, 7.54000, 0.0, 0.90408, -1.28454, 0.0]),
        "UG": np.array([9.45326, 4.02612, 0.0, 0.90408, -1.28454, 0.0]),
        "AG": np.array([9.13434, 5.59292, 0.0, 1.07873, -1.14181, 0.0]),
        "GA": np.array([9.09762, 8.63178, 0.0, 1.07873, -1.14181, 0.0]),
        "AC": np.array([6.75337, 7.83193, 0.0, 0.92261, -1.27130, 0.0]),
        "CA": np.array([9.53980, 3.99242, 0.0, 0.92261, -1.27130, 0.0]),
        "CU": np.array([7.72288, 4.10471, 0.0, 0.85747, -1.31611, 0.0]),
        "UC": np.array([6.87503, 5.40606, 0.0, 0.85747, -1.31611, 0.0]),
        "Aa": np.array([10.3229, 5.95118, 0.0, 0.97632, -1.23053, 0.0]),
        "aA": np.array([8.14233, 8.69945, 0.0, 0.97632, -1.23053, 0.0]),
        "Cc": np.array([9.34386, 2.92242, 0.0, 0.96078, -1.24270, 0.0]),
        "cC": np.array([5.18078, 8.30709, 0.0, 0.96078, -1.24270, 0.0]),
        "Uu": np.array([5.98207, 6.37077, 0.0, 0.83475, -1.33064, 0.0]),
        "uU": np.array([8.33919, 2.61339, 0.0, 0.83475, -1.33064, 0.0])
    },
    "tWW": {
        "AU": np.array([]),
        "UA": np.array([]),
        "GC": np.array([]),
        "CG": np.array([]),
        "GU": np.array([]),
        "UG": np.array([]),
        "AG": np.array([]),
        "GA": np.array([]),
        "AC": np.array([]),
        "CA": np.array([]),
        "CU": np.array([]),
        "UC": np.array([]),
        "Aa": np.array([]),
        "aA": np.array([]),
        "Cc": np.array([]),
        "cC": np.array([]),
        "Gg": np.array([]),
        "gG": np.array([]),
        "Uu": np.array([]),
        "uU": np.array([])
    }
}

# invcov matrixes taken from excel sheet
invcov_matrix = {
    "cWW": {
        "AU": np.array([
            [13.55826, 8.02138, 0.0404, -29.53008, -22.9821, -1.23507],
            [8.02138, 12.96972, -0.07211, -64.40897, -43.15129, -4.01616],
            [0.0404, -0.07211, 5.13112, -20.74689, 30.80657, -20.17125],
            [-29.53008, -64.40897, -20.74689, 562.59037, 94.67343, 110.46312],
            [-22.9821, -43.15129, 30.80657, 94.67343, 472.04396, -106.17188],
            [-1.23507, -4.01616, -20.17125, 110.46312, -106.17188, 147.62622]
        ]),
        "UA": np.array([
            [11.27767, 4.70054, 1.13706, -14.1379, -4.74048, -2.7195],
            [4.70054, 8.8458, 0.42297, -47.80563, -27.71306, -4.02268],
            [1.13706, 0.42297, 4.56653, -33.06255, 44.59317, -18.36334],
            [-14.1379, -47.80563, -33.06255, 648.28216, -165.02073, 158.74867],
            [-4.74048, -27.71306, 44.59317, -165.02073, 691.92657, -172.03707],
            [-2.7195, -4.02268, -18.36334, 158.74867, -172.03707, 138.34679]
        ]),
        "GC": np.array([
            [16.76655, 11.19735, -0.67311, -31.86366, -27.65246, 0.96855],
            [11.19735, 16.84836, -0.33814, -69.93124, -43.68352, -2.37928],
            [-0.67311, -0.33814, 5.70653, -18.26342, 36.72818, -17.9357],
            [-31.86366, -69.93124, -18.26342, 594.22984, 122.57701, 101.80008],
            [-27.65246, -43.68352, 36.72818, 122.57701, 540.17747, -131.96106],
            [0.96855, -2.37928, -17.9357, 101.80008, -131.96106, 138.33946]
        ]),
        "CG": np.array([
            [12.75304, 4.20645, 1.3332, 3.35052, 3.13658, 12.9201],
            [4.20645, 4.94252, 0.32632, -11.77944, -8.13286, 2.91874],
            [1.3332, 0.32632, 0.94621, 2.3077, 2.06322, -2.35375],
            [3.35052, -11.77944, 2.3077, 210.22211, 149.09845, 10.61418],
            [3.13658, -8.13286, 2.06322, 149.09845, 107.29988, 9.94896],
            [12.9201, 2.91874, -2.35375, 10.61418, 9.94896, 137.63921]
        ]),
        "GU": np.array([
            [6.14983, 4.53381, -1.66521, -9.31333, -13.147, 1.15824],
            [4.53381, 9.47389, -1.79479, -34.74102, -27.13904, -1.34473],
            [-1.66521, -1.79479, 5.09925, -17.41295, 30.14442, -14.02503],
            [-9.31333, -34.74102, -17.41295, 338.73393, -28.67187, 75.81462],
            [-13.147, -27.13904, 30.14442, -28.67187, 363.12726, -77.10039],
            [1.15824, -1.34473, -14.02503, 75.81462, -77.10039, 118.08205]
        ]),
        "UG": np.array([
            [8.01659, 2.99509, -0.5014, -8.43183, -25.75277, 3.24744],
            [2.99509, 4.64048, 0.17589, -27.50495, -19.72279, -2.00215],
            [-0.5014, 0.17589, 4.32772, -34.44641, 40.21057, -26.82234],
            [-8.43183, -27.50495, -34.44641, 540.62496, -236.33652, 213.66433],
            [-25.75277, -19.72279, 40.21057, -236.33652, 594.13522, -248.4766],
            [3.24744, -2.00215, -26.82234, 213.66433, -248.4766, 238.97191]
        ]),
        "AG": np.array([
            [2.44971, 1.01804, 0.96861, -7.12068, 10.01609, -4.132],
            [1.01804, 3.62772, 0.46036, -23.19345, -11.03765, -3.1873],
            [0.96861, 0.46036, 1.74525, -14.74863, 12.54138, -6.99841],
            [-7.12068, -23.19345, -14.74863, 284.13513, -15.86877, 64.69728],
            [10.01609, -11.03765, 12.54138, -15.86877, 215.98294, -36.4202],
            [-4.132, -3.1873, -6.99841, 64.69728, -36.4202, 61.64508]
        ]),
        "GA": np.array([
            [2.92186, 1.15512, 0.71438, -2.91434, 2.57092, -1.23552],
            [1.15512, 3.98699, 0.80988, -21.77229, -11.59781, -5.17632],
            [0.71438, 0.80988, 1.38047, -14.58653, 7.20147, -4.37459],
            [-2.91434, -21.77229, -14.58653, 271.12344, -19.37818, 71.03563],
            [2.57092, -11.59781, 7.20147, -19.37818, 149.56117, -10.56494],
            [-1.23552, -5.17632, -4.37459, 71.03563, -10.56494, 44.89276]
        ]),
        "AC": np.array([
            [2.90346, 2.48774, -1.00589, 0.145, 0.05514, 3.98995],
            [2.48774, 2.83322, -1.21548, -1.09971, -0.87473, 4.64892],
            [-1.00589, -1.21548, 1.2515, 0.25224, 0.52475, -2.66625],
            [0.145, -1.09971, 0.25224, 35.76014, 24.55745, 9.61951],
            [0.05514, -0.87473, 0.52475, 24.55745, 18.1417, 9.70055],
            [3.98995, 4.64892, -2.66625, 9.61951, 9.70055, 104.16673]
        ]),
        "CA": np.array([
            [3.93613, 1.27629, -0.60271, 3.36515, -20.18644, 4.67083],
            [1.27629, 0.95065, -0.321, -3.77926, -8.23695, 2.02022],
            [-0.60271, -0.321, 1.78797, -12.4956, 17.10423, -12.42359],
            [3.36515, -3.77926, -12.4956, 221.43384, -117.14727, 97.42466],
            [-20.18644, -8.23695, 17.10423, -117.14727, 265.40114, -107.49323],
            [4.67083, 2.02022, -12.42359, 97.42466, -107.49323, 128.80284]
        ]),
        "CU": np.array([
            [2.27126, 0.08956, -0.17983, -1.46106, -0.31135, 3.10152],
            [0.08956, 0.47483, -0.00457, -1.76338, -0.97067, -0.25292],
            [-0.17983, -0.00457, 0.58186, -0.49464, 0.22372, -2.9957],
            [-1.46106, -1.76338, -0.49464, 24.97706, 12.96167, -0.31092],
            [-0.31135, -0.97067, 0.22372, 12.96167, 8.53667, -0.53641],
            [3.10152, -0.25292, -2.9957, -0.31092, -0.53641, 50.41059]
        ]),
        "UC": np.array([
            [5.04792, 2.51892, -0.66581, -15.64968, -7.96024, -1.80229],
            [2.51892, 2.84934, -0.59506, -15.35832, -13.03778, -0.65301],
            [-0.66581, -0.59506, 1.71249, -5.42764, 10.68248, -7.21742],
            [-15.64968, -15.35832, -5.42764, 164.15558, 7.17199, 45.78259],
            [-7.96024, -13.03778, 10.68248, 7.17199, 158.93999, -48.00957],
            [-1.80229, -0.65301, -7.21742, 45.78259, -48.00957, 59.36237]]),
        "Aa": np.array([
            [2.03576, 0.63584, 0.24659, -4.29611, -0.16063, -3.43854],
            [0.63584, 2.68935, 0.35883, -23.32, -5.04096, -1.80848],
            [0.24659, 0.35883, 1.13002, -8.82944, 8.21718, -6.71275],
            [-4.29611, -23.32, -8.82944, 262.21994, 5.11539, 37.50712],
            [-0.16063, -5.04096, 8.21718, 5.11539, 119.49465, -61.62453],
            [-3.43854, -1.80848, -6.71275, 37.50712, -61.62453, 98.86961]
        ]),
        "aA": np.array([
            [1.63391, 0.42262, 0.34362, -0.83917, 0.24422, 3.32071],
            [0.42262, 1.65174, -0.09366, -6.3054, -7.87257, -1.34914],
            [0.34362, -0.09366, 1.80225, -9.0799, 17.03299, -3.56722],
            [-0.83917, -6.3054, -9.0799, 107.0075, -57.21662, 31.11755],
            [0.24422, -7.87257, 17.03299, -57.21662, 220.03034, -27.90796],
            [3.32071, -1.34914, -3.56722, 31.11755, -27.90796, 45.60939]
        ]),
        "Cc": np.array([
            [10.45331, 4.3869, -2.19694, -41.58137, -36.39587, 2.08609],
            [4.3869, 3.49538, -0.72658, -24.74877, -17.40935, -0.16083],
            [-2.19694, -0.72658, 2.7756, -7.54539, 23.61091, -16.71211],
            [-41.58137, -24.74877, -7.54539, 371.62124, -19.51775, 96.43066],
            [-36.39587, -17.40935, 23.61091, -19.51775, 348.50302, -122.09668],
            [2.08609, -0.16083, -16.71211, 96.43066, -122.09668, 185.88577]
        ]),
        "cC": np.array([
            [8.27004, 4.68337, -2.88004, -4.59606, -15.31663, 0.23793],
            [4.68337, 4.8292, -2.78258, -13.49597, -21.88676, 3.18351],
            [-2.88004, -2.78258, 4.66576, -7.87747, 24.56692, -14.36258],
            [-4.59606, -13.49597, -7.87747, 174.16767, -37.44976, 38.8961],
            [-15.31663, -21.88676, 24.56692, -37.44976, 298.778, -72.07486],
            [0.23793, 3.18351, -14.36258, 38.8961, -72.07486, 128.88083]
        ]),
        "Uu": np.array([
            [9.35543, 3.9644, 0.50826, -13.10381, 2.44425, -2.97145],
            [3.9644, 6.39323, -0.21584, -27.01554, -13.76121, -2.0349],
            [0.50826, -0.21584, 2.51723, -7.20798, 11.64702, -4.46665],
            [-13.10381, -27.01554, -7.20798, 190.46394, 20.31856, 19.82061],
            [2.44425, -13.76121, 11.64702, 20.31856, 210.32551, 7.46356],
            [-2.97145, -2.0349, -4.46665, 19.82061, 7.46356, 50.22318]
        ]),
        "uU": np.array([
            [10.03653, 3.66121, 0.97648, -28.22909, -24.55052, -1.45626],
            [3.66121, 5.54721, 0.50136, -33.97719, -26.17112, -3.89809],
            [0.97648, 0.50136, 1.95707, -11.81432, 10.89475, -6.76246],
            [-28.22909, -33.97719, -11.81432, 353.87718, 63.32121, 67.66717],
            [-24.55052, -26.17112, 10.89475, 63.32121, 303.66708, -41.24858],
            [-1.45626, -3.89809, -6.76246, 67.66717, -41.24858, 56.65976]
        ])
    },
    "tWW": {}
}


def compute_logq(quaternion, pair_type: str, class_type: str):
    sqrt_sum_of_squares = sqrt(
        quaternion[1] ** 2 + quaternion[2] ** 2 + quaternion[3] ** 2
    )
    log = np.array([0.0, 0.0, 0.0])
    if fabs(avg_vector[class_type][pair_type][3]) > fabs(avg_vector[class_type][pair_type][5]):
        if quaternion[1] > 0:
            log[0] = (
                    quaternion[1]
                    * asin(sqrt_sum_of_squares)
                    / sqrt_sum_of_squares
            )
            log[1] = (
                    quaternion[2]
                    * asin(sqrt_sum_of_squares)
                    / sqrt_sum_of_squares
            )
            log[2] = (
                    quaternion[3]
                    * asin(sqrt_sum_of_squares)
                    / sqrt_sum_of_squares
            )
        else:
            log[0] = (
                    quaternion[1]
                    * (asin(sqrt_sum_of_squares) - pi)
                    / sqrt_sum_of_squares
            )
            log[1] = (
                    quaternion[2]
                    * (asin(sqrt_sum_of_squares) - pi)
                    / sqrt_sum_of_squares
            )
            log[2] = (
                    quaternion[3]
                    * (asin(sqrt_sum_of_squares) - pi)
                    / sqrt_sum_of_squares
            )
    else:
        if quaternion[3] > 0:
            log[0] = (
                    quaternion[1]
                    * asin(sqrt_sum_of_squares)
                    / sqrt_sum_of_squares
            )
            log[1] = (
                    quaternion[2]
                    * asin(sqrt_sum_of_squares)
                    / sqrt_sum_of_squares
            )
            log[2] = (
                    quaternion[3]
                    * asin(sqrt_sum_of_squares)
                    / sqrt_sum_of_squares
            )
        else:
            log[0] = (
                    quaternion[1]
                    * (asin(sqrt_sum_of_squares) - pi)
                    / sqrt_sum_of_squares
            )
            log[1] = (
                    quaternion[2]
                    * (asin(sqrt_sum_of_squares) - pi)
                    / sqrt_sum_of_squares
            )
            log[2] = (
                    quaternion[3]
                    * (asin(sqrt_sum_of_squares) - pi)
                    / sqrt_sum_of_squares
            )
    return log


def kabsch_rotate(P, Q):
    """
    Kabsch algorithm for optimal rotation
    P: actual coordinates
    Q: reference coordinates
    """
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    R = np.dot(V, W)
    return R


def get_nucleotides(model):
    # read all residues
    nucleotides = []
    beg = 0
    if len(model.atom) == 0:
        return nucleotides
    i: int
    for i in range(len(model.atom)):
        if (
                model.atom[i].resi != model.atom[i - 1].resi
                or model.atom[i].resn != model.atom[i - 1].resn
        ):
            nucleotides.append(
                nucleotide(
                    model,
                    beg,
                    i,
                    model.atom[i - 1].resi,
                    model.atom[i - 1].resn,
                    model.atom[i - 1].chain,
                )
            )
            beg = i
    nucleotides.append(
        nucleotide(
            model,
            beg,
            len(model.atom),
            model.atom[beg].resi,
            model.atom[beg].resn,
            model.atom[beg].chain,
        )
    )

    # check nucleotides and return filtered list of nucleotides only
    [nuc.check_nucleotide() for nuc in nucleotides]
    return [nucl for nucl in nucleotides if nucl.nuc != "None"]


def HB(don, h, acc, r_cut=3.5, a_cut=50.0):
    # read points and vectors
    p1, p2, p3 = np.array(don), np.array(h), np.array(acc)
    b1 = p2 - p1
    b2 = p3 - p1

    # do vector algebra
    r = np.linalg.norm(b2)
    b1 /= np.linalg.norm(b1)
    b2 /= np.linalg.norm(b2)
    a = np.degrees(np.arccos(np.dot(b1, b2)))

    # report results
    if r < r_cut and a < a_cut:
        return True
    else:
        return False


def cWW(nt1, nt2):
    if nt1.nuc == "A" and nt2.nuc == "A":
        pass
    if nt1.nuc == "A" and nt2.nuc == "C":
        pass
    if nt1.nuc == "A" and nt2.nuc == "G":
        pass
    if nt1.nuc == "A" and (nt2.nuc == "U" or nt2.nuc == "T"):
        if HB(nt1.atom_name["N6"].coord, nt1.RH[2], nt2.atom_name["O4"].coord) and HB(
                nt2.atom_name["N3"].coord, nt2.YH[0], nt1.atom_name["N1"].coord
        ):
            return True
    if nt1.nuc == "C" and nt2.nuc == "A":
        pass
    if nt1.nuc == "C" and nt2.nuc == "C":
        pass
    if nt1.nuc == "C" and nt2.nuc == "G":
        if (
                HB(nt1.atom_name["N4"].coord, nt1.YH[1], nt2.atom_name["O6"].coord)
                and HB(nt2.atom_name["N1"].coord, nt2.RH[3], nt1.atom_name["N3"].coord)
                and HB(nt2.atom_name["N2"].coord, nt2.RH[4], nt1.atom_name["O2"].coord)
        ):
            return True
    if nt1.nuc == "C" and (nt2.nuc == "U" or nt2.nuc == "T"):
        pass
    if nt1.nuc == "G" and nt2.nuc == "A":
        pass
    if nt1.nuc == "G" and nt2.nuc == "C":
        if (
                HB(nt2.atom_name["N4"].coord, nt2.YH[1], nt1.atom_name["O6"].coord)
                and HB(nt1.atom_name["N1"].coord, nt1.RH[3], nt2.atom_name["N3"].coord)
                and HB(nt1.atom_name["N2"].coord, nt1.RH[4], nt2.atom_name["O2"].coord)
        ):
            return True
    if nt1.nuc == "G" and nt2.nuc == "G":
        pass
    if nt1.nuc == "G" and (nt2.nuc == "U" or nt2.nuc == "T"):
        pass
    if (nt1.nuc == "U" or nt1.nuc == "T") and nt2.nuc == "A":
        if HB(nt2.atom_name["N6"].coord, nt2.RH[2], nt1.atom_name["O4"].coord) and HB(
                nt1.atom_name["N3"].coord, nt1.YH[0], nt2.atom_name["N1"].coord
        ):
            return True
    if (nt1.nuc == "U" or nt1.nuc == "T") and nt2.nuc == "C":
        pass
    if (nt1.nuc == "U" or nt1.nuc == "T") and nt2.nuc == "G":
        pass
    if (nt1.nuc == "U" or nt1.nuc == "T") and (nt2.nuc == "U" or nt2.nuc == "T"):
        pass
    return False


def backbone_53(nt1, nt2, min_val=1.0, max_val=2.0):
    if "P" not in nt2.atom_name:
        return False
    dist_sq = sum(
        (nt1.atom_name["O3'"].coord[i] - nt2.atom_name["P"].coord[i])
        * (nt1.atom_name["O3'"].coord[i] - nt2.atom_name["P"].coord[i])
        for i in range(3)
    )
    if min_val * min_val <= dist_sq <= max_val * max_val:
        return True
    else:
        return False


########################################################################################################
#
# class nucleotide
#
########################################################################################################


class nucleotide:
    def __init__(self, model, start, end, id, name, chain: str):
        self.atom = np.array(model.atom[start:end])
        self.id = id
        self.name = name
        self.nuc = "None"
        self.chain = chain
        self.atom_name = {}
        for at in self.atom:
            if at.name not in self.atom_name:
                self.atom_name[at.name] = at
        self.origin = np.array([0.0, 0.0, 0.0])
        self.basis = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        self.ex = self.basis[0]
        self.ey = self.basis[1]
        self.ez = self.basis[2]
        self.RH = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ])
        self.YH = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ])

    def print_full(self):
        for at in self.atom:
            print(at.name, " ", at.resn, " ", at.resi, "total len: ", len(self.atom))
        print("TER")

    def print(self):
        print("residue ", self.id, " type: ", self.nuc, " chain: ", self.chain)

    def print_vectors(self):
        print(f"ex: {self.ex}\ney: {self.ey}\nez: {self.ez}")

    def check_bond(self, name1, name2, min_val=1.0, max_val=2.0):
        dist_sq = 0.0
        for i in range(3):
            dist_sq += (
                               self.atom_name[name1].coord[i] - self.atom_name[name2].coord[i]
                       ) * (self.atom_name[name1].coord[i] - self.atom_name[name2].coord[i])
        if min_val * min_val <= dist_sq <= max_val * max_val:
            return True
        else:
            return False

    def check_nucleotide(self):
        # nucleotides topology
        G_names = ["G", "RG", "DG", "G3", "RG3", "DG3", "G5", "RG5", "DG5"]
        A_names = ["A", "RA", "DA", "A3", "RA3", "DA3", "A5", "RA5", "DA5"]
        C_names = ["C", "RC", "DC", "C3", "RC3", "DC3", "C5", "RC5", "DC5"]
        U_names = ["U", "RU", "U3", "RU3", "U5", "RU5"]
        T_names = ["T", "DT", "T3", "DT3", "T5", "DT5"]
        G_atoms = [
            "O5'",
            "C5'",
            "C4'",
            "C3'",
            "C2'",
            "C1'",
            "O4'",
            "O3'",
            "N9",
            "C8",
            "N7",
            "C5",
            "C6",
            "N1",
            "C2",
            "N3",
            "C4",
            "O6",
            "N2",
        ]
        A_atoms = [
            "O5'",
            "C5'",
            "C4'",
            "C3'",
            "C2'",
            "C1'",
            "O4'",
            "O3'",
            "N9",
            "C8",
            "N7",
            "C5",
            "C6",
            "N1",
            "C2",
            "N3",
            "C4",
            "N6",
        ]
        C_atoms = [
            "O5'",
            "C5'",
            "C4'",
            "C3'",
            "C2'",
            "C1'",
            "O4'",
            "O3'",
            "N1",
            "C6",
            "C5",
            "C4",
            "N3",
            "C2",
            "O2",
            "N4",
        ]
        U_atoms = [
            "O5'",
            "C5'",
            "C4'",
            "C3'",
            "C2'",
            "C1'",
            "O4'",
            "O3'",
            "N1",
            "C6",
            "C5",
            "C4",
            "N3",
            "C2",
            "O2",
            "O4",
        ]
        T_atoms = [
            "O5'",
            "C5'",
            "C4'",
            "C3'",
            "C2'",
            "C1'",
            "O4'",
            "O3'",
            "N1",
            "C6",
            "C5",
            "C4",
            "N3",
            "C2",
            "O2",
            "O4",
            "C7",
        ]
        G_base = [
            "C1'",
            "N9",
            "C8",
            "N7",
            "C5",
            "C6",
            "N1",
            "C2",
            "N3",
            "C4",
            "O6",
            "N2",
        ]
        A_base = ["C1'", "N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4", "N6"]
        C_base = ["C1'", "N1", "C6", "C5", "C4", "N3", "C2", "O2", "N4"]
        U_base = ["C1'", "N1", "C6", "C5", "C4", "N3", "C2", "O2", "O4"]
        T_base = ["C1'", "N1", "C6", "C5", "C4", "N3", "C2", "O2", "O4", "C7"]
        G_bonds = [
            ("C1'", "N9"),
            ("N9", "C8"),
            ("C8", "N7"),
            ("N7", "C5"),
            ("C5", "C6"),
            ("C6", "N1"),
            ("N1", "C2"),
            ("C2", "N3"),
            ("N3", "C4"),
            ("C4", "C5"),
            ("C4", "N9"),
            ("C6", "O6"),
            ("C2", "N2"),
        ]
        A_bonds = [
            ("C1'", "N9"),
            ("N9", "C8"),
            ("C8", "N7"),
            ("N7", "C5"),
            ("C5", "C6"),
            ("C6", "N1"),
            ("N1", "C2"),
            ("C2", "N3"),
            ("N3", "C4"),
            ("C4", "C5"),
            ("C4", "N9"),
            ("C6", "N6"),
        ]
        C_bonds = [
            ("C1'", "N1"),
            ("N1", "C6"),
            ("C6", "C5"),
            ("C5", "C4"),
            ("C4", "N3"),
            ("N3", "C2"),
            ("C2", "N1"),
            ("C2", "O2"),
            ("C4", "N4"),
        ]
        U_bonds = [
            ("C1'", "N1"),
            ("N1", "C6"),
            ("C6", "C5"),
            ("C5", "C4"),
            ("C4", "N3"),
            ("N3", "C2"),
            ("C2", "N1"),
            ("C2", "O2"),
            ("C4", "O4"),
        ]
        T_bonds = [
            ("C1'", "N1"),
            ("N1", "C6"),
            ("C6", "C5"),
            ("C5", "C4"),
            ("C4", "N3"),
            ("N3", "C2"),
            ("C2", "N1"),
            ("C2", "O2"),
            ("C4", "O4"),
            ("C5", "C7"),
        ]
        sugar_bonds = [
            ("O5'", "C5'"),
            ("C5'", "C4'"),
            ("C4'", "C3'"),
            ("C3'", "C2'"),
            ("C2'", "C1'"),
            ("C1'", "O4'"),
            ("O4'", "C4'"),
            ("C3'", "O3'"),
        ]

        # data curation
        if self.name in G_names:
            if (
                    all(names in self.atom_name for names in G_atoms)
                    and all(self.check_bond(name1, name2) for name1, name2 in G_bonds)
                    and all(self.check_bond(name1, name2) for name1, name2 in sugar_bonds)
            ):
                base = np.array(
                    [
                        self.atom_name[name].coord
                        for name in G_base
                        if name in self.atom_name
                    ]
                )
                self.nuc = "G"
        elif self.name in A_names:
            if (
                    all(names in self.atom_name for names in A_atoms)
                    and all(self.check_bond(name1, name2) for name1, name2 in A_bonds)
                    and all(self.check_bond(name1, name2) for name1, name2 in sugar_bonds)
            ):
                base = np.array(
                    [
                        self.atom_name[name].coord
                        for name in A_base
                        if name in self.atom_name
                    ]
                )
                self.nuc = "A"
        elif self.name in C_names:
            if (
                    all(names in self.atom_name for names in C_atoms)
                    and all(self.check_bond(name1, name2) for name1, name2 in C_bonds)
                    and all(self.check_bond(name1, name2) for name1, name2 in sugar_bonds)
            ):
                base = np.array(
                    [
                        self.atom_name[name].coord
                        for name in C_base
                        if name in self.atom_name
                    ]
                )
                self.nuc = "C"
        elif self.name in U_names:
            if (
                    all(names in self.atom_name for names in U_atoms)
                    and all(self.check_bond(name1, name2) for name1, name2 in U_bonds)
                    and all(self.check_bond(name1, name2) for name1, name2 in sugar_bonds)
            ):
                base = np.array(
                    [
                        self.atom_name[name].coord
                        for name in U_base
                        if name in self.atom_name
                    ]
                )
                self.nuc = "U"
        elif self.name in T_names:
            if (
                    all(names in self.atom_name for names in T_atoms)
                    and all(self.check_bond(name1, name2) for name1, name2 in T_bonds)
                    and all(self.check_bond(name1, name2) for name1, name2 in sugar_bonds)
            ):
                base = np.array(
                    [
                        self.atom_name[name].coord
                        for name in T_base
                        if name in self.atom_name
                    ]
                )
                self.nuc = "T"

        # add H sites
        self.add_Hsites()
        # nuclease basis
        if self.nuc != "None":
            self.origin[:] = self.atom_name["C1'"].coord
            ref_coords_nuc = ref_coords[self.nuc]

            actual_coords = np.array([self.atom_name[atom].coord for atom in ref_coords_nuc.keys()])
            reference_coords = np.array(list(ref_coords_nuc.values()))

            actual_centered = actual_coords - np.mean(actual_coords, axis=0)
            reference_centered = reference_coords - np.mean(reference_coords, axis=0)

            R = kabsch_rotate(actual_centered, reference_centered)

            self.ex[:] = R[:, 0]
            self.ey[:] = R[:, 1]
            self.ez[:] = R[:, 2]

    def get_cgo_basis(self, w=0.06, l=0.75, h=0.25, k=1.618):
        from pymol import cgo

        d = w * k
        tl = l + h
        o = self.origin
        x = self.ex
        y = self.ey
        z = self.ez
        cgo_obj = [
            cgo.CYLINDER,
            o[0],
            o[1],
            o[2],
            o[0] + l * x[0],
            o[1] + l * x[1],
            o[2] + l * x[2],
            w,
            1.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            cgo.CYLINDER,
            o[0],
            o[1],
            o[2],
            o[0] + l * y[0],
            o[1] + l * y[1],
            o[2] + l * y[2],
            w,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.0,
            cgo.CYLINDER,
            o[0],
            o[1],
            o[2],
            o[0] + l * z[0],
            o[1] + l * z[1],
            o[2] + l * z[2],
            w,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            cgo.CONE,
            o[0] + l * x[0],
            o[1] + l * x[1],
            o[2] + l * x[2],
            o[0] + tl * x[0],
            o[1] + tl * x[1],
            o[2] + tl * x[2],
            d,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.0,
            cgo.CONE,
            o[0] + l * y[0],
            o[1] + l * y[1],
            o[2] + l * y[2],
            o[0] + tl * y[0],
            o[1] + tl * y[1],
            o[2] + tl * y[2],
            d,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.0,
            1.0,
            1.0,
            cgo.CONE,
            o[0] + l * z[0],
            o[1] + l * z[1],
            o[2] + l * z[2],
            o[0] + tl * z[0],
            o[1] + tl * z[1],
            o[2] + tl * z[2],
            d,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.0,
            1.0,
        ]
        return cgo_obj

    def dihedral(self, at1, at2, at3, at4):
        import numpy as np

        # read points and vectors
        p1, p2, p3, p4 = (
            np.array(self.atom_name[at1].coord),
            np.array(self.atom_name[at2].coord),
            np.array(self.atom_name[at3].coord),
            np.array(self.atom_name[at4].coord),
        )
        b1 = p1 - p2
        b2 = p3 - p2
        b3 = p4 - p3
        b2 /= np.linalg.norm(b2)
        b1 /= np.linalg.norm(b1)
        b3 /= np.linalg.norm(b3)

        # do vector algebra
        n1 = np.cross(b1, b2)
        n2 = np.cross(b3, b2)
        m1 = np.cross(b2, n1)
        x = np.dot(m1, n2)
        y = np.dot(n1, n2)

        # calculate dihedral
        dihedral_angle = np.arctan2(x, y)
        dihedral_angle_degrees = np.degrees(dihedral_angle)
        return dihedral_angle_degrees

    def chi(self):
        if self.nuc == "G" or self.nuc == "A":
            return self.dihedral("O4'", "C1'", "N9", "C4")
        elif self.nuc == "C" or self.nuc == "U" or self.nuc == "T":
            return self.dihedral("O4'", "C1'", "N1", "C2")

    def pucker(self):
        import numpy as np

        tau0 = self.dihedral("C4'", "O4'", "C1'", "C2'")
        tau1 = self.dihedral("O4'", "C1'", "C2'", "C3'")
        tau2 = self.dihedral("C1'", "C2'", "C3'", "C4'")
        tau3 = self.dihedral("C2'", "C3'", "C4'", "O4'")
        tau4 = self.dihedral("C3'", "C4'", "O4'", "C1'")
        P = np.arctan(
            ((tau4 + tau1) - (tau3 + tau0))
            / (2 * tau2 * (np.sin(np.radians(36)) + np.sin(np.radians(72))))
        )
        return np.degrees(P)

    def add_sp2_Hsite(self, at1, at2, at3, r=1.0):

        # read points and vectors
        p1, p2, p3 = (
            np.array(self.atom_name[at1].coord),
            np.array(self.atom_name[at2].coord),
            np.array(self.atom_name[at3].coord),
        )
        b1 = p2 - p1
        b2 = p2 - p3
        b2 /= np.linalg.norm(b2)
        b1 /= np.linalg.norm(b1)

        # do vector algebra
        b3 = b1 + b2
        b3 /= np.linalg.norm(b3)
        b3 *= r
        result = p2 + b3

        return result.tolist()

    def add_sp3_Hsite(self, at1, at2, at3, at4, r=1.0, a=120.0):

        # read points and vectors
        p1, p2, p3, p4 = (
            np.array(self.atom_name[at1].coord),
            np.array(self.atom_name[at2].coord),
            np.array(self.atom_name[at3].coord),
            np.array(self.atom_name[at4].coord),
        )
        b1 = p1 - p2
        b2 = p3 - p2
        b3 = p2 - p4
        b1 /= np.linalg.norm(b1)
        b2 /= np.linalg.norm(b2)
        b3 /= np.linalg.norm(b3)

        # do vector algebra
        d = self.dihedral(at1, at2, at4, at3)
        if d < 0.0:
            d = -90.0 + d / 2
        else:
            d = 90.0 + d / 2
        n1 = np.cross(b1, b3)
        n1 /= np.linalg.norm(n1)
        m1 = np.cross(b3, n1)
        h1 = r * np.cos(np.radians(a)) * b3 + r * np.sin(np.radians(a)) * (
                np.cos(np.radians(d)) * m1 + np.sin(np.radians(d)) * n1
        )
        result = p4 + h1

        return result.tolist()

    def add_Hsites(self):
        if self.nuc == "G":
            self.RH[0][:] = self.add_sp2_Hsite("N9", "C8", "N7")
            self.RH[3][:] = self.add_sp2_Hsite("C6", "N1", "C2")
            self.RH[5][:] = self.add_sp3_Hsite("N1", "C2", "N3", "N2")
            self.RH[4][:] = self.add_sp3_Hsite("N3", "C2", "N1", "N2")
        elif self.nuc == "A":
            self.RH[0][:] = self.add_sp2_Hsite("N9", "C8", "N7")
            self.RH[2][:] = self.add_sp3_Hsite("C5", "C6", "N1", "N6")
            self.RH[1][:] = self.add_sp3_Hsite("N1", "C6", "C5", "N6")
            self.RH[3][:] = self.add_sp2_Hsite("C6", "N1", "C2")
            self.RH[6][:] = self.add_sp2_Hsite("N1", "C2", "N3")
        elif self.nuc == "C":
            self.YH[0][:] = self.add_sp2_Hsite("C2", "N3", "C4")
            self.YH[2][:] = self.add_sp3_Hsite("N3", "C4", "C5", "N4")
            self.YH[1][:] = self.add_sp3_Hsite("C5", "C4", "N3", "N4")
            self.YH[3][:] = self.add_sp2_Hsite("C4", "C5", "C6")
            self.YH[4][:] = self.add_sp2_Hsite("C5", "C6", "N1")
        elif self.nuc == "U":
            self.YH[0][:] = self.add_sp2_Hsite("C2", "N3", "C4")
            self.YH[3][:] = self.add_sp2_Hsite("C4", "C5", "C6")
            self.YH[4][:] = self.add_sp2_Hsite("C5", "C6", "N1")
        elif self.nuc == "T":
            self.YH[0][:] = self.add_sp2_Hsite("C2", "N3", "C4")
            self.YH[4][:] = self.add_sp2_Hsite("C5", "C6", "N1")

    def get_cgo_Hsites(self):
        from pymol import cgo

        cgo_obj = [cgo.COLOR, 0.0, 0.0, 1.0]
        if self.nuc == "G" or self.nuc == "A":
            for hs in self.RH:
                if any(coord != 0.0 for coord in hs):
                    cgo_obj.extend([cgo.SPHERE, hs[0], hs[1], hs[2], 0.2])
        elif self.nuc == "C" or self.nuc == "U" or self.nuc == "T":
            for hs in self.YH:
                if any(coord != 0.0 for coord in hs):
                    cgo_obj.extend([cgo.SPHERE, hs[0], hs[1], hs[2], 0.2])
        return cgo_obj


########################################################################################################
#
# class mutual_pose
#
########################################################################################################


class mutual_pose:
    def __init__(
            self, translation=None, quaternion=None, nucleotide1=None, nucleotide2=None
    ):
        self.classification = None
        if translation is not None and quaternion is not None:
            self.translation = translation
            self.rev_translation = None
            self.quaternion = quaternion
            self.rev_quaternion = None
        elif nucleotide1 is not None and nucleotide2 is not None:
            self.translation = self.calculate_translation(nucleotide1, nucleotide2)
            self.rev_translation = self.calculate_translation(nucleotide2, nucleotide1)
            self.quaternion = self.calculate_quaternion(nucleotide1, nucleotide2)
            self.rev_quaternion = self.calculate_quaternion(nucleotide2, nucleotide1)
            self.nucleotide1_nuc = nucleotide1.nuc
            self.nucleotide2_nuc = nucleotide2.nuc
            self.nucleotide1_chain = nucleotide1.chain
            self.nucleotide2_chain = nucleotide2.chain
            self.nucleotide1_id = nucleotide1.id
            self.nucleotide2_id = nucleotide2.id
        else:
            raise ValueError(
                "Error in mutual_pose constructor: either translation and quaternion or two nucleotides."
            )

    def get_classification(self):
        return self.classification

    def set_classification(self, t_class: str):
        self.classification = t_class

    def calculate_translation(self, nucleotide1, nucleotide2):

        trans = np.array(
            [nucleotide2.origin[i] - nucleotide1.origin[i] for i in range(3)]
        )
        basis = np.array(nucleotide1.basis)
        result = np.dot(basis, trans)
        for coord in result:
            if fabs(coord) > 15.0 or isclose(coord, 0, abs_tol=1e-10):
                return None
        return result

    def calculate_quaternion(self, nucleotide1, nucleotide2):

        # Extract basis set vectors from nucleotide objects
        basis1 = np.array(nucleotide1.basis)
        basis2 = np.array(nucleotide2.basis)

        # Allocate rot. matrix R
        R = np.zeros((3, 3))

        # Calculate R
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    R[i][j] += basis1[i][k] * basis2[j][k]

        # Calculate quaternion from rot. matrix
        q_r = np.sqrt(1 + R[0][0] + R[1][1] + R[2][2]) / 2
        q_x = (R[2][1] - R[1][2]) / (4 * q_r)
        q_y = (R[0][2] - R[2][0]) / (4 * q_r)
        q_z = (R[1][0] - R[0][1]) / (4 * q_r)

        norm = sqrt(q_r ** 2 + q_x ** 2 + q_y ** 2 + q_z ** 2)
        q_r /= norm
        q_x /= norm
        q_y /= norm
        q_z /= norm

        # Return quaternion
        return np.array([q_r, q_x, q_y, q_z])

    def print(self, text):
        if self.translation is None:
            print("Error: translation None.")
        if self.quaternion is None:
            print("Error: quaternion None.")
        print(
            text,
            "translation: ",
            self.translation[0],
            self.translation[1],
            self.translation[2],
            ", quaternion: ",
            self.quaternion[0],
            self.quaternion[1],
            self.quaternion[2],
            self.quaternion[3],
        )

    def find_mahalanobis_distance(self, logq, logq_rev, inv_cov_matrix_t, inv_cov_matrix_rev_t, avg_vector_t,
                                  avg_vector_rev_t):

        real_vectors_rev = np.array(
            [self.rev_translation[0], self.rev_translation[1], self.rev_translation[2], logq_rev[0], logq_rev[1],
             logq_rev[2]])

        real_vectors = np.array(
            [self.translation[0], self.translation[1], self.translation[2], logq[0], logq[1], logq[2]])

        vector_diff = real_vectors - avg_vector_t
        vector_diff_rev = real_vectors_rev - avg_vector_rev_t
        result = np.dot(vector_diff, np.dot(inv_cov_matrix_t, vector_diff.T))
        result_rev = np.dot(vector_diff_rev, np.dot(inv_cov_matrix_rev_t, vector_diff_rev.T))
        return result, result_rev

    def find_pair(self, nucleotide1: nucleotide, nucleotide2: nucleotide):
        chisq_cutoff = 22.45774448  # 99.9 % cutoff
        base_pair = nucleotide1.nuc + nucleotide2.nuc
        if self.translation is None or self.quaternion is None or self.rev_translation is None or self.rev_quaternion is None:
            print("Error: translation or quaternion isn't calculated!.")
            return False
        for class_type in classification_combinations.get(base_pair):
            if class_type != 'cWW':
                continue
            if avg_vector.get(class_type, {}).get(base_pair) is None and \
                    avg_vector[class_type][
                        base_pair[0] + base_pair[1].lower()] is None:
                continue
            if base_pair[0] == base_pair[1]:
                base_pair = base_pair[0] + base_pair[1].lower()
            logq = compute_logq(self.quaternion, base_pair, class_type)
            logq_rev = compute_logq(self.rev_quaternion, base_pair[::-1], class_type)
            inv_cov_matrix_p = np.array(invcov_matrix.get(class_type, {}).get(base_pair))
            inv_cov_matrix_rev_p = np.array(invcov_matrix.get(class_type, {}).get(base_pair[::-1]))
            avg_vector_p = np.array(avg_vector[class_type][base_pair])
            avg_vector_rev_p = np.array(avg_vector[class_type][base_pair[::-1]])
            mahalanobis_distance = self.find_mahalanobis_distance(logq, logq_rev, inv_cov_matrix_p,
                                                                  inv_cov_matrix_rev_p, avg_vector_p, avg_vector_rev_p)
            if mahalanobis_distance is None:
                continue

            if mahalanobis_distance[0] < chisq_cutoff or mahalanobis_distance[1] < chisq_cutoff:
                self.set_classification(class_type)
                print(base_pair, mahalanobis_distance[0], mahalanobis_distance[1])
                return True
        return False

    def get_mutual_pose_cgo(self, ref_nucleotide, w=0.06, l=1.5, h=0.5, k=1.618):
        from pymol import cgo
        import numpy as np

        d = w * k
        tran = np.array(self.translation)
        norm_tr = np.linalg.norm(tran)
        if norm_tr > 2 * (l + h):
            q = 1.0 - h / norm_tr
        else:
            q = 0.75
        origin = ref_nucleotide.origin
        qq = np.array(self.quaternion[1:])
        basis = np.array(ref_nucleotide.basis)
        rot_axis = np.dot(qq, basis).tolist()
        trans_vec = np.dot(tran, basis).tolist()
        cgo_obj = [
            cgo.CYLINDER,
            origin[0],
            origin[1],
            origin[2],
            origin[0] + q * trans_vec[0],
            origin[1] + q * trans_vec[1],
            origin[2] + q * trans_vec[2],
            w,
            1.0,
            1.0,
            0.0,
            1.0,
            1.0,
            0.0,
            cgo.CONE,
            origin[0] + q * trans_vec[0],
            origin[1] + q * trans_vec[1],
            origin[2] + q * trans_vec[2],
            origin[0] + trans_vec[0],
            origin[1] + trans_vec[1],
            origin[2] + trans_vec[2],
            d,
            0.0,
            1.0,
            1.0,
            0.0,
            1.0,
            1.0,
            0.0,
            1.0,
            1.0,
            cgo.CYLINDER,
            origin[0],
            origin[1],
            origin[2],
            origin[0] + 0.75 * rot_axis[0],
            origin[1] + 0.75 * rot_axis[1],
            origin[2] + 0.75 * rot_axis[2],
            w,
            0.0,
            1.0,
            1.0,
            0.0,
            1.0,
            1.0,
            cgo.CONE,
            origin[0] + 0.75 * rot_axis[0],
            origin[1] + 0.75 * rot_axis[1],
            origin[2] + 0.75 * rot_axis[2],
            origin[0] + rot_axis[0],
            origin[1] + rot_axis[1],
            origin[2] + rot_axis[2],
            d,
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            1.0,
            1.0,
            1.0,
            1.0,
        ]
        return cgo_obj

    def get_rot_tran_axis(self, ref_nucleotide):
        # provide rotational translation axis, starting and ending points of translation vector on axis and rotation angle
        # return list origin[x,y,z], end_point[x,y,z], angle

        # read vectors and transform them to canonical basis set
        basis = np.array(ref_nucleotide.basis)
        tran = np.array(self.translation)
        origin = np.array(ref_nucleotide.origin)
        qq = np.array(self.quaternion[1:])
        qq /= np.linalg.norm(qq)
        rot_axis = np.dot(qq, basis)
        trans_vec = np.dot(tran, basis)
        rot_angle_rad = 2 * np.arccos(self.quaternion[0])

        # do vector algebra
        vert_tran_dist = np.dot(rot_axis, trans_vec)
        trans_lateral = trans_vec - vert_tran_dist * rot_axis
        lateral_dist = np.linalg.norm(trans_lateral)
        trans_lateral /= lateral_dist
        to_rot_axis = np.cross(rot_axis, trans_lateral)
        to_rot_axis_dist = lateral_dist * np.sqrt(
            1 / (2 * (1 - np.cos(rot_angle_rad))) - 0.25
        )
        origin_to_axis = (
                0.5 * lateral_dist * trans_lateral + to_rot_axis_dist * to_rot_axis
        )
        start = origin + origin_to_axis
        end = start + vert_tran_dist * rot_axis

        # return values
        return [
            start[0],
            start[1],
            start[2],
            end[0],
            end[1],
            end[2],
            np.degrees(rot_angle_rad),
        ]


class base_pair_info:
    def __init__(
            self, first_strand_name, second_strand_name, first_strand_id, second_strand_id, first_strand_nuc,
            second_strand_nuc, lw_classification=None
    ):
        self.first_strand_name = first_strand_name
        self.second_strand_name = second_strand_name
        self.first_strand_id = first_strand_id
        self.second_strand_id = second_strand_id
        self.first_strand_nuc = first_strand_nuc
        self.second_strand_nuc = second_strand_nuc
        if lw_classification is None:
            self.lw_class = ''
        else:
            self.lw_class = self.parse_classification(lw_classification)
        if self.lw_class is None:
            raise ValueError

    def print(self):
        print(
            f"First strand: {self.first_strand_name} {self.first_strand_id} {self.first_strand_nuc}\n" +
            f"Second strand: {self.second_strand_name} {self.second_strand_id} {self.second_strand_nuc}\n" +
            f"Classification: {self.lw_class}"
        )

    def parse_classification(self, lw_classification):
        orientation = ' '
        if lw_classification[:3] == '-/-' or lw_classification[:3] == '+/+':
            return 'cWW'
        elif lw_classification[-3:] == 'cis':
            orientation = 'c'
        elif lw_classification[-4:] == 'tran':
            orientation = 't'
        else:
            return None

        lw = np.array(['', ''])
        if lw_classification[0] == 'W':
            lw[0] = 'W'
        elif lw_classification[0] == 'H':
            lw[0] = 'H'
        elif lw_classification[0] == 'S':
            lw[0] = 'S'
        else:
            return None

        if lw_classification[2] == 'W':
            lw[1] = 'W'
        elif lw_classification[2] == 'H':
            lw[1] = 'H'
        elif lw_classification[2] == 'S':
            lw[1] = 'S'
        else:
            return None

        return orientation + lw[0] + lw[1]
