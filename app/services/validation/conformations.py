from typing import List, Dict, Any

class ConformationManager:
    @staticmethod
    def detect_alt_locs(pdb_content: str) -> List[Dict[str, Any]]:
        """Identifikuje rezidua s alternativními konformacemi pro frontend."""
        res_map = {}
        for line in pdb_content.splitlines():
            if line.startswith(("ATOM", "HETATM")) and len(line) > 16:
                alt_id = line[16].strip()
                if alt_id:
                    res_name, res_id, chain = line[17:20].strip(), line[22:26].strip(), line[21].strip()
                    key = f"{chain}-{res_id}-{res_name}"
                    if key not in res_map:
                        res_map[key] = {"chain": chain, "res_id": res_id, "res_name": res_name, "variants": set()}
                    res_map[key]["variants"].add(alt_id)
        return [v for v in res_map.values() if len(v["variants"]) > 1]

    @staticmethod
    def filter_pdb_by_selection(pdb_content: str, selections: Dict[str, str]) -> str:
        """Vyfiltruje PDB a ponechá pouze vybrané konformace reziduí."""
        output = []
        for line in pdb_content.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                alt_id = line[16].strip()
                if not alt_id:
                    output.append(line)
                    continue
                key = f"{line[21].strip()}-{line[22:26].strip()}-{line[17:20].strip()}"
                if key in selections:
                    if alt_id == selections[key]:
                        output.append(line[:16] + " " + line[17:])
                elif alt_id == 'A':
                    output.append(line[:16] + " " + line[17:])
            else:
                output.append(line)
        return "\n".join(output)