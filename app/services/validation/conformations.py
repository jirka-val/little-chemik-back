from typing import List, Dict, Any

class ConformationManager:
    @staticmethod
    def detect_alt_locs(pdb_content: str) -> List[Dict[str, Any]]:
        res_map = {}
        for line in pdb_content.splitlines():
            if line.startswith(("ATOM", "HETATM")) and len(line) > 16:
                alt_id = line[16].strip()
                if alt_id:
                    key = f"{line[21].strip()}-{line[22:26].strip()}-{line[17:20].strip()}"
                    if key not in res_map:
                        res_map[key] = {
                            "chain": line[21].strip(),
                            "res_id": line[22:26].strip(),
                            "res_name": line[17:20].strip(),
                            "variants": set()
                        }
                    res_map[key]["variants"].add(alt_id)
        return [v for v in res_map.values() if len(v["variants"]) > 1]

    @staticmethod
    def filter_pdb_by_selection(pdb_content: str, selections: Dict[str, str]) -> str:
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