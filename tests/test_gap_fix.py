#!/usr/bin/env python
"""
Test: Kdy se vytváří GAP tokeny
"""
from app.services.analysis_service import build_sequence_tokens

# Simulace PDB s přejmenovaným reziduum
pdb_content = """REMARK Test PDB
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C
ATOM      3  N   GAL A   2       2.000   2.000   2.000  1.00  0.00           N
ATOM      4  CA  GAL A   2       3.000   3.000   3.000  1.00  0.00           C
ATOM      5  N   ALA A   3       4.000   4.000   4.000  1.00  0.00           N
ATOM      6  CA  ALA A   3       5.000   5.000   5.000  1.00  0.00           C
END
"""

result = build_sequence_tokens(pdb_content, chain="A", fill_gaps=True)

print("TOKENY PRO ŘETĚZEC A:")
for token in result["chains"]["A"]["tokens"]:
    print(f"  Pos {token['position']:2}: resseq={str(token['resseq']):5} name={token['pdb_resname']:4} is_gap={token['is_gap']} known={token['known']}")

print(f"\nPočet GAP tokenů: {sum(1 for t in result['chains']['A']['tokens'] if t['is_gap'])}")
print(f"Očekávaný počet: 0 (protože všechna rezidua existují v PDB)")

