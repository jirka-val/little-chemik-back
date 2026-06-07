#!/usr/bin/env python
"""
Test: Skutečný GAP - když reziduum OPRAVDU chybí
"""
from app.services.analysis_service import build_sequence_tokens

# PDB bez ALA-2 (skutečná mezera!)
pdb_content = """REMARK Test PDB
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C
ATOM      3  N   ALA A   3       4.000   4.000   4.000  1.00  0.00           N
ATOM      4  CA  ALA A   3       5.000   5.000   5.000  1.00  0.00           C
END
"""

result = build_sequence_tokens(pdb_content, chain="A", fill_gaps=True)

print("TOKENY PRO ŘETĚZEC A (s OPRAVDOVÝM gapem):")
for token in result["chains"]["A"]["tokens"]:
    if token['is_gap']:
        print(f"  Pos {token['position']:2}: resseq={str(token['resseq']):5}  GAP TOKEN ")
    else:
        print(f"  Pos {token['position']:2}: resseq={str(token['resseq']):5} name={token['pdb_resname']:4}")

print(f"\nPočet GAP tokenů: {sum(1 for t in result['chains']['A']['tokens'] if t['is_gap'])}")
print(f"Očekávaný počet: 1 (mezi ALA-1 a ALA-3)")

