
# CLASS AMBER RESIDUE

import math
import bisect

from adams4sims_processing_library.utils import alias

class residue:

 def __init__(self,topology,i):
  first = topology["RESIDUE_POINTER"][i]-1
  if i != len(topology["RESIDUE_LABEL"])-1:
   last = topology["RESIDUE_POINTER"][i+1]-1
  else:
   last = topology["POINTERS"][0]
  self.res_name = topology["RESIDUE_LABEL"][i]
  self.size = last - first
  self.atom_name = topology["ATOM_NAME"][first:last]
  self.charge = topology["CHARGE"][first:last]
  self.mass = topology["MASS"][first:last]
  self.R = []
  self.eps = []
  for n in range(first,last):
   self.R.append(topology["NONBONDED_DIAGONAL_R"][topology["ATOM_TYPE_INDEX"][n]-1])
   self.eps.append(topology["NONBONDED_DIAGONAL_EPS"][topology["ATOM_TYPE_INDEX"][n]-1])
  self.bonds_h, self.bonds_res_h = self.create_bond_dict("BONDS_H",topology,i)
  self.bonds_a, self.bonds_res_a = self.create_bond_dict("BONDS_A",topology,i)
  self.angles_h, self.angles_res_h = self.create_angle_dict("ANGLES_H",topology,i)
  self.angles_a, self.angles_res_a = self.create_angle_dict("ANGLES_A",topology,i)
  self.dihedrals_h, self.dihedrals_res_h = self.create_dihedral_dict("DIHEDRALS_H",topology,i)
  self.dihedrals_a, self.dihedrals_res_a = self.create_dihedral_dict("DIHEDRALS_A",topology,i)
  self.impropers_h, self.impropers_res_h = self.create_dihedral_dict("IMPROPERS_H",topology,i)
  self.impropers_a, self.impropers_res_a = self.create_dihedral_dict("IMPROPERS_A",topology,i)

 def __str__(self):
  print("RESIDUE: ",self.res_name)
  for i in range(self.size):
   print("{} mass {:.2f} charge {:.6f} R {:.4f} eps {:.4f}".format(self.atom_name[i],self.mass[i],self.charge[i]/18.2223,self.R[i],self.eps[i]))
  if len(self.bonds_h) + len(self.bonds_a) != 0:
   print("BONDS:")
  for i in self.bonds_h:
   print("{}-{} {:.1f} {:.3f}".format(i[0],i[1],self.bonds_h[i][0],self.bonds_h[i][1]))
  for i in self.bonds_a:
   print("{}-{} {:.1f} {:.3f}".format(i[0],i[1],self.bonds_a[i][0],self.bonds_a[i][1]))
  if len(self.angles_h) + len(self.angles_a) != 0:
   print("ANGLES:")
  for i in self.angles_h:
   print("{}-{}-{} {:.1f} {:.2f}".format(i[0],i[1],i[2],self.angles_h[i][0],self.angles_h[i][1]))
  for i in self.angles_a:
   print("{}-{}-{} {:.1f} {:.2f}".format(i[0],i[1],i[2],self.angles_a[i][0],self.angles_a[i][1]))
  if len(self.dihedrals_h) + len(self.dihedrals_a) != 0:
   print("DIHEDRALS:")
  for i in self.dihedrals_h:
   print("{}-{}-{}-{}".format(i[0],i[1],i[2],i[3]))
   for j in self.dihedrals_h[i]:
    print("{:.2f} {:.1f} {:.0f}".format(j[0],j[1],j[2]))
  for i in self.dihedrals_a:
   print("{}-{}-{}-{}".format(i[0],i[1],i[2],i[3]))
   for j in self.dihedrals_a[i]:
    print("{:.2f} {:.1f} {:.0f}".format(j[0],j[1],j[2]))
  if len(self.impropers_h) + len(self.impropers_a) != 0:
   print("IMPROPERS:")
  for i in self.impropers_h:
   print("{}-{}-{}-{}".format(i[0],i[1],i[2],i[3]))
   for j in self.impropers_h[i]:
    print("{:.2f} {:.1f} {:.0f}".format(j[0],j[1],j[2]))
  for i in self.impropers_a:
   print("{}-{}-{}-{}".format(i[0],i[1],i[2],i[3]))
   for j in self.impropers_a[i]:
    print("{:.2f} {:.1f} {:.0f}".format(j[0],j[1],j[2]))
  return ""

 def compare(self,topology,i):
  first = topology["RESIDUE_POINTER"][i]-1
  if i != len(topology["RESIDUE_LABEL"])-1:
   last = topology["RESIDUE_POINTER"][i+1]-1
  else:
   last = topology["POINTERS"][0]
  if self.res_name != topology["RESIDUE_LABEL"][i]:
   return False
  if self.size != last - first:
   return False
  if self.atom_name != topology["ATOM_NAME"][first:last]:
   return False
  if self.charge != topology["CHARGE"][first:last]:
   return False
  if self.mass != topology["MASS"][first:last]:
   return False
  R = []
  eps = []
  for n in range(first,last):
   R.append(topology["NONBONDED_DIAGONAL_R"][topology["ATOM_TYPE_INDEX"][n]-1])
   eps.append(topology["NONBONDED_DIAGONAL_EPS"][topology["ATOM_TYPE_INDEX"][n]-1])
  if self.R != R:
   return False
  if self.eps != eps:
   return False
  if self.bonds_h != self.create_bond_dict("BONDS_H",topology,i)[0]:
   return False
  if self.bonds_a != self.create_bond_dict("BONDS_A",topology,i)[0]:
   return False
  if self.angles_h != self.create_angle_dict("ANGLES_H",topology,i)[0]:
   return False
  if self.angles_a != self.create_angle_dict("ANGLES_A",topology,i)[0]:
   return False
  if self.dihedrals_h != self.create_dihedral_dict("DIHEDRALS_H",topology,i)[0]:
   return False
  if self.dihedrals_a != self.create_dihedral_dict("DIHEDRALS_A",topology,i)[0]:
   return False
  if self.impropers_h != self.create_dihedral_dict("IMPROPERS_H",topology,i)[0]:
   return False
  if self.impropers_a != self.create_dihedral_dict("IMPROPERS_A",topology,i)[0]:
   return False
  return True

 def create_bond_dict(self,type,topology,i):
  bond_dict = {}
  bond_res_dict = {}
  for bond in topology[type][i]:
   at1, at2, ind = bond
   res1 = topology["RESIDUE_LABEL"][bisect.bisect_right(topology["RESIDUE_POINTER"],at1+1)-1]
   res2 = topology["RESIDUE_LABEL"][bisect.bisect_right(topology["RESIDUE_POINTER"],at2+1)-1]
   seq_at1 = (topology["ATOM_NAME"][at1],topology["ATOM_NAME"][at2])
   seq_at2 = (topology["ATOM_NAME"][at2],topology["ATOM_NAME"][at1])
   seq_res1 = (res1,res2)
   seq_res2 = (res2,res1)
   (key, value) = min((seq_at1,seq_res1),(seq_at2,seq_res2))
   bond_dict[key] = [topology["BOND_FORCE_CONSTANT"][ind-1],topology["BOND_EQUIL_VALUE"][ind-1]]
   bond_res_dict[key] = value
  return bond_dict, bond_res_dict

 def create_angle_dict(self,type,topology,i):
  angle_dict = {}
  angle_res_dict = {}
  for angle in topology[type][i]:
   at1, at2, at3, ind = angle
   res1 = topology["RESIDUE_LABEL"][bisect.bisect_right(topology["RESIDUE_POINTER"],at1+1)-1]
   res2 = topology["RESIDUE_LABEL"][bisect.bisect_right(topology["RESIDUE_POINTER"],at2+1)-1]
   res3 = topology["RESIDUE_LABEL"][bisect.bisect_right(topology["RESIDUE_POINTER"],at3+1)-1]
   seq_at1 = (topology["ATOM_NAME"][at1],topology["ATOM_NAME"][at2],topology["ATOM_NAME"][at3])
   seq_at2 = (topology["ATOM_NAME"][at3],topology["ATOM_NAME"][at2],topology["ATOM_NAME"][at1])
   seq_res1 = (res1,res2,res3)
   seq_res2 = (res3,res2,res1)
   (key,value) = min((seq_at1,seq_res1),(seq_at2,seq_res2))
   angle_dict[key] = [topology["ANGLE_FORCE_CONSTANT"][ind-1],math.degrees(topology["ANGLE_EQUIL_VALUE"][ind-1])]
   angle_res_dict[key] = value
  return angle_dict, angle_res_dict

 def create_dihedral_dict(self,type,topology,i):
  dihedral_dict = {}
  dihedral_res_dict = {}
  for dihedral in topology[type][i]:
   at1, at2, at3, at4, ind = dihedral
   res1 = topology["RESIDUE_LABEL"][bisect.bisect_right(topology["RESIDUE_POINTER"],at1+1)-1]
   res2 = topology["RESIDUE_LABEL"][bisect.bisect_right(topology["RESIDUE_POINTER"],at2+1)-1]
   res3 = topology["RESIDUE_LABEL"][bisect.bisect_right(topology["RESIDUE_POINTER"],at3+1)-1]
   res4 = topology["RESIDUE_LABEL"][bisect.bisect_right(topology["RESIDUE_POINTER"],at4+1)-1]
   seq_at1 = (topology["ATOM_NAME"][at1],topology["ATOM_NAME"][at2],topology["ATOM_NAME"][at3],topology["ATOM_NAME"][at4])
   seq_at2 = (topology["ATOM_NAME"][at4],topology["ATOM_NAME"][at3],topology["ATOM_NAME"][at2],topology["ATOM_NAME"][at1])
   seq_at3 = (topology["ATOM_NAME"][at1],topology["ATOM_NAME"][at4],topology["ATOM_NAME"][at3],topology["ATOM_NAME"][at2])
   seq_at4 = (topology["ATOM_NAME"][at2],topology["ATOM_NAME"][at1],topology["ATOM_NAME"][at3],topology["ATOM_NAME"][at4])
   seq_at5 = (topology["ATOM_NAME"][at2],topology["ATOM_NAME"][at4],topology["ATOM_NAME"][at3],topology["ATOM_NAME"][at1])
   seq_at6 = (topology["ATOM_NAME"][at4],topology["ATOM_NAME"][at2],topology["ATOM_NAME"][at3],topology["ATOM_NAME"][at1])
   seq_at7 = (topology["ATOM_NAME"][at4],topology["ATOM_NAME"][at1],topology["ATOM_NAME"][at3],topology["ATOM_NAME"][at2])
   seq_res1 = (res1,res2,res3,res4)
   seq_res2 = (res4,res3,res2,res1)
   seq_res3 = (res1,res4,res3,res2)
   seq_res4 = (res2,res1,res3,res4)
   seq_res5 = (res2,res4,res3,res1)
   seq_res6 = (res4,res2,res3,res1)
   seq_res7 = (res4,res1,res3,res2)
   if type.startswith("DIHEDRALS"):
    (key,value) = min((seq_at1,seq_res1),(seq_at2,seq_res2))
   else:
    (key,value) = min((seq_at1,seq_res1),(seq_at3,seq_res3),(seq_at4,seq_res4),(seq_at5,seq_res5),(seq_at6,seq_res6),(seq_at7,seq_res7))
   new_entry = [topology["DIHEDRAL_FORCE_CONSTANT"][ind-1],math.degrees(topology["DIHEDRAL_PHASE"][ind-1]),topology["DIHEDRAL_PERIODICITY"][ind-1]]
   if key not in dihedral_dict:
    dihedral_dict[key] = [new_entry]
    dihedral_res_dict[key] = value
   elif not any(entry[2] == new_entry[2] for entry in dihedral_dict[key]):
    dihedral_dict[key].append(new_entry)
  for key in dihedral_dict:
   dihedral_dict[key].sort(key=lambda x: x[2])
  return dihedral_dict, dihedral_res_dict

 def amber2gromacs_units(self):
  R_coeff = 2**(-1.0/6.0)/5.0
  E_coeff = 4.184
  Q_coeff = 18.2223
  self.R = [x*R_coeff for x in self.R]
  self.eps = [x*E_coeff for x in self.eps]
  self.charge = [x/Q_coeff for x in self.charge]
  self.bonds_h = {key: [value[0]*E_coeff*200.0, value[1]*0.1] for key, value in self.bonds_h.items()}
  self.bonds_a = {key: [value[0]*E_coeff*200.0, value[1]*0.1] for key, value in self.bonds_a.items()}
  self.angles_h = {key: [value[0]*E_coeff*2.0, value[1]] for key, value in self.angles_h.items()}
  self.angles_a = {key: [value[0]*E_coeff*2.0, value[1]] for key, value in self.angles_a.items()}
  self.dihedrals_h = {key: [[sublist[0]*E_coeff, sublist[1], sublist[2]] for sublist in value] for key, value in self.dihedrals_h.items()}
  self.dihedrals_a = {key: [[sublist[0]*E_coeff, sublist[1], sublist[2]] for sublist in value] for key, value in self.dihedrals_a.items()}
  self.impropers_h = {key: [[sublist[0]*E_coeff, sublist[1], sublist[2]] for sublist in value] for key, value in self.impropers_h.items()}
  self.impropers_a = {key: [[sublist[0]*E_coeff, sublist[1], sublist[2]] for sublist in value] for key, value in self.impropers_a.items()}

 def checkFF(self,ff,mol):
  # CHECK CHARGES
  if any(abs(q_res - q_ff) > 0.00001 for q_res, q_ff in zip(self.charge,ff.units[alias.resn_alias(self.res_name)]["atoms"]["charge"])):
   return False
  # CHECK R
  if any(abs(R_res - R_ff) > 0.000001 for R_res, R_ff in zip(self.R,ff.units[alias.resn_alias(self.res_name)]["atoms"]["R"])):
   return False
  # CHECK Eps
  if any(abs(eps_res - eps_ff) > 0.00001 for eps_res, eps_ff in zip(self.eps,ff.units[alias.resn_alias(self.res_name)]["atoms"]["eps"])):
   return False
  # CHECK BONDS
  for bond in self.bonds_h:
   bond_type = tuple(ff.find_atom_type(alias.name_alias(mol,name),alias.resn_alias(resn)) for name,resn in zip(bond,self.bonds_res_h[bond]))
   key = min((bond_type[0],bond_type[1]),(bond_type[1],bond_type[0]))
   if abs(self.bonds_h[bond][0]-ff.b["bondtypes"][key][0]) > 0.01 or abs(self.bonds_h[bond][1]-ff.b["bondtypes"][key][1]) > 0.0001:
    return False
  for bond in self.bonds_a:
   bond_type = tuple(ff.find_atom_type(alias.name_alias(mol,name),alias.resn_alias(resn)) for name,resn in zip(bond,self.bonds_res_a[bond]))
   key = min((bond_type[0],bond_type[1]),(bond_type[1],bond_type[0]))
   if abs(self.bonds_a[bond][0]-ff.b["bondtypes"][key][0]) > 0.01 or abs(self.bonds_a[bond][1]-ff.b["bondtypes"][key][1]) > 0.0001:
    return False
  # CHECK ANGLES
  for angle in self.angles_h:
   angle_type = tuple(ff.find_atom_type(alias.name_alias(mol,name),alias.resn_alias(resn)) for name,resn in zip(angle,self.angles_res_h[angle]))
   key = min((angle_type[0],angle_type[1],angle_type[2]),(angle_type[2],angle_type[1],angle_type[0]))
   if abs(self.angles_h[angle][0]-ff.b["angletypes"][key][0]) > 0.001 or abs(self.angles_h[angle][1]-ff.b["angletypes"][key][1]) > 0.001:
    return False
  for angle in self.angles_a:
   angle_type = tuple(ff.find_atom_type(alias.name_alias(mol,name),alias.resn_alias(resn)) for name,resn in zip(angle,self.angles_res_a[angle]))
   key = min((angle_type[0],angle_type[1],angle_type[2]),(angle_type[2],angle_type[1],angle_type[0]))
   if abs(self.angles_a[angle][0]-ff.b["angletypes"][key][0]) > 0.001 or abs(self.angles_a[angle][1]-ff.b["angletypes"][key][1]) > 0.001:
    return False
  # CHECK DIHEDRALS
  for dihedral in self.dihedrals_h:
   dihedral_type = tuple(ff.find_atom_type(alias.name_alias(mol,name),alias.resn_alias(resn)) for name,resn in zip(dihedral,self.dihedrals_res_h[dihedral]))
   key = min((dihedral_type[0],dihedral_type[1],dihedral_type[2],dihedral_type[3]),(dihedral_type[3],dihedral_type[2],dihedral_type[1],dihedral_type[0]))
   ff_filtered = [row for row in ff.b["dihedraltypes"][key] if row[0] != 0.0]
   res_filtered = [row for row in self.dihedrals_h[dihedral] if row[0] != 0.0]
   if len(ff_filtered) != len(res_filtered):
    return False
#   if not all(abs(ff_per[0] - res_per[0])<0.00001 and abs(ff_per[1] - res_per[1])<0.1 and ff_per[2]==res_per[2] for ff_per, res_per in zip(ff_filtered,res_filtered)):
   if not all(((abs(ff_per[0] - res_per[0])<0.00001 and abs(ff_per[1] - res_per[1])<0.1) or (abs(ff_per[0] + res_per[0])<0.00001 and abs((ff_per[1] - res_per[1]) % 360 - 180)<0.1)) and ff_per[2]==res_per[2] for ff_per, res_per in zip(ff_filtered,res_filtered)):
    return False
  for dihedral in self.dihedrals_a:
   dihedral_type = tuple(ff.find_atom_type(alias.name_alias(mol,name),alias.resn_alias(resn)) for name,resn in zip(dihedral,self.dihedrals_res_a[dihedral]))
   key = min((dihedral_type[0],dihedral_type[1],dihedral_type[2],dihedral_type[3]),(dihedral_type[3],dihedral_type[2],dihedral_type[1],dihedral_type[0]))
   ff_filtered = [row for row in ff.b["dihedraltypes"][key] if row[0] != 0.0]
   res_filtered = [row for row in self.dihedrals_a[dihedral] if row[0] != 0.0]
   if len(ff_filtered) != len(res_filtered):
    return False
#   if not all(abs(ff_per[0] - res_per[0])<0.00001 and abs(ff_per[1] - res_per[1])<0.1 and ff_per[2]==res_per[2] for ff_per, res_per in zip(ff_filtered,res_filtered)):
   if not all(((abs(ff_per[0] - res_per[0])<0.00001 and abs(ff_per[1] - res_per[1])<0.1) or (abs(ff_per[0] + res_per[0])<0.00001 and abs((ff_per[1] - res_per[1]) % 360 - 180)<0.1)) and ff_per[2]==res_per[2] for ff_per, res_per in zip(ff_filtered,res_filtered)):
    return False
  # CHECK IMPROPERS
  for improper in self.impropers_h:
   it = tuple(ff.find_atom_type(alias.name_alias(mol,name),alias.resn_alias(resn)) for name,resn in zip(improper,self.impropers_res_h[improper]))
   key = min((it[0],it[1],it[2],it[3]),(it[1],it[0],it[2],it[3]),(it[0],it[3],it[2],it[1]),(it[3],it[0],it[2],it[1]),(it[1],it[3],it[2],it[0]),(it[3],it[1],it[2],it[0]))
   ff_filtered = [row for row in ff.b["impropertypes"][key] if row[0] != 0.0]
   res_filtered = [row for row in self.impropers_h[improper] if row[0] != 0.0]
   if len(ff_filtered) != len(res_filtered):
    return False
#   if not all(abs(ff_per[0] - res_per[0])<0.00001 and abs(ff_per[1] - res_per[1])<0.1 and ff_per[2]==res_per[2] for ff_per, res_per in zip(ff_filtered,res_filtered)):
   if not all(((abs(ff_per[0] - res_per[0])<0.00001 and abs(ff_per[1] - res_per[1])<0.1) or (abs(ff_per[0] + res_per[0])<0.00001 and abs((ff_per[1] - res_per[1]) % 360 - 180)<0.1)) and ff_per[2]==res_per[2] for ff_per, res_per in zip(ff_filtered,res_filtered)):
    return False
  for improper in self.impropers_a:
   it = tuple(ff.find_atom_type(alias.name_alias(mol,name),alias.resn_alias(resn)) for name,resn in zip(improper,self.impropers_res_a[improper]))
   key = min((it[0],it[1],it[2],it[3]),(it[1],it[0],it[2],it[3]),(it[0],it[3],it[2],it[1]),(it[3],it[0],it[2],it[1]),(it[1],it[3],it[2],it[0]),(it[3],it[1],it[2],it[0]))
   ff_filtered = [row for row in ff.b["impropertypes"][key] if row[0] != 0.0]
   res_filtered = [row for row in self.impropers_a[improper] if row[0] != 0.0]
   if len(ff_filtered) != len(res_filtered):
    return False
#   if not all(abs(ff_per[0] - res_per[0])<0.00001 and abs(ff_per[1] - res_per[1])<0.1 and ff_per[2]==res_per[2] for ff_per, res_per in zip(ff_filtered,res_filtered)):
   if not all(((abs(ff_per[0] - res_per[0])<0.00001 and abs(ff_per[1] - res_per[1])<0.1) or (abs(ff_per[0] + res_per[0])<0.00001 and abs((ff_per[1] - res_per[1]) % 360 - 180)<0.1)) and ff_per[2]==res_per[2] for ff_per, res_per in zip(ff_filtered,res_filtered)):
    return False
  return True

