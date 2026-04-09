
# FUNCTIONS FOR HANDLING AMBER TOPOLOGY FILES

import re
import math
import bisect
from itertools import combinations
from collections import deque
from datetime import datetime

from app.utils.adams4sims_processing_library.utils.conversion_dict import get_conversion_dict
from app.utils.adams4sims_processing_library.utils import alias

# Function to parse line in AMBER topology
def parse_line(line, format_spec):
 data_type = format_spec[0]
 if '.' in format_spec:
  item_length = int(format_spec[1:].split('.')[0])
 else:
  item_length = int(format_spec[1:])

 items = [line[i:i+item_length] for i in range(0, len(line), item_length)]

 if data_type == 'I':
  return [int(item) for item in items]
 elif data_type == 'a':
  return [item.strip() for item in items]
 elif data_type == 'E':
  return [float(item) for item in items]
 else:
  raise ValueError(f"Unsupported format: {data_type}")

# Function to read block of data
def block_of_data(file_name, name):
 try:
  with open(file_name, 'r') as file:
   found_flag = False
   format_spec = None
   data = []

   for line in file:
    if found_flag:
     if line.startswith("%FLAG"):
      break
     if format_spec:
      parsed_data = parse_line(line.rstrip("\n"), format_spec)
      data.extend(parsed_data)
    elif line.startswith(f"%FLAG {name}"):
     found_flag = True
     # Read the next line for the format
     format_line = next(file).strip()
     match = re.match(r"%FORMAT\((\d+)([aIE]\d+(\.\d+)?)\)", format_line)
     if match:
      num_items, format_spec = match.groups()[0], match.groups()[1]
      format_spec = f"{format_spec}"
     else:
      raise ValueError("Invalid format specification")

   if not found_flag:
    print(f"No section starting with '%FLAG {name}' found in the file.")
    return []
   return data
 except FileNotFoundError:
  print(f"The file '{file_name}' does not exist.")
  return []
 except Exception as e:
  print(f"An error occurred: {e}")
  return []

# Function reading complete topology
def read_AMBER_topology(file_name: str) -> dict:
 """
     Reads and parses an AMBER topology file.

     Parameters
     ----------
     file_name : str
         The path to the AMBER topology file to be read.

     Returns
     -------
     dict
         A dictionary containing the parsed topology data. The keys represent
         different sections of the topology file, and the values are the corresponding
         data arrays.

     Raises
     ------
     FileNotFoundError
         If the specified file does not exist.
     ValueError
         If the format specification in the topology file is invalid.
     Exception
         If any other error occurs during the parsing process.
     """
 try:
  with open(file_name, 'r') as file:
   format_spec = None
   flag = None
   topology = {}
   data = []

   for line in file:
    if line.startswith("%FLAG"):
     if flag:
      topology[flag] = data
      data = []
     flag = line.split()[1]
     format_line = next(file).strip()
     match = re.match(r"%FORMAT\((\d+)([aIE]\d+(\.\d+)?)\)", format_line)
     if match:
      num_items, format_spec = match.groups()[0], match.groups()[1]
      format_spec = f"{format_spec}"
     else:
      raise ValueError("Invalid format specification")
    elif format_spec:
     parsed_data = parse_line(line.rstrip("\n"), format_spec)
     data.extend(parsed_data)
   if flag not in topology:
    topology[flag] = data
   check_AMBER_topology(topology)
   return topology
 except FileNotFoundError:
  print(f"The file '{file_name}' does not exist.")
  return {}
 except Exception as e:
  print(f"An error occurred: {e}")
  return {}

# Function for checking AMBER topology
def check_AMBER_topology(topology):
 checks = [
    (["ATOM_NAME","CHARGE","ATOMIC_NUMBER","MASS","ATOM_TYPE_INDEX","NUMBER_EXCLUDED_ATOMS","AMBER_ATOM_TYPE","TREE_CHAIN_CLASSIFICATION","JOIN_ARRAY","IROTAT","RADII","SCREEN"],topology["POINTERS"][0],"atoms"),
    (["RESIDUE_LABEL","RESIDUE_POINTER"],topology["POINTERS"][11],"residues"),
    (["NONBONDED_PARM_INDEX"],topology["POINTERS"][1]**2,"nonbond pairs"),
    (["LENNARD_JONES_ACOEF","LENNARD_JONES_BCOEF"],topology["POINTERS"][1]*(topology["POINTERS"][1]+1)//2,"Lennard-Jones parameters"),
    (["BOND_FORCE_CONSTANT","BOND_EQUIL_VALUE"],topology["POINTERS"][15],"bond parameters"),
    (["ANGLE_FORCE_CONSTANT","ANGLE_EQUIL_VALUE"],topology["POINTERS"][16],"angle parameters"),
    (["DIHEDRAL_FORCE_CONSTANT","DIHEDRAL_PERIODICITY","DIHEDRAL_PHASE","SCEE_SCALE_FACTOR","SCNB_SCALE_FACTOR"],topology["POINTERS"][17],"dihedral parameters"),
    (["BONDS_INC_HYDROGEN"],3*topology["POINTERS"][2],"bonds with hydrogen"),
    (["BONDS_WITHOUT_HYDROGEN"],3*topology["POINTERS"][12],"bonds without hydrogen"),
    (["ANGLES_INC_HYDROGEN"],4*topology["POINTERS"][4],"angles with hydrogen"),
    (["ANGLES_WITHOUT_HYDROGEN"],4*topology["POINTERS"][13],"angles without hydrogen"),
    (["DIHEDRALS_INC_HYDROGEN"],5*topology["POINTERS"][6],"dihedrals with hydrogen"),
    (["DIHEDRALS_WITHOUT_HYDROGEN"],5*topology["POINTERS"][14],"dihedrals without hydrogen"),
    (["SOLTY"],topology["POINTERS"][18],"solty parameters"),
    (["EXCLUDED_ATOMS_LIST"],topology["POINTERS"][10],"excluded atom list"),
    (["HBOND_ACOEF","HBOND_BCOEF","HBCUT"],topology["POINTERS"][19],"HB parameters")
 ]

 if len(topology["POINTERS"])<31:
  raise Exception("Inconsistent number of POINTERS in AMBER topology file.")
 for keys, expected_length, description in checks:
  for key in keys:
   if len(topology[key]) != expected_length:
    raise Exception(f"Inconsistent number of {description} in AMBER topology file for {key}.")

# Function for adding atomic nonbonded parameters and NBfix
def add_atomic_nonbond_parm(topology):
 ntypes = topology["POINTERS"][1]
 nonbonded_diagonal_index = []
 nonbonded_diagoval_R = []
 nonbonded_diagoval_eps = []
 nonbonded_NBfix = []

 for i in range(ntypes):
  nonbonded_diagonal_index.append(topology["NONBONDED_PARM_INDEX"][ntypes*i+i]-1)
  if topology["LENNARD_JONES_ACOEF"][nonbonded_diagonal_index[i]] == 0.0 or topology["LENNARD_JONES_BCOEF"][nonbonded_diagonal_index[i]] == 0.0:
   nonbonded_diagoval_R.append(0.0)
   nonbonded_diagoval_eps.append(0.0)
  else:
   nonbonded_diagoval_R.append(0.5*(2*topology["LENNARD_JONES_ACOEF"][nonbonded_diagonal_index[i]]/topology["LENNARD_JONES_BCOEF"][nonbonded_diagonal_index[i]])**(1.0/6.0))
   nonbonded_diagoval_eps.append(0.25*topology["LENNARD_JONES_BCOEF"][nonbonded_diagonal_index[i]]*topology["LENNARD_JONES_BCOEF"][nonbonded_diagonal_index[i]]/topology["LENNARD_JONES_ACOEF"][nonbonded_diagonal_index[i]])

 for i in range(ntypes):
  for j in range(i):
   index = topology["NONBONDED_PARM_INDEX"][ntypes*i+j]-1
   if index != topology["NONBONDED_PARM_INDEX"][ntypes*j+i]-1:
    raise Exception("Nonbonded matrix in topology file is not symmetric.")
   if index >=0:
    if topology["LENNARD_JONES_ACOEF"][index] == 0.0 or topology["LENNARD_JONES_BCOEF"][index] == 0.0:
     if nonbonded_diagoval_R[i] != 0.0 and nonbonded_diagoval_R[j] != 0.0 and nonbonded_diagoval_eps[i] != 0.0 and nonbonded_diagoval_eps[j] != 0.0:
      nonbonded_NBfix.append((i,j))
    else:
     R = (2*topology["LENNARD_JONES_ACOEF"][index]/topology["LENNARD_JONES_BCOEF"][index])**(1.0/6.0)
     eps = 0.25*topology["LENNARD_JONES_BCOEF"][index]*topology["LENNARD_JONES_BCOEF"][index]/topology["LENNARD_JONES_ACOEF"][index]
     if abs(R-nonbonded_diagoval_R[i]-nonbonded_diagoval_R[j])>0.0001 or abs(math.sqrt(nonbonded_diagoval_eps[i]*nonbonded_diagoval_eps[j])-eps)>0.0001:
      nonbonded_NBfix.append((i,j))
 topology["NONBONDED_DIAGONAL_R"] = nonbonded_diagoval_R
 topology["NONBONDED_DIAGONAL_EPS"] = nonbonded_diagoval_eps
 topology["NONBONDED_NBFIX"] = nonbonded_NBfix

# Function for the splitting of the bonded parameters into residues
def split_bonded_parm_into_residues(topology):
 topology["BONDS_H"] = [[] for _ in range(topology["POINTERS"][11])]
 for i in range(0,len(topology["BONDS_INC_HYDROGEN"]),3):
  at1 = topology["BONDS_INC_HYDROGEN"][i]//3
  at2 = topology["BONDS_INC_HYDROGEN"][i+1]//3
  ind = topology["BONDS_INC_HYDROGEN"][i+2]
  res1 = bisect.bisect_right(topology["RESIDUE_POINTER"],at1+1)-1
  res2 = bisect.bisect_right(topology["RESIDUE_POINTER"],at2+1)-1
  topology["BONDS_H"][max(res1,res2)].append([at1, at2, ind])
  #topology["BONDS_H"][res1].append([at1, at2, ind])
  #if res2 != res1:
  # topology["BONDS_H"][res2].append([at1, at2, ind])
 topology["BONDS_A"] = [[] for _ in range(topology["POINTERS"][11])]
 for i in range(0,len(topology["BONDS_WITHOUT_HYDROGEN"]),3):
  at1 = topology["BONDS_WITHOUT_HYDROGEN"][i]//3
  at2 = topology["BONDS_WITHOUT_HYDROGEN"][i+1]//3
  ind = topology["BONDS_WITHOUT_HYDROGEN"][i+2]
  res1 = bisect.bisect_right(topology["RESIDUE_POINTER"],at1+1)-1
  res2 = bisect.bisect_right(topology["RESIDUE_POINTER"],at2+1)-1
  topology["BONDS_A"][max(res1,res2)].append([at1, at2, ind])
  #topology["BONDS_A"][res1].append([at1, at2, ind])
  #if res2 != res1:
  # topology["BONDS_A"][res2].append([at1, at2, ind])
 topology["ANGLES_H"] = [[] for _ in range(topology["POINTERS"][11])]
 for i in range(0,len(topology["ANGLES_INC_HYDROGEN"]),4):
  at1 = topology["ANGLES_INC_HYDROGEN"][i]//3
  at2 = topology["ANGLES_INC_HYDROGEN"][i+1]//3
  at3 = topology["ANGLES_INC_HYDROGEN"][i+2]//3
  ind = topology["ANGLES_INC_HYDROGEN"][i+3]
  res1 = bisect.bisect_right(topology["RESIDUE_POINTER"],at1+1)-1
  res2 = bisect.bisect_right(topology["RESIDUE_POINTER"],at2+1)-1
  res3 = bisect.bisect_right(topology["RESIDUE_POINTER"],at3+1)-1
  topology["ANGLES_H"][res2].append([at1, at2, at3, ind])
  #if res2 != res1:
  # topology["ANGLES_H"][res2].append([at1, at2, at3, ind])
 topology["ANGLES_A"] = [[] for _ in range(topology["POINTERS"][11])]
 for i in range(0,len(topology["ANGLES_WITHOUT_HYDROGEN"]),4):
  at1 = topology["ANGLES_WITHOUT_HYDROGEN"][i]//3
  at2 = topology["ANGLES_WITHOUT_HYDROGEN"][i+1]//3
  at3 = topology["ANGLES_WITHOUT_HYDROGEN"][i+2]//3
  ind = topology["ANGLES_WITHOUT_HYDROGEN"][i+3]
  res1 = bisect.bisect_right(topology["RESIDUE_POINTER"],at1+1)-1
  res2 = bisect.bisect_right(topology["RESIDUE_POINTER"],at2+1)-1
  res3 = bisect.bisect_right(topology["RESIDUE_POINTER"],at3+1)-1
  topology["ANGLES_A"][res2].append([at1, at2, at3, ind])
  #if res2 != res1: 
  # topology["ANGLES_A"][res2].append([at1, at2, at3, ind])
 topology["DIHEDRALS_H"] = [[] for _ in range(topology["POINTERS"][11])]
 topology["IMPROPERS_H"] = [[] for _ in range(topology["POINTERS"][11])]
 for i in range(0,len(topology["DIHEDRALS_INC_HYDROGEN"]),5):
  at1 = topology["DIHEDRALS_INC_HYDROGEN"][i]//3
  at2 = topology["DIHEDRALS_INC_HYDROGEN"][i+1]//3
  at3 = abs(topology["DIHEDRALS_INC_HYDROGEN"][i+2]//3)
  at4 = topology["DIHEDRALS_INC_HYDROGEN"][i+3]//3
  ind = topology["DIHEDRALS_INC_HYDROGEN"][i+4]
  res1 = bisect.bisect_right(topology["RESIDUE_POINTER"],at1+1)-1
  res2 = bisect.bisect_right(topology["RESIDUE_POINTER"],at2+1)-1
  res3 = bisect.bisect_right(topology["RESIDUE_POINTER"],at3+1)-1
  res4 = bisect.bisect_right(topology["RESIDUE_POINTER"],abs(at4)+1)-1
  if at4 >= 0: 
   topology["DIHEDRALS_H"][max(res2,res3)].append([at1, at2, at3, at4, ind])
   #topology["DIHEDRALS_H"][res2].append([at1, at2, at3, at4, ind])
   #if res3 != res2:
   # topology["DIHEDRALS_H"][res3].append([at1, at2, at3, at4, ind])
  else:
   topology["IMPROPERS_H"][res3].append([at1, at2, at3, -at4, ind])
   #topology["IMPROPERS_H"][res1].append([at1, at2, at3, -at4, ind])
   #if res4 != res1:
   # topology["IMPROPERS_H"][res4].append([at1, at2, at3, -at4, ind])
   #if res2 != res1 and res2 != res4:
   # topology["IMPROPERS_H"][res2].append([at1, at2, at3, -at4, ind])
 topology["DIHEDRALS_A"] = [[] for _ in range(topology["POINTERS"][11])]
 topology["IMPROPERS_A"] = [[] for _ in range(topology["POINTERS"][11])]
 for i in range(0,len(topology["DIHEDRALS_WITHOUT_HYDROGEN"]),5):
  at1 = topology["DIHEDRALS_WITHOUT_HYDROGEN"][i]//3
  at2 = topology["DIHEDRALS_WITHOUT_HYDROGEN"][i+1]//3
  at3 = abs(topology["DIHEDRALS_WITHOUT_HYDROGEN"][i+2]//3)
  at4 = topology["DIHEDRALS_WITHOUT_HYDROGEN"][i+3]//3
  ind = topology["DIHEDRALS_WITHOUT_HYDROGEN"][i+4]
  res1 = bisect.bisect_right(topology["RESIDUE_POINTER"],at1+1)-1
  res2 = bisect.bisect_right(topology["RESIDUE_POINTER"],at2+1)-1
  res3 = bisect.bisect_right(topology["RESIDUE_POINTER"],at3+1)-1
  res4 = bisect.bisect_right(topology["RESIDUE_POINTER"],abs(at4)+1)-1
  # proline exception
  #resn1 = alias.resn_alias(topology["RESIDUE_LABEL"][res1]).strip()
  #resn2 = alias.resn_alias(topology["RESIDUE_LABEL"][res2]).strip()
  #resn3 = alias.resn_alias(topology["RESIDUE_LABEL"][res3]).strip()
  #resn4 = alias.resn_alias(topology["RESIDUE_LABEL"][res4]).strip()
  #atmn1 = alias.name_alias(None,topology["ATOM_NAME"][at1]).strip()
  #atmn2 = alias.name_alias(None,topology["ATOM_NAME"][at2]).strip()
  #atmn3 = alias.name_alias(None,topology["ATOM_NAME"][at3]).strip()
  #atmn4 = alias.name_alias(None,topology["ATOM_NAME"][abs(at4)]).strip()
  if at4 >= 0:
   #topology["DIHEDRALS_A"][res2].append([at1, at2, at3, at4, ind])
   #if res3 != res2:
   # topology["DIHEDRALS_A"][res3].append([at1, at2, at3, at4, ind])
   #target_residues = [res2]
   #if res3 != res2:
   # target_residues.append(res3)
   #if res3 != res2 and ((resn1 == "PRO" and atmn1 == "CD") or (resn4 == "PRO" and atmn4 == "CD")):
   # pro_targets = []
   # if resn2 == "PRO":
   #  pro_targets.append(res2)
   # if resn3 == "PRO":
   #  pro_targets.append(res3)
   # if len(pro_targets) == 1:
   #  target_residues = pro_targets
   #for r in target_residues:
   # topology["DIHEDRALS_A"][r].append([at1, at2, at3, at4, ind])
   topology["DIHEDRALS_A"][max(res2,res3)].append([at1, at2, at3, at4, ind])
  else:
   #topology["IMPROPERS_A"][res1].append([at1, at2, at3, -at4, ind])
   #if res4 != res1:
   # topology["IMPROPERS_A"][res4].append([at1, at2, at3, -at4, ind])
   #if res2 != res1 and res2 != res4:
   # topology["IMPROPERS_A"][res2].append([at1, at2, at3, -at4, ind])
   #target_residues = [res1]
   #if res4 != res1:
   # target_residues.append(res4)
   #if res2 != res1 and res2 != res4:
   # target_residues.append(res2)
   #if (resn3 == "PRO" and atmn3 == "N") and (
   # (resn1 == "PRO" and atmn1 == "CD") or
   # (resn2 == "PRO" and atmn2 == "CD") or
   # (resn4 == "PRO" and atmn4 == "CD")
   #):
   # target_residues = [res3]
   #for r in target_residues:
   # topology["IMPROPERS_A"][r].append([at1, at2, at3, -at4, ind])
   topology["IMPROPERS_A"][res3].append([at1, at2, at3, -at4, ind])

def write_flag(file,title,format):
 flag_str = "%FLAG " + title
 format_str = "%FORMAT(" + format + ")"
 file.write(f"{flag_str:80s}" + "\n")
 file.write(f"{format_str:80s}" + "\n")

def write_int(file,values):
 for i, val in enumerate(values):
  file.write(f"{val:8d}")
  if (i + 1) % 10 == 0:
   file.write("\n")
 if len(values) % 10 != 0:
  file.write("\n")

def write_str(file,values):
 for i, val in enumerate(values):
  file.write(f"{val:4s}")
  if (i + 1) % 20 == 0:
   file.write("\n")
 if len(values) % 20 != 0:
  file.write("\n")

def write_dbl(file,values):
 for i, val in enumerate(values):
  file.write(f"{val:16.8E}") 
  if (i + 1) % 5 == 0:
   file.write("\n")
 if len(values) % 5 != 0:
  file.write("\n")

def write_AMBER_topology(file_name, topology):
 with open(file_name, 'w') as file:
  now = datetime.now()
  date_str = now.strftime("%m/%d/%y  %H:%M:%S")
  line = f'%VERSION  VERSION_STAMP = V0001.000  DATE = {date_str}'
  file.write(line + "\n")
  int_f = "10I8"
  str_f = "20a4"
  dbl_f = "5E16.8"
  write_flag(file,"TITLE",str_f)
  file.write("Exa4Mind_topology\n")
  write_flag(file,"POINTERS",int_f)
  write_int(file,topology["POINTERS"])
  write_flag(file,"ATOM_NAME",str_f)
  write_str(file,topology["ATOM_NAME"]) 
  write_flag(file,"CHARGE",dbl_f)
  write_dbl(file,topology["CHARGE"])
  write_flag(file,"ATOMIC_NUMBER",int_f)
  write_int(file,topology["ATOMIC_NUMBER"])
  write_flag(file,"MASS",dbl_f)
  write_dbl(file,topology["MASS"])
  write_flag(file,"ATOM_TYPE_INDEX",int_f)
  write_int(file,topology["ATOM_TYPE_INDEX"])
  write_flag(file,"NUMBER_EXCLUDED_ATOMS",int_f)
  write_int(file,topology["NUMBER_EXCLUDED_ATOMS"])
  write_flag(file,"NONBONDED_PARM_INDEX",int_f)
  write_int(file,topology["NONBONDED_PARM_INDEX"])
  write_flag(file,"RESIDUE_LABEL",str_f)
  write_str(file,topology["RESIDUE_LABEL"])
  write_flag(file,"RESIDUE_POINTER",int_f)
  write_int(file,topology["RESIDUE_POINTER"])
  write_flag(file,"BOND_FORCE_CONSTANT",dbl_f)
  write_dbl(file,topology["BOND_FORCE_CONSTANT"])
  write_flag(file,"BOND_EQUIL_VALUE",dbl_f)
  write_dbl(file,topology["BOND_EQUIL_VALUE"])
  write_flag(file,"ANGLE_FORCE_CONSTANT",dbl_f)
  write_dbl(file,topology["ANGLE_FORCE_CONSTANT"])
  write_flag(file,"ANGLE_EQUIL_VALUE",dbl_f)
  write_dbl(file,topology["ANGLE_EQUIL_VALUE"])
  write_flag(file,"DIHEDRAL_FORCE_CONSTANT",dbl_f)
  write_dbl(file,topology["DIHEDRAL_FORCE_CONSTANT"])
  write_flag(file,"DIHEDRAL_PERIODICITY",dbl_f)
  write_dbl(file,topology["DIHEDRAL_PERIODICITY"])
  write_flag(file,"DIHEDRAL_PHASE",dbl_f)
  write_dbl(file,topology["DIHEDRAL_PHASE"])
  write_flag(file,"SCEE_SCALE_FACTOR",dbl_f)
  write_dbl(file,topology["SCEE_SCALE_FACTOR"])
  write_flag(file,"SCNB_SCALE_FACTOR",dbl_f)
  write_dbl(file,topology["SCNB_SCALE_FACTOR"])
  write_flag(file,"SOLTY",dbl_f)
  write_dbl(file,topology["SOLTY"])
  write_flag(file,"LENNARD_JONES_ACOEF",dbl_f)
  write_dbl(file,topology["LENNARD_JONES_ACOEF"])
  write_flag(file,"LENNARD_JONES_BCOEF",dbl_f)
  write_dbl(file,topology["LENNARD_JONES_BCOEF"])
  write_flag(file,"BONDS_INC_HYDROGEN",int_f)
  write_int(file,topology["BONDS_INC_HYDROGEN"])
  write_flag(file,"BONDS_WITHOUT_HYDROGEN",int_f)
  write_int(file,topology["BONDS_WITHOUT_HYDROGEN"])
  write_flag(file,"ANGLES_INC_HYDROGEN",int_f)
  write_int(file,topology["ANGLES_INC_HYDROGEN"])
  write_flag(file,"ANGLES_WITHOUT_HYDROGEN",int_f)
  write_int(file,topology["ANGLES_WITHOUT_HYDROGEN"])
  write_flag(file,"DIHEDRALS_INC_HYDROGEN",int_f)
  write_int(file,topology["DIHEDRALS_INC_HYDROGEN"])
  write_flag(file,"DIHEDRALS_WITHOUT_HYDROGEN",int_f)
  write_int(file,topology["DIHEDRALS_WITHOUT_HYDROGEN"])
  write_flag(file,"EXCLUDED_ATOMS_LIST",int_f)
  write_int(file,topology["EXCLUDED_ATOMS_LIST"])
  write_flag(file,"HBOND_ACOEF",dbl_f)
  write_dbl(file,topology["HBOND_ACOEF"])
  write_flag(file,"HBOND_BCOEF",dbl_f)
  write_dbl(file,topology["HBOND_BCOEF"])
  write_flag(file,"HBCUT",dbl_f)
  write_dbl(file,topology["HBCUT"])
  write_flag(file,"AMBER_ATOM_TYPE",str_f)
  write_str(file,topology["AMBER_ATOM_TYPE"])
  write_flag(file,"TREE_CHAIN_CLASSIFICATION",str_f)
  write_str(file,topology["TREE_CHAIN_CLASSIFICATION"])
  write_flag(file,"JOIN_ARRAY",int_f)
  write_int(file,topology["JOIN_ARRAY"])
  write_flag(file,"IROTAT",int_f)
  write_int(file,topology["IROTAT"])
  if "SOLVENT_POINTERS" in topology:
   write_flag(file,"SOLVENT_POINTERS","3I8")
   write_int(file,topology["SOLVENT_POINTERS"])
  if "ATOMS_PER_MOLECULE" in topology:
   write_flag(file,"ATOMS_PER_MOLECULE",int_f)
   write_int(file,topology["ATOMS_PER_MOLECULE"])
  if "BOX_DIMENSIONS" in topology:
   write_flag(file,"BOX_DIMENSIONS",dbl_f)
   write_dbl(file,topology["BOX_DIMENSIONS"])
  if "RADIUS_SET" in topology:
   write_flag(file,"RADIUS_SET","1a80")
   write_str(file,topology["RADIUS_SET"])
  write_flag(file,"RADII",dbl_f)
  write_dbl(file,topology["RADII"])
  write_flag(file,"SCREEN",dbl_f)
  write_dbl(file,topology["SCREEN"])
  if "IPOL" in topology:
   write_flag(file,"IPOL","1I8")
   write_int(file,topology["IPOL"])

def strip(topology, ntwprt):
 if (ntwprt+1) in topology["RESIDUE_POINTER"]:
  nres = topology["RESIDUE_POINTER"].index(ntwprt+1)
  topology["POINTERS"][11] = nres
 elif ntwprt >= topology["POINTERS"][0] or ntwprt == 0:
  return
 else:
  raise exception("Error: stripping the topology within the residue is not possible.")
 topology["POINTERS"][0] = ntwprt
 topology["ATOM_NAME"] = topology["ATOM_NAME"][:ntwprt]
 topology["CHARGE"] = topology["CHARGE"][:ntwprt]
 topology["ATOMIC_NUMBER"] = topology["ATOMIC_NUMBER"][:ntwprt]
 topology["MASS"] = topology["MASS"][:ntwprt]
 topology["ATOM_TYPE_INDEX"] = topology["ATOM_TYPE_INDEX"][:ntwprt]
 ntypes = topology["POINTERS"][1]
 topology["POINTERS"][1] = max(topology["ATOM_TYPE_INDEX"])
 topology["NUMBER_EXCLUDED_ATOMS"] = topology["NUMBER_EXCLUDED_ATOMS"][:ntwprt]
 NNB = sum(topology["NUMBER_EXCLUDED_ATOMS"])
 topology["POINTERS"][10] = NNB
 if ntypes != topology["POINTERS"][1]:
  new_parm_ind = []
  m = topology["POINTERS"][1]
  for i in range(m):
   new_parm_ind.extend(topology["NONBONDED_PARM_INDEX"][(i*ntypes):(i*ntypes + m)])
  topology["NONBONDED_PARM_INDEX"] = new_parm_ind
  max_pair_typ = max(topology["NONBONDED_PARM_INDEX"])
  topology["LENNARD_JONES_ACOEF"] = topology["LENNARD_JONES_ACOEF"][:max_pair_typ]
  topology["LENNARD_JONES_BCOEF"] = topology["LENNARD_JONES_BCOEF"][:max_pair_typ]
 topology["RESIDUE_LABEL"] = topology["RESIDUE_LABEL"][:nres]
 topology["RESIDUE_POINTER"] = topology["RESIDUE_POINTER"][:nres]
 ind_lim = 3*ntwprt
 new_bnd_h = []
 for i in range(0,len(topology["BONDS_INC_HYDROGEN"]),3):
  if topology["BONDS_INC_HYDROGEN"][i] < ind_lim and topology["BONDS_INC_HYDROGEN"][i+1] < ind_lim:
   new_bnd_h.extend([topology["BONDS_INC_HYDROGEN"][i], topology["BONDS_INC_HYDROGEN"][i+1], topology["BONDS_INC_HYDROGEN"][i+2]])
 topology["BONDS_INC_HYDROGEN"] = new_bnd_h
 topology["POINTERS"][2] = len(topology["BONDS_INC_HYDROGEN"])//3
 new_bnd_a = []
 for i in range(0,len(topology["BONDS_WITHOUT_HYDROGEN"]),3):
  if topology["BONDS_WITHOUT_HYDROGEN"][i] < ind_lim and topology["BONDS_WITHOUT_HYDROGEN"][i+1] < ind_lim:
   new_bnd_a.extend([topology["BONDS_WITHOUT_HYDROGEN"][i], topology["BONDS_WITHOUT_HYDROGEN"][i+1], topology["BONDS_WITHOUT_HYDROGEN"][i+2]])
 topology["BONDS_WITHOUT_HYDROGEN"] = new_bnd_a
 topology["POINTERS"][3] = len(topology["BONDS_WITHOUT_HYDROGEN"])//3
 topology["POINTERS"][12] = topology["POINTERS"][3]
 new_ang_h = []
 for i in range(0,len(topology["ANGLES_INC_HYDROGEN"]),4):
  if topology["ANGLES_INC_HYDROGEN"][i] < ind_lim and topology["ANGLES_INC_HYDROGEN"][i+1] < ind_lim and topology["ANGLES_INC_HYDROGEN"][i+2] < ind_lim:
   new_ang_h.extend([topology["ANGLES_INC_HYDROGEN"][i], topology["ANGLES_INC_HYDROGEN"][i+1], topology["ANGLES_INC_HYDROGEN"][i+2], topology["ANGLES_INC_HYDROGEN"][i+3]])
 topology["ANGLES_INC_HYDROGEN"] = new_ang_h
 topology["POINTERS"][4] = len(topology["ANGLES_INC_HYDROGEN"])//4
 new_ang_a = []
 for i in range(0,len(topology["ANGLES_WITHOUT_HYDROGEN"]),4):
  if topology["ANGLES_WITHOUT_HYDROGEN"][i] < ind_lim and topology["ANGLES_WITHOUT_HYDROGEN"][i+1] < ind_lim and topology["ANGLES_WITHOUT_HYDROGEN"][i+2] < ind_lim:
   new_ang_a.extend([topology["ANGLES_WITHOUT_HYDROGEN"][i], topology["ANGLES_WITHOUT_HYDROGEN"][i+1], topology["ANGLES_WITHOUT_HYDROGEN"][i+2], topology["ANGLES_WITHOUT_HYDROGEN"][i+3]])
 topology["ANGLES_WITHOUT_HYDROGEN"] = new_ang_a
 topology["POINTERS"][5] = len(topology["ANGLES_WITHOUT_HYDROGEN"])//4
 topology["POINTERS"][13] = topology["POINTERS"][5]
 new_dih_h = []
 for i in range(0,len(topology["DIHEDRALS_INC_HYDROGEN"]),5):
  if topology["DIHEDRALS_INC_HYDROGEN"][i] < ind_lim and topology["DIHEDRALS_INC_HYDROGEN"][i+1] < ind_lim and abs(topology["DIHEDRALS_INC_HYDROGEN"][i+2]) < ind_lim and abs(topology["DIHEDRALS_INC_HYDROGEN"][i+3]) < ind_lim:
   new_dih_h.extend([topology["DIHEDRALS_INC_HYDROGEN"][i], topology["DIHEDRALS_INC_HYDROGEN"][i+1], topology["DIHEDRALS_INC_HYDROGEN"][i+2], topology["DIHEDRALS_INC_HYDROGEN"][i+3], topology["DIHEDRALS_INC_HYDROGEN"][i+4]])
 topology["DIHEDRALS_INC_HYDROGEN"] = new_dih_h
 topology["POINTERS"][6] = len(topology["DIHEDRALS_INC_HYDROGEN"])//5
 new_dih_a = []
 for i in range(0,len(topology["DIHEDRALS_WITHOUT_HYDROGEN"]),5):
  if topology["DIHEDRALS_WITHOUT_HYDROGEN"][i] < ind_lim and topology["DIHEDRALS_WITHOUT_HYDROGEN"][i+1] < ind_lim and abs(topology["DIHEDRALS_WITHOUT_HYDROGEN"][i+2]) < ind_lim and abs(topology["DIHEDRALS_WITHOUT_HYDROGEN"][i+3]) < ind_lim:
   new_dih_a.extend([topology["DIHEDRALS_WITHOUT_HYDROGEN"][i], topology["DIHEDRALS_WITHOUT_HYDROGEN"][i+1], topology["DIHEDRALS_WITHOUT_HYDROGEN"][i+2], topology["DIHEDRALS_WITHOUT_HYDROGEN"][i+3], topology["DIHEDRALS_WITHOUT_HYDROGEN"][i+4]])
 topology["DIHEDRALS_WITHOUT_HYDROGEN"] = new_dih_a
 topology["POINTERS"][7] = len(topology["DIHEDRALS_WITHOUT_HYDROGEN"])//5
 topology["POINTERS"][14] = topology["POINTERS"][7]
 topology["POINTERS"][15] = max(topology["BONDS_INC_HYDROGEN"][2::3] + topology["BONDS_WITHOUT_HYDROGEN"][2::3])
 topology["BOND_FORCE_CONSTANT"] = topology["BOND_FORCE_CONSTANT"][:topology["POINTERS"][15]]
 topology["BOND_EQUIL_VALUE"] = topology["BOND_EQUIL_VALUE"][:topology["POINTERS"][15]]
 topology["POINTERS"][16] = max(topology["ANGLES_INC_HYDROGEN"][3::4] + topology["ANGLES_WITHOUT_HYDROGEN"][3::4])
 topology["ANGLE_FORCE_CONSTANT"] = topology["ANGLE_FORCE_CONSTANT"][:topology["POINTERS"][16]]
 topology["ANGLE_EQUIL_VALUE"] = topology["ANGLE_EQUIL_VALUE"][:topology["POINTERS"][16]]
 topology["POINTERS"][17] = max(topology["DIHEDRALS_INC_HYDROGEN"][4::5] + topology["DIHEDRALS_WITHOUT_HYDROGEN"][4::5])
 topology["DIHEDRAL_FORCE_CONSTANT"] = topology["DIHEDRAL_FORCE_CONSTANT"][:topology["POINTERS"][17]]
 topology["DIHEDRAL_PERIODICITY"] = topology["DIHEDRAL_PERIODICITY"][:topology["POINTERS"][17]]
 topology["DIHEDRAL_PHASE"] = topology["DIHEDRAL_PHASE"][:topology["POINTERS"][17]]
 topology["SCEE_SCALE_FACTOR"] = topology["SCEE_SCALE_FACTOR"][:topology["POINTERS"][17]]
 topology["SCNB_SCALE_FACTOR"] = topology["SCNB_SCALE_FACTOR"][:topology["POINTERS"][17]]
 topology["EXCLUDED_ATOMS_LIST"] = topology["EXCLUDED_ATOMS_LIST"][:NNB]
 # SOLTY currently unused, no need to modify it
 NPHB = min(topology["NONBONDED_PARM_INDEX"])
 if NPHB < 0:
  topology["POINTERS"][19] = -NPHB
 else:
  topology["POINTERS"][19] = 0
 topology["HBOND_ACOEF"] = topology["HBOND_ACOEF"][:topology["POINTERS"][19]]
 topology["HBOND_BCOEF"] = topology["HBOND_BCOEF"][:topology["POINTERS"][19]]
 topology["HBCUT"] = topology["HBCUT"][:topology["POINTERS"][19]]
 topology["AMBER_ATOM_TYPE"] = topology["AMBER_ATOM_TYPE"][:ntwprt]
 topology["TREE_CHAIN_CLASSIFICATION"] = topology["TREE_CHAIN_CLASSIFICATION"][:ntwprt]
 topology["JOIN_ARRAY"] = topology["JOIN_ARRAY"][:ntwprt]
 topology["IROTAT"] = topology["IROTAT"][:ntwprt]
 topology["RADII"] = topology["RADII"][:ntwprt]
 topology["SCREEN"] = topology["SCREEN"][:ntwprt]
 if topology["POINTERS"][27] != 0:
  cumulative_sum = 0
  for i, value in enumerate(topology["ATOMS_PER_MOLECULE"]):
   cumulative_sum += value
   if cumulative_sum >= ntwprt:
    NSPM = i+1
    break
  topology["ATOMS_PER_MOLECULE"] = topology["ATOMS_PER_MOLECULE"][:NSPM]
  if topology["SOLVENT_POINTERS"][0] > nres:
   topology["SOLVENT_POINTERS"][0] = nres
  topology["SOLVENT_POINTERS"][1] = NSPM
  if topology["SOLVENT_POINTERS"][2] > NSPM:
   topology["SOLVENT_POINTERS"][2] = NSPM+1
 atoms_per_res = []
 for i in range(1,len(topology["RESIDUE_POINTER"])):
  atoms_per_res.append(topology["RESIDUE_POINTER"][i]-topology["RESIDUE_POINTER"][i-1])
 atoms_per_res.append(topology["POINTERS"][0]-topology["RESIDUE_POINTER"][-1])
 topology["POINTERS"][28] = max(atoms_per_res)

def read_mol_type(topology):
 molecule_type = []
 molecule_type_per_res = []
 molecule_types_by_atom_types = [[] for _ in range(topology["POINTERS"][1])]
 converted_residue_labels = []
 convert_dict = get_conversion_dict()
 for i in range(topology["POINTERS"][11]):
  beg = topology["RESIDUE_POINTER"][i] - 1
  end = topology["RESIDUE_POINTER"][i+1] - 1 if i < (topology["POINTERS"][11] - 1) else topology["POINTERS"][0]
  resn = alias.resn_alias(topology["RESIDUE_LABEL"][i])
  atom_names = [alias.name_alias(None,name) for name in topology["ATOM_NAME"][beg:end]]
  target_keys = set(atom_names)
  for mol_type, residues in convert_dict.items():
   if resn not in residues:
    continue
   candidate_resnames = [resn]
   term_list = residues[resn].get("terminal",[])
   for alt in term_list:
    if alt in residues:
     candidate_resnames.append(alt)
   for cand in candidate_resnames:
    atoms = residues[cand]["atom"]
    if set(atoms.keys()) == target_keys:
     for atom_type_index in topology["ATOM_TYPE_INDEX"][beg:end]:
      if mol_type not in molecule_types_by_atom_types[atom_type_index-1]:
       molecule_types_by_atom_types[atom_type_index-1].append(mol_type)
     molecule_type.extend([mol_type] * len(atom_names))
     molecule_type_per_res.append(mol_type)
     converted_residue_labels.append(cand)
     break
   else:
    continue
   break
  else:
   molecule_type.extend(["unk"] * len(atom_names))
   molecule_type_per_res.append("unk")
   converted_residue_labels.append(resn)
 for k, lst in enumerate(molecule_types_by_atom_types):
  if not lst:
   molecule_types_by_atom_types[k] = ["unk"]
 topology["MOLECULE_TYPE"] = molecule_type
 topology["MOLECULE_TYPE_PER_RES"] = molecule_type_per_res
 topology["MOLECULE_TYPES_BY_ATOM_TYPES"] = molecule_types_by_atom_types
 topology["CONVERTED_RESIDUE_LABEL"] = converted_residue_labels

def correct_VdWmatrix_by_molecule_types(topology):
 # check topology nonbonded information
 required = (
  "NONBONDED_DIAGONAL_R",
  "NONBONDED_DIAGONAL_EPS",
  "NONBONDED_NBFIX",
  "MOLECULE_TYPE",
  "MOLECULE_TYPES_BY_ATOM_TYPES",
 )
 missing = [k for k in required if k not in topology]
 if missing:
  raise Exception(f"Missing keys in topology: {missing} – run nonbonded preprocessing first.")

 has_inter_nbfix = False
 for i, j in topology["NONBONDED_NBFIX"]:
  mols_i = set(topology["MOLECULE_TYPES_BY_ATOM_TYPES"][i])
  mols_j = set(topology["MOLECULE_TYPES_BY_ATOM_TYPES"][j])
  if mols_i != mols_j:
   has_inter_nbfix = True
   break
 if not has_inter_nbfix:
  return True

 VdW_tol = 1e-3
 zero_tol = 1e-7
 
 # build NBfix 0/1 matrix and VdW matrix
 ntyp = topology["POINTERS"][1]
 M = [[0]*ntyp for _ in range(ntyp)]
 Am = [[0.0]*ntyp for _ in range(ntyp)]
 Bm = [[0.0]*ntyp for _ in range(ntyp)]
 VdW_R = [[0.0]*ntyp for _ in range(ntyp)]
 VdW_eps = [[0.0]*ntyp for _ in range(ntyp)]
 for i, j in topology["NONBONDED_NBFIX"]:
  M[i][j] = 1
  M[j][i] = 1
 for i in range(ntyp):
  for j in range(ntyp):
   index = topology["NONBONDED_PARM_INDEX"][ntyp*i+j]
   if index > 0:
    A = topology["LENNARD_JONES_ACOEF"][index-1]
    B = topology["LENNARD_JONES_BCOEF"][index-1]
   elif index < 0:
    A = topology["HBOND_ACOEF"][-index-1]
    B = topology["HBOND_BCOEF"][-index-1]
   else:
    raise Exception("NONBONDED_PARM_INDEX matrix contains index equal to zero")
   Am[i][j] = A
   Bm[i][j] = B
   if abs(A) < zero_tol or abs(B) < zero_tol:
    M[i][j] = 2
    VdW_R[i][j] = 0.0
    VdW_eps[i][j] = 0.0
   else:
    VdW_R[i][j] = (2*A/B)**(1.0/6.0)
    VdW_eps[i][j] = 0.25*B*B/A

 # find candidates for VdW matrix correction
 candidates = []
 for i in range(ntyp):
  mols_i = topology["MOLECULE_TYPES_BY_ATOM_TYPES"][i]
  outside_js = [j for j in range(ntyp) if j != i and any(mol_j not in mols_i for mol_j in topology["MOLECULE_TYPES_BY_ATOM_TYPES"][j])]
  if outside_js and all(M[i][j] == 1 or M[i][j] == 2 for j in outside_js) and not all(M[i][j] == 2 for j in range(ntyp)):
   candidates.append(i)

 # find VdW matrix correction
 nonbonded_diagonal_R = topology["NONBONDED_DIAGONAL_R"].copy()
 nonbonded_diagonal_eps = topology["NONBONDED_DIAGONAL_EPS"].copy()
 for i in candidates:
  mols_i = topology["MOLECULE_TYPES_BY_ATOM_TYPES"][i]
  outside_js = [j for j in range(ntyp) if j != i and any(mol_j not in mols_i for mol_j in topology["MOLECULE_TYPES_BY_ATOM_TYPES"][j]) and j not in candidates and M[i][j] != 2]
  R_can = []
  eps_can = []
  for j in outside_js:
    R_can.append(VdW_R[i][j]-topology["NONBONDED_DIAGONAL_R"][j])
    eps_can.append(VdW_eps[i][j]*VdW_eps[i][j]/topology["NONBONDED_DIAGONAL_EPS"][j])
  R_tol = max(R_can, default=topology["NONBONDED_DIAGONAL_R"][i]) - min(R_can, default=topology["NONBONDED_DIAGONAL_R"][i])
  R_est = (max(R_can, default=topology["NONBONDED_DIAGONAL_R"][i]) + min(R_can, default=topology["NONBONDED_DIAGONAL_R"][i]))/2.0
  eps_tol = max(eps_can, default=topology["NONBONDED_DIAGONAL_EPS"][i]) - min(eps_can, default=topology["NONBONDED_DIAGONAL_EPS"][i])
  eps_est = (max(eps_can, default=topology["NONBONDED_DIAGONAL_EPS"][i]) + min(eps_can, default=topology["NONBONDED_DIAGONAL_EPS"][i]))/2.0
  #print(R_can)
  #print(eps_can)
  #print(R_est,R_tol,eps_est,eps_tol)
  if abs(R_tol) > VdW_tol or abs(eps_tol) > VdW_tol:
   R_est = topology["NONBONDED_DIAGONAL_R"][i]
   eps_est = topology["NONBONDED_DIAGONAL_EPS"][i]
  nonbonded_diagonal_R[i] = R_est
  nonbonded_diagonal_eps[i] = eps_est

 # Update NBfix matrix
 N = [[0]*ntyp for _ in range(ntyp)]
 new_nbfix_pairs = []
 for i in range(ntyp):
  for j in range(ntyp):
   Rij = nonbonded_diagonal_R[i] + nonbonded_diagonal_R[j]
   epsij = math.sqrt(nonbonded_diagonal_eps[i]*nonbonded_diagonal_eps[j])
   if abs(VdW_R[i][j]) < zero_tol or abs(VdW_eps[i][j]) < zero_tol:
    N[i][j] = 2
   elif abs(Rij - VdW_R[i][j]) > VdW_tol or abs(epsij - VdW_eps[i][j]) > VdW_tol:
    mols_i = set(topology["MOLECULE_TYPES_BY_ATOM_TYPES"][i])
    mols_j = set(topology["MOLECULE_TYPES_BY_ATOM_TYPES"][j])
    if mols_i != mols_j:
     return False
    N[i][j] = 1
    if i <= j:
     new_nbfix_pairs.append((i, j))
   else:
    N[i][j] = 0

 topology["NONBONDED_DIAGONAL_R"]   = nonbonded_diagonal_R
 topology["NONBONDED_DIAGONAL_EPS"] = nonbonded_diagonal_eps
 topology["NONBONDED_NBFIX"]        = new_nbfix_pairs

 return True

def create_AMBER_topology(mol):

    # initial setting
    top = {}
    bondii_screen_by_Z = {
        1: 0.85,   # H
        6: 0.72,   # C
        7: 0.79,   # N
        8: 0.85,   # O
        9: 0.88,   # F
        15: 0.86,  # P
        16: 0.96,  # S
    }
    bondii_radii_by_Z = {
        1: 1.2,    # H (will be modified for mbondii)
        6: 1.7,    # C (ALL-ATOM, ignore UA heuristic)
        7: 1.55,   # N
        8: 1.5,    # O
        9: 1.5,    # F
        14: 2.1,   # Si
        15: 1.85,  # P
        16: 1.8,   # S
        17: 1.7,   # Cl
    }
    tree_map = {
        0: 'M',
        1: 'E',
        2: 'S',
        3: 'B',
        4: '3',
        5: '4',
        6: '5',
        7: '6'
    }


    # initialize topology lists
    top['POINTERS'] = [0]*31
    top['ATOM_NAME'] = []
    top['CHARGE'] = []
    top['ATOMIC_NUMBER'] = []
    top['MASS'] = []
    top['ATOM_TYPE_INDEX'] = []
    top['NUMBER_EXCLUDED_ATOMS'] = []
    top['NONBONDED_PARM_INDEX'] = []
    top['RESIDUE_LABEL'] = [] 
    top['RESIDUE_POINTER'] = [] 
    top['BOND_FORCE_CONSTANT'] = []
    top['BOND_EQUIL_VALUE'] = []
    top['ANGLE_FORCE_CONSTANT'] = []
    top['ANGLE_EQUIL_VALUE'] = []
    top['DIHEDRAL_FORCE_CONSTANT'] = []
    top['DIHEDRAL_PERIODICITY'] = []
    top['DIHEDRAL_PHASE'] = []
    top['SCEE_SCALE_FACTOR'] = []
    top['SCNB_SCALE_FACTOR'] = []
    top['SOLTY'] = []
    top['LENNARD_JONES_ACOEF'] = []
    top['LENNARD_JONES_BCOEF'] = []
    top['BONDS_INC_HYDROGEN'] = []
    top['BONDS_WITHOUT_HYDROGEN'] = []
    top['ANGLES_INC_HYDROGEN'] = []
    top['ANGLES_WITHOUT_HYDROGEN'] = []
    top['DIHEDRALS_INC_HYDROGEN'] = []
    top['DIHEDRALS_WITHOUT_HYDROGEN'] = []
    top['EXCLUDED_ATOMS_LIST'] = []
    top['HBOND_ACOEF'] = []
    top['HBOND_BCOEF'] = []
    top['HBCUT'] = []
    top['AMBER_ATOM_TYPE'] = []
    top['TREE_CHAIN_CLASSIFICATION'] = []
    top['JOIN_ARRAY'] = []
    top['IROTAT'] = []
    top['RADIUS_SET'] = (['modi','fied',' Bon','di r','adii',' (mb','ondi',')   '] + ['    ']*12)
    top['RADII'] = []
    top['SCREEN'] = []
    top['IPOL'] = [0]

    # auxiliary settings
    new_types = []
    connectivity = {}
    mol_types = []
    VdW = []
    bond_parms = []
    angle_parms = []
    dihedral_parms = []
    previous_unit = None
    previous_subsum = 0
    old_types = []
    natom = 0
    nmxrs = 0
    numextra = 0
    atoms_in_res = []
    atoms_in_mol = []
    if 'box' in mol:
        ifbox = 1
    else:
        ifbox = 0
    nspsol = 0
    iptres = 0
    solute = True

    # atomic parameters
    for res in mol['residues']:
        res_name = res['resn']
        mol_type = res['mol_type']
        force_field = mol['force_field_data'][mol_type]
        unit = force_field.units[res_name]
        top['RESIDUE_LABEL'].append(res_name)
        top['RESIDUE_POINTER'].append(natom+1)
        res_size = len(unit['atoms']['name'])
        atoms_in_res.append(res_size)
        if nmxrs < res_size:
            nmxrs = res_size
        for at, name in enumerate(unit['atoms']['name']):
            type = unit['atoms']['type'][at]
            charge = unit['atoms']['charge'][at]
            at_num = unit['atoms']['at_num'][at]
            if at_num == 0:
                numextra += 1
            mass = unit['atoms']['mass'][at]
            R = unit['atoms']['R'][at]
            eps = unit['atoms']['eps'][at]
            if force_field.types[type].strip().upper() == 'HW':
                hw_type = 1
            else:
                hw_type = 0
            if (R, eps, hw_type) not in VdW:
                VdW.append((R,eps,hw_type))
            mol_types.append(mol_type)
            top['ATOM_NAME'].append(name)
            top['CHARGE'].append(charge*18.2223)
            top['ATOMIC_NUMBER'].append(at_num)
            top['MASS'].append(mass)
            top['ATOM_TYPE_INDEX'].append(VdW.index((R,eps,hw_type))+1)
            new_types.append(type)
            top['AMBER_ATOM_TYPE'].append(force_field.types[type])
            if force_field.types[type] not in old_types:
                old_types.append(force_field.types[type])
        connected_to_previous = False
        for bond in reversed(unit['bonds']):
            if bond[0][0] == '-':
                if previous_unit is None:
                    raise Exception("First residue is not terminal residue and has connectivity to previous.")
                a = previous_subsum + previous_unit['atoms']['name'].index(bond[0][1:])
                b = natom + unit['atoms']['name'].index(bond[1])
                connected_to_previous = True
                #print(f"{bond[0][1:]}-{bond[1]} : {a+1} {b+1} : {new_types[a]}-{new_types[b]}")
            else:
                a = natom + unit['atoms']['name'].index(bond[0])
                b = natom + unit['atoms']['name'].index(bond[1])
                #print(f"{bond[0]}-{bond[1]} : {a+1} {b+1} : {new_types[a]}-{new_types[b]}")
            if a not in connectivity:
                connectivity[a] = []
            if b not in connectivity:
                connectivity[b] = []
            connectivity[a].append(b)
            connectivity[b].append(a)
            key = min((new_types[a],new_types[b]),(new_types[b],new_types[a]))
            #print(force_field.b['bondtypes'][key])
            bond_p = force_field.b['bondtypes'][key]
            if bond_p not in bond_parms:
                #print("appending bond parm list")
                bond_parms.append(bond_p)
            if top['ATOMIC_NUMBER'][a] == 1 or top['ATOMIC_NUMBER'][b] == 1: # bonds including hydrogen
                top['BONDS_INC_HYDROGEN'].extend([3*a, 3*b, bond_parms.index(bond_p)])
            else:
                top['BONDS_WITHOUT_HYDROGEN'].extend([3*a, 3*b, bond_parms.index(bond_p)])
        if connected_to_previous:
            atoms_in_mol[-1] += res_size
        else:
            atoms_in_mol.append(res_size)
        if solute:
            if mol_type.startswith('W'):
                nspsol = len(atoms_in_mol)
                solute = False
            else:
                iptres += 1
        previous_unit = unit
        previous_subsum = natom
        natom += len(unit['atoms']['name'])

    # bonded parameters
    for b, neighbors in connectivity.items():  
        force_field = mol['force_field_data'][mol_types[b]]
        for a, c in combinations(neighbors,2):
            #print(f"{a+1}-{b+1}-{c+1}")
            key = min((new_types[a],new_types[b],new_types[c]),(new_types[c],new_types[b],new_types[a]))
            #print(force_field.b['angletypes'][key])
            angle_p = force_field.b['angletypes'][key]
            if angle_p not in angle_parms:
                angle_parms.append(angle_p)
            if any(top['ATOMIC_NUMBER'][i] == 1 for i in [a,b,c]):
                top['ANGLES_INC_HYDROGEN'].extend([3*a, 3*b, 3*c, angle_parms.index(angle_p)])
            else:
                top['ANGLES_WITHOUT_HYDROGEN'].extend([3*a, 3*b, 3*c, angle_parms.index(angle_p)])
        for c in (nt for nt in neighbors if nt > b):
            for a in (nt for nt in neighbors if nt != c):
                for d in (nt for nt in connectivity[c] if nt != b and nt != a):
                    #print(f"{a+1}-{b+1}-{c+1}-{d+1}")
                    key = min((new_types[a],new_types[b],new_types[c],new_types[d]),(new_types[d],new_types[c],new_types[b],new_types[a]))
                    #print(force_field.b['dihedraltypes'][key])
                    dihedral_p = force_field.b['dihedraltypes'][key]
                    pn_fact = 1
                    for dihe_term in dihedral_p:
                        if abs(dihe_term[0]) < 1e-7:
                            continue
                        if dihe_term not in dihedral_parms:
                            dihedral_parms.append(dihe_term)
                        if any(top['ATOMIC_NUMBER'][i] == 1 for i in [a,b,c,d]):
                            top['DIHEDRALS_INC_HYDROGEN'].extend([3*a, 3*b, 3*c*pn_fact, 3*d, dihedral_parms.index(dihe_term)])
                        else:
                            top['DIHEDRALS_WITHOUT_HYDROGEN'].extend([3*a, 3*b, 3*c*pn_fact, 3*d, dihedral_parms.index(dihe_term)])
                        pn_fact = -1
    for ic0, neighbors in connectivity.items():
        force_field = mol['force_field_data'][mol_types[ic0]]
        if len(neighbors) == 3: # impropers
            ia0, ib0, id0 = neighbors[0], neighbors[1], neighbors[2]
            key_candidates = (
                ((ia0,ib0,ic0,id0),(new_types[ia0],new_types[ib0],new_types[ic0],new_types[id0])),
                ((ia0,id0,ic0,ib0),(new_types[ia0],new_types[id0],new_types[ic0],new_types[ib0])),
                ((ib0,ia0,ic0,id0),(new_types[ib0],new_types[ia0],new_types[ic0],new_types[id0])),
                ((ib0,id0,ic0,ia0),(new_types[ib0],new_types[id0],new_types[ic0],new_types[ia0])),
                ((id0,ia0,ic0,ib0),(new_types[id0],new_types[ia0],new_types[ic0],new_types[ib0])),
                ((id0,ib0,ic0,ia0),(new_types[id0],new_types[ib0],new_types[ic0],new_types[ia0]))
            )
            for (ia,ib,ic,id), key in key_candidates:
                if key in force_field.b['impropertypes']:
                    dihedral_p = force_field.b['impropertypes'][key]
                    break
            else:
                raise Exception(
                    f"Improper dihedral type not found for any permutation: "
                    f"{[k for _, k in key_candidates]} in {force_field.b['impropertypes']}"
                )
            #print(f"{ia+1}-{ib+1}-{ic+1}-{id+1} {key}")
            for dihe_term in dihedral_p:
                if abs(dihe_term[0]) < 1e-7:
                    continue
                if dihe_term not in dihedral_parms:
                    dihedral_parms.append(dihe_term)
                if any(top['ATOMIC_NUMBER'][i] == 1 for i in [ia,ib,ic,id]):
                    top['DIHEDRALS_INC_HYDROGEN'].extend([3*ia, 3*ib, -3*ic, -3*id, dihedral_parms.index(dihe_term)])
                else:
                    top['DIHEDRALS_WITHOUT_HYDROGEN'].extend([3*ia, 3*ib, -3*ic, -3*id, dihedral_parms.index(dihe_term)])

    # excluded list and GB parameters
    for at in range(natom):
        # excluded list code
        if at not in connectivity:
            top['NUMBER_EXCLUDED_ATOMS'].append(1)
            top['EXCLUDED_ATOMS_LIST'].append(0)
            continue
        visited = {at}
        q = deque([(at,0)])
        out = set()
        while q:
            node, dist = q.popleft()
            if dist == 3:
                continue
            for nb in connectivity.get(node, []):
                if nb in visited:
                    continue
                visited.add(nb)
                nd = dist + 1
                out.add(nb)
                q.append((nb,nd))
        out_sorted = sorted(nb + 1 for nb in out if nb > at)
        top['NUMBER_EXCLUDED_ATOMS'].append(len(out_sorted))
        top['EXCLUDED_ATOMS_LIST'].extend(out_sorted)
        # GB radii and screen code
        Z = top['ATOMIC_NUMBER'][at]
        screen = bondii_screen_by_Z.get(Z, 0.80)
        radius = bondii_radii_by_Z.get(Z, 1.5)
        if Z == 1:
            if top['AMBER_ATOM_TYPE'][at].strip().upper() == 'HW':
                radius = 0.8
            else:
                neighbors = connectivity.get(at,[])
                if neighbors:
                    nb0 = neighbors[0]
                    Znb = top['ATOMIC_NUMBER'][nb0]
                    if Znb == 6:  # H bonded on C
                        radius = 1.3
                    elif Znb in (8, 16): # H on O or S
                        radius = 0.8
                    elif Znb == 7: # H on N
                        radius = 1.3
        top['RADII'].append(radius)
        top['SCREEN'].append(screen)

    # parameter settings
    for parm in bond_parms:
        top['BOND_FORCE_CONSTANT'].append(parm[0]/4.184/100/2)
        top['BOND_EQUIL_VALUE'].append(parm[1]*10)
    for parm in angle_parms:
        top['ANGLE_FORCE_CONSTANT'].append(parm[0]/4.184/2)
        top['ANGLE_EQUIL_VALUE'].append(math.radians(parm[1]))
    for parm in dihedral_parms:
        top['DIHEDRAL_FORCE_CONSTANT'].append(parm[0]/4.184)
        top['DIHEDRAL_PERIODICITY'].append(parm[2])
        top['DIHEDRAL_PHASE'].append(math.radians(parm[1]))
        top['SCEE_SCALE_FACTOR'].append(1.2)
        top['SCNB_SCALE_FACTOR'].append(2.0)
    top['SOLTY'] = [0.0]*len(old_types)
    top['JOIN_ARRAY'] = [0]*natom
    top['IROTAT'] = [0]*natom

    # Van der Waals
    VdW_types = len(VdW)
    hb_ind = 0
    for i in range(VdW_types):
        VdW[i] = (VdW[i][0]*(2**(1/6))*5,VdW[i][1]/4.184,VdW[i][2])
    for i in range(VdW_types):
        #print(f"R={VdW[i][0]*(2**(1/6))*5}, eps={VdW[i][1]/4.184}")
        #print(f"R={VdW[i][0]}, eps={VdW[i][1]}")
        for j in range(VdW_types):
            ind_square_mat = VdW_types*i + j
            m = max(i,j)
            n = min(i,j)
            ind_triang_mat = m * (m + 1) // 2 + n + 1
            eps = math.sqrt(VdW[i][1]*VdW[j][1])
            R   = VdW[i][0]+VdW[j][0]
            R2  = R*R
            R6  = R2*R2*R2
            R12 = R6*R6
            A = eps*R12
            B = 2*eps*R6
            if VdW[i][2] == 1 and VdW[j][2] == 1: # water HW-HW interaction
                hb_ind += 1
                top['NONBONDED_PARM_INDEX'].append(-hb_ind)
                if i >= j:
                    top['LENNARD_JONES_ACOEF'].append(0.0)
                    top['LENNARD_JONES_BCOEF'].append(0.0)
                    top['HBOND_ACOEF'].append(A)
                    top['HBOND_BCOEF'].append(B)
                    top['HBCUT'].append(0.0)
            else:
                top['NONBONDED_PARM_INDEX'].append(ind_triang_mat)
                if i >= j:
                    top['LENNARD_JONES_ACOEF'].append(A)
                    top['LENNARD_JONES_BCOEF'].append(B)

    # tree chain classification
    for at in range(natom):
        neighbors = connectivity.get(at, [])
        nb_len = len(neighbors)
        if mol_types[at].startswith('W'):
            val = 'BLA'
        elif len(new_types[at]) > 1 and new_types[at][1] == 'B' and nb_len != 1:
            val = 'M'
        else:
            val = tree_map.get(nb_len, 'X')
        top['TREE_CHAIN_CLASSIFICATION'].append(val)

    # box info
    if ifbox == 1:
        top['SOLVENT_POINTERS'] = [iptres, len(atoms_in_mol), nspsol]
        top['ATOMS_PER_MOLECULE'] = atoms_in_mol
        top['BOX_DIMENSIONS'] = [90.0, mol['box'][0], mol['box'][1], mol['box'][2]]

    # pointers
    top['POINTERS'][0] = natom
    top['POINTERS'][1] = VdW_types
    top['POINTERS'][2] = len(top['BONDS_INC_HYDROGEN'])
    top['POINTERS'][3] = len(top['BONDS_WITHOUT_HYDROGEN'])
    top['POINTERS'][4] = len(top['ANGLES_INC_HYDROGEN'])
    top['POINTERS'][5] = len(top['ANGLES_WITHOUT_HYDROGEN'])
    top['POINTERS'][6] = len(top['DIHEDRALS_INC_HYDROGEN'])
    top['POINTERS'][7] = len(top['DIHEDRALS_WITHOUT_HYDROGEN'])
    top['POINTERS'][8] = 0
    top['POINTERS'][9] = 0
    top['POINTERS'][10] = len(top['EXCLUDED_ATOMS_LIST'])
    top['POINTERS'][11] = len(top['RESIDUE_LABEL'])
    top['POINTERS'][12] = top['POINTERS'][3]
    top['POINTERS'][13] = top['POINTERS'][5]
    top['POINTERS'][14] = top['POINTERS'][7]
    top['POINTERS'][15] = len(top['BOND_FORCE_CONSTANT'])
    top['POINTERS'][16] = len(top['ANGLE_FORCE_CONSTANT'])
    top['POINTERS'][17] = len(top['DIHEDRAL_FORCE_CONSTANT'])
    top['POINTERS'][18] = len(old_types)
    top['POINTERS'][19] = len(top['HBOND_ACOEF'])
    top['POINTERS'][20] = 0
    top['POINTERS'][21] = 0
    top['POINTERS'][22] = 0
    top['POINTERS'][23] = 0
    top['POINTERS'][24] = 0
    top['POINTERS'][25] = 0
    top['POINTERS'][26] = 0
    top['POINTERS'][27] = ifbox
    top['POINTERS'][28] = nmxrs
    top['POINTERS'][29] = 0
    top['POINTERS'][30] = numextra

    return top

