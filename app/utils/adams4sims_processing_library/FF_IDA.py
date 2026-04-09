from adams4sims_processing_library.utils.loader_helper import get_lines


# CLASS FOR HANDLING IDA FF TYPE

class ff:

 def __init__(self, rtp_file, nb_itp_file, b_itp_file, atp_file = None):
  self.units = {}
  self.nb = {}
  self.b  = {}
  self.types = {}
  self.read_rtp(rtp_file)
  self.read_nonbonded(nb_itp_file)
  self.read_bonded(b_itp_file)
  if atp_file is not None:
    self.read_atp(atp_file)

 def find_atom_type(self,name,resn):
  try:
   i = self.units[resn]["atoms"]["name"].index(name)
   return self.units[resn]["atoms"]["type"][i]
  except Exception as e:
   print(f"Cannot find type for atom {name} in residue {resn}, exception code {e}")
   return None

 def read_atp(self,atp_file):
  for line in get_lines(atp_file):
    if not line or line.startswith(";") or line.startswith("#") or len(line.split()) == 0:
     continue
    field = line.split()
    if len(field) < 4:
     raise Exception("Old atom types missing in UM3F atp file.")
    self.types[field[0]] = field[3]

 def read_rtp(self,rtp_file):
  section = None
  for line in get_lines(rtp_file):
    if not line or line.startswith(";") or line.startswith("#") or len(line.split()) == 0:
     continue
    if line.startswith("[") and line.split()[2] == "]":
     residue = line.split()[1]
     self.units[residue] = {"atoms":{"name":[], "type":[], "charge":[]}, "bonds":[], "impropers":[]}
    elif line.split()[0] == "[":
     section = line.split()[1]
    elif section == "atoms":
     field = line.split()
     if len(field) < 4:
      raise Exception("Incorrect number of filed in atoms section of rtp file.")
     self.units[residue]["atoms"]["name"].append(field[0])
     self.units[residue]["atoms"]["type"].append(field[1])
     self.units[residue]["atoms"]["charge"].append(float(field[2]))
    elif section == "bonds":
     field = line.split()
     if len(field) < 2:
      raise Exception("Incorrect number of filed in bonds section of rtp file.")
     self.units[residue]["bonds"].append(field)
    elif section == "impropers":
     field = line.split()
     if len(field) < 4:
      raise Exception("Incorrect number of filed in impropers section of rtp file.")
     self.units[residue]["impropers"].append(field)

 def read_nonbonded(self,nb_itp_file):
  # READ NONBONDED.ITP
  for line in get_lines(nb_itp_file):
    if not line or line.startswith(";") or line.startswith("#") or len(line.split()) == 0:
     continue
    if line.startswith("["):
     if line.split()[2] != "]" or line.split()[1] != "atomtypes":
      raise Exception("Incorrect section in nonbonded file.")
    else:
     field = line.split()
     if len(field) < 7:
      print(field)
      raise Exception("Incorrect number of filed in atomtypes section of nonbonded itp file.")
     self.nb[field[0]] = [int(field[1]),float(field[2]),float(field[5]),float(field[6])]
  # PASS RELEVANT DATA TO SELF.UNITS
  for res in self.units:
   self.units[res]["atoms"]["mass"] = []
   self.units[res]["atoms"]["R"] = []
   self.units[res]["atoms"]["eps"] = []
   self.units[res]["atoms"]["at_num"] = []
   for at in self.units[res]["atoms"]["type"]:
    self.units[res]["atoms"]["at_num"].append(self.nb[at][0])
    self.units[res]["atoms"]["mass"].append(self.nb[at][1])
    self.units[res]["atoms"]["R"].append(self.nb[at][2])
    self.units[res]["atoms"]["eps"].append(self.nb[at][3])

 def read_bonded(self,b_itp_file):
  # READ BONDED.ITP
  section = None
  self.b = {"bondtypes":{},"angletypes":{},"dihedraltypes":{},"impropertypes":{}}
  for line in get_lines(b_itp_file):
    if not line or line.startswith(";") or line.startswith("#") or len(line.split()) == 0:
     continue
    if line.startswith("["):
     if line.split()[2] != "]" or (line.split()[1] != "bondtypes" and line.split()[1] != "angletypes" and line.split()[1] != "dihedraltypes"):
      raise Exception("Incorrect section in nonbonded file.")
     section = line.split()[1]
    elif section == "bondtypes":
     field = line.split()
     if len(field) < 5:
      print(field)
      raise Exception("Incorrect number of filed in bondtypes section of bonded itp file.")
     at1 = field[0]
     at2 = field[1]
     key = min((at1,at2),(at2,at1))
     self.b["bondtypes"][key] = [float(field[4]),float(field[3])]
    elif section == "angletypes":
     field = line.split()
     if len(field) < 6:
      print(field)
      raise Exception("Incorrect number of filed in angletypes section of bonded itp file.")
     at1 = field[0]
     at2 = field[1]
     at3 = field[2]
     key = min((at1,at2,at3),(at3,at2,at1))
     self.b["angletypes"][key] = [float(field[5]),float(field[4])]
    elif section == "dihedraltypes":
     field = line.split()
     if len(field) < 8:
      print(field)
      raise Exception("Incorrect number of filed in dihedraltypes section of bonded itp file.")
     at1 = field[0] 
     at2 = field[1]
     at3 = field[2]
     at4 = field[3]
     # proper dihedral
     if field[4] == "9":
      key = min((at1,at2,at3,at4),(at4,at3,at2,at1))
      if key not in self.b["dihedraltypes"]:
       self.b["dihedraltypes"][key] = [[float(field[6]),float(field[5]),float(field[7])]]
      else:
       self.b["dihedraltypes"][key].append([float(field[6]),float(field[5]),float(field[7])])
     # improper dihedral
     elif field[4] == "4":
      key = min((at1,at2,at3,at4),(at2,at1,at3,at4),(at2,at4,at3,at1),(at4,at2,at3,at1),(at1,at4,at3,at2),(at4,at1,at3,at2))
      if key not in self.b["impropertypes"]:
       self.b["impropertypes"][key] = [[float(field[6]),float(field[5]),float(field[7])]]
      else:
       self.b["impropertypes"][key].append([float(field[6]),float(field[5]),float(field[7])])
     else:
      raise Exception(f"Unsupported dihedral type in bonded itp file: {field[4]}.")
    else:
     raise Exception("Some data in bonded.itp out of any section.")


