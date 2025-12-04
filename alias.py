
name_alias_dict ={
      'O5*': "O5'",
      'C5*': "C5'",
      'C4*': "C4'",
      'O4*': "O4'",
      "O1'": "O4'",
      'C3*': "C3'",
      'O3*': "O3'",
      'C2*': "C2'",
      'O2*': "O2'",
      'C1*': "C1'",
      'C5M': 'C7',
      'H1*': "H1'",
      'H2*1': "H2'",
      "H2'1": "H2'",
      'H2*2': "H2''",
      "H2'2": "H2''",
      "'H2'": "H2''",
      'H3*': "H3'",
      'H4*': "H4'",
      'H5*1': "H5'",
      "H5'1": "H5'",
      'H5*2': "H5''",
      "'H5'": "H5''",
      "H5'2": "H5''",
      "HO'2": "HO2'",
      'H5T': "HO5'",
      'H3T': "HO3'",
      'OA': 'OP1',
      'O1P': 'OP1',
      'OB': 'OP2',
      'O2P': 'OP2',
      'O3P': 'OP3'
}

#name_alias_dict ={
#     'R': {
#      'O5*': "O5'", 
#      'C5*': "C5'", 
#      'C4*': "C4'", 
#      'O4*': "O4'", 
#      "O1'": "O4'", 
#      'C3*': "C3'", 
#      'O3*': "O3'", 
#      'C2*': "C2'", 
#      'O2*': "O2'", 
#      'C1*': "C1'", 
#      'C7': 'C5M', 
#      'C5M': 'C7', 
#      'H1*': "H1'", 
#      'H2*1': "H2'", 
#      "H2'1": "H2'", 
#      'H2*2': "H2''", 
#      "H2'2": "H2''", 
#      'H3*': "H3'", 
#      'H4*': "H4'", 
#      'H5*1': "H5'", 
#      "H5'1": "H5'", 
#      'H5*2': "H5''", 
#      "H5'2": "H5''", 
#      "HO'2": "HO2'", 
#      'H5T': "HO5'", 
#      'H3T': "HO3'", 
#      'OA': 'OP1', 
#      'O1P': 'OP1', 
#      'OB': 'OP2', 
#      'O2P': 'OP2', 
#      'O3P': 'OP3'}, 
#     'D': {
#      'O5*': "O5'", 
#      'C5*': "C5'", 
#      'C4*': "C4'", 
#      'O4*': "O4'", 
#      "O1'": "O4'", 
#      'C3*': "C3'", 
#      'O3*': "O3'", 
#      'C2*': "C2'", 
#      'O2*': "O2'", 
#      'C1*': "C1'", 
#      'C7': 'C5M', 
#      'C5M': 'C7', 
#      'H1*': "H1'", 
#      'H2*1': "H2'", 
#      "H2'1": "H2'", 
#      "'H2'": "H2''", 
#      'H2*2': "H2''", 
#      "H2'2": "H2''", 
#      'H3*': "H3'", 
#      'H4*': "H4'", 
#      'H5*1': "H5'", 
#      "H5'1": "H5'", 
#      'H5*2': "H5''", 
#      "H5'2": "H5''", 
#      "'H5'": "H5''", 
#      "HO'2": "HO2'", 
#      'H5T': "HO5'", 
#      'H3T': "HO3'", 
#      'OA': 'OP1', 
#      'O1P': 'OP1', 
#      'OB': 'OP2', 
#      'O2P': 'OP2'},
#     'P' : {}}

resn_alias_dict ={
     'AN': 'RAN',
     'A5': 'RA5',
     'A3': 'RA3', 
     'A' : 'RA', 
     'CN': 'RCN', 
     'C5': 'RC5', 
     'C3': 'RC3', 
     'C' : 'RC', 
     'GN': 'RGN', 
     'G5': 'RG5', 
     'G3': 'RG3', 
     'G' : 'RG', 
     'UN': 'RUN', 
     'U5': 'RU5', 
     'U3': 'RU3', 
     'U': 'RU'}

def resn_alias(res):
 if res in resn_alias_dict:
  return resn_alias_dict[res]
 else:
  return res

def name_alias(mol, name):
 if name in name_alias_dict:
  return name_alias_dict[name]
 else:
  return name

#def name_alias(mol,name):
# if mol in name_alias_dict:
#  if name in name_alias_dict[mol]:
#   return name_alias_dict[mol][name]
#  else:
#   return name
# else:
#  return name



# def return_dict():
  
#     alias_dict ={'R': {
#                         'O5*': "O5'", 
#                         'C5*': "C5'", 
#                         'C4*': "C4'", 
#                         'O4*': "O4'", 
#                         "O1'": "O4'", 
#                         'C3*': "C3'", 
#                         'O3*': "O3'", 
#                         'C2*': "C2'", 
#                         'O2*': "O2'", 
#                         'C1*': "C1'", 
#                         'C7': 'C5M', 
#                         'C5M': 'C7', 
#                         'H1*': "H1'", 
#                         'H2*1': "H2'", 
#                         "H2'1": "H2'", 
#                         'H2*2': "H2''", 
#                         "H2'2": "H2''", 
#                         'H3*': "H3'", 
#                         'H4*': "H4'", 
#                         'H5*1': "H5'", 
#                         "H5'1": "H5'", 
#                         'H5*2': "H5''", 
#                         "H5'2": "H5''", 
#                         "HO'2": "HO2'", 
#                         'H5T': "HO5'", 
#                         'H3T': "HO3'", 
#                         'OA': 'OP1', 
#                         'O1P': 'OP1', 
#                         'OB': 'OP2', 
#                         'O2P': 'OP2', 
#                         'O3P': 'OP3'}, 
#                  'D': {
#                        'O5*': "O5'", 
#                        'C5*': "C5'", 
#                        'C4*': "C4'", 
#                        'O4*': "O4'", 
#                        "O1'": "O4'", 
#                        'C3*': "C3'", 
#                        'O3*': "O3'", 
#                        'C2*': "C2'", 
#                        'O2*': "O2'", 
#                        'C1*': "C1'", 
#                        'C7': 'C5M', 
#                        'C5M': 'C7', 
#                        'H1*': "H1'", 
#                        'H2*1': "H2'", 
#                        "H2'1": "H2'", 
#                        "'H2'": "H2''", 
#                        'H2*2': "H2''", 
#                        "H2'2": "H2''", 
#                        'H3*': "H3'", 
#                        'H4*': "H4'", 
#                        'H5*1': "H5'", 
#                        "H5'1": "H5'", 
#                        'H5*2': "H5''", 
#                        "H5'2": "H5''", 
#                        "'H5'": "H5''", 
#                        "HO'2": "HO2'", 
#                        'H5T': "HO5'", 
#                        'H3T': "HO3'", 
#                        'OA': 'OP1', 
#                        'O1P': 'OP1', 
#                        'OB': 'OP2', 
#                        'O2P': 'OP2'},
#                 'P' : {}}

#     return alias_dict

# def return_residue(res):


#     d = {'AN': 'RAN', 
#          'A5' : 'RA5',
#          'A3': 'RA3', 
#          'A': 'RA', 
#          'CN': 'RCN', 
#          'C5': 'RC5', 
#          'C3': 'RC3', 
#          'C': 'RC', 
#          'GN': 'RGN', 
#          'G5': 'RG5', 
#          'G3': 'RG3', 
#          'G': 'RG', 
#          'UN': 'RUN', 
#          'U5': 'RU5', 
#          'U3': 'RU3', 
#          'U': 'RU'}


#     if res in d:
     
#         return d[res]
    
#     else:
#         return res










