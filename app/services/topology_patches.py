def apply_topology_patches(ff_instance, mol_type: str):
    """
    Univerzální patch pro ionty z nabídky.
    Zajišťuje, že AMBER_topology.py vždy najde potřebná fyzikální data.
    """
    if not hasattr(ff_instance, 'b'): ff_instance.b = {}
    if not hasattr(ff_instance, 'units'): ff_instance.units = {}
    if not hasattr(ff_instance, 'types'): ff_instance.types = {}

    # --- ČÁST A: KOMPLETNÍ DATABÁZE IONTOŮ (z tvé UI nabídky) ---
    # Hodnoty R a eps jsou standardní Amber parametry (Joung-Cheatham)
    standard_ions = {
        # Kationty (Cations)
        'Li+': {'real': 'LI', 'atom': 'LI', 'charge': 1.0, 'mass': 6.941, 'at_num': 3, 'R': 0.0816, 'eps': 0.0765},
        'Na+': {'real': 'NA', 'atom': 'NA', 'charge': 1.0, 'mass': 22.990, 'at_num': 11, 'R': 0.1369, 'eps': 0.3558},
        'K+': {'real': 'K', 'atom': 'K', 'charge': 1.0, 'mass': 39.098, 'at_num': 19, 'R': 0.1705, 'eps': 0.3558},
        'Rb+': {'real': 'RB', 'atom': 'RB', 'charge': 1.0, 'mass': 85.468, 'at_num': 37, 'R': 0.1813, 'eps': 0.6276},
        'Cs+': {'real': 'CS', 'atom': 'CS', 'charge': 1.0, 'mass': 132.905, 'at_num': 55, 'R': 0.1977, 'eps': 0.8786},
        'Mg2+': {'real': 'MG', 'atom': 'MG', 'charge': 2.0, 'mass': 24.305, 'at_num': 12, 'R': 0.1180, 'eps': 0.0400},
        'Ca2+': {'real': 'CA', 'atom': 'CA', 'charge': 2.0, 'mass': 40.078, 'at_num': 20, 'R': 0.1470, 'eps': 0.4500},
        'Zn2+': {'real': 'ZN', 'atom': 'ZN', 'charge': 2.0, 'mass': 65.380, 'at_num': 30, 'R': 0.0926, 'eps': 0.0523},

        # Anionty (Anions)
        'F-': {'real': 'F', 'atom': 'F', 'charge': -1.0, 'mass': 18.998, 'at_num': 9, 'R': 0.2303, 'eps': 0.0031},
        'Cl-': {'real': 'CL', 'atom': 'CL', 'charge': -1.0, 'mass': 35.453, 'at_num': 17, 'R': 0.2513, 'eps': 0.1272},
        'Br-': {'real': 'BR', 'atom': 'BR', 'charge': -1.0, 'mass': 79.904, 'at_num': 35, 'R': 0.2470, 'eps': 0.3765},
        'I-': {'real': 'I', 'atom': 'I', 'charge': -1.0, 'mass': 126.904, 'at_num': 53, 'R': 0.2860, 'eps': 0.1674},
    }

    # Propojíme aliasy a injektujeme chybějící data
    for alias, data in standard_ions.items():
        # 1. Pokud silové pole iont vůbec nezná (časté u vody) -> Injekce
        if alias not in ff_instance.units and data['real'] not in ff_instance.units:
            ff_instance.units[alias] = {
                'atoms': {
                    'name': [data['atom']],
                    'type': [data['atom']],
                    'charge': [data['charge']],
                    'mass': [data['mass']],
                    'at_num': [data['at_num']],
                    'R': [data['R']],
                    'eps': [data['eps']]
                },
                'bonds': []
            }
            ff_instance.types[data['atom']] = data['atom']

        # 2. Pokud zná jen jeden název (třeba NA), ale hledáme alias (Na+) -> Bridge
        elif data['real'] in ff_instance.units and alias not in ff_instance.units:
            ff_instance.units[alias] = ff_instance.units[data['real']]
            # Zajistíme registraci typu
            atom_type = ff_instance.units[data['real']]['atoms']['type'][0]
            ff_instance.types[atom_type] = atom_type

    return ff_instance