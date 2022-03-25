import re

def get_positions(block):
    block = repr(block)
    atoms = []
    atom_lis = []
    pattern = re.compile(r'\s[^0-9\sM]\s') #find atom types in text block

    matches = pattern.finditer(block)

    for match in matches:
        atom_lis.append(match.group()[1])
        atoms.append({match.group()[1]: {'x': None, 'y': None, 'z': None}})


    pattern = re.compile(r'-?\d+\.\d{4}') #find positions in text

    matches = pattern.finditer(block)
    match_lis = []
    for match in matches:
        match_lis.append(match.group())

    for i in range(len(atoms)): #assing positions to atoms
        atoms[i][atom_lis[i]]['x'] = float(match_lis[i * 3])
        atoms[i][atom_lis[i]]['y'] = float(match_lis[i * 3 + 1])
        atoms[i][atom_lis[i]]['z'] = float(match_lis[i * 3 + 2])

    pattern = re.compile(r'[^0\D][\s|\d]{2}\s?\d+\s{2}[^0\D]') #find bonds in text
    matches = pattern.finditer(repr(block))
    match_lis = []
    for match in matches:
        match_lis.append(match.group())

    bonds = []
    i = 0
    for match in match_lis: #find which atom the bond starts and ends at
        pattern = re.compile(r'\d+')
        matches = pattern.finditer(match)
        match_lis = []
        for match in matches:
            match_lis.append(match.group())

        origin = int(match_lis[0]) - 1
        destination = int(match_lis[1]) - 1
        start = atoms[int(origin)][atom_lis[origin]] #get position of starting atom
        end = atoms[int(destination)][atom_lis[destination]] #gert position of ending atom
        num_bonds = int(match_lis[2])
        bonds.append({'start': start, 'end': end, 'num bonds': num_bonds})
        i += 1

    return atoms, atom_lis, bonds