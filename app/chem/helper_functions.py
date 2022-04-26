"""Copyright (c) 2022 VIKTOR B.V.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

VIKTOR B.V. PROVIDES THIS SOFTWARE ON AN "AS IS" BASIS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import re

def get_positions(block):
    """get specifice substrings with information from the large molblock string"""
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
    match_lis = [match.group() for match in matches]

    for count in enumerate(atoms): #assing positions to atoms
        count = count[0]
        atoms[count][atom_lis[count]]['x'] = float(match_lis[count * 3])
        atoms[count][atom_lis[count]]['y'] = float(match_lis[count * 3 + 1])
        atoms[count][atom_lis[count]]['z'] = float(match_lis[count * 3 + 2])

    pattern = re.compile(r'[^0\D][\s|\d]{2}\s?\d+\s{2}[^0\D]') #find bonds in text
    matches = pattern.finditer(repr(block))
    match_lis = [match.group() for match in matches]

    bonds = []
    i = 0
    for hit in match_lis: #find which atom the bond starts and ends at
        pattern = re.compile(r'\d+')
        matches = pattern.finditer(hit)
        match_lis = []
        for match in matches:
            string = match.group()
            if len(string) == 5:
                match_lis.append(string[0:2])
                match_lis.append(string[2:])
            else:
                match_lis.append(match.group())

        origin = int(match_lis[0]) - 1
        destination = int(match_lis[1]) - 1
        start = atoms[int(origin)][atom_lis[origin]] #get position of starting atom
        end = atoms[int(destination)][atom_lis[destination]] #gert position of ending atom
        num_bonds = int(match_lis[2])
        bonds.append({'start': start, 'end': end, 'num bonds': num_bonds})
        i += 1

    return atoms, atom_lis, bonds
