import pandas as pd
from viktor import Color
from viktor.core import ViktorController
from viktor.views import DataResult, DataGroup, DataView, DataItem, GeometryResult,GeometryView, SVGView, SVGResult
from viktor.geometry import Point, Sphere, Material, CircularExtrusion, Line

import chemlib
import pandas
import rdkit
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

from .helper_functions import get_positions


from .parametrization import ChemistryParametrization


class ChemistryController(ViktorController):
    label = 'Chemistry'
    parametrization = ChemistryParametrization

    @DataView('Element properties', duration_guess=1)
    def get_properties(self, params, **kwargs):
        element = chemlib.chemistry.Element(params.elements.input)

        res = DataGroup(DataItem(label='Atomic number', value=element.AtomicNumber),
                        DataItem(label='Atomic mass', value=element.AtomicMass),
                        DataItem(label='Atomic radius', value=element.AtomicRadius),
                        DataItem(label='density', value=element.Density),
                        DataItem(label='melting point', value=element.MeltingPoint),
                        DataItem(label='boiling point', value=element.BoilingPoint))

        return DataResult(res)
    def to_chemlib_compound(self,strings):
        strings = list(filter(lambda val: val != '' and val != None, strings))
        # reactants = list(filter(lambda val: val != None, reactants))

        compounds = []
        for string in strings:
            compounds.append(chemlib.Compound(string))

        return compounds

    @DataView('Reaction', duration_guess=1)
    def get_reaction(self, params, **kwargs):
        reactants = [params.reactions.reactants[0]['first'],
                     params.reactions.reactants[0]['second'],
                     params.reactions.reactants[0]['third'],
                     params.reactions.reactants[0]['fourth']]
        products = [params.reactions.products[0]['first'],
                    params.reactions.products[0]['second'],
                     params.reactions.products[0]['third'],
                    params.reactions.products[0]['fourth']]

        reactants = self.to_chemlib_compound(reactants)
        products = self.to_chemlib_compound(products)

        reaction = chemlib.Reaction(reactants, products)
        reaction.balance()
        reaction = reaction.formula

        result = DataGroup(DataItem(label='Balanced reaction', value=reaction))

        return DataResult(result)

    @SVGView('SMILE', duration_guess=1)
    def get_SMILE(self, params, **kwargs):
        molecule = Chem.MolFromSmiles(params.smile.smile)
        fig = Draw.MolToFile(molecule, 'molecule_drawing.svg')

        return SVGResult.from_path('molecule_drawing.svg')

    def get_positions(self, block):
        block = block.split('\n')
        amounts = block[3].split(' ')
        amounts = list(filter(lambda val: val != '', amounts))
        num_molecules = int(amounts[0])
        num_bonds = int(amounts[1])
        bonds = block[-(num_bonds + 2):-2]
        atoms = block[4:-(num_bonds + 2)]
        mol_list = []

        temporary = []
        for line in atoms:
            line = line[3:-37]
            mol_list.append(line[-1])
            line = line[:-2]
            line = line.replace('   ', '    ')
            line = line.split('    ')
            temporary.append(line)
        atoms = temporary

        temporary = []
        for el in atoms:
            inter = []
            for str in el:
                str = float(str)
                inter.append(str)
            temporary.append(inter)
        atoms = temporary

        temporary = []
        for bond in bonds:
            bond = bond.replace('  ', ' ')
            bond = bond[1:]
            bond = bond.split(' ')
            temporary.append([atoms[int(bond[0]) - 1], atoms[int(bond[1]) - 1], int(bond[2])])
        bonds = temporary

        return atoms,mol_list, bonds

    @GeometryView('3D', duration_guess=2)
    def get_Geom(self, params, **kwargs):
        double_bond_separation = 0.2
        molecule = Chem.MolFromSmiles(params.smile.smile)
        molecule = Chem.AddHs(molecule)
        AllChem.EmbedMolecule(molecule)
        block = Chem.MolToMolBlock(molecule)

        atoms,atom_list, bonds= get_positions(block)

        geometries = []
        i = 0
        for dict in atoms:
            atom = atom_list[i]
            point = Point(dict[atom]['x'],dict[atom]['y'],dict[atom]['z'])
            if atom == 'H':
                geometries.append(Sphere(point, 0.2, material=Material('white',color= Color.white())))
            elif atom == 'O':
                geometries.append(Sphere(point, 0.5, material=Material('red', color=Color(200,20,20))))
            elif atom == 'N':
                geometries.append(Sphere(point, 0.5, material=Material('blue', color=Color(20,20,200))))
            else:
                geometries.append(Sphere(point, 0.5, material=Material('black', color = Color.viktor_black())))
            i += 1

        for bond in bonds:
            if bond['num bonds'] == 2:
                shifts = np.array([bond['start']['x'] - bond['end']['x'],
                                   bond['start']['y'] - bond['end']['y'],
                                   bond['start']['z'] - bond['end']['z']])
                shifts = (shifts / sum(abs(shifts)))*double_bond_separation

                point_start_1 = Point(bond['start']['x'] + shifts[1],
                                      bond['start']['y'] + shifts[2],
                                      bond['start']['z'] + shifts[0])
                point_end_1 = Point(bond['end']['x'] + shifts[1],
                                    bond['end']['y'] + shifts[2],
                                    bond['end']['z'] + shifts[0])
                left_line_1 = Line(point_start_1, point_end_1)
                cylinder_1 = CircularExtrusion(diameter=0.1, line=left_line_1,
                                             material=Material('black', color=Color.viktor_black()))
                geometries.append(cylinder_1)

                point_start_2 = Point(bond['start']['x'] - shifts[1],
                                      bond['start']['y'] - shifts[2],
                                      bond['start']['z'] - shifts[0])
                point_end_2 = Point(bond['end']['x'] - shifts[1],
                                    bond['end']['y'] - shifts[2],
                                    bond['end']['z'] - shifts[0])
                left_line_2 = Line(point_start_2, point_end_2)
                cylinder_2 = CircularExtrusion(diameter=0.1, line=left_line_2,
                                               material=Material('black', color=Color.viktor_black()))
                geometries.append(cylinder_2)
            else:
                point_start = Point(bond['start']['x'], bond['start']['y'], bond['start']['z'])
                point_end = Point(bond['end']['x'], bond['end']['y'], bond['end']['z'])
                left_line = Line(point_start, point_end)
                cylinder = CircularExtrusion(diameter=0.2, line=left_line,material=Material('black', color = Color.viktor_black()))
                geometries.append(cylinder)

        return GeometryResult(geometries)

