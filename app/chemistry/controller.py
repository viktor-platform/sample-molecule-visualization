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


from .parametrization import chemistryParametrization


class chemistryController(ViktorController):
    label = 'chemistry'
    parametrization = chemistryParametrization

    @DataView('Element properties', duration_guess=1)
    def get_properties(self, params, **kwargs):
        element = chemlib.chemistry.Element(params.elements.input)
        num = element.AtomicNumber
        mass = element.AtomicMass
        radius = element.AtomicRadius
        density = element.Density
        melt = element.MeltingPoint
        boil = element.BoilingPoint

        res = DataGroup(DataItem(label='Atomic number', value=num),
                        DataItem(label='Atomic mass', value=mass),
                        DataItem(label='Atomic radius', value=radius),
                        DataItem(label='density', value=density),
                        DataItem(label='melting point', value=melt),
                        DataItem(label='boiling point', value=boil))

        return DataResult(res)

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

        reactants = list(filter(lambda val: val != '', reactants))
        reactants = list(filter(lambda val: val != None, reactants))
        products = list(filter(lambda val: val != '', products))
        products = list(filter(lambda val: val != None, products))

        compounds = []
        for reactant in reactants:
            compounds.append(chemlib.Compound(reactant))
        reactants = compounds

        compounds = []
        for product in products:
            compounds.append(chemlib.Compound(product))
        products = compounds

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

    @GeometryView('3D', duration_guess=1)
    def get_Geom(self, params, **kwargs):
        double_bond_separation = 0.2
        molecule = Chem.MolFromSmiles(params.smile.smile)
        molecule = Chem.AddHs(molecule)
        AllChem.EmbedMolecule(molecule)
        block = Chem.MolToMolBlock(molecule)

        block,mol_list, bonds = self.get_positions(block)

        geometries = []
        i = 0
        for positions in block:
            point = Point(positions[0], positions[1], positions[2])
            if mol_list[i] == 'H':
                geometries.append(Sphere(point, 0.2, material=Material('white',color= Color.white())))
            elif mol_list[i] == 'O':
                geometries.append(Sphere(point, 0.5, material=Material('red', color=Color(200,20,20))))
            elif mol_list[i] == 'N':
                geometries.append(Sphere(point, 0.5, material=Material('blue', color=Color(20,20,200))))
            else:
                geometries.append(Sphere(point, 0.5, material=Material('black', color = Color.viktor_black())))
            i += 1

        for bond in bonds:
            if bond[2] == 2:
                shifts = np.array([bonds[0][0][0] - bonds[0][1][0], bonds[0][0][1] - bonds[0][1][1], bonds[0][0][2] - bonds[0][1][2]])
                shifts = (shifts / sum(abs(shifts)))*double_bond_separation
                point_start_1 = Point(bond[0][0] + shifts[0], bond[0][1] + shifts[1], bond[0][2] + shifts[2])
                point_end_1 = Point(bond[1][0] + shifts[0], bond[1][1] + shifts[1], bond[1][2] + shifts[2])
                left_line_1 = Line(point_start_1, point_end_1)
                cylinder_1 = CircularExtrusion(diameter=0.1, line=left_line_1,
                                             material=Material('black', color=Color.viktor_black()))
                geometries.append(cylinder_1)

                point_start_2 = Point(bond[0][0] - shifts[0], bond[0][1] - shifts[1], bond[0][2] - shifts[2])
                point_end_2 = Point(bond[1][0] - shifts[0], bond[1][1] - shifts[1], bond[1][2] - shifts[2])
                left_line_2 = Line(point_start_2, point_end_2)
                cylinder_2 = CircularExtrusion(diameter=0.1, line=left_line_2,
                                               material=Material('black', color=Color.viktor_black()))
                geometries.append(cylinder_2)
            else:
                point_start = Point(bond[0][0], bond[0][1], bond[0][2])
                point_end = Point(bond[1][0], bond[1][1], bond[1][2])
                left_line = Line(point_start, point_end)
                cylinder = CircularExtrusion(diameter=0.2, line=left_line,material=Material('black', color = Color.viktor_black()))
                geometries.append(cylinder)

        return GeometryResult(geometries)

