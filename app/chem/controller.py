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

import chemlib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from viktor import Color
from viktor.core import ViktorController
from viktor.geometry import CircularExtrusion
from viktor.geometry import Line
from viktor.geometry import Material
from viktor.geometry import Point
from viktor.geometry import Sphere
from viktor.views import DataGroup
from viktor.views import DataItem
from viktor.views import DataResult
from viktor.views import DataView
from viktor.views import GeometryResult
from viktor.views import GeometryView
from viktor.views import SVGResult
from viktor.views import SVGView

from .helper_functions import get_positions
from .parametrization import ChemistryParametrization


class ChemistryController(ViktorController):
    """for visualization and calculation"""
    label = 'Chemistry'
    parametrization = ChemistryParametrization

    @DataView('Element properties', duration_guess=1)
    def get_properties(self, params, **kwargs):
        """view for getting element properties"""
        element = chemlib.chemistry.Element(params.elements.element)

        res = DataGroup(DataItem(label='Atomic number', value=element.AtomicNumber),
                        DataItem(label='Atomic mass', value=element.AtomicMass),
                        DataItem(label='Atomic radius', value=element.AtomicRadius),
                        DataItem(label='density', value=element.Density),
                        DataItem(label='melting point', value=element.MeltingPoint),
                        DataItem(label='boiling point', value=element.BoilingPoint))

        return DataResult(res)

    def to_chemlib_compound(self,strings):
        """make chemical formulas from strings"""
        strings = list(filter(lambda val: val not in ('',None), strings))

        compounds = [chemlib.Compound(string) for string in strings]

        return compounds

    @DataView('Reaction', duration_guess=1)
    def get_reaction(self, params, **kwargs):
        """for constructing a reaction with the given in/outputs"""
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
    def get_smile_2d(self, params, **kwargs):
        """for visualizing a drawing corresponding to a smile"""
        molecule = Chem.MolFromSmiles(params.smile.smile)
        Draw.MolToFile(molecule, 'molecule_drawing.svg')

        return SVGResult.from_path('molecule_drawing.svg')

    @GeometryView('3D', duration_guess=2)
    def get_3d(self, params, **kwargs):
        """to make a 3D representation corresponding to a smile"""
        double_bond_separation = 0.2
        molecule = Chem.MolFromSmiles(params.smile.smile)
        molecule = Chem.AddHs(molecule)
        AllChem.EmbedMolecule(molecule)
        block = Chem.MolToMolBlock(molecule)
        atoms,atom_list, bonds = get_positions(block)

        geometries = []
        i = 0
        for atom_dict in atoms:
            atom = atom_list[i]
            point = Point(atom_dict[atom]['x'],atom_dict[atom]['y'],atom_dict[atom]['z'])
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
                cylinder = CircularExtrusion(diameter=0.2, line=left_line,
                                             material=Material('black', color=Color.viktor_black()))
                geometries.append(cylinder)

        return GeometryResult(geometries)
