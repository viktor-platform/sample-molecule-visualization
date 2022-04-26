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

from viktor.parametrization import Parametrization
from viktor.parametrization import Section
from viktor.parametrization import Table
from viktor.parametrization import TextField

class ChemistryParametrization(Parametrization):
    """For providing the right input fields to users"""
    elements = Section('Element input')
    elements.element = TextField('Which element')

    reactions = Section('Reaction input')
    reactions.reactants = Table('reactions')
    reactions.reactants.first = TextField('reactant 1')
    reactions.reactants.second = TextField('reactant 2')
    reactions.reactants.third = TextField('reactant 3')
    reactions.reactants.fourth = TextField('reactant 4')
    reactions.products = Table('products')
    reactions.products.first = TextField('product 1')
    reactions.products.second = TextField('product 2')
    reactions.products.third = TextField('product 3')
    reactions.products.fourth = TextField('product 4')

    smile = Section('SMILE input')
    smile.smile = TextField('SMILE', flex= 100)
