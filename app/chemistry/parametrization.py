from viktor.parametrization import Parametrization, Section, NumberField, TextField, Table, ToggleButton, Lookup, IsNotEqual


class ChemistryParametrization(Parametrization):
    elements = Section('Element input')
    elements.input = TextField('Which element')

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

