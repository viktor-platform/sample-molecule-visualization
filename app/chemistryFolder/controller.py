from viktor.core import ViktorController

class ChemistryFolderController(ViktorController):
    label = 'Chemistry Folder'
    children = ['Chemistry']
    show_children_as = 'Cards'

viktor_convert_entity_field = True