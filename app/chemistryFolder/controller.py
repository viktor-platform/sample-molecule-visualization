from viktor.core import ViktorController

class chemistryFolderController(ViktorController):
    label = 'chemistry Folder'
    children = ['chemistry']
    show_children_as = 'Cards'

    viktor_convert_entity_field = True