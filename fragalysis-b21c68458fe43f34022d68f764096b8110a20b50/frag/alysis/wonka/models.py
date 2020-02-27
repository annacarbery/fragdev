import uuid

class Object(object):

    def __init__(self, value_array, object_desc):
        """
        An Object will have a type and a value arrays
        :param value_array: e.g. x,y,z coords - what is clustered
        :param object_desc: e.g. "Water" or "H-bond acceptor"
        """
        self.value_array = value_array
        self.object_desc = object_desc
        self.uuid = str(uuid.uuid4())
        self.cluster = -1

class Owner(object):

    def __init__(self, object_list, title):
        """
        An Owner will own multiple objects. E.g. PDB 4CUP owns waters.
        :param object_list: the list of objects it owns
        :param title: the title of the object
        """
        self.object_list = object_list
        self.title = title
