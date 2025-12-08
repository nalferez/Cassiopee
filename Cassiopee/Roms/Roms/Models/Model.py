# The Model class
import Converter.PyTree as C
import Converter.Mpi as Cmpi

class Model:
    def __init__(self, name, type):
        # model name
        if name[-4:] == ".mod": name = name[:-4]
        if name[-5:] == ".mod/": name = name[:-5]
        self.name = name
        # model type
        self.type = type
        # model file
        self.fileName = name+'.mod'
        # pytree link
        self.db = None # link to db
        self.ref = None # link to ref name in db
        self.variables = None # variable list
        self.filter = None # sub indices in global matrix
        # matrix
        self.A = None # distributed
        self.W = None # distributed

    def load(self):
        return None

    def save(self):
        return None

    def build(self):
        return None

    def instantiate(self):
        return None