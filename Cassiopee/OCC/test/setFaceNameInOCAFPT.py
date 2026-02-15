# - setFaceNameInOCAF (pyTree) -
import OCC.PyTree as OCC
import Converter.PyTree as C

hook = OCC.createEmptyCAD()
OCC._addSphere(hook, (0,0,0), 1.)
OCC._setFaceNameInOCAF(hook, [1], 'exterior')

#OCC.printOCAF(hook)

t = OCC.meshAll(hook, hmax=0.1)
OCC._addOCAFCompoundNames(hook, t)
C.convertPyTree2File(t, 'out.cgns')
