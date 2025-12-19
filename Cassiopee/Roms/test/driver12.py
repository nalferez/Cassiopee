# driver: parametric sketch for naca4 digits
# + volume meshing
import Roms.Driver as D

# Thickness
xx = D.Scalar('xx')
xx.range = [10., 14., 1.]

# max camber
M = D.Scalar('M')
M.range = [0., 5., 1.]

# camber position
P = D.Scalar('P')
P.range = [0., 5., 1.]

# Contrainte
D.Lt(P, M)

# Create points for exterior domain
P1 = D.Point('P1', (-10,-10,0))

P2 = D.Point('P2', (-10,10,0))

P3 = D.Point('P3', (10,10,0))

P4 = D.Point('P4', (10,-10,0))

# Create lines
line1 = D.Line('line1', P1, P2)
line2 = D.Line('line2', P2, P3)
line3 = D.Line('line3', P3, P4)
line4 = D.Line('line4', P4, P1)

# create profile
naca1 = D.Naca('naca1', M, P, xx)

# Create sketch
sketch1 = D.Sketch('sketch1', [naca1], h=[0.01,0.01,0.01])

sketch2 = D.Sketch('sketch2', [line1,line2,line3,line4], h=[1.,1.,1.])

# Create volume mesh
vol1 = D.Volume2D('vol1', [sketch1, sketch2], orders=[-1,-1])

# solve
D.DRIVER.solve()

#D.DRIVER.instantiate({'M':0., 'P': 0., 'xx': 12.})
#sketch1.writeCAD('out.step')
#exit(0)

#refMesh = sketch1.mesh(0.1,0.1,0.1)

#D.DRIVER.instantiate({'hauteur': 0.8})
#m = sketch1.rmesh2(refMesh,0.1,0.1,0.1)
#m = sketch1.rmesh(refMesh)
#m = sketch1.pmesh(refMesh)
#import Converter
#Converter.convertArrays2File(m, 'out.plt')

#refMesh = sketch1.Mesh(0.1,0.1,0.1)
#D.DRIVER.instantiate({'hauteur': 0.8})
#m = sketch1.Rmesh(refMesh,0.1,0.1,0.1)
#m = sketch1.Rmesh2(refMesh) # FAIL
#import Converter.PyTree as C
#C.convertPyTree2File(m, 'out.cgns')
#exit(0)

#m = sketch1.rmesh2(m)

#import Converter.PyTree as C
#C.convertPyTree2File(m, 'out.cgns')

# direct mesh
#refMesh = vol1.SetReferenceMesh(0.1,0.1,0.1)
#C.convertPyTree2File(refMesh, "out.cgns")

# indirect mesh
#D.DRIVER.instantiate({'hauteur': 0.7})
#dmesh = vol1.Dmesh()

#import CPlot.PyTree as CPlot
import CPlot
import Converter
import time
pt = D.DRIVER.walkDOE()
while pt is not None:
    m = vol1.mesh()
    #Converter.convertArrays2File(m, 'out.plt')
    #m = sketch1.mesh(0.01,0.01,0.01)
    CPlot.display(m)
    pt = D.DRIVER.walkDOE()
    time.sleep(0.5)