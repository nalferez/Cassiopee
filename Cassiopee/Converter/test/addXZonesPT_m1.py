# - addXZones (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# Cree le fichier test
if Cmpi.rank == 0:
    a = G.cart((0,0,0), (1,1,1), (10,10,10))
    b = G.cart((9,0,0), (1,1,1), (10,10,10))
    t = C.newPyTree(['Base',a,b])
    C.convertPyTree2File(t, LOCAL+'/test.cgns')
Cmpi.barrier()

# Relit des zones par procs
t = Cmpi.convertFile2SkeletonTree(LOCAL+'/test.cgns')
(t, dic) = Distributor2.distribute(t, NProc=Cmpi.size, algorithm='fast')
t = Cmpi.readZones(t, LOCAL+'/test.cgns', rank=Cmpi.rank)
Cmpi._convert2PartialTree(t)

# Cree le bbox tree
tb = Cmpi.createBBoxTree(t)
# Cree le graph
graph = Cmpi.computeGraph(tb)
# Add X Zones
t = Cmpi.addXZones(t, graph)
if Cmpi.rank == 0:
    #Internal.printTree(t)
    test.testT(t, 1)
