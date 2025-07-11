# - extractionIBM a la paroi (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import Post.PyTree as P
import Dist2Walls.PyTree as DTW
import Transform.PyTree as T
import Initiator.PyTree as I
import Converter.Internal as Internal
import Connector.ToolboxIBM as IBM
import KCore.test as test
import numpy
test.TOLERANCE = 1.e-8

N = 41
a = G.cart((0,0,0),(1./(N-1),1./(N-1),1./(N-1)),(N,N,N))
xm = 0.5*N/(N-1)
s = D.sphere((xm,xm,xm),0.1,N=20)
s = C.convertArray2Tetra(s); s = G.close(s)
t = C.newPyTree(['Base', a])

# Blanking
bodies = [[s]]
BM = numpy.array([[1]],Internal.E_NpyInt)
t = X.blankCells(t,bodies,BM,blankingType='center_in')
X._setHoleInterpolatedPoints(t,depth=-2)
# Dist2Walls
DTW._distance2Walls(t,[s],type='ortho',loc='centers',signed=1)
t = C.center2Node(t,'centers:TurbulentDistance')
# Gradient de distance localise en centres => normales
t = P.computeGrad(t, 'TurbulentDistance')
C._initVars(t,"centers:Density",1.)
C._initVars(t,"centers:VelocityX",0.2)
C._initVars(t,"centers:VelocityY",0.)
C._initVars(t,"centers:VelocityZ",0.)
C._initVars(t,"centers:Temperature",1.)
tc = C.node2Center(t)

tb = C.newPyTree(['Base', s])
C._addState(tb, 'EquationDimension',3)
C._addState(tb, 'GoverningEquations', 'NSTurbulent')

X._setIBCData(t, tc, loc='centers', storage='inverse', bcType=0)

#test avec arbre tc compact
vars=['Density','VelocityX','VelocityY','VelocityZ','Temperature']

X.miseAPlatDonorTree__(t, tc, graph=None)
# attention compact=0 car t n est pas compacte
X._setInterpTransfers(t,tc, bcType=0,varType=2,variablesIBC=vars,compact=0,compactD=1)
z = IBM.extractIBMWallFields(tc, tb=tb)
test.testT(z,1)
#
Internal._rmNodesFromName(t,"Parameter_int")
Internal._rmNodesFromName(t,"Parameter_real")
tc = C.node2Center(t)
X._setIBCData(t, tc, loc='centers', storage='inverse', bcType=3)
X.miseAPlatDonorTree__(t, tc, graph=None)
X._setInterpTransfers(t, tc, bcType=3, varType=2,variablesIBC=vars,compact=0,compactD=1)
tb_out = IBM.extractIBMWallFields(tc, tb=tb)
test.testT(tb_out,2)
