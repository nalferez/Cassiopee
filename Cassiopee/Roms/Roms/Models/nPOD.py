# POD subclass - numpy implementation
from . import POD
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
import Compressor.PyTree as Compressor
import Converter.Internal as Internal
import numpy

class nPOD( POD.POD ):
    def __init__(self, name, type):
        super().__init__(name, type)

    def set(self, A, W=None, db=None, ref=None, variables=[]):
        self.A = A
        self.W = W
        self.db = db
        self.ref = ref
        self.variables = variables

    def savePhi(self):
        """Save Phi base modes."""
        if Cmpi.rank > 0: return None
        if self.Phi is None:
            raise ValueError("save: Phi modeds are not built.")
        t = C.newPyTree(['POD'])
        node = ["phi", self.Phi, [], 'DataArray_t']
        Compressor._packNode(node)
        t[2][1][2].append(node)
        C.convertPyTree2File(t, self.fileName, format='bin_hdf')
        return None

    # build Phi from matrix M keeping K modes with svd
    def buildPhi(self, K=-1):
        """Build phi modes with SVD."""
        if self.A is None:
            raise ValueError("build: A matrix is not present.")
        # build modes
        Phi, S, Vt = numpy.linalg.svd(self.A, full_matrices=False)
        # energy of each modes
        #energy = S**2 / numpy.sum(S**2)
        if K > 0: self.K = K
        else: self.K = Phi.shape[1]
        self.Phi = Phi[:, 0:self.K]
        return Phi, S, Vt

    def buildAndSaveCoords(self):
        """Build snapshots coords."""
        if self.A is None:
            raise ValueError("build: A matrix is not present.")
        if self.Phi is None:
            raise ValueError("build: Phi modeds are not built.")

        ncols = self.A.shape[1]
        self.coords = numpy.empty( (ncols, self.K), dtype=numpy.float64 )
        coords = numpy.empty( (self.K), dtype=numpy.float64 )
        for j in range(ncols):
            m = self.A[:,j].ravel('k')
            for i in range(self.K):
                c0 = numpy.dot(self.Phi[:,i], m)
                coords[i] = c0
            node = ["%05d"%j, coords, [], 'DataArray_t']
            Compressor._packNode(node)
            #print(node, flush=True)
            Filter.writeNodesFromPaths(self.fileName, 'CGNSTree/POD', node, format='bin_hdf')
            self.coords[j,:] = coords[:]

    def readCoords(self, j):
        node = Filter.readNodesFromPaths(self.fileName, ['CGNSTree/POD/%05d'%j])[0]
        Compressor._unpackNode(node)
        coord = node[1]
        return coord

    def instantiate(self, coords):
        """Return field vector from coords."""
        m = self.Phi @ coords
        return m

    def buildTree(self, coords):
        """Rebuild tree from field vector."""
        if self.db is None:
            raise ValueError("instantiate: db is missing.")
        if self.Phi is None:
            raise ValueError("instantiate: Phi modes are missing.")
        m = self.Phi @ coords
        cgnsName = self.db.dirName+'/%s'%self.ref+'.cgns'
        tref = C.convertFile2PyTree(cgnsName)
        for v in self.variables:
            v = v.split(':')
            if len(v) == 2:
                loc = v[0]; v = v[1]
            else: loc = 'nodes'; v = v[0]
            for z in Internal.getZones(tref):
                npts = C.getNPts(z)
                ncells = C.getNCells(z)
                pt0 = 0
                if loc == 'nodes':
                    FS = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
                    if FS is None:
                        FS = Internal.newFlowSolution(Internal.__FlowSolutionNodes__, 'Vertex', parent=z)
                    b = numpy.zeros( (npts, 1), dtype=numpy.float64)
                    Internal.newDataArray(v, value=b, parent=FS)
                    b[:,0] = m[pt0:pt0+npts]
                    pt0 += npts
                else:
                    FC = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
                    if FC is None:
                        FC = Internal.newFlowSolution(Internal.__FlowSolutionCenters__, 'CellCenter', parent=z)
                    b = numpy.zeros( (ncells, 1), dtype=numpy.float64)
                    Internal.newDataArray(v, value=b, parent=FC)
                    b[:,0] = m[pt0:pt0+ncells]
                    pt0 += ncells
        return tref

    def fetchTree(self, point):
        """Rebuild pyTree for given point."""
        from scipy.interpolate import RBFInterpolator
        q = self.db.query()
        xo = self.db.fetchPointVector(q)
        yo = self.coords
        # Create the RBF interpolator
        rbfInterpolator = RBFInterpolator(xo, yo, kernel='linear')
        # Define the input point
        inputPoint = numpy.empty( (1,xo.shape[1]), dtype=numpy.float64 )
        for c, i in enumerate(self.db.parameters):
            inputPoint[0,c] = point[i]
        # Evaluate the interpolator at the input point
        coords = rbfInterpolator(inputPoint)
        coords = coords.ravel('k')
        t = self.buildTree(coords)
        return t
