# - model -
import Roms.DB.DataBase as DataBase
import Roms.Models.nPOD as nPOD
import Converter.PyTree as C

db = DataBase.DataBase('NACA1.db', mode='r')
q = db.query()
A = db.fetchMatrix(q, variables=['centers:Density'])

mod = nPOD.nPOD('density', type='npod')
mod.set(A, db=db, ref="ref1", variables=['centers:Density'])

# Build phi
mod.buildPhi()
mod.savePhi()
mod.buildAndSaveCoords()

# instantiate an existing snapshot
coords = mod.readCoords(0)
t = mod.buildTree(coords)
#C.convertPyTree2File(t, 'out.cgns')

# interpolate snapshot
t = mod.fetchTree({'Mach':0.75, 'AoA': -1.38})
C.convertPyTree2File(t, 'out.cgns')
