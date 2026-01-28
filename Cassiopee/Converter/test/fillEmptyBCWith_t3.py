# - fillEmptyBC (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# BE
a = G.cartHexa((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
# ajoute en tant que GC
C._fillEmptyBCWith(a, 'wallv', 'BCWallViscous')
test.testT(a, 1)

# ME to complete
