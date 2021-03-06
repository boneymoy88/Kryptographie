import project_path
from Tocas import Polynomring, PolynomringElement, Z, GanzzahlRestklassenring, GanzzahlRestklassenringElement
from Extension.Projekt import Berlekamp_Massey as FR
from Extension import endlicher_koerper
from Extension.PolynomRestklassenring import PolynomRestklassenring, PolynomRestklassenringElement
import Extension.polynomring_extension
import Poly
import Modular
import Algorithms

F = GanzzahlRestklassenring(7)
FX = Polynomring(F)
g = PolynomringElement([1,0,1],FX)
FX_g = PolynomRestklassenring(g)

sequenz = [1,4,2,4,5,2,2,1,3,4]
print(g.irreduzibel())
print(FX_g.ist_endlicher_koerper())
C = PolynomRestklassenringElement([1], FX_g)
print(FX_g)
#FR.berlekamp_massey(sequenz,FX_g)
field = Modular.Modular(101)
polyRing = Poly.Poly(field)
B = [field.one()]
print(B)
print(polyRing.timesXPower(B,6))
