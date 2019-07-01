import copy
import math

from Tocas import Polynomring, PolynomringElement

def _polynomring_ExtGGT(a: PolynomringElement, b: PolynomringElement):
    if not (isinstance(a, PolynomringElement) and isinstance(b, PolynomringElement)):
        raise TypeError('Argumente nicht vom Typ PolynomringElement')

    if a.basisring != b.basisring:
        raise ValueError(
            'PolynomringElement haben nicht den gleichen Grundring')

    u, v = a.ring.eins, a.ring.null
    s, t = a.ring.null, a.ring.eins
    while b != a.ring.null:
        q = a // b
        a, b = b, a - q * b
        u, s = s, u - q * s
        v, t = t, v - q * t
    return a.ring.eins, u // a, v // a

def _polynomring_ext_ggt(self: Polynomring, a: PolynomringElement, b: PolynomringElement):
    if not (isinstance(a, PolynomringElement) and isinstance(b, PolynomringElement)):
        raise TypeError('Argumente nicht vom Typ PolynomringElement')

    if a.ring != self or b.ring != self:
        raise TypeError('PolynomringElement nicht im Polynomring')

    return Polynomring.ExtGGT(a, b)

def _polynomring_ist_endlicher_koerper(_):
    return False


Polynomring.ExtGGT = staticmethod(_polynomring_ExtGGT)
Polynomring.ext_ggt = _polynomring_ext_ggt
Polynomring.ist_endlicher_koerper = _polynomring_ist_endlicher_koerper


def _polynom_div(a: PolynomringElement, b: PolynomringElement):
    """Polynomdivision ohne Rest"""
    a_koeffizienten = copy.copy(a.koeffizienten.koeffizienten)
    b_koeffizienten = b.koeffizienten.koeffizienten
    q_koeffizienten = []

    while a_koeffizienten and len(a_koeffizienten) >= len(b_koeffizienten):
        t = a_koeffizienten[-1] / b_koeffizienten[-1]
        q_koeffizienten = [t] + q_koeffizienten

        for i in range(len(b_koeffizienten)):
            a_koeffizienten[len(a_koeffizienten) -
                            len(b_koeffizienten) + i] -= t * b_koeffizienten[i]

        del a_koeffizienten[-1]

    q = PolynomringElement(q_koeffizienten, a.ring)
    r = PolynomringElement(a_koeffizienten, a.ring)

    return q, r


def _polynom_floordiv(a: PolynomringElement, b: PolynomringElement):
    return _polynom_div(a, b)[0]


def _polynom_mod(a: PolynomringElement, b: PolynomringElement):
    return _polynom_div(a, b)[1]


PolynomringElement.__floordiv__ = _polynom_floordiv
PolynomringElement.__mod__ = _polynom_mod


def _polynom_ExtGGT(a: PolynomringElement, b: PolynomringElement, norm=True):
    if not (isinstance(a, PolynomringElement) and isinstance(b, PolynomringElement)):
        raise TypeError('Argumente nicht vom Typ PolynomringElement')

    if a.basisring != b.basisring:
        raise ValueError(
            'PolynomringElement haben nicht den gleichen Grundring')

    u, v = a.ring.eins, a.ring.null
    s, t = a.ring.null, a.ring.eins
    while b != a.ring.null:
        q = a // b
        a, b = b, a - q * b
        u, s = s, u - q * s
        v, t = t, v - q * t

    if norm:
        a, u, v = a.ring.eins, u // a, v // a
    return a, u, v


def _polynom_ext_ggt(self: Polynomring, a: PolynomringElement, b: PolynomringElement, norm=True):
    if not (isinstance(a, PolynomringElement) and isinstance(b, PolynomringElement)):
        raise TypeError('Argumente nicht vom Typ PolynomringElement')

    if a.ring != self or b.ring != self:
        raise TypeError('PolynomringElement nicht im Polynomring')

    return Polynomring.ExtGGT(a, b, norm)


Polynomring.ExtGGT = staticmethod(_polynom_ExtGGT)
Polynomring.ext_ggt = _polynom_ext_ggt



def _primzahlen(n, primzahl_test):
    divisors = [d for d in range(2, n+1) if n % d == 0]
    return [d for d in divisors if primzahl_test(d)]


def _polynom_irreduzibel(f: PolynomringElement, primzahl_test=None):
    import Extension.Primzahl as primzahl
    from Extension.PolynomRestklassenring import PolynomRestklassenring, PolynomRestklassenringElement

    if not primzahl_test:
        primzahl_test = primzahl.miller_rabin
    
    if not f.basisring.ist_endlicher_koerper():
        raise TypeError('Bassisring muss endlicher Körper sein')
    if f.grad == 0:
        return False

    KX_f = PolynomRestklassenring(f)
    var = PolynomRestklassenringElement(f.ring.variable, KX_f)

    n = [int(f.grad / p) for p in _primzahlen(f.grad, primzahl_test)]

    if isinstance(f.basisring, PolynomRestklassenring):
        q = f.basisring.modulus.basisring.modulus ** f.basisring.modulus.grad
    else:
        q = f.basisring.modulus

    for i in range(0, len(n)):
        h = (((var ** q) ** n[i]) - var)
        g, _, _ = Polynomring.ExtGGT(f, h.wert, False)
        if g.grad != 0:
            return False

    g = ((var ** (q ** f.grad)) - var)
    if g == var.ring.null:
        return True
    return False

PolynomringElement.__floordiv__ = _polynom_floordiv
PolynomringElement.__mod__ = _polynom_mod
PolynomringElement.irreduzibel = _polynom_irreduzibel
