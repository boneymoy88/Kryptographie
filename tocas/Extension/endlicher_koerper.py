import math
import random
import decimal

from Tocas import GanzzahlRestklassenring, Polynomring, PolynomringElement, Ring

import Extension.polynomring_extension
import Extension.ganzzahlrestklassenring_extension
from Extension import Primzahl
from Extension.PolynomRestklassenring import PolynomRestklassenring


def EndlicherKoerper(p, n=None, variablenname='x', primzahl_test=Primzahl.miller_rabin):
    if isinstance(p, Ring):
        if p.ist_endlicher_koerper():
            R = p
        else:
            raise TypeError(
                'Erstes Argument muss ein endlicher Körper sein')
    else:
        if not n:
            q = p
            for exp in range(2, int(math.log(q, 2)) + 1):
                basis = decimal.Decimal(
                    q) ** (decimal.Decimal(1)/decimal.Decimal(exp))
                if basis % 1 == 0 and primzahl_test(int(basis)):
                    p = int(basis)
                    n = exp
                    break
            else:
                if primzahl_test(p):
                    n = 1
                else:
                    raise TypeError(
                        'Erstes Argument muss Primzahl oder Primpotenz sein')
        elif not primzahl_test(p):
            raise TypeError('p ist keine Primzahl')
        R = GanzzahlRestklassenring(p)
    if n < 1:
        raise TypeError('n ist keine natürliche Zahl')

    if n == 1:
        return R

    RX = Polynomring(R, variablenname)

    koeffizienten = [RX.basisring.null] * (n + 1)

    random.seed(0)
    while True:
        koeffizienten[random.randint(0, n)] = RX.basisring.zufaellig()

        if koeffizienten[-1] == RX.basisring.null:
            koeffizienten[-1] += RX.basisring.eins

        f = PolynomringElement(koeffizienten, RX)
        if f.irreduzibel(primzahl_test):
            return PolynomRestklassenring(f)