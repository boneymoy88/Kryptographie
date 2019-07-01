from Tocas import GanzzahlRestklassenring

from Extension import Primzahl


def _ganzzahlrestklassenring_ist_endlicher_koerper(k, primzahl_test=Primzahl.miller_rabin):
    return primzahl_test(k.modulus)


GanzzahlRestklassenring.ist_endlicher_koerper = _ganzzahlrestklassenring_ist_endlicher_koerper