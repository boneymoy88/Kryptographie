from Tocas import Ganzzahlring


def _ganzzahlring_ist_endlicher_koerper(_):
    return False


Ganzzahlring.ist_endlicher_koerper = _ganzzahlring_ist_endlicher_koerper
