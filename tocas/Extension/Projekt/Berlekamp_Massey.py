from Tocas import Polynomring, PolynomringElement, Z
from Extension import endlicher_koerper
from Extension.PolynomRestklassenring import PolynomRestklassenring, PolynomRestklassenringElement
import Extension.Primzahl as Primzahl
import time

                        
R_X = Polynomring(Z)


def diskrepanz(polynom,seq,N,K):
    """Berechnet Diskrepanz der erechneten Sequenz und gegbenen Sequenz

    Args:
        polynom (PolynomringElement): Polynom das vorgegeben Sequenz errechnen soll
        seq                   (list): Die zu berechnende Sequenz
        N                      (int): Bisher betrachtete laenge der Sequenz
        K              (Endlicher K): Elemente aus Sequenz sind aus diesem endlichen Koerper
    Returns:
        integer: Diskrepanz zwischen tatsaechlicher Sequenz und durch Polynom errechnete Sequenz
    
    """
    total = K.null
    for i in range(polynom.grad + 1):
        d = polynom.koeffizienten[i] * seq[N-i]
        total = total + d
    return total

def reverseKoeffizienten(C):
    """Kehrt Koeffizientenliste um
    
    C = PolynomringElement([1,1,0,1],K) := 1 + x + x**3
    >>> reverseKoeffizienten(C)
    ---> [1,0,1,1]   := x + x**2 + x**3
    
    Args: 
        C (PolynomringElement): PolynomringElement, dessen Koeffizienten umgekehrt werden sollen
    Return:
        PolynomringElement: Gibt PolynomringElement mit Koeffizienten in umgekehrter Reihenfolge zurueck
    """
    anzahlKoeffizienten = C.koeffizienten.laenge
    reverseList = [0]*anzahlKoeffizienten
    for i in range(anzahlKoeffizienten):
        reverseList[anzahlKoeffizienten - 1 - i] = C.koeffizienten[i]
    return PolynomringElement(reverseList, C.ring)


def scalarMassey(seq,K):
    """Berlekamp Massey Algorithmus findet kleinstes Polynom welches gegebene Sequenz aus endlichem Koerper berechnet
    
    Args: 
        seq                  (list): Die zu berechnende Sequenz
        K (GanzzahlRestklassenring): Endlicher Koerper des zu Berechnenden Polynoms 
    
    Returns:
        PolynomringElement: Zurueckgegebenes Polynom berechnet gegebene Sequenz

    Raises:
        RuntimeError: Wenn K kein endlicher Koerper
                      Wenn K nicht uerber Primzahl
    """
    
    tStart = time.time() #: Laufzeitberechnung intialisieren
    
    if (not K.ist_endlicher_koerper()): # Testen ob K endlicher_koerper 
        raise RuntimeError("Berlekamp Massey nur auf endlichen Koerpern.")
    
    # Initialisierungen von
    # Polynomring: polyring
    # Polynomring Elemente: C,B
    
    P = Polynomring(K)

    if (not Primzahl.miller_rabin(P.basisring.modulus)):
        raise RuntimeError("Polynomring nicht ueber Primzahl")

    C = PolynomringElement([K.eins],P)
    B = PolynomringElement([K.eins],P)  # Rueckkopplungspolynom, bei letzten laengenaenderung des LFSR
    x = 1                               # Anzahl Iterationen seit letzter laengenaenderung des LSFR
    L = 0                               # aktuelle Laenge des LFSR
    b = K.eins                          # wert Diskrepanz der letzten laengenaenderung des LFSR
    N = 0                               # Aktuell betrachtete Elemente aus Sequenz

    
    for N in range(len(seq)):     # Durchlaufe alle Elemente der gegebenen Sequenz

        d = diskrepanz(C,seq,N,K) # Berechnet Diskrepanz der erechneten Sequenz und gegbenen Sequenz
        if (d == K.null):         # Falls die Diskrepanz null betraegt wurde Sequenz korrekt durch das Polynom C beschrieben
            x = x + 1             # somit kann mit dem naechsten Element fortgefuehrt werden
        else:                     # Falls nicht, muss das Polynom C angepasst werden  
            if 2*L > N:           # Groesse des aktuellen linearen Feedback Register reicht aus
                                  # falls L mehr als halb so gross wie aktuell betrachtete seq-laenge
                mB = B * (d/b).wert                       # (d/b) * B = mB
                hilfKoeff = mB * P.variable**x         
                C = C - PolynomringElement(hilfKoeff, P)  # C = C - (d/b) * B * P.variable**x
                x = x + 1                                 # keine laengenaenderung --> x erhoehen

            else:                 # Lineares Feedback Register zu klein und muss erhoeht werden

                mB = B * (d/b).wert 
                B = C
                hilfKoeff = mB * P.variable**x
                C = C - PolynomringElement(hilfKoeff, P) 
                L = N + 1 - L                            # Anpassung der laenge des LFSR
                b = d                                    # Bei Anpassung alte Diskrepanz in b
                x = 1                                    # zuruecksetzen der Iterationenen seit L vergroessert wurde
                
    C = reverseKoeffizienten(C)
    return C, L
    #return time.time() - tStart, P.basisring.modulus



