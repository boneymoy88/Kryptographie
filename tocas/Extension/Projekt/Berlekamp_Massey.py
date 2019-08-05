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

def skalarPolyMultiplikation(polynom,skalar):
    """Hilfsfunktion um Multiiplikation von Polynom und Skalar durchzufuehren
    
    Args:
        polynom (PolynomringElement): Operand fuer Multiplikation
        skalar             (integer): Operand fuer Multiplikation
    Returns:
        PolynomringElement: Ergebnis der Multiplikation von Polynom und Skalar     

    """
    l = polynom.koeffizienten.laenge
    helpKoeff = [0]*l
    for i in range(0,l):
        helpKoeff[i] = polynom.koeffizienten[i] * skalar
        if(i + 1 == l):
            break
    return PolynomringElement(helpKoeff,polynom.ring)

def createKoeffList2(C):
    """Gibt Liste von Koeffizienten eines gegebenen Polynoms zurueck
    
    Args:
        C (PolynomringElement): Polynom aus welchem die Koeffizienten zu extrahieren sind
    Returns:
        list: Liste der Koeffizienten des gegebene Polynoms
    """
    helpKoeff = [0] * (C.koeffizienten.laenge)
    for i in range(0, C.koeffizienten.laenge):
        if(C.koeffizienten[i] != 0):
            helpKoeff[i] = C.koeffizienten[i]
    return helpKoeff

# 
def potenzieren(C, x, K):
    """ Potenziert Koeffizienten des Polynoms C um x
    
    C = PolynomringElement([1,0,1],K) := 1 + x**2
    >>> potenzieren(C,x,K)
    x = 1 ---> [0,1,0,1]   := x + x**3
    x = 2 ---> [0,0,1,0,1] := x**2 + x**4

    Args: 
        C (PolynomringElement): Zu potenzirenedes Polynom
        x (integer): Erhoehung der Potenz der Koeffizienten des Polynoms um x
    Returns:
        list: Potenzierte Liste von Koeffizienten

    """
    koeffizienten = createKoeffList2(C)
    zeros = [K.null for i in range(x)]
    zeros[len(zeros):] = koeffizienten
    return zeros

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
    B = PolynomringElement([K.eins],P)
    x = 1 # Iterator
    L = 0 # Anzahl angenommener Fehler 
    b = K.eins
    N = 0 # Aktuell betrachtete Elemente aus Sequenz

    
    for N in range(len(seq)):     # Durchlaufe alle Elemente der gegebenen Sequenz

        d = diskrepanz(C,seq,N,K) # Berechnet Diskrepanz der erechneten Sequenz und gegbenen Sequenz
        if (d == K.null):         # Falls die Diskrepanz null betraegt wurde Sequenz korrekt berechnet durch das Polynom C beschrieben
            x = x + 1             # somit kann mit dem naechsten Elemetn fortgefuehrt werden
        else:                     # Falls nicht, muss das Polynom C angepasst werden  
            if 2*L > N: 
                mB = skalarPolyMultiplikation(B,(d/b))
                hilfkoeff = potenzieren(mB, x, K)
                C = C - PolynomringElement(hilfkoeff, P)
                x = x + 1
            else:
                mB = skalarPolyMultiplikation(B,(d/b))
                B = C
                hilfkoeff = potenzieren(mB, x, K)
                C = C - PolynomringElement(hilfkoeff, P) 
                L = N + 1 - L
                b = d
                x = 1
    C = reverseKoeffizienten(C)
    return time.time() - tStart, P.basisring.modulus


















def xor_list(list1, list2):
    temp = list1
    for i in range(0, len(list1)):
        if ((list1[i] == 1 or list2[i] == 1) and list1[i] != list2[i]):
            temp[i] = 1
        else : temp[i] = 0
    return temp

def createKoeffList(f,N):
    list1 = [0] * N 
    for i in range(0, f.koeffizienten.laenge - 1):
        list1[i] = f.koeffizienten[i]
    return list1

def Berlekamp_Massey_algorithm_binary(sequence):
    N = len(sequence)
    s = sequence[:]

    for k in range(N):
        if s[k] == 1:
            break
    var = [0] * (k + 2)
    var[k+1] = 1
    var[0] = 1
    f = PolynomringElement(var,R_X)
    
    l = k + 1 
    g = PolynomringElement([1],R_X)
    a = k 
    b = 0

    for n in range(k + 1, N):
        d = 0
        #print("#Koeffizienten in f:")
        #print(f)
        for i in range(0, f.koeffizienten.laenge):
            if(f.koeffizienten[i] != 0):
                d ^= s[i + n - l]
       
        if d == 0:
            b += 1
        else:
            if 2 * l > n:
                #print("Koeffizienten in g: ")
                #print(g)
                help_koeff = [0] * N
                for i in range(0, g.koeffizienten.laenge):
                    if(g.koeffizienten[i] != 0):
                        help_koeff[a - b + i] = 1
                t = createKoeffList(f, N)

                koeffizienten = xor_list(help_koeff,t)
                f = PolynomringElement(koeffizienten, R_X)
                b += 1
            else:
                temp = f
                
                help_koeff = [0] * N
            
                #print(b-a)
                #print("Koeffizienten von f vor ^ g")
                #print(f)
                #print(g)
                
                for i in range(0, f.koeffizienten.laenge):
                    if(f.koeffizienten[i] != 0 ):
                        help_koeff[b - a + i] = 1
                t1 = createKoeffList(g, N)
                koeffizienten = xor_list(help_koeff, t1)

                f = PolynomringElement(koeffizienten,R_X)
                #print("Koeffizienten von f nach ^ g")
                #print(f)

                l = n + 1 - l
                g = temp
                a = b
                b = n - l + 1

    return (f, l)
