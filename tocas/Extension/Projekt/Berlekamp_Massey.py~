from Tocas import Polynomring, PolynomringElement, Z
from Extension import endlicher_koerper
from Extension.PolynomRestklassenring import PolynomRestklassenring, PolynomRestklassenringElement
import Extension.Primzahl as Primzahl
import time

                        
R_X = Polynomring(Z)

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

# Berechnet Diskrepanz der erechneten Sequenz und gegbenen Sequenz
def diskrepanz(polynom,seq,N,K):
    total = K.null
    for i in range(polynom.grad + 1):
        d = polynom.koeffizienten[i] * seq[N-i]
        total = total + d
    return total

# Hilfsfunktion um Multiiplikation von Polynom und Skalar durchzufuehrenf
def skalarPolyMultiplikation(polynom,skalar,Ring):
    l = polynom.koeffizienten.laenge
    helpKoeff = [0]*l
    for i in range(0,l):
        helpKoeff[i] = polynom.koeffizienten[i] * skalar
        if(i + 1 == l):
            break
    return PolynomringElement(helpKoeff,Ring)

# Gibt Liste von Koeffizienten eines gegebenen Polynoms zurueck
def createKoeffList2(C,N):
    helpKoeff = [0] * (N+1)
    for i in range(0, C.koeffizienten.laenge):
        if(C.koeffizienten[i] != 0):
            helpKoeff[i] = C.koeffizienten[i]
    return helpKoeff

# Potenziert Koeffizienten des Polynoms C um x
def potenzieren(C, N, x, K):
    koeffizienten = createKoeffList2(C,N)
    zeros = [K.null for i in range(x)]
    zeros[len(zeros):] = koeffizienten
    return zeros


# Berlekamp Massey Algorithmus findet kleinstes Polynom welches gegebene Sequenz
# aus endlichem Koerper berechnet
#
# seq: Sequenz von Elementen aus
# K: endlicher Koerper 
def scalarMassey(seq,K):

    #
    # Laufzeitberechnung intialisieren
    #
    tStart = time.time()
    #
    # Testen ob K endlicher_koerper 
    #
    if (not K.ist_endlicher_koerper()):
        raise RuntimeError("Berlekamp Massey nur auf endlichen Koerpern.")
    #
    # Initialisierungen von
    # Polynomring: polyring
    # Polynomring Elemente: C,B
    # 
    #
    P = Polynomring(K)

    if (not Primzahl.miller_rabin(P.basisring.modulus)):
        raise RuntimeError("Polynomring nicht ueber Primzahl")

    C = PolynomringElement([K.eins],P)
    B = PolynomringElement([K.eins],P)
    x = 1
    L = 0
    b = K.eins
    N = 0

    #
    # Durchlaufe alle Elemente der gegebenen Sequenz
    #
    for N in range(len(seq)):
        # Berechnet Diskrepanz der erechneten Sequenz und gegbenen Sequenz
        d = diskrepanz(C,seq,N,K)
        # Falls die Diskrepanz null betraegt wurde Sequenz korrekt berechnet durch das Polynom C beschrieben
        # somit kann mit dem naechsten Elemetn fortgefuehrt werden
        if (d == K.null):
            x = x + 1
        # Falls nicht, muss das Polynom C angepasst werden  
        else:
            if 2*L > N:
                mB = skalarPolyMultiplikation(B,(d/b),P)# B * (d/b)
                hilfkoeff = potenzieren(mB, N, x, K)
                C = C - PolynomringElement(hilfkoeff, P)
                x = x + 1
            else:
                mB = skalarPolyMultiplikation(B,(d/b),P)
                B = C
                hilfkoeff = potenzieren(mB, N, x, K)
                # Potenzieren in dem 0 an anfang der Liste eingesetzt wird
                C = C - PolynomringElement(hilfkoeff, P) 
                L = N + 1 - L
                b = d
                x = 1
    C = reverseKoeffizienten(C)
    #print("Laufzeit in ms: ", time.time() - tStart)
    return time.time() - tStart, P.basisring.modulus

def reverseKoeffizienten(C):
    anzahlKoeffizienten = C.koeffizienten.laenge
    reverseList = [0]*anzahlKoeffizienten
    for i in range(anzahlKoeffizienten):
        reverseList[anzahlKoeffizienten - 1 - i] = C.koeffizienten[i]
    return PolynomringElement(reverseList, C.ring)
        
