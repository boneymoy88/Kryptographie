from Tocas import Polynomring, PolynomringElement, Z
from Extension import endlicher_koerper
from Extension.PolynomRestklassenring import PolynomRestklassenring, PolynomRestklassenringElement

                        
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
def disc(poly,seq,n,polyRing,field):
    total = field.null
    #print(n)
    #print(poly.koeffizienten[0] * seq[n])
    for i in range(poly.grad + 1):
        d = poly.koeffizienten[i] * seq[n-i]
        #print(type(field.eins))
        #print(type(d))
        #print(type(total))
        total = total + d
    #print(total)
    return total

# Hilfsfunktion um Multiiplikation von Polynom und Skalar durchzufuehren
def skalarPolyMultiplikation(F,d,Ring):
    helpkoeff = [0]*F.koeffizienten.laenge
    #print(F.koeffizienten.laenge)
    for i in range(0,F.koeffizienten.laenge):
        #print("koeff")
        #print(F.koeffizienten[i] * d)
        helpkoeff[i] = F.koeffizienten[i] * d
        if(i + 1 == F.koeffizienten.laenge):
            break
    return PolynomringElement(helpkoeff,Ring)

# Gibt Liste von Koeffizienten eines gegebenen Polynoms zurueck
def createKoeffList2(C,N):
    list1 = [0] * (N+1)
    #print(N)
    #print(f.koeffizienten[0])
    for i in range(0, C.koeffizienten.laenge):
        #print(f.koeffizienten[i])
        if(C.koeffizienten[i] != 0):
            #print(f.koeffizienten[i])
            list1[i] = C.koeffizienten[i]
    return list1

# Potenziert Koeffizienten des Polynoms C um x
def potenzieren(C, N, x, field):
    koeffizienten = createKoeffList2(C,N)
    zeros = [field.null for i in range(x)]
    zeros[len(zeros):] = koeffizienten
    return zeros


# Berlekamp Massey Algorithmus findet kleinstes Polynom welches gegebene Sequenz
# aus endlichem Koerper berechnet
# seq: Sequenz von Elementen aus
# field: endlicher Koerper 
def scalarMassey(seq,field):
    # test if field is endlicher_koerper
    polyRing = Polynomring(field)
    C = PolynomringElement([field.eins],polyRing)
    B = PolynomringElement([field.eins],polyRing)
    x = 1
    L = 0
    b = field.eins
    N = 0
    # Durchlaufe alle Elemente der gegebenen Sequenz
    for N in range(len(seq)):
        # Berechnet Diskrepanz der erechneten Sequenz und gegbenen Sequenz
        d = disc(C,seq,N,polyRing,field)
        # Falls die Diskrepanz null betraegt wurde Sequenz korrekt berechnet durch das Polynom C beschrieben
        # somit kann mit dem naechsten Elemetn fortgefuehrt werden
        if (d == field.null):
            x = x+1
        # Falls nicht, muss das Polynom C angepasst werden  
        else:
            if 2*L > N:
                mB = skalarPolyMultiplikation(B,(d/b),polyRing)# B * (d/b)
                hilfkoeff = potenzieren(mB, N, x, field)
                C = C - PolynomringElement(hilfkoeff, polyRing)
                x = x + 1
            else:
                mB = skalarPolyMultiplikation(B,(d/b),polyRing)
                B = C
                hilfkoeff = potenzieren(mB, N, x, field)
                # Potenzieren in dem 0 an anfang der Liste eingesetzt wird
                C = C - PolynomringElement(hilfkoeff, polyRing) 
                L = N + 1 - L
                b = d
                x = 1
    #C = PolynomringElement( ,polyRing)
    return C

