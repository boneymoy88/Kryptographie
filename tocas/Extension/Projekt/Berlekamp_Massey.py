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

def Berlekamp_Massey_algorithm_old(sequence):

    N = len(sequence)
    s = sequence[:]

    for k in range(N):
        if s[k] == 1:
            break
    f = set([k + 1, 0])  # use a set to denote polynomial
    l = k + 1
    g = set([0]) 
            # output the polynomial
    def print_poly(polynomial):
        result = ''
        lis = sorted(polynomial, reverse=True)
        for i in lis:
            if i == 0:
                result += '1'
            else:
                result += 'x^%s' % str(i)

            if i != lis[-1]:
                result += ' + '

        return result

    a = k
    b = 0
    for n in range(k + 1, N):
        d = 0
        #print("#Koeffizienten in f: ")
        #print(print_poly(f))
        for ele in f:
            d ^= s[ele + n - l]

        if d == 0:
            b += 1
        else:
            if 2 * l > n:
                #print("in g: ")
                #print(print_poly(g))
                f ^= set([a - b + ele for ele in g])
                b += 1
            else:
                #print(b-a)
                temp = f.copy()
                #print([b - a + ele for ele in f])
                #print("Koeffizienten in f vor ^g") 
                #print(print_poly(f))
                #print(print_poly(g))
                f = set([b - a + ele for ele in f])
                #print(print_poly(f))
                f ^= g
                #print("Koeffizienten in f nach ^g")
                #print(print_poly(set([ele for ele in f])))
                l = n + 1 - l
                g = temp
                a = b
                b = n - l + 1

    return (print_poly(f), l)

def berlekamp_massey(s,K):

    
    if (not K.ist_endlicher_koerper()):
        print("Kein endlicher Koerper!")
        return 0
    else:
        print("Ist endlichr Koerper!")

    print(K)
    C = PolynomRestklassenringElement([1],K)
    B = PolynomRestklassenringElement([1],K)
    n = 0
    L = 0
    m = 1
    b = 1 #field K?
    N = len(s)
    
    for n in range(0, N-1):
        helpSum = 0
        for k in range(0,L):
            print("helpSum")
            print(k)
            helpSum += C.wert.koeffizienten[k]*s[n-k]
        d = s[n] + helpSum
        print(d)
        if (d == 1):
            m = m + 1
        elif (2 * L <= n):
            T = C
            print(C)
            C = C - d * b**(-1) * K.erzeuger**m * B
            print(C)
            L = n + 1 - L
            B = T
            b = d
            m = 1
        else:
            if(b != 0):
                C = C - d * b**(-1) * K.erzeuger**m * B
            else:
                C = C
            m = m + 1
    return (C, L)

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

def createKoeffList2(f,N):
    list1 = [0] * (N+1)
    #print(N)
    #print(f.koeffizienten[0])
    for i in range(0, f.koeffizienten.laenge):
        #print(f.koeffizienten[i])
        if(f.koeffizienten[i] != 0):
            #print(f.koeffizienten[i])
            list1[i] = f.koeffizienten[i]
    return list1

def scalarMassey(seq,field):
    #print(type(field.eins))
    polyRing = Polynomring(field)
    C = PolynomringElement([field.eins],polyRing)
    B = PolynomringElement([field.eins],polyRing)
    x = 1
    L = 0
    b = field.eins
    N = 0
    for N in range(len(seq)):
        d = disc(C,seq,N,polyRing,field)
        #print(N)
        if (d == field.null):
            x = x+1
        else:
            if 2*L > N:
                mB = skalarPolyMultiplikation(B,(d/b),polyRing)# B * (d/b)
                print("test")
                #print(mB)
                #print(C)
                hilfkoeff = createKoeffList2(mB,N)   
                print(hilfkoeff)
                zeros=[field.null for i in range(x)]
                zeros[len(zeros):] = hilfkoeff
                hilfkoeff = zeros
                #hilfkoeff = hilfkoeff[:N+1]
                print(hilfkoeff)
                C = C - PolynomringElement(hilfkoeff, polyRing)
                #polyRing.sub(C,polyRing.timesXPower(mB,x))
                #print(C)
                x = x + 1
            else:
                print("mod")
                mB = skalarPolyMultiplikation(B,(d/b),polyRing)
                #print(mB)
                B = C
                hilfkoeff = createKoeffList2(mB,N)   
                print(hilfkoeff)
                print(x)
                zeros = [field.null for i in range(x)]
                print(zeros)
                # Potenzieren in dem 0 an anfang der Liste eingesetzt wird
                zeros[len(zeros):] = hilfkoeff
                hilfkoeff = zeros
                #hilfkoeff = hilfkoeff[:N+1]
                print(hilfkoeff)
                C = C - PolynomringElement(hilfkoeff, polyRing) #polyRing.sub(C,polyRing.timesXPower(mB,x))
                #print(C)
                L = N + 1 - L
                b = d
                x = 1
    #C = PolynomringElement( ,polyRing)
    return C

    
#input binary sequence s
#output polynom
#def berlekamp_massey(s):
#    N = len(s)
#    # erste auftretende 1 in k merken
#    for k in range(N):
#        if s[k] == 1:
#            break
#    
#    C = PolynomringElement([1], R_X)
#    T = PolynomringElement([1], R_X)
#    L = 0
#    m = -1
#    B = PolynomringElement([1], R_X)
#    N = 0
#    #calculate discrepancy
#    while(n < N):
#        d = s[n] + sumOfSequence(s, L, m, N)
#        if (d == 1):
#            T = C
#            C = C + B * R_X.variable**(N-m)
#            if (L <= (N/2)):
#                L = N + 1 - L
#                m = N
#                B = T
#        else:
#            N += 1
#    return C

#def sumOfSequence(s, L, m, N):
#    discrepancy = 0
#    for i in range(0, m - 1):
#        discrepancy += C.koeffizienten[i] * s[N - 1 - i]
#    return discrepancy
