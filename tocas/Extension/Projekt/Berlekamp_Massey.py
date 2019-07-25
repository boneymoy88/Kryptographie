from Tocas import Polynomring, PolynomringElement, Z
from Extension import endlicher_koerper 

                        
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
    for i in range(0, f.koeffizienten.laenge):
        list1[i] = f.koeffizienten[i]
    return list1



def Berlekamp_Massey_algorithm_tocas(sequence):
    N = len(sequence)
    s = sequence[:]

    for k in range(N):
        if s[k] == 1:
            break
    var = [0] * (k + 2)
    var[k+1] = 1
    var[0] = 1
    f = PolynomringElement(var,R_X)
    
    l = k + 1 # l=3
    g = PolynomringElement([1],R_X)
    a = k #a=2
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
                t = createKoeffList(f, N)
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

def Berlekamp_Massey_algorithm(sequence):

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
