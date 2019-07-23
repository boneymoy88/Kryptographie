from Tocas import Polynomring, PolynomringElement, Z
from Extension import endlicher_koerper 

                        
R_X = Polynomring(Z)
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
        print("#Koeffizienten:")
        for i in range(0, f.koeffizienten.laenge):
            if(f.koeffizienten[i] == 1):
                print(i)
                d ^= s[i + n - l]
       
        if d == 0:
            b += 1
        else:
            if 2 * l > n:
                print("in g ")
                print([i for i in range(0, g.koeffizienten.laenge)])
                koeff = [a - b + i for i in range(0, g.koeffizienten.laenge)]
                f += PolynomringElement(koeff, R_X)
                b += 1
            else:
                temp = f
                print([i for i in range(0, f.koeffizienten.laenge)])
                koeff = [b - a + i for i in range(0, f.koeffizienten.laenge)]
                f = PolynomringElement(koeff,R_X) + g
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
    
    g = set([0])
    a = k
    b = 0
    for n in range(k + 1, N):
        d = 0
        print("#Koeffizienten:")
        for ele in f:
            print(ele)
            d ^= s[ele + n - l]

        if d == 0:
            b += 1
        else:
            if 2 * l > n:
                print("in g ")
                print([ele for ele in g])
                f ^= set([a - b + ele for ele in g])
                b += 1
            else:
                temp = f.copy()
                print([ele for ele in f])
                f = set([b - a + ele for ele in f]) ^ g
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
