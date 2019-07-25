# -*- coding=utf-8 -*-

import project_path
from Tocas import Polynomring, PolynomringElement, Z
from Extension.Projekt import Berlekamp_Massey as LFSR

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
R_X = Polynomring(Z)

a = 2
b = 2
f = PolynomringElement([1,0,0,1], R_X)
g = PolynomringElement([1], R_X)
if (f.koeffizienten.laenge > g.koeffizienten.laenge):
    help_koeff = [0] * f.koeffizienten.laenge
else: help_koeff = [0] * g.koeffizienten.laenge

for i in range(0, g.koeffizienten.laenge):
    help_koeff[a - b + i] = 1

def xor_list(list1, list2):
    temp = list1
    for i in range(0, len(list1)):
        if ((list1[i] == 1 or list2[i] == 1) and list1[i] != list2[i]):
            temp[i] = 1
        else : temp[i] = 0
    return temp

def createKoeffList(f):
    list1 = [0] * f.koeffizienten.laenge 
    for i in range(0, f.koeffizienten.laenge):
        list1[i] = f.koeffizienten[i]
    return list1

s = 0,1,1,1,0,1,1,0,1,1,0
print(LFSR.Berlekamp_Massey_algorithm(s)) 
print(LFSR.Berlekamp_Massey_algorithm_tocas(s))


