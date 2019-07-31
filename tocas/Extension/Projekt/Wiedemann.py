import numpy as np
from Tocas import Polynomring, PolynomringElement, Z
from Extension import endlicher_koerper
#Wiedemann algorithm 

R_X = Polynomring(Z)

def wiedemann(A, b):
    #nonsingular square case (n x n)
    listOfMatrices = []
    n = len(A)
    for i in range(0,2*n - 1):
        listOfMatrices.append(A^i * b)

    k = 0
    g0 = PolynomringElement([1], R_X)
    
    
    return 0
#compute A^i * b for i = 0,1,...,2n-1

