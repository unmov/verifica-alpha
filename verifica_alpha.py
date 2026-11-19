"""Verifica della localizzazione di alpha - Misure di Hausdorff e applicazioni (PoliTO)"""
from fractions import Fraction as F
from math import comb
from itertools import product

def add(z, w): return (z[0] + w[0], z[1] + w[1])
def sub(z, w): return (z[0] - w[0], z[1] - w[1])
def mul(z, w): return (z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0])
def scale(a, z): return (a*z[0], a*z[1])
def norm2(z): return z[0]**2 + z[1]**2        
def upper(z): return abs(z[0]) + abs(z[1])        # maggiorante di |z|
def lower(z): return max(abs(z[0]), abs(z[1]))    # minorante di |z|

def power(z, k):
    r = (F(1), F(0))
    for _ in range(k):
        r = mul(r, z)
    return r

def evaluate(coeff, z):
    r = (F(0), F(0))
    for k, a in enumerate(coeff):
        r = add(r, scale(a, power(z, k)))
    return r

# P(z) = 2z^3 - (1 + z + ... + z^12)
a = [F(-1)] * 13
a[3] = F(1)
c = (F(-418, 1000), F(551, 1000))
eps = F(1, 1000)

# Coefficienti di P(c+h) = sum_k A_k h^k (formula di Taylor).
A = [evaluate([a[j] * comb(j, k) for j in range(k, 13)], c) for k in range(13)]

# Ipotesi del teorema di Rouche' su |h| = eps, con f(h) = A_1 h e g = P(c+h) - f.
resto = upper(A[0]) + sum(upper(A[k]) * eps**k for k in range(2, 13))
assert resto < lower(A[1]) * eps
print("P ha esattamente una radice nel disco B(c, 1/1000).")
