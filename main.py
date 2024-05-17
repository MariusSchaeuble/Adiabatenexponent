import matplotlib
import numpy
import numpy as np
import sympy as sym
from Helpers import identifier, isCharacter
import math
from numpy import matrix, array, mean, std, max, linspace, ones, sin, cos, tan, arctan, pi, sqrt, exp, arcsin, arccos, arctan2, sinh, cosh
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, title, savefig, errorbar, grid
import scipy.optimize as opt
from GPII import *
from math import sqrt
pi = math.pi

matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)

def gauss(term):
    ids = identifier(term)
    symbols = []
    for str1 in ids:
        symbols.append(sym.sympify(str1))
    termSymbol = sym.sympify(term)
    values = []
    for ident in ids:
        exec("values.append(" + ident + ")")

    derivatives = []
    i = 0
    while i < len(symbols):
        r = sym.diff(termSymbol, symbols[i])
        j = 0
        while j < len(symbols):
            # exec('r.evalf(subs={symbols[j]: ' + values[j] + '})')
            r = r.evalf(subs={symbols[j]: values[j]})
            j += 1
        derivatives.append(r.evalf())
        i += 1
    i = 0
    while i < len(derivatives):
        exec("derivatives[i] *= sigma_" + ids[i])
        i = i + 1
    res = 0
    for z in derivatives:
        res += z ** 2
    return math.sqrt(res)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def linear(x, a, b):
    return a*x+b



#1. teil
W1 = matrix("""
18.6    4.1;
16.3    3.6;
14.7    3.6;
23.1    5;
19.5    4.2;
19.2    4.2;
14.5    3.4;
16.6    3.4;
18.5    3.7 ;
12.2    2.4;
12      2.4;
16.9    3.6;
13.5    2.9;
24.5    5.2;
22.6    4.9;
19.6    4.1;
18      4.2;
25.8    5.5;
23.7    5.1;
16.6    3.6;
22.2    5;
21.9    5.1;
14.7    3.4;
23.1    5;
20.1    4.4;
22.9    5.3;
19.2    4.4;
22      5.3;
21.9    4.8;
21      5;
22.6    4.8
""")# h1, h2


adiabat1 = ones(len(W1))
for i in range(len(W1)):
    adiabat1[i] = W1[i, 0]/(W1[i, 0] - W1[i, 1])




#2. teil
#L: länge volumen in cm
L_offset=6
m_Schwingkoerper=11.67/1000
sigma_m_Schwingkoerper = 0.01/1000
#d: durchmesser schwingkörper
d=1.65/100
sigma_d = 0.25/1000

#atmosphärendrucl literaturwert

#CO2
W_CO2 = matrix("""
0    6044;
3    7324;
6    8468;
9    9456;
12   10329;
18   11478;
22   12665;
28   14045
""")# eingestellt Länge Volumenkörper in cm , zeit für 30 schwingungen in ms
#L_offset+L_matrix=L_ges


V = ones(8)
Tq = ones(8)
sigma_Tq = ones(len(Tq))
for i in range(len(V)):
    V[i] = W_CO2[i, 0]/100*pi*(10.3/200)**2

    t = W_CO2[i, 1]
    sigma_t = 500
    Tq[i] = (t/30000)**2
    sigma_Tq[i] = gauss("(t/30000)**2")


errorbar(V, Tq, sigma_Tq, None, 'x', label="Measured Time, CO2", linewidth=3)
optimizedParameters, s = opt.curve_fit(linear, V, Tq)
plt.plot(V, linear(V, *optimizedParameters), label="fit")
a_std, b_std = np.sqrt(np.diag(s))

#plt.title(r"$\alpha$", fontsize=25)
plt.xlabel(r"Volumina in m^3", fontsize=25)
plt.ylabel(r"T^2 in s^2", fontsize=25)
legend(fontsize=15)
grid()
plt.tight_layout()
savefig('co2.png')
show()

m = optimizedParameters[0]
sigma_m = a_std

p = 1005*100
sigma_p = 10*100

gamma_CO2 = 64*m_Schwingkoerper/(m*p*d**4)
sigma_gamma_CO2 = gauss("64*m_Schwingkoerper/(m*p*d**4)")



W_Argon = matrix("""
0    5393;
4    7032;
8    8356;
12   9456;
16   10115;
18   10818;
24   12329;
28   13577
""")

V = ones(8)
Tq = ones(8)
sigma_Tq = ones(len(Tq))
for i in range(len(V)):
    V[i] = W_Argon[i, 0]/100*pi*(10.3/200)**2

    t = W_Argon[i, 1]
    sigma_t = 500
    Tq[i] = (t/30000)**2
    sigma_Tq[i] = gauss("(t/30000)**2")


errorbar(V, Tq, sigma_Tq, None, 'x', label="Measured Time, Argon", linewidth=3)
optimizedParameters, s = opt.curve_fit(linear, V, Tq)
plt.plot(V, linear(V, *optimizedParameters), label="fit")
a_std, b_std = np.sqrt(np.diag(s))

#plt.title(r"$\alpha$", fontsize=25)
plt.xlabel(r"Volumina in m^3", fontsize=25)
plt.ylabel(r"T^2 in s^2", fontsize=25)
legend(fontsize=15)
grid()
plt.tight_layout()
savefig('Argon.png')
show()

m = optimizedParameters[0]
sigma_m = a_std

p = 1005*100
sigma_p = 10*100

gamma_Argon = 64*m_Schwingkoerper/(m*p*d**4)
sigma_gamma_Argon = gauss("64*m_Schwingkoerper/(m*p*d**4)")

W_Stick = matrix("""
0     5874;
4     7515;
8     8823;
12    9894;
16    10520;

24    12866;
28    13650
""")

# 20    112026; unrealistisch

V = ones(7)
Tq = ones(7)
sigma_Tq = ones(len(Tq))
for i in range(len(V)):
    V[i] = W_Stick[i, 0]/100*pi*(10.3/200)**2

    t = W_Stick[i, 1]
    sigma_t = 500
    Tq[i] = (t/30000)**2
    sigma_Tq[i] = gauss("(t/30000)**2")


errorbar(V, Tq, sigma_Tq, None, 'x', label="Measured Time, Stickstoff", linewidth=3)
optimizedParameters, s = opt.curve_fit(linear, V, Tq)
plt.plot(V, linear(V, *optimizedParameters), label="fit")
a_std, b_std = np.sqrt(np.diag(s))

#plt.title(r"$\alpha$", fontsize=25)
plt.xlabel(r"Volumina in m^3", fontsize=25)
plt.ylabel(r"T^2 in s^2", fontsize=25)
legend(fontsize=15)
grid()
plt.tight_layout()
savefig('Stick.png')
show()

m = optimizedParameters[0]
sigma_m = a_std

p = 1005*100
sigma_p = 10*100

gamma_Stick = 64*m_Schwingkoerper/(m*p*d**4)
sigma_gamma_Stick = gauss("64*m_Schwingkoerper/(m*p*d**4)")


#Messwerte drucken
latexTable(UC(W1[:, 0], 'u'), UC(W1[:, 1], 'u'))
latexTable(UC(W_CO2[:, 0], 'cm'), UC(W_CO2[:, 1], 'ms'))
latexTable(UC(W_Argon[:, 0], 'cm'), UC(W_Argon[:, 1], 'ms'))
latexTable(UC(W_Stick[:, 0], 'cm'), UC(W_Stick[:, 1], 'ms'))