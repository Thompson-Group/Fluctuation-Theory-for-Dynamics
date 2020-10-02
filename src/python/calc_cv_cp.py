import numpy as np


p = 1.0 # bar
Na = 6.0221409e23
conv = 2.39005736e-29 * Na # barA^3 -> kcal/mol
conv = 1.43932e-5
conv2 = 1000 # kcal/mol -> cal/mol
R = 1.9872036e-3 # Gas Constant (kcal/mol/K)
R = R*1000 # -> cal/mol/K
T = 298.15 # K
N = 343 # Molecules
div = (R*T*T)
E_file = "e_init.out"
V_file = "vol_init.out"

E = np.genfromtxt(E_file, usecols=0)
V = np.genfromtxt(V_file, usecols=0)

E = [x for x in E]

E_av = np.average(E)
V_av = np.average(V)

print("<E>: %s" % E_av)
print("<V>: %s" % V_av)

dE = [x - E_av for x in E]
dV = [x - V_av for x in V]

dE_av = np.average(dE)
dV_av = np.average(dV)

print("<dE>: %s" % dE_av)
print("<dV>: %s" % dV_av)

dE2 = [x**2.0 for x in dE]
dV2 = [x**2.0 for x in dV]
dEdV = [x*y for x,y in zip(dE,dV)]

dE2_av = np.average(dE2)
dV2_av = np.average(dV2)
dEdV_av = np.average(dEdV)
print("<dE2>: %s" % dE2_av)
print("<dV2>: %s" % dV2_av)
print("<dEdV>: %s" % dEdV_av)

term1 = dE2_av*conv2**2.0                              # (cal/mol/K)^2
term2 = (p**2.0)*dV2_av*(conv**2.0)*conv2**2.0         # bar^2A^6 -> (kcal/mol/K)^2 ->(cal/mol/K)^2
term3 = 2*p*dEdV_av*conv*conv2**2.0                    # barA3(kcal/mol/K) -> (kcal/mol/K)^2->(cal/mol/K)^2
dH2_av = (term1+term2+term3)
print("<dH2>: %s" % dH2_av)
print("terms: %s %s %s" % (term1, term2, term3))
print("terms: %s %s %s" % (term1/div, term2/div, term3/div))
print("R (cal/mol/K): %s" % R)

Cv = term1/div/N
Cp = dH2_av/div/N
print("Cv (cal mol^-1 K^-1): %s" % Cv)
print("Cp (cal mol^-1 K^-1): %s" % Cp)
