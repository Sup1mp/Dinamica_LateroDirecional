from Dinamica_LateroDirecional import Dinamica_LateroDirecional
from Aeronave import Aircraft, Finn, Wing
import matplotlib.pyplot as plt

# Derivadas do McDonnell F-4C Phantom, body axes reference
ro = 0.3809     # densidade do ar kg/m^3

m = 17642       # kg
Ix = 33898      # kgm^2
Iz = 189496     # kgm^2
Ixz = 2952      # kgm^2

S = 49.239      # area da asa m^2
b = 11.787      # envergadura da asa m

V0 = 178        # velocidade m/s
theta_e = 9.4   # graus

asa = Wing(S, b, 1, 1, [1, 1])
ev = Finn(1, 1, 1, 1, [1,1], 1, 1, 1)

aviao = Aircraft(m, Ix, Ixz, Iz, V0, theta_e)
aviao.ad_mass(ro, asa.S, asa.b)

ld = Dinamica_LateroDirecional(aviao, asa, ev)

# derivadas
Y = {
    'v': -0.5974,
    'p': 0,
    'r': 0,
    'e': -0.0159,
    'c': 0.1193
}
L = {
    'v': -0.1048,
    'p': -0.1164,
    'r': 0.0455,
    'e': 0.0454,
    'c': 0.0086
}
N = {
    'v': 0.0987,
    'p': -0.0045,
    'r': -0.1132,
    'e': 0.00084,
    'c': -0.0741
}

A = ld.matriz_A (Y['v'], L['v'], N['v'], Y['p'], L['p'], N['p'], Y['r'], L['r'], N['r'])
B = ld.matriz_B (Y['e'], L['e'], N['e'], Y['c'], L['c'], N['c'])

# print(f"{A}\n\n{B}")

delta, N = ld.matriz_G()

# print(f'\n{delta}\n\n{N}')

ld.step()

plt.show()