from pandas import DataFrame
import numpy as np
import math
import matplotlib.pyplot as plt

from Util import erro_dataframe
from Aeronave import *

# exemplo 13.1 COOK (1 velocidade)
V0 = 25.75          # m/s
alpha_e = 2.577     # deg

ro = 0.3809         # densidade do ar kg/m^3
T = 20              # temperatura do ar °C
n = 11              # numero de elementos

# dados do exemplo
r = np.array([
[25.750, -0.053, -0.402, 0.097, -0.505, 0.012, 0.055, -0.036, -0.024, 0.0057, -0.053, -0.235, 0, 0.068, 0, 0.173]
])
real = DataFrame(np.round(r, 4), columns=['V0', 'Lv', 'Lp', 'Lr', 'Le', 'Lc', 'Nv', 'Np', 'Nr', 'Ne', 'Nc', 'Yv', 'Yp', 'Yr', 'Ye', 'Yc'])
#==================================================================================================
# geometry
# asa
w = Wing(
    S = 12.7,
    b = 15,
    mac = 0.835,
    c12 = [1.04, 0.63]
)
w.set_angles(
    T = 2,
    V_c4 = -0.8,
    V_LE = 0,
    inc = 9
)
# aileron
ai = Aileron(
    S = 0.38*w.b/2 * 0.23,
    c = 0.23,
    y1 = 0.56*w.b/2,
    y2 = 0.94*w.b/2
)
ai.set_CLa(4.136)

y = np.linspace(0, w.b/2, n)
cy = 1.04 - 0.82*y/15

# EV
f = Fin(
    S = 0.96,
    b = 1.25,
    c12 = [0.878, 0.405]
)
f.set_angles(
    V_c4 = 16
)
# leme
ru = Rudder(
    S = 0.42,
    c = 0
)
ru.set_CLa(f.AR, 5.405)

h = np.linspace(0, f.b, n)
ch = f.chord(h)

# EH
t = Tail(1, 1, 1, [1, 1])

# profundor
el = Elevator(1, 1)

# fuselagem
b = Body(
    Sl = 3.88,
    h = 0.88
)

#==================================================================================================
a = Aircraft(w, f, t, b, V = V0)
a.set_mass(
    mass = 318,
    Ix = 1370,
    Iz = 1770,
    Ixz = -4.1
)
a.adim_mass(
    ro = ro
)
a.set_angles(
    alpha = alpha_e,
    theta = alpha_e
)
a.set_fin(
    lf = 4.595 + 5*w.mac/12,
    Lf = 4.595,
    hf = 0.118 + 1.25/2
)
a.set_control(
    aileron = ai,
    rudder = ru,
    elevator = el
)
#==================================================================================================
# calculos
a.estimate_Coefs(
    k = 1,
    T = 20,
    ro = ro
)
CDe = 0.013 + 1.13*a.get_CL_eq()**2/(math.pi*w.AR)

a.derivatives(
    dCL_day = a.w.CLa,
    dCD_day = a.w.CDa,
    dCL_dah = a.f.CLa,
    dCDy_de = a.w.CDa,
    CDy = CDe,
    CLy = a.get_CL_eq(),
    cy = cy,
    ch = ch
)
deri = a.get_derivatives()  # organiza as derivadas em dataframe

# print para comparar as derivadas
print(f"Referência:\n{real}")
print(f"Calculado:\n{deri}")
print(f"Erro:\n{erro_dataframe(real, deri)}")

a.dim_derivatives(ro = ro)  # dimensionaliza derivadas

#==================================================================================================
# DINAMICA
from Dinamica_LateroDirecional import Dinamica_LateroDirecional
din = Dinamica_LateroDirecional(a)

w, c = din.aprox_freq()             # frequências natural e amortecida
tr = np.round(din.aprox_Tr(), 3)    # tempo de rolagem
ts = np.round(din.aprox_Ts(), 3)    # tempo de espiral

print(f"omega_d: {np.round(w, 3)}\nzeta_d: {np.round(c, 3)}\nTr: {tr}s\nTs: {ts}s")

mA, mB = din.A_B()    # matrizes A e B
delta, mN = din.G()  # polinomio caracteristico e matriz N

print(f"Matriz A:\n{mA}\nMatriz B:\n{mB}\ndelta: {delta}")

din.step()      # plot 
din.root_map()  # plot das frequências

#==================================================================================================
# NORMA
from Norma import MILF8587C
norm = MILF8587C(
    Class = 1,
    Category = "B"
)
# calcula o level do voo
l_dutch = norm.dutch_roll(w, c)
l_roll = norm.roll_mode(tr)
l_spiral = norm.spiral_stability(ts)

print(f"\nlevel da norma MIL-F-8587C\nDutch: {l_dutch}\nRoll: {l_roll}\nSpiral: {l_spiral}")

plt.show()
