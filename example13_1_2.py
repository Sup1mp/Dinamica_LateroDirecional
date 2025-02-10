from pandas import DataFrame
import Aeronave as A
import numpy as np
import math
import matplotlib.pyplot as plt
import Util.Util as Util

# exemplo 13.1 COOK (todas as velocidades)
V0 = np.array([round(i*0.5144444, 3) for i in range(35, 90, 5)])
alpha_e = np.array([9.209, 6.161, 4.072, 2.577, 1.471, 0.63, -0.025, -0.544, -0.963, -1.306, -1.59])

ro = 0.3809     # densidade do ar kg/m^3
# ro = 1.2754     # densidade do ar kg/m^3 para 20°C
T = 20          # temperatura do ar °C
n = 11          # numero de elementos

# dados do exemplo
r = np.array([
[18.025, -0.042, -0.403, 0.187, -0.505, 0.006, 0.056, -0.074, -0.028, 0.0120, -0.054, -0.235, 0, 0.070, 0, 0.173],
[20.600, -0.047, -0.402, 0.145, -0.505, 0.009, 0.056, -0.055, -0.026, 0.0089, -0.054, -0.235, 0, 0.069, 0, 0.173],
[23.175, -0.050, -0.402, 0.117, -0.505, 0.011, 0.056, -0.044, -0.025, 0.0070, -0.054, -0.235, 0, 0.069, 0, 0.173],
[25.750, -0.053, -0.402, 0.097, -0.505, 0.012, 0.055, -0.036, -0.024, 0.0057, -0.053, -0.235, 0, 0.068, 0, 0.173],
[28.325, -0.054, -0.401, 0.082, -0.505, 0.013, 0.055, -0.030, -0.023, 0.0047, -0.053, -0.235, 0, 0.068, 0, 0.173],
[30.900, -0.056, -0.401, 0.070, -0.505, 0.014, 0.055, -0.025, -0.023, 0.0039, -0.053, -0.235, 0, 0.068, 0, 0.173],
[33.475, -0.057, -0.401, 0.061, -0.505, 0.015, 0.055, -0.022, -0.023, 0.0035, -0.053, -0.235, 0, 0.068, 0, 0.173],
[36.050, -0.058, -0.401, 0.054, -0.505, 0.015, 0.055, -0.019, -0.023, 0.0029, -0.052, -0.235, 0, 0.067, 0, 0.173],
[38.625, -0.058, -0.401, 0.048, -0.505, 0.016, 0.054, -0.017, -0.022, 0.0025, -0.052, -0.235, 0, 0.067, 0, 0.173],
[41.200, -0.059, -0.401, 0.044, -0.505, 0.016, 0.054, -0.015, -0.022, 0.0022, -0.052, -0.235, 0, 0.067, 0, 0.173],
[43.775, -0.059, -0.401, 0.040, -0.505, 0.016, 0.054, -0.013, -0.022, 0.0020, -0.052, -0.235, 0, 0.067, 0, 0.173]
])
real = DataFrame(np.round(r, 4), columns=['V0', 'Lv', 'Lp', 'Lr', 'Le', 'Lc', 'Nv', 'Np', 'Nr', 'Ne', 'Nc', 'Yv', 'Yp', 'Yr', 'Ye', 'Yc'])

#====================================================================
# geometry

# asa
w = A.Wing(
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
ai = A.Aileron(
    S = 0.38*w.b/2 * 0.23,
    c = 0.23,
    y1 = 0.56*w.b/2,
    y2 = 0.94*w.b/2
)
ai.set_CLa(4.136)

y = np.linspace(0, w.b/2, n)
cy = 1.04 - 0.82*y/15

# EV
f = A.Fin(
    S = 0.96,
    b = 1.25,
    c12 = [0.878, 0.405]
)
f.set_angles(
    V_c4 = 16
)
# leme
ru = A.Rudder(
    S = 0.42,
    c = 0
)
ru.set_CLa(f.AR, 5.405)

h = np.linspace(0, f.b, n)
ch = f.chord(h)

# EH
t = A.Tail(1, 1, 1, [1, 1])

# profundor
el = A.Elevator(1, 1)

# fuselagem
b = A.Body(
    Sl = 3.88,
    h = 0.88
)

#====================================================================
a = A.Aircraft(w, f, t, b, V = V0)
a.set_mass(
    mass = 318,
    Ix = 1370,
    Iz = 1770,
    Ixz = -4.1
)
# a.adim_mass(
#     ro = ro
# )
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
#====================================================================
# calculos
a.estimate_Coefs(
    k = 1,
    T = 20,
    ro = ro
)
CDe = 0.013 + 1.13*a.get_CL_eq()**2/(math.pi*w.AR)

a.derivatives(
    dCL_day = np.array([a.w.CLa for _ in range (len(V0))]).transpose(),
    dCD_day = np.array([a.w.CDa for _ in range (len(V0))]).transpose(),
    dCL_dah = np.array([a.f.CLa for _ in range (len(V0))]).transpose(),
    dCDy_de = np.array([a.w.CDa for _ in range (len(V0))]).transpose(),
    CDy = CDe,
    CLy = a.get_CL_eq(),
    cy = cy,
    ch = ch
)
deri = a.get_derivatives()
err = Util.erro_dataframe(real, deri)
# print(f"Referência:\n{real}")
print(f"Calculado:\n{deri}")
print(f"Erro:\n{err}")
sts = np.round(err.describe(), 4)
print(f"Stats:\n{sts.loc['mean']}\nMean sum: {sts.loc['mean'].sum()}")

# DINAMICA
from Dinamica_LateroDirecional import Dinamica_LateroDirecional
din = Dinamica_LateroDirecional(a)

w, c = din.aprox_freq()             # frequências natural e amortecida
tr = np.round(din.aprox_Tr(), 3)    # tempo de rolagem
ts = np.round(din.aprox_Ts(), 3)    # tempo de espiral

print(f"omega_d: {np.round(w, 3)}\nzeta_d: {np.round(c, 3)}\nTr: {tr}s\nTs: {ts}s")

# Norma
from Norma import MILF8587C
norm = MILF8587C(
    Class = 1,
    Category = "B"
)
l_dutch = []
l_roll = []
l_spiral = []
for n in range(len(V0)):
    l_dutch.append(norm.dutch_roll(w[n], c[n]))
    l_roll.append(norm.roll_mode(tr[n]))
    l_spiral.append(norm.spiral_stability(ts[n]))

print(f"\nlevel da norma MIL-F-8587C\nDutch: {l_dutch}\nRoll: {l_roll}\nSpiral: {l_spiral}")
