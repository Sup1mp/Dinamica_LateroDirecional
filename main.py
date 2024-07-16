from Dinamica_LateroDirecional import Dinamica_LateroDirecional
from pandas import DataFrame
import Aeronave as A
import numpy as np
import math
import matplotlib.pyplot as plt

V0 = [round(i*0.5144444, 3) for i in range(35, 90, 5)]
alpha_e = [9.209, 6.161, 4.072, 2.577, 1.471, 0.63, -0.025, -0.544, -0.963, -1.306, -1.59]

ro = 0.3809     # densidade do ar kg/m^3
# ro = 1.2754     # densidade do ar kg/m^3 para 20°C
T = 20          # temperatura do ar °C
n = 11          # numero de elementos

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
a = A.Aircraft(w, f, t, b, V = V0[5])
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
    alpha = alpha_e[5],
    theta = alpha_e[5]
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
    dCL_day = a.w.CLa,
    dCD_day = a.w.CDa,
    dCL_dah = a.f.CLa,
    dCDy_de = a.w.CDa,
    CDy = CDe,
    CLy = a.get_CL_eq(),
    cy = cy,
    ch = ch
)
a.dim_derivatives(
    ro = ro
)

din = Dinamica_LateroDirecional(a)
din.A_B()
din.G()
din.step()

print(f"Tr: {round(din.aprox_Tr(), 3)}s")
print(f"Ts: {round(din.aprox_Ts(), 3)}s")
w, c = din.aprox_freq()
print(f"omega_d: {round(w, 3)}\nzeta_d: {round(c, 3)}")

plt.show()
