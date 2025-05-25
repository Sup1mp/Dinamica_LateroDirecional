import Exemples as ex
import Util.util as util
import matplotlib.pyplot as plt
import numpy as np

from Dinamica_LateroDirecional import Dinamica_LateroDirecional
from Norma import MILF8587C

modo = 1

# Exemplo utilizado
a, real = ex.Dart_T51_Sailplane(modo = modo)

# comparação de erro
util.compara_derivadas(real, a.get_derivatives())

# DINAMICA

din = Dinamica_LateroDirecional(a)

w, c = din.aprox_freq()             # frequências natural e amortecida
tr = np.round(din.aprox_Tr(), 3)    # tempo de rolagem
ts = np.round(din.aprox_Ts(), 3)    # tempo de espiral

print(f"omega_d: {np.round(w, 3)}\nzeta_d: {np.round(c, 3)}\nTr: {tr}s\nTs: {ts}s")

# Norma

norm = MILF8587C(
    Class = 1,
    Category = "B"
)
if modo == 1:
    l_dutch = norm.dutch_roll(w, c)
    l_roll = norm.roll_mode(tr)
    l_spiral = norm.spiral_stability(ts)

elif modo == 2:
    l_dutch = []
    l_roll = []
    l_spiral = []
    for n in range(len(w)):
        l_dutch.append(norm.dutch_roll(w[n], c[n]))
        l_roll.append(norm.roll_mode(tr[n]))
        l_spiral.append(norm.spiral_stability(ts[n]))

print(f"\nlevel da norma MIL-F-8587C\nDutch: {l_dutch}\nRoll: {l_roll}\nSpiral: {l_spiral}")

plt.show()
