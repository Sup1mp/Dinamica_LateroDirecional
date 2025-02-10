from Exemples import Dart_T51_Sailplane
import Util.Util as util
import numpy as np

# Exemplo utilizado
a, real = Dart_T51_Sailplane(modo = 2)

# comparação de erro
util.compara_derivadas(real, a.get_derivatives())

# DINAMICA
from Dinamica_LateroDirecional import Dinamica_LateroDirecional
din = Dinamica_LateroDirecional(a)

w, c = din.aprox_freq()             # frequências natural e amortecida
tr = np.round(din.aprox_Tr(), 3)    # tempo de rolagem
ts = np.round(din.aprox_Ts(), 3)    # tempo de espiral

print(f"omega_d: {np.round(w, 3)}\nzeta_d: {np.round(c, 3)}\nTr: {tr}s\nTs: {ts}s")

# # Norma
# from Norma import MILF8587C
# norm = MILF8587C(
#     Class = 1,
#     Category = "B"
# )
# l_dutch = []
# l_roll = []
# l_spiral = []
# for n in range(len(V0)):
#     l_dutch.append(norm.dutch_roll(w[n], c[n]))
#     l_roll.append(norm.roll_mode(tr[n]))
#     l_spiral.append(norm.spiral_stability(ts[n]))

# print(f"\nlevel da norma MIL-F-8587C\nDutch: {l_dutch}\nRoll: {l_roll}\nSpiral: {l_spiral}")
