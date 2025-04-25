import numpy as np
import pandas as pd
from Util.unidades import ft2m
from Util.util import mach, compara_derivadas, getAtmosphere
from Aeronave import *

def Boeing_747_100 (modo: int = 1):
     
    if modo == 1:
        V0 = 157.886        # m/s
        alt = ft2m(20000)   # m
        T, P, ro, s = getAtmosphere(alt)

    elif modo == 2:
        V0 = [67.3608, 157.886, 265.481]    # m/s
        alt = [0, ft2m(20000), ft2m(40000)] # m
        T, P, ro, s = getAtmosphere(alt)
    else:
        raise ValueError("modo desconhecido, colocar somente 1 ou 2")
    
    n = 3
    # real = pd.DataFrame(np.round(r, 4), columns=['V0', 'Lv', 'Lp', 'Lr', 'Le', 'Lc', 'Nv', 'Np', 'Nr', 'Ne', 'Nc', 'Yv', 'Yp', 'Yr', 'Ye', 'Yc'])
    #=============================================================================================
    w = Wing(
        S = 510.9667,             # m^2
        b = 59.643264,            # m
        mac = 8.324088,           # m
        c12 = [48.995, 13.3],     # m
        th = 0.113
    )
    w.set_angles(
        T = 8,        # deg
        V_c4 = 51.5,  # deg
        V_LE = 48,    # deg
        inc = 2.5     # deg
    )
    y1 = np.linspace(3.403092, 12.7635, n//2)
    cy1 = w.c12[0] + (9.460992 - w.c12[0])*y1/12.7635
    y2 = np.linspace(12.7635, w.b/2, n//2+1 if n%2 != 0 else n/2)
    cy2 = 9.460992 + (w.c12[1] - 9.460992)*y2/(w.b/2)
    cy = np.append(cy1, cy2)

    al = Aileron(
        S = 7.8837,       # m^2
        c = 1.2192,       # m
        y1 = 21.56,       # m
        y2 = 28.0263      # m
    )
    #=============================================================================================
    f = Fin(
        S = 74.539,       # m
        b = 9.7139,       # m
        c12 = [11.7165, 3.6302],  # m
        th = 0.113    # aproxima como mesmo perfil da asa
    )
    f.set_angles(
        V_c4 = 45,    # deg
        V_LE = 40     # deg
    )

    h = np.linspace(0, f.b, n)
    ch = f.chord(h)
    
    ru = Rudder(
        S = 22.4522,    # m^2
        c = (1.26492 + 3.959352)/2
    )
    #=============================================================================================
    t = Tail(1, 1, 1)
    el = Elevator(1, 1)
    #=============================================================================================
    b = Body(
        Sl = 305.0581,    # m^2
        h = 6.43128       # m
    )
    #=============================================================================================
    a = Aircraft(w, f, t, b, V0)
    a.set_control(al, el, ru)
    a.set_angles(
        alpha = 0,
        theta = 2
    )

    a.set_fin(
        lf = 30.323282,
        Lf = 29.629608,
        hf = 8.330184
    )
    a.set_mass(
        ro = ro,
        mass = 288756.903,
        Ix = 24675886.69,
        Iz = 67384152.12,
        Ixz = 1315143.412
    )

    a.estimate_Coefs(
        k = 1,
        ro = ro,
        M = mach(V0, alt)
    )

    a.derivatives(
        dCL_day = a.w.CLa if modo == 1 else np.array([a.w.CLa for _ in range (len(V0))]).transpose(),
        dCD_day = a.w.CDa if modo == 1 else np.array([a.w.CDa for _ in range (len(V0))]).transpose(),
        dCL_dah = a.f.CLa if modo == 1 else np.array([a.f.CLa for _ in range (len(V0))]).transpose(),
        dCDy_de = a.f.CDa if modo == 1 else np.array([a.w.CDa for _ in range (len(V0))]).transpose(),
        CDy = 0.04,
        CLy = a.get_CL_eq(),
        cy = cy,
        ch = ch
    )

    return a

def Dart_T51_Sailplane(modo: int = 1):
        '''
        Dados do exemplo 13.1 COOK:\n
        modo = 1 : Apenas 1 velocidade\n
        modo = 2 : Todas as velocidades
        '''
        if modo == 1:
            # apenas dados para 1 velocidade
            V0 = 25.75          # m/s
            alpha_e = 2.577     # deg
            r = np.array([
            [25.750, -0.053, -0.402, 0.097, -0.505, 0.012, 0.055, -0.036, -0.024, 0.0057, -0.053, -0.235, 0, 0.068, 0, 0.173]
            ])

        elif modo == 2:
            # todos os dados dado no exemplo
            V0 = np.array([round(i*0.5144444, 3) for i in range(35, 90, 5)])
            alpha_e = [9.209, 6.161, 4.072, 2.577, 1.471, 0.63, -0.025, -0.544, -0.963, -1.306, -1.59]
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
        else:
            raise ValueError("modo desconhecido, colocar somente 1 ou 2")
        
        alt = 20000 # m
        # alt = np.linspace(0, 10000, 11)
        T, P, ro, S = getAtmosphere(alt)
        # ro = 0.3809     # densidade do ar kg/m^3
        # ro = 1.2754     # densidade do ar kg/m^3 para 20°C
        n = 11          # numero de elementos

        # derivadas e velocidade
        
        real = pd.DataFrame(np.round(r, 4), columns=['V0', 'Lv', 'Lp', 'Lr', 'Le', 'Lc', 'Nv', 'Np', 'Nr', 'Ne', 'Nc', 'Yv', 'Yp', 'Yr', 'Ye', 'Yc'])

        # asa
        w = Wing(
            S = 12.7,           # m^2
            b = 15,             # m
            mac = 0.835,        # m
            c12 = [1.04, 0.63],
            th = 0.18
        )

        w.set_CL(
              CL0 = 0.3875, # 1/rad
              CLa = 5.55    # 1/rad
        )
        w.set_CD(
             CD0= 0.013,
             CDa = 1.13
        )
        w.set_angles(
            T = 2,          # deg
            V_c4 = -0.8,    # deg
            V_LE = 0,       # deg
            inc = 9         # deg
        )

        # aileron
        ai = Aileron(
            S = 0.38*w.b/2 * 0.23,  # m^2
            c = 0.23,               # m
            y1 = 0.56*w.b/2,        # m
            y2 = 0.94*w.b/2         # m
        )
        ai.set_CLa(
              CLa = 4.136           # 1/rad
        )

        y = np.linspace(0, w.b/2, n)
        cy = (w.c12[1] - w.c12[0])*y/(w.b/2) + w.c12[0]

        # EV=======================================================================================
        f = Fin(
            S = 0.96,               # m^2
            b = 1.25,               # m
            c12 = [0.878, 0.405],
            th = 0.15
        )
        f.set_CL(
              CL0 = 0,
              CLa = 3.68
        )

        f.set_angles(
            V_c4 = 16               # deg
        )
        # leme
        ru = Rudder(
            S = 0.42,               # m^2
            c = 0.336               # m
        )
        ru.set_CLa(
              CLa = 3.68            # 1/rad
        )

        h = np.linspace(0, f.b, n)
        ch = f.chord(h)

        # EH=======================================================================================
        t = Tail(
              S = 1.14,             # m^2
              b = 2.61,             # m
              mac = 0.439,          # m
              c12 = [0.574, 0.304],
              th = 0.15
        )

        # profundor
        el = Elevator(1, 1)

        # fuselagem
        b = Body(
            Sl = 3.88,
            h = 0.88
        )

        a = Aircraft(w, f, t, b, V0 = V0)
        a.set_control(      # coloca as superfícies de controle
            aileron = ai,
            rudder = ru,
            elevator = el
        )
        a.set_mass(
            ro = ro,        # kg/m^3
            mass = 318,     # kg
            Ix = 1370,
            Iz = 1770,
            Ixz = -4.1
        )
        a.set_angles(
            alpha = alpha_e,
            theta = alpha_e
        )
        a.set_fin(          # define as dimensões e posiçao do leme
            lf = 4.595 + 5*w.mac/12,
            Lf = 4.595,
            hf = 0.118 + 1.25/2
        )
        # calculos=================================================================================
        a.estimate_Coefs(
            k = 1,
            ro = ro,
            M = mach(V0, alt)
        )

        CDe = 0.013 + 1.13*a.get_CL_eq()**2/(math.pi*w.AR)

        a.derivatives(
            dCL_day = a.w.CLa if modo == 1 else np.array([a.w.CLa for _ in range (len(V0))]).transpose(),
            dCD_day = a.w.CDa if modo == 1 else np.array([a.w.CDa for _ in range (len(V0))]).transpose(),
            dCL_dah = a.f.CLa if modo == 1 else np.array([a.f.CLa for _ in range (len(V0))]).transpose(),
            dCDy_de = a.f.CDa if modo == 1 else np.array([a.w.CDa for _ in range (len(V0))]).transpose(),
            CDy = CDe,
            CLy = a.get_CL_eq(),
            cy = cy,
            ch = ch
        )

        return a, real

def DR_1_AeroBat (modo: int = 1):
     w = Wing(
          S = 118
     )
     return

#=============================================================================== 
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    # EXEMPLO DART
    a, real = Dart_T51_Sailplane(modo = 1)
    # a = Boeing_747_100(modo=1)

    # comparação de erro
    # compara_derivadas(real, a.get_derivatives())

    # DINAMICA
    from Dinamica_LateroDirecional import Dinamica_LateroDirecional
    din = Dinamica_LateroDirecional(a)

    w, c = din.aprox_freq()             # frequências natural e amortecida
    tr = np.round(din.aprox_Tr(), 3)    # tempo de rolagem
    ts = np.round(din.aprox_Ts(), 3)    # tempo de espiral

    print(f"\nomega_d: {np.round(w, 3)}\nzeta_d: {np.round(c, 3)}\nTr: {tr} s\nTs: {ts} s")
    
    A, B = din.A_B()        # calculo matrizes A e B
    delta, N = din.G()      # calculo matriz G e N

    # Matrizes
    # print(f"A: {np.round(A, 3)}")
    # print(f"B: {np.round(B, 3)}")
    # print(f"N: {np.round(N, 3)}")

    # print(f"Delta: {np.round(delta, 3)}")

    din.step()
    din.root_map()

    ro = ft2m(getAtmosphere(20000)[2])
    # Lf, Sf = din.get_Lf_Sf(ro, w, c)
    # print(f"Lf: {Lf}\nSf: {Sf}")
    Sf = din.get_Sf_for(ro, w, c)
    print(f"\nSf real: {a.f.S}\nSf calc: {Sf}")

    # # Norma
    # from Norma import MILF8587C
    # norm = MILF8587C(
    #     Class = 1,
    #     Category = "B"
    # )
    # l_dutch = []
    # l_roll = []
    # l_spiral = []
    # for n in range(len(w)):
    #     l_dutch.append(norm.dutch_roll(w[n], c[n]))
    #     l_roll.append(norm.roll_mode(tr[n]))
    #     l_spiral.append(norm.spiral_stability(ts[n]))

    # print(f"\nlevel da norma MIL-F-8587C\nDutch: {l_dutch}\nRoll: {l_roll}\nSpiral: {l_spiral}")
    
    plt.show()