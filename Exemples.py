import numpy as np
import pandas as pd
from Aeronave import *

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
            alpha_e = np.array([9.209, 6.161, 4.072, 2.577, 1.471, 0.63, -0.025, -0.544, -0.963, -1.306, -1.59])
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

        ro = 0.3809     # densidade do ar kg/m^3
        # ro = 1.2754     # densidade do ar kg/m^3 para 20°C
        T = 20          # temperatura do ar °C
        n = 11          # numero de elementos

        # derivadas e velocidade
        
        pd.DataFrame(np.round(r, 4), columns=['V0', 'Lv', 'Lp', 'Lr', 'Le', 'Lc', 'Nv', 'Np', 'Nr', 'Ne', 'Nc', 'Yv', 'Yp', 'Yr', 'Ye', 'Yc'])

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

        # EV=======================================================================================
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

        # EH=======================================================================================
        t = Tail(1, 1, 1, [1, 1])

        # profundor
        el = Elevator(1, 1)

        # fuselagem
        b = Body(
            Sl = 3.88,
            h = 0.88
        )

        a = Aircraft(w, f, t, b, V = V0)
        a.set_mass(
            ro = ro,
            mass = 318,
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
        a.set_control(      # coloca as superfícies de controle
            aileron = ai,
            rudder = ru,
            elevator = el
        )
        # calculos=================================================================================
        a.estimate_Coefs(
            k = 1,
            T = 20,
            ro = ro
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


if __name__ == "__main__":
        pass