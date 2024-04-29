import numpy as np
from scipy import signal
from math import sin, cos, radians, sqrt
import matplotlib.pyplot as plt

class Dinamica_LateroDirecional:
    def __init__(self, m:float, Ix:float, Iz:float, Ixz:float, ro:float, V0:float, S:float, b:float, theta_e:float):
        '''
        Classe da Dinâmica Latero-Direcional de aeronaves, se inicializa com:
            m : massa (kg)
            Ix : momento de inércia Ix eixo aeronautico (kg*m^2)
            Iz : momento de inércia Iz eixo aeronautico(kg*m^2)
            Ixz : momento de inércia Ixz eixo aeronautico(kg*m^2)
            ro : densidade do ar (kg/m^3)
            V0 : velocidade da aeronave (m/s)
            S : área da asa (m^2)
            b : envergadura da asa (m)
            theta_e : ângulo de equilíbrio entre horizonte e direção de voo (deg)
        '''

        self.S = S                          # área da asa (m^2)
        self.b = b                          # envergadura da asa (m)
        self.V0 = V0                        # velocidade (m/s)
        self.theta_e = radians(theta_e)     # angulo de equilibrio entre direção de voo e o horizonte (deg)

        # adimensionalizadores
        ad1 = 0.5*ro*V0*S
        ad2 = ad1*b

        self.m = m/ad1                      # massa adimensional
        self.Ix = Ix/ad2                    # momento de inercia Ix adimensional
        self.Iz = Iz/ad2                    # momento de inercia Iz adimensional
        self.Ixz = Ixz/ad2                  # momento de inercia Ixz adimensional

        return

    def matriz_A (self, Yv:float, Lv:float, Nv:float, Yp:float, Lp:float, Np:float, Yr:float, Lr:float, Nr:float):
        '''
        Calcula a matriz de estabilidade A (COOK et. al, 2013) utilizando as derivadas adimensionais
        '''

        a1 = self.Ix*self.Iz - (self.Ixz**2)   # simplificação
        g = 9.81    # gravidade

        # voo reto e simétrico
        Ue = self.V0*cos(self.theta_e)
        We = self.V0*sin(self.theta_e)

        # matriz de estabilidade
        self.A = np.array([
            [Yv/self.m, (Yp*self.b + self.m*We)/self.m, (Yr*self.b - self.m*Ue)/self.m, g*cos(self.theta_e), g*sin(self.theta_e)],
            [(self.Ixz*Nv + self.Iz*Lv)/a1, (self.Ixz*Np + self.Iz*Lp)*self.b/a1, (self.Ixz*Nr + self.Iz*Lr)*self.b/a1, 0, 0],
            [(self.Ix*Nv + self.Ixz*Lv)/a1, (self.Ix*Np + self.Ixz*Lp)*self.b/a1, (self.Ix*Nr + self.Ixz*Lr)*self.b/a1, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0]
        ])

        return self.A

    def matriz_B (self, Ye:float, Le:float, Ne:float, Yc:float, Lc:float, Nc:float):
        '''
        Calcula a matriz de controle B (COOK et. al, 2013) utilizando as derivadas adimensionais
        '''

        a1 = self.Ix*self.Iz - (self.Ixz**2)   # simplificação

        # matriz de controle
        self.B = self.V0 * np.array([
            [Ye/self.m, Yc/self.m],
            [(self.Ixz*Ne + self.Iz*Le)/a1, (self.Ixz*Nc + self.Iz*Lc)/a1],
            [(self.Ix*Ne + self.Ixz*Le)/a1, (self.Ix*Nc + self.Ixz*Lc)/a1],
            [0, 0],
            [0, 0]
        ])

        return self.B
    
    def matriz_G (self):

        yv = self.A[0,0]
        yp = self.A[0,1]
        yr = self.A[0,2]
        yphi = self.A[0,3]
        ypsi = self.A[0,4]

        lv = self.A[1,0]
        lp = self.A[1,1]
        lr = self.A[1,2]

        nv = self.A[2,0]
        npp = self.A[2,1]
        nr = self.A[2,2]

        ye = self.B[0,0]
        yc = self.B[0,1]

        le = self.B[1,0]
        lc = self.B[1,1]

        ne = self.B[2,0]
        nc = self.B[2,1]

        self.delta = np.array([
            1,
            -(lp + nr + yv),
            (lp*nr - lr*npp) + (nr*yv - nv*yr) + (lp*yv - lv*yp),
            lv*(nr*yp - npp*yr - yphi) + nv*(lp*yr - lr*yp - ypsi) + yv*(lr*npp - lp*nr),
            lv*(nr*yphi - npp*ypsi) + nv*(lp*ypsi - lr*yphi),
            0
        ])

        Nv_e = np.array([
            ye,
            le*yp + ne*yr - ye*(lp + nr),
            le*(npp*yr - nr*yp + yphi) + ne*(lr*yp - lp*yr + ypsi) + ye*(lp*nr - lr*npp),
            le*(npp*ypsi - nr*yphi) + ne*(lr*yphi - lp*ypsi),
            0
        ])

        Nv_c = np.array([
            yc,
            lc*yp + nc*yr - yc*(lp + nr),
            lc*(npp*yr - nr*yp + yphi) + nc*(lr*yp - lp*yr + ypsi) + yc*(lp*nr - lr*npp),
            lc*(npp*ypsi - nr*yphi) + nc*(lr*yphi - lp*ypsi),
            0
        ])

        Np_e = np.array([
            le,
            -le*(nr + yv) + ne*lr + ye*lv,
            le*(nr*yv - nv*yr) + ne*(lv*yr - lr*yv) + ye*(lr*nv - lv*nr),
            -le*nv*ypsi + ne*lv*ypsi,
            0
        ])

        Np_c = np.array([
            lc,
            -lc*(nr + yv) + nc*lr + yc*lv,
            lc*(nr*yv - nv*yr) + nc*(lv*yr - lr*yv) + yc*(lr*nv - lv*nr),
            -lc*nv*ypsi + nc*lv*ypsi,
            0
        ])

        Nphi_e = np.array([
            le,
            -le*(nr + yv) + ne*lr + yc*lv,
            le*(nr*yv - nv*yr) + ne*(lv*yr - lr*yv) + yc*(lr*nv - lv*nr),
            -le*nv*ypsi + ne*lv*ypsi
        ])

        Nphi_c = np.array([
            lc,
            -lc*(nr + yv) + nc*lr + yc*lv,
            lc*(nr*yv - nv*yr) + nc*(lv*yr - lr*yv) + yc*(lr*nv - lv*nr),
            -lc*nv*ypsi + nc*lv*ypsi
        ])

        Nr_e = np.array([
            ne,
            le*npp - ne*(lp + yv) + ye*nv,
            le*(nv*yp - npp*yv) + ne*(lp*yv - lv*yp) + ye*(lv*npp - lp*nv),
            le*nv*yphi - ne*lv*yphi,
            0
        ])

        Nr_c = np.array([
            nc,
            lc*npp - nc*(lp + yv) + yc*nv,
            lc*(nv*yp - npp*yv) + nc*(lp*yv - lv*yp) + yc*(lv*npp - lp*nv),
            lc*nv*yphi - nc*lv*yphi,
            0
        ])

        Npsi_e = np.array([
            ne,
            le*npp - ne*(lp + yv) + ye*nv,
            le*(nv*yp - npp*yv) + ne*(lp*yv - lv*yp) + ye*(lv*npp - lp*nv),
            le*nv*yphi - ne*lv*yphi
        ])

        Npsi_c = np.array([
            nc,
            lc*npp - nc*(lp + yv) + yc*nv,
            lc*(nv*yp - npp*yv) + nc*(lp*yv - lv*yp) + yc*(lv*npp - lp*nv),
            lc*nv*yphi - nc*lv*yphi
        ])

        self.N = [
            Nv_e, Nv_c,
            Np_e, Np_c,
            Nr_e, Nr_c,
            Nphi_e, Nphi_c,
            Npsi_e, Npsi_c
        ]


        return self.delta, self.N
    
    def aprox_Tr (self, Ix_d, Lp_d):
        '''
        Retorna a constante de tempo Tr do Roll Mode aproximada
            Ix_d : momento de inércia Ix eixo aeronautico (kg*m^2)
            Lp_d : derivada dimensional Lp
        '''
        self.Tr_ap = -Ix_d/Lp_d
        
        return self.Tr_ap
    
    def aprox_Ts (self, V0, Lv, Nv, Lp, Np, Lr, Nr):
        '''
        Retorna a constante de tempo Ts do Spiral Mode aproximada
            V0 : velocidade da aeronave (m/s)
        '''
        self.Ts_ap = -V0/9.81 * (Lv*Np - Lp*Nv)/(Lr*Nv - Lv*Nr)

        return self.Ts_ap
    
    def aprox_freq (self, Iz_d, m_d, Yv_d, Nv_d, Nr_d):
        '''
        Retorna as frequências naturais omega_d e de amortecimento zeta_d aproximadas para o Dutch Roll
            Iz_d : momento de inércia Iz eixo aeronautico (kg*m^2)
            m_d : massa da aeronave (kg)
            Yv_d : derivada dimensional Yv
            Nv_d : derivada dimensional Nv
            Nr_d : derivada dimensional Nr
        '''

        self.wd_ap = sqrt(V0*Nv_d/Iz_d)                         # frequencia natural omega_d
        self.cd_ap = -(Nr_d/Iz_d + Yv_d/m_d)/(2*self.wd_ap)     # frequencia de amortecimento zeta_d

        return self.wd_ap, self.cd_ap
    
    def step(self):
        '''
        Apresenta o gráfico para a resposta em step do sistema
        '''
        
        fig, ax = plt.subplots(5, 2)

        titles = ['$v_{\epsilon}$', '$v_{\zeta}$',
                  '$r_{\epsilon}$', '$r_{\zeta}$',
                  '$p_{\epsilon}$', '$p_{\zeta}$',
                  '$\phi_{\epsilon}$', '$\phi_{\zeta}$',
                  '$\psi_{\epsilon}$', '$\psi_{\zeta}$'
        ]

        for i in range(len(self.N)):

            t, y = signal.step(signal.lti(self.N[i], self.delta))   # step response

            ax[i//2][i%2].plot(t, y, 'k')
            ax[i//2][i%2].set(xlabel = 't(s)', ylabel=titles[i])
            ax[i//2][i%2].grid()
            
        fig.tight_layout()

        return

if __name__ == "__main__":
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

    ld = Dinamica_LateroDirecional(m, Ix, Iz, Ixz, ro, V0, S, b, theta_e)

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