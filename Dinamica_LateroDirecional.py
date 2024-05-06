import numpy as np
from scipy import signal
from math import sin, cos, radians, sqrt
import matplotlib.pyplot as plt

from Aeronave import Aircraft, Finn, Wing

class Dinamica_LateroDirecional:
    def __init__ (self, Aviao : Aircraft):
        '''
        Classe da Dinâmica Latero-Direcional de aeronaves
        '''
        self.a = Aviao
        self.w = Aviao.w
        self.f = Aviao.f

        return

    def matriz_A (self, Yv:float, Lv:float, Nv:float, Yp:float, Lp:float, Np:float, Yr:float, Lr:float, Nr:float):
        '''
        Calcula a matriz de estabilidade A (COOK et. al, 2013) utilizando as derivadas adimensionais
        '''

        a1 = self.a.Ix*self.a.Iz - (self.a.Ixz**2)   # simplificação
        g = 9.81    # gravidade

        # voo reto e simétrico
        Ue = self.a.V0*cos(self.a.theta_e)
        We = self.a.V0*sin(self.a.theta_e)

        # matriz de estabilidade
        self.A = np.array([
            [Yv/self.a.m, (Yp*self.w.b + self.a.m*We)/self.a.m, (Yr*self.w.b - self.a.m*Ue)/self.a.m, g*cos(self.a.theta_e), g*sin(self.a.theta_e)],
            [(self.a.Ixz*Nv + self.a.Iz*Lv)/a1, (self.a.Ixz*Np + self.a.Iz*Lp)*self.w.b/a1, (self.a.Ixz*Nr + self.a.Iz*Lr)*self.w.b/a1, 0, 0],
            [(self.a.Ix*Nv + self.a.Ixz*Lv)/a1, (self.a.Ix*Np + self.a.Ixz*Lp)*self.w.b/a1, (self.a.Ix*Nr + self.a.Ixz*Lr)*self.w.b/a1, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0]
        ])

        return self.A

    def matriz_B (self, Ye:float, Le:float, Ne:float, Yc:float, Lc:float, Nc:float):
        '''
        Calcula a matriz de controle B (COOK et. al, 2013) utilizando as derivadas adimensionais
        '''

        a1 = self.a.Ix*self.a.Iz - (self.a.Ixz**2)   # simplificação

        # matriz de controle
        self.B = self.a.V0 * np.array([
            [Ye/self.a.m, Yc/self.a.m],
            [(self.a.Ixz*Ne + self.a.Iz*Le)/a1, (self.a.Ixz*Nc + self.a.Iz*Lc)/a1],
            [(self.a.Ix*Ne + self.a.Ixz*Le)/a1, (self.a.Ix*Nc + self.a.Ixz*Lc)/a1],
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

        # polinomio caracteristico 
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

        # matrix N
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
            Ix_d : momento de inércia a.Ix eixo aeronautico (kg*m^2)
            Lp_d : derivada dimensional Lp
        '''
        self.Tr_ap = -Ix_d/Lp_d
        
        return self.Tr_ap
    
    def aprox_Ts (self, V0, Lv, Nv, Lp, Np, Lr, Nr):
        '''
        Retorna a constante de tempo Ts do Spiral Mode aproximada
            a.V0 : velocidade da aeronave (m/s)
        '''
        self.Ts_ap = -V0/9.81 * (Lv*Np - Lp*Nv)/(Lr*Nv - Lv*Nr)

        return self.Ts_ap
    
    def aprox_freq (self, Iz_d, m_d, Yv_d, Nv_d, Nr_d):
        '''
        Retorna as frequências naturais omega_d e de amortecimento zeta_d aproximadas para o Dutch Roll
            Iz_d : momento de inércia a.Iz eixo aeronautico (kg*m^2)
            m_d : massa da aeronave (kg)
            Yv_d : derivada dimensional Yv
            Nv_d : derivada dimensional Nv
            Nr_d : derivada dimensional Nr
        '''

        self.wd_ap = sqrt(self.a.V0*Nv_d/Iz_d)                         # frequencia natural omega_d
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
    import main