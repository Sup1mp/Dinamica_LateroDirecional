import numpy as np
from scipy import signal
import control as ct
import matplotlib.pyplot as plt

from Aeronave import Aircraft

class Dinamica_LateroDirecional:
    def __init__ (self, aero: Aircraft):
        '''
        Classe da Dinâmica Latero-Direcional de aeronaves que utiliza as equações e métodos presentes no COOK et. al (2013)
        '''
        self.aero = aero

        self.z = np.zeros((self.aero._len_velocities)) if self.aero._many_velocities else 0  # vetor de 0
        self.o = np.ones((self.aero._len_velocities))  if self.aero._many_velocities else 1  # vetor de 1

        return

    def A_B (self):
        '''
        Calcula as matrizes "A" e "B" utilizando as derivadas adimensionais do apêndice 8 do mesmo livro 
        '''
        a1 = self.aero.Ix1*self.aero.Iz1 - (self.aero.Ixz1**2)   # pra facilitar
        g = 9.81    # gravidade

        # voo reto e simétrico permite essa simplificação
        Ue = self.aero.V0*np.cos(self.aero.theta)
        We = self.aero.V0*np.sin(self.aero.theta)

        # matriz de estabilidade
        self.A = np.array([
            [self.aero.Yv/self.aero.m1, (self.aero.Yp*self.aero.w.b + self.aero.m1*We)/self.aero.m1, (self.aero.Yr*self.aero.w.b - self.aero.m1*Ue)/self.aero.m1, g*np.cos(self.aero.theta), g*np.sin(self.aero.theta)],
            [(self.aero.Ixz1*self.aero.Nv + self.aero.Iz1*self.aero.Lv)/a1, (self.aero.Ixz1*self.aero.Np + self.aero.Iz1*self.aero.Lp)*self.aero.w.b/a1, (self.aero.Ixz1*self.aero.Nr + self.aero.Iz1*self.aero.Lr)*self.aero.w.b/a1, self.z, self.z],
            [(self.aero.Ix1*self.aero.Nv + self.aero.Ixz1*self.aero.Lv)/a1, (self.aero.Ix1*self.aero.Np + self.aero.Ixz1*self.aero.Lp)*self.aero.w.b/a1, (self.aero.Ix1*self.aero.Nr + self.aero.Ixz1*self.aero.Lr)*self.aero.w.b/a1, self.z, self.z],
            [self.z, self.o, self.z, self.z, self.z],
            [self.z, self.z, self.o, self.z, self.z]
        ])

        # matriz de controle
        self.B = self.aero.V0 * np.array([
            [self.aero.Ye/self.aero.m1, self.aero.Yc/self.aero.m1],
            [(self.aero.Ixz1*self.aero.Ne + self.aero.Iz1*self.aero.Le)/a1, (self.aero.Ixz1*self.aero.Nc + self.aero.Iz1*self.aero.Lc)/a1],
            [(self.aero.Ix1*self.aero.Ne + self.aero.Ixz1*self.aero.Le)/a1, (self.aero.Ix1*self.aero.Nc + self.aero.Ixz1*self.aero.Lc)/a1],
            [self.z, self.z],
            [self.z, self.z]
        ])

        return self.A, self.B
    
    def G (self):
        '''
        Retorna a equação caracteristica (delta) e a matriz de polinomios (N) para cada variavel utilizando
        as equações 
        '''
        if self.aero._many_velocities:
            yv = self.A[0,0,:]
            yp = self.A[0,1,:]
            yr = self.A[0,2,:]
            yphi = self.A[0,3,:]
            ypsi = self.A[0,4,:]

            lv = self.A[1,0,:]
            lp = self.A[1,1,:]
            lr = self.A[1,2,:]

            nv = self.A[2,0,:]
            npp = self.A[2,1,:]
            nr = self.A[2,2,:]

            ye = self.B[0,0,:]
            yc = self.B[0,1,:]

            le = self.B[1,0,:]
            lc = self.B[1,1,:]

            ne = self.B[2,0,:]
            nc = self.B[2,1,:]
        else:
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
            self.o,
            -(lp + nr + yv),
            (lp*nr - lr*npp) + (nr*yv - nv*yr) + (lp*yv - lv*yp),
            lv*(nr*yp - npp*yr - yphi) + nv*(lp*yr - lr*yp - ypsi) + yv*(lr*npp - lp*nr),
            lv*(nr*yphi - npp*ypsi) + nv*(lp*ypsi - lr*yphi),
            self.z
        ]).transpose()*180/np.pi    # rad => deg

        Nv_e = np.array([
            ye,
            le*yp + ne*yr - ye*(lp + nr),
            le*(npp*yr - nr*yp + yphi) + ne*(lr*yp - lp*yr + ypsi) + ye*(lp*nr - lr*npp),
            le*(npp*ypsi - nr*yphi) + ne*(lr*yphi - lp*ypsi),
            self.z
        ])

        Nv_c = np.array([
            yc,
            lc*yp + nc*yr - yc*(lp + nr),
            lc*(npp*yr - nr*yp + yphi) + nc*(lr*yp - lp*yr + ypsi) + yc*(lp*nr - lr*npp),
            lc*(npp*ypsi - nr*yphi) + nc*(lr*yphi - lp*ypsi),
            self.z
        ])

        Np_e = np.array([
            le,
            -le*(nr + yv) + ne*lr + ye*lv,
            le*(nr*yv - nv*yr) + ne*(lv*yr - lr*yv) + ye*(lr*nv - lv*nr),
            -le*nv*ypsi + ne*lv*ypsi,
            self.z
        ])*180/np.pi    # rad/s => deg/s

        Np_c = np.array([
            lc,
            -lc*(nr + yv) + nc*lr + yc*lv,
            lc*(nr*yv - nv*yr) + nc*(lv*yr - lr*yv) + yc*(lr*nv - lv*nr),
            -lc*nv*ypsi + nc*lv*ypsi,
            self.z
        ])*180/np.pi    # rad/s => deg/s

        Nphi_e = np.array([
            self.z,
            le,
            -le*(nr + yv) + ne*lr + yc*lv,
            le*(nr*yv - nv*yr) + ne*(lv*yr - lr*yv) + yc*(lr*nv - lv*nr),
            -le*nv*ypsi + ne*lv*ypsi
        ])*180/np.pi    # rad/s => deg/s

        Nphi_c = np.array([
            self.z,
            lc,
            -lc*(nr + yv) + nc*lr + yc*lv,
            lc*(nr*yv - nv*yr) + nc*(lv*yr - lr*yv) + yc*(lr*nv - lv*nr),
            -lc*nv*ypsi + nc*lv*ypsi
        ])*180/np.pi    # rad/s => deg/s

        Nr_e = np.array([
            ne,
            le*npp - ne*(lp + yv) + ye*nv,
            le*(nv*yp - npp*yv) + ne*(lp*yv - lv*yp) + ye*(lv*npp - lp*nv),
            le*nv*yphi - ne*lv*yphi,
            self.z
        ])*180/np.pi    # rad => deg

        Nr_c = np.array([
            nc,
            lc*npp - nc*(lp + yv) + yc*nv,
            lc*(nv*yp - npp*yv) + nc*(lp*yv - lv*yp) + yc*(lv*npp - lp*nv),
            lc*nv*yphi - nc*lv*yphi,
            self.z
        ])*180/np.pi    # rad => deg

        Npsi_e = np.array([
            self.z,
            ne,
            le*npp - ne*(lp + yv) + ye*nv,
            le*(nv*yp - npp*yv) + ne*(lp*yv - lv*yp) + ye*(lv*npp - lp*nv),
            le*nv*yphi - ne*lv*yphi
        ])*180/np.pi    # rad => deg

        Npsi_c = np.array([
            self.z,
            nc,
            lc*npp - nc*(lp + yv) + yc*nv,
            lc*(nv*yp - npp*yv) + nc*(lp*yv - lv*yp) + yc*(lv*npp - lp*nv),
            lc*nv*yphi - nc*lv*yphi
        ])*180/np.pi    # rad => deg

        # matrix N
        self.N = [
            Nv_e.transpose(), Nv_c.transpose(),
            Np_e.transpose(), Np_c.transpose(),
            Nr_e.transpose(), Nr_c.transpose(),
            Nphi_e.transpose(), Nphi_c.transpose(),
            Npsi_e.transpose(), Npsi_c.transpose()
        ]

        return self.delta, self.N

    # def g (self):
    #     import sympy as sym
    #     from sympy import pprint
    #     from sympy.abc import s

    #     def gg (A, B):
    #         d = s*sym.eye(5) - A
    #         f = sym.Poly(sym.simplify(d.det()))
            
    #         pprint(f)
    #         delt = sym.Poly(d.det(), domain='ZZ')
    #         n = d.adjoint()*sym.Matrix(B)
    #         return delt, n
        
    #     shape = np.shape(self.A)
    #     if len(shape) == 3:
    #         for i in range(shape[2]):
    #             d, n = gg(self.A[:, :, i], self.B[:, :, i])
    #     else:
    #         self.delta, self.N = gg(self.A, self.B)

    #     return np.round(self.delta, 4), self.N

    def aprox_Tr (self):
        '''
        Retorna a constante de tempo Tr do Roll Mode aproximada
        '''
        return -self.aero.Ix/self.aero.Lp1
    
    def aprox_Ts (self):
        '''
        Retorna a constante de tempo Ts do Spiral Mode aproximada
        '''
        return -self.aero.V0/9.81 * (self.aero.Lv*self.aero.Np - self.aero.Lp*self.aero.Nv)/(self.aero.Lr*self.aero.Nv - self.aero.Lv*self.aero.Nr)
    
    def aprox_freq (self):
        '''
        Retorna as frequências naturais omega_d e de amortecimento zeta_d aproximadas para o Dutch Roll
        '''
        wd_ap = np.sqrt(self.aero.V0*self.aero.Nv1/self.aero.Iz)                            # frequencia natural omega_d
        cd_ap = -(self.aero.Nr1/self.aero.Iz + self.aero.Yv1/self.aero.m)/(2*wd_ap)     # frequencia de amortecimento zeta_d

        return wd_ap, cd_ap
    
    def step(self, tmax = 10):
        '''
        Apresenta o gráfico para a resposta em step do sistema\n
            tmax : tempo máximo analisado (s)
        '''
        fig, ax = plt.subplots(5, 2)

        # organização dos titulos na ordem em que aparecem
        titles = ['$v_{\epsilon}$ (m/s)', '$v_{\zeta}$ (m/s)',
                  '$r_{\epsilon}$ (deg/s)', '$r_{\zeta}$ (deg/s)',
                  '$p_{\epsilon}$ (deg/s)', '$p_{\zeta}$ (deg/s)',
                  '$\phi_{\epsilon}$ (deg)', '$\phi_{\zeta}$ (deg)',
                  '$\psi_{\epsilon}$ (deg)', '$\psi_{\zeta}$ (deg)'
        ]

        for i in range(10):
            for j in range(self.aero._len_velocities):
                
                # step response
                if self.aero._many_velocities:
                    # t, y = signal.step(signal.TransferFunction(self.N[i][j], self.delta[j, :]), T=np.linspace(0, tmax, tmax*10))
                    result = ct.step_response(ct.tf(self.N[i][j], self.delta[j, :]), T=np.linspace(0, tmax, tmax*10))
                else:
                    # t, y = signal.step(signal.TransferFunction(self.N[i], self.delta), T=np.linspace(0, tmax, tmax*10))
                    result = ct.step_response(ct.tf(self.N[i], self.delta), T=np.linspace(0, tmax, tmax*10))

                # coloca as legendas e o nome de cada gráfico
                # ax[i//2][i%2].plot(t, y)
                ax[i//2][i%2].plot(result.time, result.outputs)
            ax[i//2][i%2].set(ylabel=titles[i])
            ax[i//2][i%2].grid()
            ax[i//2][i%2].set_xlim([0, tmax])
            # ax[i//2][i%2].legend([f"{self.aero.V0[jj]} m/s" for jj in range(len(self.aero.V))])
        ax[4][0].set(xlabel="t (s)")
        ax[4][1].set(xlabel="t (s)")

        fig.tight_layout()  # ajusta o tamanho

        return
    
    def root_map (self):
        '''
        Apresenta o gráfico dos polos e zeros do polinômio característico e retorna as frequências naturais e amortecidas
        '''
        # raizes do polinomio caracteristico
        r = np.array([np.roots(self.delta[i]) for i in range(self.aero._len_velocities)]) if self.aero._many_velocities else np.roots(self.delta)
        
        wd = np.sqrt(r.real**2 + r.imag**2)         # frequência natural
        cd = r.real / np.where(wd == 0, 1, wd)      # frequência de amortecimento

        plt.figure()
        if self.aero._many_velocities:
            for i in range(self.aero._len_velocities):
                plt.scatter(r[i].real, r[i].imag)
        else:
            plt.scatter(r.real, r.imag)

        plt.grid()
        plt.xlabel('$\sigma$ (rad/s)')
        plt.ylabel('$j \gamma$ (rad/s)')

        # legenda do gráfico
        plt.legend([f"{self.aero.V0[jj]} m/s" if self.aero._many_velocities else f"{self.aero.V0} m/s" for jj in range(self.aero._len_velocities)])

        return wd, cd

    # def get_Lf_Sf (self, ro, wd, cd):
    #     '''
    #     Retorna uma distância Lf e uma área Sf da empenagem vertical baseado nas frequências naturais\n
    #         ro : densidade do ar
    #         wd : frequência natural
    #         cd : frequência de amortecimento
    #     '''
    #     # baseado nas equações de aproximação das frequencias
    #     A = 2*(wd**2)*self.aero.Iz/(ro*(self.aero.V0**2)*self.aero.f.CLa)
    #     B = -self.aero.Iz/self.aero.f.CLa * (4*wd*cd*self.aero.m/(ro*self.aero.V0) + self.aero.b.Sl*self.aero.b.CDl)
    #     delta = np.sqrt(B**2 + 4*A**2*self.aero.m*self.aero.Iz)

    #     # distancia do CG até o CA da EV
    #     Lf = np.array([
    #         (B + delta)/(2*A*self.aero.m),
    #         (B - delta)/(2*A*self.aero.m)
    #     ])

    #     # área da EV
    #     Sf = np.array([
    #         (2*A**2*self.aero.m)/(B + delta),
    #         (2*A**2*self.aero.m)/(B - delta)
    #     ])
    #     return np.round(Lf, 3), np.round(Sf, 3)
    
    def get_Sf_for (self, ro, wd, cd):
        '''
        Retorna uma área Sf da empenagem vertical baseado nas frequências naturais\n
            ro : densidade do ar
            wd : frequência natural
            cd : frequência de amortecimento
        '''

        n = 4*cd*wd*self.aero.m*self.aero.Iz + \
            ro*self.aero.V0*(self.aero.b.Sl*self.aero.b.CDl*self.aero.Iz - \
                             self.aero.w.b*self.aero.w.c12[1]*self.aero.m*self.aero.w.get_CD(self.aero.alpha))
        
        d = ro*self.aero.V0*self.aero.f.CLa*(self.aero.Lf*self.aero.lf * self.aero.m + self.aero.Iz)

        return n/d


if __name__ == "__main__":
    import main