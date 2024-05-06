import math

class Wing:
    def __init__(self, S:float, b:float, mac:float, inc:float, c12:list):
        '''
        Wing object, parâmetros geométricos:
            S : área
            b : envergadura
            mac : corda média aerodinamica
            inc : incidência (deg)
            c12 : cordas na raiz e na ponta [raiz, ponta]
        '''
        self. S = S
        self.b = b
        self.mac = mac
        self.inc = math.radians(inc)
        self.c12 = c12

        self.CA = self.mac/4                        # posição do centro aerodinamico
        self.AR = self.mac * self.b / (self.S**2)   # alongamento
        self.lbd = self.c12[1] / self.c12[0]         # afilamento

        return

    def coefs (self, CL, CD, Cma):
        '''
        Coeficientes aerodinâmicos:
            CL : coef sustentação
            CD : coef arrasto
            Cma : coef de momento entorno do centro aerodinamico
        '''
        self.CD = CD
        self.CL = CL
        self.Cma = Cma

        return

    def angulos (self, T, V):
        '''
        Ângulos:
            T : ang diedro
            V : ang enflexamento
        '''
        self.T = T      # angulo de diedro
        self.V = V      # angulo de enflexamento

        return
#=======================================================================================================
class Empennage (Wing):
    def __init__(self, S, b, mac, inc, c12: list, l, h):
        '''
        Empennage object, parâmetros geométricos:
            l : distância do CG até o CA da empenagem em x
            h : distância do CG até o CA da empenagem em z
        '''
        super().__init__(S, b, mac, inc, c12)

        self.l = l
        self.h = h

        return
    
    def vol_cauda (self, Sw, bw):
        '''
        Calcula o Volume de Cauda Vc da empenagem:
            Sw : área da asa
            bw : envergadura da asa
        '''
        self.Vc = self.l * self.S / (Sw * bw)

        return
#=======================================================================================================
class Finn(Empennage):
    def __init__(self, S, b, mac, inc, c12: list, l, h, k):
        '''
        Empenagem Vertical, parâmetros geométricos:
            S : área
            b : envergadura
            mac : corda média aerodinamica
            inc : incidência (deg)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            l : distância do CG até o CA da empenagem em x
            h : distância do CG até o CA da empenagem em z
            k : quantidade de EV's
        '''
        super().__init__(S, b, mac, inc, c12, l, h)

        self.k = k
        
        return
#=======================================================================================================
class Tail (Empennage):
    def __init__(self, S, b, mac, inc, c12: list, l, h):
        '''
        Empenagem Horizontal, parâmetros geométricos:
            S : área
            b : envergadura
            mac : corda média aerodinamica
            inc : incidência (deg)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            l : distância do CG até o CA da empenagem em x
            h : distância do CG até o CA da empenagem em z
        '''
        super().__init__(S, b, mac, inc, c12, l, h)
        return
#=======================================================================================================
class Aircraft:
    def __init__(self,
                 m : float,
                 Ix : float,
                 Ixz : float,
                 Iz : float,
                 V0 : float,
                 ro : float,
                 theta_e : float,
                 Asa : Wing,
                 EV : Finn,
                 # EH : Tail
        ):
        '''
        Aeronave:
            m : massa (kg)
            Ix : momento de inércia Ix eixo aeronautico (kg*m^2)
            Iz : momento de inércia Iz eixo aeronautico(kg*m^2)
            Ixz : momento de inércia Ixz eixo aeronautico(kg*m^2)
            ro : densidade do ar (kg/m^3)
            V0 : velocidade da aeronave (m/s)
            theta_e : ângulo de equilíbrio entre horizonte e direção de voo (deg)
        '''
        self.m = m
        self.V0 = V0
        self.theta_e = math.radians(theta_e)

        self.w = Asa
        self.f = EV
        # self.t = EH

        self.f.vol_cauda(self.w.S, self.w.b)
        # self.t.vol_cauda(self.w.S, self.w.b)

        # adimensionalizadores
        ad1 = 0.5*ro*V0*Asa.S
        ad2 = ad1*Asa.b

        self.m = m/ad1                      # massa adimensional
        self.Ix = Ix/ad2                    # momento de inercia Ix adimensional
        self.Iz = Iz/ad2                    # momento de inercia Iz adimensional
        self.Ixz = Ixz/ad2                  # momento de inercia Ixz adimensional
        
        return
#=======================================================================================================
if __name__ == "__main__":
    import main