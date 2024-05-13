import math

def Tau (S_super, S_control):
    '''
    Retorna o valor do parâmetros TAU para superfícies de controle
        S_super : área total da empenagem
        S_control : área da superficie de controle (leme ou profundor)
    '''
    x = S_control / S_super

    # equação da figura 2.20 do NELSON
    return 12.267*(x**5) - 25.8*(x**4) + 21.477*(x**3) - 9.6446*(x**2) + 3.26*x + 0.01

#=======================================================================================================
class AeroSurface:
    def __init__(self, S: float, b: float, mac: float, inc: float, c12: list = None):
        '''
        Surface object, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            inc : ângulo de incidência (deg)
            c12 : cordas na raiz e na ponta [raiz, ponta]
        '''
        self.S = S
        self.b = b
        self.mac = mac
        self.inc = math.radians(inc)
        self.c12 = c12 if c12 != None else [mac, mac]

        self.x_CA = self.mac/4                          # posição do centro aerodinamico
        self.AR = self.mac * self.b / (self.S**2)       # alongamento
        self.lbd = self.c12[1] / self.c12[0]            # afilamento

        return
    
    def set_CL (self, CL0, CLa):
        '''
        Coeficientes aerodinâmicos:
            CLa : derivada do coef sustentação em relação a alpha
            CL0 : coef sustentação em alpha = 0
        '''
        self.CLa = CLa
        self.CL0 = CL0

        return
    
    def set_CD (self, CD0, CDa):
        '''
        Coeficientes aerodinâmicos:
            CDa : derivada do coef de arrasto em relação a alpha
            CD0 : coef de arrasto em alpha = 0
        '''
        self.CDa = CDa
        self.CD0 = CD0

        return

    def angulos (self, T: float, V: float):
        '''
        Ângulos:
            T : ang diedro (deg)
            V : ang enflexamento (deg)
        '''
        self.T = math.radians(T)    # angulo de diedro
        self.V = math.radians(V)    # angulo de enflexamento

        return
#=======================================================================================================
class Empennage (AeroSurface):
    def __init__(self, S: float, b: float, mac: float, inc: float, c12: list, l: float, h: float):
        '''
        Empennage object, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            inc : ângulo de incidência (deg)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            l : distância do CG até o CA da empenagem em x (m)
            h : distância do CG até o CA da empenagem em z (m)
        '''
        super().__init__(S, b, mac, inc, c12)

        self.l = l
        self.h = h

        return
    
    def vol_cauda (self, Sw: float, bw: float):
        '''
        Calcula o Volume de Cauda Vc da empenagem:
            Sw : área da asa (m)
            bw : envergadura da asa (m)
        '''
        self.Vc = self.l * self.S / (Sw * bw)

        return
#=======================================================================================================
class Wing(AeroSurface):
    def __init__(self, S: float, b: float, mac: float, inc: float, Cm_CA: float, c12: list = None):
        '''
        Surface object, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            inc : ângulo de incidência (deg)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            Cm_CA : coef de momento em torno do CA
        '''
        super().__init__(S, b, mac, inc, c12)

        self.Cm_CA = Cm_CA

        return
    
    def Cm_0 (self, x_CG):
        '''
        Retorna a contribuição da asa para o coeficiente de momento para alpha = 0
        '''
        self.Cm0 = self.Cm_CA + (x_CG - self.x_CA)/self.mac * (self.CL0 + self.CLa*self.inc)

        return self.Cm0
    
    def Cm_a (self, x_CG):
        '''
        Retorna a contribuição da asa para o coeficiente de momento em funçãod e alpha
        '''
        self.Cma = self.CLa * (x_CG - self.x_CA)/self.mac

        return self.Cma
#=======================================================================================================
class Finn(Empennage):
    def __init__(self, S: float, b: float, mac: float, inc: float, c12: list, l: float, h: float, k: int):
        '''
        Empenagem Vertical, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            inc : ângulo de incidência (deg)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            l : distância do CG até o CA da empenagem em x (m)
            h : distância do CG até o CA da empenagem em z (m)
            k : quantidade de EV's
        '''
        super().__init__(S, b, mac, inc, c12, l, h)

        self.k = k
        
        return
    
    def Cl_b (self):
        '''
        Retorna a contribuição da EV para o coeficiente de momento lateral em função de beta (sideslip)
        Depende do Volume de cauda Vc e do coef CL alpha
        '''
        self.Clb = - self.Vc * self.CLa * self.h / self.l
        # equação 3.39 do COOK

        return self.Clb
    
    def Cn_b (self, kn, krl, Vf, ETA):
        '''
        Retorna o coeficiente de momento angular de guinada
            kn : fator empírico de interferência asa-fuselagem, função direta da  geometria da fuselagem
            krl : fator empírico do numero de Reynolds  da fuselagem
            Vf : "volume de cauda fuselagem" 
        '''
        self.Cnb = - kn * krl * Vf + self.Vc * self.CLa * ETA
        # equação 5.91 TAPERA

        return self.Cnb
    
    def Cn_dr (self, ETA, S_rudder):
        '''
        Retorna o coeficiente de momento direcional do leme em função da deflexão delta
            ETA : eficiencia de cauda
            S_rudder : área total do leme
        '''
        self.Cndr = - self.k * ETA * self.Vc * self.CLa * Tau(self.S, S_rudder)

        return self.Cndr
#=======================================================================================================
class Tail (Empennage):
    def __init__(self, S: float, b: float, mac: float, inc: float, c12: list, l: float, h: float):
        '''
        Empenagem Horizontal, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            inc : ângulo de incidência (deg)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            l : distância do CG até o CA da empenagem em x (m)
            h : distância do CG até o CA da empenagem em z (m)
        '''
        super().__init__(S, b, mac, inc, c12, l, h)
        return
    
    def Cm_a (self, ETA, de_da):
        '''
        Retorna a contribuição da EH para o coeficiente de momento em função de alpha
            ETA : eficiencia de cauda
            de_da : derivada do downwash em função de alpha
        '''
        self.Cma = - ETA * self.Vc * self.CLa * (1 - de_da)

        return self.Cma
    
    def Cm_0 (self, ETA, e_0):
        '''
        Retorna a contribuição da EH para o coeficiente de momento em alpha = 0
        '''
        self.Cm0 = - ETA * self.Vc * (self.CL0 + self.CLa * (self.inc - e_0))

        return self.Cm0
    
    def Cm_de (self, ETA, S_elevator):
        '''
        Retorna o coef de momento do profundor em função da deflexão delta
        '''
        self.Cmde = - ETA * Tau(self.S, S_elevator) * self.Vc * self.CLa

        return self.Cmde
#=======================================================================================================
class Aircraft:
    def __init__(self, m: float, Ix: float, Ixz: float, Iz: float, V0: float, theta_e: float):
        '''
        Aeronave:
            m : massa (kg)
            Ix : momento de inércia Ix eixo aeronautico (kg*m^2)
            Iz : momento de inércia Iz eixo aeronautico (kg*m^2)
            Ixz : momento de inércia Ixz eixo aeronautico (kg*m^2)
            ro : densidade do ar (kg/m^3)
            V0 : velocidade da aeronave (m/s)
            theta_e : ângulo de equilíbrio entre horizonte e direção de voo (deg)
        '''
        self.m = m
        self.Ix = Ix
        self.Ixz = Ixz
        self.Iz = Iz

        self.V0 = V0
        self.theta_e = math.radians(theta_e)

    def ad_mass (self, ro: float, Sw: float, bw: float):
        '''
        Admensionaliza a massa e os momentos de inércia em relação á asa:
            ro : densidade do ar (kg/m^3)
            Sw : área da asa (m^2)
            bw : envergadura da asa (m)
        '''

        # adimensionalizadores
        ad1 = 0.5*ro*self.V0*Sw
        ad2 = ad1*bw

        self.m1 = self.m/ad1                      # massa adimensional
        self.Ix1 = self.Ix/ad2                    # momento de inercia Ix adimensional
        self.Iz1 = self.Iz/ad2                    # momento de inercia Iz adimensional
        self.Ixz1 = self.Ixz/ad2                  # momento de inercia Ixz adimensional
        
        return
#=======================================================================================================
if __name__ == "__main__":
    import main