import math
import numpy as np
from pandas import DataFrame
import Util

#=======================================================================================================
class AeroSurface:
    def __init__(self, S: float, b: float, mac: float, c12: list = None):
        '''
        Surface object, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            c12 : cordas na raiz e na ponta [raiz, ponta]
        '''
        self.S = S
        self.b = b
        self.mac = mac
        self.c12 = c12 if c12 != None else [mac, mac]

        self.x_CA = self.mac/4                          # posição do centro aerodinamico
        self.AR = (self.b**2) / self.S                  # alongamento
        self.lbd = self.c12[1] / self.c12[0]            # afilamento

        self.set_angles()   # ângulos = 0

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
    
    def get_CD (self, alpha):
        # Oswald's factor
        e = 1/(1.05 + 0.007*math.pi*self.AR)    # Obert 5 <= AR <= 25
        # e = 1.78*(1 - 0.045*self.AR**0.68) - 0.64   # internet
        
        # pode ser deduzida a partir da equação 3.1 do "Methods for estimating stability and control\
        #  derivatives of conventional subsonic airplanes" de Jan Roskam
        return self.CD0 + (self.get_CL(alpha)**2 - self.CL0**2)/(math.pi*self.AR*e)
    
    def get_CL (self, alpha):
        return self.CL0 + self.CLa * (alpha + self.inc)

    def set_angles (self, T: float = 0, V_c4: float = 0, V_LE: float = 0, inc: float = 0):
        '''
        Ângulos:
            T : ang diedro (deg)
            V_c4 : ang enflexamento a 1/4 da corda (deg)
            V_LE : ang de enflexamento no bordo de ataque (deg)
            inc : ang de incidência (deg)
        '''
        self.T = math.radians(T)            # angulo de diedro
        self.V_c4 = math.radians(V_c4)      # angulo de enflexamento
        self.V_LE = math.radians(V_LE)      # angulo de enflexamento no bordo de ataque
        self.inc = math.radians(inc)        # angulo de incidência

        return
    
    def estimate_CLa (self, k: float, V_LE: float, M: float):
        '''
        Estima o valor de CLa com base em métodos paramétricos presentes no "Methods for estimating stability and control derivatives\
        of conventional subsonic airplanes" de Jan Roskam:
            k : ratio of actual average wing section lift curve slope, CLa to 2pi
            V_LE : ângulo de enflexamento no bordo de ataque (rad)
            M : numero de mach
        '''
        # equação 2.3
        tg_V_c2 = np.tan(V_LE) - ((1 - self.lbd)/(1 + self.lbd))*2/self.AR

        # equação 3.8 + 3.9 modificada
        self.CLa = 2*np.pi*self.AR/(2 + np.sqrt((1 - M**2 + tg_V_c2**2)*(self.AR/k)**2 + 4))
        return
    
    def estimate_CDa (self, ro, V, m):
        '''
        Estima o valor de CDa com base em métodos paramétricos presentes no "Methods for estimating stability and control derivatives\
        of conventional subsonic airplanes" de Jan Roskam:
            ro : densidade do ar (kg/m^3)
            V : velocidade (m/s)
            m : massa da aeronave (kg)
        '''
        # Oswald's factor
        e = 1/(1.05 + 0.007*math.pi*self.AR)    # Obert
        # e = 1.78*(1 - 0.045*self.AR**0.68) - 0.64   # internet

        # equação 3.2 + 3.3 + 3.4 modificada
        self.CDa = 4*m*self.CLa/(ro*V**2 * self.S*math.pi*self.AR*e)
        return
#=======================================================================================================    
class ControlSurface:
    def __init__(self, S: float, c: float):
        '''
        Control object, parâmetros geométricos
            S : área (m^2)
            c : corda (m)
        '''
        self.S = S
        self.c = c
        return
    
    def set_CLa (self, CLa):
        '''
        Coef de sustentaçãoe em função de alpha
        '''
        self.CLa = CLa
        return
    
    def TAU (self, S):
        '''
        Retorna o valor do parâmetros TAU para superfícies de controle
            S : área total da superfície onde o controle está
        '''
        x = self.S / S
        # equação da figura 2.20 do NELSON
        return 12.267*(x**5) - 25.8*(x**4) + 21.477*(x**3) - 9.6446*(x**2) + 3.26*x + 0.01
#=======================================================================================================
class Rudder (ControlSurface):
    def __init__(self, S: float, c: float):
        '''
        Define as dimensioes básicas do leme:
            S : área (m^2)
            c : corda (m)
        '''
        super().__init__(S, c)
        return
    
    def set_CLa (self, AR_f, CLa, f = 1):
        '''
        Correção de CLa do leme devido ao efeito geométrico da EV\
        usando o método descrito em Abbot and Von Doenhoff (1959)
            AR_f : alongamento da EV
            CLa : coef de sustentação em função de alpha
            f : fator de correção empirico
        '''
        self.CLa = f*CLa/(1 + CLa/(math.pi*AR_f))
        # equação 13.254 COOK
        return
#=======================================================================================================
class Elevator (ControlSurface):
    def __init__(self, S: float, c: float):
        '''
        Define as dimensioes básicas do profundor:
            S : área (m^2)
            c : corda (m)
        '''
        super().__init__(S, c)
        return
#=======================================================================================================
class Aileron (ControlSurface):
    def __init__(self, S: float, c: float, y1: float, y2: float):
        '''
        Define as dimensioes básicas do aileron
            S : área (m^2)
            c : corda
            y1 : distância perpendicular do eixo x até o começo do aileron
            y2 : distância perpendicular do eixo x até o final do aileron
        '''
        super().__init__(S, c)
        self.y1 = y1
        self.y2 = y2
        return
#=======================================================================================================
class Wing (AeroSurface):
    def __init__(self, S: float, b: float, mac: float, c12: list = None, Cm_CA: float = 0):
        '''
        Surface object, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            Cm_CA : coef de momento em torno do CA da asa
        '''
        super().__init__(S, b, mac, c12)

        self.Cm_CA = Cm_CA

        return
    
    def downwash (self):
        '''
        Retorna derivada de_da do angulo de downwash por alpha
        '''
        # equação 2.23 NELSON
        return 2*self.CLa / (math.pi * self.AR)
#=======================================================================================================
class Fin (AeroSurface):
    def __init__(self, S: float, b: float, c12: list = None, k: int = 1):
        '''
        Empenagem Vertical, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            k : quantidade de EV's
        '''
        super().__init__(S, b, (c12[0] + c12[1])/2, c12)
        
        self.V_c4 = math.atan(3*(c12[0] - c12[1])/(4*b))            # ângulo de enflexamento em 1/4 de corda
        self.V_LE = math.atan((self.c12[0] - self.c12[1])/self.b)   # ângulo de enflexamento no bordo de ataque
        self.k = k
        
        return
    
    def effective_AR (self, AR_B_AR, AR_HB_AR_B, KH):
        '''
        Calcula o AR effetivo da EV, obtido no Methods for estimating stability and control derivatives\
         of conventional subsonic airplanes by ROSKAM
            AR_B_AR : razão do AR na preseça do corpo com o AR da EV isolada (figura 7.5)
            AR_HB_AR_B : razão do AR na presença da EH e do corpo com o AR só na presença do corpo (figura 7.6)
            KH : fator que conta com o tamanho relativo da EH e da EV (figura 7.7)
        '''
        # equação 7.6 do Methods for estimating stability and control derivatives of conventional\
        #  subsonic airplanes by ROSKAM
        return AR_B_AR * self.AR * (1 + KH* (AR_HB_AR_B - 1))
    
    def chord(self, h):
        '''
        Retorna a corda local na coordenada h de uma EV alinhada no bordo de fuga
        '''
        return self.c12[0] + (self.c12[1] - self.c12[0])*h/self.b
#=======================================================================================================
class Tail (AeroSurface):
    def __init__(self, S: float, b: float, mac: float, c12: list):
        '''
        Empenagem Horizontal, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            c12 : cordas na raiz e na ponta [raiz, ponta]
        '''
        super().__init__(S, b, mac, c12)

        return
#=======================================================================================================
class Body:
    def __init__(self, Sl: float, h: float):
        '''
        Fuselagem
            Sl : área lateral projetada da fuselagem
            h : altura
        '''
        self.Sl = Sl
        self.h = h

        self.CDl = -(0.00714 + 0.674*h**2/Sl)   # estimation of sideforce due to sideslip for body
        return
#=======================================================================================================
class Aircraft:
    def __init__(self, wing: Wing, fin: Fin, tail: Tail, body: Body, V):
        '''
        Aircraft:
            wing : asa obj
            fin : empenagem vertical obj
            tail : empenagem horizontal obj
            body : fuselagem obj
            V : velocidade (m/s)
        '''
        self.w = wing
        self.f = fin
        self.t = tail
        self.b = body
        self.V = V

        self._many_velocities = True if type(self.V) == np.ndarray else False
        self._len_velocities = len(self.V) if self._many_velocities else 1

        return
    
    def derivatives (self, dCL_day, dCD_day, dCL_dah, dCDy_de, CDy, CLy, cy: list, ch: list):
        '''
        Calcula as derivadas adimensionais de estabilidade presentes no apêndice 8 do livro
            dCL_day : derivada do coef de sustentação em função de alpha na coordenada y
            dCD_day : derivada do coef de arrasto em função de alpha na coordenada y
            dCL_dah : derivada do coef de sustentação em função de alpha na coordenada h
            dCDy_de : derivada do coef de arrasto em função de xi (deflexão do aileron) na coordenada y
            CDy : coef de arrasto local na coordenada y
            CLy : coef de sustentação local na coordenada y
            cy : corda local na coordenada y
            ch : corda local na coordenada h
        '''
        def Lv_int1 ():
            return cy * dCL_day * self.w.T * y          # Wing with Dihedral
        
        def Lv_int2 ():
            return cy * y                               # wing with aft sweep
        
        def Yp_int ():
            return dCL_dah * ch * h                     # fin contribution
        
        def Lp_int ():
            return (dCL_day + CDy)*cy * y**2            # wing contribuiton
        
        def Np_int ():
            return (CLy - dCD_day)*cy * y**2            # wing contribution
        
        def Lr_int ():
            return CLy * cy * y**2                      # wing contribution
        
        def Nr_int ():
            return CDy * cy * y**2                      # wing contribution
        
        def Le_int ():
            return cy * ya
        
        def Ne_int ():
            return dCDy_de * cy * ya

        s = self.w.b/2               # pra facilitar

        n1 = len(cy)    # size of wing chord data
        n2 = len(ch)    # size of fin chord data

        y = np.linspace(0, s, n1)
        ya = np.linspace(self.a.y1, self.a.y2, n1)
        h = np.linspace(0, self.f.b, n2)

        # derivadas de "v" (sideslip)==============================================================
        # Side force
        self.Yv = (self.b.Sl*self.b.CDl - self.f.S*self.f.CLa)/self.w.S

        # Rolling moment
        self.Lv = -Util.trapezoidal(Lv_int1(), 0, s, n1)/(self.w.S * s)\
                - 2*self.CLe*math.tan(self.w.V_c4)*Util.trapezoidal(Lv_int2(), 0, s, n1)/(self.w.S * s)\
                - self.f.CLa * self.Vv * (self.hf/self.Lf)

        # Yawing moment
        self.Nv = self.f.CLa * self.Vv    #??

        # derivadas de "p" (roll rate)=============================================================
        # Side force 
        self.Yp = -Util.trapezoidal(Yp_int(), 0, self.f.b, n2)/(self.w.S * self.w.b)

        # Rolling moment
        self.Lp =  -Util.trapezoidal(Lp_int(), 0, s, n1)/(2*self.w.S * s**2)

        # Yawing moment
        self.Np = -Util.trapezoidal(Np_int(), 0, s, n1)/(2*self.w.S * s**2)

        # derivadas de "r" (yaw rate)==============================================================
        # Side force
        self.Yr = self.Vv*self.f.CLa      #??

        # Rolling moment
        self.Lr = Util.trapezoidal(Lr_int(), 0, s, n1)/(self.w.S * s**2)\
                + self.f.CLa*self.Vv*(self.hf/self.w.b)

        # Yawing moment
        self.Nr = -Util.trapezoidal(Nr_int(), 0, s, n1)/(self.w.S * s**2)\
                - self.f.CLa*self.Vv*(self.Lf/self.w.b)
        
        # derivadas de "e" (aileron)===============================================================
        # Side force
        self.Ye = 0

        # Rolling moment
        self.Le = -self.a.CLa*Util.trapezoidal(Le_int(), self.a.y1, self.a.y2, n1)/(self.w.S * s)
        
        # Yawing moment
        self.Ne = Util.trapezoidal(Ne_int(), self.a.y1, self.a.y2, n1)/(self.w.S * s)

        # derivadas de "c" (rudder)================================================================
        # Side force
        self.Yc = self.f.S*self.r.CLa/self.w.S

        # Rolling moment
        self.Lc = self.Vv*self.r.CLa*(self.hf/self.Lf)

        # Yawing moment
        self.Nc = -self.Vv*self.r.CLa

        # Dimensionalização =======================================================================
        # side force
        self.Yv1 = self.Yv * self._ad1
        self.Yp1 = self.Yp * self._ad2
        self.Yr1 = self.Yr * self._ad2
        self.Ye1 = self.Ye * self._ad4
        self.Yc1 = self.Yc * self._ad4

        # rolling moment
        self.Lv1 = self.Lv * self._ad2
        self.Lp1 = self.Lp * self._ad3
        self.Lr1 = self.Lr * self._ad3
        self.Le1 = self.Le * self._ad5
        self.Lc1 = self.Lc * self._ad5

        # yawing moment
        self.Nv1 = self.Nv * self._ad2
        self.Np1 = self.Np * self._ad3
        self.Nr1 = self.Nr * self._ad3
        self.Ne1 = self.Ne * self._ad5
        self.Nc1 = self.Nc * self._ad5

        return

    def set_angles (self, alpha, theta):
        '''
        Ângulos:
            alpha : ang de ataque (deg)
            theta : ang de arfagem (deg)
        '''
        self.alpha = np.radians(alpha)
        self.theta = np.radians(theta)
        return
    
    def set_control (self, aileron: Aileron, elevator: Elevator, rudder : Rudder):
        '''
        Recebe as superfícies de controle da aeronave
            aileron : aileron obj
            elevator : profundor obj
            rudder : leme obj
        '''
        self.a = aileron
        self.e = elevator
        self.r = rudder
        return
    
    def set_fin (self, lf: float, Lf: float, hf: float):
        '''
        Informa a posição da EV
            lf : distância do CA da asa até o CA da EV em x (m)
            Lf : distância do CG até o CA da EV em x (m)
            hf : distância do eixo x até o CA da EV em z (m)
        '''
        self.lf = lf
        self.Lf = Lf
        self.hf = hf

        # volume de cauda vertical
        self.Vv = self.Lf * self.f.S / (self.w.S * self.w.b)

        # # volume de cauda modificado para EV que possibilita empiricamente ajuste dinâmico
        # self.Vv = self.f.S*(self.f.mac+0.7*self.z*math.tan(self.V))/(Sw*bw)
        # # equação 13.244 do COOK
        return
    
    def set_mass (self, ro:float, mass: float, Ix: float, Iz: float, Ixz: float):
        '''
        Informa a massa e os momentos de inércia
            ro : densidade do ar (kg/m^3)
            mass : massa da aeronave (kg)
            Ix : momento de inércia Ix eixo aeronautico (kg*m^2)
            Iz : momento de inércia Iz eixo aeronautico (kg*m^2)
            Ixz : momento de inércia Ixz eixo aeronautico (kg*m^2)
        '''
        # dimensionais
        self.m = mass
        self.Ix = Ix
        self.Iz = Iz
        self.Ixz = Ixz

        # adimensionalizadores
        self._ad1 = 0.5*ro*self.V*self.w.S
        self._ad2 = self._ad1*self.w.b
        self._ad3 = self._ad2*self.w.b
        self._ad4 = self._ad1*self.V
        self._ad5 = self._ad4*self.w.b

        self.m1 = self.m/self._ad1        # massa adimensional
        self.Ix1 = self.Ix/self._ad2      # momento de inercia Ix adimensional
        self.Iz1 = self.Iz/self._ad2      # momento de inercia Iz adimensional
        self.Ixz1 = self.Ixz/self._ad2    # momento de inercia Ixz adimensional

        return
    
    def set_tail (self, lt: float, Lt: float, ht: float):
        '''
        Informa a posição da EH
            lt : distância do CA da asa até o CA da EH em x (m)
            Lt : distância do CG até o CA da EH em x (m)
            ht : distância do eixo x até o CA da EH em z (m)
        '''
        self.lt = lt
        self.Lt = Lt
        self.ht = ht

        # volume de caudda horizontal
        self.Vh = self.Lt * self.t.S / (self.w.S * self.w.mac)
        return
    
    def estimate_Coefs (self, k: float, T: float, ro : float):
        '''
        Estima os valores dos coeficientes com base em métodos paramétricos presentes no "Methods for estimating stability and control derivatives\
        of conventional subsonic airplanes" de Jan Roskam:
            k : ratio of actual average wing section lift curve slope, CLa to 2pi
            T : temperatura (°C)
            ro : densidade do ar (kg/m^3)
        '''
        # CLa total da aeronave
        self.w.estimate_CLa(k, self.w.V_LE, Util.mach(self.V, T))
        self.f.estimate_CLa(k, self.f.V_LE, Util.mach(self.V, T))
        # self.t.estimate_CLa(k, self.t.V_LE, Util.mach(self.V, T))

        # CDa total da aeronave
        self.w.estimate_CDa(ro, self.V, self.m)
        self.f.estimate_CDa(ro, self.V, self.m)
        # self.t.estimate_CDa(ro, self.V, self.m)

        self.CLe = 2*self.m/(ro*self.V**2 * self.w.S)
        # equação 3.3 Roskam
        
        return
    
    def get_CL_eq (self):
        '''
        Retorna CL de equilíbrio da aeronave
        '''
        return self.CLe
    
    def get_CLa (self):
        '''
        Retorna CLa da aeronave
        '''
        return self.w.CLa + self.t.CLa
    
    def get_CDa (self):
        '''
        Retorna CDa da aeronave
        '''
        return self.w.CDa + self.f.CDa + self.t.CDa
    
    def get_derivatives (self):
        '''
        Retorna um dataframe com a velocidade e as derivadas
        '''
        # lists with names and orders
        name = ['V0', 'Lv', 'Lp', 'Lr', 'Le', 'Lc',\
                'Nv', 'Np', 'Nr', 'Ne', 'Nc',\
                'Yv', 'Yp', 'Yr', 'Ye', 'Yc']
        order = [self.V, self.Lv, self.Lp, self.Lr, self.Le, self.Lc,\
                self.Nv, self.Np, self.Nr, self.Ne, self.Nc,\
                self.Yv, self.Yp, self.Yr, self.Ye, self.Yc]
        
        if self._many_velocities:
            # make sure that the size stays homogeneous even with float inputs
            deri = np.array(list(map(
                lambda x: x if type(x) == np.ndarray else np.array([x for _ in range(len(self.V))]),
                order))).transpose()
        else:
            deri = np.array([order])

        return DataFrame(np.round(deri, 4), columns= name)
    
    def curve (self, phi):
        '''
        Retorna o tempo para completar uma curva e a taxa de curvatura
            phi : angulo de rolagem da aeronave (deg)
        '''
        t_c = 2*np.pi*self.V/(9.81*np.tan(np.radians(phi)))
        # equação 2.25 do COOK

        turn_rate = 2*np.pi/t_c
        # equação 2.26 do COOK

        return t_c, turn_rate
#=======================================================================================================
if __name__ == "__main__":
    import main