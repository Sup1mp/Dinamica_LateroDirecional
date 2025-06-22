import math
import numpy as np
from pandas import DataFrame
from Util.util import trapezoidal
from Util import cFit

#=======================================================================================================
class AeroSurface:
    def __init__(self, S: float, b: float, mac: float, c12: list = None, th: float = None):
        '''
        Surface object, parâmetros geométricos:\n
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            th : esperssura do perfil em porcentagem
        '''
        self.S = S
        self.b = b
        self.mac = mac
        self.c12 = c12 if c12 != None else [mac, mac]
        self.th = th*self.mac if th != None else th

        self.__params__()   # calcula os parâmetros principais de um perfil aerodinâmico

        self.set_angles(0,0,0,0)   # ângulos = 0
        self.__coefs__()    # inicializa coeficientes nulos

        return
    
    def __coefs__(self):
        # inicializa coeficientes como nulos
        self.CL0 = self.CLa = None
        self.CD0 = self.CDa = None
        self._noCDset = True   # indica se CD foi colocado pelo usuário ou não
        self._noCLset = True   # indica se CL foi colocado pelo usuário ou não
        return
    
    def __params__(self):
        # calcula os parâmetros principais de um perfil aerodinâmico
        self.x_CA = self.mac/4                          # posição do centro aerodinamico
        self.AR = (self.b**2) / self.S                  # alongamento
        self.lbd = self.c12[1] / self.c12[0]            # afilamento
        return

    def set_CL (self, CL0, CLa):
        '''
        Coeficientes aerodinâmicos:\n
            CLa : derivada do coef sustentação em relação a alpha (1/deg)
            CL0 : coef sustentação em alpha = 0 (1/deg)
        '''
        self.CLa = CLa*180/np.pi
        self.CL0 = CL0*180/np.pi

        self._noCLset = False   # indica que o CL foi setado pelo usuário

        return
    
    def set_CD (self, CD0, CDa):
        '''
        Coeficientes aerodinâmicos:\n
            CDa : derivada do coef de arrasto em relação a alpha (1/deg)
            CD0 : coef de arrasto em alpha = 0  (1/deg)
        '''
        self.CDa = CDa*180/np.pi
        self.CD0 = CD0*180/np.pi

        self._noCDset = False   # indica que o CD foi setado pelo usuário

        return
    
    def get_CD (self, alpha):
        '''
        Retorna CD (1/rad) para ângulo de ataque alpha (deg)
        '''
        # # Oswald's factor
        # e = 1/(1.05 + 0.007*math.pi*self.AR)    # Obert 5 <= AR <= 25
        # # e = 1.78*(1 - 0.045*self.AR**0.68) - 0.64   # internet
        
        # # pode ser deduzida a partir da equação 3.1 do "Methods for estimating stability and control\
        # #  derivatives of conventional subsonic airplanes" de Jan Roskam
        # return self.CD0 + ((self.CLa * (alpha + self.inc))**2)/(math.pi*self.AR*e)
        return self.CD0 + self.CDa*(np.radians(alpha) + self.inc)
    
    def get_CL (self, alpha):
        '''
        Retorna CL (1/rad) para ângulo de ataque alpha (deg)
        '''
        return self.CL0 + self.CLa * (np.radians(alpha) + self.inc)

    def set_angles (self, T: float = None, V_c4: float = None, V_LE: float = None, inc: float = None):
        '''
        Ângulos da superfície:\n
            T : ang diedro (deg)
            V_c4 : ang enflexamento a 1/4 da corda (deg)
            V_LE : ang de enflexamento no bordo de ataque (deg)
            inc : ang de incidência (deg)
        '''
        # só modifica o ângulo se desejado, caso contrário permanece com o valor que já tinha
        self.T = math.radians(T) if T != None else self.T               # angulo de diedro
        self.V_c4 = math.radians(V_c4) if V_c4 != None else self.V_c4   # angulo de enflexamento
        self.V_LE = math.radians(V_LE) if V_LE != None else self.V_LE   # angulo de enflexamento no bordo de ataque
        self.inc = math.radians(inc) if inc != None else self.inc       # angulo de incidência

        return
    
    def estimate_CLa (self, K: float = 0.9, M = 0):
        '''
        Estima o valor de CLa com base em métodos paramétricos presentes no "Methods for estimating stability and control derivatives\
        of conventional subsonic airplanes" de Jan Roskam:\n
            K : fator obtido na figura B.1,1a do Ektins
            M : numero de mach
        '''
        if self._noCLset:
            # equação 2.3
            tg_V_c2 = np.tan(self.V_LE) - ((1 - self.lbd)/(1 + self.lbd))*2/self.AR

            # kappa com CLa teórico
            k = 1.05*K*cFit.Cla_t(self.th/self.mac)/(2*np.pi)
            # ratio of actual average wing section lift curve slope, CLa to 2pi

            # equação 3.8 + 3.9 modificada
            self.CLa = 2*np.pi*self.AR/(2 + np.sqrt((1 - M**2 + tg_V_c2**2)*(self.AR/k)**2 + 4))
        return
    
    def estimate_CDa (self, ro, V0, m):
        '''
        Estima o valor de CDa com base em métodos paramétricos presentes no "Methods for estimating stability and control derivatives\
        of conventional subsonic airplanes" de Jan Roskam:\n
            ro : densidade do ar (kg/m^3)
            V0 : velocidade (m/s)
            m : massa da aeronave (kg)
        '''
        if self._noCDset:
            # Oswald's factor
            e = 1/(1.05 + 0.007*math.pi*self.AR)    # Obert
            # e = 1.78*(1 - 0.045*self.AR**0.68) - 0.64   # internet

            # equação 3.2 + 3.3 + 3.4 modificada
            self.CDa = 4*m*self.CLa/(ro*V0**2 * self.S*math.pi*self.AR*e)
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

        self._noCLdset = True
        self._noCLaset = True
        return
    
    def set_CLd (self, CLd):
        '''
        Coeficientes aerodinâmicos:\n
            CLd : derivada do coef de arrasto em relação a delta (1/deg)
        '''
        self.CLd = CLd*180/np.pi
        self._noCLdset = False
        return
    
    def set_CLa (self, CLa):
        '''
        Coeficientes aerodinâmicos:\n
            CLa : derivada do coef sustentação em relação a alpha (1/deg)
        '''
        self.CLa = CLa*180/np.pi
        self._noCLaset = False
        return
    
    def estimate_CLa (self, surf, M = 0):
        '''
        Coef de sustentação em função de alpha:\n
            surf : superfície aerodinâmica onde o controle se encontra
            M : numero de mach
        '''
        if self._noCLaset:
            beta = np.sqrt(1 - M**2)    # Prandlt-Glauert compressibility factor
            
            K = 0.9
            # kappa com CLa teórico
            k = 1.05*K*cFit.Cla_t(surf.th/surf.mac)/(2*np.pi)
            # Apêndice B.1 do Ektins

            tg_V_c2 = np.tan(surf.V_LE) - ((1 - surf.lbd)/(1 + surf.lbd))*2/surf.AR
            # equação 2.3 do "Methods for estimating stability and control derivatives of conventional subsonic airplanes" de Jan Roskam
            
            self.CLa = 2*np.pi*self.S/(2 + np.sqrt((self.S*beta/k)**2 * (1 + (tg_V_c2/beta)**2) + 4 ))
            # equação da figura B.1,2 para voo subsonico
        
        return
    
    def estimate_CLd (self, surf, M = 0):
        '''
        Coef de sustentaçãoe em função da deflexão delta:\n
            surf : superfície aerodinâmica onde o controle se encontra
            M : numero de mach
        '''
        if self._noCLdset:
            t_c = surf.th/surf.mac      # razão da espessura pela corda do perfil aerodinâmico
            cf_c = self.c/surf.mac      # razão entre a corda da superfície de controle e a corda do perfil em que está

            K1 = cFit.K1(cf_c, surf.AR)     # fator Flap-chord
            # Figura B.2,2 do Etkins

            K = 0.9
            Cla = 1.05*K*cFit.Cla_t(t_c)/ np.sqrt(1 - M**2)
            # Equação B.1,1 do Etkins

            Cld = cFit.Cld(t_c, cf_c, self.CLa/Cla)     # Eficiencia de controle para fluxo incompressível de duas dimensões
            # Figura B.2,1 do Etkins

            self.CLd = Cld * (surf.CLa/self.CLa) * K1
            # apêndice B.2 do Etkins
        return
    
    def TAU (self, S):
        '''
        Retorna o valor do parâmetros TAU para superfícies de controle:\n
            S : área total da superfície onde o controle está
        '''
        x = self.S / S
        # equação da figura 2.20 do NELSON
        return 12.267*(x**5) - 25.8*(x**4) + 21.477*(x**3) - 9.6446*(x**2) + 3.26*x + 0.01
#=======================================================================================================
class Rudder (ControlSurface):
    def __init__(self, S: float, c: float):
        '''
        Define as dimensioes básicas do leme:\n
            S : área (m^2)
            c : corda média (m)
        '''
        super().__init__(S, c)
        return
    
    def estimate_CLd (self, surf, M = 0, f = 1):
        '''
        Coef de sustentaçãoe em função da deflexão delta do leme\
        corrigido pelo efeito geométrico da EV usando o método descrito\
        em Abbot and Von Doenhoff (1959):\n
            surf : superfície aerodinâmica onde o controle se encontra
            M : numero de mach
            f : fator de correção empirico
        '''
        super().estimate_CLd(surf, M)

        self.CLd = f*self.CLd/(1 + self.CLd/(math.pi*surf.AR))
        # equação 13.254 COOK
        return
#=======================================================================================================
class Elevator (ControlSurface):
    def __init__(self, S: float, c: float):
        '''
        Define as dimensioes básicas do profundor:\n
            S : área (m^2)
            c : corda (m)
        '''
        super().__init__(S, c)
        return
#=======================================================================================================
class Aileron (ControlSurface):
    def __init__(self, S: float, c: float, y1: float, y2: float):
        '''
        Define as dimensioes básicas do aileron:\n
            S : área (m^2)
            c : corda
            y1 : distância perpendicular do eixo x até o começo do aileron
            y2 : distância perpendicular do eixo x até o final do aileron
        '''
        super().__init__(S, c)
        self.y1 = y1
        self.y2 = y2
        return
    
    def estimate_CLd(self, surf, M = 0):
        '''
        Coef de sustentaçãoe em função da deflexão delta:\n
            surf : superfície aerodinâmica onde o controle se encontra
            M : numero de mach
        '''
        if self._noCLdset:
            super().estimate_CLd(surf, M)

            # Fator de envergadura para flaps
            K2 = cFit.K2(self.y2, surf.b, surf.lbd) - cFit.K2(self.y1, surf.b, surf.lbd)
            # Figura B.2,3 do Etkins

            self.CLd = self.CLd * K2
            # equação do apêndice B.2 do Etikins
        return
#=======================================================================================================
class Wing (AeroSurface):
    def __init__(self, S: float, b: float, mac: float, c12: list = None, th: float = None, Cm_CA: float = 0):
        '''
        Surface object, parâmetros geométricos:\n
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            th : esperssura do perfil em porcentagem
            Cm_CA : coef de momento em torno do CA da asa
        '''
        super().__init__(S, b, mac, c12, th)

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
    def __init__(self, S: float, b: float, c12: list, th: float = None, k: int = 1):
        '''
        Empenagem Vertical, parâmetros geométricos:\n
            S : área (m^2)
            b : envergadura (m)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            th : esperssura do perfil em porcentagem
            k : quantidade de EV's
        '''
        super().__init__(S, b, (c12[0] + c12[1])/2, c12, th)
        
        self.V_c4 = math.atan(3*(c12[0] - c12[1])/(4*b))            # ângulo de enflexamento em 1/4 de corda
        self.V_LE = math.atan((self.c12[0] - self.c12[1])/self.b)   # ângulo de enflexamento no bordo de ataque
        self.k = k
        
        return
    
    def effective_AR (self, AR_B_AR, AR_HB_AR_B, KH):
        '''
        Calcula o AR effetivo da EV, obtido no Methods for estimating stability and control derivatives\
         of conventional subsonic airplanes by ROSKAM:\n
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
    def __init__(self, S: float, b: float, mac: float, c12: list = None, th: float = None):
        '''
        Empenagem Horizontal, parâmetros geométricos:\n
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            th : esperssura do perfil em porcentagem
        '''
        super().__init__(S, b, mac, c12, th)

        return
#=======================================================================================================
class Body:
    def __init__(self, Sl: float, h: float):
        '''
        Fuselagem:\n
            Sl : área lateral projetada da fuselagem
            h : altura
        '''
        self.Sl = Sl
        self.h = h

        self.CDl = -(0.00714 + 0.674*h**2/Sl)   # estimation of sideforce due to sideslip for body
        return
#=======================================================================================================
class Aircraft:
    def __init__(self, wing: Wing, fin: Fin, tail: Tail, body: Body, V0: float|list):
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
        self.V0 = V0 if type(V0) == float else np.array(V0)

        self._many_velocities = True if type(self.V0) != float else False
        self._len_velocities = len(self.V0) if self._many_velocities else 1

        self.set_angles() # inicializa os ângulos como zero

        return
    
    def derivatives (self, dCL_day, dCD_day, dCL_dah, dCDy_de, CDy, CLy, cy: list, ch: list):
        '''
        Calcula as derivadas adimensionais de estabilidade presentes no apêndice 8 do livro do COOK:\n
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
            return cy * y_a
        
        def Ne_int ():
            return dCDy_de * cy * y_a

        s = self.w.b/2               # pra facilitar

        n1 = len(cy)    # size of wing chord data
        n2 = len(ch)    # size of fin chord data

        y = np.linspace(0, s, n1)
        h = np.linspace(0, self.f.b, n2)

        # cy_a = np.array([cy[i] if self.a.y1 <= y[i] <= self.a.y2 else 0 for i in range(n1)])
        y_a = np.array([y[i] if self.a.y1 <= y[i] <= self.a.y2 else 0 for i in range(n1)]) if n1 > 3 else y

        # derivadas de "v" (sideslip)==============================================================
        # Side force
        self.Yv = (self.b.Sl*self.b.CDl - self.f.S*self.f.CLa)/self.w.S

        # Rolling moment
        self.Lv = -trapezoidal(Lv_int1(), 0, s, n1)/(self.w.S * s)\
                - 2*self.CLe*math.tan(self.w.V_c4)*trapezoidal(Lv_int2(), 0, s, n1)/(self.w.S * s)\
                - self.f.CLa * self.Vv * (self.hf/self.Lf)

        # Yawing moment
        self.Nv = self.f.CLa * self.Vv    #??

        # derivadas de "p" (roll rate)=============================================================
        # Side force 
        self.Yp = -trapezoidal(Yp_int(), 0, self.f.b, n2)/(self.w.S * self.w.b)

        # Rolling moment
        self.Lp =  -trapezoidal(Lp_int(), 0, s, n1)/(2*self.w.S * s**2)

        # Yawing moment
        self.Np = -trapezoidal(Np_int(), 0, s, n1)/(2*self.w.S * s**2)

        # derivadas de "r" (yaw rate)==============================================================
        # Side force
        self.Yr = self.Vv*self.f.CLa      #??

        # Rolling moment
        self.Lr = trapezoidal(Lr_int(), 0, s, n1)/(self.w.S * s**2)\
                + self.f.CLa*self.Vv*(self.hf/self.w.b)

        # Yawing moment
        self.Nr = -trapezoidal(Nr_int(), 0, s, n1)/(self.w.S * s**2)\
                - self.f.CLa*self.Vv*(self.Lf/self.w.b)
        
        # derivadas de "e" (aileron)===============================================================
        # Side force
        self.Ye = 0

        # Rolling moment
        self.Le = -self.a.CLd*trapezoidal(Le_int(), self.a.y1, self.a.y2, n1)/(self.w.S * s)
        
        # Yawing moment
        self.Ne = trapezoidal(Ne_int(), self.a.y1, self.a.y2, n1)/(self.w.S * s)

        # derivadas de "c" (rudder)================================================================
        # Side force
        self.Yc = self.f.S*self.r.CLd/self.w.S

        # Rolling moment
        self.Lc = self.Vv*self.r.CLd*(self.hf/self.Lf)

        # Yawing moment
        self.Nc = -self.Vv*self.r.CLd

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

    def set_angles (self, alpha: float | list = 0, theta: float | list = 0):
        '''
        Ângulos da aeronave:\n
            alpha : ang. de ataque (deg)
            theta : ang. de arfagem (deg)
        '''
        if self._many_velocities:
            # tranforma angulos de float para listas com o mesmo tamanho de V0 usando o mesmo valor de float
            self.alpha = np.radians(alpha) if type(alpha) != float and type(alpha) != int else np.radians([alpha for _ in range(self._len_velocities)])
            self.theta = np.radians(theta) if type(theta) != float and type(theta) != int else np.radians([theta for _ in range(self._len_velocities)])
            
            if type(alpha) != float and type(alpha) != int:
                if len(alpha) != self._len_velocities :
                    # caso tamanho da lista de ângulos não seja o mesmo de velocidades, raise error
                    raise Exception("Alpha size don't mathce V0")
                
            if type(theta) != float and type(theta) != int:
                if len(theta) != self._len_velocities:
                    raise Exception("Theta size don't mathce V0")
        else:
            self.alpha = np.radians(alpha)
            self.theta = np.radians(theta)
        return
    
    def set_control (self, aileron: Aileron, elevator: Elevator, rudder : Rudder):
        '''
        Recebe as superfícies de controle da aeronave:\n
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
        Informa a posição da EV:\n
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
    
    def set_mass (self, ro:float | list, mass: float, Ix: float, Iz: float, Ixz: float):
        '''
        Informa a massa e os momentos de inércia:\n
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
        self._ad1 = 0.5*ro*self.V0*self.w.S
        self._ad2 = self._ad1*self.w.b
        self._ad3 = self._ad2*self.w.b
        self._ad4 = self._ad1*self.V0
        self._ad5 = self._ad4*self.w.b

        self.m1 = self.m/self._ad1        # massa adimensional
        self.Ix1 = self.Ix/self._ad2      # momento de inercia Ix adimensional
        self.Iz1 = self.Iz/self._ad2      # momento de inercia Iz adimensional
        self.Ixz1 = self.Ixz/self._ad2    # momento de inercia Ixz adimensional

        # coeficiente de sustentação de equilibrio da aeronave
        self.CLe = 2*self.m/(ro*self.V0**2 * self.w.S)
        # equação 3.3 Roskam

        return
    
    def set_tail (self, lt: float, Lt: float, ht: float):
        '''
        Informa a posição da EH:\n
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
    
    def estimate_CLa (self, k: float, M = 0):
        '''
        Estima os valores de CLa para todas as superfícies, até as de controle:\n
            k : ratio of actual average wing section lift curve slope, CLa to 2pi
            M : numero de mach
        '''

        # CLa total da aeronave
        self.w.estimate_CLa(k, M)      # Asa
        self.f.estimate_CLa(k, M)      # EV
        # self.t.estimate_CLa(k, M)      # EH

        self.a.estimate_CLa(self.w, M)      # Aileron
        self.r.estimate_CLa(self.f, M)      # Leme
        # self.e.estimate_CLa(self.t, M)      # Profundor
        return

    def estimate_CLd (self, M = 0):
        '''
        Estima os valores de CLd (Efetividade de Controle) para todas as superfícies de controle:\n
            M : numero de mach
        '''
        self.a.estimate_CLd(self.w, M)     # aileron
        self.r.estimate_CLd(self.f, M)     # leme
        # self.e.estimate_CLd(self.t)     # profundor
        return
    
    def estimate_CDa (self, ro):
        '''
        Estima os valores de CDa para todas as superfícies:\n
            ro : densidade do ar (kg/m^3)
        '''
        # CDa total da aeronave
        self.w.estimate_CDa(ro, self.V0, self.m)
        self.f.estimate_CDa(ro, self.V0, self.m)
        # self.t.estimate_CDa(ro, self.V0, self.m)
        return

    def estimate_Coefs (self, k: float, ro : float, M = 0):
        '''
        Estima os valores dos coeficientes com base em métodos paramétricos presentes no "Methods for estimating stability and control derivatives\
        of conventional subsonic airplanes" de Jan Roskam:\n
            k : ratio of actual average wing section lift curve slope, CLa to 2pi
            ro : densidade do ar (kg/m^3)
            M : numero de mach
        '''
        self.estimate_CLa(k, M)

        self.estimate_CLd(M)
        
        self.estimate_CDa(ro)
        
        return
    
    def get_CL_eq (self):
        '''
        Retorna CL de equilíbrio da aeronave CLe
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
        order = [self.V0, self.Lv, self.Lp, self.Lr, self.Le, self.Lc,\
                self.Nv, self.Np, self.Nr, self.Ne, self.Nc,\
                self.Yv, self.Yp, self.Yr, self.Ye, self.Yc]
        
        if self._many_velocities:
            # make sure that the size stays homogeneous even with float inputs
            deri = np.array(list(map(
                lambda x: x if type(x) == np.ndarray else np.array([x for _ in range(len(self.V0))]),
                order))).transpose()
        else:
            deri = np.array([order])

        return DataFrame(np.round(deri, 4), columns= name)
    
    def curve (self, phi):
        '''
        Retorna o tempo para completar uma curva e a taxa de curvatura:\n
            phi : angulo de rolagem da aeronave (deg)
        '''
        t_c = 2*np.pi*self.V0/(9.81*np.tan(np.radians(phi)))
        # equação 2.25 do COOK

        turn_rate = 2*np.pi/t_c
        # equação 2.26 do COOK

        return t_c, turn_rate
    
    def new_Fin (self, new_fin: Fin, new_lf = None, new_Lf = None, new_hf = None):
        '''
        Troca a EV da aeronave, os dados não modificados não serão alterados:\n
            new_fin : novo objeto Fin
            new_lf : nova distância do CA da asa até o CA da EV no eixo x (m)
            new_Lf : nova distância do CG até o CA da EV no eixo x (m)
            new_hf : nova altura do eixo x até o CA da EV no eixo z (m)
        '''

        lf = new_lf if new_lf != None else self.lf
        Lf = new_Lf if new_Lf != None else self.Lf
        hf = new_hf if new_hf != None else self.hf

        self.f = new_fin
        self.set_fin(lf, Lf, hf)

        return
#=======================================================================================================
if __name__ == "__main__":
    import main