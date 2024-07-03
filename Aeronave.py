import math
import numpy as np
from pandas import DataFrame

def Tau (S_super, S_control):
    '''
    Retorna o valor do parâmetros TAU para superfícies de controle
        S_super : área total da empenagem
        S_control : área da superficie de controle (leme ou profundor)
    '''
    x = S_control / S_super

    # equação da figura 2.20 do NELSON
    return 12.267*(x**5) - 25.8*(x**4) + 21.477*(x**3) - 9.6446*(x**2) + 3.26*x + 0.01

def remove_space (var):
    return [var[j] for j in range(len(var)) if var[j] != '']

def get_xflr5_table (raw, index):
    values = []

    for i in range(index+1, len(raw)):
        var = raw[i].split(' ')     # filter for "words/numbers"

        if len(var) > 1:    # checks for empty spaces
            try:
                # try to get float
                values.append([float(var[j]) for j in range(len(var)) if var[j] != ''])
            except:
                # not float, ignores and move on
                pass
        else:
            # empty space = end of table
            break

    return np.array(values)
    
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
        self.AR = (self.b**2) / self.S                  # alongamento
        self.lbd = self.c12[1] / self.c12[0]            # afilamento

        self.polar = {}     # dados relacionados à polar
        self.span = {}      # dados ao longo da envergadura

        self.set_angles(T = 0, V_c4 = 0)   # Ang Diedro = Ang Enflexamento = 0

        return
    
    def read_polar (self, filename):
        '''
        Coleta os dados dos arquivos de polar do xflr5
        '''
        with open(filename, "r") as file:
            raw = file.read().split('\n')   # reads file

        for i in range(len(raw)):

            var = raw[i].split(' ')     # filter for "words/numbers"

            # coleta as variaveis assossiadas
            if "alpha" in raw[i]:
                self.polar['polar'] = DataFrame(get_xflr5_table(raw, i), columns=remove_space(var))

            # coleta a velocidade
            if "speed" in raw[i]:
                self.polar["V0"] = float(var[-2])

        # obtem CL
        fin = self.polar['polar']['CL'].argmax()
        ini = self.polar['polar']['CL'].argmin()

        self.CLa = round((self.polar['polar'].iloc[fin]['CL'] - self.polar['polar'].iloc[ini]['CL']) / (self.polar['polar'].iloc[fin]['alpha'] - self.polar['polar'].iloc[ini]['alpha']), 6)
        self.CL0 = self.polar['polar'][self.polar['polar']['alpha'] == 0]['CL'][0]

        # obtem CD
        fin = self.polar['polar']['CD'].argmax()
        ini = self.polar['polar']['CD'].argmin()
        
        self.CDa = round((self.polar['polar'].iloc[fin]['CD'] - self.polar['polar'].iloc[ini]['CD']) / (self.polar['polar'].iloc[fin]['alpha'] - self.polar['polar'].iloc[ini]['alpha']), 6)
        self.CD0 = self.polar['polar'][self.polar['polar']['alpha'] == 0]['CD'][0]

        return
    
    def read_span (self, filename):
        '''
        coleta dados do arquivo de asa do xflr5
        '''
        with open(filename, "r") as file:
            raw = file.read().split('\n')   # reads file

        for i in range(len(raw)):
            var = raw[i].split(' ')

            if '=' in raw[i]:
                # coleta avulsos
                var = remove_space(var)
                for j in range(len(var)):
                    if var[j] == '=':
                        try:
                            self.span[var[j - 1]] = float(var[j + 1])
                        except:
                            self.span[var[j - 1]] = float(var[j + 1][:-2])

            if "y-span" in raw[i]:
                # coleta os dados ao longo da envergadura
                self.span['span'] = DataFrame(get_xflr5_table(raw, i), columns=remove_space(var))
            
            if "Panel" in raw[i]:
                # coleta dados do centro do pressão em cada painel
                self.span['Cp'] = DataFrame(get_xflr5_table(raw, i), columns=remove_space(var)).iloc[:, 1:]
        
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
        e = 1/(1.05 + 0.007*math.pi*self.AR)    # Obert
        # e = 1.78*(1 - 0.045*self.AR**0.68) - 0.64   # internet
        
        # pode ser deduzida a partir da equação 3.1 do "Methods for estimating stability and control\
        #  derivatives of conventional subsonic airplanes" de Jan Roskam
        return self.CD0 + (self.get_CL(alpha)**2 - self.CL0**2)/(math.pi*self.AR*e)
    
    def get_CL (self, alpha):
        return self.CL0 + self.CLa * (alpha + self.inc)

    def set_angles (self, T: float, V_c4: float):
        '''
        Ângulos:
            T : ang diedro (deg)
            V_c4 : ang enflexamento (deg)
        '''
        self.T = math.radians(T)            # angulo de diedro
        self.V_c4 = math.radians(V_c4)      # angulo de enflexamento

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
    
    def estimate_CLe (self, ro, V, m):
        '''
        Estima o valor do CL de equilibrio com base em métodos paramétricos presentes no "Methods for estimating stability and control derivatives\
        of conventional subsonic airplanes" de Jan Roskam:
            ro : densidade do ar (kg/m^3)
            V : velocidade (m/s)
            m : massa da aeronave (kg)
        '''
        # equação 3.3
        self.CLe = 2*m/(ro*V**2 * self.S)
        
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
class Empennage (AeroSurface):
    def __init__(self, S: float, b: float, mac: float, inc: float, c12: list, l: float, L: float, h: float):
        '''
        Empennage object, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            inc : ângulo de incidência (deg)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            l : distância do CA até o CA da empenagem em x (m)
            L : distância do CG até o CA da empenagem em x (m)
            h : distância do eixo x até o CA da empenagem em z (m)
        '''
        super().__init__(S, b, mac, inc, c12)

        self.l = l
        self.L = L
        self.h = h

        return
#=======================================================================================================
class Rudder (AeroSurface):
    def __init__(self):
        return
    
    def correct_CLa (self, f, AR_f):
        '''
        Correção de CLa do leme devido ao efeito geométrico da EV\
        usando o método descrito em Abbot and Von Doenhoff (1959)
            f : fator de correção empirico
            AR_f : alongamento da EV
        '''
        self.CLa = f*self.CLa/(1 + self.CLa/(math.pi*AR_f))
        # equação 13.254 COOK
        return
#=======================================================================================================
class Elevator (AeroSurface):
    def __init__(self):
        return
#=======================================================================================================
class Aileron (AeroSurface):
    def __init__(self):
        return
    
    def set_geometry (self, b: float, c: float, y1: float, y2: float = None):
        '''
        Define as dimensioes básicas do aileron
            b : envergadura
            c : corda
            y1 : distância perpendicular do eixo x até o começo do aileron
            y2 : distância perpendicular do eixo x até o final do aileron
        '''
        self.b = b
        self.c = c
        self.y1 = y1
        self.y2 = y2 if y2 != None else y1 + b
        
        return
#=======================================================================================================
class Wing (AeroSurface):
    def __init__(self, S: float, b: float, mac: float, inc: float, c12: list = None, Cm_CA: float = 0):
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
        self.a = Aileron()      # aileron

        return
    
    def downwash (self):
        '''
        Calcula a derivada de_da do angulo de downwash por alpha
        '''
        self.de_da = 2*self.CLa / (math.pi * self.AR)
        # equação 2.23 NELSON

        return self.de_da
    
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
    
    def Cn_b (self, kn, krl, SFs, lF):
        '''
        Contribuição da asa para o coeficiente de momento angular de guinada
            kn : coeficiente empirico de interferência asa-fuselagem (pode ser obtido pela figura 2.28 do NELSON)
            krl : fator de correção em função do reynolds da fuselagem
            SFs : área lateral projetada da fuselagem
            lF : comprimento da fuselagem
        '''
        self.Cnb = - kn * krl * SFs * lF /(self.S * self.b)     # per deg
        # equação 2.74 do NELSON

        return self.Cnb
#=======================================================================================================
class Fin(Empennage):
    def __init__(self, b: float, c12: list, l: float, L: float, h: float, S: float = None, k: int = 1):
        '''
        Empenagem Vertical, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            l : distância do CA até o CA da empenagem em x (m)
            L : distância do CG até o CA da empenagem em x (m)
            h : distância do eixo x até o CA da empenagem em z (m)
            k : quantidade de EV's
        '''
        S = S if S!= None else (c12[0] + c12[1])*b*k/2  # area de trapézio

        super().__init__(S, b, (c12[0] + c12[1])/2, 0, c12, l, L, h)
        
        self.V_c4 = math.atan(3*(c12[0] - c12[1])/(4*b))    # alinhamento no bordo de fuga
        self.k = k
        self.r = Rudder()       # leme
        
        return
    
    def vol_cauda (self, Sw: float, bw: float):
        '''
        Calcula o Volume de Cauda Vc da empenagem:
            Sw : área da asa (m)
            bw : envergadura da asa (m)
        '''
        self.Vc = self.L * self.S / (Sw * bw)

        # # volume de cauda modificado para EV que possibilita empiricamente ajuste dinâmico
        # self.Vc = self.S*(self.mac+0.7*self.z*math.tan(self.V))/(Sw*bw)
        # # equação 13.244 do COOK

        return
    
    def effective_AR (self, AR_B_AR, AR_HB_AR_B, KH):
        '''
        Calcula o AR effetivo da EV, obtido no Methods for estimating stability and control derivatives\
         of conventional subsonic airplanes by ROSKAM
            AR_B_AR : razão do AR na preseça do corpo com o AR da EV isolada (figura 7.5)
            AR_HB_AR_B : razão do AR na presença da EH e do corpo com o AR só na presença do corpo (figura 7.6)
            KH : fator que conta com o tamanho relativo da EH e da EV (figura 7.7)
        '''
        self.AR = AR_B_AR * self.AR * (1 + KH* (AR_HB_AR_B - 1))
        # equação 7.6 do Methods for estimating stability and control derivatives of conventional\
        #  subsonic airplanes by ROSKAM
        
        return

    def Cl_b (self):
        '''
        Retorna a contribuição da EV para o coeficiente de momento lateral em função de beta (sideslip)
        Depende do Volume de cauda Vc e do coef CL alpha
        '''
        self.Clb = - self.Vc * self.CLa * self.h / self.l
        # equação 3.39 do COOK

        return self.Clb
    
    def Cn_b (self, asa : Wing, zw: float, d: float):
        '''
        Retorna o coeficiente de momento angular de guinada
            zw : distancia paralela ao eixo z da corda da asa até a linha de centro da fuselagem
            d : profundidade da fuselagem
        '''
                
        self.Cnb = self.Vc * self.CLa * (0.724 + 3.06*(self.S/(asa.S*(1 + math.cos(asa.V_c4)))) + 0.4*zw/d + 0.009 * asa.AR)
        # equação 2.80 + 2.81 do NELSON

        return self.Cnb
    
    def Cn_dr (self, ETA):
        '''
        Retorna o coeficiente de momento direcional do leme em função da deflexão delta
            ETA : eficiencia de cauda
            S_rudder : área total do leme
        '''
        self.Cndr = - self.k * ETA * self.Vc * self.CLa * Tau(self.S, self.r.S)

        return self.Cndr
#=======================================================================================================
class Tail (Empennage):
    def __init__(self, S: float, b: float, mac: float, inc: float, c12: list, l: float, L: float, h: float):
        '''
        Empenagem Horizontal, parâmetros geométricos:
            S : área (m^2)
            b : envergadura (m)
            mac : corda média aerodinamica (m)
            inc : ângulo de incidência (deg)
            c12 : cordas na raiz e na ponta [raiz, ponta]
            l : distância do CA até o CA da empenagem em x (m)
            L : distância do CG até o CA da empenagem em x (m)
            h : distância do eixo x até o CA da empenagem em z (m)
        '''
        super().__init__(S, b, mac, inc, c12, l, L, h)

        self.e = Elevator()     # profundor

        return
    
    def vol_cauda (self, Sw: float, MACw: float):
        '''
        Calcula o Volume de Cauda Vc da empenagem:
            Sw : área da asa (m)
            MACw : corda media aerodinamica da asa (m)
        '''
        self.Vc = self.L * self.S / (Sw * MACw)

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
    
    def Cm_de (self, ETA):
        '''
        Retorna o coef de momento do profundor em função da deflexão delta
        '''
        self.Cmde = - ETA * Tau(self.S, self.e.S) * self.Vc * self.CLa

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
            V0 : velocidade da aeronave (m/s)
            theta_e : ângulo de equilíbrio entre horizonte e direção de voo (deg)
        '''
        self.m = m
        self.Ix = Ix
        self.Ixz = Ixz
        self.Iz = Iz

        self.V0 = V0
        self.theta_e = np.radians(theta_e)

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
    
    def curve (self, phi):
        '''
        Retorna o tempo para completar uma curva e a taxa de curvatura
            phi : angulo de rolagem da aeronave
        '''
        self.t_c = 2*math.pi*self.V0/(9.81*math.tan(phi))
        # equação 2.25 do COOK

        self.turn_rate = 2*math.pi/self.t_c
        # equação 2.26 do COOK

        return self.t_c, self.turn_rate
#=======================================================================================================
if __name__ == "__main__":
    import main