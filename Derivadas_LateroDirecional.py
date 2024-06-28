from math import tan
import numpy as np
from Aeronave import Wing, Fin

def weddle (yi, a, b):
    '''
    Weddle's Rule para calculo de integral definida entre "a" e "b" com n = 5 com 7 termos (5 a mais)
        yi : valores de f(xi) igualmente espaçados entre "a" e "b" de modo que yi[0] = f(a) e yi[6] = f(b)
    '''

    # Weddle's Rule, ESDU 85046
    return (3*(b - a)/50)*(yi[0] + 5*yi[1] + yi[2] + 6*yi[3] + yi[4] + 5*yi[5] + yi[6])

def trapezoidal (yi: list, a: float , b: float, n: int):
    '''
    Trapezoidal Rule para calculo de integral definida entre "a" e "b" com precisão de n termos igualmente espaçados
        yi : valores de f(xi) igualmente espaçados entre "a" e "b" de modo que yi[0] = f(a) e yi[-1] = f(b)
    '''
    if len(yi) == n:
        return (b - a)/(2*n) * (yi[0] + yi[-1] + 2*np.sum(yi[1:-1]))
    raise ValueError("Size of list yi not match with n")


class Derivadas_LateroDirecional:
    def __init__(self, Asa: Wing, EV: Fin):
        '''
        Classe das derivadas dinâmicas latero-direcionais obtidas do COOK et. al (2013)
            Asa : asa da aeronave
            EV : empenagem vertical da aeronave
        '''
        self.w = Asa
        self.f = EV

        self.f.vol_cauda(self.w.S, self.w.b)    # garante que o volume de cauda vertical foi calculado

        return

    def estabilidade (self, dCL_day, dCD_day, dCL_dah, CDy, CLy, cy: list, ch: list, y_B: float, S_B: float):
        '''
        Calcula as derivadas adimensionais de estabilidade presentes no apêndice 8 do livro
            dCL_day : derivada do coef de sustentação em função de alpha na coordenada y
            dCD_day : derivada do coef de arrasto em função de alpha na coordenada y
            dCL_dah : derivada do coef de sustentação em função de alpha na coordenada h
            CDy : coef de arrasto local na coordenada y
            CLy : coef de sustentação local na coordenada y
            cy : corda local na coordenada y
            ch : corda local na coordenada h
            y_B : coef de arrasto lateral total
            S_B : área lateral projetada da fuselagem
        '''
        def Lv_int1 ():
            return cy * dCL_day * self.w.T * y        # Wing with Dihedral
        
        def Lv_int2 ():
            return cy * y                           # wing with aft sweep
        
        def Yp_int ():
            return dCL_dah * ch * h                 # fin contribution
        
        def Lp_int ():
            return (dCL_day + CDy)*cy * y**2        # wing contribuiton
        
        def Np_int ():
            return (CLy - dCD_day)*cy * y**2        # wing contribution
        
        def Lr_int ():
            return CLy * cy * y**2                  # wing contribution
        
        def Nr_int ():
            return CDy * cy * y**2                  # wing contribution

        s = self.w.b/2  # pra facilitar

        n1 = len(cy)    # size of wing chord data
        n2 = len(ch)    # size of fin chord data

        y = np.linspace(0, s, n1)
        h = np.linspace(0, self.f.b, n2)

        # derivadas de "v" (sideslip)==============================================================
        # Side force
        self.Yv = (S_B*y_B - self.f.S*self.f.CLa)/self.w.S

        # Rolling moment
        self.Lv = -trapezoidal(Lv_int1(), 0, s, n1)/(self.w.S * s)\
                - 2*self.w.CLe*tan(self.w.V)*trapezoidal(Lv_int2(), 0, s, n1)/(self.w.S *s)\
                - self.f.CLa * self.f.Vc * (self.f.h/self.f.l)

        # Yawing moment
        self.Nv = self.f.CLa * self.f.Vc

        # derivadas de "p" (roll rate)=============================================================
        # Side force 
        self.Yp = -trapezoidal(Yp_int(), 0, self.f.b, n2)/(self.w.S * self.w.b)

        # Rolling moment
        self.Lp =  -trapezoidal(Lp_int(), 0, s, n1)/(2*self.w.S*s**2)

        # Yawing moment
        self.Np = -trapezoidal(Np_int(), 0, s, n1)/(2*self.w.S*s**2)

        # derivadas de "r" (yaw rate)==============================================================
        # Side force
        self.Yr = self.f.Vc*self.f.CLa

        # Rolling moment
        self.Lr = trapezoidal(Lr_int(), 0, s, n1)/(self.w.S*s**2)\
                + self.f.CLa*self.f.Vc*(self.f.h/self.w.b)

        # Yawing moment
        self.Nr = -trapezoidal(Nr_int(), 0, s, n1)/(self.w.S*s**2)\
                - self.f.CLa*self.f.Vc*(self.f.l/self.w.b)

        return
    
    def controle (self, cy, dCDy_de):
        '''
        Calcula as derivadas adimensionais de controle presentes no apêndice 8 do livro
            cy : corda local na coordenada y
            dCDy_de : derivada do coef de arrasto em função de xi (deflexão do aileron) na coordenada y
        '''
        def Le_int ():
            return cy * y
        
        def Ne_int ():
            return dCDy_de * cy * y
        
        s = self.w.b/2  # pra facilitar

        n1 = len(cy)
        y = np.linspace(self.w.a.y1, self.w.a.y2, n1)

        # derivadas de "e" (aileron)===============================================================
        # Side force
        self.Ye = 0

        # Rolling moment
        self.Le = -self.w.a.CLa*trapezoidal(Le_int(), self.w.a.y1, self.w.a.y2, n1)/(self.w.S*s)
        
        # Yawing moment
        self.Ne = trapezoidal(Ne_int(), self.w.a.y1, self.w.a.y2, n1)/(self.w.S*s)

        # derivadas de "c" (rudder)================================================================
        # Side force
        self.Yc = self.f.S/self.f.r.CLa

        # Rolling moment
        self.Lc = self.f.Vc*self.f.r.CLa*(self.f.h/self.f.l)

        # Yawing moment
        self.Nc = -self.f.Vc*self.f.r.CLa

        return
    
    # def estabilidade_american (self, V0, ro, m):

    #     ad1 = ro*V0*self.w.S*self.w.b
    #     ad2 = ro*V0*self.w.S*self.w.b**2

    #     # derivadas de "v" (sideslip)==============================================================
    #     self.Yv = ad1*Cyb/(4*m)

    #     self.Yp = ad1*Cyp/(4*m)

    #     self.Yr = ad1*Cyr/(4*m)

    #     # derivadas de "p" (roll rate)=============================================================      
    #     self.Lv = ad1*Clb/(2*Ix)

    #     self.Lp = ad2*Clp/(4*Ix)

    #     self.Lr = ad2*Clr/(4*Ix)

    #     # derivadas de "r" (yaw rate)==============================================================
    #     self.Nv = ad1*Cnb/(2*Iz)

    #     self.Np = ad2*Cnp/(4*Iz)

    #     self.Nr = ad2*Cnr/(4*Iz)

    #     return
    
    # def controle_american (self, V0, ro, m):

    #     ad = ro*V0**2*self.w.S*self.w.b
        
    #     # derivadas de "e" (aileron)===============================================================
    #     self.Ye = ad*Cyde/(2*m)

    #     self.Le = ad*Clde/(2*Ix)

    #     self.Ne = ad*Cnde/(2*Iz)
       
    #     # derivadas de "c" (rudder)================================================================
    #     self.Yc = ad*Cydr/(2*m)

    #     self.Lc = ad*Cldr/(2*Ix)

    #     self.Nc = ad*Cndr/(2*Iz)

    #     return

    def show (self):
        print(f"\nYv: {self.Yv}\nYp: {self.Yp}\nYr: {self.Yr}\nYe: {self.Ye}\nYc: {self.Yc}\n")
        print(f"Nv: {self.Nv}\nNp: {self.Np}\nNr: {self.Nr}\nNe: {self.Ne}\nNc: {self.Nc}\n")
        print(f"Lv: {self.Lv}\nLp: {self.Lp}\nLr: {self.Lr}\nLe: {self.Le}\nLc: {self.Lc}\n")
    
if __name__ == "__main__":
    import main