from math import tan, cos
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


class Derivadas_LateroDirecional_eng:
    def __init__(self):
        '''
        Classe das derivadas dinâmicas latero-direcionais obtidas do COOK et. al (2013)
        '''
        return

    def estabilidade (self, w: Wing, f: Fin, dCL_day, dCD_day, dCL_dah, CDy, CLy, cy: list, ch: list, y_B: float, S_B: float):
        '''
        Calcula as derivadas adimensionais de estabilidade presentes no apêndice 8 do livro
            w : asa 
            f : empenagem vertical
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
            return cy * dCL_day * w.T * y        # Wing with Dihedral
        
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

        f.vol_cauda(w.S, w.b)   # garante que o volume de cauda vertical foi calculado
        s = w.b/2               # pra facilitar

        n1 = len(cy)    # size of wing chord data
        n2 = len(ch)    # size of fin chord data

        y = np.linspace(0, s, n1)
        h = np.linspace(0, f.b, n2)

        # derivadas de "v" (sideslip)==============================================================
        # Side force
        self.Yv = (S_B*y_B - f.S*f.CLa)/w.S

        # Rolling moment
        self.Lv = -trapezoidal(Lv_int1(), 0, s, n1)/(w.S * s)\
                - 2*w.CLe*tan(w.V_c4)*trapezoidal(Lv_int2(), 0, s, n1)/(w.S *s)\
                - f.CLa * f.Vc * (f.h/f.l)

        # Yawing moment
        self.Nv = f.CLa * f.Vc

        # derivadas de "p" (roll rate)=============================================================
        # Side force 
        self.Yp = -trapezoidal(Yp_int(), 0, f.b, n2)/(w.S * w.b)

        # Rolling moment
        self.Lp =  -trapezoidal(Lp_int(), 0, s, n1)/(2*w.S*s**2)

        # Yawing moment
        self.Np = -trapezoidal(Np_int(), 0, s, n1)/(2*w.S*s**2)

        # derivadas de "r" (yaw rate)==============================================================
        # Side force
        self.Yr = f.Vc*f.CLa

        # Rolling moment
        self.Lr = trapezoidal(Lr_int(), 0, s, n1)/(w.S*s**2)\
                + f.CLa*f.Vc*(f.h/w.b)

        # Yawing moment
        self.Nr = -trapezoidal(Nr_int(), 0, s, n1)/(w.S*s**2)\
                - f.CLa*f.Vc*(f.l/w.b)

        return
    
    def controle (self, w: Wing, f: Fin, cy, dCDy_de):
        '''
        Calcula as derivadas adimensionais de controle presentes no apêndice 8 do livro
            w : asa
            f : empenagem vertical
            cy : corda local na coordenada y
            dCDy_de : derivada do coef de arrasto em função de xi (deflexão do aileron) na coordenada y
        '''
        def Le_int ():
            return cy * y
        
        def Ne_int ():
            return dCDy_de * cy * y
        
        s = w.b/2  # pra facilitar

        n1 = len(cy)
        y = np.linspace(w.a.y1, w.a.y2, n1)

        # derivadas de "e" (aileron)===============================================================
        # Side force
        self.Ye = 0

        # Rolling moment
        self.Le = -w.a.CLa*trapezoidal(Le_int(), w.a.y1, w.a.y2, n1)/(w.S*s)
        
        # Yawing moment
        self.Ne = trapezoidal(Ne_int(), w.a.y1, w.a.y2, n1)/(w.S*s)

        # derivadas de "c" (rudder)================================================================
        # Side force
        self.Yc = f.S/f.r.CLa

        # Rolling moment
        self.Lc = f.Vc*f.r.CLa*(f.h/f.l)

        # Yawing moment
        self.Nc = -f.Vc*f.r.CLa

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


# class Derivadas_LateroDirecional_ame:
#     def __init__(self):
#         '''
#         Classe das derivadas dinâmicas latero-direcionais obtidas no Methods for estimating\
#         stability and control derivatives of conventional subsonic airplanes by ROSKAM
#         '''
#         return
    
#     def estabilidade (self, w: Wing, f: Fin):
#         '''
#         Calcula as derivadas adimensionais de estabilidade
#             So : área transversal da fuselagem no ponto Xo (figura 7.2)
#             Zw : distância vertical de 1/4 da raiz da asa até a linha de centro da fuselagem
#             d : altura maxima na intersecção wing-body
#             S_B : área média da fuselagem
#             l_B : comprimento da fuselagem
#         '''
#         f.vol_cauda(w.S, w.b)    # garante que o volume de cauda vertical foi calculado

#         ETA = 0.724 + 3.06*f.S/(w.S*(1+cos(w.V_c4))) + 0.4*Zw/d + 0.009*w.AR
        
#         a = 2*Zw/d
#         Ki = 0.7*a/0.6 if a >= 0 else -1.85*a   # fator de interferência wing-body para derivada Cyb estimado pela figura 7.1
        

#         Cyb_W = -0.0001 * abs(w.T) * 57.3   # wing
#         Cyb_B = -2*Ki * (So/w.S)            # fuselage (body)
#         Cyb_F = -f.CLa*ETA*f.S/w.S        # fin (vertical)

#         self.Cyb = Cyb_W + Cyb_B + Cyb_F    # variation of side force coef with sideslip angle
        
#         Clb_WB = -1.2*np.sqrt(w.AR)*2*Zw*np.sqrt(S_B/0.7854)/(w.b**2)
#         Clb_F = Cyb_F * (f.h*np.cos(alpha_e) - f.L*np.sin(alpha_e))/w.b

#         self.Clb = Clb_WB + Clb_F   # variation of rolling moment coef with sideslip angle

#         Cnb_B = -57.3*Kn*S_B*l_B/(w.S*w.b)
#         Cnb_F = - Cyb_F *(f.L*np.cos(alpha_e) + f.h*np.sin(alpha_e))/w.b

#         self.Cnb = Cnb_B + Cnb_F    # variation of yawing moment coef with sideslip angle

#         Cyp_F = 2*Clb_F
#         self.Cyp = Cyp_F            # variation of side force coef with roll rate


#         return
    
#     def controle (self):

#         return
        

if __name__ == "__main__":
    import main