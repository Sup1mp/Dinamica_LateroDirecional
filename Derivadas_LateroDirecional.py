from scipy import quad
from math import tan
from Aeronave import Wing, Finn

def weddle (a, b, yi):
    '''
    Weddle's Rule para calculo de integral definida entre "a" e "b" com n = 5 com 7 termos (5 a mais)
        yi : valores de f(xi) igualmente espaçados entre "a" e "b" de modo que yi[0] = f(a) e yi[6] = f(b) 
    '''

    # Weddle's Rule, ESDU 85046
    return (3*(b - a)/50)*(yi[0] + 5*yi[1] + yi[2] + 6*yi[3] + yi[4] + 5*yi[5] + yi[6])


class Derivadas_LateroDirecional:
    def __init__(self, Asa: Wing, EV: Finn):
        '''
        Classe das derivadas dinâmicas latero-direcionais obtidas do COOK et. al (2013)
            Asa : asa da aeronave
            EV : empenagem vertical da aeronave
        '''
        self.w = Asa
        self.f = EV

        self.f.vol_cauda(self.w.S, self.w.b)    # garante que o volume de cauda vertical foi calculado

        return

    def estabilidade (self, dCL_day: list, dCD_day: list, dCL_dah: list, CDy: list, CLy: list, cy: list, ch: list, y_B: float, S_B: float):
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
        def Lv_int1 (y):
            return cy * dCL_day * self.T * y
        
        def Lv_int2 (y):
            return cy * y
        
        def Yp_int (h):
            return dCL_dah * ch * h
        
        def Lp_int (y):
            return (dCL_day + CDy)*cy*y**2 
        
        def Np_int (y):
            return (CLy - dCD_day)*cy*y**2
        
        def Lr_int (y):
            return CLy*cy*y**2
        
        def Nr_int (y):
            return CDy*cy*y**2

        s = self.w.b/2  # pra facilitar

        # derivadas de "v"
        self.Yv = (S_B*y_B - self.f.S*self.f.CLa)/self.w.S

        self.Lv = -quad(Lv_int1, 0, s)[0]/(self.w.S * s) - 2*self.w.CL*tan(self.w.V)*quad(Lv_int2, 0, s)[0]/(self.w.S *s) - self.f.CLa * self.f.Vc * (self.f.h/self.f.l)

        self.Nv = self.f.CLa * self.f.Vc

        # derivadas de "p"
        self.Yp = -quad(Yp_int, 0, self.f.b)[0]/(self.w.S * self.w.b)

        self.Lp =  -quad(Lp_int, 0, s)[0]/(2*self.w.S*s**2)

        self.Np = -quad(Np_int, 0, s)[0]/(2*self.w.S*s**2)

        # derivadas de "r"
        self.Yr = self.f.Vc*self.f.CLa

        self.Lr = quad(Lr_int, 0, s)[0]/(self.w.S*s**2) + self.f.CLa*self.f.Vc*(self.f.h/self.w.b)

        self.Nr = -quad(Nr_int, 0, s)[0]/(self.w.S*s**2) - self.f.CLa*self.f.Vc*(self.f.l/self.w.b)

        return
    
    def controle (self, cy, dCDy_de):
        '''
        Calcula as derivadas adimensionais de controle presentes no apêndice 8 do livro
            cy : corda local na coordenada y
            dCDy_de :
        '''
        def Le_int (y):
            return cy*y
        
        def Ne_int (y):
            return dCDy_de*cy*y
        
        s = self.w.b/2  # pra facilitar

        # derivadas de "e"
        self.Ye = 0

        self.Le = -self.w.a.CLa * quad(Le_int, self.w.a.y1, self.w.a.y2)[0]/(self.w.S*s)

        self.Ne = quad(Ne_int, self.w.a.y1, self.w.a.y2)[0]/(self.w.S*s)

        # derivadas de "c"
        self.Yc = self.f.S/self.f.r.CLa

        self.Lc = self.f.Vc*self.f.r.CLa*(self.f.h/self.f.l)

        self.Nc = -self.f.Vc*self.f.r.CLa

        return
    
if __name__ == "__main__":
    import main