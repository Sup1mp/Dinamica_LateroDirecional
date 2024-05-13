from scipy import quad
from math import tan
from Aeronave import Wing, Finn

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

    def estabilidade (self, ay: list, cy: list, ah: list, ch: list, CDy: list, CLy: list, dCD_day: list, y_B: float):
        '''
        Calcula as derivadas adimensionais de estabilidade presentes no apêndice 8 do livro
            ay : derivada de sustentação local na coordenada y
            cy : corda local na coordenada y
            ah : derivada de sustentação local na coordenada h
            ch : corda local na coordenada h
            CDy : coef de arrasto na coordenada y
            CLy : coef de sustentação na coordenada y
            dCD_day :
            y_B : coef de arrasto lateral total
        '''
        def Lv_int1 (y):
            return cy * ay * self.T * y
        
        def Lv_int2 (y):
            return cy * y
        
        def Yp_int (h):
            return ah * ch * h
        
        def Lp_int (y):
            return (ay + CDy)*cy*y**2 
        
        def Np_int (y):
            return (CLy - dCD_day)*cy*y**2
        
        def Lr_int (y):
            return CLy*cy*y**2
        
        def Nr_int (y):
            return CDy*cy*y**2

        s = self.w.b/2  # pra facilitar

        # derivadas de "v"
        self.Yv = (self.S_B*y_B - self.f.S*self.f.CLa)/self.w.S

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