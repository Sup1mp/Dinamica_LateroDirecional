from scipy import quad
from math import tan
from Aeronave import Wing, Finn

class Derivadas_LateroDirecional:
    def __init__(self, Asa:Wing, Empenagem:Finn):

        self.w = Asa
        self.f = Empenagem
        self.f.vol_cauda(self.w.S, self.w.b)

        return

    def estabilidade (self,
                      ay,
                      cy,
                      ah,
                      ch,
                      CDy,
                      CLy,
                      dCD_day,
                      a1_F,
                      H_F,
                      y_B,
                      lambda_1_4
                      ):

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

        s = self.w.b/2

        self.Yv = (self.S_B*y_B - self.f.S*a1_F)/self.w.S

        self.Lv = -quad(Lv_int1, 0, s)[0]/(self.w.S * s) - 2*self.w.CL*tan(lambda_1_4)*quad(Lv_int2, 0, s)[0]/(self.w.S *s) - a1_F * self.f.Vc * (self.f.h/self.f.l)

        self.Nv = a1_F * self.f.Vc

        self.Yp = -quad(Yp_int, 0, H_F)[0]/(self.w.S * self.w.b)

        self.Lp =  -quad(Lp_int, 0, s)[0]/(2*self.w.S*s**2)

        self.Np = -quad(Np_int, 0, s)[0]/(2*self.w.S*s**2)

        self.Yr = self.f.Vc*a1_F

        self.Lr = quad(Lr_int, 0, s)[0]/(self.w.S*s**2) + a1_F*self.f.Vc*(self.f.h/self.w.b)

        self.Nr = -quad(Nr_int, 0, s)[0]/(self.w.S*s**2) - a1_F*self.f.Vc*(self.f.l/self.w.b)

        return
    
    def controle (self,
                  cy,
                  dCDy_de,
                  a2_A,
                  a2_R,
                  y1,
                  y2
                  ):

        def Le_int (y):
            return cy*y
        
        def Ne_int (y):
            return dCDy_de*cy*y
        
        s = self.w.b/2

        self.Ye = 0

        self.Le = -a2_A * quad(Le_int, y1, y2)[0]/(self.w.S*s)

        self.Ne = quad(Ne_int, y1, y2)[0]/(self.w.S*s)

        self.Yc = self.f.S/a2_R

        self.Lc = self.f.Vc*a2_R*(self.f.h/self.f.l)

        self.Nc = -self.f.Vc*a2_R

        return
    
if __name__ == "__main__":
    import main