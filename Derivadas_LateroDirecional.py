from scipy import quad
from math import tan

class Derivadas_LateroDirecional:
    def __init__(
            self,
            ro,
            m,
            Ix,
            Iz,
            Ixz,
            V0,
            S,
            b,
        ):
        
        self.S = S
        self.b = b
        self.V_F = V_F

        # adimensionalizadores
        ad1 = 0.5*ro*V0*S
        ad2 = ad1*b

        self.m = m/ad1                      # massa adimensional
        self.Ix = Ix/ad2                    # momento de inercia Ix adimensional
        self.Iz = Iz/ad2                    # momento de inercia Iz adimensional
        self.Ixz = Ixz/ad2                  # momento de inercia Ixz adimensional

        return

    # def estabilidade (self):

    #     def wing_dihedral (y):
    #         return cy * ay * T * y
        
    #     def wing_sweep (y):
    #         return cy * y
        
    #     def finn_roll (h):
    #         return ah * ch * h

    #     s = self.b/2

    #     self.Yv = (self.S_B*y_B - self.S_F*a1_F)/self.S

    #     self.Lv = -quad(wing_dihedral, 0, s)[0]/(self.S * s) - 2*self.C_L*tan(lambda_1_4)*quad(wing_sweep, 0, s)/(self.S *s) - a1_F * V_F * (h_F/l_F)

    #     self.Nv = a1_F * V_F

    #     self.Yp = -quad(finn_roll, 0, H_F)/(self.S * self.b)