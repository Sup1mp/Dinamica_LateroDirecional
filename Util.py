import math
import numpy as np

def mach(V, T):
    '''
    Return the mach number for a certain temperature
        V : velocity (m/s)
        T : temperature (°C)
    '''
    return V/math.sqrt(1.4*8.31*(T + 273.15)/0.02897)

def calc_erro (real, aprox):
    return np.round(abs((aprox - real)/real), 3)

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

if __name__ == "__main__":
    pass