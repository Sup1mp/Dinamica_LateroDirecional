import numpy as np
from pandas import DataFrame

def remove_space (var):
    '''
    Remove spaces and return only the data on a string
    '''
    return [var[j] for j in range(len(var)) if var[j] != '']

def mach(V, T):
    '''
    Return the mach number for a certain temperature
        V : velocity (m/s)
        T : temperature (°C)
    '''
    return V/np.sqrt(1.4*8.31*(T + 273.15)/0.02897)

def erro (real, aprox):
    return np.where(real!=0, np.round(abs((aprox - real)/real), 3), np.nan)

def erro_dataframe (real: DataFrame, aprox: DataFrame):
    return abs((aprox - real)/real.replace(0, 1))

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
    match len(yi.shape):
        case 1:
            return (b - a)/(2*n) * (yi[0] + yi[-1] + 2*np.sum(yi[1:-1]))
        case 2:
            return (b - a)/(2*n) * (yi[:,0] + yi[:,-1] + 2*np.sum(yi[:,1:-1], axis=1))

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

def read_polar (filename):
    '''
    Coleta os dados dos arquivos de polar do xflr5
    '''
    with open(filename, "r") as file:
        raw = file.read().split('\n')   # reads file

    polar = {}  # dados relacionados à polar

    for i in range(len(raw)):

        var = raw[i].split(' ')     # filter for "words/numbers"

        # coleta as variaveis assossiadas
        if "alpha" in raw[i]:
            polar['polar'] = DataFrame(get_xflr5_table(raw, i), columns=remove_space(var))

        # coleta a velocidade
        if "speed" in raw[i]:
            polar["V0"] = float(var[-2])
    return polar

def read_span (filename):
    '''
    coleta dados do arquivo de asa do xflr5
    '''
    with open(filename, "r") as file:
        raw = file.read().split('\n')   # reads file
    
    span = {}   # dados ao longo da envergadura

    for i in range(len(raw)):
        var = raw[i].split(' ')

        if '=' in raw[i]:
            # coleta avulsos
            var = remove_space(var)
            for j in range(len(var)):
                if var[j] == '=':
                    try:
                        span[var[j - 1]] = float(var[j + 1])
                    except:
                        span[var[j - 1]] = float(var[j + 1][:-2])

        if "y-span" in raw[i]:
            # coleta os dados ao longo da envergadura
            span['span'] = DataFrame(get_xflr5_table(raw, i), columns=remove_space(var))
        
        if "Panel" in raw[i]:
            # coleta dados do centro do pressão em cada painel
            span['Cp'] = DataFrame(get_xflr5_table(raw, i), columns=remove_space(var)).iloc[:, 1:]
    
    return span

if __name__ == "__main__":
    pass