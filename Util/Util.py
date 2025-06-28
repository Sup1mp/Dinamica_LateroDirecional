import numpy as np
from pandas import DataFrame
import json

def remove_space (var):
    '''
    Remove spaces and return only the data on a string
    '''
    return [var[j] for j in range(len(var)) if var[j] != '']

def mach(V0, alt):
    '''
    Return the mach number for a certain temperature
        V0 : velocity (m/s)
        alt : altitude (m)
    '''
    # return V0/np.sqrt(1.4*8.31*(T + 273.15)/0.02897)
    return V0/getAtmosphere(alt)[3]

def erro (real, aprox):
    return np.round(np.where(real!=0, abs((aprox - real)/real), abs((aprox - real))), 3)

def erro_dataframe (real: DataFrame, aprox: DataFrame):
    return np.round(abs((aprox - real)/real.replace(0, 1)), 3)

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

def compara_derivadas (real, calc):
    '''
    Compara o erro entre as derivadas
    '''
    err = erro_dataframe(real, calc)

    print(f"Referência:\n{real}")
    print(f"Calculado:\n{calc}")
    print(f"Erro Percentual:\n{err}")
    
    st = np.round(err.describe(), 4)
    print(f"Média Erros:\n{st.loc['mean']}\nMédia Total: {st.loc['mean'].sum()/len(st.loc['mean'])}")
    return

def getAtmosphere (altitude):
    '''
    Retorna dados de Temperatura, Pressão, Densidade e Velocidade do Som respectivamente para uma data altitude:\n
        altitude : altitude (m)
    OBS: Limite de altitude de 18900 m
    '''
    # from scipy.interpolate import interp

    # pega dados atmosféricos
    data = DataFrame(get_json("Atmospheric_data"))
    alt = data['Altitude(m)']

    temp = np.interp(altitude, alt, data['Temperature(K)'])         # interpolação linear
    pres = np.interp(altitude, alt, data['Pressure(N/m^2)'])        # interpolação linear
    dens = np.interp(altitude, alt, data['Density(kg/m^3)'])        # interpolação linear
    sound = np.interp(altitude, alt, data['Speed of Sound(m/s)'])   # interpolação linear

    return temp, pres, dens, sound

def get_json (filename):
    '''
    Retorna dados de json do caminho "Util/Data/filename"
    '''
    with open(f"Util/Data/{filename}.json", 'r') as file:
        data = json.load(file)  # lê dados do arquivo
    return data

def save_json(data, filename):
    '''
    Salva dados em json no caminho "Util/Data/filename"
    '''
    with open(f"Util/Data/{filename}.json", '+w') as file:
        json.dump(data, file)   # salva dados no arquivo
    return

def oswald (lbd, AR, T, V_c4, M):

    # contabildade do efeito de couple entre lbd e V_c4
    d_lbd = -0.357 + 0.45*np.exp(-0.0375*V_c4)

    # fit de curva por Nita e Scholz
    f_lbd = np.polyval(np.array([0.0524, -0.15, 0.1659, -0.0706, 0.0119]), lbd-d_lbd)

    # fator de correção para velocidade de Mach > 0.3 
    k_e = np.where(M > 0.3, -0.001521*((M/0.3 - 1)**10.82) + 1, 1)

    # fator de correção para diedro
    k_T = (1 + (1/np.cos(T) - 1)/2.83)**2

    # Equação do fator de Oswald dado por Hörner e trabalhada por Nita e Scholz
    return (1/(1+f_lbd*AR)) * 0.114  * k_T * k_e
    

if __name__ == "__main__":
    # # erro test
    # dn_1 = np.random.random((10))
    # dn_2 = np.array([i for i in range(0, 10)])
    # dn1 = 12.5
    # dn2 = 10
    # print(f"erro_1: {erro(dn_2, dn_1)}")
    # print(f"erro1: {erro(dn2, dn1)}")

    # # dataframe error test
    # ddf_1 = DataFrame(dn_1)
    # ddf_2 = DataFrame(dn_2)
    # print(f"erro_dataframe: {erro_dataframe(ddf_2, ddf_1)}")

    # t, p, d, s = getAtmosphere([19000, 200])
    # print(f"temp: {t} K \npress: {p} N/m^2\ndens: {d} kg/m^3\nsound: {s} m/s")

    # V0 = 200
    # alt = 5000
    # print(f"Para {V0} m/s a {alt} m de altitude:\nmach: {mach(V0, alt)}")
    print(np.exp(1))