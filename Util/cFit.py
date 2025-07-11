import numpy as np
from scipy.interpolate import interpn
import pandas as pd

def ag_Cl (cf_c):
    '''
    razão de mudança do ângulo de zero-lift com a deflexão do flap\n
        cf_c : razão entre a corda da superfície de controle e a corda do perfil em que está
    '''
    # r² : 0.999171266752709
    a = -0.32121329*np.log(22.64310540900852*cf_c + 1.0601873765142054)
    return np.where(a >= -1, a, -1)

def Cld (t_c, cf_c, Cla_Cla_t):
    '''
    Razão Cld corrigido por Cld teórico (por radianos) obtido do gráfico da figura B.2,1 do Etkins\n
        t_c : razão da espessura pela corda do perfil aerodinâmico
        cf_c : razão entre a corda da superfície de controle e a corda do perfil em que está
        Cla_Cla_t : razão Cla corrigido por Cla teórico do perfil
    '''
    if type(cf_c) != np.ndarray and type(Cla_Cla_t) != np.ndarray:
        cf_c = max(0.1, min(cf_c, 0.5))
        Cla_Cla_t = max(0.7, min(Cla_Cla_t, 0.98))

    # Cld_098 = np.array([ 0.31489457, -0.26165478, -0.08271601,  0.13248927,  0.94637163])
    # # r² : 0.9752489305008382

    # Cld_09 = np.array([-4.14329494,  5.22107478, -2.60883147,  0.76749809,  0.75179737])
    # # r² : 0.9962246748286544

    # Cld_084 = np.array([-7.44014254, 10.11408234, -5.05865866,  1.31736998,  0.60376918])
    # # r² : 0.9963943169798939

    # Cld_078 = np.array([-6.57936124,  8.70387513, -4.42133854,  1.31452314,  0.47042494])
    # # r² : 0.9981353589984262

    # Cld_07 = np.array([-0.41396058,  0.98766589, -0.89591411,  0.72072819,  0.31629528])
    # # r² : 0.9965485499970037

    a = np.polyval(np.array([-1106.09364955,  3141.89348546, -2912.71952857,   878.40801028]), Cla_Cla_t)
    # r² : 0.9910946477006675

    b = np.polyval(np.array([ 1535.09229464, -4323.4329261 ,  3979.33427803, -1192.70645244]), Cla_Cla_t)
    # r² : 0.9813068617955021

    c = np.polyval(np.array([ -762.49666592,  2132.61954378, -1950.89867755,   581.34060383]), Cla_Cla_t)
    # r² : 0.9774636821755441

    d = np.polyval(np.array([ 153.60554985, -429.01078602,  390.47437052, -115.09150056]), Cla_Cla_t)
    # r² : 0.9889530781715157

    e = np.polyval(np.array([ -5.90609821,  16.21850973, -12.37903544,   3.0605189 ]), Cla_Cla_t)
    # r² : 0.999981611406111

    # Cld corrigido
    return np.polyval(np.array([a, b, c, d, e]), cf_c) * Cld_t(t_c, cf_c)

def Cld_t(t_c, cf_c):
    '''
    Cld teórico (por radianos) obtido do gráfico da figura B.2,1 do Etkins\n
        t_c : razão da espessura pela corda do perfil aerodinâmico
        cf_c : razão entre a corda da superfície de controle e a corda do perfil em que está
    '''
    if type(cf_c) != np.ndarray and type(t_c) != np.ndarray:
        cf_c = max(0.1, min(cf_c, 0.5))
        t_c = max(0, min(t_c, 0.15))

    # Cldt_015 = np.array([ 13.78624775, -21.56185942,  17.18739997,   1.02728285])
    # # r² : 0.9946907350538143

    # Cldt_01 = np.array([ 0.63493929, -7.2730938 , 11.65382782,  1.58511417])
    # # r² : 0.9993434000757871

    # Cldt_006 = np.array([ 19.29148536, -28.13476663,  18.81863419,   0.6930432 ])
    # # r² : 0.9974711743165151

    # Cldt_0 = np.array([  8.84566959, -16.55524344,  13.92124889,   1.21889171])
    # # r² : 0.9835579723985991

    a = np.polyval(np.array([ 9.67332847e+04, -2.18824314e+04,  1.13880299e+03,  8.84566959e+00]), t_c)
    # r² : 1.0

    b = np.polyval(np.array([-1.07436861e+05,  2.43352365e+04, -1.26633354e+03, -1.65552434e+01]), t_c)
    # r² : 1.0

    c = np.polyval(np.array([ 3.88489278e+04, -8.82326092e+03,  4.71162604e+02,  1.39212489e+01]), t_c)
    # r² : 1.0

    d = np.polyval(np.array([-4.54946112e+03,  1.03857294e+03, -5.47004582e+01,  1.21889171e+00]), t_c)
    # r² : 1.0

    # Cld Teórico
    return np.polyval(np.array([a, b, c, d]), cf_c)

def Cla_t (t_c):
    '''
    Cla teórico (por radiano) obtido do gráfico da figura B.1,1(b) do Etkins
        t_c : razão da espessura pela corda do perfil aerodinâmico
    '''
    # r² : 0.988574903403187
    return np.polyval(np.array([4.91957246, 6.29415932]), t_c)

def K1 (cf_c, Aw):
    '''
    Fator de corda do flap K1 obtido no gráfico da figura B.2,2 do Etkins\n
        cf_c : razão entre a corda da superfície de controle e a corda do perfil em que está
        Aw : Aspect Ratio (alongamento) da asa
    '''
    try:
        from util import get_json
    except:
        from Util.util import get_json
    

    data = pd.DataFrame(get_json('K1_data'))
    ag = np.sort(np.unique(data['ag_Cl']))[::-1]       # valores unicos ag_Cl
    aw = np.unique(data['Aw'])          # valores unicos Aw
    stru = (len(ag), len(aw))           # estrutura
    points = (ag, aw)
    k1 = np.reshape(data['K1'], stru)   # estrutura resposta como os valores

    # força e limita AR maximo para 10 (limite da interpolação)
    return interpn(points, k1, (min(max(-0.9,ag_Cl(cf_c)), -0.1), min(Aw, 10)))[0]

def K2 (y, bw, lbd):
    '''
    Fator de envergadura do flap K2 obtido na figura B.2,3 do Etkins\n
        y : posição da superfície de controle em y
        bw : envergadura da asa
        lbd : afilamento da asa
    '''
    # K2_0 = np.array([-0.30326651, -0.16742802,  1.46001205,  0.02070015])
    # # r² : 0.926855137757906

    # K2_05 = np.array([-0.27384951, -0.07201788,  1.35541055,  0.00868353])
    # # r² : 0.9707734775429683

    # K2_1 = np.array([-3.84493899e-01,  1.68648997e-01,  1.21840169e+00,  3.97510259e-04])
    # # r² : 0.9623390915632927

    a = np.polyval(np.array([-0.28012278,  0.19889539, -0.30326651]), lbd)
    # r² : 1.0

    b = np.polyval(np.array([ 0.29051347,  0.04556354, -0.16742802]), lbd)
    # r² : 1.0

    c = np.polyval(np.array([-0.06481472, -0.17679564,  1.46001205]), lbd)
    # r² : 1.0

    d = np.polyval(np.array([ 0.0074612 , -0.02776384,  0.02070015]), lbd)
    # r² : 1.0

    return np.polyval(np.array([a, b, c, d]), 2*y/bw)

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # ag_Cl

    # cf_c = np.linspace(0, 1)
    # plt.plot(cf_c, ag_Cl(cf_c), 'k')
    # plt.ylabel('($\\alpha_{\delta})_{C_{l}}$')
    # plt.xlabel('$c_{f}/c$')
    
    #==============================================================================================
    # Cld

    # t_c = np.linspace(0, 0.15)
    # cf_c = np.linspace(0.1, 0.5)
    # Cla_Cla_t = np.linspace(0.7, 0.98, 5)
    
    # for i in Cla_Cla_t:
    #     plt.plot(cf_c, Cld(t_c, cf_c, i)/Cld_t(t_c, cf_c))

    # plt.ylabel('($C_{l_{\delta}}/(C_{l_{\delta}})_{teorico}$')
    # plt.xlabel('$c_{f}/c$')
    # plt.legend(['$C_{l_{\\alpha}}/(C_{l_{\\alpha}})_{teorico}$ = '+str(round(i, 2)) for i in Cla_Cla_t])

    #==============================================================================================
    # Cld_t

    # t_c = np.linspace(0, 0.15, 3)
    # cf_c = np.linspace(0.1, 0.5)

    # for t in t_c:
    #     plt.plot(cf_c, Cld_t(t, cf_c))
    
    # plt.ylabel('$(C_{l_{\delta}})_{teorico}$')
    # plt.xlabel('$c_{f}/c$')
    # plt.legend(['t/c = '+str(round(i, 2)) for i in t_c])

    #==============================================================================================
    # Cla_t

    # t_c = np.linspace(0, 0.2)

    # plt.plot(t_c, Cla_t(t_c), 'k')
    # plt.xlabel('$t/c$')
    # plt.ylabel('$(C_{l_{\\alpha}})_{teorico}$')
    #==============================================================================================
    K1

    # cf_c = np.linspace(0.01, 0.78, 8)
    cf_c = np.array([0.02, 0.05, 0.12, 0.3, 0.4, 0.5, 0.78])
    Aw = np.linspace(0.5, 10)
    ag = []
    
    for c in cf_c:
        x = []
        y = []
        for a in Aw:
            x.append(a)
            y.append(K1(c, a))
        plt.plot(x, y)
        ag.append(max(-0.9,ag_Cl(c)))

    plt.ylabel("K1")
    plt.xlabel('$A_{w}$')
    plt.xlim([0, 10])
    plt.ylim([1, 2])
    plt.legend(['$(\\alpha_{\delta})_{C_{l}}$ = '+str(np.round(i, 2)) for i in ag])
    #==============================================================================================
    # K2

    # lbd = np.linspace(0, 1, 5)
    # bw = 4.5
    # y = np.linspace(0, bw/2)
    
    # for l in lbd:
    #     plt.plot(2*y/bw, K2(y, bw, l))
    
    # plt.xlabel("$\eta$")
    # plt.ylabel("K2")
    # plt.legend(['$\eta$ = '+str(round(i, 2)) for i in lbd])

    #==============================================================================================
    plt.grid()
    plt.show()