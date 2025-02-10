import numpy as np

def Cld_t (t_c, cf_c):
    '''
    Cld teórico (por radiano) obtido do gráfico da figura B.2,1(a) do Etkins
        t_c : razão da espessura pela corda do perfil aerodinâmico
        cf_c : razão entre a corda da superfície de controle e a corda do perfil em que está
    '''

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

    return np.polyval(np.array([a, b, c, d]), cf_c)

def Cld_Cld_t (cf_c, Cla_Cla_t):
    '''
    Razão Cld corrigido por Cld teórico obtido do gráfico da figura B.2,1(b) do Etkins
        cf_c : razão entre a corda da superfície de controle e a corda do perfil em que está
        Cla_Cla_t : razão Cla corrigido por Cla teórico do perfil
    '''

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

    a = np.polyval(np.array([-1106.09364955,  3141.89348546, -2912.71952857,   878.40801028]), cf_c)
    # r² : 0.9910946477006675

    b = np.polyval(np.array([ 1535.09229464, -4323.4329261 ,  3979.33427803, -1192.70645244]), cf_c)
    # r² : 0.9813068617955021

    c = np.polyval(np.array([ -762.49666592,  2132.61954378, -1950.89867755,   581.34060383]), cf_c)
    # r² : 0.9774636821755441

    d = np.polyval(np.array([ 153.60554985, -429.01078602,  390.47437052, -115.09150056]), cf_c)
    # r² : 0.9889530781715157

    e = np.polyval(np.array([ -5.90609821,  16.21850973, -12.37903544,   3.0605189 ]), cf_c)
    # r² : 0.999981611406111

    return np.polyval(np.array([a, b, c, d, e]), Cla_Cla_t)

def Cla_t (t_c):
    '''
    Cla teórico (por radiano) obtido do gráfico da figura B.1,1(b) do Etkins
        t_c : razão da espessura pela corda do perfil aerodinâmico
    '''
    # r² : 0.988574903403187
    return np.polyval(np.array([4.91957246, 6.29415932]), t_c)

def K1 (cf_c, Aw):
    '''
    Fator de corda do flap K1 obtido no gráfico da figura B.2,2 do Etkins
        cf_c : razão entre a corda da superfície de controle e a corda do perfil em que está
        Aw : Aspect Ratio (alongamento) da asa
    '''

    # K1_01 = 12.447125825584799/(2.17970246308371*np.log(12.530365460259038*Aw - 0.6443923831392444))
    # # r² : 0.9993129506627653

    # K1_02 = 9.055286295555177/(1.129418809869344*np.log(94.28126561344699*Aw - 9.507320283506893)) - 0.05871595813845618  
    # # r² : 0.9989776325285039

    # K1_04 = -3.7812801386437993/(-6461.822979877319*np.log(0.0007576742238701205*Aw + 1.0006737613491827)) + 0.9831439334251642
    # # r² : 0.9994912385691715

    # K1_06 = -21001.919110059305/(-32458.11604649426*np.log(3.6280251154866745*Aw + 2.3404851955535007)) + 0.8427593803422163
    # # r² : 0.9990433000732162

    # razão de mudança do ângulo de zero-lift com a deflexão do flap
    ag_Cl = -0.32121329*np.log(22.64310540900852*cf_c + 1.0601873765142054)
    # r² : 0.999171266752709

    a = np.polyval(np.array([-524430.76867469, 367000.65661605, -73423.96157289, 4209.26748563]), ag_Cl)
    # r² : 1.0

    b = np.polyval(np.array([-2.72971789e+05, 8.33993884e+04, -5.92229414e+03, 3.33870209e+01]), ag_Cl)
    # r² : 1.0

    c = np.polyval(np.array([ 11040.4379927, -12024.67839898, 3652.08186174, -243.47147472]), ag_Cl)
    # r² : 1.0

    d = np.polyval(np.array([-1170.33289345, 1289.79718951, -393.64513332, 26.99248195]), ag_Cl)
    # r² : 1.0

    e = np.polyval(np.array([-68.19917138, 67.06095009, -15.93150261, 0.99073993]), ag_Cl)
    # r² : 1.0

    return a/(b*np.log(c*Aw + d)) + e

def K2 (y1, y2, b, lbd):
    '''
    Fator de envergadura do flap K2 obtido na figura B.2,3 do Etkins
        y1 : posição inical da superfície de controle
        y2 : posição final da suerpfície de controle
        b : envergadura da asa
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

    eta1 = 2*y1/b   # eta inicial
    eta2 = 2*y2/b   # eta final

    return np.polyval(np.array([a, b, c, d]), eta2) - np.polyval(np.array([a, b, c, d]), eta1)

if __name__ == "__main__":
    pass