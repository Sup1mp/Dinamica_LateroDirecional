import numpy as np

def Cld (t_c, cf_c, Cla_Cla_t):
    '''
    Razão Cld corrigido por Cld teórico (por radianos) obtido do gráfico da figura B.2,1 do Etkins\n
        t_c : razão da espessura pela corda do perfil aerodinâmico\n
        cf_c : razão entre a corda da superfície de controle e a corda do perfil em que está\n
        Cla_Cla_t : razão Cla corrigido por Cla teórico do perfil
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

    # Cld Teórico
    Cld_t = np.polyval(np.array([a, b, c, d]), cf_c)

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

    # Cld corrigido
    return np.polyval(np.array([a, b, c, d, e]), Cla_Cla_t) * Cld_t

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

    # razão de mudança do ângulo de zero-lift com a deflexão do flap
    ag_Cl = -0.32121329*np.log(22.64310540900852*cf_c + 1.0601873765142054)
    # r² : 0.999171266752709

    # aproximação das curvas de K1
    # K1_01 = 0.011888507746894404/(1.5611044248944355*np.log(0.005821683478021238*Aw + 1.81342892061381) - 0.9247049739637537) + 0.9725943999683583
    # r² : 0.9997067471885825

    # K1_02 = 1.3031342670941815/(0.16242847912458203*np.log(0.8855599843840107*Aw - 0.08945724221775199) + 0.7585216968844564) - 0.059124386335223274
    # r² : 0.998977632595784

    # K1_04 = 0.006249156356691067/(5.320730945367886*np.log(0.0013515378229472428*Aw + 0.8899486489969634) + 0.6275424485379149) + 0.9828504663433438
    # r² : 0.9994912322507276

    # K1_06 = -187.76114396563466/(-290.18142990992874*np.log(3.8429762534987417*Aw + 2.4791506797989658) + 16.702016429820777) + 0.8427592723308773
    # r² : 0.9990433000731884

    # K1_07 = 55.145405798083026/(68.5134149431511*np.log(378.9483693658721*Aw + 131.95059480464863) - 241.87871341874964) + 0.8432261325608006
    # r² : 0.998934087667205

    # K1_08 = 1.386825583915899/(4.914383310766874*np.log(1.6986551221860975*Aw + 1.7991618334176516) + 1.899349056200815) + 0.9159200529327468
    # r² : 0.9986224938225269

    # K1_09 = 0.4499060868722573/(0.031959913406766176*np.log(2.331413923458465*Aw + 0.5465037936517184) + 0.6078786494289624) + 0.34631570365359404
    # r² : 0.996164683462444

    a = np.polyval(np.array([-15.91546612,  16.31609497,  -1.40180908]), ag_Cl)
    # r² : 0.9856930415969143

    b = np.polyval(np.array([-140.4479337 ,  185.78799018,  -59.89180115,    5.83285257]), ag_Cl)
    # r² : 1.0

    c = np.polyval(np.array([ 22.16902481, -35.0176765 ,  17.75085422,  -1.441256  ]), ag_Cl)
    # r² : 1.0

    d = np.polyval(np.array([-67.59078383, 106.03066234, -46.10670546,   5.43138363]), ag_Cl)
    # r² : 1.0

    e = np.polyval(np.array([-21.24783627,  23.15839958,  -3.02600008]), ag_Cl)
    # r² : 0.9999941238154317

    f = np.polyval(np.array([-34.39889194,  54.89915532, -24.37901202,   2.89590294]), ag_Cl)
    # r² : 1.0

    # força e limita AR maximo para 10 (limite da interpolação)
    return a/(b*np.log(abs(c*np.where(Aw > 10, 10, Aw) + d)) + e) + f

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
    import matplotlib.pyplot as plt

    # Cld
    # t_c = np.linspace(0, 0.15)
    # cf_c = np.linspace(0.1, 1)
    # Cla_Cla_t = np.linspace(0.7, 0.98)

    # plt.plot(cf_c, Cld(t_c, cf_c, Cla_Cla_t))
    # plt.grid()
    # plt.show()
    #==============================================================================================
    # Cla_t
    # t_c = np.linspace(0, 0.2)

    # plt.plot(t_c, Cla_t(t_c))
    # plt.grid()
    # plt.show()
    #==============================================================================================
    # K1
    cf_c = np.linspace(0, 1)
    Aw = np.linspace(0, 10)
    
    plt.plot(Aw, K1(cf_c, Aw))
    plt.grid()
    plt.show()
    #==============================================================================================
