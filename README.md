# Dinâmica Latero-Direcional
Este projeto constitui no meu TCC em Engenharia Mecânica e possui a finalidade de auxiliar no calculo da estabilidade dinâmica latero-direcional dos aeromodelos utilizados na competição AERO SAE Brasil

## Bibliotecas utilizadas
  - Numpy
  - Sympy
  - pandas
  - matplotlib

## Lista de demais variáveis
  - alpha = $\alpha$ : ângulo de ataque (°)
  - AR : alongamento da superfície
  - b : envergadura da suprefície ($m$)
  - c : corda da superfície (no geral de controle) ($m$)
  - c12 : lista com cordas na raiz e na ponta respectivamente ($m$)
  - inc : ângulo de incidência da superfície (°)
  - lbd = $\lambda$ : afilamento da superfície
  - M : velocidade de mach da aeronave
  - m : massa ($kg$)
  - mac : corda média aerodinâmica ($m$)
  - ro = $\rho$ : densidade do ar ($kg/m^3$)
  - S : área da superfície ($m^2$)
  - T = $\Gamma$ : ângulo de diedro da asa (°)
  - th : espessura do perfil aerodinâmico
  - V0 = $V_0$ : velocidade da aeronave ($m/s$)
  - V_c4 = $\Lambda_{c/4}$ : ângulo de enflexamento na posição de 1/4 de corda (°)
  - V_LE = $\Lambda_{LE}$ : ângulo de enflexamento no bordo de ataque (°)

## Lista de Coeficientes
  - CLa = $C_{L_{\alpha}}$ : Coeficiente de sustentação em função do ângulo de ataque $\alpha$
  - CL0 = $C_{L_0}$ : Coeficiente de sustentação em $\alpha = 0$
  - CDa = $C_{D_{\alpha}}$ : Coeficiente de arrasto em função de $\alpha$
  - CD0 = $C_{D_0}$ : Coeficiente de arrasto em $\alpha = 0$
  - CLd = $C_{L_{\delta}}$ : Coeficiente de sustentação em função da deflexão da superfície de comando $\delta$
  - CLe = $C_{L_e}$ : Coeficiente de sustentação de equilíbrio da aeronave

## Referências
  - Cook, M. V., Flight dynamics principles: a linear systems approach to aircraft stability and control, 3rd ed, 2013
  - Etkins B., Reid L. D., Dynamics of Flight: Stabilty and Control, 3rd ed, 1996 
  - Nelson R. C., Flight stability and automatic control, 1989
  - MIL-F-8785C, U.S Military, 1980

# Documentação
## **Aeronave.AeroSurface (S: _float_, b: _float_, mac: _float_, c12: _list_ = None, th: _float_ = None)**
Classe básica de `Superfícies Aerodinâmicas` (Asa, Empenagem Horizontal e Vertical)
|Métodos | Descrição | Inputs|
|---|---|---|
|set_CL | Define os valores dos coeficentes de sustentação da superfície, se não for chamada coeficientes são definidos como nulos| CL0, CLa|
|set_CD | Define os valores dos coeficentes de arrasto da superfície, se não for chamada coeficientes são definidos como nulos | CD0, CDa|
|get_CD | Retorna o valor de CD aproximado com base nas equações de Jan Roskan para o ângulo de ataque alpha desejado | alpha|
|get_CL | Retorna o valor de CL aproximado linearmente para o ângulo de ataque desejado | alpha|
|set_angles| Define os valores de ângulos importântes para a superfície, são convertidos para radianos | T=0, V_c4=0, V_LE=0, inc=0|
|estimate_CLa | Estima o valor de CLa com base nas equações de Jan Roskan e do Ektins | K=0.9, M=0|
|estimate_CDa | Estima o valor de CDa com base nas equações de Jan Roskan | ro, V0, m|


## **Aeronave.Aircraft (wing: _Wing_, fin: _Fin_, tail: _Tail_, body: _Body_, V0: _float/list_)**
Classe que representa a `aeronave` completa
|Métodos | Descrição | Inputs|
|---|---|---|
|derivatives | Calcula as derivadas aerodinâmicas das superefícies | dCL_day, dCD_day, dCL_dah, dCDy_de, CDy, CLy, cy: _list_, ch: _list_|
|set_angles | | alpha: _float/list_ = 0, theta: _float/list_ = 0|
|set_control | | aileron: _Aileron_, elevator: _Elevator_, rudder : _Rudder_|
|set_fin | | lf: _float_, Lf: _float_, hf: _float_|
|set_mass | | ro: _float/list_, mass: _float_, Ix: _float_, Iz: _float_, Ixz: _float_|
|set_tail | | lt: _float_, Lt: _float_, ht: _float_|
|estimate_CLa | | k: _float_, M = 0|
|estimate_CLd | | M = 0|
|estimate_CDa | | ro: _float_|
|estimate_Coefs | | k: _float_, ro : _float_, M = 0|
|get_CL_eq | | |
|get_CLa | | |
|get_CDa | | |
|get_derivatives | | |
|curve | Retorna o tempo para completar uma curva e a taxa de curvatura respectivamente | phi|


## **Aeronave.Aileron (S: _float_, c: _float_, y1: _float_, y2: _float_)**
Classe filha de ControlSurface que representa o `Aileron`, superfície de controle da Asa, responsável pela rolagem


## **Aeronave.Body (Sl: _float_, h: _float_)**
Classe que representa a `Fuselagem`


## **Aeronave.ControlSurface (S: _float_, c: _float_)**
Classe básica para as `Superfíces de Comando` (Aileron, Leme e Profundor)
|Métodos | Descrição | Inputs|
|---|---|---|
|set_CLd | Define o Valor do coeficiente de sustentação em função da deflexão da superfície | CLd|
|set_CLa | Define o valor do coeficiente de sustentação em função de alpha | CLa|
|estimate_Cla | Estima o valor de Cla, coeficiente de sustentação teórico 2D, da superfície de controle com base nas equações do Jan Roskan e do Ektins | surf, M=0|
|estimate_CLd | Estima do valor de CLd da superfície de controle com base nas equações do Ektins | surf, M=0|
|TAU | Retorna o valor de $\tau$ com base no gráfico encontrado no Nelson | S|


## **Aeronave.Elevator (S: _float_, c: _float_)**
Classe filha de ControlSurface que representa o `Profundor`, superfície de controle da EH, responsável pela cabragem


## **Aeronave.Fin (S: _float_, b: _float_, c12: _list_, th: _float_ = None, k: _int_ = 1)**
Classe filha de AeroSurface para representar a `Empenagem Vertical` (EV)
|Métodos | Descrição | Inputs|
|---|---|---|
|effective_AR | Retorna o valor de AR efetivo da EV considderando influência da EH (não funciona ainda) | AR_B_AR, AR_HB_AR, KH|
|chord | Retorna a corda local na coordenada h da EV, considerada que a EV possui formato trapezoidal alinhado no bordo de fuga | h|


## **Aeronave.Rudder (S: _float_, c: _float_)**
Classe filha de ControlSurface que representa o `Leme`, superfície de controle da EV, responsável pela guinada


## **Aeronave.Tail (S: _float_, b: _float_, mac: _float_, c12: _list_ = None, th: _float_ = None)**
Classe filha de AeroSurface para representar a `Empenagem Horizontal` (EH)


## **Aeronave.Wing (S: _float_, b: _float_, mac: _float_, c12: _list_ = None, th: _float_ = None, Cm_CA: _float_ = 0)**
Classe filha de AeroSurface para representar a `Asa`
|Métodos | Descrição | Inputs|
|---|---|---|
|downwash | Retorna o valor do DownWash teórico baseado no Nelson | |

## **Util.util.mach(V0: _float/list_, alt: _float/list_)**





