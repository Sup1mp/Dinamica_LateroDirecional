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
## AeroSurface (S, b, mac, c12 = None, th = None)
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

## Wing (S, b, mac, c12 = None, th = None, Cm_CA = 0)
Classe filha de AeroSurface para representar a `Asa`

|Métodos | Descrição | Inputs|
|---|---|---|
|downwash | Retorna o valor do DownWash teórico baseado no Nelson | |

## Tail (S, b, mac, c12 = None, th = None)
Classe filha de AeroSurface para representar a `Empenagem Horizontal` (EH)

## Fin (S, b, c12, th = None, k = 1)
Classe filha de AeroSurface para representar a `Empenagem Vertical` (EV)

|Métodos | Descrição | Inputs|
|---|---|---|
|effective_AR | Retorna o valor de AR efetivo da EV considderando influência da EH (não funciona ainda) | AR_B_AR, AR_HB_AR, KH|
|chord | Retorna a corda local na coordenada h da EV, considerada que a EV possui formato trapezoidal alinhado no bordo de fuga | h|

## ControlSurface (S, c)
Classe básica para as `Superfíces de Comando` (Aileron, Leme e Profundor)

|Métodos | Descrição | Inputs|
|---|---|---|
|set_CLd | Define o Valor do coeficiente de sustentação em função da deflexão da superfície | CLd|
|set_CLa | Define o valor do coeficiente de sustentação em função de alpha | CLa|
|estimate_Cla | Estima o valor de Cla, coeficiente de sustentação teórico 2D, da superfície de controle com base nas equações do Jan Roskan e do Ektins | surf, M=0|
|estimate_CLd | Estima do valor de CLd da superfície de controle com base nas equações do Ektins | surf, M=0|
|TAU | Retorna o valor de $\tau$ com base no gráfico encontrado no Nelson | S|

## Aileron (S, c, y1, y2)
Classe filha de ControlSurface que representa o `Aileron`, superfície de controle da Asa, responsável pela rolagem

## Elevator (S, c)
Classe filha de ControlSurface que representa o `Profundor`, superfície de controle da EH, responsável pela cabragem

## Rudder (S, c)
Classe filha de ControlSurface que representa o `Leme`, superfície de controle da EV, responsável pela guinada

## Body (Sl, h)
Classe que representa a `Fuselagem`

## Aircraft (wing, fin, tail, body, V0)
Classe que representa a aeronave como um todo
