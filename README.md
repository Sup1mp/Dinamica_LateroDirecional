## Dinâmica Latero-Direcional
Este projeto constitui no meu TCC em Engenharia Mecânica e possui a finalidade de auxiliar no calculo da estabilidade dinâmica latero-direcional dos aeromodelos utilizados na competição AERO SAE Brasil

## Bibliotecas utilizadas
  - Numpy
  - Sympy
  - pandas
  - matplotlib

## Lista de demais variáveis
  - S : área da superfície
  - b : envergadura da suprefície
  - c : corda da superfície (no geral de controle)
  - mac : corda média aerodinâmica
  - AR : alongamento da superfície
  - lbd = $\lambda$ : afilamento da superfície

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
