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
  - c : corda da superfície (no geral de controle) ($m%)
  - c12 : lista com cordas na raiz e na ponta respectivamente ($m$)
  - inc : ângulo de incidência da superfície (°)
  - lbd = $\lambda$ : afilamento da superfície
  - M : velocidade de mach da aeronave
  - m : massa ($kg$)
  - mac : corda média aerodinâmica ($m$)
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
Classe básica de superfícies aerodinâmicas (Asa, Empenagem Horizontal e Vertical)

### AeroSurface.set_CL (CL0, CLa)
Define os valores dos coeficentes de sustentação da superfície, se não for chamada coeficientes são definidos como nulos

### AeroSurface.set_CD (CD0, CDa)
Define os valores dos coeficentes de arrasto da superfície, se não for chamada coeficientes são definidos como nulos

### AeroSurface.get_CD (alpha)
Retorna o valor de CD aproximado com base nas equações de Jan Roskan para o ângulo de ataque alpha desejado

### AeroSurface.get_CL (alpha)
Retorna o valor de CL aproximado linearmente para o ângulo de ataque desejado

### AeroSurface.set_angles (T = 0, V_c4 = 0, V_LE = 0, inc = 0)
Define os valores de ângulos importântes para a superfície, são convertidos para radianos

### AeroSurface.estimate_CLa (K = 0.9, M = 0)
Estima o valor de CLa com base nas equações de Jan Roskan e do Ektins

### AeroSurface.estima_CDa (ro, V0, m)
Estima o valor de CDa com base nas equações de Jan Roskan

## ControlSurface (S, c)
Classe básica para as superfíces de comando (Aileron, Leme e Profundor)

### ControlSurface.set_CLd (CLd)
Define o Valor do coeficiente de sustentação em função da deflexão da superfície

### ControlSurface.set_CLa (CLa)
Define o valor do coeficiente de sustentação em função de alpha

### ControlSurface.estimate_Cla (surf, M = 0)
Estima o valor de CLa da superfície de controle com base nas equações do Jan Roskan e do Ektins

### ControlSurface.estimate_CLd (surf, M = 0)
Estima do valor de CLd da superfície de controle com base nas equações do Ektins
