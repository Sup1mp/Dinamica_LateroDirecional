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
  - b : envergadura da superfície ($m$)
  - c : corda da superfície (no geral de controle) ($m$)
  - c12 : lista com cordas na raiz e na ponta respectivamente ($m$)
  - inc : ângulo de incidência da superfície (°)
  - lf : distância, no eixo x, entre os 𝐶_𝐴’s da asa e da EV [m]
  - Lf : distância, no eixo x, entre o 𝐶_𝐺 e o 𝐶_𝐴 da EV [m]
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
  - Vv : volume de cauda vertical
  - Vh : volume de cauda horizontal
  - y1 e y2 : distância, no eixo y, entre a linha de centro da fuselagem e o começo e final do aileron respectivamente


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
## `Aeronave.AeroSurface (S: float, b: float, mac: float, c12: list = None, th: float = None)`
**Descrição:**
Classe base para as _Superfícies Aerodinâmicas_ (Asa, Empenagem Horizontal e Vertical) que possui diversos métodos comuns entre elas.

**Inputs:**
  - S: Área da superfície $[m]$
  - b: Envergadura da superfície $[m]$
  - mac: Corda Média aerodinâmica da superfície $[m]$
  - c12: Lista com as cordas na Raiz e na Ponta da superfície, respectivamente. Caso não seja especificado, possuirá os valores de `mac` $[m]$
  - th: Espessura percentual do perfil aerodinãmico


### `Aeronave.AeroSurface.estimate_CDa (ro, V0, m)`
**Descrição:**
Estima o valor de CDa com base nas equações de Jan Roskan

**Inputs:**
  - ro: Densidade do Ar $[kg/m^3]$
  - V0: Velocidade da aeronave $[m/s]$
  - m: Massa da aeronave  $[kg]$

**Outputs:**
  - CDa: Coeficiente de Arrasto em função de $\alpha$ estimado

  
### `Aeronave.AeroSurface.estimate_CLa (K = 0.9, M = 0)`
**Descrição:**
Estima o valor de CLa com base nas equações de Jan Roskan e do Ektins

**Inputs:**
  - K:
  - M: Número de Mach da aeronave, para voos com velocidades muito baixas pode ser aproximado para 0

**Outputs:**
  - CLa: Coeficiente de Sustentação em função de $\alpha$ estimado


### `Aeronave.AeroSurface.get_CD (alpha)`
**Descrição:**
Retorna o valor de CD aproximado com base nas equações de Jan Roskan para o ângulo de ataque alpha desejado

**Inputs:**
  - alpha: Ângulo de ataque desejado da superfície $[deg]$

**Outputs:**
  - CD: Coeficiente de Arrasto para o alpha desejado


### `Aeronave.AeroSurface.get_CL (alpha)`
**Descrição:**
Retorna o valor de CL aproximado linearmente para o ângulo de ataque desejado

**Inputs:**
  - alpha: Ângulo de ataque desejado da superfície $[deg]$

**Outputs:**
  - CD: Coeficiente de Arrasto para o alpha desejado


### `Aeronave.AeroSurface.get_angles (T = 0, V_c4 = 0, V_LE = 0, inc = 0)`
**Descrição:**
Define os valores de ângulos importântes para a superfície, são convertidos para radianos

**Inputs:**
  - T: Ângulo de Diedro da superfície $[deg]$
  - V_c4: Ângulo de Enflexamento em 1/4 da corda da superfície $[deg]$
  - V_LE: Ângulo de Enflexamento no bordo de ataque da superfície $[deg]$
  - inc: Ângulo de incidÊncia da superfície $[deg]$


### `Aeronave.AeroSurface.set_CD (CD0, CDa)`
**Descrição:**
Define os valores dos coeficentes de arrasto da superfície, se não for chamada coeficientes são definidos como nulos.

**Inputs:**
  - CD0: Coeficiente de Arrasto para ângulo de ataque igual a 0 ($\alpha = 0$)
  - CDa: Coeficiente de Arrasto em função do ângulo de ataque ($\alpha$) $[1/deg]$


### `Aeronave.AeroSurface.set_CL (CL0 = None, CLa = None)`
**Descrição:**
Define os valores dos coeficentes de sustentação da superfície, se não for chamada coeficientes são definidos como nulos.

**Inputs:**
  - CL0: Coeficiente de Sustentação para ângulo de ataque igual a 0 ($\alpha = 0$)
  - CLa: Coeficiente de Sustentação em função do ângulo de ataque ($\alpha$) $[1/deg]$


## `Aeronave.Aircraft (wing: Wing, fin: Fin, tail: Tail, body: Body, V0: float | list)`
**Descrição:**
Classe que representa a `aeronave` completa

### `Aeronave.Aircraft.curve (phi)`
**Descrição:**
Retorna o tempo para completar uma curva e a taxa de curvatura respectivamente

### `Aeronave.Aircraft.derivatives (dCL_day, dCD_day, dCL_dah, dCDy_de, CDy, CLy, cy: list, ch: list)`
**Descrição:**
Calcula as derivadas aerodinâmicas da aeronave

### `Aeronave.Aircraft.estimate_CDa (ro: float)`

### `Aeronave.Aircraft.estimate_CLa (k: float, M = 0)`

### `Aeronave.Aircraft.estimate_Cld (M = 0)`

### `Aeronave.Aircraft.estimate_Coefs (k: float, ro : float, M = 0)`

### `Aeronave.Aircraft.get_CDa (ro: float)`

### `Aeronave.Aircraft.get_CL_eq ()`

### `Aeronave.Aircraft.get_CLa ()`

### `Aeronave.Aircraft.get_derivatives ()`

### `Aeronave.Aircraft.set_angles (alpha: float | list = 0, theta: float | list = 0)`

### `Aeronave.Aircraft.set_control (aileron: Aileron, elevator: Elevator, rudder : Rudder)`

### `Aeronave.Aircraft.set_fin (lf: float, Lf: float, hf: float)`

### `Aeronave.Aircraft.set_mass (ro: float | list, mass: float, Ix: float, Iz: float, Ixz: float)`

### `Aeronave.Aircraft.set_tail (lt: float, Lt: float, ht: float)`


## `Aeronave.Aileron (S: float, c: float, y1: float, y2: float)`
**Descrição:**
Classe filha de `ControlSurface` que representa o _Aileron_, superfície de controle da Asa, responsável pela rolagem

**Inputs:**
  - S: Área do aileron $[m^2]$
  - c: Cordar do aileron $[m]$
  - y1: Distância, no eixo y, da primeira corda do aileron até a linha de centro da aeronave $[m]$
  - y1: Distância, no eixo y, da última corda do aileron até a linha de centro da aeronave $[m]$


## `Aeronave.Body (Sl: float, h: float)`
**Descrição:**
Classe que representa a _Fuselagem_ da aeronave.

**Inputs:**
  - Sl: Área lateral da fuselagem $[m^2]$
  - h: Altura da fuselagem $[m]$


## `Aeronave.ControlSurface (S: float, c: float)`
**Descrição:**
Classe básica para as _Superfíces de Comando_ (Aileron, Leme e Profundor)

**Inputs:**
  - S: Área da superfície $[m^2]$
  - b: Envergadura da superfície $[m]$

### `Aeronave.ControlSurface.estimate_CLd (surf: Fin | Tail | Wing, M: float = 0)`
**Descrição:**
Estima do valor de CLd da superfície de controle com base nas equações do Ektins

**Inputs:**
  - surf: Superfície à qual o controle se refere
  - M: Número de MACH

**Outputs:**
  - CLd: Coeficiente de sustentação em função do ângulo de deflexão ($\delta$) estimado $[1/deg]$

### `Aeronave.ControlSurface.estimate_Cla (surf: Fin | Tail | Wing, M: float = 0)`
**Descrição:**
Estima o valor de Cla, coeficiente de sustentação teórico 2D, da superfície de controle com base nas equações do Jan Roskan e do Ektins

**Inputs:**
  - surf: Superfície à qual o controle se refere
  - M: Número de MACH

**Outputs:**
  - Cla: Coeficiente de sustentação teórico 2D estimado $[1/deg]$

### `Aeronave.ControlSurface.set_CLa (CLa)`
**Descrição:**
Define o valor do coeficiente de sustentação em função de alpha

**Inputs:**
  - CLa: Coeficiente de Sustentação em função de $\alpha$ desejado $[1/deg]$

### `Aeronave.ControlSurface.set_CLd (CLd)`
**Descrição:**
Define o Valor do coeficiente de sustentação em função da deflexão da superfície

**Inputs:**
  - CLd: Coeficiente de Sustentação em função da deflexão $\delta$ desejado $[1/deg]$

### `Aeronave.ControlSurface.TAU (S: float)`
**Descrição:**
Retorna o valor de $\tau$ com base no gráfico encontrado no Nelson

**Inputs:**
  - S: Área da superfície referente ao controle.
  

## `Aeronave.Elevator (S: float, c: float)`
Classe filha de `ControlSurface` que representa o _Profundor_, superfície de controle da EH, responsável pela cabragem


## `Aeronave.Fin (S: float, b: float, c12: list, th: float = None, k: int = 1)`
Classe filha de `AeroSurface` para representar a _Empenagem Vertical_ (EV)

**Inputs:**
  - k: Número de EV's presentes na aeronave


### `Aeronave.Fin.chord (h: float)`
**Descrição:**
Retorna a corda local na coordenada h da EV, considerada que a EV possui formato trapezoidal alinhado no bordo de fuga

**Inputs:**
  - h: Altura da corda local ao longo da envergadura, ou seja, _h_ varia de 0 a _b_ $[m]$


### `Aeronave.Fin.effective_AR (AR_B_AR, AR_HB_AR, KH)`
**Descrição:**
Retorna o valor de AR efetivo da EV considderando influência da EH (não funciona ainda)


## `Aeronave.Rudder (S: float, c: float)`
Classe filha de `ControlSurface` que representa o _Leme_, superfície de controle da EV, responsável pela guinada


## `Aeronave.Tail (S: float, b: float, mac: float, c12: list = None, th: float = None)`
**Descrição:**
Classe filha de `AeroSurface` para representar a _Empenagem Horizontal_ (EH)


## `Aeronave.Wing (S: _float_, b: float, mac: float, c12: list = None, th: float = None, Cm_CA: float = 0)`
**Descrição:**
Classe filha de `AeroSurface` para representar a _Asa_

**Inputs:**
  - Cm_CA: Coeficiente de momento da asa em torno do Centro Aerodinâmico


### `Aeronave.Wing.downwash ()`
**Descrição:**
Retorna o valor do DownWash teórico baseado no Nelson


## `Util.util.mach(V0: float | list, alt: float | list)`
**Descrição:**
Calcula o número de MACH de uma dada velocidade para uma certa altitude, utilizando interpolações e aproximações de dados atmosféricos.

**Inputs:**
  - V0: Velocidade desejada $[m/s]$
  - alt: Altitude desejada $[m]$





