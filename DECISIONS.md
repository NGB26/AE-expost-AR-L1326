# Informe de decisiones — Refactor del modelo SIVQRD (Argentina COVID-19)

Este documento acompaña al refactor de `scripts/modelo_sivqrd_argentina.R`
(script unico) hacia la arquitectura modular en `R/*.R` + `scripts/run.R` +
`tests/testthat/`. Sigue el procedimiento pedido: diagnostico con linea e
impacto, propuesta y justificacion, implementacion, resultados de ejecucion,
y que no se pudo verificar.

Todos los numeros de este informe salen de una corrida real en
`outputs/` (datos reales de OWID, cache con hash en `data/cache/`), no de
datos inventados. El archivo `outputs/sessionInfo.txt` fija las versiones
exactas de R y paquetes usadas.

## Indice

1. [Diagnostico (Prioridad 0)](#1-diagnostico-prioridad-0)
2. [Matematica (Prioridad 1)](#2-matematica-prioridad-1)
3. [Datos (Prioridad 2)](#3-datos-prioridad-2)
4. [Calibracion (Prioridad 3)](#4-calibracion-prioridad-3)
5. [Arquitectura (Prioridad 4)](#5-arquitectura-prioridad-4)
6. [Que NO se pudo verificar](#6-que-no-se-pudo-verificar)
7. [Lista de cambios](#7-lista-de-cambios)

---

## 1. Diagnostico (Prioridad 0)

### 1.1 Bloque duplicado y llaves huerfanas (lineas ~235-267 del script roto en `main`)

**Sintoma:** `Rscript -e 'parse("modelo_sivqrd_argentina.R")'` fallaba con
`unexpected '}'` en la linea 267. El archivo no se podia ni parsear.

**Causa:** un merge manual (commit directo a `main`, fuera de un PR) combino
una edicion hecha sobre una version VIEJA del script (que todavia leia
parametros de Excel) con el PR que ya habia reescrito esas mismas secciones.
El resultado dejo, entre la definicion de `ajustar_betas()` y su uso, un
fragmento sobrante de la version anterior de `calcular_error_tramo()`
(variable `comparacion` sin definir, `parametros_fijos` referenciando un
objeto que ya no existia como global, y una llave de cierre sin su apertura
correspondiente).

**Impacto:** el script no corria en absoluto. Explica literalmente el "no
veo que se pueda adaptar facilmente": no era una cuestion de legibilidad,
el archivo estaba roto a nivel de sintaxis.

**Resolucion:** eliminado. La logica de `ajustar_betas()` quedo en una unica
definicion limpia, ahora en `R/calibrate.R::calibrar_secuencial_alternante()`
(ver seccion 4).

### 1.2 `obtener_parametro(parametros_excel, ...)` y `ruta_excel_parametros`

**Sintoma:** aun si el error de sintaxis de 1.1 no existiera, las lineas
que llamaban a `obtener_parametro(parametros_excel, "...")` para `s0_frac`,
`i0_personas`, `dias_infecciosos`, `ifr_base`, `rezago_semanas`, `ve_muerte`,
`sigma` y `ve_contagio` habrian fallado con `object 'parametros_excel' not
found`: ni `library(readxl)` ni la construccion de `parametros_excel` existian
en el archivo (se habian eliminado en el PR anterior, que saco la lectura de
Excel por pedido explicito).

**Impacto:** el script fallaria en la primera linea de la seccion de
parametros, antes de llegar a ninguna simulacion.

**Resolucion:** `R/config.R` declara todos los parametros a mano (ver
seccion 5), sin ninguna dependencia de Excel. `ruta_excel_parametros`
(que ademas apuntaba a un path personal de Windows,
`C:/Users/NICOLASGA/OneDrive - ...`) se elimino: no tenia ningun uso real.

---

## 2. Matematica (Prioridad 1)

### 2.1 Conservacion `s+i+iv+v+q+r+d = 1`

Demostrada algebraicamente en `R/model.R` (comentario sobre
`derivadas_sivqrd()`): sumando las 7 derivadas, los terminos de nuevas
infecciones se cancelan entre `S`/`I` y entre `V`/`Iv`, el termino de
vacunacion se cancela entre `S`/`V`, y los terminos de salida de `I`, `Q` e
`Iv` se cancelan exactamente contra sus entradas a `R` y `D`. La suma de
derivadas es identicamente 0 para cualquier valor de los parametros.
Verificado tambien por simulacion en `tests/testthat/test-model.R`
("la poblacion total se conserva"): para una corrida de 20 semanas,
`|suma - 1| < 1e-6` en todas las semanas.

### 2.2 Positividad

No se agrego ningun piso artificial (`pmax(x, 0)`) a los compartimentos: en
un sistema donde cada termino de SALIDA de un compartimento es proporcional
a ese mismo compartimento (como este), la frontera `x=0` es invariante —si
un compartimento llega a 0, su derivada en ese instante es `>= 0`— asi que
no hace falta forzarlo. Se verifico (no se asumio) con una corrida de 20
semanas con beta y theta altos (`test-model.R`, "los compartimentos fisicos
no son negativos"): todos los valores quedan `>= -1e-8` (la tolerancia
cubre el error numerico de punto flotante, no una violacion real).

### 2.3 Parametrizacion de gamma/delta (y la correccion que se revirtio)

El pedido fue `gamma=(1-IFR)*7/D` y `delta=IFR*7/D`. Esa es exactamente la
formula que ya estaba en la version previa al ultimo commit directo a
`main`; en ese commit alguien la habia cambiado a
`delta = (IFR/(1-IFR)) * gamma_total` (y la version vacunada analoga).

**Por que se revierte:** con `gamma_total = 7/dias_infecciosos` fijo, la
UNICA combinacion que cumple a la vez (a) que `gamma+delta` reproduzca la
duracion total del periodo infeccioso `D` y (b) que la probabilidad de
muerte por riesgos competitivos, `delta/(gamma+delta)`, sea exactamente el
IFR, es `delta = IFR*gamma_total`, `gamma = (1-IFR)*gamma_total` (demostracion
algebraica completa en el comentario de
`R/model.R::derivar_parametros_biologicos()`). La formula alternativa
(`IFR/(1-IFR)`) rompe la condicion (a): `gamma+delta` deja de ser
`7/dias_infecciosos`, asi que "D dias de duracion del periodo infeccioso"
deja de ser cierto.

**Verificacion por simulacion** (no solo algebra), en
`test-model.R::"el IFR se recupera en un proceso competitivo simple"`: se
simula una poblacion que arranca enteramente en `I` (sin `S`, sin nuevas
infecciones), se la deja evolucionar 200 semanas (tiempo de sobra para que
`I` se vacie), y se mide `D/(D+R)` al final. Para IFR de prueba de 1%, 5% y
20%, el valor recuperado coincide con el IFR usado hasta `1e-4` de
tolerancia. Misma logica para vacunados (`gamma_v`, `delta_v`, con
`ifr_vacunados = ifr_base*(1-ve_muerte)`).

Esta misma derivacion se aplica **igual** a vacunados y no vacunados: no hay
un tratamiento especial, solo que `ifr_vacunados` reemplaza a `ifr_base` y
`dias_infecciosos_vac` a `dias_infecciosos`.

**INTERPRETACION DEL IFR que esto asume** (para que quede explicita, como
pide el enunciado): se asume que, mientras alguien esta en `I` (o `Iv`),
compite entre dos salidas de **hazard constante**: recuperarse o morir. Bajo
hazards constantes, la fraccion que termina muriendo es
`delta/(gamma+delta)` (formula estandar de riesgos competitivos
exponenciales) — y se define `delta`, `gamma` para que esa fraccion sea
exactamente el IFR. Es una simplificacion: en la realidad el riesgo de
muerte no es un hazard constante a lo largo de la enfermedad (es mas alto
mas tarde en el curso clinico), pero es la simplificacion estandar de un
modelo compartimental de un solo compartimento `I` sin sub-etapas.

### 2.4 Rezago de muertes vs. tiempo-a-muerte implicito en la EDO — **NO RESUELTO, requiere tu confirmacion**

Esto es lo mas importante que el pedido de auditoria encontro y que **no**
se puede resolver sin mas informacion.

El modelo tiene DOS fuentes de demora entre que alguien se infecta y que su
muerte aparece en la serie que se usa para calibrar:

1. **Implicita en la EDO:** el tiempo que alguien pasa en `I` antes de
   transicionar a `D` tiene media `1/gamma_total = dias_infecciosos/7`
   semanas (10 dias en el valor actual de `R/config.R`).
2. **Explicita, `rezago_semanas`:** un `dplyr::lag()` adicional de 3 semanas
   que se aplica a las muertes que emite el modelo ANTES de compararlas con
   las muertes reportadas por OWID (`R/metrics.R::alinear_muertes_con_rezago()`).

El Excel original (`SIRD_Argentina_PCR.xlsx`, hoja `Parametros`, fila
`Rezago infeccion -> muerte`) documenta esos 3 semanas como "19 dias segun
meta-analisis, redondeado a 3 semanas por el paso del modelo", **sin citar
el estudio**. Esto importa porque, segun que mida ese "19 dias":

- Si mide **tiempo biologico infeccion-a-muerte** (p.ej. Verity et al. 2020,
  Lancet Infect Dis, que reporta ~18 dias de onset-a-muerte): ese plazo YA
  incluye, en gran parte, los 10 dias que el modelo ya cuenta como
  `dias_infecciosos`. Aplicar el rezago ENCIMA de eso duplicaria parte del
  tiempo: el modelo terminaria representando ~10+21=31 dias de infeccion a
  muerte reportada, cuando la literatura dice ~19-21.
- Si mide **demora administrativa de reporte** (el tiempo entre que la
  muerte ocurre y que se carga en el sistema de vigilancia, tipicamente
  mas corto y un concepto DISTINTO del tiempo clinico), entonces el rezago
  es correctamente ADITIVO al tiempo biologico que ya modela la EDO, y no
  hay doble conteo.

**No tengo forma de saber cual de las dos cosas mide el "19 dias" del Excel
original sin el estudio citado.** El valor y el comportamiento (aditivo) se
mantuvieron IGUAL que en el script anterior, por continuidad, pero quedan
marcados en `R/config.R` (campo `rezago_semanas`) como pendientes de
confirmar, con esta misma explicacion. Si se puede identificar la fuente
exacta, el arreglo es inmediato (cambiar `rezago_semanas`, o si el 19 dias
ya incluye el periodo infeccioso, usar algo como `max(0, round((19 -
dias_infecciosos)/7))` en vez del valor fijo actual).

### 2.5 R0 y Re via matriz de proxima generacion (NGM)

Derivacion completa en el comentario de `R/model.R::calcular_r0()`. Los
compartimentos "infectados" (los unicos que generan nuevas infecciones) son
`i` e `iv`. Las matrices de la NGM, evaluadas en un punto `(s*, v*)`:

```
F = [ beta*s*     beta*s*   ]      V = [ theta+gamma+delta      0           ]
    [ beta_v*v*   beta_v*v* ]          [ 0                      gamma_v+delta_v ]
```

(las dos columnas de `F` son iguales porque tanto los infectados de origen
`S` como los de origen `V` contagian a traves del mismo pool `i+iv`). Esto
hace que `F` sea **rango 1** (producto externo de `(beta*s*, beta_v*v*)` y
`(1,1)`), asi que `F*V^-1` tiene un unico autovalor no nulo, igual a su
traza:

```
R = beta*s*/(theta+gamma+delta) + beta_v*v*/(gamma_v+delta_v)
```

que es EXACTAMENTE la ecuacion (31) del documento original. `R0` (ecuacion
30) es el caso particular `s*=1, v*=0` (nadie vacunado todavia). Verificado
numericamente (no solo algebraicamente) en `test-model.R`: para 9
combinaciones de `(s,v)`, el autovalor dominante de la NGM calculado via
`eigen()` coincide con la formula cerrada hasta `1e-10`.
`calcular_re_ngm_numerico()` queda en `R/model.R` como la verificacion
independiente, documentada como redundante con `calcular_re()` (mas lenta,
no se usa en el pipeline principal).

### 2.6 `tasa_contagio_semanal`: de `beta*s_prom` a incidencia integrada

**Problema con `beta*s_prom`:** no es una tasa de nuevas infecciones. Le
faltaba multiplicar por `(i+iv)` (la fuerza de infeccion depende de cuanta
gente esta efectivamente infectada, no solo de `beta` y `s`), y ademas
ignoraba por completo las infecciones de origen vacunado
(`beta_v*v*(i+iv)`, ecuacion 25/26).

**Solucion implementada:** se agrego un OCTAVO estado a la EDO,
`inf_acum` (`R/model.R::derivadas_sivqrd()`), puramente contable (no es
parte de la poblacion fisica, no entra en la suma de conservacion de la
seccion 2.1), con derivada `d(inf_acum)/dt = beta*s*(i+iv) +
beta_v*v*(i+iv)` — la suma de AMBOS flujos de nuevas infecciones. La
incidencia semanal exacta se obtiene igual que las muertes, con `diff()`
sobre ese acumulado (`R/simulate.R::simular_serie()`), y
`R/metrics.R::calcular_tasa_contagio_por_tramo()` promedia esa incidencia
(ya en fraccion de poblacion) por tramo. Es una integracion EXACTA del flujo
de nuevas infecciones, no una aproximacion con el valor promedio de `s`.

---

## 3. Datos (Prioridad 2)

Implementado en `R/data.R`. Dos problemas reales aparecieron al testear esto
contra los datos verdaderos de OWID (no son hipoteticos):

### 3.1 La ultima semana de la serie tenia 1 solo dia de datos

`fecha_fin_serie = 2022-06-05` es un domingo; la "semana 118" original
cubria solo ese unico dia (84 muertes reportadas ESE dia) pero se trataba
como si fueran las muertes de una semana entera — subestimando por
construccion la ultima semana del tramo Omicron. `agregar_owid_semanal()`
ahora cuenta `dias_observados` por semana y, por defecto
(`incluir_semanas_incompletas = FALSE`), DESCARTA las semanas con menos de 7
dias de la serie usable (con advertencia explicita). Se verifico que, en la
serie real 2020-2022, la UNICA semana afectada es la ultima (la primera
tiene sus 7 dias completos). Efecto practico: la serie usable paso de 118 a
117 semanas, y el tramo Omicron de `91-117` en vez de `91-118`.

### 3.2 Valores negativos y NA en `new_cases`/`new_deaths`

`validar_owid_diario()` los detecta y reporta como advertencia (no frena la
corrida: una correccion retroactiva de OWID es legitima, no es un error de
datos). **Verificado sobre los datos reales**: en el rango
`2020-03-08` a `2022-06-05` para Argentina, **no hay ningun valor negativo
ni NA** en `new_cases`/`new_deaths` (820 dias, 0 NA, 0 negativos). Si
aparecieran, `agregar_owid_semanal()` suma los negativos tal cual en el
total semanal (es lo correcto: la correccion debe descontarse de alguna
semana) y solo interviene si el TOTAL semanal da negativo (lo lleva a 0, con
advertencia), porque un conteo semanal negativo no tiene significado
epidemiologico utilizable por el modelo.

### 3.3 Cache: el hash del objeto en memoria no es estable

`digest::digest()` sobre un `data.frame` en memoria hashea su representacion
interna (tipos exactos, atributos), que NO sobrevive un roundtrip
`write_csv()`/`read_csv()` aunque los VALORES sean identicos. La primera
version de la cache fallaba su propio chequeo de integridad en la SEGUNDA
corrida (falso positivo de "cache corrupta"). Se cambio a hashear el
ARCHIVO en disco (`hash_archivo()`, sha256 sobre los bytes), que es estable
por construccion. Test de regresion en `test-data.R`.

### 3.4 sigma (modelo) vs. "vacunados completos" de OWID (dato): significados distintos

`sigma` es un **hazard semanal CONSTANTE**, calibrado para que, aplicado
desde `fecha_vacunacion` hasta el final de la serie, reproduzca la
cobertura FINAL observada
(`-ln(1-cobertura_final)/semanas_de_vacunacion`, en
`calcular_sigma_desde_owid()`). La campana real tuvo un ritmo variable
(arranca lento, acelera, se aplana), que un hazard constante por
construccion NO puede reproducir semana a semana.
`verificar_significado_vacunacion()` deja esto explicito comparando, semana
a semana, la cobertura observada (OWID) contra la modelada (`v(t)/N` del
modelo), en vez de asumir tacitamente que son la misma cosa — ver
`outputs/grafico_cobertura_vacunacion.png` y
`outputs/verificacion_cobertura_vacunacion.csv`: coinciden en el nivel
final (por construccion) pero divergen bastante en la forma del camino.

### 3.5 Snapshot reproducible

`R/data.R::obtener_datos_owid_diario()` cachea en `data/cache/` (CSV +
`.meta.json` con URL fuente, timestamp UTC, hash sha256 y numero de filas) y
usa esa cache en corridas posteriores salvo `forzar_descarga=TRUE`. El
directorio `data/cache/` queda fuera de git (ver `.gitignore`): es
regenerable y pesaria mas de lo razonable para versionar. El `.meta.json`
SI deja trazabilidad completa de que version de los datos se uso.

---

## 4. Calibracion (Prioridad 3)

Numeros de esta seccion: corrida real, `outputs/metricas_globales.csv`,
`outputs/comparacion_metodos_beta.csv`, `outputs/multistart_theta.csv`,
`outputs/sensibilidad_oat.csv`, `outputs/bootstrap_theta.csv`,
`outputs/perfil_identificabilidad_dispo.{csv,png}`.

### 4.1 Baseline (secuencial-alternante) vs. ajuste conjunto global

| Metodo | MSRE-log | MSE nivel | MSE log | beta (5 tramos) | theta |
|---|---|---|---|---|---|
| Secuencial-alternante (baseline) | **0.1612** | 9.46M | 4.94 | 2.00 / 1.63 / 0.94 / 2.00 / 2.00 | 0.088 |
| Conjunto global (`optim`, L-BFGS-B) | **0.0858** | 0.90M | 1.65 | 2.00 / 2.00 / 1.55 / 1.90 / 2.00 | 0.771 |

El ajuste conjunto, arrancando DESDE el optimo del baseline, converge
(`convergence=0` de `optim`) a un error **47% menor** en las 3 metricas, con
un `theta` casi 9 veces mayor y betas distintos en 3 de los 5 tramos (ver
`comparacion_metodos_beta.csv`). Esto es evidencia directa de que el
baseline secuencial, aunque es el metodo confirmado como punto de partida
del proyecto, **no encuentra el optimo global**: encuentra un optimo de un
procedimiento GREEDY (cada tramo resuelve su propio problema local, sin
anticipar como afecta eso al ajuste de los tramos siguientes via el
rezago). Costo: el baseline corre en ~9s; el conjunto, en esta maquina,
tarda **~3.5 minutos** (80 evaluaciones de una funcion que corre 5
integraciones de EDO cada una).

### 4.2 Identificabilidad entre beta y theta

**Perfil de perdida** (`outputs/perfil_identificabilidad_dispo.png`, tramo
DISPO, grilla de 12x12 con los demas beta fijos en el optimo del baseline):
el mapa de calor muestra una **cresta diagonal** clara — la banda de error
bajo se desplaza hacia beta mas alto a medida que theta sube. Es la firma
visual esperada de no-identificabilidad: `theta` saca gente de `I` hacia
`Q` (dejan de contagiar) de forma parecida a como bajar `beta` reduce
directamente los contagios, asi que distintas combinaciones `(beta,theta)`
explican curvas de muertes casi iguales. (1 de las 144 celdas, en la
esquina de beta y theta simultaneamente maximos, no convergio
numericamente por rigidez del sistema y quedo en NA explicito, no en un
numero enganoso — ver `R/simulate.R`.)

**Multi-start** (4 puntos de partida para theta): el resultado final
depende fuertemente de donde arranca:

| theta inicial | theta final | error |
|---|---|---|
| 0 | 0.00004 | 0.244 |
| 0.2 | 0.323 | 0.135 |
| 0.5 | 0.774 | **0.095** |
| 1 | 1.135 | 0.131 |

Cuatro arranques, cuatro `theta` finales sustancialmente distintos y
ninguno coincide con el `theta=0.088` del baseline con `theta_inicial=0.1`
(que es el que usa el pipeline por defecto). El error tambien varia 2.5x
entre arranques. Esto confirma, por fuera del perfil de perdida, que el
problema de calibracion tiene multiples optimos locales razonables.

**Sensibilidad OAT:** `theta` es MUY sensible a los parametros fijos,
particularmente:
- `dias_infecciosos`: 10 -> theta=0.088; 7 -> theta=0.384; 14 -> theta=0.020.
- `ifr_base`: 0.006 -> theta=0.088; 0.004 -> theta=0.466; 0.009 -> theta=0.371.
- `ve_contagio`: 0.586 -> theta=0.088; 0.68 -> theta=0.772.

Es otra forma del mismo problema de identificabilidad: como `theta` no esta
anclado a ningun dato independiente de las muertes (no hay una serie
observada de "gente en cuarentena" contra la cual calibrarlo), absorbe
cualquier error de especificacion de los demas parametros.

**Bootstrap de theta** (20 repeticiones, residuos de log1p resampleados,
betas del baseline fijos — ver limitaciones abajo): IC empirico 95% =
**[0.078, 0.861]**. El intervalo es ANCHO (cubre casi todo el rango
`[theta_min, theta_max]=[0,2]` relevante) y el punto estimado del baseline
(0.088) queda cerca del extremo inferior, no del centro.

**Conclusion sobre identificabilidad:** con la informacion disponible (solo
una serie semanal de muertes, sin datos independientes de ocupacion de
cuarentena/aislamiento), `theta` no esta bien identificado. Los resultados
de `theta` de este modelo (en cualquiera de los dos metodos) no deberian
interpretarse como una estimacion precisa de "la" tasa de cuarentena real;
son, como maximo, un valor que ayuda a que el modelo ajuste mejor la curva
de muertes, dentro de una familia amplia de valores casi-equivalentes. Esto
es una limitacion del DISEÑO del problema (que dato se usa para calibrar),
no un bug de la implementacion.

**Advertencia sobre el bootstrap:** no es un intervalo de confianza formal.
Asume residuos intercambiables entre semanas (no lo son del todo: hay
autocorrelacion serial en una curva epidemica) y fija los beta del baseline
(no propaga su incertidumbre). Se documenta en el codigo
(`R/calibrate.R::bootstrap_theta()`) como una cota INFERIOR de la
incertidumbre real, no como el intervalo definitivo.

### 4.3 Sin inferencia causal

Ningun resultado de este informe (beta por tramo, theta, R0/Re) se
interpreta como el efecto causal de ninguna politica o del prestamo: el
modelo no tiene una estrategia de identificacion causal (no hay grupo de
control, variacion exogena, ni diseño cuasi-experimental). `beta(t)` captura
TODO lo que cambio en cada periodo (politica, estacionalidad, variantes,
comportamiento), no un efecto aislado. Esto ya estaba documentado en el
`LEEME` del Excel original ("NO estima el efecto causal del prestamo BID") y
se mantiene como principio para este refactor.

---

## 5. Arquitectura (Prioridad 4)

```
R/
  config.R     Parametros (valor, unidad, fuente) + validar_config()
  data.R       OWID: descarga, cache+hash, validacion, agregacion semanal, sigma
  model.R      EDO (derivadas_sivqrd), parametros biologicos derivados, R0/Re (NGM)
  simulate.R   simular_tramo, simular_serie, ajustar_betas_secuencial (puras)
  metrics.R    MSRE-log, MSE nivel/log, residuos, tasas por tramo
  calibrate.R  secuencial-alternante, conjunto global, perfil, multistart,
               sensibilidad OAT, bootstrap
scripts/
  run.R        Orquestacion end-to-end -> outputs/
tests/
  testthat.R               punto de entrada
  testthat/helper-fixtures.R   datos y config sinteticos (sin internet)
  testthat/test-*.R        92 expectativas, ver seccion de ejecucion
data/
  cache/       cache de OWID (CSV + meta.json), gitignored
  raw/         Excel original (solo referencia historica, ya no se lee)
outputs/       CSV, PNG, RDS, sessionInfo.txt de la ultima corrida de run.R
```

Sin estado global: toda funcion en `R/*.R` recibe sus insumos como
argumentos explicitos (ver firmas de funciones). La duplicacion entre
`ajustar_betas()` y `ajustar_betas_con_beta_fijo()` del script original se
elimino: `R/simulate.R::simular_serie()` (forward, sin optimizar) y
`R/simulate.R::ajustar_betas_secuencial()` (con `optimize()` por tramo)
comparten la logica de recorrer tramos; `R/calibrate.R` combina ambas segun
haga falta fijar beta o theta.

`renv`: no se agrego (no hay acceso a un espejo de CRAN con versiones
congeladas en este entorno para generar un lockfile real); `outputs/sessionInfo.txt`
cumple el mismo objetivo de trazabilidad de versiones para cada corrida.
Ver README.md para instalar manualmente las mismas versiones si hace falta
una reproduccion exacta en otra maquina.

---

## 6. Que NO se pudo verificar

- **Las ecuaciones (23)-(31) en si:** se tomaron literalmente de
  `Modelo_SIR_para_PCR.docx` (extraido con pandoc en una sesion anterior de
  este mismo proyecto). No se verificaron contra una segunda fuente
  independiente (un libro de texto, un paper que derive el mismo modelo):
  la UNICA fuente es ese documento. Si existe una version mas reciente o
  corregida, no se detectaria.
- **La fuente exacta del rezago de 3 semanas** (seccion 2.4): el Excel
  original cita "19 dias segun meta-analisis" sin nombrar el estudio. No se
  pudo determinar si mide tiempo biologico o demora de reporte.
- **Los odds ratio de Rearte et al. 2022** (usados para `ve_contagio`, ya
  completados en un paso anterior de este proyecto): se obtuvieron de
  resultados de busqueda web, no del PDF original del paper (Lancet;
  `thelancet.com`, `pmc.ncbi.nlm.nih.gov` y `researchgate.net` estan
  bloqueados por la politica de red de este entorno — se confirmo con
  `curl` que la conexion es rechazada por el proxy). Los valores puntuales
  (OR=0.36/0.32/0.56) aparecieron consistentes en dos busquedas
  independientes, pero no hay forma de confirmar los intervalos de
  confianza ni la tabla completa sin acceso al texto completo.
- **El mix de vacunas aplicadas en Argentina** usado para ponderar
  `ve_contagio`: la cifra de Chequeado es de **septiembre de 2021**
  (aplicaciones), no de la campana completa hasta 2022; se uso por ser la
  mejor fuente de DOSIS APLICADAS (no solo recibidas) encontrada, pero no
  es exactamente el periodo completo del modelo.
- **Corrida completa de `scripts/run.R` con TODOS los diagnosticos activados
  en un unico proceso en background**: en este entorno, el proceso en
  background se interrumpio sin error dos veces (aparentemente por
  contencion de recursos al correr en paralelo con otros procesos de R de
  esta misma sesion), aunque el sistema de notificaciones lo reporto como
  "completado". La corrida que SI se uso para los numeros de este informe
  se hizo en foreground, sin otra carga concurrente, y termino
  correctamente (`outputs/` completo, "Listo." al final). Si se corre esto
  en CI o en otra maquina, correr `scripts/run.R` en foreground (o con
  suficientes recursos dedicados) hasta confirmarlo una vez.

---

## 7. Lista de cambios

Resumen (no un diff linea por linea; usar `git diff` sobre el commit de este
refactor para el detalle completo):

- **Nuevo:** `R/config.R`, `R/data.R`, `R/model.R`, `R/simulate.R`,
  `R/metrics.R`, `R/calibrate.R`, `scripts/run.R`, `tests/testthat.R`,
  `tests/testthat/*.R` (7 archivos), `DECISIONS.md` (este archivo), este
  README actualizado.
- **Eliminado del script original:** el bloque duplicado/roto, toda
  referencia a Excel (`library(readxl)`, `obtener_parametro()`,
  `parametros_excel`, `ruta_excel_parametros`), la logica de
  `ajustar_betas()`/`ajustar_betas_con_beta_fijo()` duplicada.
- **Corregido:** formula de `delta`/`delta_v` (vuelta a `IFR*gamma_total`),
  hash de cache (de objeto en memoria a archivo en disco), semana final
  incompleta tratada como semana completa.
- **Agregado (antes no existia):** octavo estado `inf_acum` (incidencia
  integrada), R0/Re derivados y verificados via NGM, calibracion conjunta
  global como alternativa comparable, perfil de identificabilidad,
  multi-start, sensibilidad OAT, bootstrap de theta, validacion de datos
  (columnas/tipos/NA/negativos/semanas incompletas), cache con hash y
  metadata, comparacion explicita sigma-vs-dato-observado, 92 tests
  automatizados.
- **`scripts/modelo_sivqrd_argentina.R` (el script unico anterior):** queda
  en el repo sin tocar, como referencia; el pipeline vigente es
  `scripts/run.R` + `R/*.R`.
