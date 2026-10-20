# AE expost AR-L1326

This project is to build the economic analysis of the project AR-L1326
(COVID-19 prototype). The AE will be an ex-post cost-benefit analysis.

## Modelo epidemiologico SIVQRD

Implementa en R el modelo SIVQRD descrito en `Modelo_SIR_para_PCR.docx`
(ecuaciones 23-31), calibrado contra la serie semanal de muertes de Our
World in Data (OWID) para Argentina (2020-2023). No depende de Excel: todos
los parametros se declaran a mano en `R/config.R`, con su unidad y su fuente
bibliografica (o la razon de no tenerla).

### Especificacion de la calibracion

- **4 fechas de corte HARDCODEADAS** (en `R/config.R`) que definen **5
  segmentos**: comienzo de ASPO, comienzo de DISPO, comienzo de vacunacion y
  Omicron.
- **La UNICA variable de ajuste es `beta`** (un valor por segmento).
- **Funcion de perdida: RMSE de `log(1+muertes)`** del modelo vs. OWID
  (`R/metrics.R::calcular_rmse_log`).
- **Ajuste GLOBAL y SIMULTANEO**: los 5 `beta` se optimizan a la vez para
  minimizar el RMSE de toda la serie (`R/calibrate.R::calibrar_beta_global`),
  no segmento por segmento. El ajuste secuencial queda solo como comparacion.
- **Parametros FIJOS de literatura**: efectividad vacunal frente al contagio
  (`ve_contagio`), IFR, y de ahi `gamma`/`delta` (via la duracion del periodo
  infeccioso de Kissler et al.). **`theta` tambien es fijo** (no se calibra):
  el efecto del confinamiento lo absorbe el `beta` propio de cada segmento.
- **`sigma` se calcula de los datos** reales de vacunacion de OWID, no de
  literatura (`R/data.R::calcular_sigma_desde_owid`).

Ademas, `R/rates.R` arma una **estructura de tasas IMPLICITAS por semana**
(beta, sigma, delta, Re...) reconstruidas de los datos crudos al estilo de una
planilla de Markov, como diagnostico y punto de partida de la calibracion.

**Para entender las decisiones detras de este diseño** (que bugs se
encontraron y donde, por que se corrigio cada formula, los resultados de
comparar el metodo de calibracion secuencial contra uno conjunto, el
analisis de identificabilidad entre beta y theta, y que cosas quedaron sin
poder verificar), ver **[DECISIONS.md](DECISIONS.md)**. Este README es solo
la guia de instalacion y ejecucion.

### Estructura del repositorio

```
R/
  config.R     Todos los parametros del modelo (valor, unidad, fuente) + validacion
  data.R       Descarga/cache/validacion de OWID, agregacion semanal, calculo de sigma
  model.R      Ecuaciones diferenciales (SIVQRD), parametros biologicos derivados, R0/Re
  rates.R      Tasas IMPLICITAS por semana (reconstruccion tipo Markov-Excel) + resumen por tramo
  simulate.R   Integracion numerica (deSolve) y ajuste secuencial de beta, sin estado global
  metrics.R    Funcion de perdida principal (RMSE-log) + MSRE y metricas de diagnostico
  calibrate.R  Calibracion de beta (solo beta): global simultaneo (principal), secuencial
               (comparacion), sensibilidad a theta, sensibilidad OAT, multi-start, bootstrap
scripts/
  run.R                           Orquesta todo el pipeline y escribe outputs/
  modelo_sivqrd_argentina.R       Version anterior (script unico), queda de referencia historica
tests/
  testthat.R                      Punto de entrada de los tests
  testthat/helper-fixtures.R      Datos y config sinteticos (sin internet)
  testthat/test-*.R               Tests sobre R/*.R
data/
  raw/     Excel original (SIRD_Argentina_PCR.xlsx), solo como referencia historica
  cache/   Cache local de OWID (se genera sola al correr; no versionada en git)
outputs/   CSV, PNG, RDS y sessionInfo.txt de la ultima corrida de scripts/run.R
DECISIONS.md   Informe de decisiones: diagnostico, matematica, datos, calibracion
```

### Requisitos

- R >= 4.3 (probado con R 4.3.3).
- Paquetes: `tidyverse`, `deSolve`, `scales`, `jsonlite`, `digest`,
  `testthat` (>= 3.0). Versiones con las que se probo este refactor (ver
  tambien `outputs/sessionInfo.txt` tras correr `scripts/run.R`):

  | Paquete | Version |
  |---|---|
  | tidyverse | 2.0.0 |
  | deSolve | 1.40 |
  | scales | 1.3.0 |
  | jsonlite | 1.8.8 |
  | digest | 0.6.34 |
  | testthat | 3.2.1 |

- Acceso a internet la PRIMERA vez que se corre `scripts/run.R` (para
  descargar OWID); las corridas siguientes usan la cache local en
  `data/cache/` salvo que se pida `forzar_descarga_owid = TRUE`. Los tests
  (`tests/testthat/`) usan datos sinteticos y **no** necesitan internet.

#### Instalar dependencias

Sin `renv` (no se congelo un lockfile en este refactor: ver
[DECISIONS.md, seccion 5](DECISIONS.md#5-arquitectura-prioridad-4) para por
que). Instalar a mano:

```r
install.packages(c("tidyverse", "deSolve", "scales", "jsonlite", "digest", "testthat"))
```

o, en Debian/Ubuntu con paquetes del sistema (mas rapido, sin compilar):

```bash
sudo apt-get install -y r-cran-tidyverse r-cran-desolve r-cran-scales \
  r-cran-jsonlite r-cran-digest r-cran-testthat
```

### Correr los tests

Desde la raiz del repo:

```bash
Rscript tests/testthat.R
```

No requiere internet (usa `tests/testthat/helper-fixtures.R`, datos
sinteticos deterministicos). Deberia terminar sin fallos, mostrando algo
como:

```
config: ........
data: .............
integracion-sintetica: ..............
metrics: ..........
model: ...........................
simulate: ....................
```

### Correr el pipeline completo

Desde la raiz del repo:

```bash
Rscript scripts/run.R
```

Descarga (o usa la cache de) OWID, calibra el modelo, y escribe todo en
`outputs/`: tablas (`*.csv`), graficos (`*.png`), los objetos de resultado
completos (`*.rds`, para seguir explorando en una sesion de R con
`readRDS()`) y `sessionInfo.txt`.

**Tiempo esperado:** con la serie 2020-2023 (~199 semanas), cada calibracion
global tarda ~1-3 minutos (el ODE corre sobre toda la serie en cada evaluacion
del optimizador). Con los defaults (secuencial + global + sensibilidad a
theta), una corrida tarda ~7-10 minutos. Los diagnosticos que corren VARIAS
calibraciones estan **apagados** por defecto; prenderlos en el bloque
`opciones` al principio de `scripts/run.R` para el analisis completo:

```r
opciones <- list(
  forzar_descarga_owid      = FALSE,
  correr_secuencial         = TRUE,   # ajuste secuencial: comparacion + warm start (~10-20s)
  correr_multistart         = FALSE,  # robustez del optimo global (~5 calibraciones)
  correr_sensibilidad_theta = TRUE,   # barrido del theta FIJO (~4 calibraciones)
  correr_sensibilidad_oat   = FALSE,  # OAT sobre parametros de literatura (~11 calibraciones, lento)
  correr_bootstrap          = FALSE,  # bootstrap de los beta (muy lento)
  n_rep_bootstrap           = 20
)
```

**Nota de este entorno:** al correr `scripts/run.R` con TODOS los
diagnosticos en un proceso en segundo plano (background) en paralelo con
otra carga de trabajo, el proceso se interrumpio dos veces sin mensaje de
error (ver [DECISIONS.md, seccion 6](DECISIONS.md#6-que-no-se-pudo-verificar)).
Corrido en primer plano, sin otra carga concurrente, termino bien las dos
veces que se probo. Si al correrlo `outputs/` queda incompleto sin ningun
error en pantalla, es probable que sea esto: re-correr en primer plano, con
los recursos de la maquina dedicados a este proceso.

### Editar parametros o fechas de corte

Todo en `R/config.R`, funcion `construir_config()`. Cada parametro tiene su
`valor`, `unidad` y `fuente` (o un `# EDITAR`/`# CONFIRMAR` si es un valor
de trabajo sin cita bibliografica, o si hay una ambiguedad sin resolver —
ver [DECISIONS.md, seccion 2.4](DECISIONS.md#24-rezago-de-muertes-vs-tiempo-a-muerte-implicito-en-la-edo--no-resuelto-requiere-tu-confirmacion)
sobre `rezago_semanas` en particular). Despues de editar, `Rscript
scripts/run.R` valida la config de entrada (`validar_config()`, falla
explicito si algo quedo en `NA` o fuera de rango) antes de correr nada.

Dos parametros a revisar en particular: **`theta`** (fijo, no calibrado:
valor de trabajo `0.3`, marcado `# CONFIRMAR`; el barrido de sensibilidad en
`outputs/sensibilidad_theta.csv` muestra cuanto cambia el ajuste con el) e
**`ifr_base`** (IFR base, `# CONFIRMAR` si se prefiere uno especifico de
Argentina). La clinica (`dias_infecciosos`, `dias_infecciosos_vac`) sale de
Kissler et al.; `ve_contagio` de Rearte et al. (Lancet 2022, Argentina).
