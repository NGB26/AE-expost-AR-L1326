# AE expost AR-L1326

This project is to build the economic analysis of the project AR-L1326
(COVID-19 prototype). The AE will be an ex-post cost-benefit analysis.

## Modelo epidemiologico SIVQRD

Implementa en R el modelo SIVQRD descrito en `Modelo_SIR_para_PCR.docx`
(ecuaciones 23-31), calibrado contra la serie semanal de muertes de Our
World in Data (OWID) para Argentina. No depende de Excel: todos los
parametros se declaran a mano en `R/config.R`, con su unidad y su fuente
bibliografica (o la razon de no tenerla).

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
  simulate.R   Integracion numerica (deSolve) y ajuste secuencial de beta, sin estado global
  metrics.R    Funcion de perdida (MSRE-log) y metricas de diagnostico
  calibrate.R  Calibracion: baseline secuencial-alternante, alternativa conjunta,
               identificabilidad (perfil, multi-start), sensibilidad OAT, bootstrap
scripts/
  run.R                           Orquesta todo el pipeline y escribe outputs/
  modelo_sivqrd_argentina.R       Version anterior (script unico), queda de referencia
tests/
  testthat.R                      Punto de entrada de los tests
  testthat/helper-fixtures.R      Datos y config sinteticos (sin internet)
  testthat/test-*.R               92 expectativas sobre R/*.R
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

**Tiempo esperado:** con todos los diagnosticos activados (que es el
default), entre 5 y 10 minutos — el paso mas caro es la calibracion
conjunta global (~3-4 minutos). Para una corrida rapida (~15-30 segundos,
solo el baseline), editar el bloque `opciones` al principio de
`scripts/run.R` y poner en `FALSE` los pasos que no se necesiten:

```r
opciones <- list(
  forzar_descarga_owid   = FALSE,
  correr_calibracion_conjunta = FALSE,      # ~3-4 min
  correr_perfil_identificabilidad = FALSE,  # ~3-10s
  correr_multistart      = FALSE,           # ~30s
  correr_sensibilidad    = FALSE,           # ~1-2 min
  correr_bootstrap       = FALSE,           # ~20-60s
  n_rep_bootstrap        = 20
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
