# ============================================================================
# MODELO EPIDEMIOLOGICO SIVQRD - ARGENTINA COVID-19
# Insumo para el analisis economico ex post del PCR AR-L1326 (BID)
#
# Implementa en R el modelo SIVQRD descrito en "Modelo_SIR_para_PCR.docx"
# (ecuaciones 23-29). Ya no depende del Excel "SIRD_Argentina_PCR.xlsx": los
# parametros biologicos y las fechas de corte de beta se declaran a mano en
# este mismo script (Bloque 2.1 y 2.2), y la unica fuente de datos externa es
# OWID (de ahi tambien se toma la poblacion N).
#
# Fuente de calibracion: muertes semanales de Our World in Data (OWID).
#
# Metodo de ajuste: proceso ITERATIVO de dos pasos (Bloque 3).
#   (A) Con theta fijo, se ajustan los 5 beta de forma secuencial por tramo:
#       cada beta minimiza el error dentro de SU tramo unicamente, usando el
#       estado final del tramo anterior como condicion inicial.
#   (B) Con esos 5 beta ya fijos, se ajusta theta (un unico valor para toda
#       la serie) para minimizar el error GLOBAL.
# (A) y (B) se repiten hasta que el error global deje de mejorar (o se
# alcance el maximo de iteraciones). Los demas parametros (gamma, delta,
# sigma, beta_v, gamma_v, delta_v) son fijos y no se calibran.
#
# FUNCION DE PERDIDA: MSRE (Mean Squared Relative Error) semanal, calculado
# sobre log1p(muertes) en vez de sobre el nivel (ver Bloque 3 para el porque),
# tanto para beta como para theta.
# ============================================================================

library(tidyverse)   # dplyr, tidyr, purrr, ggplot2
library(lubridate)   # manejo de fechas
library(deSolve)     # integracion numerica de las ecuaciones diferenciales
library(scales)      # formato de ejes en los graficos

# ----------------------------------------------------------------------------
# BLOQUE 0 - CONSTANTES DE SERIE Y FECHAS DE CORTE DE BETA (EDITAR ACA)
# ----------------------------------------------------------------------------

fecha_inicio_serie <- as.Date("2020-03-08")   # semana 1
fecha_fin_serie    <- as.Date("2022-06-05")   # ultima semana de la serie

# Fechas de corte de beta: definen el inicio de cada uno de los 5 tramos del
# modelo SIVQRD. Van aca (y no mas abajo, junto a tramos_def) porque sigma
# (Bloque 2.1) tambien las necesita.
fecha_aspo       <- as.Date("2020-03-16")   # EDITAR: inicio del ASPO
fecha_dispo      <- as.Date("2020-06-01")   # EDITAR: inicio del DISPO
fecha_vacunacion <- as.Date("2020-12-29")   # EDITAR: inicio de la vacunacion
fecha_omicron    <- as.Date("2021-12-01")   # EDITAR: aparicion de Omicron

# Numero de semana (entero) que CONTIENE una fecha dada, con semana 1 =
# fecha_inicio_serie. Se usa tanto para agregar los datos diarios a semanal
# como para ubicar las fechas de corte de beta dentro de esa misma grilla.
semana_de_fecha <- function(fecha) as.integer(floor(as.numeric(fecha - fecha_inicio_serie) / 7)) + 1L

# ----------------------------------------------------------------------------
# BLOQUE 1 - IMPORTACION DE DATOS (OWID)
# ----------------------------------------------------------------------------
# Archivo historico completo (congelado; Argentina dejo de reportar en 2022).
# https://github.com/owid/covid-19-data/blob/master/public/data/owid-covid-data.csv

url_owid <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"

owid_argentina_diario <- read_csv(url_owid, show_col_types = FALSE) %>%
  filter(iso_code == "ARG") %>%
  select(date, new_cases, new_deaths, new_tests, positive_rate, people_fully_vaccinated, population) %>%
  filter(date >= fecha_inicio_serie, date <= fecha_fin_serie) %>%
  mutate(semana = semana_de_fecha(date))

datos_owid_semanal <- owid_argentina_diario %>%
  group_by(semana) %>%
  summarise(
    fecha                    = min(date),
    casos_owid               = sum(new_cases, na.rm = TRUE),
    muertes_owid             = sum(new_deaths, na.rm = TRUE),
    tests_owid               = sum(new_tests, na.rm = TRUE),
    positividad_owid         = mean(positive_rate, na.rm = TRUE),
    vacunados_completos_owid = last(na.omit(people_fully_vaccinated)),
    poblacion_owid            = mean(population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(semana)

# ----------------------------------------------------------------------------
# BLOQUE 2 - FUNCIONES Y MODELO SIVQRD (ecuaciones diferenciales)
# ----------------------------------------------------------------------------

## 2.1 Parametros del modelo (EDITAR ACA) ------------------------------------
# N se extrae de los datos de OWID (columna "population", constante en toda
# la serie). El resto de los parametros biologicos y de atribucion se declara
# a mano: completar cada uno con su fuente bibliografica antes de usar el
# modelo para resultados finales.

n_poblacion <- mean(datos_owid_semanal$poblacion_owid)

s0_frac          <- 1       # EDITAR: fraccion susceptible inicial (S0/N)
i0_personas      <- 100     # EDITAR: infecciosos iniciales (I0)
dias_infecciosos <- 10      # EDITAR: duracion del periodo infeccioso no vacunado, en dias (D)
ifr_base         <- 0.006   # EDITAR: IFR base pre-vacunacion
rezago_semanas   <- 3       # EDITAR: rezago infeccion -> muerte, en semanas
ve_muerte        <- 0.7     # EDITAR: efectividad vacunal contra la muerte

dias_infecciosos_vac <- dias_infecciosos   # EDITAR si la literatura da una duracion distinta para Iv

# sigma: tasa semanal de vacunacion de susceptibles (S -> V). No es un
# parametro biologico de la literatura sino el RITMO de la propia campana de
# vacunacion argentina, asi que se extrae de los datos de OWID en vez de
# hardcodearse: es el hazard semanal CONSTANTE que, aplicado desde el inicio
# de la vacunacion (fecha_vacunacion) hasta el final de la serie, reproduce
# la cobertura de "vacunados completos" finalmente observada. Si
# V(T)/N = 1 - exp(-sigma*T), con T en semanas desde el inicio de la
# vacunacion, entonces sigma = -ln(1 - V(T)/N) / T.
semanas_de_vacunacion <- max(datos_owid_semanal$semana) - semana_de_fecha(fecha_vacunacion)
cobertura_final_owid  <- tail(na.omit(datos_owid_semanal$vacunados_completos_owid), 1) / n_poblacion
sigma <- -log(1 - cobertura_final_owid) / semanas_de_vacunacion

# ve_contagio: efectividad vacunal contra el CONTAGIO (no contra la muerte;
# define beta_v = beta*(1-ve_contagio)). Se construye igual que ve_muerte en
# el Excel original: promedio ponderado por el mix de vacunas efectivamente
# aplicadas en Argentina.
#
# Odds ratios de infeccion a las 2 dosis, de Rearte A, et al. "Effectiveness
# of rAd26-rAd5, ChAdOx1 nCoV-19, and BBIBP-CorV vaccines for risk of
# infection with SARS-CoV-2 and death due to COVID-19 in people older than
# 60 years in Argentina: a test-negative, case-control, and retrospective
# longitudinal study". Lancet. 2022;399(10331):1254-1264.
# https://doi.org/10.1016/S0140-6736(22)00011-3
#   rAd26-rAd5 (Sputnik V):        OR = 0.36 -> VE = 0.64
#   ChAdOx1 nCoV-19 (AstraZeneca):  OR = 0.32 -> VE = 0.68
#   BBIBP-CorV (Sinopharm):         OR = 0.56 -> VE = 0.44
#
# Mix de dosis APLICADAS en Argentina a sep-2021, segun Chequeado con datos
# del Monitor Publico de Vacunacion (AstraZeneca 33%, Sinopharm 32%,
# Sputnik V 29%, Moderna 5%):
# https://chequeado.com/hilando-fino/cuales-son-las-vacunas-contra-el-coronavirus-mas-aplicadas-en-la-argentina/
# Moderna no tiene OR de infeccion en Rearte et al.: se renormalizan los
# pesos sobre las 3 vacunas con dato (94% del total aplicado a esa fecha).
ve_infeccion_por_vacuna <- c(sputnik_v = 1 - 0.36, astrazeneca = 1 - 0.32, sinopharm = 1 - 0.56)
peso_por_vacuna         <- c(sputnik_v = 0.29,     astrazeneca = 0.33,     sinopharm = 0.32)
peso_por_vacuna         <- peso_por_vacuna / sum(peso_por_vacuna)
ve_contagio <- sum(ve_infeccion_por_vacuna * peso_por_vacuna)   # ~0.586

# theta (tasa de cuarentena) ya NO se fija a mano: se calibra en el Bloque 3
# junto con beta (ver proceso iterativo). Aca solo se declara el punto de
# partida y el rango de busqueda del optimizador.
theta_inicial <- 0.1   # EDITAR si se quiere otro punto de partida para la calibracion
theta_min     <- 0     # EDITAR: cota inferior de busqueda de theta
theta_max     <- 2     # EDITAR: cota superior de busqueda de theta

# gamma y delta se derivan de la tasa de salida total del compartimento I
# (gamma_total = 7/D, en semanas) partida segun el IFR, de forma que ambos
# queden como los dos parametros independientes que define el modelo SIVQRD
# (ecuaciones 23-29). Con gamma = (1-IFR)*gamma_total y delta = IFR*gamma_total
# se cumplen a la vez las dos condiciones que definen a D y al IFR:
# gamma+delta = gamma_total = 7/D (la duracion total del periodo infeccioso
# es D) y delta/(gamma+delta) = IFR (la probabilidad de morir, dado que te
# infectaste, es el IFR). Cualquier otra combinacion que no sea esta rompe
# una de las dos igualdades.
gamma_total <- 7 / dias_infecciosos
delta       <- ifr_base * gamma_total        # tasa de mortalidad de I  (parametro delta, ec. 24 y 29)
gamma       <- (1 - ifr_base) * gamma_total  # tasa de recuperacion de I (parametro gamma, ec. 24 y 28)

ifr_vacunados   <- ifr_base * (1 - ve_muerte)
gamma_total_vac <- 7 / dias_infecciosos_vac
delta_v         <- ifr_vacunados * gamma_total_vac         # tasa de mortalidad de Iv (ec. 25 y 29)
gamma_v         <- (1 - ifr_vacunados) * gamma_total_vac   # tasa de recuperacion de Iv (ec. 25 y 28)

## 2.2 Tramos de beta ---------------------------------------------------------
# Las fechas de corte (fecha_aspo, fecha_dispo, fecha_vacunacion,
# fecha_omicron) se declaran en el Bloque 0, no aca: sigma (arriba) tambien
# las necesita.

# cortes_semana[k] = primera semana (entera) del tramo k, en la misma grilla
# semanal que los datos (ver semana_de_fecha en el Bloque 0).
cortes_semana <- c(1,
                    semana_de_fecha(fecha_aspo),
                    semana_de_fecha(fecha_dispo),
                    semana_de_fecha(fecha_vacunacion),
                    semana_de_fecha(fecha_omicron),
                    max(datos_owid_semanal$semana) + 1)

tramos_def <- tibble(
  tramo      = 1:5,
  nombre     = c("Pre-intervencion", "ASPO", "DISPO", "Vacunacion", "Omicron"),
  semana_ini = cortes_semana[1:5],
  semana_fin = cortes_semana[2:6] - 1
)

## 2.3 Ecuaciones diferenciales del SIVQRD (ecuaciones 23 a 29) -------------
# Estado: fracciones de poblacion s, i, iv, v, q, r, d (s+i+iv+v+q+r+d = 1).
# theta y sigma se "activan" recien desde el hito de politica correspondiente
# (no hubo cuarentena antes del ASPO ni vacunas antes de diciembre de 2020):
# esto es un supuesto de implementacion, no parte de las ecuaciones en si, y
# queda documentado aqui para que pueda discutirse.

derivadas_sivqrd <- function(t, estado, parametros) {
  with(as.list(c(estado, parametros)), {

    beta_t   <- beta_fun(t)
    beta_v_t <- beta_t * (1 - ve_contagio)
    theta_t  <- if (t >= semana_inicio_aspo) theta else 0
    sigma_t  <- if (t >= semana_inicio_vacunacion) sigma else 0

    d_s  <- -beta_t * s * (i + iv) - sigma_t * s                     # Ecuacion (23)
    d_i  <-  beta_t * s * (i + iv) - i * (theta_t + delta + gamma)   # Ecuacion (24)
    d_iv <-  beta_v_t * v * (i + iv) - iv * (delta_v + gamma_v)      # Ecuacion (25)
    d_v  <-  sigma_t * s - beta_v_t * v * (i + iv)                   # Ecuacion (26)
    d_q  <-  theta_t * i - (gamma + delta) * q                       # Ecuacion (27)
    d_r  <-  gamma * (i + q) + gamma_v * iv                          # Ecuacion (28)
    d_d  <-  delta * (i + q) + delta_v * iv                          # Ecuacion (29)

    list(c(s = d_s, i = d_i, iv = d_iv, v = d_v, q = d_q, r = d_r, d = d_d))
  })
}

## 2.4 Simulacion de un tramo con beta constante -----------------------------

# Arma la lista de parametros fijos para un theta dado (el resto no cambia
# durante la calibracion). Se usa tanto para ajustar beta (theta fijo) como
# para ajustar theta (beta fijo) en el Bloque 3.
construir_parametros_fijos <- function(theta_valor) {
  list(
    theta = theta_valor, sigma = sigma, ve_contagio = ve_contagio,
    gamma = gamma, delta = delta, gamma_v = gamma_v, delta_v = delta_v,
    semana_inicio_aspo = semana_de_fecha(fecha_aspo),
    semana_inicio_vacunacion = semana_de_fecha(fecha_vacunacion)
  )
}

simular_tramo <- function(estado_inicial, semana_ini, semana_fin, beta_valor, parametros_fijos) {
  # Se pide un punto extra al final (semana_fin + 1): representa el estado al
  # INICIO de la semana siguiente (= fin de la ultima semana del tramo), y es
  # lo que permite computar el flujo de muertes de esa ultima semana con
  # diff(). Ese punto extra se descarta luego para graficos/promedios (ver
  # Bloque 3 y 4), para no duplicarlo con la primera semana del tramo que le
  # sigue.
  tiempos <- seq(semana_ini, semana_fin + 1, by = 1)
  parms <- c(parametros_fijos, list(beta_fun = function(t) beta_valor))
  ode(y = estado_inicial, times = tiempos, func = derivadas_sivqrd, parms = parms, method = "lsoda") %>%
    as.data.frame() %>%
    as_tibble() %>%
    rename(semana = time)
}

# R0 (ec. 30) y Re(t) (ec. 31), para reportar en el Bloque 4.
calcular_r0 <- function(beta_valor, theta_valor) beta_valor / (theta_valor + gamma + delta)

calcular_re <- function(beta_valor, theta_valor, s_prom, v_prom) {
  (beta_valor * s_prom) / (theta_valor + gamma + delta) +
    (beta_valor * (1 - ve_contagio) * v_prom) / (gamma_v + delta_v)
}

# ----------------------------------------------------------------------------
# BLOQUE 3 - AJUSTE Y CALIBRACION DEL MODELO (proceso iterativo beta <-> theta)
# ----------------------------------------------------------------------------

estado_inicial <- c(
  s = s0_frac - i0_personas / n_poblacion,
  i = i0_personas / n_poblacion,
  iv = 0, v = 0, q = 0, r = 0, d = 0
)

# Error de ajuste (MSRE, Mean Squared Relative Error, agregacion semanal) de
# una serie de "emisiones" (muertes nuevas del modelo, antes del rezago) ya
# simulada, contra las muertes semanales de OWID. Se reutiliza tanto para
# calibrar beta y theta como para reportar el error global en el Bloque 4.
#
# El MSRE se calcula sobre log1p(muertes), no sobre el nivel: el error
# relativo en NIVEL esta acotado por abajo en -1 (subestimar) pero no tiene
# techo por arriba (sobreestimar), asi que en una dinamica exponencial como
# esta el optimizador tiende a elegir beta casi nulo con tal de no arriesgar
# una sobreestimacion explosiva (se probo: el tramo mas largo terminaba en el
# piso del rango de busqueda y el modelo caia a ~0 el resto de la serie). El
# error relativo sobre log1p(.) no tiene esa asimetria y mantiene el espiritu
# del MSRE (pondera por igual semanas con pocas y con muchas muertes).
#
# MSRE = promedio( ((log1p(muertes_modelo) - log1p(muertes_owid)) / log1p(muertes_owid))^2 ).
# Las semanas con 0 muertes observadas se excluyen (log1p(0) = 0, el error
# relativo no esta definido ahi); en la serie de Argentina 2020-2022 no hay
# ninguna, pero se cubre el caso igual.
#
# El rezago (dplyr::lag) se calcula SIEMPRE sobre la serie COMPLETA de
# emisiones recibida (arrange() la ordena primero), para que las primeras
# semanas de un tramo tomen el rezago del historial de emisiones de los
# tramos anteriores. semanas_objetivo filtra DESPUES de aplicar el rezago, no
# antes, para no perder ese historial al pedir el error de un tramo puntual.
calcular_error_global <- function(emisiones, semanas_objetivo = NULL) {
  comparacion <- emisiones %>%
    arrange(semana) %>%
    mutate(muertes_comparables = dplyr::lag(muertes_emitidas, n = rezago_semanas, default = 0)) %>%
    left_join(datos_owid_semanal %>% select(semana, muertes_owid), by = "semana")
  if (!is.null(semanas_objetivo)) {
    comparacion <- comparacion %>% filter(semana %in% semanas_objetivo)
  }
  log_modelo <- log1p(comparacion$muertes_comparables)
  log_owid   <- log1p(comparacion$muertes_owid)
  error_relativo <- (log_modelo - log_owid) / log_owid
  error_relativo[comparacion$muertes_owid == 0] <- NA
  mean(error_relativo^2, na.rm = TRUE)
}

# (A) Ajusta los 5 beta de forma secuencial por tramo, para un theta fijo:
# cada beta minimiza el error DENTRO de su propio tramo, usando el estado
# final del tramo anterior como condicion inicial del siguiente.
ajustar_betas <- function(theta_valor) {
  parametros_fijos <- construir_parametros_fijos(theta_valor)
  estado_actual     <- estado_inicial
  emisiones         <- tibble(semana = numeric(), muertes_emitidas = numeric())
  resultados_sim    <- list()
  betas_ajustados   <- numeric(nrow(tramos_def))

  for (k in seq_len(nrow(tramos_def))) {

    semana_ini_k  <- tramos_def$semana_ini[k]
    semana_fin_k  <- tramos_def$semana_fin[k]
    semanas_tramo <- seq(semana_ini_k, semana_fin_k, by = 1)

    calcular_error_tramo <- function(beta_k) {
      sim <- simular_tramo(estado_actual, semana_ini_k, semana_fin_k, beta_k, parametros_fijos)
      muertes_nuevas <- diff(sim$d) * n_poblacion
      emisiones_candidatas <- tibble(semana = semanas_tramo, muertes_emitidas = muertes_nuevas)
      calcular_error_global(bind_rows(emisiones, emisiones_candidatas), semanas_objetivo = semanas_tramo)
    }

    ajuste <- optimize(calcular_error_tramo, interval = c(0.01, 2), tol = 1e-5)
    betas_ajustados[k] <- ajuste$minimum

    sim_final <- simular_tramo(estado_actual, semana_ini_k, semana_fin_k, betas_ajustados[k], parametros_fijos)
    muertes_nuevas_final <- diff(sim_final$d) * n_poblacion

    emisiones <- bind_rows(emisiones, tibble(semana = semanas_tramo, muertes_emitidas = muertes_nuevas_final)) %>%
      arrange(semana)

    # Se descarta el punto extra (semana_fin_k + 1): es solo el estado
    # inicial del tramo siguiente, no una semana propia de este tramo.
    resultados_sim[[k]] <- sim_final %>%
      filter(semana <= semana_fin_k) %>%
      mutate(tramo = k, nombre_tramo = tramos_def$nombre[k])

    estado_actual <- sim_final %>%
      slice_tail(n = 1) %>%
      select(s, i, iv, v, q, r, d) %>%
      unlist()
  }

  list(betas = betas_ajustados, emisiones = emisiones, simulacion = bind_rows(resultados_sim))
}

# (B) Para un vector de beta YA fijo, corre la serie completa con un theta
# dado y devuelve su error global. Se usa para calibrar theta en el paso (B).
evaluar_theta <- function(theta_candidato, betas_fijos) {
  calcular_error_global(ajustar_betas_con_beta_fijo(theta_candidato, betas_fijos)$emisiones)
}

# Corre la serie completa (los 5 tramos, uno atras del otro) con un vector de
# beta ya fijo (sin optimizar) y un theta dado. Es la contraparte de
# ajustar_betas() cuando lo que se quiere variar es theta, no beta.
ajustar_betas_con_beta_fijo <- function(theta_valor, betas_fijos) {
  parametros_fijos <- construir_parametros_fijos(theta_valor)
  estado_actual     <- estado_inicial
  emisiones         <- tibble(semana = numeric(), muertes_emitidas = numeric())
  resultados_sim    <- list()

  for (k in seq_len(nrow(tramos_def))) {
    semana_ini_k  <- tramos_def$semana_ini[k]
    semana_fin_k  <- tramos_def$semana_fin[k]
    semanas_tramo <- seq(semana_ini_k, semana_fin_k, by = 1)

    sim <- simular_tramo(estado_actual, semana_ini_k, semana_fin_k, betas_fijos[k], parametros_fijos)
    muertes_nuevas <- diff(sim$d) * n_poblacion

    emisiones <- bind_rows(emisiones, tibble(semana = semanas_tramo, muertes_emitidas = muertes_nuevas)) %>%
      arrange(semana)

    resultados_sim[[k]] <- sim %>%
      filter(semana <= semana_fin_k) %>%
      mutate(tramo = k, nombre_tramo = tramos_def$nombre[k])

    estado_actual <- sim %>% slice_tail(n = 1) %>% select(s, i, iv, v, q, r, d) %>% unlist()
  }

  list(emisiones = emisiones, simulacion = bind_rows(resultados_sim))
}

# Proceso iterativo: (A) ajustar beta -> (B) ajustar theta -> (A) volver a
# ajustar beta -> ... hasta que el error global deje de mejorar.
#
# (A) minimiza el error DENTRO DE CADA TRAMO y (B) minimiza el error GLOBAL:
# son dos objetivos distintos, asi que alternar entre ambos no garantiza que
# cada pasada mejore el ajuste global (puede quedar mejor por tramo y peor en
# conjunto). Por eso se guarda el mejor par (beta, theta) visto hasta ahora
# segun el error global, y el proceso corta apenas una pasada deja de
# mejorarlo (con la tolerancia de convergencia definida abajo) en vez de
# quedarse con el resultado de la ultima iteracion sin mas.
max_iteraciones_externas <- 5      # EDITAR si se quieren mas/menos pasadas
tolerancia_convergencia  <- 1e-6

theta_actual     <- theta_inicial
mejor_theta      <- theta_inicial
mejor_betas      <- NULL
mejor_error      <- Inf
historial_ajuste <- tibble(iteracion = integer(), theta = double(), error_global = double())

for (iteracion in seq_len(max_iteraciones_externas)) {

  ajuste_beta    <- ajustar_betas(theta_actual)          # (A) beta con theta fijo
  betas_actuales <- ajuste_beta$betas

  ajuste_theta <- optimize(
    evaluar_theta, interval = c(theta_min, theta_max), tol = 1e-4,
    betas_fijos = betas_actuales
  )                                                       # (B) theta con beta fijo

  historial_ajuste <- bind_rows(historial_ajuste, tibble(
    iteracion = iteracion, theta = ajuste_theta$minimum, error_global = ajuste_theta$objective
  ))

  if (ajuste_theta$objective < mejor_error - tolerancia_convergencia) {
    mejor_error  <- ajuste_theta$objective
    mejor_theta  <- ajuste_theta$minimum
    mejor_betas  <- betas_actuales
    theta_actual <- ajuste_theta$minimum
  } else {
    break
  }
}

# Se recalcula la simulacion final con el mejor par (beta, theta) encontrado
# (sin volver a optimizar beta), para que quede exactamente consistente con
# mejor_error y con lo que se reporta en el Bloque 4.
ajuste_final        <- ajustar_betas_con_beta_fijo(mejor_theta, mejor_betas)
betas_ajustados      <- mejor_betas
emisiones_global     <- ajuste_final$emisiones
simulacion_completa  <- ajuste_final$simulacion
theta                <- mejor_theta   # valor final de theta, ya calibrado

# ----------------------------------------------------------------------------
# BLOQUE 4 - RESULTADOS DEL MODELO
# ----------------------------------------------------------------------------

## 4.1 Tabla resumen por tramo: beta, R0, Re, tasa de contagio y mortalidad -
# La tasa de contagio (nuevas infecciones semanales / N) y la tasa de
# mortalidad (muertes semanales / N) se calculan directamente a partir de la
# simulacion y de las muertes emitidas por tramo.

tasas_por_tramo <- simulacion_completa %>%
  group_by(tramo) %>%
  summarise(s_prom = mean(s), v_prom = mean(v), .groups = "drop") %>%
  left_join(
    emisiones_global %>%
      left_join(simulacion_completa %>% select(semana, tramo), by = "semana") %>%
      group_by(tramo) %>%
      summarise(tasa_mortalidad_semanal = mean(muertes_emitidas) / n_poblacion, .groups = "drop"),
    by = "tramo"
  )

resumen_tramos <- tramos_def %>%
  mutate(
    beta = betas_ajustados,
    theta_vigente = if_else(nombre == "Pre-intervencion", 0, theta)
  ) %>%
  left_join(tasas_por_tramo, by = "tramo") %>%
  mutate(
    r0 = calcular_r0(beta, theta_vigente),
    re = calcular_re(beta, theta_vigente, s_prom, v_prom),
    tasa_contagio_semanal = beta * s_prom
  ) %>%
  select(tramo, nombre, semana_ini, semana_fin, beta, r0, re,
         tasa_contagio_semanal, tasa_mortalidad_semanal)

mse_global <- {
  comparacion_global <- emisiones_global %>%
    mutate(muertes_comparables = dplyr::lag(muertes_emitidas, n = rezago_semanas, default = 0)) %>%
    left_join(datos_owid_semanal %>% select(semana, muertes_owid), by = "semana")
  tibble(
    # msre es la funcion de perdida que efectivamente se minimizo (debe
    # coincidir con mejor_error). mse_nivel y mse_log son solo diagnostico.
    msre      = calcular_error_global(emisiones_global),
    mse_nivel = mean((comparacion_global$muertes_comparables - comparacion_global$muertes_owid)^2, na.rm = TRUE),
    mse_log   = mean((log1p(comparacion_global$muertes_comparables) - log1p(comparacion_global$muertes_owid))^2, na.rm = TRUE)
  )
}

cat("Theta calibrado:", theta, " (punto de partida:", theta_inicial, ")\n")
print(historial_ajuste)
print(resumen_tramos)
print(mse_global)

## 4.2 Grafico: evolucion semanal de los compartimentos ---------------------

grafico_compartimentos <- simulacion_completa %>%
  mutate(fecha = fecha_inicio_serie + (semana - 1) * 7) %>%
  pivot_longer(cols = c(s, i, iv, v, q, r, d), names_to = "compartimento", values_to = "fraccion") %>%
  mutate(compartimento = factor(compartimento,
                                  levels = c("s", "i", "iv", "v", "q", "r", "d"),
                                  labels = c("Susceptibles", "Infectados", "Infectados vacunados",
                                              "Vacunados", "Cuarentena", "Recuperados", "Fallecidos"))) %>%
  ggplot(aes(x = fecha, y = fraccion, color = compartimento)) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(title = "Modelo SIVQRD - Argentina", subtitle = "Evolucion semanal de los compartimentos",
       x = NULL, y = "Fraccion de la poblacion", color = NULL) +
  theme_minimal()

## 4.3 Grafico: muertes semanales reales vs. predichas por el modelo --------

fechas_quiebre <- c(fecha_aspo, fecha_dispo, fecha_vacunacion, fecha_omicron)

grafico_muertes <- emisiones_global %>%
  mutate(muertes_modelo = dplyr::lag(muertes_emitidas, n = rezago_semanas, default = 0)) %>%
  left_join(datos_owid_semanal %>% select(semana, muertes_owid), by = "semana") %>%
  mutate(fecha = fecha_inicio_serie + (semana - 1) * 7) %>%
  select(fecha, muertes_modelo, muertes_owid) %>%
  pivot_longer(cols = c(muertes_modelo, muertes_owid), names_to = "serie", values_to = "muertes") %>%
  mutate(serie = recode(serie, muertes_modelo = "Modelo SIVQRD", muertes_owid = "OWID (observado)")) %>%
  ggplot(aes(x = fecha, y = muertes, color = serie)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = fechas_quiebre, linetype = "dashed", color = "grey50") +
  scale_x_date(date_labels = "%b-%Y") +
  labs(title = "Muertes semanales: modelo vs. observado", x = NULL, y = "Muertes por semana", color = NULL) +
  theme_minimal()

grafico_compartimentos
grafico_muertes
