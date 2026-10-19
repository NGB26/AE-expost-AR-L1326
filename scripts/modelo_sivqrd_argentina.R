# ============================================================================
# MODELO EPIDEMIOLOGICO SIVQRD - ARGENTINA COVID-19
# Insumo para el analisis economico ex post del PCR AR-L1326 (BID)
#
# Adapta el libro "SIRD_Argentina_PCR.xlsx" (modelo SIRD simple, ajuste semanal
# por formulas de planilla) a un modelo SIVQRD completo, resuelto e integrado
# en R. Las ecuaciones siguen la numeracion de "Modelo_SIR_para_PCR.docx".
#
# Fuente de calibracion: muertes semanales de Our World in Data (OWID), igual
# que en la planilla original. La serie oficial de datos.gob.ar se importa
# como serie de contraste/validacion (no interviene en el ajuste de beta).
#
# Metodo de ajuste: secuencial por tramo. Para cada uno de los 5 tramos de
# beta se corre el modelo con el estado final del tramo anterior como
# condicion inicial, y se busca el beta que minimiza el error dentro de ESE
# tramo unicamente. Los parametros restantes (gamma, delta, theta, sigma,
# beta_v, gamma_v, delta_v) son fijos y no se calibran.
#
# FUNCION DE PERDIDA: MSE(log1p(muertes modelo), log1p(muertes OWID)) sobre
# la agregacion semanal. El MSE en nivel se sigue reportando en el Bloque 4
# como diagnostico, pero no interviene en el ajuste de beta.
#
# NOTA SOBRE DATOS.GOB.AR: esta sesion no tuvo acceso de red al host
# datos.gob.ar para confirmar en vivo los nombres de columnas del dataset
# "Casos registrados de COVID-19 en la Republica Argentina"
# (https://www.datos.gob.ar/dataset/salud-covid-19-casos-registrados-republica-argentina).
# Los nombres usados abajo (clasificacion_resumen, fecha_apertura,
# fecha_fallecimiento, fallecido) son los del diccionario de datos publicado
# historicamente por el Ministerio de Salud para ese dataset; verificar contra
# el diccionario vigente antes de la primera corrida.
# ============================================================================

library(tidyverse)   # dplyr, tidyr, purrr, ggplot2
library(readxl)      # lectura de Parametros y Tramos_Beta desde el Excel
library(lubridate)   # manejo de fechas
library(deSolve)     # integracion numerica de las ecuaciones diferenciales
library(scales)      # formato de ejes en los graficos

# ----------------------------------------------------------------------------
# BLOQUE 0 - RUTAS Y CONSTANTES DE SERIE
# ----------------------------------------------------------------------------

ruta_excel_parametros <- "data/raw/SIRD_Argentina_PCR.xlsx"

fecha_inicio_serie <- as.Date("2020-03-08")   # semana 1, igual que Datos_OWID
fecha_fin_serie    <- as.Date("2022-06-05")   # ultima semana cubierta por el Excel

# Numero de semana (entero) que CONTIENE una fecha dada, con semana 1 =
# fecha_inicio_serie. Se usa tanto para agregar los datos diarios a semanal
# como para ubicar las fechas de quiebre de politica dentro de esa misma
# grilla semanal (mismo criterio que la columna "Semana ini" de Tramos_Beta).
semana_de_fecha <- function(fecha) as.integer(floor(as.numeric(fecha - fecha_inicio_serie) / 7)) + 1L

# ----------------------------------------------------------------------------
# BLOQUE 1 - IMPORTACION DE DATOS
# ----------------------------------------------------------------------------

## 1.1 Parametros fijos y fechas de tramo, leidos desde el Excel adjunto -----
# Se lee por ETIQUETA (no por celda fija) para poder seguir editando el Excel
# sin romper el script, igual que el libro original invita a hacer en la hoja
# "Parametros" ("AQUI se edita").

obtener_parametro <- function(tabla, etiqueta) {
  fila <- tabla %>% filter(parametro == etiqueta)
  if (nrow(fila) != 1) {
    stop(sprintf("Parametro '%s' no encontrado o ambiguo (%d coincidencias) en la hoja Parametros.",
                  etiqueta, nrow(fila)))
  }
  fila$valor[[1]]
}

parametros_excel <- read_excel(ruta_excel_parametros, sheet = "Parametros", range = "B4:G30") %>%
  set_names(c("parametro", "valor", "min", "max", "unidad", "fuente")) %>%
  filter(!is.na(valor))

tramos_excel <- read_excel(ruta_excel_parametros, sheet = "Tramos_Beta", range = "B5:I14") %>%
  set_names(c("id", "hito", "fecha_inicio", "semana_ini", "semana_fin", "beta_excel", "r0_excel", "fuente"))

# Fechas de quiebre de politica: ASPO, DISPO, inicio de vacunacion y Omicron.
# Son las mismas 4 fechas que ya estan documentadas en Tramos_Beta; aqui se
# usan para construir los 5 TRAMOS DE BETA del modelo SIVQRD (en vez de los 9
# tramos de la planilla original).
fecha_aspo       <- tramos_excel %>% filter(hito == "Cierre fronteras + ASPO (Dto 297/20)") %>% pull(fecha_inicio) %>% as.Date()
fecha_dispo      <- tramos_excel %>% filter(hito == "DISPO / flexibilizacion") %>% pull(fecha_inicio) %>% as.Date()
fecha_vacunacion <- tramos_excel %>% filter(hito == "Inicio vacunacion") %>% pull(fecha_inicio) %>% as.Date()
fecha_omicron    <- tramos_excel %>% filter(hito == "Omicron") %>% pull(fecha_inicio) %>% as.Date()

## 1.2 Our World in Data (OWID) - serie de calibracion ----------------------
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
    fecha                   = min(date),
    casos_owid               = sum(new_cases, na.rm = TRUE),
    muertes_owid             = sum(new_deaths, na.rm = TRUE),
    tests_owid                = sum(new_tests, na.rm = TRUE),
    positividad_owid          = mean(positive_rate, na.rm = TRUE),
    vacunados_completos_owid = last(na.omit(people_fully_vaccinated)),
    .groups = "drop"
  ) %>%
  arrange(semana)

## 1.3 datos.gob.ar - serie oficial argentina, para contraste ---------------
# Dataset: "Casos registrados de COVID-19 en la Republica Argentina"
# https://www.datos.gob.ar/dataset/salud-covid-19-casos-registrados-republica-argentina
# Se usa la API CKAN para resolver la URL del recurso vigente en lugar de un
# link fijo, porque la planilla original advierte que la frecuencia y el
# nombre del archivo cambiaron durante la pandemia.

url_api_dataset <- "https://datos.gob.ar/api/3/action/package_show?id=salud-covid-19-casos-registrados-republica-argentina"

resolver_url_recurso_datos_gob_ar <- function() {
  metadata <- jsonlite::fromJSON(url_api_dataset)
  recursos <- metadata$result$resources
  # Se prioriza el recurso identificado como archivo historico de "Casos".
  candidato <- recursos %>%
    as_tibble() %>%
    filter(str_detect(tolower(name), "caso")) %>%
    slice(1)
  candidato$url
}

url_datos_gob_ar <- resolver_url_recurso_datos_gob_ar()

datos_argentina_oficial_diario <- read_csv(
  url_datos_gob_ar,
  col_select = c(clasificacion_resumen, fecha_apertura, fallecido, fecha_fallecimiento),
  show_col_types = FALSE,
  locale = locale(encoding = "UTF-8")
) %>%
  filter(clasificacion_resumen == "Confirmado")

casos_oficial_semanal <- datos_argentina_oficial_diario %>%
  mutate(fecha_apertura = as.Date(fecha_apertura)) %>%
  filter(fecha_apertura >= fecha_inicio_serie, fecha_apertura <= fecha_fin_serie) %>%
  mutate(semana = semana_de_fecha(fecha_apertura)) %>%
  count(semana, name = "casos_oficial_ar")

muertes_oficial_semanal <- datos_argentina_oficial_diario %>%
  filter(fallecido == "SI") %>%
  mutate(fecha_fallecimiento = as.Date(fecha_fallecimiento)) %>%
  filter(fecha_fallecimiento >= fecha_inicio_serie, fecha_fallecimiento <= fecha_fin_serie) %>%
  mutate(semana = semana_de_fecha(fecha_fallecimiento)) %>%
  count(semana, name = "muertes_oficial_ar")

## 1.4 Serie semanal unica, para calibracion y grafico de contraste ---------

datos_semanales <- datos_owid_semanal %>%
  left_join(casos_oficial_semanal, by = "semana") %>%
  left_join(muertes_oficial_semanal, by = "semana") %>%
  mutate(across(c(casos_oficial_ar, muertes_oficial_ar), ~ replace_na(.x, 0)))

# ----------------------------------------------------------------------------
# BLOQUE 2 - FUNCIONES Y MODELO SIVQRD (ecuaciones diferenciales)
# ----------------------------------------------------------------------------

## 2.1 Parametros fijos del modelo -------------------------------------------
# N, S0/N, I0, duracion infecciosa e IFR vienen del Excel. theta, sigma,
# beta_v (via ve_contagio), gamma_v y delta_v no estan en el Excel original
# (el modelo SIRD simple no los usa): quedan como constantes "a completar con
# la literatura", exactamente como las demas filas en blanco de la hoja
# Parametros. NO son valores calibrados.

n_poblacion     <- obtener_parametro(parametros_excel, "Poblacion Argentina (N)")
s0_frac         <- obtener_parametro(parametros_excel, "Fraccion susceptible inicial (S0/N)")
i0_personas     <- obtener_parametro(parametros_excel, "Infecciosos iniciales (I0)")
dias_infecciosos <- obtener_parametro(parametros_excel, "Duracion del periodo infeccioso (D)")
ifr_base        <- obtener_parametro(parametros_excel, "IFR base pre-vacunacion")
rezago_semanas  <- obtener_parametro(parametros_excel, "Rezago infeccion -> muerte")
ve_muerte       <- obtener_parametro(parametros_excel, "Efectividad vacunal contra muerte")

# Parametros adicionales del SIVQRD: EDITAR con la fuente bibliografica
# correspondiente antes de usar el modelo para resultados finales.
theta       <- NA_real_   # EDITAR: tasa semanal de entrada a cuarentena de los infectados (I -> Q)
sigma       <- NA_real_   # EDITAR: tasa semanal de vacunacion de susceptibles (S -> V)
ve_contagio <- NA_real_   # EDITAR: efectividad vacunal contra el contagio (define beta_v = beta*(1-ve_contagio))
dias_infecciosos_vac <- dias_infecciosos  # EDITAR si la literatura da una duracion distinta para Iv

# gamma y delta se derivan de la tasa de salida total del compartimento I
# (gamma_total = 7/D, en semanas) partida segun el IFR, igual que hace la
# planilla (hoja Modelo_SIRD, columna F) pero ahora como dos parametros
# independientes, tal como los define el modelo SIVQRD (ecuaciones 23-29).
gamma_total       <- 7 / dias_infecciosos
delta             <- ifr_base * gamma_total        # tasa de mortalidad de I  (parametro delta, ec. 24 y 29)
gamma             <- (1 - ifr_base) * gamma_total   # tasa de recuperacion de I (parametro gamma, ec. 24 y 28)

ifr_vacunados     <- ifr_base * (1 - ve_muerte)
gamma_total_vac   <- 7 / dias_infecciosos_vac
delta_v           <- ifr_vacunados * gamma_total_vac   # tasa de mortalidad de Iv (ec. 25 y 29)
gamma_v           <- (1 - ifr_vacunados) * gamma_total_vac  # tasa de recuperacion de Iv (ec. 25 y 28)

## 2.2 Tramos de beta (5 tramos, segun hitos de politica) -------------------
# semana_de_fecha() ya esta definida en el Bloque 0 (misma grilla semanal que
# los datos). cortes_semana[k] = primera semana (entera) del tramo k.

cortes_semana <- c(1,
                    semana_de_fecha(fecha_aspo),
                    semana_de_fecha(fecha_dispo),
                    semana_de_fecha(fecha_vacunacion),
                    semana_de_fecha(fecha_omicron),
                    max(datos_semanales$semana) + 1)

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

parametros_fijos <- list(
  theta = theta, sigma = sigma, ve_contagio = ve_contagio,
  gamma = gamma, delta = delta, gamma_v = gamma_v, delta_v = delta_v,
  semana_inicio_aspo = semana_de_fecha(fecha_aspo),
  semana_inicio_vacunacion = semana_de_fecha(fecha_vacunacion)
)

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
# BLOQUE 3 - AJUSTE Y CALIBRACION DEL MODELO
# ----------------------------------------------------------------------------
# Variable de ajuste: beta, un valor por tramo (5 en total). El resto de los
# parametros es fijo (Bloque 2). Objetivo por tramo: minimizar
#   MSE(log1p(muertes modelo con rezago), log1p(muertes OWID))
# sobre la agregacion SEMANAL (no diaria). El ajuste es secuencial: beta_1 se
# fija usando solo las semanas del tramo 1; el estado final de ese tramo pasa
# como condicion inicial del tramo 2, y asi sucesivamente.

estado_inicial <- c(
  s = s0_frac - i0_personas / n_poblacion,
  i = i0_personas / n_poblacion,
  iv = 0, v = 0, q = 0, r = 0, d = 0
)

estado_actual  <- estado_inicial
emisiones_global <- tibble(semana = numeric(), muertes_emitidas = numeric())
resultados_sim   <- list()
betas_ajustados  <- numeric(nrow(tramos_def))

for (k in seq_len(nrow(tramos_def))) {

  semana_ini_k  <- tramos_def$semana_ini[k]
  semana_fin_k  <- tramos_def$semana_fin[k]
  semanas_tramo <- seq(semana_ini_k, semana_fin_k, by = 1)

  calcular_error_tramo <- function(beta_k) {
    sim <- simular_tramo(estado_actual, semana_ini_k, semana_fin_k, beta_k, parametros_fijos)

    # Muertes nuevas emitidas por el modelo en la semana w = incremento del
    # compartimento d (acumulado) entre el inicio de w y el inicio de w+1,
    # multiplicado por la poblacion. sim tiene semana_fin_k - semana_ini_k + 2
    # filas (incluye el punto extra al final), por eso diff() da exactamente
    # un valor por cada semana de semanas_tramo.
    muertes_nuevas <- diff(sim$d) * n_poblacion
    emisiones_candidatas <- tibble(semana = semanas_tramo, muertes_emitidas = muertes_nuevas)

    emisiones_comparables <- bind_rows(emisiones_global, emisiones_candidatas) %>%
      arrange(semana) %>%
      mutate(muertes_comparables = dplyr::lag(muertes_emitidas, n = rezago_semanas, default = 0)) %>%
      filter(semana %in% semanas_tramo) %>%
      left_join(datos_semanales %>% select(semana, muertes_owid), by = "semana")

    mean((log1p(emisiones_comparables$muertes_comparables) - log1p(emisiones_comparables$muertes_owid))^2, na.rm = TRUE)
  }

  ajuste <- optimize(calcular_error_tramo, interval = c(0.01, 2), tol = 1e-5)
  betas_ajustados[k] <- ajuste$minimum

  sim_final <- simular_tramo(estado_actual, semana_ini_k, semana_fin_k, betas_ajustados[k], parametros_fijos)
  muertes_nuevas_final <- diff(sim_final$d) * n_poblacion

  emisiones_global <- bind_rows(emisiones_global, tibble(semana = semanas_tramo, muertes_emitidas = muertes_nuevas_final)) %>%
    arrange(semana)

  # Se descarta el punto extra (semana_fin_k + 1) de los resultados que se
  # grafican: es solo el estado inicial del tramo siguiente.
  resultados_sim[[k]] <- sim_final %>%
    filter(semana <= semana_fin_k) %>%
    mutate(tramo = k, nombre_tramo = tramos_def$nombre[k])

  estado_actual <- sim_final %>%
    slice_tail(n = 1) %>%
    select(s, i, iv, v, q, r, d) %>%
    unlist()
}

simulacion_completa <- bind_rows(resultados_sim)

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
    left_join(datos_semanales %>% select(semana, muertes_owid), by = "semana")
  tibble(
    mse_nivel = mean((comparacion_global$muertes_comparables - comparacion_global$muertes_owid)^2, na.rm = TRUE),
    mse_log   = mean((log1p(comparacion_global$muertes_comparables) - log1p(comparacion_global$muertes_owid))^2, na.rm = TRUE)
  )
}

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
  left_join(datos_semanales %>% select(semana, muertes_owid), by = "semana") %>%
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
