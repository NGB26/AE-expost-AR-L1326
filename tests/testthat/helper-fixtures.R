# ==============================================================================
# Fixtures sinteticas, deterministicas, SIN dependencia de internet. Se usan
# en todos los tests de tests/testthat/.
# ==============================================================================

#' Config de prueba: serie corta (20 semanas), 2 tramos, parametros simples
#' (no son valores de literatura, son valores de prueba redondos).
fixture_valores <- function() {
  list(
    fecha_inicio_serie = as.Date("2021-01-03"),
    fecha_fin_serie    = as.Date("2021-05-23"),  # 20 semanas
    fecha_aspo         = as.Date("2021-02-14"),  # separa tramo 1 de tramo 2 (semana 7)
    fecha_dispo        = as.Date("2021-02-14"),
    fecha_vacunacion   = as.Date("2021-03-21"),  # semana 12
    fecha_omicron      = as.Date("2021-05-23"),  # ultima semana: sin tramo 5 util
    s0_frac = 1, i0_personas = 50,
    dias_infecciosos = 11, dias_infecciosos_vac = 8.7,
    ifr_base = 0.01, ve_muerte = 0.7, ve_contagio = 0.5,
    rezago_semanas = 2,
    sigma = 0.05,
    theta = 0.3,
    beta_min = 0.01, beta_max = 2,
    maxit_optimizador = 200
  )
}

#' n_poblacion de prueba: chico, para que las simulaciones sean rapidas y
#' las cantidades absolutas (muertes, etc.) queden en un rango razonable.
fixture_n_poblacion <- function() 1e6

#' tramos_def de prueba: 2 tramos sobre una serie de 20 semanas.
fixture_tramos_def <- function() {
  tibble::tibble(
    tramo = 1:2, nombre = c("Tramo1", "Tramo2"),
    semana_ini = c(1, 8), semana_fin = c(7, 20)
  )
}

#' Datos semanales sinteticos: una curva de muertes determinista (sube y
#' baja, como una ola), sin ceros, para que MSRE-log este siempre definido.
fixture_datos_semanal <- function(n_semanas = 20, n_poblacion = fixture_n_poblacion()) {
  semana <- 1:n_semanas
  pico <- 10
  # Curva log-normal discreta simple: sube, pica en pico, baja. +1 para que
  # nunca de exactamente 0 (MSRE-log no esta definido en 0).
  muertes <- round(50 * exp(-((semana - pico)^2) / (2 * 4^2)) + 1)
  vacunados <- pmin(round(cumsum(ifelse(semana >= 12, 2000, 0))), n_poblacion * 0.5)

  tibble::tibble(
    semana = semana,
    fecha = as.Date("2021-01-03") + (semana - 1) * 7,
    dias_observados = 7L,
    casos_owid = muertes * 20L,
    muertes_owid = muertes,
    tests_owid = muertes * 200L,
    positividad_owid = 0.1,
    vacunados_completos_owid = vacunados,
    poblacion_owid = n_poblacion,
    semana_completa = TRUE
  )
}

#' Parametros biologicos de prueba, ya derivados.
fixture_biologicos <- function(valores = fixture_valores()) {
  derivar_parametros_biologicos(valores)
}

#' Estado inicial de prueba.
fixture_estado_inicial <- function(valores = fixture_valores(), n_poblacion = fixture_n_poblacion()) {
  construir_estado_inicial(valores$s0_frac, valores$i0_personas, n_poblacion)
}

#' parametros_fijos de prueba, para un theta dado.
fixture_parametros_fijos <- function(theta_valor, valores = fixture_valores(), biologicos = fixture_biologicos(valores)) {
  construir_parametros_fijos(
    theta_valor, valores$sigma, valores$ve_contagio, biologicos$gamma, biologicos$delta,
    biologicos$gamma_v, biologicos$delta_v,
    semana_de_fecha(valores$fecha_aspo, valores$fecha_inicio_serie),
    semana_de_fecha(valores$fecha_vacunacion, valores$fecha_inicio_serie)
  )
}
