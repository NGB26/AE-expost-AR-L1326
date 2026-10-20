# ==============================================================================
# CALIBRATE - Ajuste de los beta por segmento (SOLO beta es variable libre)
# ==============================================================================
# Especificacion (confirmada con el usuario):
#   * Hay 5 segmentos definidos por 4 fechas de corte HARDCODEADAS: comienzo de
#     ASPO, comienzo de DISPO, comienzo de vacunacion, Omicron (ver R/config.R
#     y scripts/run.R). Cada segmento tiene UN beta propio.
#   * La UNICA variable de ajuste es beta (un valor por segmento). theta, IFR,
#     gamma, delta y la efectividad vacunal (ve_contagio) estan FIJOS (de
#     literatura); sigma se calcula de los datos (R/data.R).
#   * La funcion de perdida es el RMSE de log(1+muertes) del modelo vs. OWID
#     (R/metrics.R::calcular_rmse_log).
#   * El metodo PRINCIPAL (calibrar_beta_global) minimiza ese RMSE de forma
#     GLOBAL y SIMULTANEA sobre toda la serie, moviendo los 5 beta a la vez
#     (no segmento por segmento). calibrar_secuencial() queda como comparacion.
#
# Todas las funciones reciben sus insumos explicitos (estado_inicial,
# tramos_def, parametros biologicos, theta, datos_semanal, n_poblacion, cotas)
# y no leen variables globales.
# ==============================================================================

library(dplyr)
library(purrr)
library(tidyr)

#' Arma la lista de parametros fijos del modelo para un theta dado. Durante la
#' calibracion NADA de esto cambia (theta es fijo; sigma, ve_contagio, gamma,
#' delta, gamma_v, delta_v tampoco): lo unico que se mueve son los beta, que
#' NO van en esta lista sino como argumento aparte de simular_serie().
construir_parametros_fijos <- function(theta_valor, sigma, ve_contagio, gamma, delta, gamma_v, delta_v,
                                        semana_inicio_aspo, semana_inicio_vacunacion) {
  list(
    theta = theta_valor, sigma = sigma, ve_contagio = ve_contagio,
    gamma = gamma, delta = delta, gamma_v = gamma_v, delta_v = delta_v,
    semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion
  )
}

#' Estado inicial del modelo (s, i, iv, v, q, r, d, inf_acum), a partir de la
#' fraccion susceptible inicial y los infecciosos iniciales.
construir_estado_inicial <- function(s0_frac, i0_personas, n_poblacion) {
  c(
    s = s0_frac - i0_personas / n_poblacion,
    i = i0_personas / n_poblacion,
    iv = 0, v = 0, q = 0, r = 0, d = 0, inf_acum = 0
  )
}

# ------------------------------------------------------------------------------
# METODO PRINCIPAL: ajuste CONJUNTO y GLOBAL de los beta (theta fijo)
# ------------------------------------------------------------------------------

#' Calibra los beta (uno por segmento) de forma CONJUNTA y SIMULTANEA,
#' minimizando el RMSE de log(1+muertes) de TODA la serie de una sola vez.
#'
#' Este es el metodo que pide la especificacion: a diferencia de un ajuste
#' secuencial (cada beta minimiza el error de SU segmento, con el estado
#' heredado del segmento anterior), aca optim() mueve los 5 beta a la vez para
#' que el error GLOBAL sea minimo. Importa porque el estado al final de un
#' segmento (cuantos susceptibles quedan, cuantos vacunados) depende de los
#' beta de TODOS los segmentos anteriores: el optimo global puede sacrificar
#' un poco el ajuste de un segmento para mejorar el de los que vienen despues,
#' cosa que el secuencial no puede hacer.
#'
#' theta NO es variable de ajuste: entra fijo via `theta`. Idem sigma,
#' ve_contagio y los parametros biologicos (gamma, delta, gamma_v, delta_v).
#'
#' @param betas_iniciales Punto de partida de optim() (vector de largo
#'   nrow(tramos_def)). Si es NULL, arranca en el centro-bajo del rango. En el
#'   pipeline se lo pasa el resultado del ajuste secuencial, que es un punto de
#'   partida razonable y barato de calcular (practica estandar).
#' @param maxit Tope de iteraciones de optim() (L-BFGS-B, respeta cotas).
#' @return Lista con betas, theta (el fijo, por trazabilidad), emisiones,
#'   incidencia, simulacion, estado_final, mejor_error (RMSE-log global),
#'   convergencia (codigo de optim(): 0 = ok), n_evaluaciones y metodo.
calibrar_beta_global <- function(estado_inicial, tramos_def, biologicos, theta, sigma, ve_contagio,
                                 semana_inicio_aspo, semana_inicio_vacunacion,
                                 datos_semanal, n_poblacion, rezago_semanas,
                                 beta_min, beta_max, betas_iniciales = NULL, maxit = 500) {

  n_tramos <- nrow(tramos_def)
  parametros_fijos <- construir_parametros_fijos(
    theta, sigma, ve_contagio, biologicos$gamma, biologicos$delta,
    biologicos$gamma_v, biologicos$delta_v, semana_inicio_aspo, semana_inicio_vacunacion
  )

  n_evaluaciones <- 0
  # PENALIZACION_ODE: si para un vector de beta el ODE no integra (region de
  # parametros numericamente infactible; ver R/simulate.R::simular_tramo), se
  # devuelve un valor finito y grande en vez de dejar que el error corte la
  # calibracion. Asi optim() se aleja de esa region en vez de morir ahi.
  PENALIZACION_ODE <- 1e6
  objetivo_global <- function(betas) {
    n_evaluaciones <<- n_evaluaciones + 1
    err <- tryCatch({
      sim <- simular_serie(estado_inicial, tramos_def, betas, parametros_fijos)
      # RMSE-log sobre TODA la serie (semanas_objetivo = NULL): el "global y
      # simultaneo" de la especificacion.
      calcular_rmse_log(sim$emisiones, datos_semanal, n_poblacion, rezago_semanas)
    }, error = function(e) NA_real_)
    if (!is.finite(err)) PENALIZACION_ODE else err
  }

  if (is.null(betas_iniciales)) {
    betas_iniciales <- rep(mean(c(beta_min, beta_max)) / 2, n_tramos)
  }
  stopifnot(length(betas_iniciales) == n_tramos)

  ajuste <- optim(
    par = betas_iniciales, fn = objetivo_global, method = "L-BFGS-B",
    lower = rep(beta_min, n_tramos), upper = rep(beta_max, n_tramos),
    control = list(maxit = maxit, factr = 1e7)
  )

  betas_finales   <- ajuste$par
  resultado_final <- simular_serie(estado_inicial, tramos_def, betas_finales, parametros_fijos)

  list(
    betas = betas_finales, theta = theta,
    emisiones = resultado_final$emisiones, incidencia = resultado_final$incidencia,
    simulacion = resultado_final$simulacion, estado_final = resultado_final$estado_final,
    mejor_error = ajuste$value, convergencia = ajuste$convergence,
    n_evaluaciones = ajuste$counts[["function"]], metodo = "beta_global_simultaneo"
  )
}

#' Multi-start del ajuste global: corre calibrar_beta_global() desde varios
#' puntos de partida (uno "central" determinista y n_aleatorios sorteados
#' uniformemente en [beta_min, beta_max]) y se queda con el de menor RMSE-log.
#' Si todos convergen al mismo optimo, el resultado es robusto al arranque; si
#' no, el optimo global depende del punto de partida y hay que reportarlo.
#'
#' @return Lista con `mejor` (el resultado de calibrar_beta_global del mejor
#'   arranque) y `tabla` (tibble: arranque, error, betas en list-col) para
#'   auditar la dispersion entre starts.
calibrar_beta_global_multistart <- function(estado_inicial, tramos_def, biologicos, theta, sigma, ve_contagio,
                                             semana_inicio_aspo, semana_inicio_vacunacion,
                                             datos_semanal, n_poblacion, rezago_semanas,
                                             beta_min, beta_max, betas_centrales = NULL,
                                             n_aleatorios = 4, maxit = 500, semilla = 1) {
  set.seed(semilla)
  n_tramos <- nrow(tramos_def)
  if (is.null(betas_centrales)) betas_centrales <- rep(mean(c(beta_min, beta_max)) / 2, n_tramos)

  arranques <- c(
    list(betas_centrales),
    purrr::map(seq_len(n_aleatorios), function(i) runif(n_tramos, beta_min, beta_max))
  )

  resultados <- purrr::imap(arranques, function(b0, idx) {
    r <- calibrar_beta_global(
      estado_inicial, tramos_def, biologicos, theta, sigma, ve_contagio,
      semana_inicio_aspo, semana_inicio_vacunacion, datos_semanal, n_poblacion, rezago_semanas,
      beta_min, beta_max, betas_iniciales = b0, maxit = maxit
    )
    list(idx = idx, resultado = r)
  })

  tabla <- purrr::map_dfr(resultados, function(x) tibble(
    arranque = if (x$idx == 1) "central" else sprintf("aleatorio_%d", x$idx - 1),
    error = x$resultado$mejor_error, betas = list(x$resultado$betas)
  ))

  mejor_idx <- which.min(purrr::map_dbl(resultados, ~ .x$resultado$mejor_error))
  list(mejor = resultados[[mejor_idx]]$resultado, tabla = tabla)
}

# ------------------------------------------------------------------------------
# COMPARACION: ajuste secuencial por segmento (theta fijo), mismo loss
# ------------------------------------------------------------------------------

#' Ajuste SECUENCIAL: cada beta minimiza el RMSE-log de SU propio segmento, en
#' orden, heredando el estado del segmento anterior (ver
#' R/simulate.R::ajustar_betas_secuencial). Se incluye SOLO como punto de
#' comparacion contra calibrar_beta_global(): sirve para mostrar cuanto se
#' gana (si algo) resolviendo el problema global en vez de segmento a segmento,
#' y es un buen punto de partida barato para el optimizador global. theta fijo.
#'
#' @return Lista con betas, theta, emisiones, incidencia, simulacion,
#'   estado_final, mejor_error (RMSE-log GLOBAL de la serie resultante) y metodo.
calibrar_secuencial <- function(estado_inicial, tramos_def, biologicos, theta, sigma, ve_contagio,
                                semana_inicio_aspo, semana_inicio_vacunacion,
                                datos_semanal, n_poblacion, rezago_semanas, beta_min, beta_max) {

  parametros_fijos <- construir_parametros_fijos(
    theta, sigma, ve_contagio, biologicos$gamma, biologicos$delta,
    biologicos$gamma_v, biologicos$delta_v, semana_inicio_aspo, semana_inicio_vacunacion
  )

  funcion_error <- function(emisiones, semanas_objetivo) {
    calcular_rmse_log(emisiones, datos_semanal, n_poblacion, rezago_semanas, semanas_objetivo)
  }

  r <- ajustar_betas_secuencial(estado_inicial, tramos_def, parametros_fijos, funcion_error, beta_min, beta_max)
  error_global <- calcular_rmse_log(r$emisiones, datos_semanal, n_poblacion, rezago_semanas)

  list(
    betas = r$betas, theta = theta,
    emisiones = r$emisiones, incidencia = r$incidencia, simulacion = r$simulacion,
    estado_final = r$estado_final, mejor_error = error_global, metodo = "secuencial_beta"
  )
}

# ------------------------------------------------------------------------------
# Sensibilidad al theta FIJO (theta no se calibra: hay que mostrar cuanto
# dependen beta y el ajuste del valor que se le asigno)
# ------------------------------------------------------------------------------

#' Como theta esta FIJO (no se calibra), su valor es un supuesto. Esta funcion
#' barre una grilla de theta: para cada theta candidato recalibra los 5 beta de
#' forma global y reporta el RMSE-log alcanzado y los beta resultantes. Si el
#' error y los beta cambian poco en un rango amplio de theta, el resultado es
#' robusto al supuesto; si cambian mucho, theta importa y conviene justificar
#' mejor su valor (o, si el usuario lo habilita, pasarlo a variable libre).
#'
#' @param thetas Vector de valores de theta a probar.
#' @return Tibble(theta, rmse_log, beta_1..beta_k, convergencia).
analisis_sensibilidad_theta <- function(thetas, estado_inicial, tramos_def, biologicos, sigma, ve_contagio,
                                        semana_inicio_aspo, semana_inicio_vacunacion,
                                        datos_semanal, n_poblacion, rezago_semanas,
                                        beta_min, beta_max, betas_iniciales = NULL, maxit = 300) {
  purrr::map_dfr(thetas, function(theta_val) {
    r <- calibrar_beta_global(
      estado_inicial, tramos_def, biologicos, theta_val, sigma, ve_contagio,
      semana_inicio_aspo, semana_inicio_vacunacion, datos_semanal, n_poblacion, rezago_semanas,
      beta_min, beta_max, betas_iniciales = betas_iniciales, maxit = maxit
    )
    betas_nombrados <- setNames(as.list(r$betas), paste0("beta_", seq_along(r$betas)))
    tibble(theta = theta_val, rmse_log = r$mejor_error, convergencia = r$convergencia) %>%
      bind_cols(as_tibble(betas_nombrados))
  })
}

# ------------------------------------------------------------------------------
# Sensibilidad OAT (one-at-a-time) a los parametros biologicos/de literatura
# ------------------------------------------------------------------------------

#' Perturba UN parametro de literatura a la vez (dias_infecciosos, ifr_base,
#' rezago_semanas, ve_muerte, ve_contagio, dias_infecciosos_vac o theta) a un
#' valor alternativo, re-deriva lo que dependa de el (gamma/delta via
#' derivar_parametros_biologicos) y vuelve a correr la calibracion GLOBAL de
#' los beta, para ver cuanto cambian los beta, R0 y el error. No es un analisis
#' de sensibilidad formal (Sobol, FAST): es un barrido OAT simple, mas rapido y
#' facil de interpretar para un informe tecnico. Documentar las conclusiones
#' con esa cautela.
#'
#' @param cfg_valores Lista plana de aplanar_config(cfg) (incluye theta fijo).
#' @param parametro Nombre del parametro a perturbar (debe existir en cfg_valores).
#' @param valores_alternativos Valores alternativos a probar (ademas del base).
#' @return Tibble(parametro, valor, rmse_log, beta_1..beta_k).
analisis_sensibilidad <- function(cfg_valores, parametro, valores_alternativos,
                                  tramos_def, datos_semanal, n_poblacion,
                                  semana_inicio_aspo, semana_inicio_vacunacion,
                                  beta_min, beta_max, maxit = 300) {

  valores_a_probar <- unique(c(cfg_valores[[parametro]], valores_alternativos))

  purrr::map_dfr(valores_a_probar, function(valor_alt) {
    cfg_mod <- cfg_valores
    cfg_mod[[parametro]] <- valor_alt
    biologicos <- derivar_parametros_biologicos(cfg_mod)
    estado_inicial <- construir_estado_inicial(cfg_mod$s0_frac, cfg_mod$i0_personas, n_poblacion)

    r <- calibrar_beta_global(
      estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
      theta = cfg_mod$theta, sigma = cfg_mod$sigma, ve_contagio = cfg_mod$ve_contagio,
      semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
      datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = cfg_mod$rezago_semanas,
      beta_min = beta_min, beta_max = beta_max, maxit = maxit
    )

    betas_nombrados <- setNames(as.list(r$betas), paste0("beta_", seq_along(r$betas)))
    tibble(parametro = parametro, valor = valor_alt, rmse_log = r$mejor_error) %>%
      bind_cols(as_tibble(betas_nombrados))
  })
}

# ------------------------------------------------------------------------------
# Bootstrap de residuos para los beta (opcional, acotado en costo)
# ------------------------------------------------------------------------------

#' Intervalo de incertidumbre EMPIRICO para los beta, via bootstrap no
#' parametrico de los residuos semanales (en log1p) del mejor ajuste: para
#' cada repeticion, se resamplea CON REPOSICION el vector de residuos, se
#' arma una serie de muertes sinteticas = modelo + residuo resampleado, y se
#' recalibra el vector completo de beta (global) contra esa serie sintetica.
#'
#' ADVERTENCIA: no es un IC en el sentido formal de un modelo con verosimilitud
#' explicita. Asume que los residuos son intercambiables entre semanas (no lo
#' son del todo: hay autocorrelacion serial en una curva epidemica). Es una
#' cota INFERIOR de la incertidumbre real; usar con esa salvedad en el informe.
#' Cada repeticion corre una calibracion global completa, asi que es cara: n_rep
#' bajo por defecto.
#'
#' @return Matriz n_rep x n_tramos con los beta de cada repeticion bootstrap.
bootstrap_betas <- function(estado_inicial, tramos_def, biologicos, theta, sigma, ve_contagio,
                            semana_inicio_aspo, semana_inicio_vacunacion,
                            datos_semanal, n_poblacion, rezago_semanas,
                            betas_punto, beta_min, beta_max, n_rep = 20, semilla = 1) {

  set.seed(semilla)
  n_tramos <- nrow(tramos_def)

  parametros_fijos <- construir_parametros_fijos(
    theta, sigma, ve_contagio, biologicos$gamma, biologicos$delta,
    biologicos$gamma_v, biologicos$delta_v, semana_inicio_aspo, semana_inicio_vacunacion
  )
  sim_base    <- simular_serie(estado_inicial, tramos_def, betas_punto, parametros_fijos)
  comparacion <- alinear_muertes_con_rezago(sim_base$emisiones, datos_semanal, n_poblacion, rezago_semanas)
  residuos_log <- log1p(comparacion$muertes_owid) - log1p(comparacion$muertes_modelo)

  resultados <- matrix(NA_real_, nrow = n_rep, ncol = n_tramos)
  for (i in seq_len(n_rep)) {
    residuos_resampleados <- sample(residuos_log, length(residuos_log), replace = TRUE)
    muertes_sinteticas <- pmax(round(expm1(log1p(comparacion$muertes_modelo) + residuos_resampleados)), 0)

    datos_sinteticos <- datos_semanal
    # alinear por semana (comparacion viene ordenada por semana, igual que datos_semanal)
    datos_sinteticos$muertes_owid <- muertes_sinteticas[match(datos_sinteticos$semana, comparacion$semana)]

    r <- calibrar_beta_global(
      estado_inicial, tramos_def, biologicos, theta, sigma, ve_contagio,
      semana_inicio_aspo, semana_inicio_vacunacion, datos_sinteticos, n_poblacion, rezago_semanas,
      beta_min, beta_max, betas_iniciales = betas_punto, maxit = 200
    )
    resultados[i, ] <- r$betas
  }
  colnames(resultados) <- paste0("beta_", seq_len(n_tramos))
  resultados
}
