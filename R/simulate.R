# ==============================================================================
# SIMULATE - Integracion numerica del SIVQRD, sin estado global
# ==============================================================================
# Todas las funciones son puras: reciben sus insumos como argumentos y
# devuelven su resultado, sin leer ni modificar variables globales. Esto
# reemplaza a ajustar_betas() y ajustar_betas_con_beta_fijo() del script
# original (que repetian casi linea por linea el mismo bucle sobre tramos,
# una optimizando beta y la otra sin optimizar): aca hay una sola funcion de
# simulacion hacia adelante (simular_serie) y una sola funcion de ajuste
# secuencial (ajustar_betas_secuencial), y la calibracion (R/calibrate.R)
# combina ambas en vez de duplicar el bucle.
# ==============================================================================

library(deSolve)
library(dplyr)
library(tibble)

#' Simula un tramo con beta CONSTANTE, usando deSolve::ode() (metodo lsoda,
#' adaptive step: mas preciso que un Euler de paso semanal como el de la
#' planilla original).
#'
#' @param estado_inicial Vector c(s, i, iv, v, q, r, d, inf_acum).
#' @param semana_ini,semana_fin Semanas (enteras) que cubre el tramo.
#' @param beta_valor Beta constante para todo el tramo.
#' @param parametros_fijos Lista para derivadas_sivqrd (sin beta_fun).
#' @return Tibble con columnas semana, s, i, iv, v, q, r, d, inf_acum. Tiene
#'   (semana_fin - semana_ini + 2) filas: incluye un punto EXTRA en
#'   semana_fin+1 (estado al inicio de la semana siguiente), necesario para
#'   poder calcular con diff() el flujo semanal de la ULTIMA semana del
#'   tramo. simular_serie() descarta ese punto extra al concatenar tramos
#'   (para no duplicarlo con el primero del tramo siguiente).
simular_tramo <- function(estado_inicial, semana_ini, semana_fin, beta_valor, parametros_fijos) {
  tiempos <- seq(semana_ini, semana_fin + 1, by = 1)
  parms <- c(parametros_fijos, list(beta_fun = function(t) beta_valor))
  # maxsteps subido del default de lsoda (5000): en combinaciones EXTREMAS de
  # beta y theta (cerca de sus cotas superiores a la vez, como las que prueban
  # los barridos de sensibilidad sobre la serie larga 2020-2023) el sistema se
  # vuelve rigido. Se intenta primero lsoda (rapido, auto-switch); si lsoda
  # falla o devuelve no-finitos, se reintenta con un solver stiff explicito
  # (BDF via lsode, "method=bdf"). Solo si AMBOS fallan se levanta un error
  # explicito, que los optimizadores de R/calibrate.R atrapan y convierten en
  # una penalizacion finita (para no cortar la calibracion por una region de
  # parametros infactible). Las ecuaciones con beta/theta razonables integran
  # bien con lsoda; esto es una red de seguridad para el barrido del optimo.
  integrar <- function(metodo) {
    ode(y = estado_inicial, times = tiempos, func = derivadas_sivqrd, parms = parms,
        method = metodo, maxsteps = 100000)
  }
  salida <- tryCatch(integrar("lsoda"), error = function(e) NULL)
  if (is.null(salida) || any(!is.finite(salida[, -1]))) {
    salida <- tryCatch(integrar("bdf"), error = function(e) NULL)  # reintento stiff
  }
  if (is.null(salida) || any(!is.finite(salida[, -1]))) {
    stop(sprintf(
      "simular_tramo(): la integracion no convergio entre semana %s y %s (beta=%.4f) ni con lsoda ni con bdf: el sistema quedo numericamente rigido/inestable en ese rango de parametros.",
      semana_ini, semana_fin, beta_valor
    ), call. = FALSE)
  }
  salida %>%
    as.data.frame() %>%
    as_tibble() %>%
    rename(semana = time)
}

#' Simula la serie COMPLETA (todos los tramos, uno atras del otro) con un
#' vector de beta YA FIJO (sin optimizar). Es la pieza comun que tanto la
#' calibracion secuencial como la evaluacion de theta (con beta fijo) usan.
#'
#' @param estado_inicial Estado al comienzo del primer tramo.
#' @param tramos_def Tibble con columnas tramo, semana_ini, semana_fin.
#' @param betas Vector de beta, uno por fila de tramos_def (mismo orden).
#' @param parametros_fijos Lista para derivadas_sivqrd (sin beta_fun; ya
#'   debe traer el theta que corresponda).
#' @return Lista con:
#'   - emisiones: tibble(semana, muertes_emitidas) -- flujo CRUDO de muertes
#'     por semana (antes del rezago de reporte).
#'   - incidencia: tibble(semana, nuevas_infecciones) -- flujo de NUEVAS
#'     infecciones por semana (de origen susceptible + vacunado).
#'   - simulacion: tibble con la trayectoria semanal de los 7 compartimentos
#'     fisicos (s..d) por tramo, SIN el punto extra de cada simular_tramo().
#'   - estado_final: estado al cierre del ultimo tramo.
simular_serie <- function(estado_inicial, tramos_def, betas, parametros_fijos) {
  stopifnot(length(betas) == nrow(tramos_def))

  estado_actual  <- estado_inicial
  emisiones      <- tibble(semana = numeric(), muertes_emitidas = numeric())
  incidencia     <- tibble(semana = numeric(), nuevas_infecciones = numeric())
  resultados_sim <- vector("list", nrow(tramos_def))

  for (k in seq_len(nrow(tramos_def))) {
    semana_ini_k  <- tramos_def$semana_ini[k]
    semana_fin_k  <- tramos_def$semana_fin[k]
    semanas_tramo <- seq(semana_ini_k, semana_fin_k, by = 1)

    sim <- simular_tramo(estado_actual, semana_ini_k, semana_fin_k, betas[k], parametros_fijos)

    # Flujo semanal = incremento del acumulado entre el inicio de la semana
    # w y el inicio de w+1 (ver nota en simular_tramo sobre el punto extra).
    # Nota de unidades: estas son FRACCIONES de poblacion por semana; se
    # multiplican por n_poblacion en R/calibrate.R y R/metrics.R, no aca.
    emisiones  <- bind_rows(emisiones,  tibble(semana = semanas_tramo, muertes_emitidas = diff(sim$d))) %>% arrange(semana)
    incidencia <- bind_rows(incidencia, tibble(semana = semanas_tramo, nuevas_infecciones = diff(sim$inf_acum))) %>% arrange(semana)

    resultados_sim[[k]] <- sim %>%
      filter(semana <= semana_fin_k) %>%
      mutate(tramo = k, nombre_tramo = tramos_def$nombre[k])

    estado_actual <- sim %>%
      slice_tail(n = 1) %>%
      select(s, i, iv, v, q, r, d, inf_acum) %>%
      unlist()
  }

  list(
    emisiones  = emisiones,
    incidencia = incidencia,
    simulacion = bind_rows(resultados_sim),
    estado_final = estado_actual
  )
}

#' Ajusta un beta POR TRAMO, de forma secuencial: para cada tramo (en orden),
#' busca con optimize() el beta que minimiza funcion_error() DENTRO de ese
#' tramo, usando el estado final del tramo anterior como condicion inicial
#' del siguiente. Es el metodo baseline del proyecto (confirmado con el
#' usuario): cada beta resuelve un problema LOCAL (su propio tramo), no el
#' error global de toda la serie -- ver R/calibrate.R para la discusion de
#' identificabilidad y la alternativa de ajuste conjunto.
#'
#' @param estado_inicial Estado al comienzo del primer tramo.
#' @param tramos_def Tibble con columnas tramo, semana_ini, semana_fin.
#' @param parametros_fijos Lista para derivadas_sivqrd (sin beta_fun; ya con
#'   el theta vigente).
#' @param funcion_error Function(emisiones_historial, semanas_objetivo) ->
#'   numero. Recibe el historial COMPLETO de emisiones (tramos ya resueltos +
#'   el candidato actual) para que el rezago de reporte pueda mirar hacia
#'   atras del tramo (ver R/metrics.R).
#' @param beta_min,beta_max Cotas de busqueda de beta.
#' @param tol Tolerancia de optimize().
#' @return Lista con betas (vector), emisiones, incidencia, simulacion,
#'   estado_final -- mismo formato que simular_serie().
ajustar_betas_secuencial <- function(estado_inicial, tramos_def, parametros_fijos,
                                      funcion_error, beta_min, beta_max, tol = 1e-5) {

  estado_actual   <- estado_inicial
  emisiones       <- tibble(semana = numeric(), muertes_emitidas = numeric())
  incidencia      <- tibble(semana = numeric(), nuevas_infecciones = numeric())
  resultados_sim  <- vector("list", nrow(tramos_def))
  betas_ajustados <- numeric(nrow(tramos_def))

  for (k in seq_len(nrow(tramos_def))) {

    semana_ini_k  <- tramos_def$semana_ini[k]
    semana_fin_k  <- tramos_def$semana_fin[k]
    semanas_tramo <- seq(semana_ini_k, semana_fin_k, by = 1)

    calcular_error_tramo <- function(beta_k) {
      # tryCatch -> penalizacion finita si el ODE no integra para este beta
      # (ver nota en R/calibrate.R y en simular_tramo), para no cortar el
      # ajuste secuencial en una region de parametros infactible.
      err <- tryCatch({
        sim <- simular_tramo(estado_actual, semana_ini_k, semana_fin_k, beta_k, parametros_fijos)
        candidatas <- tibble(semana = semanas_tramo, muertes_emitidas = diff(sim$d))
        funcion_error(bind_rows(emisiones, candidatas), semanas_objetivo = semanas_tramo)
      }, error = function(e) NA_real_)
      if (!is.finite(err)) 1e6 else err
    }

    ajuste <- optimize(calcular_error_tramo, interval = c(beta_min, beta_max), tol = tol)
    betas_ajustados[k] <- ajuste$minimum

    sim_final <- simular_tramo(estado_actual, semana_ini_k, semana_fin_k, betas_ajustados[k], parametros_fijos)

    emisiones  <- bind_rows(emisiones,  tibble(semana = semanas_tramo, muertes_emitidas = diff(sim_final$d))) %>% arrange(semana)
    incidencia <- bind_rows(incidencia, tibble(semana = semanas_tramo, nuevas_infecciones = diff(sim_final$inf_acum))) %>% arrange(semana)

    resultados_sim[[k]] <- sim_final %>%
      filter(semana <= semana_fin_k) %>%
      mutate(tramo = k, nombre_tramo = tramos_def$nombre[k])

    estado_actual <- sim_final %>%
      slice_tail(n = 1) %>%
      select(s, i, iv, v, q, r, d, inf_acum) %>%
      unlist()
  }

  list(
    betas = betas_ajustados, emisiones = emisiones, incidencia = incidencia,
    simulacion = bind_rows(resultados_sim), estado_final = estado_actual
  )
}
