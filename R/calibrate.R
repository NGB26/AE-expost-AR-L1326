# ==============================================================================
# CALIBRATE - Ajuste de beta y theta: baseline, alternativa conjunta,
# identificabilidad, sensibilidad y (opcional) bootstrap.
# ==============================================================================
# Todas las funciones reciben sus insumos explicitos (estado_inicial,
# tramos_def, parametros biologicos, datos_semanal, n_poblacion, cotas) y no
# leen variables globales.
# ==============================================================================

library(dplyr)
library(purrr)
library(tidyr)

#' Arma la lista de parametros fijos para un theta dado (el resto -- sigma,
#' ve_contagio, gamma, delta, gamma_v, delta_v -- no cambia durante la
#' calibracion).
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
# Baseline: secuencial por tramo, alternando con theta
# ------------------------------------------------------------------------------

#' Metodo BASELINE del proyecto (confirmado con el usuario): para un theta
#' fijo, ajusta los 5 beta de forma secuencial por tramo (cada uno minimiza
#' el error DENTRO de su propio tramo); despues, con esos beta ya fijos,
#' ajusta theta (un unico valor para minimizar el error GLOBAL); se repite
#' hasta que el error global deje de mejorar.
#'
#' Beta minimiza un objetivo LOCAL (por tramo) y theta uno GLOBAL: son dos
#' objetivos distintos, asi que alternar entre ambos no garantiza que cada
#' pasada mejore el ajuste global -- se verifico empiricamente (ver
#' DECISIONS.md). Por eso se guarda el mejor par (beta, theta) visto segun el
#' error global, y el proceso corta apenas una pasada deja de mejorarlo.
#'
#' @return Lista con betas, theta, emisiones, incidencia, simulacion,
#'   historial_ajuste (tibble iteracion/theta/error_global), mejor_error,
#'   n_evaluaciones_ode (cuantas veces se corrio simular_tramo en total, como
#'   medida de costo computacional).
calibrar_secuencial_alternante <- function(estado_inicial, tramos_def, biologicos, sigma, ve_contagio,
                                            semana_inicio_aspo, semana_inicio_vacunacion,
                                            datos_semanal, n_poblacion, rezago_semanas,
                                            theta_inicial, theta_min, theta_max,
                                            beta_min, beta_max,
                                            max_iteraciones = 5, tolerancia = 1e-6) {

  n_evaluaciones_ode <- 0

  funcion_error <- function(emisiones, semanas_objetivo) {
    calcular_msre_log(emisiones, datos_semanal, n_poblacion, rezago_semanas, semanas_objetivo)
  }

  ajustar_beta_dado_theta <- function(theta_valor) {
    parametros_fijos <- construir_parametros_fijos(
      theta_valor, sigma, ve_contagio, biologicos$gamma, biologicos$delta,
      biologicos$gamma_v, biologicos$delta_v, semana_inicio_aspo, semana_inicio_vacunacion
    )
    resultado <- ajustar_betas_secuencial(estado_inicial, tramos_def, parametros_fijos, funcion_error, beta_min, beta_max)
    n_evaluaciones_ode <<- n_evaluaciones_ode + nrow(tramos_def) * 20  # ~20 evals/optimize(), orden de magnitud
    resultado
  }

  evaluar_error_global_theta <- function(theta_candidato, betas_fijos) {
    parametros_fijos <- construir_parametros_fijos(
      theta_candidato, sigma, ve_contagio, biologicos$gamma, biologicos$delta,
      biologicos$gamma_v, biologicos$delta_v, semana_inicio_aspo, semana_inicio_vacunacion
    )
    sim <- simular_serie(estado_inicial, tramos_def, betas_fijos, parametros_fijos)
    n_evaluaciones_ode <<- n_evaluaciones_ode + nrow(tramos_def)
    calcular_msre_log(sim$emisiones, datos_semanal, n_poblacion, rezago_semanas)
  }

  theta_actual     <- theta_inicial
  mejor_theta      <- theta_inicial
  mejor_betas      <- NULL
  mejor_error      <- Inf
  historial_ajuste <- tibble(iteracion = integer(), theta = double(), error_global = double())

  for (iteracion in seq_len(max_iteraciones)) {

    ajuste_beta    <- ajustar_beta_dado_theta(theta_actual)
    betas_actuales <- ajuste_beta$betas

    ajuste_theta <- optimize(
      evaluar_error_global_theta, interval = c(theta_min, theta_max), tol = 1e-4,
      betas_fijos = betas_actuales
    )

    historial_ajuste <- bind_rows(historial_ajuste, tibble(
      iteracion = iteracion, theta = ajuste_theta$minimum, error_global = ajuste_theta$objective
    ))

    if (ajuste_theta$objective < mejor_error - tolerancia) {
      mejor_error  <- ajuste_theta$objective
      mejor_theta  <- ajuste_theta$minimum
      mejor_betas  <- betas_actuales
      theta_actual <- ajuste_theta$minimum
    } else {
      break
    }
  }

  parametros_finales <- construir_parametros_fijos(
    mejor_theta, sigma, ve_contagio, biologicos$gamma, biologicos$delta,
    biologicos$gamma_v, biologicos$delta_v, semana_inicio_aspo, semana_inicio_vacunacion
  )
  resultado_final <- simular_serie(estado_inicial, tramos_def, mejor_betas, parametros_finales)

  list(
    betas = mejor_betas, theta = mejor_theta,
    emisiones = resultado_final$emisiones, incidencia = resultado_final$incidencia,
    simulacion = resultado_final$simulacion,
    historial_ajuste = historial_ajuste, mejor_error = mejor_error,
    n_evaluaciones_ode = n_evaluaciones_ode, metodo = "secuencial_alternante"
  )
}

# ------------------------------------------------------------------------------
# Alternativa: ajuste CONJUNTO global de los 5 beta + theta
# ------------------------------------------------------------------------------

#' Ajusta los 5 beta Y theta de forma CONJUNTA (optim(), L-BFGS-B, con cotas
#' explicitas), minimizando el MSRE-log GLOBAL de toda la serie de una sola
#' vez -- a diferencia del baseline, que optimiza cada beta contra el error
#' de SU tramo. Sirve de punto de comparacion para la pregunta de
#' identificabilidad: si el optimo conjunto da una combinacion (beta,theta)
#' MUY distinta a la del baseline con un error global similar o mejor, es
#' evidencia de que la superficie de perdida tiene una cresta ancha (varias
#' combinaciones casi-equivalentes) y que el baseline secuencial es solo UNA
#' solucion razonable entre varias, no "la" solucion.
#'
#' @param valores_iniciales Vector c(beta_1..beta_5, theta) para arrancar
#'   optim(); si es NULL, se usa el resultado del baseline secuencial como
#'   punto de partida (practica estandar: partir de una solucion razonable).
#' @return Lista con la misma forma que calibrar_secuencial_alternante(), mas
#'   `convergencia` (codigo de optim(), 0 = convergencia reportada por el
#'   optimizador) y `n_evaluaciones` (contador real de optim(), no estimado).
calibrar_conjunto_global <- function(estado_inicial, tramos_def, biologicos, sigma, ve_contagio,
                                      semana_inicio_aspo, semana_inicio_vacunacion,
                                      datos_semanal, n_poblacion, rezago_semanas,
                                      beta_min, beta_max, theta_min, theta_max,
                                      valores_iniciales = NULL) {

  n_tramos <- nrow(tramos_def)

  evaluar <- function(theta_params) {
    betas <- theta_params[1:n_tramos]
    theta_valor <- theta_params[n_tramos + 1]
    parametros_fijos <- construir_parametros_fijos(
      theta_valor, sigma, ve_contagio, biologicos$gamma, biologicos$delta,
      biologicos$gamma_v, biologicos$delta_v, semana_inicio_aspo, semana_inicio_vacunacion
    )
    sim <- simular_serie(estado_inicial, tramos_def, betas, parametros_fijos)
    calcular_msre_log(sim$emisiones, datos_semanal, n_poblacion, rezago_semanas)
  }

  if (is.null(valores_iniciales)) {
    valores_iniciales <- c(rep(mean(c(beta_min, beta_max)) / 2, n_tramos), mean(c(theta_min, theta_max)))
  }

  ajuste <- optim(
    par = valores_iniciales, fn = evaluar, method = "L-BFGS-B",
    lower = c(rep(beta_min, n_tramos), theta_min),
    upper = c(rep(beta_max, n_tramos), theta_max),
    control = list(factr = 1e7, maxit = 200)
  )

  betas_finales <- ajuste$par[1:n_tramos]
  theta_final   <- ajuste$par[n_tramos + 1]
  parametros_finales <- construir_parametros_fijos(
    theta_final, sigma, ve_contagio, biologicos$gamma, biologicos$delta,
    biologicos$gamma_v, biologicos$delta_v, semana_inicio_aspo, semana_inicio_vacunacion
  )
  resultado_final <- simular_serie(estado_inicial, tramos_def, betas_finales, parametros_finales)

  list(
    betas = betas_finales, theta = theta_final,
    emisiones = resultado_final$emisiones, incidencia = resultado_final$incidencia,
    simulacion = resultado_final$simulacion,
    mejor_error = ajuste$value, convergencia = ajuste$convergence,
    n_evaluaciones = ajuste$counts[["function"]], metodo = "conjunto_global"
  )
}

# ------------------------------------------------------------------------------
# Identificabilidad: perfil de perdida en la grilla (beta, theta)
# ------------------------------------------------------------------------------

#' Genera una grilla de error GLOBAL en funcion de (beta_tramo, theta), para
#' un tramo elegido, con los demas beta fijos en un valor de referencia
#' (tipicamente el resultado del baseline). Una "cresta" diagonal en el mapa
#' de calor (muchas combinaciones beta-theta con error similar) es la firma
#' visual de no-identificabilidad entre ambos parametros: theta reduce la
#' fuerza de infeccion efectiva sacando gente de I hacia Q, de forma
#' parecida a bajar beta directamente, asi que distintas combinaciones
#' pueden explicar casi la misma curva de muertes.
#'
#' @param tramo_objetivo Indice de fila en tramos_def cuyo beta se barre.
#' @param betas_referencia Vector de beta (uno por tramo) a usar para los
#'   tramos que NO son tramo_objetivo (p.ej., el resultado del baseline).
#' @param rango_beta,rango_theta Vectores numericos con los valores a probar.
#' @return Tibble(beta, theta, error_global).
generar_perfil_identificabilidad <- function(estado_inicial, tramos_def, biologicos, sigma, ve_contagio,
                                              semana_inicio_aspo, semana_inicio_vacunacion,
                                              datos_semanal, n_poblacion, rezago_semanas,
                                              tramo_objetivo, betas_referencia, rango_beta, rango_theta) {

  grilla <- expand_grid(beta = rango_beta, theta = rango_theta)

  # Combinaciones extremas (beta y theta a la vez cerca de sus cotas
  # superiores) pueden volver el sistema numericamente rigido para lsoda
  # (ver maxsteps en R/simulate.R::simular_tramo). Esas celdas se reportan
  # como NA_real_, NO como un numero de error grande: tratarlas como "error
  # alto" las confundiria con evidencia real de mal ajuste, cuando en
  # realidad es una falla de convergencia numerica.
  n_fallos <- 0
  grilla$error_global <- purrr::map2_dbl(grilla$beta, grilla$theta, function(beta_val, theta_val) {
    betas <- betas_referencia
    betas[tramo_objetivo] <- beta_val
    parametros_fijos <- construir_parametros_fijos(
      theta_val, sigma, ve_contagio, biologicos$gamma, biologicos$delta,
      biologicos$gamma_v, biologicos$delta_v, semana_inicio_aspo, semana_inicio_vacunacion
    )
    resultado <- tryCatch({
      sim <- simular_serie(estado_inicial, tramos_def, betas, parametros_fijos)
      calcular_msre_log(sim$emisiones, datos_semanal, n_poblacion, rezago_semanas)
    }, error = function(e) NA_real_)
    if (is.na(resultado)) n_fallos <<- n_fallos + 1
    resultado
  })

  if (n_fallos > 0) {
    warning(sprintf(
      "generar_perfil_identificabilidad(): %d de %d celdas de la grilla no convergieron numericamente y quedaron en NA (tipicamente beta y theta simultaneamente cerca de sus cotas superiores).",
      n_fallos, nrow(grilla)
    ), call. = FALSE)
  }

  grilla
}

# ------------------------------------------------------------------------------
# Multi-start: robustez del optimo del metodo secuencial-alternante
# ------------------------------------------------------------------------------

#' Repite calibrar_secuencial_alternante() desde varios theta_inicial
#' distintos. Si todos convergen al mismo (o casi el mismo) mejor_error y
#' theta final, el optimo encontrado es robusto al punto de partida; si no,
#' hay multiples optimos locales y el resultado del baseline depende de
#' donde arranco theta_inicial.
#'
#' @param thetas_iniciales Vector de puntos de partida a probar.
#' @return Tibble(theta_inicial, theta_final, mejor_error, betas [list-col]).
calibrar_multistart <- function(thetas_iniciales, ...) {
  resultados <- purrr::map(thetas_iniciales, function(theta_ini) {
    r <- calibrar_secuencial_alternante(theta_inicial = theta_ini, ...)
    tibble(theta_inicial = theta_ini, theta_final = r$theta, mejor_error = r$mejor_error, betas = list(r$betas))
  })
  bind_rows(resultados)
}

# ------------------------------------------------------------------------------
# Sensibilidad OAT (one-at-a-time)
# ------------------------------------------------------------------------------

#' Perturba UN parametro biologico a la vez (dias_infecciosos, ifr_base,
#' rezago_semanas, ve_muerte, ve_contagio o sigma) a un valor alternativo,
#' re-deriva lo que dependa de el (gamma/delta via derivar_parametros_biologicos)
#' y vuelve a correr la calibracion secuencial-alternante completa, para ver
#' cuanto cambian beta, theta, R0 y el error de ajuste. No es un analisis de
#' sensibilidad formal (Sobol, FAST): es un barrido OAT simple, mas rapido y
#' mas facil de interpretar para un informe tecnico. Documentar cualquier
#' conclusion causal con esta cautela.
#'
#' @param cfg_valores Lista plana de aplanar_config(cfg).
#' @param parametro Nombre del parametro a perturbar (debe existir en cfg_valores).
#' @param valores_alternativos Vector de valores alternativos a probar
#'   (ademas del valor base, que siempre se incluye).
#' @return Tibble(parametro, valor, mejor_error, theta, r0_tramo1, ...).
analisis_sensibilidad <- function(cfg_valores, parametro, valores_alternativos,
                                   tramos_def, datos_semanal, n_poblacion,
                                   semana_inicio_aspo, semana_inicio_vacunacion, ...) {

  valores_a_probar <- unique(c(cfg_valores[[parametro]], valores_alternativos))

  resultados <- purrr::map(valores_a_probar, function(valor_alt) {
    cfg_mod <- cfg_valores
    cfg_mod[[parametro]] <- valor_alt
    biologicos <- derivar_parametros_biologicos(cfg_mod)
    estado_inicial <- construir_estado_inicial(cfg_mod$s0_frac, cfg_mod$i0_personas, n_poblacion)

    r <- calibrar_secuencial_alternante(
      estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
      sigma = cfg_mod$sigma, ve_contagio = cfg_mod$ve_contagio,
      semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
      datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = cfg_mod$rezago_semanas,
      theta_inicial = cfg_mod$theta_inicial, theta_min = cfg_mod$theta_min, theta_max = cfg_mod$theta_max,
      beta_min = cfg_mod$beta_min, beta_max = cfg_mod$beta_max, ...
    )

    tibble(
      parametro = parametro, valor = valor_alt, mejor_error = r$mejor_error, theta = r$theta,
      beta_1 = r$betas[1], beta_2 = r$betas[2], beta_3 = r$betas[3], beta_4 = r$betas[4], beta_5 = r$betas[5]
    )
  })

  bind_rows(resultados)
}

# ------------------------------------------------------------------------------
# Bootstrap de residuos (opcional, acotado en costo)
# ------------------------------------------------------------------------------

#' Intervalo de incertidumbre EMPIRICO para theta, via bootstrap no
#' parametrico de los residuos semanales (en log1p) del mejor ajuste: para
#' cada repeticion, se resamplea CON REPOSICION el vector de residuos, se
#' construye una serie de muertes sinteticas = modelo + residuo resampleado,
#' y se recalibra theta (con los beta del ajuste original fijos, por costo
#' computacional) contra esa serie sintetica.
#'
#' ADVERTENCIA SOBRE QUE TAN DEFENDIBLE ES ESTO: no es un intervalo de
#' confianza en el sentido formal de un modelo estadistico con una
#' verosimilitud explicita -- es una aproximacion empirica que asume que los
#' residuos son intercambiables entre semanas (no lo son del todo: hay
#' autocorrelacion serial evidente en una curva epidemica) y que fijar los
#' beta originales no subestima la incertidumbre total. Se ofrece como una
#' cota INFERIOR de la incertidumbre real, no como el intervalo definitivo.
#' Usar con esa salvedad explicita en cualquier informe.
#'
#' @param n_rep Numero de repeticiones bootstrap. Cada una corre una
#'   optimizacion de theta (rapida, 1 dimension): n_rep=30 por defecto para
#'   mantener el costo bajo; subir si se dispone de mas tiempo de computo.
#' @return Vector numerico de longitud n_rep con los theta bootstrap.
bootstrap_theta <- function(estado_inicial, tramos_def, biologicos, sigma, ve_contagio,
                             semana_inicio_aspo, semana_inicio_vacunacion,
                             datos_semanal, n_poblacion, rezago_semanas,
                             betas_fijos, theta_min, theta_max, n_rep = 30, semilla = 1) {

  set.seed(semilla)

  parametros_base <- construir_parametros_fijos(
    mean(c(theta_min, theta_max)), sigma, ve_contagio, biologicos$gamma, biologicos$delta,
    biologicos$gamma_v, biologicos$delta_v, semana_inicio_aspo, semana_inicio_vacunacion
  )
  sim_base <- simular_serie(estado_inicial, tramos_def, betas_fijos, parametros_base)
  comparacion <- alinear_muertes_con_rezago(sim_base$emisiones, datos_semanal, n_poblacion, rezago_semanas)
  residuos_log <- log1p(comparacion$muertes_owid) - log1p(comparacion$muertes_modelo)

  evaluar_theta_sobre_serie_sintetica <- function(theta_val, muertes_sinteticas) {
    parametros_fijos <- construir_parametros_fijos(
      theta_val, sigma, ve_contagio, biologicos$gamma, biologicos$delta,
      biologicos$gamma_v, biologicos$delta_v, semana_inicio_aspo, semana_inicio_vacunacion
    )
    sim <- simular_serie(estado_inicial, tramos_def, betas_fijos, parametros_fijos)
    pred <- alinear_muertes_con_rezago(sim$emisiones, datos_semanal, n_poblacion, rezago_semanas)
    error_rel <- (log1p(pred$muertes_modelo) - log1p(muertes_sinteticas)) / log1p(muertes_sinteticas)
    error_rel[muertes_sinteticas == 0] <- NA
    mean(error_rel^2, na.rm = TRUE)
  }

  vapply(seq_len(n_rep), function(i) {
    residuos_resampleados <- sample(residuos_log, length(residuos_log), replace = TRUE)
    muertes_sinteticas <- pmax(expm1(log1p(comparacion$muertes_owid) - residuos_log + residuos_resampleados), 0)
    ajuste <- optimize(
      evaluar_theta_sobre_serie_sintetica, interval = c(theta_min, theta_max), tol = 1e-3,
      muertes_sinteticas = muertes_sinteticas
    )
    ajuste$minimum
  }, numeric(1))
}
