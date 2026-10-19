# ==============================================================================
# METRICS - Funciones de error y de resumen, puras y testeables
# ==============================================================================
# Todas reciben sus insumos como argumentos (datos_semanales, n_poblacion,
# rezago_semanas, etc.) en vez de leerlos de variables globales, para poder
# testearlas de forma aislada con datos sinteticos.
# ==============================================================================

library(dplyr)
library(tidyr)

#' Convierte el flujo de emisiones (FRACCION de poblacion por semana, como lo
#' devuelve R/simulate.R) a muertes ABSOLUTAS, aplica el rezago de reporte y
#' las deja alineadas por semana con las muertes observadas.
#'
#' El rezago se aplica como dplyr::lag() sobre la serie COMPLETA de
#' emisiones (ordenada por semana) ANTES de filtrar a semanas_objetivo: asi,
#' las primeras semanas de un tramo toman el rezago del historial de
#' emisiones de los tramos anteriores en vez de perderlo. Filtrar antes de
#' aplicar el lag fue un bug real en una version anterior del script (ver
#' DECISIONS.md).
#'
#' @param emisiones Tibble(semana, muertes_emitidas) en FRACCION de n_poblacion.
#' @param datos_semanales Tibble con columnas semana, muertes_owid (conteo).
#' @param n_poblacion Poblacion total, para pasar de fraccion a conteo.
#' @param rezago_semanas Semanas de rezago de reporte (ver R/config.R).
#' @param semanas_objetivo Vector de semanas a las que restringir el
#'   resultado DESPUES de aplicar el rezago (NULL = todas).
#' @return Tibble(semana, muertes_modelo, muertes_owid).
alinear_muertes_con_rezago <- function(emisiones, datos_semanales, n_poblacion, rezago_semanas, semanas_objetivo = NULL) {
  comparacion <- emisiones %>%
    arrange(semana) %>%
    mutate(muertes_modelo = dplyr::lag(muertes_emitidas * n_poblacion, n = rezago_semanas, default = 0)) %>%
    select(semana, muertes_modelo) %>%
    left_join(datos_semanales %>% select(semana, muertes_owid), by = "semana")

  if (!is.null(semanas_objetivo)) {
    comparacion <- comparacion %>% filter(semana %in% semanas_objetivo)
  }
  comparacion
}

#' Funcion de perdida MSRE (Mean Squared Relative Error), calculada sobre
#' log1p(muertes) en vez de sobre el nivel.
#'
#' POR QUE SOBRE log1p Y NO SOBRE EL NIVEL: el error relativo en nivel,
#' (modelo-obs)/obs, esta acotado por abajo en -1 (subestimar como maximo
#' predice 0) pero NO tiene techo por arriba (sobreestimar puede dar
#' cualquier valor positivo). En una dinamica de crecimiento exponencial
#' como una ola epidemica, eso vuelve al MSRE en nivel fuertemente asimetrico:
#' el optimizador prefiere un beta casi nulo (prediccion ~0, error acotado en
#' -1) antes que arriesgar una sobreestimacion explosiva (error sin techo).
#' Se verifico empiricamente: con MSRE en nivel, el tramo mas largo de la
#' serie calibraba en el piso del rango de busqueda de beta y el modelo
#' predicho colapsaba a ~0 el resto de la serie (ver DECISIONS.md). Tomando
#' el error relativo sobre log1p(muertes) se mantiene el espiritu del MSRE
#' (ponderar por igual semanas con pocas y con muchas muertes) sin esa
#' asimetria.
#'
#' Las semanas con 0 muertes observadas se excluyen (log1p(0)=0: el error
#' relativo no esta definido ahi); en la serie real de Argentina 2020-2022
#' no hay ninguna semana en 0, pero la funcion cubre el caso para datos
#' sinteticos o de otro pais/periodo.
#'
#' @return Un numero (promedio de los errores relativos al cuadrado).
calcular_msre_log <- function(emisiones, datos_semanales, n_poblacion, rezago_semanas, semanas_objetivo = NULL) {
  comparacion <- alinear_muertes_con_rezago(emisiones, datos_semanales, n_poblacion, rezago_semanas, semanas_objetivo)
  log_modelo <- log1p(comparacion$muertes_modelo)
  log_owid   <- log1p(comparacion$muertes_owid)
  error_relativo <- (log_modelo - log_owid) / log_owid
  error_relativo[comparacion$muertes_owid == 0] <- NA
  mean(error_relativo^2, na.rm = TRUE)
}

#' MSE en nivel (diagnostico, NO es la funcion de perdida que se minimiza).
calcular_mse_nivel <- function(emisiones, datos_semanales, n_poblacion, rezago_semanas) {
  comparacion <- alinear_muertes_con_rezago(emisiones, datos_semanales, n_poblacion, rezago_semanas)
  mean((comparacion$muertes_modelo - comparacion$muertes_owid)^2, na.rm = TRUE)
}

#' MSE de log1p (diagnostico, NO es la funcion de perdida que se minimiza:
#' se diferencia de calcular_msre_log() en que esta no divide por log1p(obs)).
calcular_mse_log <- function(emisiones, datos_semanales, n_poblacion, rezago_semanas) {
  comparacion <- alinear_muertes_con_rezago(emisiones, datos_semanales, n_poblacion, rezago_semanas)
  mean((log1p(comparacion$muertes_modelo) - log1p(comparacion$muertes_owid))^2, na.rm = TRUE)
}

#' Residuos semana a semana (modelo - observado, en nivel y en log1p), para
#' graficar o auditar el ajuste por tramo.
calcular_residuos <- function(emisiones, datos_semanales, n_poblacion, rezago_semanas, tramos_def = NULL) {
  comparacion <- alinear_muertes_con_rezago(emisiones, datos_semanales, n_poblacion, rezago_semanas) %>%
    mutate(
      residuo_nivel = muertes_modelo - muertes_owid,
      residuo_log   = log1p(muertes_modelo) - log1p(muertes_owid)
    )
  if (!is.null(tramos_def)) {
    comparacion <- comparacion %>%
      mutate(tramo = purrr::map_int(semana, function(s) {
        fila <- which(s >= tramos_def$semana_ini & s <= tramos_def$semana_fin)
        if (length(fila) == 0) NA_integer_ else tramos_def$tramo[fila[1]]
      }))
  }
  comparacion
}

#' Convierte la incidencia (FRACCION de poblacion por semana, de
#' R/simulate.R) a nuevas infecciones ABSOLUTAS por semana.
calcular_incidencia_absoluta <- function(incidencia, n_poblacion) {
  incidencia %>% mutate(nuevas_infecciones_abs = nuevas_infecciones * n_poblacion)
}

#' Tasa de mortalidad semanal promedio por tramo (muertes/N), a partir de las
#' emisiones CRUDAS (sin rezago: es una propiedad del modelo, no de cuando se
#' reporta la muerte).
calcular_tasa_mortalidad_por_tramo <- function(emisiones, n_poblacion, tramos_def) {
  emisiones %>%
    mutate(tramo = purrr::map_int(semana, function(s) {
      fila <- which(s >= tramos_def$semana_ini & s <= tramos_def$semana_fin)
      if (length(fila) == 0) NA_integer_ else tramos_def$tramo[fila[1]]
    })) %>%
    group_by(tramo) %>%
    summarise(tasa_mortalidad_semanal = mean(muertes_emitidas), .groups = "drop")
}

#' Tasa de contagio (incidencia) semanal promedio por tramo (nuevas
#' infecciones/N), a partir de la incidencia INTEGRADA (ver R/model.R: suma
#' de nuevas infecciones de origen susceptible y de origen vacunado), NO de
#' beta*s_prom (que subestimaba al ignorar las infecciones en vacunados y no
#' era, en rigor, una tasa de incidencia sino un termino de la EDO).
calcular_tasa_contagio_por_tramo <- function(incidencia, tramos_def) {
  incidencia %>%
    mutate(tramo = purrr::map_int(semana, function(s) {
      fila <- which(s >= tramos_def$semana_ini & s <= tramos_def$semana_fin)
      if (length(fila) == 0) NA_integer_ else tramos_def$tramo[fila[1]]
    })) %>%
    group_by(tramo) %>%
    summarise(tasa_contagio_semanal = mean(nuevas_infecciones), .groups = "drop")
}
