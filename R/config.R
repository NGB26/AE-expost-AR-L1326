# ==============================================================================
# CONFIG - Parametros del modelo SIVQRD (Argentina COVID-19)
# ==============================================================================
# Fuente unica de verdad para todos los parametros del modelo. No se lee de
# Excel: cada valor se declara aca, con su unidad, su fuente bibliografica (o
# la razon de no tenerla) y, cuando aplica, un rango plausible para la
# validacion. construir_config() arma la lista; validar_config() la audita y
# frena con un error explicito ante cualquier parametro sin valor o fuera de
# rango -- preferible a que el modelo corra en silencio con un NA o un
# numero sin sentido.
#
# Las ecuaciones del modelo (23 a 31) estan tomadas de "Modelo_SIR_para_PCR.docx"
# (seccion "4. Modelo extendido SIVQRD"). Ese documento es la UNICA fuente de
# las ecuaciones: no se verifico contra una segunda fuente independiente (ver
# DECISIONS.md, seccion "Que no se pudo verificar").
# ==============================================================================

#' Arma la configuracion completa del modelo SIVQRD.
#'
#' @return Una lista con, para cada parametro, su `valor`, `unidad` y
#'   `fuente`. Los parametros marcados con fuente "EDITAR: ..." son los que
#'   el usuario debe completar o revisar antes de usar el modelo para
#'   resultados finales.
construir_config <- function() {

  parametro <- function(valor, unidad, fuente) list(valor = valor, unidad = unidad, fuente = fuente)

  list(

    # -- Serie temporal -------------------------------------------------------
    fecha_inicio_serie = parametro(
      as.Date("2020-03-08"), "fecha",
      "Primer caso confirmado en Argentina reportado por OWID; semana 1 del modelo."
    ),
    fecha_fin_serie = parametro(
      as.Date("2022-06-05"), "fecha",
      "Ultima semana con datos homogeneos en la serie historica de OWID."
    ),

    # -- Fechas de corte de beta (5 tramos) ------------------------------------
    # APROXIMACION SEMANAL: cada fecha se ubica en la semana ENTERA que la
    # contiene (ver semana_de_fecha() en R/model.R); si la politica empezo a
    # mitad de semana, esa semana completa queda del lado del tramo nuevo.
    # Mismo criterio que usaba la columna "Semana ini" del Excel original
    # (verificado: reproduce exactamente semana 2 para ASPO, 13 para DISPO,
    # 43 para inicio de vacunacion y 91 para Omicron).
    fecha_aspo = parametro(
      as.Date("2020-03-16"), "fecha",
      "Decreto 297/20 (ASPO). Fecha exacta de decreto; el modelo la aproxima a la semana que la contiene."
    ),
    fecha_dispo = parametro(
      as.Date("2020-06-01"), "fecha",
      "Inicio de DISPO / flexibilizacion. EDITAR si se dispone de la fecha de decreto exacta por jurisdiccion."
    ),
    fecha_vacunacion = parametro(
      as.Date("2020-12-29"), "fecha",
      "Inicio de la campana de vacunacion en Argentina (OWID: primera dosis con cobertura > 0)."
    ),
    fecha_omicron = parametro(
      as.Date("2021-12-01"), "fecha",
      "Aparicion de Omicron en Argentina. EDITAR si se dispone de la fecha de deteccion de la primera secuencia."
    ),

    # -- Condiciones iniciales y demografia ------------------------------------
    s0_frac = parametro(
      1, "fraccion",
      "Se asume toda la poblacion susceptible al inicio de la serie (sin inmunidad natural previa)."
    ),
    i0_personas = parametro(
      100, "personas",
      "Semilla de infecciosos al 08-mar-2020. EDITAR: valor de trabajo, no tiene cita bibliografica especifica."
    ),

    # -- Curso clinico no vacunado ----------------------------------------------
    dias_infecciosos = parametro(
      10, "dias",
      "Duracion del periodo infeccioso (recuperacion + muerte combinadas). Consenso post-pandemia citado en el Excel original (SIRD_Argentina_PCR.xlsx, hoja Parametros)."
    ),
    ifr_base = parametro(
      0.006, "fraccion",
      "IFR (infection fatality rate) base pre-vacunacion. Meta-analisis de seroprevalencia citado en el Excel original."
    ),

    # -- Rezago de reporte ------------------------------------------------------
    # OJO: ver DECISIONS.md, seccion "Rezago de muertes vs. tiempo-a-muerte
    # implicito en la EDO", antes de cambiar este valor. No se verifico si el
    # "19 dias segun meta-analisis" original mide tiempo biologico
    # infeccion-a-muerte (que se solaparia con dias_infecciosos) o demora
    # administrativa de reporte (que es aditiva). Se mantiene el valor y el
    # comportamiento del script original (aditivo) por continuidad, pero
    # queda marcado como pendiente de confirmar.
    rezago_semanas = parametro(
      3, "semanas",
      "EDITAR/CONFIRMAR: Excel original cita '19 dias segun meta-analisis, redondeado a 3 semanas', sin nombrar el estudio. No se verifico si ese plazo ya incluye el periodo infeccioso del modelo (ver DECISIONS.md)."
    ),

    # -- Vacunacion: eficacia --------------------------------------------------
    ve_muerte = parametro(
      0.7, "fraccion",
      "Efectividad vacunal contra la muerte. Lancet Reg Health Am 2023;19:100438 (promedio ponderado del mix argentino), citado en el Excel original."
    ),
    ve_contagio = parametro(
      0.5859574, "fraccion",
      paste(
        "Promedio ponderado por el mix de vacunas aplicadas en Argentina (AstraZeneca 33%, Sinopharm 32%,",
        "Sputnik V 29%, normalizado a 94%), con odds ratios de infeccion a 2 dosis de Rearte A, et al.",
        "'Effectiveness of rAd26-rAd5, ChAdOx1 nCoV-19, and BBIBP-CorV vaccines for risk of infection with",
        "SARS-CoV-2 and death due to COVID-19 in people older than 60 years in Argentina: a test-negative,",
        "case-control, and retrospective longitudinal study'. Lancet. 2022;399(10331):1254-1264.",
        "https://doi.org/10.1016/S0140-6736(22)00011-3",
        "(OR informados en busqueda web, no verificados contra el PDF original: ver DECISIONS.md)."
      )
    ),
    dias_infecciosos_vac = parametro(
      10, "dias",
      "Duracion del periodo infeccioso en vacunados. EDITAR si la literatura da un valor distinto al no vacunado (dias_infecciosos); por ahora se usa el mismo."
    ),

    # -- Vacunacion: ritmo (sigma) ----------------------------------------------
    # Se calcula en R/data.R::calcular_sigma_desde_owid(), NO se hardcodea:
    # es el ritmo observado de la campana, no una constante biologica. Este
    # placeholder documenta la formula para que quede trazable en config.
    sigma = parametro(
      NA_real_, "1/semana",
      "Se calcula en tiempo de ejecucion a partir de la cobertura de 'vacunados completos' de OWID (ver R/data.R::calcular_sigma_desde_owid). NO es un valor de literatura: es el hazard semanal constante que reproduce la cobertura final observada."
    ),

    # -- Cuarentena (theta) -------------------------------------------------
    # theta YA NO es un parametro fijo: se calibra junto con beta (Bloque de
    # calibracion). Aca solo se fija el punto de partida y las cotas de
    # busqueda del optimizador.
    theta_inicial = parametro(0.1, "1/semana", "Punto de partida para la calibracion. No es un valor estimado."),
    theta_min      = parametro(0,   "1/semana", "Cota inferior de busqueda de theta."),
    theta_max      = parametro(2,   "1/semana", "Cota superior de busqueda de theta."),

    # -- Calibracion ------------------------------------------------------------
    beta_min = parametro(0.01, "1/semana", "Cota inferior de busqueda de beta, igual que en el Excel original."),
    beta_max = parametro(2,    "1/semana", "Cota superior de busqueda de beta, igual que en el Excel original."),
    max_iteraciones_externas = parametro(5, "iteraciones", "Tope de pasadas beta<->theta del proceso iterativo."),
    tolerancia_convergencia  = parametro(1e-6, "adimensional", "Mejora minima de error_global para seguir iterando.")
  )
}

#' Extrae solo los `valor` de la configuracion, para pasar a las funciones del
#' modelo sin tener que acceder a `cfg$parametro$valor` en cada lugar.
#'
#' @param cfg Lista devuelta por construir_config().
#' @return Lista plana nombre -> valor.
aplanar_config <- function(cfg) {
  setNames(lapply(cfg, function(p) p$valor), names(cfg))
}

#' Valida la configuracion: frena con un error explicito si falta un valor o
#' si un valor esta fuera de un rango plausible. Se llama SIEMPRE antes de
#' correr el modelo (scripts/run.R lo hace de entrada).
#'
#' @param cfg Lista devuelta por construir_config().
#' @param exigir_sigma Si es TRUE (default), sigma debe tener un valor no-NA;
#'   poner en FALSE solo si se va a calcular sigma DESPUES de validar (caso
#'   de uso normal: se valida de nuevo despues de calcularlo).
validar_config <- function(cfg, exigir_sigma = TRUE) {

  valores <- aplanar_config(cfg)

  rangos <- list(
    s0_frac                  = c(0, 1),
    ifr_base                 = c(0, 1),
    ve_muerte                = c(0, 1),
    ve_contagio               = c(0, 1),
    dias_infecciosos          = c(1, 60),
    dias_infecciosos_vac       = c(1, 60),
    rezago_semanas            = c(0, 10),
    i0_personas                = c(0, Inf),
    theta_min                  = c(0, Inf),
    theta_max                  = c(0, Inf),
    beta_min                   = c(0, Inf),
    beta_max                   = c(0, Inf),
    sigma                      = c(0, 10)
  )

  parametros_a_exigir <- setdiff(names(cfg), if (exigir_sigma) character(0) else "sigma")

  for (nombre in parametros_a_exigir) {
    valor <- valores[[nombre]]
    if (length(valor) == 0 || all(is.na(valor))) {
      stop(sprintf(
        "config invalida: el parametro '%s' no tiene valor (NA). Fuente registrada: '%s'. Completarlo en R/config.R antes de correr el modelo.",
        nombre, cfg[[nombre]]$fuente
      ), call. = FALSE)
    }
  }

  for (nombre in intersect(names(rangos), parametros_a_exigir)) {
    valor <- valores[[nombre]]
    rango <- rangos[[nombre]]
    if (!is.na(valor) && (valor < rango[1] || valor > rango[2])) {
      stop(sprintf(
        "config invalida: '%s' = %s esta fuera del rango plausible [%s, %s].",
        nombre, valor, rango[1], rango[2]
      ), call. = FALSE)
    }
  }

  if (valores$theta_min >= valores$theta_max) {
    stop("config invalida: theta_min debe ser menor que theta_max.", call. = FALSE)
  }
  if (valores$beta_min >= valores$beta_max) {
    stop("config invalida: beta_min debe ser menor que beta_max.", call. = FALSE)
  }
  if (valores$fecha_inicio_serie >= valores$fecha_fin_serie) {
    stop("config invalida: fecha_inicio_serie debe ser anterior a fecha_fin_serie.", call. = FALSE)
  }
  fechas_corte <- c(valores$fecha_aspo, valores$fecha_dispo, valores$fecha_vacunacion, valores$fecha_omicron)
  if (is.unsorted(fechas_corte, strictly = TRUE)) {
    stop("config invalida: fecha_aspo < fecha_dispo < fecha_vacunacion < fecha_omicron no se cumple.", call. = FALSE)
  }

  invisible(TRUE)
}
