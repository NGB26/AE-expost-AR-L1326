# ==============================================================================
# DATA - Ingesta de OWID: validacion, correcciones, cache y agregacion semanal
# ==============================================================================

library(dplyr)
library(readr)
library(tibble)
library(jsonlite)
library(digest)

COLUMNAS_OWID_REQUERIDAS <- c(
  "iso_code", "date", "new_cases", "new_deaths", "new_tests",
  "positive_rate", "people_fully_vaccinated", "population"
)

# ------------------------------------------------------------------------------
# Descarga y cache
# ------------------------------------------------------------------------------

#' Descarga el CSV historico completo de OWID y lo filtra a Argentina y las
#' columnas que el modelo necesita (el archivo completo pesa ~100MB; el
#' filtrado se guarda en cache, no el original).
#'
#' @param url URL del CSV. Default: el historico "congelado" de GitHub
#'   (owid/covid-19-data dejo de actualizar sus extracts en 2024; Argentina
#'   dejo de reportar diariamente en 2022, asi que para el periodo de este
#'   modelo es, de hecho, un snapshot estable).
#' @return Tibble diario de Argentina con las columnas de COLUMNAS_OWID_REQUERIDAS.
descargar_owid_argentina <- function(url = "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv") {
  datos <- read_csv(url, show_col_types = FALSE)
  faltantes <- setdiff(COLUMNAS_OWID_REQUERIDAS, names(datos))
  if (length(faltantes) > 0) {
    stop(sprintf(
      "El CSV de OWID descargado de '%s' no tiene las columnas esperadas: %s. El esquema del archivo pudo haber cambiado.",
      url, paste(faltantes, collapse = ", ")
    ), call. = FALSE)
  }
  datos %>%
    filter(iso_code == "ARG") %>%
    select(all_of(COLUMNAS_OWID_REQUERIDAS))
}

#' Hash reproducible (sha256) de un ARCHIVO en disco (no del objeto R en
#' memoria: digest::digest() sobre un data.frame hashea su representacion
#' interna -- tipos exactos, atributos -- que NO sobrevive un roundtrip
#' write_csv()/read_csv() aunque los VALORES queden identicos. Hasheando el
#' archivo de una vez evita falsos negativos de "la cache no coincide".
hash_archivo <- function(ruta) {
  digest::digest(file = ruta, algo = "sha256")
}

#' Guarda datos + metadata (fecha de descarga, URL fuente, hash, n filas) en
#' dir_cache, para poder reproducir la corrida sin volver a pegarle a la red
#' y para poder auditar DE DONDE salio el archivo que se uso.
guardar_cache_owid <- function(datos, url_fuente, dir_cache = "data/cache") {
  dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
  ruta_csv  <- file.path(dir_cache, "owid_argentina_diario.csv")
  ruta_meta <- file.path(dir_cache, "owid_argentina_diario.meta.json")

  write_csv(datos, ruta_csv)
  metadata <- list(
    url_fuente = url_fuente,
    descargado_en = format(Sys.time(), tz = "UTC", usetz = TRUE),
    n_filas = nrow(datos),
    hash_sha256 = hash_archivo(ruta_csv),
    columnas = names(datos)
  )
  write_json(metadata, ruta_meta, auto_unbox = TRUE, pretty = TRUE)
  invisible(metadata)
}

#' Lee la cache si existe y es consistente (el hash guardado coincide con el
#' hash del CSV tal como quedo en disco). Devuelve NULL si no hay cache o si
#' quedo corrupta/editada a mano sin actualizar el meta.json.
leer_cache_owid <- function(dir_cache = "data/cache") {
  ruta_csv  <- file.path(dir_cache, "owid_argentina_diario.csv")
  ruta_meta <- file.path(dir_cache, "owid_argentina_diario.meta.json")
  if (!file.exists(ruta_csv) || !file.exists(ruta_meta)) return(NULL)

  metadata <- fromJSON(ruta_meta)

  if (hash_archivo(ruta_csv) != metadata$hash_sha256) {
    warning(sprintf(
      "La cache en '%s' no coincide con su hash guardado (meta.json): se ignora y se vuelve a descargar.",
      ruta_csv
    ), call. = FALSE)
    return(NULL)
  }
  datos <- read_csv(ruta_csv, show_col_types = FALSE)
  list(datos = datos, metadata = metadata)
}

#' Punto de entrada: devuelve los datos diarios de Argentina, usando la
#' cache local si existe (y es consistente) salvo que forzar_descarga=TRUE.
#'
#' @return Lista con `datos` (tibble diario) y `metadata` (procedencia:
#'   descargado ahora, o leido de cache con su fecha/hash originales).
obtener_datos_owid_diario <- function(forzar_descarga = FALSE, dir_cache = "data/cache",
                                        url = "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv") {
  if (!forzar_descarga) {
    cache <- leer_cache_owid(dir_cache)
    if (!is.null(cache)) {
      cache$metadata$origen <- "cache"
      return(cache)
    }
  }
  datos <- descargar_owid_argentina(url)
  metadata <- guardar_cache_owid(datos, url, dir_cache)
  metadata$origen <- "descarga"
  list(datos = datos, metadata = metadata)
}

# ------------------------------------------------------------------------------
# Validacion estructural
# ------------------------------------------------------------------------------

#' Audita los datos diarios: columnas, tipos, NA, valores negativos. No
#' modifica los datos (eso lo hace agregar_owid_semanal); solo reporta.
#'
#' Sobre los valores negativos: OWID puede tener new_cases/new_deaths
#' negativos cuando un pais corrige retroactivamente un conteo acumulado
#' sobreestimado. Sumarlos TAL CUAL en la agregacion semanal es lo correcto
#' (la correccion debe descontarse de alguna semana); lo que no tiene sentido
#' epidemiologico es que el TOTAL semanal de muertes de una semana entera de
#' negativo, eso se trata en agregar_owid_semanal().
#'
#' @return Lista con `ok` (TRUE/FALSE) y `advertencias` (character vector,
#'   vacio si no hay nada que reportar).
validar_owid_diario <- function(datos) {
  advertencias <- character(0)

  faltantes <- setdiff(COLUMNAS_OWID_REQUERIDAS, names(datos))
  if (length(faltantes) > 0) {
    stop(sprintf("Faltan columnas en los datos de OWID: %s", paste(faltantes, collapse = ", ")), call. = FALSE)
  }

  if (!inherits(datos$date, "Date")) {
    stop("La columna 'date' no es de tipo Date.", call. = FALSE)
  }
  for (col in c("new_cases", "new_deaths", "new_tests", "positive_rate", "people_fully_vaccinated", "population")) {
    if (!is.numeric(datos[[col]])) {
      stop(sprintf("La columna '%s' no es numerica (es %s).", col, class(datos[[col]])[1]), call. = FALSE)
    }
  }

  n_na_cases  <- sum(is.na(datos$new_cases))
  n_na_deaths <- sum(is.na(datos$new_deaths))
  if (n_na_cases > 0)  advertencias <- c(advertencias, sprintf("new_cases tiene %d NA.", n_na_cases))
  if (n_na_deaths > 0) advertencias <- c(advertencias, sprintf("new_deaths tiene %d NA.", n_na_deaths))

  n_neg_cases  <- sum(datos$new_cases < 0, na.rm = TRUE)
  n_neg_deaths <- sum(datos$new_deaths < 0, na.rm = TRUE)
  if (n_neg_cases > 0)  advertencias <- c(advertencias, sprintf("new_cases tiene %d dias con valor negativo (correccion retroactiva de OWID).", n_neg_cases))
  if (n_neg_deaths > 0) advertencias <- c(advertencias, sprintf("new_deaths tiene %d dias con valor negativo (correccion retroactiva de OWID).", n_neg_deaths))

  poblaciones_distintas <- length(unique(na.omit(datos$population)))
  if (poblaciones_distintas > 1) {
    advertencias <- c(advertencias, sprintf("La columna 'population' tiene %d valores distintos (se esperaba uno solo).", poblaciones_distintas))
  }

  fechas_duplicadas <- sum(duplicated(datos$date))
  if (fechas_duplicadas > 0) {
    stop(sprintf("Hay %d fechas duplicadas en los datos diarios.", fechas_duplicadas), call. = FALSE)
  }

  list(ok = length(advertencias) == 0, advertencias = advertencias)
}

# ------------------------------------------------------------------------------
# Agregacion semanal
# ------------------------------------------------------------------------------

#' Agrega los datos diarios a semanal, con semana 1 = fecha_inicio_serie
#' (igual criterio que semana_de_fecha() en R/model.R). Detecta semanas
#' incompletas (con menos de 7 fechas distintas) y, por defecto, las
#' DESCARTA de la serie usable en vez de dejarlas contar como una semana
#' entera. Esto importa: con fecha_fin_serie = 2022-06-05 (un domingo), la
#' "semana 118" original tenia un solo dia (el 5 de junio) con 84 muertes
#' registradas ESE DIA, pero quedaba en la tabla como si fueran las muertes
#' de una semana completa -- subestimando por construccion la ultima semana
#' del ultimo tramo de calibracion. Se verifico que la unica semana afectada
#' en la serie real 2020-2022 es la ultima (la primera, con
#' fecha_inicio_serie, tiene sus 7 dias completos).
#'
#' Si el TOTAL semanal de muertes o casos da negativo (una correccion
#' retroactiva grande concentrada en una semana), se deja en 0 con una
#' advertencia: un conteo negativo de muertes/semana no tiene significado
#' epidemiologico utilizable por el modelo.
#'
#' @param incluir_semanas_incompletas Si es FALSE (default), las semanas con
#'   menos de 7 dias de datos se eliminan del resultado (quedan reportadas
#'   en `advertencias` igual). Poner en TRUE para conservarlas.
#' @return Lista con `datos` (tibble semanal) y `advertencias` (character).
agregar_owid_semanal <- function(datos_diario, fecha_inicio_serie, fecha_fin_serie, incluir_semanas_incompletas = FALSE) {
  advertencias <- character(0)

  datos_diario <- datos_diario %>%
    filter(date >= fecha_inicio_serie, date <= fecha_fin_serie) %>%
    mutate(semana = semana_de_fecha(date, fecha_inicio_serie))

  semanal <- datos_diario %>%
    group_by(semana) %>%
    summarise(
      fecha                    = min(date),
      dias_observados           = n_distinct(date),
      casos_owid                = sum(new_cases, na.rm = TRUE),
      muertes_owid              = sum(new_deaths, na.rm = TRUE),
      tests_owid                = sum(new_tests, na.rm = TRUE),
      positividad_owid          = mean(positive_rate, na.rm = TRUE),
      vacunados_completos_owid  = last(na.omit(people_fully_vaccinated)),
      poblacion_owid            = mean(population, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(semana) %>%
    mutate(semana_completa = dias_observados == 7)

  semanas_incompletas <- semanal %>% filter(!semana_completa) %>% pull(semana)
  if (length(semanas_incompletas) > 0) {
    if (incluir_semanas_incompletas) {
      advertencias <- c(advertencias, sprintf(
        "Semanas con menos de 7 dias de datos (se CONSERVAN por incluir_semanas_incompletas=TRUE): %s.",
        paste(semanas_incompletas, collapse = ", ")
      ))
    } else {
      advertencias <- c(advertencias, sprintf(
        "Semanas con menos de 7 dias de datos: %s. Se DESCARTAN de la serie usable (incluir_semanas_incompletas=FALSE).",
        paste(semanas_incompletas, collapse = ", ")
      ))
      semanal <- semanal %>% filter(semana_completa)
    }
  }

  n_negativas <- sum(semanal$muertes_owid < 0)
  if (n_negativas > 0) {
    advertencias <- c(advertencias, sprintf(
      "%d semana(s) con total de muertes negativo (correccion retroactiva grande): se llevan a 0.",
      n_negativas
    ))
    semanal <- semanal %>% mutate(muertes_owid = pmax(muertes_owid, 0))
  }

  list(datos = semanal, advertencias = advertencias)
}

# ------------------------------------------------------------------------------
# sigma: ritmo de vacunacion observado (NO es un parametro de literatura)
# ------------------------------------------------------------------------------

#' Calcula sigma (tasa semanal de vacunacion de susceptibles, S -> V) como el
#' hazard CONSTANTE que, aplicado desde fecha_vacunacion hasta el final de la
#' serie, reproduce la cobertura final de "vacunados completos" observada en
#' OWID: si V(T)/N = 1 - exp(-sigma*T), con T en semanas desde el inicio de
#' la vacunacion, entonces sigma = -ln(1 - V(T)/N) / T.
#'
#' OJO: esta sigma es una aproximacion de UN SOLO numero (hazard constante) a
#' un ritmo de vacunacion real que en los datos varia semana a semana
#' (arranca lento, acelera, despues se aplana). No va a reproducir la
#' TRAYECTORIA semana a semana de vacunados_completos_owid, solo su nivel
#' FINAL (ver verificar_significado_vacunacion() para la comparacion
#' explicita entre ambas series).
#'
#' @param datos_semanal Tibble con columnas semana, vacunados_completos_owid.
#' @param n_poblacion Poblacion total.
#' @param fecha_vacunacion,fecha_inicio_serie Fechas (ver R/config.R).
#' @return sigma (numero, 1/semana).
calcular_sigma_desde_owid <- function(datos_semanal, n_poblacion, fecha_vacunacion, fecha_inicio_serie) {
  semana_inicio_vac <- semana_de_fecha(fecha_vacunacion, fecha_inicio_serie)
  semana_fin        <- max(datos_semanal$semana)
  t_semanas         <- semana_fin - semana_inicio_vac
  if (t_semanas <= 0) {
    stop("fecha_vacunacion esta despues (o en) la ultima semana de la serie: no hay ventana para calcular sigma.", call. = FALSE)
  }
  cobertura_final <- tail(na.omit(datos_semanal$vacunados_completos_owid), 1) / n_poblacion
  if (is.na(cobertura_final) || cobertura_final <= 0 || cobertura_final >= 1) {
    stop(sprintf("Cobertura final de vacunacion fuera de rango (0,1): %s.", cobertura_final), call. = FALSE)
  }
  -log(1 - cobertura_final) / t_semanas
}

#' Compara, semana a semana, la cobertura de vacunados OBSERVADA (OWID,
#' dato) contra la cobertura MODELADA (V(t)/N, resultado de simular con el
#' sigma constante de calcular_sigma_desde_owid()). Son dos cosas
#' conceptualmente distintas -- una es el dato crudo de la campana real, la
#' otra es la trayectoria de un modelo de hazard constante calibrado solo a
#' matchear el nivel final -- y esta funcion lo deja explicito en una sola
#' tabla en vez de asumir que son lo mismo.
#'
#' @param datos_semanal Tibble con semana, vacunados_completos_owid.
#' @param simulacion Tibble con semana, v (del resultado de R/simulate.R).
#' @param n_poblacion Poblacion total.
#' @return Tibble(semana, cobertura_observada, cobertura_modelada, diferencia).
verificar_significado_vacunacion <- function(datos_semanal, simulacion, n_poblacion) {
  datos_semanal %>%
    select(semana, vacunados_completos_owid) %>%
    mutate(cobertura_observada = vacunados_completos_owid / n_poblacion) %>%
    inner_join(simulacion %>% select(semana, v) %>% rename(cobertura_modelada = v), by = "semana") %>%
    mutate(diferencia = cobertura_modelada - cobertura_observada) %>%
    select(semana, cobertura_observada, cobertura_modelada, diferencia)
}
