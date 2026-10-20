# ==============================================================================
# RATES - Tasas IMPLICITAS por semana (reconstruccion tipo "Markov en Excel")
# ==============================================================================
# Esta es la "estructura de datos que permite calcular las tasas (beta, sigma,
# gamma, delta, Re...) implicitas para cada semana" que pide la especificacion.
#
# Idea, igual que una planilla de Markov del SIR:
#   1. Se reconstruyen, semana a semana, los STOCKS de cada compartimento
#      (S, I, R, D, V) a partir de los flujos observados en OWID (casos nuevos,
#      muertes nuevas, vacunados acumulados), con una recursion discreta de
#      paso semanal (un Euler explicito, como el que haria una planilla).
#   2. De esos stocks y flujos se DESPEJAN las tasas que, en las ecuaciones del
#      modelo, generarian esos mismos flujos. Eso da un beta, un sigma, un
#      delta, etc. IMPLICITOS para cada semana.
#
# PARA QUE SIRVE: (a) diagnostico -- ver como vario de hecho la transmision
# (beta) semana a semana en los datos crudos, sin imponer los 5 tramos; (b)
# puntos de partida y rangos razonables para la calibracion; (c) chequear que
# los parametros fijos de literatura (gamma, delta) no sean absurdos frente a
# lo que muestran los datos.
#
# ADVERTENCIAS IMPORTANTES (la estructura las documenta, no las esconde):
#   * Se usan los CASOS REPORTADOS como proxy de las infecciones nuevas. La
#     subdeteccion (sobre todo en 2020 y en asintomaticos) hace que estas
#     infecciones esten SUBESTIMADAS; por lo tanto el beta implicito y el stock
#     de infecciosos tambien. Son tasas "sobre casos detectados", no verdad de
#     terreno. Para corregir subdeteccion habria que dividir casos por una tasa
#     de ascertainment (no disponible aca) -- queda como # EDITAR si se tiene.
#   * No se modela el solapamiento entre "vacunados" y "ya infectados/
#     recuperados": a S se le restan por separado infecciones acumuladas y
#     vacunados acumulados, lo que puede doble-contar personas y llevar S por
#     debajo de 0 en semanas tardias; se clampea a 0 y se avisa.
#   * gamma y delta del modelo son CONSTANTES biologicas; aca ademas se reporta
#     un delta y un CFR "implicitos" (muertes / stock infeccioso, muertes /
#     casos) como contraste empirico, no como insumo del modelo.
# ==============================================================================

library(dplyr)
library(tibble)

#' Reconstruye stocks semanales y despeja las tasas implicitas de cada semana.
#'
#' @param datos_semanal Tibble de agregar_owid_semanal(): columnas semana,
#'   fecha, casos_owid, muertes_owid, vacunados_completos_owid, poblacion_owid
#'   (y las demas, que no se usan aca).
#' @param n_poblacion Poblacion total (constante).
#' @param dias_infecciosos Duracion del periodo infeccioso (dias), para derivar
#'   gamma_total = 7/dias_infecciosos (1/semana), la tasa de salida de I.
#' @param ifr_base IFR de literatura, para partir gamma_total en gamma (recup.)
#'   y delta (muerte): delta = ifr*gamma_total, gamma = (1-ifr)*gamma_total.
#' @param i0_personas Semilla de infecciosos en la semana 1 (si NULL, usa los
#'   casos de la primera semana como stock inicial).
#' @return Lista con:
#'   - tasas: tibble(semana, fecha, stocks en conteo y fraccion, y tasas
#'     implicitas beta/sigma/delta/cfr/Re por semana).
#'   - advertencias: character vector (p.ej. semanas con S clampeado a 0).
construir_tasas_implicitas_semanales <- function(datos_semanal, n_poblacion,
                                                  dias_infecciosos, ifr_base,
                                                  i0_personas = NULL) {
  advertencias <- character(0)

  d <- datos_semanal %>% arrange(semana)
  n_sem <- nrow(d)

  gamma_total <- 7 / dias_infecciosos          # 1/semana: tasa total de salida de I
  delta_bio   <- ifr_base * gamma_total        # 1/semana: componente muerte
  gamma_bio   <- (1 - ifr_base) * gamma_total  # 1/semana: componente recuperacion

  # Flujos observados (proxy): infecciones nuevas = casos reportados; muertes
  # nuevas = muertes reportadas; vacunados nuevos = diferencia del acumulado.
  nuevas_infecciones <- pmax(d$casos_owid, 0)
  muertes_nuevas     <- pmax(d$muertes_owid, 0)
  vac_acum <- d$vacunados_completos_owid
  vac_acum[is.na(vac_acum)] <- 0
  vac_acum <- cummax(vac_acum)                 # el acumulado no deberia bajar; forzar monotonia
  vac_nuevos <- pmax(c(vac_acum[1], diff(vac_acum)), 0)

  # Reconstruccion de stocks (recursion semanal, Euler explicito):
  #   I_t = I_{t-1} + infecciones_t - (gamma_total) * I_{t-1}
  # El outflow de I lo fija gamma_total (literatura); las muertes observadas se
  # usan para el stock D acumulado, no para vaciar I (asi gamma_total manda la
  # duracion, como en el modelo).
  I <- numeric(n_sem); R <- numeric(n_sem); D <- numeric(n_sem)
  V <- numeric(n_sem); S <- numeric(n_sem)
  infecciones_acum <- cumsum(nuevas_infecciones)

  I[1] <- if (is.null(i0_personas)) nuevas_infecciones[1] else i0_personas
  D[1] <- muertes_nuevas[1]
  R[1] <- gamma_bio * I[1]
  V[1] <- vac_acum[1]
  S[1] <- max(n_poblacion - I[1] - R[1] - D[1] - V[1], 0)

  for (t in 2:n_sem) {
    salida_I <- gamma_total * I[t - 1]
    I[t] <- max(I[t - 1] + nuevas_infecciones[t] - salida_I, 0)
    D[t] <- D[t - 1] + muertes_nuevas[t]
    R[t] <- R[t - 1] + gamma_bio * I[t - 1]     # recuperados acumulados (resto del outflow)
    V[t] <- vac_acum[t]
    S[t] <- n_poblacion - infecciones_acum[t] - V[t]
    if (S[t] < 0) S[t] <- 0
  }

  n_clamp <- sum(S == 0 & seq_len(n_sem) > 1)
  if (n_clamp > 0) {
    advertencias <- c(advertencias, sprintf(
      "%d semana(s) con S reconstruido <= 0 (clampeado a 0): tipicamente por doble conteo vacunados/recuperados o subdeteccion; el beta implicito de esas semanas no es confiable.",
      n_clamp
    ))
  }

  s_frac <- S / n_poblacion
  i_frac <- I / n_poblacion
  v_frac <- V / n_poblacion

  # --- Despeje de las tasas implicitas (una por semana) ----------------------
  # beta_t:   infecciones_t = beta_t * s_{t-1} * i_{t-1} * N  =>
  #           beta_t = infecciones_t / (s_{t-1} * i_{t-1} * N)
  # sigma_t:  vac_nuevos_t = sigma_t * S_{t-1}  =>  sigma_t = vac_nuevos_t / S_{t-1}
  # delta_t:  muertes_t = delta_t * I_{t-1}  =>  delta_t = muertes_t / I_{t-1}
  # cfr_t:    muertes_t / casos_t  (letalidad sobre casos detectados)
  # Re_t:     beta_t * s_{t-1} / gamma_total  (proxy del numero reproductivo efectivo)
  s_prev <- c(NA, s_frac[-n_sem])
  i_prev <- c(NA, i_frac[-n_sem])
  S_prev <- c(NA, S[-n_sem])
  I_prev <- c(NA, I[-n_sem])

  denom_beta <- s_prev * i_prev * n_poblacion
  beta_implicita  <- ifelse(denom_beta > 0, nuevas_infecciones / denom_beta, NA_real_)
  sigma_implicita <- ifelse(S_prev > 0, vac_nuevos / S_prev, NA_real_)
  delta_implicita <- ifelse(I_prev > 0, muertes_nuevas / I_prev, NA_real_)
  cfr_detectada   <- ifelse(nuevas_infecciones > 0, muertes_nuevas / nuevas_infecciones, NA_real_)
  re_implicito    <- beta_implicita * s_prev / gamma_total

  tasas <- tibble(
    semana = d$semana,
    fecha  = d$fecha,
    # stocks reconstruidos (conteo)
    S = S, I = I, R = R, D = D, V = V,
    # stocks (fraccion de poblacion)
    s = s_frac, i = i_frac, r = R / n_poblacion, d = D / n_poblacion, v = v_frac,
    # flujos observados usados
    infecciones_nuevas = nuevas_infecciones,
    muertes_nuevas = muertes_nuevas,
    vacunados_nuevos = vac_nuevos,
    # tasas implicitas por semana
    beta_implicita  = beta_implicita,
    sigma_implicita = sigma_implicita,
    delta_implicita = delta_implicita,
    cfr_detectada   = cfr_detectada,
    re_implicito    = re_implicito,
    # constantes de literatura, de referencia (no varian por semana)
    gamma_total_bio = gamma_total,
    gamma_bio = gamma_bio,
    delta_bio = delta_bio
  )

  list(tasas = tasas, advertencias = advertencias)
}

#' Resume las tasas implicitas por tramo de calibracion (promedio de beta_t,
#' sigma_t, Re_t, etc. dentro de cada tramo), para comparar de un vistazo el
#' beta IMPLICITO promedio de cada tramo contra el beta CALIBRADO del modelo.
#'
#' @param tasas Tibble devuelto por construir_tasas_implicitas_semanales()$tasas.
#' @param tramos_def Tibble(tramo, nombre, semana_ini, semana_fin).
#' @return Tibble(tramo, nombre, semana_ini, semana_fin, beta_implicita_prom, ...).
resumir_tasas_por_tramo <- function(tasas, tramos_def) {
  asignar_tramo <- function(s) {
    fila <- which(s >= tramos_def$semana_ini & s <= tramos_def$semana_fin)
    if (length(fila) == 0) NA_integer_ else tramos_def$tramo[fila[1]]
  }
  tasas %>%
    mutate(tramo = vapply(semana, asignar_tramo, integer(1))) %>%
    filter(!is.na(tramo)) %>%
    group_by(tramo) %>%
    summarise(
      beta_implicita_prom  = mean(beta_implicita,  na.rm = TRUE),
      sigma_implicita_prom = mean(sigma_implicita, na.rm = TRUE),
      delta_implicita_prom = mean(delta_implicita, na.rm = TRUE),
      cfr_detectada_prom   = mean(cfr_detectada,   na.rm = TRUE),
      re_implicito_prom    = mean(re_implicito,    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(tramos_def, by = "tramo") %>%
    arrange(tramo) %>%
    select(tramo, nombre, semana_ini, semana_fin, everything())
}
