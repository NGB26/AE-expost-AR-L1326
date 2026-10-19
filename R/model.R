# ==============================================================================
# MODEL - Ecuaciones diferenciales del SIVQRD y cantidades derivadas
# ==============================================================================
# Ecuaciones (23) a (31), tomadas literalmente de "Modelo_SIR_para_PCR.docx",
# seccion 4 ("Modelo extendido SIVQRD"). Estado: fracciones de poblacion
# s, i, iv, v, q, r, d, con s+i+iv+v+q+r+d = 1 en todo momento (prueba mas
# abajo). Se agrega un OCTAVO estado, inf_acum, puramente contable (no forma
# parte de la poblacion, no entra en la suma de conservacion): acumula el
# total de NUEVAS infecciones (de origen susceptible + de origen vacunado),
# para poder calcular la incidencia semanal exacta via diff() igual que se
# hace con las muertes (ver R/simulate.R). Antes, "tasa_contagio_semanal" se
# aproximaba como beta*s_prom, que (a) no es una tasa de nuevas infecciones
# sino una tasa de salida de S por unidad de i+iv implicita, y (b) ignoraba
# por completo las infecciones de origen vacunado (beta_v*v*(i+iv)).
# ==============================================================================

#' Numero de semana (entero) que CONTIENE una fecha dada, con semana 1 =
#' fecha_inicio_serie. APROXIMACION: el modelo opera en grilla semanal, asi
#' que una fecha de corte que cae a mitad de semana queda absorbida por esa
#' semana completa (ver R/config.R, seccion "Fechas de corte de beta").
semana_de_fecha <- function(fecha, fecha_inicio_serie) {
  as.integer(floor(as.numeric(fecha - fecha_inicio_serie) / 7)) + 1L
}

#' Deriva gamma, delta, gamma_v y delta_v a partir de la duracion del
#' periodo infeccioso y el IFR (base y de vacunados).
#'
#' INTERPRETACION DEL IFR: se asume que, mientras un individuo esta en el
#' compartimento I (o Iv), compite entre dos salidas de hazard CONSTANTE:
#' recuperarse (tasa gamma) o morir (tasa delta). Bajo hazards constantes,
#' la probabilidad de que la salida sea "muerte" es delta/(gamma+delta) (la
#' formula estandar de riesgos competitivos). Para que esa probabilidad sea
#' EXACTAMENTE el IFR, y que a la vez gamma+delta reproduzca la duracion
#' total del periodo infeccioso D (gamma_total = 7/D, en 1/semana), hace
#' falta:
#'     delta = IFR * gamma_total
#'     gamma = (1 - IFR) * gamma_total
#' Verificacion algebraica: delta + gamma = IFR*gamma_total + (1-IFR)*gamma_total
#' = gamma_total (duracion correcta) y delta/(gamma+delta) = IFR*gamma_total/gamma_total
#' = IFR (letalidad correcta). Es la UNICA combinacion lineal en gamma_total que
#' cumple ambas igualdades a la vez; ver tests/testthat/test-model.R,
#' "IFR se recupera en un proceso competitivo simple" para la verificacion
#' numerica por simulacion.
#'
#' Mismo razonamiento para los vacunados, con ifr_vacunados = ifr_base*(1-ve_muerte)
#' y dias_infecciosos_vac en vez de dias_infecciosos.
derivar_parametros_biologicos <- function(valores) {
  gamma_total <- 7 / valores$dias_infecciosos
  delta <- valores$ifr_base * gamma_total
  gamma <- (1 - valores$ifr_base) * gamma_total

  ifr_vacunados <- valores$ifr_base * (1 - valores$ve_muerte)
  gamma_total_vac <- 7 / valores$dias_infecciosos_vac
  delta_v <- ifr_vacunados * gamma_total_vac
  gamma_v <- (1 - ifr_vacunados) * gamma_total_vac

  list(
    gamma_total = gamma_total, delta = delta, gamma = gamma,
    ifr_vacunados = ifr_vacunados, gamma_total_vac = gamma_total_vac,
    delta_v = delta_v, gamma_v = gamma_v
  )
}

#' Derivadas del sistema SIVQRD (ecuaciones 23 a 29), para deSolve::ode().
#'
#' @param t Tiempo, en semanas desde fecha_inicio_serie (semana 1 = t=1).
#' @param estado Vector nombrado c(s, i, iv, v, q, r, d, inf_acum).
#' @param parametros Lista con beta_fun(t), theta, sigma, ve_contagio, gamma,
#'   delta, gamma_v, delta_v, semana_inicio_aspo, semana_inicio_vacunacion.
#'   theta y sigma se "activan" recien desde el hito de politica
#'   correspondiente (no hubo cuarentena antes del ASPO ni vacunas antes de
#'   fecha_vacunacion): es un supuesto de implementacion, no parte de las
#'   ecuaciones en si.
derivadas_sivqrd <- function(t, estado, parametros) {
  with(as.list(c(estado, parametros)), {

    beta_t   <- beta_fun(t)
    beta_v_t <- beta_t * (1 - ve_contagio)
    theta_t  <- if (t >= semana_inicio_aspo) theta else 0
    sigma_t  <- if (t >= semana_inicio_vacunacion) sigma else 0

    nuevas_infecciones_s <- beta_t * s * (i + iv)      # flujo S  -> I  (termino comun a ec. 23 y 24)
    nuevas_infecciones_v <- beta_v_t * v * (i + iv)    # flujo V  -> Iv (termino comun a ec. 25 y 26)

    d_s  <- -nuevas_infecciones_s - sigma_t * s                       # Ecuacion (23)
    d_i  <-  nuevas_infecciones_s - i * (theta_t + delta + gamma)     # Ecuacion (24)
    d_iv <-  nuevas_infecciones_v - iv * (delta_v + gamma_v)          # Ecuacion (25)
    d_v  <-  sigma_t * s - nuevas_infecciones_v                       # Ecuacion (26)
    d_q  <-  theta_t * i - (gamma + delta) * q                        # Ecuacion (27)
    d_r  <-  gamma * (i + q) + gamma_v * iv                           # Ecuacion (28)
    d_d  <-  delta * (i + q) + delta_v * iv                           # Ecuacion (29)

    d_inf_acum <- nuevas_infecciones_s + nuevas_infecciones_v         # auxiliar: NO es parte de s+i+iv+v+q+r+d

    list(c(s = d_s, i = d_i, iv = d_iv, v = d_v, q = d_q, r = d_r, d = d_d, inf_acum = d_inf_acum))
  })
}

# ------------------------------------------------------------------------------
# PRUEBA DE CONSERVACION (d/dt (s+i+iv+v+q+r+d) = 0)
# ------------------------------------------------------------------------------
# Sumando las 7 derivadas (sin inf_acum):
#   -nuevas_infecciones_s - sigma_t*s                                  (d_s)
#   +nuevas_infecciones_s - i*(theta_t+delta+gamma)                    (d_i)
#   +nuevas_infecciones_v - iv*(delta_v+gamma_v)                       (d_iv)
#   +sigma_t*s - nuevas_infecciones_v                                  (d_v)
#   +theta_t*i - (gamma+delta)*q                                       (d_q)
#   +gamma*(i+q) + gamma_v*iv                                          (d_r)
#   +delta*(i+q) + delta_v*iv                                          (d_d)
#
# Los terminos "nuevas_infecciones_s" se cancelan entre d_s y d_i; los
# "nuevas_infecciones_v" entre d_v y d_iv; "sigma_t*s" entre d_s y d_v. Quedan:
#   -i*(theta_t+delta+gamma) - iv*(delta_v+gamma_v) + theta_t*i - (gamma+delta)*q
#   + gamma*(i+q) + gamma_v*iv + delta*(i+q) + delta_v*iv
# Agrupando por compartimento de origen:
#   i:  -i*theta_t - i*delta - i*gamma + theta_t*i + gamma*i + delta*i = 0
#   q:  -(gamma+delta)*q + gamma*q + delta*q = 0
#   iv: -iv*delta_v - iv*gamma_v + gamma_v*iv + delta_v*iv = 0
# Los tres grupos se anulan exactamente: d/dt(s+i+iv+v+q+r+d) = 0 para
# cualquier valor de los parametros. Verificado tambien numericamente en
# tests/testthat/test-model.R ("la poblacion total se conserva").
# ------------------------------------------------------------------------------

#' R0 (ec. 30) via la matriz de proxima generacion (NGM).
#'
#' Compartimentos "infectados" (los unicos que generan nuevas infecciones):
#' i, iv. Matriz de nuevas infecciones F (derivadas de las tasas de entrada
#' respecto de i, iv) y matriz de transiciones V (derivadas de las tasas de
#' salida, sin contar nuevas infecciones), evaluadas en un punto (s*, v*):
#'   F = [ beta*s*      beta*s*     ]      V = [ theta+gamma+delta   0           ]
#'       [ beta_v*v*    beta_v*v*   ]          [ 0                   gamma_v+delta_v ]
#' (las dos columnas de F son iguales porque tanto los infectados de origen S
#' como los de origen V contagian a traves del mismo "pool" i+iv).
#'
#' F es RANGO 1 (es el producto externo de los vectores (beta*s*, beta_v*v*)
#' y (1,1)): su unico autovalor no nulo de F*V^-1 es su traza,
#'   R = beta*s*/(theta+gamma+delta) + beta_v*v*/(gamma_v+delta_v)
#' que es exactamente la ecuacion (31) del documento. R0 (ec. 30) es el caso
#' particular s*=1, v*=0 (nadie vacunado todavia): R0 = beta/(theta+gamma+delta).
#' Ver calcular_re_ngm_numerico() para la verificacion por autovalores.
calcular_r0 <- function(beta_valor, theta_valor, gamma, delta) {
  beta_valor / (theta_valor + gamma + delta)
}

#' Re(t) (ec. 31) via la NGM, forma cerrada (ver calcular_r0() para la
#' derivacion). s_val y v_val son las fracciones de poblacion susceptible y
#' vacunada en el punto donde se evalua Re (puede ser un instante puntual o
#' un promedio temporal, segun lo que llame a esta funcion).
calcular_re <- function(beta_valor, theta_valor, s_val, v_val, gamma, delta, ve_contagio, gamma_v, delta_v) {
  beta_v_valor <- beta_valor * (1 - ve_contagio)
  (beta_valor * s_val) / (theta_valor + gamma + delta) +
    (beta_v_valor * v_val) / (gamma_v + delta_v)
}

#' Mismo Re(t), calculado NUMERICAMENTE como el radio espectral de F*V^-1
#' (autovalor dominante), para verificar por fuera de la formula cerrada que
#' calcular_re() realmente corresponde a la NGM. Se usa en
#' tests/testthat/test-model.R, no en el pipeline principal (es mas lento y
#' matematicamente redundante con calcular_re(), que ya esta demostrada).
calcular_re_ngm_numerico <- function(beta_valor, theta_valor, s_val, v_val, gamma, delta, ve_contagio, gamma_v, delta_v) {
  beta_v_valor <- beta_valor * (1 - ve_contagio)
  matriz_F <- matrix(
    c(beta_valor * s_val, beta_v_valor * v_val,
      beta_valor * s_val, beta_v_valor * v_val),
    nrow = 2, byrow = TRUE
  )
  matriz_V <- diag(c(theta_valor + gamma + delta, gamma_v + delta_v))
  ngm <- matriz_F %*% solve(matriz_V)
  max(Mod(eigen(ngm, only.values = TRUE)$values))
}
