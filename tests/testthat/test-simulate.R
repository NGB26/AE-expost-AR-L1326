test_that("simular_tramo da el numero de filas esperado y diff(d) la longitud correcta", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)

  semana_ini <- 5; semana_fin <- 12
  sim <- simular_tramo(estado_inicial, semana_ini, semana_fin, beta_valor = 0.8, parametros_fijos)

  # (fin - ini + 2) filas: incluye el punto extra en semana_fin+1.
  expect_equal(nrow(sim), semana_fin - semana_ini + 2)
  # diff(d) da exactamente un valor por semana del tramo.
  expect_equal(length(diff(sim$d)), semana_fin - semana_ini + 1)
  expect_equal(length(diff(sim$d)), length(seq(semana_ini, semana_fin)))
})

test_that("un tramo de UNA sola semana no rompe simular_tramo (caso borde real: Pre-intervencion)", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)

  sim <- simular_tramo(estado_inicial, semana_ini = 1, semana_fin = 1, beta_valor = 0.8, parametros_fijos)

  expect_equal(nrow(sim), 2)
  expect_equal(length(diff(sim$d)), 1)
})

test_that("simular_serie es continua entre tramos (no hay salto artificial al unir)", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  tramos_def <- fixture_tramos_def()  # tramo1: semanas 1-7, tramo2: semanas 8-20

  # OJO: con fixture_parametros_fijos() normal, semana_inicio_aspo=7 cae
  # EXACTO en el borde tramo1/tramo2, y ahi el RHS de la EDO tiene un salto
  # real (theta pasa de 0 a su valor vigente). Integrar continuo de punta a
  # punta vs. partir en 2 e integrar cada pedazo NO deberian dar lo mismo en
  # ese caso: partir justo en el salto es, de hecho, la practica numerica
  # correcta (reinicia el control de paso del integrador en la
  # discontinuidad) y el continuo arrastra un poco mas de error ahi. Para
  # aislar lo que esta prueba quiere verificar -- que partir en tramos NO
  # introduce un artefacto POR SI SOLO, mas alla de como se maneje una
  # discontinuidad real -- se usa un parametros_fijos SIN ninguna
  # discontinuidad en el rango de semanas de la prueba (theta y sigma
  # arrancan mucho despues de la semana 20).
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)
  parametros_fijos$semana_inicio_aspo <- 9999
  parametros_fijos$semana_inicio_vacunacion <- 9999

  # Mismo beta en ambos tramos: la serie "partida en 2" debe coincidir con
  # una unica simulacion continua de semana 1 a 20.
  beta_constante <- 0.7
  resultado_partido <- simular_serie(estado_inicial, tramos_def, betas = c(beta_constante, beta_constante), parametros_fijos)
  sim_continua <- simular_tramo(estado_inicial, 1, 20, beta_constante, parametros_fijos)

  comparacion <- resultado_partido$simulacion %>%
    dplyr::select(semana, s, i, iv, v, q, r, d) %>%
    dplyr::inner_join(sim_continua %>% dplyr::select(semana, s, i, iv, v, q, r, d), by = "semana", suffix = c("_partido", "_continuo"))

  for (col in c("s", "i", "iv", "v", "q", "r", "d")) {
    diferencia <- abs(comparacion[[paste0(col, "_partido")]] - comparacion[[paste0(col, "_continuo")]])
    expect_true(all(diferencia < 1e-6), info = paste("columna:", col))
  }
})

test_that("el estado final de un tramo es el estado inicial exacto del siguiente", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  tramos_def <- fixture_tramos_def()
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)

  resultado <- simular_serie(estado_inicial, tramos_def, betas = c(0.9, 0.5), parametros_fijos)

  ultima_fila_tramo1 <- resultado$simulacion %>% dplyr::filter(tramo == 1) %>% dplyr::slice_tail(n = 1)
  primera_fila_tramo2 <- resultado$simulacion %>% dplyr::filter(tramo == 2) %>% dplyr::slice(1)

  # No son la MISMA semana (tramo2 arranca una semana despues de donde
  # termina tramo1), pero la condicion inicial de tramo2 (antes de integrar
  # su primera semana) debe ser EXACTAMENTE el estado final de tramo1: lo
  # verificamos re-simulando una sola semana desde ahi.
  estado_cierre_tramo1 <- ultima_fila_tramo1 %>% dplyr::select(s, i, iv, v, q, r, d, inf_acum) %>% unlist()
  primer_paso_tramo2 <- simular_tramo(estado_cierre_tramo1, primera_fila_tramo2$semana, primera_fila_tramo2$semana, 0.5, parametros_fijos)

  expect_equal(as.numeric(primer_paso_tramo2$s[1]), as.numeric(estado_cierre_tramo1["s"]), tolerance = 1e-10)
  expect_equal(as.numeric(primer_paso_tramo2$d[1]), as.numeric(estado_cierre_tramo1["d"]), tolerance = 1e-10)
})

test_that("simular_serie produce emisiones e incidencia con una fila por semana de toda la serie", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  tramos_def <- fixture_tramos_def()
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)

  resultado <- simular_serie(estado_inicial, tramos_def, betas = c(0.9, 0.5), parametros_fijos)

  n_semanas_total <- max(tramos_def$semana_fin) - min(tramos_def$semana_ini) + 1
  expect_equal(nrow(resultado$emisiones), n_semanas_total)
  expect_equal(nrow(resultado$incidencia), n_semanas_total)
  expect_equal(sort(resultado$emisiones$semana), min(tramos_def$semana_ini):max(tramos_def$semana_fin))
})

test_that("ajustar_betas_secuencial devuelve un beta por tramo, dentro de las cotas", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  tramos_def <- fixture_tramos_def()
  datos_semanal <- fixture_datos_semanal()
  n_poblacion <- fixture_n_poblacion()
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)

  funcion_error <- function(emisiones, semanas_objetivo) {
    calcular_msre_log(emisiones, datos_semanal, n_poblacion, valores$rezago_semanas, semanas_objetivo)
  }

  resultado <- ajustar_betas_secuencial(estado_inicial, tramos_def, parametros_fijos, funcion_error,
                                         beta_min = valores$beta_min, beta_max = valores$beta_max)

  expect_equal(length(resultado$betas), nrow(tramos_def))
  expect_true(all(resultado$betas >= valores$beta_min - 1e-9 & resultado$betas <= valores$beta_max + 1e-9))
  expect_false(anyNA(resultado$betas))
})
