test_that("el rezago no pierde la historia del tramo anterior (regresion del bug real)", {
  # Reproduce exactamente el bug encontrado y corregido en una version
  # anterior del script: si se filtra a semanas_objetivo ANTES de aplicar el
  # rezago, las primeras `rezago_semanas` del tramo nuevo pierden el
  # historial del tramo anterior y el lag se rellena con el default (0) en
  # vez de con las emisiones reales de las ultimas semanas del tramo previo.
  rezago <- 2
  emisiones <- tibble::tibble(semana = 1:10, muertes_emitidas = (1:10) / 1e6)  # fraccion de poblacion
  n_poblacion <- 1e6
  datos_semanal <- tibble::tibble(semana = 1:10, muertes_owid = 1:10)

  # Semana 8 (en el "tramo nuevo") deberia tomar el rezago de la semana
  # 8 - rezago = 6, cuyo valor crudo es muertes_emitidas=6/1e6 -> 6 muertes.
  resultado_correcto <- alinear_muertes_con_rezago(emisiones, datos_semanal, n_poblacion, rezago, semanas_objetivo = 8:10)
  valor_en_semana_8 <- resultado_correcto$muertes_modelo[resultado_correcto$semana == 8]

  expect_equal(valor_en_semana_8, 6)  # NO debe ser 0 (que es lo que daria el bug)

  # Version BUGUEADA para contraste expreso: filtrar antes de aplicar el lag.
  emisiones_filtradas_antes <- emisiones %>% dplyr::filter(semana %in% 8:10)
  resultado_bugueado <- alinear_muertes_con_rezago(emisiones_filtradas_antes, datos_semanal, n_poblacion, rezago)
  valor_bugueado_semana_8 <- resultado_bugueado$muertes_modelo[resultado_bugueado$semana == 8]

  expect_equal(valor_bugueado_semana_8, 0)  # el comportamiento bugueado, documentado para que quede claro el contraste
  expect_false(isTRUE(all.equal(valor_en_semana_8, valor_bugueado_semana_8)))
})

test_that("calcular_msre_log no produce NA cuando los datos no tienen semanas en cero", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  tramos_def <- fixture_tramos_def()
  datos_semanal <- fixture_datos_semanal()
  n_poblacion <- fixture_n_poblacion()
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)

  resultado <- simular_serie(estado_inicial, tramos_def, betas = c(0.8, 0.5), parametros_fijos)
  comparacion <- alinear_muertes_con_rezago(resultado$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas)

  expect_false(anyNA(comparacion$muertes_owid))
  expect_false(anyNA(comparacion$muertes_modelo))

  error <- calcular_msre_log(resultado$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas)
  expect_true(is.finite(error))
  expect_false(is.na(error))
})

test_that("calcular_msre_log excluye (no rompe) semanas con 0 muertes observadas", {
  emisiones <- tibble::tibble(semana = 1:5, muertes_emitidas = c(0, 0.000001, 0.000002, 0.000001, 0))
  n_poblacion <- 1e6
  datos_semanal <- tibble::tibble(semana = 1:5, muertes_owid = c(0, 1, 2, 1, 0))

  error <- calcular_msre_log(emisiones, datos_semanal, n_poblacion, rezago_semanas = 0)

  expect_true(is.finite(error))  # no debe dar NA ni Inf por las semanas en 0
})

test_that("calcular_residuos etiqueta cada semana con su tramo", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  tramos_def <- fixture_tramos_def()
  datos_semanal <- fixture_datos_semanal()
  n_poblacion <- fixture_n_poblacion()
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)

  resultado <- simular_serie(estado_inicial, tramos_def, betas = c(0.8, 0.5), parametros_fijos)
  residuos <- calcular_residuos(resultado$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas, tramos_def)

  expect_false(anyNA(residuos$tramo))
  expect_true(all(residuos$tramo %in% tramos_def$tramo))
})
