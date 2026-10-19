test_that("el pipeline completo (datos sinteticos -> modelo -> calibracion -> metricas) corre sin internet", {
  valores <- fixture_valores()
  n_poblacion <- fixture_n_poblacion()
  datos_semanal <- fixture_datos_semanal(n_poblacion = n_poblacion)
  tramos_def <- fixture_tramos_def()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores, n_poblacion)

  semana_aspo <- semana_de_fecha(valores$fecha_aspo, valores$fecha_inicio_serie)
  semana_vac  <- semana_de_fecha(valores$fecha_vacunacion, valores$fecha_inicio_serie)

  resultado <- calibrar_secuencial_alternante(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_aspo, semana_inicio_vacunacion = semana_vac,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    theta_inicial = valores$theta_inicial, theta_min = valores$theta_min, theta_max = valores$theta_max,
    beta_min = valores$beta_min, beta_max = valores$beta_max,
    max_iteraciones = valores$max_iteraciones_externas, tolerancia = valores$tolerancia_convergencia
  )

  expect_equal(length(resultado$betas), nrow(tramos_def))
  expect_false(anyNA(resultado$betas))
  expect_true(is.finite(resultado$theta))
  expect_true(is.finite(resultado$mejor_error))

  mse_nivel <- calcular_mse_nivel(resultado$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas)
  mse_log   <- calcular_mse_log(resultado$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas)
  expect_true(is.finite(mse_nivel))
  expect_true(is.finite(mse_log))

  tasas_mort <- calcular_tasa_mortalidad_por_tramo(resultado$emisiones, n_poblacion, tramos_def)
  tasas_cont <- calcular_tasa_contagio_por_tramo(resultado$incidencia, tramos_def)
  expect_equal(nrow(tasas_mort), nrow(tramos_def))
  expect_equal(nrow(tasas_cont), nrow(tramos_def))
  expect_false(anyNA(tasas_mort$tasa_mortalidad_semanal))
  expect_false(anyNA(tasas_cont$tasa_contagio_semanal))

  resumen_tramos <- tramos_def %>%
    dplyr::mutate(beta = resultado$betas) %>%
    dplyr::left_join(tasas_mort, by = "tramo") %>%
    dplyr::left_join(tasas_cont, by = "tramo")
  expect_false(anyNA(resumen_tramos))
})

test_that("calibrar_conjunto_global corre sobre datos sinteticos y mejora o iguala al baseline", {
  valores <- fixture_valores()
  n_poblacion <- fixture_n_poblacion()
  datos_semanal <- fixture_datos_semanal(n_poblacion = n_poblacion)
  tramos_def <- fixture_tramos_def()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores, n_poblacion)
  semana_aspo <- semana_de_fecha(valores$fecha_aspo, valores$fecha_inicio_serie)
  semana_vac  <- semana_de_fecha(valores$fecha_vacunacion, valores$fecha_inicio_serie)

  baseline <- calibrar_secuencial_alternante(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_aspo, semana_inicio_vacunacion = semana_vac,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    theta_inicial = valores$theta_inicial, theta_min = valores$theta_min, theta_max = valores$theta_max,
    beta_min = valores$beta_min, beta_max = valores$beta_max, max_iteraciones = 2
  )

  conjunto <- calibrar_conjunto_global(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_aspo, semana_inicio_vacunacion = semana_vac,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    beta_min = valores$beta_min, beta_max = valores$beta_max, theta_min = valores$theta_min, theta_max = valores$theta_max,
    valores_iniciales = c(baseline$betas, baseline$theta)
  )

  expect_equal(length(conjunto$betas), nrow(tramos_def))
  expect_true(is.finite(conjunto$mejor_error))
  # Arrancando DESDE el optimo del baseline, el conjunto no deberia terminar
  # peor (L-BFGS-B es, como minimo, un refinamiento local de ese punto).
  expect_lte(conjunto$mejor_error, baseline$mejor_error + 1e-6)
})
