test_that("el pipeline completo (datos sinteticos -> modelo -> calibracion global -> metricas) corre sin internet", {
  valores <- fixture_valores()
  n_poblacion <- fixture_n_poblacion()
  datos_semanal <- fixture_datos_semanal(n_poblacion = n_poblacion)
  tramos_def <- fixture_tramos_def()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores, n_poblacion)

  semana_aspo <- semana_de_fecha(valores$fecha_aspo, valores$fecha_inicio_serie)
  semana_vac  <- semana_de_fecha(valores$fecha_vacunacion, valores$fecha_inicio_serie)

  resultado <- calibrar_beta_global(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    theta = valores$theta, sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_aspo, semana_inicio_vacunacion = semana_vac,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    beta_min = valores$beta_min, beta_max = valores$beta_max, maxit = 100
  )

  expect_equal(length(resultado$betas), nrow(tramos_def))
  expect_false(anyNA(resultado$betas))
  expect_true(all(resultado$betas >= valores$beta_min - 1e-9 & resultado$betas <= valores$beta_max + 1e-9))
  expect_equal(resultado$theta, valores$theta)  # theta FIJO: entra y sale igual, no se calibra
  expect_true(is.finite(resultado$mejor_error))

  tasas_mort <- calcular_tasa_mortalidad_por_tramo(resultado$emisiones, n_poblacion, tramos_def)
  tasas_cont <- calcular_tasa_contagio_por_tramo(resultado$incidencia, tramos_def)
  expect_equal(nrow(tasas_mort), nrow(tramos_def))
  expect_equal(nrow(tasas_cont), nrow(tramos_def))
  expect_false(anyNA(tasas_mort$tasa_mortalidad_semanal))
  expect_false(anyNA(tasas_cont$tasa_contagio_semanal))
})

test_that("calibrar_beta_global iguala o mejora al secuencial (arranca desde su optimo)", {
  valores <- fixture_valores()
  n_poblacion <- fixture_n_poblacion()
  datos_semanal <- fixture_datos_semanal(n_poblacion = n_poblacion)
  tramos_def <- fixture_tramos_def()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores, n_poblacion)
  semana_aspo <- semana_de_fecha(valores$fecha_aspo, valores$fecha_inicio_serie)
  semana_vac  <- semana_de_fecha(valores$fecha_vacunacion, valores$fecha_inicio_serie)

  secuencial <- calibrar_secuencial(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    theta = valores$theta, sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_aspo, semana_inicio_vacunacion = semana_vac,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    beta_min = valores$beta_min, beta_max = valores$beta_max
  )

  global <- calibrar_beta_global(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    theta = valores$theta, sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_aspo, semana_inicio_vacunacion = semana_vac,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    beta_min = valores$beta_min, beta_max = valores$beta_max,
    betas_iniciales = secuencial$betas, maxit = 200
  )

  expect_equal(length(global$betas), nrow(tramos_def))
  expect_true(is.finite(global$mejor_error))
  # Arrancando DESDE el optimo secuencial, el global no deberia terminar peor:
  # minimiza el MISMO loss (RMSE-log) pero sobre toda la serie a la vez.
  expect_lte(global$mejor_error, secuencial$mejor_error + 1e-6)
})

test_that("las tasas implicitas semanales se calculan y resumen por tramo", {
  n_poblacion <- fixture_n_poblacion()
  datos_semanal <- fixture_datos_semanal(n_poblacion = n_poblacion)
  tramos_def <- fixture_tramos_def()

  impl <- construir_tasas_implicitas_semanales(
    datos_semanal, n_poblacion, dias_infecciosos = 11, ifr_base = 0.01, i0_personas = 50
  )

  expect_equal(nrow(impl$tasas), nrow(datos_semanal))
  # Columnas de tasas implicitas presentes
  expect_true(all(c("beta_implicita", "sigma_implicita", "delta_implicita", "re_implicito") %in% names(impl$tasas)))
  # Los stocks reconstruidos son no-negativos (S clampeado a 0 si hiciera falta)
  expect_true(all(impl$tasas$S >= 0))
  expect_true(all(impl$tasas$I >= 0))
  # La primera semana no tiene beta implicito (no hay semana previa)
  expect_true(is.na(impl$tasas$beta_implicita[1]))
  # beta implicito no-NA debe ser no-negativo
  betas_impl <- impl$tasas$beta_implicita[!is.na(impl$tasas$beta_implicita)]
  expect_true(all(betas_impl >= 0))

  resumen <- resumir_tasas_por_tramo(impl$tasas, tramos_def)
  expect_equal(nrow(resumen), nrow(tramos_def))
  expect_true(all(c("beta_implicita_prom", "re_implicito_prom") %in% names(resumen)))
})
