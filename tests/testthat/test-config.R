test_that("validar_config falla explicito cuando sigma es NA (antes de calcularlo)", {
  cfg <- construir_config()
  expect_error(validar_config(cfg, exigir_sigma = TRUE), "sigma.*no tiene valor")
})

test_that("validar_config pasa sin sigma cuando exigir_sigma=FALSE", {
  cfg <- construir_config()
  expect_true(validar_config(cfg, exigir_sigma = FALSE))
})

test_that("validar_config frena ante un parametro fuera de rango", {
  cfg <- construir_config()
  cfg$ifr_base$valor <- 1.5  # imposible: IFR > 1
  expect_error(validar_config(cfg, exigir_sigma = FALSE), "fuera del rango")
})

test_that("validar_config frena ante theta_min >= theta_max", {
  cfg <- construir_config()
  cfg$theta_min$valor <- 1
  cfg$theta_max$valor <- 1
  expect_error(validar_config(cfg, exigir_sigma = FALSE), "theta_min")
})

test_that("validar_config frena ante fechas de corte fuera de orden", {
  cfg <- construir_config()
  cfg$fecha_dispo$valor <- cfg$fecha_aspo$valor - 1  # DISPO antes que ASPO: invalido
  expect_error(validar_config(cfg, exigir_sigma = FALSE), "fecha_aspo")
})

test_that("aplanar_config devuelve una lista plana nombre->valor", {
  cfg <- construir_config()
  valores <- aplanar_config(cfg)
  expect_type(valores, "list")
  expect_equal(valores$ifr_base, cfg$ifr_base$valor)
  expect_false(is.list(valores$ifr_base))
})
