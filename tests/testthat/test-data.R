test_that("validar_owid_diario detecta columnas faltantes", {
  datos_incompletos <- tibble::tibble(iso_code = "ARG", date = as.Date("2021-01-01"), new_cases = 1)
  expect_error(validar_owid_diario(datos_incompletos), "Faltan columnas")
})

test_that("validar_owid_diario detecta tipos incorrectos", {
  datos <- tibble::tibble(
    iso_code = "ARG", date = "2021-01-01",  # character, no Date: debe fallar
    new_cases = 1, new_deaths = 1, new_tests = 1,
    positive_rate = 0.1, people_fully_vaccinated = 0, population = 1e6
  )
  expect_error(validar_owid_diario(datos), "date")
})

test_that("validar_owid_diario reporta (sin frenar) NA y negativos como advertencias", {
  datos <- tibble::tibble(
    iso_code = "ARG", date = as.Date("2021-01-01") + 0:4,
    new_cases = c(1, NA, -3, 4, 5), new_deaths = c(0, 0, 0, 0, 0),
    new_tests = c(10, 10, 10, 10, 10), positive_rate = 0.1,
    people_fully_vaccinated = 0, population = 1e6
  )
  resultado <- validar_owid_diario(datos)
  expect_false(resultado$ok)
  expect_true(any(grepl("NA", resultado$advertencias)))
  expect_true(any(grepl("negativo", resultado$advertencias)))
})

test_that("validar_owid_diario frena ante fechas duplicadas", {
  datos <- tibble::tibble(
    iso_code = "ARG", date = rep(as.Date("2021-01-01"), 2),
    new_cases = c(1, 1), new_deaths = c(0, 0), new_tests = c(1, 1),
    positive_rate = 0.1, people_fully_vaccinated = 0, population = 1e6
  )
  expect_error(validar_owid_diario(datos), "duplicadas")
})

test_that("agregar_owid_semanal descarta semanas incompletas por defecto", {
  # 10 dias: una semana completa (7 dias) + 3 dias sueltos (semana incompleta).
  datos_diario <- tibble::tibble(
    iso_code = "ARG", date = as.Date("2021-01-03") + 0:9,
    new_cases = 1, new_deaths = 1, new_tests = 10,
    positive_rate = 0.1, people_fully_vaccinated = 0, population = 1e6
  )
  resultado <- agregar_owid_semanal(datos_diario, as.Date("2021-01-03"), as.Date("2021-01-12"))

  expect_equal(nrow(resultado$datos), 1)  # solo la semana completa
  expect_true(any(grepl("menos de 7 dias", resultado$advertencias)))
})

test_that("agregar_owid_semanal conserva semanas incompletas si se pide explicitamente", {
  datos_diario <- tibble::tibble(
    iso_code = "ARG", date = as.Date("2021-01-03") + 0:9,
    new_cases = 1, new_deaths = 1, new_tests = 10,
    positive_rate = 0.1, people_fully_vaccinated = 0, population = 1e6
  )
  resultado <- agregar_owid_semanal(datos_diario, as.Date("2021-01-03"), as.Date("2021-01-12"), incluir_semanas_incompletas = TRUE)

  expect_equal(nrow(resultado$datos), 2)
})

test_that("agregar_owid_semanal lleva a 0 (con advertencia) un total semanal negativo", {
  datos_diario <- tibble::tibble(
    iso_code = "ARG", date = as.Date("2021-01-03") + 0:6,
    new_cases = 1, new_deaths = c(-100, 1, 1, 1, 1, 1, 1),  # correccion retroactiva grande
    new_tests = 10, positive_rate = 0.1, people_fully_vaccinated = 0, population = 1e6
  )
  resultado <- agregar_owid_semanal(datos_diario, as.Date("2021-01-03"), as.Date("2021-01-09"))

  expect_equal(resultado$datos$muertes_owid[1], 0)
  expect_true(any(grepl("negativo", resultado$advertencias)))
})

test_that("calcular_sigma_desde_owid reproduce la cobertura final por construccion", {
  n_poblacion <- 1e6
  fecha_inicio <- as.Date("2021-01-03")
  fecha_vacunacion <- as.Date("2021-02-07")  # semana 6
  datos_semanal <- fixture_datos_semanal(n_poblacion = n_poblacion)

  sigma <- calcular_sigma_desde_owid(datos_semanal, n_poblacion, fecha_vacunacion, fecha_inicio)
  t_semanas <- max(datos_semanal$semana) - semana_de_fecha(fecha_vacunacion, fecha_inicio)
  cobertura_reconstruida <- 1 - exp(-sigma * t_semanas)
  cobertura_real <- tail(datos_semanal$vacunados_completos_owid, 1) / n_poblacion

  expect_equal(cobertura_reconstruida, cobertura_real, tolerance = 1e-9)
})

test_that("hash_archivo es estable ante un roundtrip write_csv/read_csv (regresion del bug de cache)", {
  datos <- tibble::tibble(a = 1:5, b = letters[1:5])
  ruta <- tempfile(fileext = ".csv")
  on.exit(unlink(ruta))
  readr::write_csv(datos, ruta)
  hash1 <- hash_archivo(ruta)
  # Releer y reescribir no deberia cambiar el archivo si los datos son los mismos.
  datos_releidos <- readr::read_csv(ruta, show_col_types = FALSE)
  readr::write_csv(datos_releidos, ruta)
  hash2 <- hash_archivo(ruta)
  expect_equal(hash1, hash2)
})
