test_that("la poblacion total se conserva (s+i+iv+v+q+r+d = 1)", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)

  sim <- simular_tramo(estado_inicial, 1, 20, beta_valor = 0.8, parametros_fijos)
  suma <- rowSums(sim[, c("s", "i", "iv", "v", "q", "r", "d")])

  expect_true(all(abs(suma - 1) < 1e-6))
})

test_that("los compartimentos fisicos no son negativos", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  parametros_fijos <- fixture_parametros_fijos(0.5, valores, biologicos)

  sim <- simular_tramo(estado_inicial, 1, 20, beta_valor = 1.5, parametros_fijos)
  valores_fisicos <- as.matrix(sim[, c("s", "i", "iv", "v", "q", "r", "d")])

  expect_true(all(valores_fisicos >= -1e-8))
})

test_that("beta=0 no produce nuevas infecciones", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  # sigma=0 tambien, para aislar el efecto de beta (si no, s baja igual por vacunacion)
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)
  parametros_fijos$sigma <- 0

  sim <- simular_tramo(estado_inicial, 1, 15, beta_valor = 0, parametros_fijos)

  # Sin nuevas infecciones, inf_acum no se mueve del cero inicial.
  expect_true(all(abs(sim$inf_acum) < 1e-10))
  # s se mantiene constante (nadie se infecta ni se vacuna).
  expect_true(all(abs(sim$s - sim$s[1]) < 1e-10))
})

test_that("sigma=0 evita que nadie se vacune (v se mantiene en 0)", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)
  parametros_fijos$sigma <- 0

  sim <- simular_tramo(estado_inicial, 1, 20, beta_valor = 0.8, parametros_fijos)

  expect_true(all(abs(sim$v) < 1e-10))
  expect_true(all(abs(sim$iv) < 1e-10))  # sin v, tampoco puede haber infectados vacunados
})

test_that("theta=0 evita que nadie entre en cuarentena (q se mantiene en 0)", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  parametros_fijos <- fixture_parametros_fijos(0, valores, biologicos)
  # theta=0 en construir_parametros_fijos ya alcanza porque el gating por
  # semana_inicio_aspo no puede subir el theta vigente por encima del valor
  # base (if (t>=aspo) theta else 0): con theta=0, siempre es 0.

  sim <- simular_tramo(estado_inicial, 1, 20, beta_valor = 0.8, parametros_fijos)

  expect_true(all(abs(sim$q) < 1e-10))
})

test_that("delta=0 y delta_v=0 evitan toda muerte (d se mantiene en 0)", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)
  estado_inicial <- fixture_estado_inicial(valores)
  parametros_fijos <- fixture_parametros_fijos(0.3, valores, biologicos)
  parametros_fijos$delta <- 0
  parametros_fijos$delta_v <- 0

  sim <- simular_tramo(estado_inicial, 1, 20, beta_valor = 0.8, parametros_fijos)

  expect_true(all(abs(sim$d) < 1e-10))
})

test_that("el IFR se recupera en un proceso competitivo simple (sin reinfeccion)", {
  # Sin nuevas infecciones (beta=0) ni cuarentena (theta=0), los i0_personas
  # iniciales en I solo pueden salir por dos caminos: recuperarse (gamma) o
  # morir (delta). Dado tiempo suficiente, la fraccion que termina en D sobre
  # el total que salio de I (R+D) debe converger al IFR usado para derivar
  # gamma y delta -- es la verificacion por SIMULACION de la cuenta algebraica
  # documentada en R/model.R::derivar_parametros_biologicos().
  for (ifr_prueba in c(0.01, 0.05, 0.2)) {
    valores <- fixture_valores()
    valores$ifr_base <- ifr_prueba
    biologicos <- derivar_parametros_biologicos(valores)

    parametros_fijos <- construir_parametros_fijos(
      theta_valor = 0, sigma = 0, ve_contagio = valores$ve_contagio,
      gamma = biologicos$gamma, delta = biologicos$delta,
      gamma_v = biologicos$gamma_v, delta_v = biologicos$delta_v,
      semana_inicio_aspo = 999999, semana_inicio_vacunacion = 999999  # nunca se activan
    )

    estado_inicial <- c(s = 0, i = 1, iv = 0, v = 0, q = 0, r = 0, d = 0, inf_acum = 0)
    # 200 semanas: de sobra para que gamma_total ~ 0.7/semana vacie I.
    sim <- simular_tramo(estado_inicial, 1, 200, beta_valor = 0, parametros_fijos)
    fila_final <- sim[nrow(sim), ]

    expect_lt(fila_final$i, 1e-6)  # I debe haberse vaciado
    ifr_observado <- fila_final$d / (fila_final$d + fila_final$r)
    expect_equal(ifr_observado, ifr_prueba, tolerance = 1e-4)
  }
})

test_that("calcular_re() (forma cerrada) coincide con el autovalor de la NGM", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)

  for (s_val in c(1, 0.6, 0.2)) {
    for (v_val in c(0, 0.3, 0.7)) {
      re_cerrada <- calcular_re(0.8, 0.2, s_val, v_val, biologicos$gamma, biologicos$delta,
                                 valores$ve_contagio, biologicos$gamma_v, biologicos$delta_v)
      re_numerica <- calcular_re_ngm_numerico(0.8, 0.2, s_val, v_val, biologicos$gamma, biologicos$delta,
                                               valores$ve_contagio, biologicos$gamma_v, biologicos$delta_v)
      expect_equal(re_cerrada, re_numerica, tolerance = 1e-10)
    }
  }
})

test_that("R0 es el caso particular de Re con s=1, v=0", {
  valores <- fixture_valores()
  biologicos <- fixture_biologicos(valores)

  r0 <- calcular_r0(0.8, 0.2, biologicos$gamma, biologicos$delta)
  re_en_dfe <- calcular_re(0.8, 0.2, s_val = 1, v_val = 0, biologicos$gamma, biologicos$delta,
                            valores$ve_contagio, biologicos$gamma_v, biologicos$delta_v)

  expect_equal(r0, re_en_dfe, tolerance = 1e-12)
})

test_that("derivar_parametros_biologicos cumple gamma+delta = 7/dias_infecciosos", {
  valores <- fixture_valores()
  biologicos <- derivar_parametros_biologicos(valores)

  expect_equal(biologicos$gamma + biologicos$delta, 7 / valores$dias_infecciosos, tolerance = 1e-12)
  expect_equal(biologicos$gamma_v + biologicos$delta_v, 7 / valores$dias_infecciosos_vac, tolerance = 1e-12)
  expect_equal(biologicos$delta / (biologicos$gamma + biologicos$delta), valores$ifr_base, tolerance = 1e-12)
})
