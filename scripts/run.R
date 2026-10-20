# ==============================================================================
# RUN - Orquestacion end-to-end del modelo SIVQRD - Argentina COVID-19
# ==============================================================================
# Pipeline completo: datos -> tasas implicitas -> modelo -> calibracion ->
# resultados, y guarda todo en outputs/ (CSV, RDS, PNG, sessionInfo). Correr:
#   Rscript scripts/run.R
# desde la raiz del repo (ver README.md para dependencias y como correr tests).
#
# CALIBRACION (especificacion del usuario): la UNICA variable de ajuste es beta
# (un valor por segmento). theta, IFR, gamma, delta y ve_contagio estan FIJOS
# (literatura); sigma se calcula de los datos. El metodo PRINCIPAL ajusta los 5
# beta de forma GLOBAL y SIMULTANEA minimizando el RMSE de log(1+muertes).
#
# Los pasos de diagnostico mas costosos estan detras de flags, apagables para
# una corrida rapida.
# ==============================================================================

# NOTA DE COSTO: con la serie 2020-2023 (~199 semanas) cada calibracion global
# tarda ~1-3 min (el ODE corre sobre toda la serie en cada evaluacion de
# optim()). Por eso los diagnosticos que corren VARIAS calibraciones estan
# apagados por defecto. Prenderlos para el analisis completo (ver tiempos).
opciones <- list(
  forzar_descarga_owid     = FALSE,  # TRUE para ignorar la cache y re-descargar OWID
  correr_secuencial        = TRUE,   # ajuste secuencial: comparacion + warm start (~10-20s)
  correr_multistart        = FALSE,  # multi-start del ajuste global, robustez (~5 calibraciones)
  correr_sensibilidad_theta = TRUE,  # barrido del theta FIJO, clave porque theta es un supuesto (~4 calibraciones)
  correr_sensibilidad_oat  = FALSE,  # OAT sobre parametros de literatura (~11 calibraciones, lento)
  correr_bootstrap         = FALSE,  # bootstrap de los beta (1 calibracion global por rep, muy lento)
  n_rep_bootstrap          = 20
)

dir_salida <- "outputs"
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)

suppressMessages({
  library(tidyverse)
  library(scales)
})

for (archivo in c("config.R", "data.R", "model.R", "rates.R", "simulate.R", "metrics.R", "calibrate.R")) {
  source(file.path("R", archivo))
}

# ------------------------------------------------------------------------------
# 1) Config
# ------------------------------------------------------------------------------
cfg <- construir_config()
validar_config(cfg, exigir_sigma = FALSE)  # sigma se calcula y se valida mas abajo
valores <- aplanar_config(cfg)

# ------------------------------------------------------------------------------
# 2) Datos (OWID: cache si existe, si no descarga; validados y agregados)
# ------------------------------------------------------------------------------
owid <- obtener_datos_owid_diario(forzar_descarga = opciones$forzar_descarga_owid)
cat(sprintf("Datos OWID: origen=%s, descargados=%s, hash=%s\n",
            owid$metadata$origen, owid$metadata$descargado_en, substr(owid$metadata$hash_sha256, 1, 12)))

validacion <- validar_owid_diario(owid$datos %>% filter(date >= valores$fecha_inicio_serie, date <= valores$fecha_fin_serie))
if (length(validacion$advertencias) > 0) {
  cat("Advertencias de validacion de datos diarios:\n"); cat(paste(" -", validacion$advertencias, collapse = "\n"), "\n")
}

agregacion <- agregar_owid_semanal(owid$datos, valores$fecha_inicio_serie, valores$fecha_fin_serie)
datos_semanal <- agregacion$datos
if (length(agregacion$advertencias) > 0) {
  cat("Advertencias de agregacion semanal:\n"); cat(paste(" -", agregacion$advertencias, collapse = "\n"), "\n")
}
write_csv(datos_semanal, file.path(dir_salida, "owid_semanal.csv"))

n_poblacion <- mean(datos_semanal$poblacion_owid)

# ------------------------------------------------------------------------------
# 3) sigma (depende de los datos), validacion final de config
# ------------------------------------------------------------------------------
valores$sigma <- calcular_sigma_desde_owid(datos_semanal, n_poblacion, valores$fecha_vacunacion, valores$fecha_inicio_serie)
cfg$sigma$valor <- valores$sigma
validar_config(cfg, exigir_sigma = TRUE)
cat(sprintf("sigma calculado desde OWID: %.5f (1/semana)\n", valores$sigma))

# ------------------------------------------------------------------------------
# 4) Parametros biologicos, tramos, estado inicial
# ------------------------------------------------------------------------------
biologicos <- derivar_parametros_biologicos(valores)

semana_inicio_aspo       <- semana_de_fecha(valores$fecha_aspo, valores$fecha_inicio_serie)
semana_inicio_dispo       <- semana_de_fecha(valores$fecha_dispo, valores$fecha_inicio_serie)
semana_inicio_vacunacion  <- semana_de_fecha(valores$fecha_vacunacion, valores$fecha_inicio_serie)
semana_inicio_omicron     <- semana_de_fecha(valores$fecha_omicron, valores$fecha_inicio_serie)
cortes_semana <- c(1, semana_inicio_aspo, semana_inicio_dispo, semana_inicio_vacunacion, semana_inicio_omicron, max(datos_semanal$semana) + 1)

tramos_def <- tibble(
  tramo = 1:5, nombre = c("Pre-intervencion", "ASPO", "DISPO", "Vacunacion", "Omicron"),
  semana_ini = cortes_semana[1:5], semana_fin = cortes_semana[2:6] - 1
)

estado_inicial <- construir_estado_inicial(valores$s0_frac, valores$i0_personas, n_poblacion)

# ------------------------------------------------------------------------------
# 5) Tasas IMPLICITAS por semana (reconstruccion tipo Markov-Excel)
# ------------------------------------------------------------------------------
cat("\nCalculando tasas implicitas por semana (estructura tipo Markov)...\n")
tasas_impl <- construir_tasas_implicitas_semanales(
  datos_semanal, n_poblacion,
  dias_infecciosos = valores$dias_infecciosos, ifr_base = valores$ifr_base,
  i0_personas = valores$i0_personas
)
if (length(tasas_impl$advertencias) > 0) {
  cat("Advertencias de tasas implicitas:\n"); cat(paste(" -", tasas_impl$advertencias, collapse = "\n"), "\n")
}
write_csv(tasas_impl$tasas, file.path(dir_salida, "tasas_implicitas_semanales.csv"))
tasas_por_tramo <- resumir_tasas_por_tramo(tasas_impl$tasas, tramos_def)
write_csv(tasas_por_tramo, file.path(dir_salida, "tasas_implicitas_por_tramo.csv"))
print(tasas_por_tramo %>% select(tramo, nombre, beta_implicita_prom, re_implicito_prom))

# ------------------------------------------------------------------------------
# 6) Calibracion
# ------------------------------------------------------------------------------
# 6a) Secuencial (comparacion + warm start barato para el global)
betas_warm <- NULL
if (opciones$correr_secuencial) {
  cat("\nCalibrando (secuencial por segmento, para comparacion/warm start)...\n")
  secuencial <- calibrar_secuencial(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    theta = valores$theta, sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    beta_min = valores$beta_min, beta_max = valores$beta_max
  )
  betas_warm <- secuencial$betas
  cat(sprintf("  betas (secuencial): %s | RMSE-log global: %.4f\n",
              paste(round(secuencial$betas, 3), collapse = ", "), secuencial$mejor_error))
}

# 6b) GLOBAL simultaneo (metodo PRINCIPAL)
cat("\nCalibrando (GLOBAL simultaneo, metodo principal: solo beta, theta fijo)...\n")
calibracion <- calibrar_beta_global(
  estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
  theta = valores$theta, sigma = valores$sigma, ve_contagio = valores$ve_contagio,
  semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
  datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
  beta_min = valores$beta_min, beta_max = valores$beta_max,
  betas_iniciales = betas_warm, maxit = valores$maxit_optimizador
)
cat(sprintf("  betas (global): %s\n", paste(round(calibracion$betas, 3), collapse = ", ")))
cat(sprintf("  theta (fijo): %.3f | RMSE-log global: %.4f | convergencia optim(): %d | evaluaciones: %d\n",
            calibracion$theta, calibracion$mejor_error, calibracion$convergencia, calibracion$n_evaluaciones))

# 6c) Multi-start (robustez del optimo al punto de partida)
if (opciones$correr_multistart) {
  cat("\nMulti-start del ajuste global...\n")
  ms <- calibrar_beta_global_multistart(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    theta = valores$theta, sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    beta_min = valores$beta_min, beta_max = valores$beta_max, betas_centrales = betas_warm,
    n_aleatorios = 4, maxit = valores$maxit_optimizador
  )
  write_csv(ms$tabla %>% select(-betas), file.path(dir_salida, "multistart_global.csv"))
  cat(sprintf("  rango de RMSE-log entre starts: [%.5f, %.5f]\n", min(ms$tabla$error), max(ms$tabla$error)))
  if (ms$mejor$mejor_error < calibracion$mejor_error - 1e-8) {
    cat("  (un arranque alternativo mejoro el optimo; se adopta el mejor del multi-start)\n")
    calibracion <- ms$mejor
  }
}

# ------------------------------------------------------------------------------
# 7) Resultados por tramo (beta calibrado vs implicito, R0, Re, tasas)
# ------------------------------------------------------------------------------
tasas_mort <- calcular_tasa_mortalidad_por_tramo(calibracion$emisiones, n_poblacion, tramos_def)
tasas_cont <- calcular_tasa_contagio_por_tramo(calibracion$incidencia, tramos_def)
promedios_sv <- calibracion$simulacion %>% group_by(tramo) %>% summarise(s_prom = mean(s), v_prom = mean(v), .groups = "drop")

resumen_tramos <- tramos_def %>%
  mutate(beta = calibracion$betas, theta_vigente = if_else(nombre == "Pre-intervencion", 0, calibracion$theta)) %>%
  left_join(promedios_sv, by = "tramo") %>%
  left_join(tasas_mort, by = "tramo") %>%
  left_join(tasas_cont, by = "tramo") %>%
  left_join(tasas_por_tramo %>% select(tramo, beta_implicita_prom), by = "tramo") %>%
  mutate(
    r0 = calcular_r0(beta, theta_vigente, biologicos$gamma, biologicos$delta),
    re = calcular_re(beta, theta_vigente, s_prom, v_prom, biologicos$gamma, biologicos$delta,
                      valores$ve_contagio, biologicos$gamma_v, biologicos$delta_v)
  ) %>%
  select(tramo, nombre, semana_ini, semana_fin, beta, beta_implicita_prom, r0, re,
         tasa_contagio_semanal, tasa_mortalidad_semanal)

print(resumen_tramos)
write_csv(resumen_tramos, file.path(dir_salida, "resumen_tramos.csv"))

metricas_globales <- tibble(
  metodo = calibracion$metodo,
  rmse_log  = calibracion$mejor_error,
  msre_log  = calcular_msre_log(calibracion$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas),
  mse_nivel = calcular_mse_nivel(calibracion$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas),
  mse_log   = calcular_mse_log(calibracion$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas)
)
if (opciones$correr_secuencial) {
  metricas_globales <- bind_rows(metricas_globales, tibble(
    metodo = secuencial$metodo, rmse_log = secuencial$mejor_error,
    msre_log  = calcular_msre_log(secuencial$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas),
    mse_nivel = calcular_mse_nivel(secuencial$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas),
    mse_log   = calcular_mse_log(secuencial$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas)
  ))
}
write_csv(metricas_globales, file.path(dir_salida, "metricas_globales.csv"))

residuos <- calcular_residuos(calibracion$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas, tramos_def)
write_csv(residuos, file.path(dir_salida, "residuos_calibracion.csv"))

verificacion_vacunacion <- verificar_significado_vacunacion(datos_semanal, calibracion$simulacion, n_poblacion)
write_csv(verificacion_vacunacion, file.path(dir_salida, "verificacion_cobertura_vacunacion.csv"))

# ------------------------------------------------------------------------------
# 8) Sensibilidad al theta FIJO (theta es un supuesto, no se calibra)
# ------------------------------------------------------------------------------
if (opciones$correr_sensibilidad_theta) {
  cat("\nBarrido de sensibilidad sobre el theta FIJO...\n")
  sens_theta <- analisis_sensibilidad_theta(
    thetas = c(0, 0.15, 0.3, 0.6),
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    beta_min = valores$beta_min, beta_max = valores$beta_max, betas_iniciales = calibracion$betas, maxit = 250
  )
  write_csv(sens_theta, file.path(dir_salida, "sensibilidad_theta.csv"))
  print(sens_theta %>% select(theta, rmse_log))
}

# ------------------------------------------------------------------------------
# 9) Sensibilidad OAT a parametros de literatura -- OPCIONAL, costosa
# ------------------------------------------------------------------------------
if (opciones$correr_sensibilidad_oat) {
  cat("\nCorriendo analisis de sensibilidad OAT...\n")
  parametros_sensibilidad <- list(
    dias_infecciosos = c(8, 14),
    ifr_base         = c(0.004, 0.009),
    rezago_semanas   = c(2, 4),
    ve_contagio      = c(0.44, 0.68),
    theta            = c(0.1, 0.5)
  )
  sensibilidad <- map_dfr(names(parametros_sensibilidad), function(p) {
    analisis_sensibilidad(
      cfg_valores = valores, parametro = p, valores_alternativos = parametros_sensibilidad[[p]],
      tramos_def = tramos_def, datos_semanal = datos_semanal, n_poblacion = n_poblacion,
      semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
      beta_min = valores$beta_min, beta_max = valores$beta_max
    )
  })
  write_csv(sensibilidad, file.path(dir_salida, "sensibilidad_oat.csv"))
  print(sensibilidad %>% select(parametro, valor, rmse_log))
}

# ------------------------------------------------------------------------------
# 10) Bootstrap de los beta -- OPCIONAL (caro)
# ------------------------------------------------------------------------------
if (opciones$correr_bootstrap) {
  cat("\nCorriendo bootstrap de los beta...\n")
  boot_betas <- bootstrap_betas(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    theta = valores$theta, sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    betas_punto = calibracion$betas, beta_min = valores$beta_min, beta_max = valores$beta_max,
    n_rep = opciones$n_rep_bootstrap
  )
  ic_betas <- as_tibble(t(apply(boot_betas, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)),
                        rownames = "beta") %>% rename(p2.5 = `2.5%`, p50 = `50%`, p97.5 = `97.5%`)
  write_csv(ic_betas, file.path(dir_salida, "bootstrap_betas_ic.csv"))
  print(ic_betas)
}

# ------------------------------------------------------------------------------
# 11) Graficos principales
# ------------------------------------------------------------------------------
fechas_quiebre <- c(valores$fecha_aspo, valores$fecha_dispo, valores$fecha_vacunacion, valores$fecha_omicron)

grafico_compartimentos <- calibracion$simulacion %>%
  mutate(fecha = valores$fecha_inicio_serie + (semana - 1) * 7) %>%
  pivot_longer(cols = c(s, i, iv, v, q, r, d), names_to = "compartimento", values_to = "fraccion") %>%
  mutate(compartimento = factor(compartimento,
                                  levels = c("s", "i", "iv", "v", "q", "r", "d"),
                                  labels = c("Susceptibles", "Infectados", "Infectados vacunados",
                                              "Vacunados", "Cuarentena", "Recuperados", "Fallecidos"))) %>%
  ggplot(aes(x = fecha, y = fraccion, color = compartimento)) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(title = "Modelo SIVQRD - Argentina", subtitle = "Evolucion semanal de los compartimentos",
       x = NULL, y = "Fraccion de la poblacion", color = NULL) +
  theme_minimal()
ggsave(file.path(dir_salida, "grafico_compartimentos.png"), grafico_compartimentos, width = 8, height = 5, dpi = 120)

comparacion_muertes <- alinear_muertes_con_rezago(calibracion$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas) %>%
  mutate(fecha = valores$fecha_inicio_serie + (semana - 1) * 7) %>%
  pivot_longer(cols = c(muertes_modelo, muertes_owid), names_to = "serie", values_to = "muertes") %>%
  mutate(serie = recode(serie, muertes_modelo = "Modelo SIVQRD", muertes_owid = "OWID (observado)"))

grafico_muertes <- ggplot(comparacion_muertes, aes(x = fecha, y = muertes, color = serie)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = fechas_quiebre, linetype = "dashed", color = "grey50") +
  scale_x_date(date_labels = "%b-%Y") +
  labs(title = "Muertes semanales: modelo vs. observado",
       subtitle = "Lineas punteadas: cortes ASPO / DISPO / vacunacion / Omicron",
       x = NULL, y = "Muertes por semana", color = NULL) +
  theme_minimal()
ggsave(file.path(dir_salida, "grafico_muertes.png"), grafico_muertes, width = 8, height = 5, dpi = 120)

# beta implicito (semanal) vs beta calibrado (escalon por tramo)
escalon_beta <- resumen_tramos %>%
  transmute(semana_ini, semana_fin, beta) %>%
  rowwise() %>% do(tibble(semana = seq(.$semana_ini, .$semana_fin), beta = .$beta)) %>% ungroup() %>%
  mutate(fecha = valores$fecha_inicio_serie + (semana - 1) * 7)

grafico_beta <- ggplot() +
  geom_line(data = tasas_impl$tasas, aes(x = fecha, y = beta_implicita, color = "beta implicito (semanal)"), linewidth = 0.5, na.rm = TRUE) +
  geom_line(data = escalon_beta, aes(x = fecha, y = beta, color = "beta calibrado (por tramo)"), linewidth = 1) +
  geom_vline(xintercept = fechas_quiebre, linetype = "dashed", color = "grey60") +
  scale_x_date(date_labels = "%b-%Y") +
  coord_cartesian(ylim = c(0, valores$beta_max)) +
  labs(title = "beta: implicito de los datos vs. calibrado por tramo",
       x = NULL, y = "beta (1/semana)", color = NULL) +
  theme_minimal()
ggsave(file.path(dir_salida, "grafico_beta_implicito_vs_calibrado.png"), grafico_beta, width = 8, height = 5, dpi = 120)

grafico_vacunacion <- verificacion_vacunacion %>%
  mutate(fecha = valores$fecha_inicio_serie + (semana - 1) * 7) %>%
  pivot_longer(cols = c(cobertura_observada, cobertura_modelada), names_to = "serie", values_to = "cobertura") %>%
  mutate(serie = recode(serie, cobertura_observada = "OWID (observado)", cobertura_modelada = "Modelo (sigma constante)")) %>%
  ggplot(aes(x = fecha, y = cobertura, color = serie)) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(title = "Cobertura de vacunacion: observada vs. modelada",
       subtitle = "sigma es un hazard CONSTANTE calibrado solo al nivel final; no reproduce el ritmo real semana a semana",
       x = NULL, y = "Cobertura (vacunados completos / poblacion)", color = NULL) +
  theme_minimal()
ggsave(file.path(dir_salida, "grafico_cobertura_vacunacion.png"), grafico_vacunacion, width = 8, height = 5, dpi = 120)

# ------------------------------------------------------------------------------
# 12) Objetos completos + sessionInfo, para reproducibilidad
# ------------------------------------------------------------------------------
saveRDS(calibracion, file.path(dir_salida, "resultado_calibracion_global.rds"))
if (opciones$correr_secuencial) saveRDS(secuencial, file.path(dir_salida, "resultado_calibracion_secuencial.rds"))
saveRDS(tasas_impl, file.path(dir_salida, "tasas_implicitas.rds"))
saveRDS(list(cfg = cfg, valores = valores, tramos_def = tramos_def, n_poblacion = n_poblacion,
             datos_semanal = datos_semanal, owid_metadata = owid$metadata),
        file.path(dir_salida, "contexto_corrida.rds"))

writeLines(capture.output(sessionInfo()), file.path(dir_salida, "sessionInfo.txt"))

cat(sprintf("\nListo. Resultados en '%s/'.\n", dir_salida))
