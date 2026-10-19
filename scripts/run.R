# ==============================================================================
# RUN - Orquestacion end-to-end del modelo SIVQRD - Argentina COVID-19
# ==============================================================================
# Corre el pipeline completo: datos -> modelo -> calibracion -> resultados, y
# guarda todo en outputs/ (CSV, RDS, PNG, sessionInfo). Para correrlo:
#   Rscript scripts/run.R
# desde la raiz del repo (ver README.md para la guia completa, dependencias
# y como correr los tests antes de esto).
#
# Los pasos de diagnostico mas costosos en tiempo de computo (calibracion
# conjunta global: ~3-4 minutos; sensibilidad con varios parametros;
# bootstrap con muchas repeticiones) estan detras de flags mas abajo,
# apagados por defecto para que una corrida de rutina sea rapida (~15-30s).
# Prenderlos para el analisis completo que pide DECISIONS.md.
# ==============================================================================

opciones <- list(
  forzar_descarga_owid   = FALSE,  # TRUE para ignorar la cache y re-descargar OWID
  correr_calibracion_conjunta = TRUE,   # comparacion vs. el baseline secuencial (~3-4 min)
  correr_perfil_identificabilidad = TRUE,  # grilla beta-theta de un tramo (rapido)
  correr_multistart      = TRUE,   # varios theta_inicial (rapido, ~15-30s)
  correr_sensibilidad    = TRUE,   # OAT sobre varios parametros (~1-2 min)
  correr_bootstrap       = TRUE,   # bootstrap de theta (~20-60s con n_rep bajo)
  n_rep_bootstrap        = 20
)

dir_salida <- "outputs"
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)

suppressMessages({
  library(tidyverse)
  library(scales)
})

for (archivo in c("config.R", "data.R", "model.R", "simulate.R", "metrics.R", "calibrate.R")) {
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
# 5) Calibracion: baseline secuencial-alternante
# ------------------------------------------------------------------------------
cat("\nCalibrando (secuencial-alternante)...\n")
baseline <- calibrar_secuencial_alternante(
  estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
  sigma = valores$sigma, ve_contagio = valores$ve_contagio,
  semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
  datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
  theta_inicial = valores$theta_inicial, theta_min = valores$theta_min, theta_max = valores$theta_max,
  beta_min = valores$beta_min, beta_max = valores$beta_max,
  max_iteraciones = valores$max_iteraciones_externas, tolerancia = valores$tolerancia_convergencia
)
cat(sprintf("  betas: %s\n", paste(round(baseline$betas, 3), collapse = ", ")))
cat(sprintf("  theta: %.4f | error (MSRE-log): %.4f\n", baseline$theta, baseline$mejor_error))

# ------------------------------------------------------------------------------
# 6) Resultados del baseline: tabla por tramo, metricas diagnosticas
# ------------------------------------------------------------------------------
tasas_mort <- calcular_tasa_mortalidad_por_tramo(baseline$emisiones, n_poblacion, tramos_def)
tasas_cont <- calcular_tasa_contagio_por_tramo(baseline$incidencia, tramos_def)
promedios_sv <- baseline$simulacion %>% group_by(tramo) %>% summarise(s_prom = mean(s), v_prom = mean(v), .groups = "drop")

resumen_tramos <- tramos_def %>%
  mutate(beta = baseline$betas, theta_vigente = if_else(nombre == "Pre-intervencion", 0, baseline$theta)) %>%
  left_join(promedios_sv, by = "tramo") %>%
  left_join(tasas_mort, by = "tramo") %>%
  left_join(tasas_cont, by = "tramo") %>%
  mutate(
    r0 = calcular_r0(beta, theta_vigente, biologicos$gamma, biologicos$delta),
    re = calcular_re(beta, theta_vigente, s_prom, v_prom, biologicos$gamma, biologicos$delta,
                      valores$ve_contagio, biologicos$gamma_v, biologicos$delta_v)
  ) %>%
  select(tramo, nombre, semana_ini, semana_fin, beta, r0, re, tasa_contagio_semanal, tasa_mortalidad_semanal)

print(resumen_tramos)
write_csv(resumen_tramos, file.path(dir_salida, "resumen_tramos.csv"))

metricas_globales <- tibble(
  metodo = "secuencial_alternante",
  msre_log = baseline$mejor_error,
  mse_nivel = calcular_mse_nivel(baseline$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas),
  mse_log  = calcular_mse_log(baseline$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas),
  n_evaluaciones_ode_aprox = baseline$n_evaluaciones_ode
)

residuos <- calcular_residuos(baseline$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas, tramos_def)
write_csv(residuos, file.path(dir_salida, "residuos_baseline.csv"))

verificacion_vacunacion <- verificar_significado_vacunacion(datos_semanal, baseline$simulacion, n_poblacion)
write_csv(verificacion_vacunacion, file.path(dir_salida, "verificacion_cobertura_vacunacion.csv"))

# ------------------------------------------------------------------------------
# 7) Calibracion conjunta global (comparacion) -- OPCIONAL, costosa
# ------------------------------------------------------------------------------
if (opciones$correr_calibracion_conjunta) {
  cat("\nCalibrando (conjunto global, puede tardar unos minutos)...\n")
  conjunto <- calibrar_conjunto_global(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    beta_min = valores$beta_min, beta_max = valores$beta_max, theta_min = valores$theta_min, theta_max = valores$theta_max,
    valores_iniciales = c(baseline$betas, baseline$theta)
  )
  cat(sprintf("  betas: %s\n", paste(round(conjunto$betas, 3), collapse = ", ")))
  cat(sprintf("  theta: %.4f | error (MSRE-log): %.4f | convergencia optim(): %d | evaluaciones: %d\n",
              conjunto$theta, conjunto$mejor_error, conjunto$convergencia, conjunto$n_evaluaciones))

  metricas_globales <- bind_rows(metricas_globales, tibble(
    metodo = "conjunto_global", msre_log = conjunto$mejor_error,
    mse_nivel = calcular_mse_nivel(conjunto$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas),
    mse_log  = calcular_mse_log(conjunto$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas),
    n_evaluaciones_ode_aprox = conjunto$n_evaluaciones * nrow(tramos_def)
  ))

  comparacion_metodos <- tibble(
    tramo = tramos_def$tramo, nombre = tramos_def$nombre,
    beta_secuencial = baseline$betas, beta_conjunto = conjunto$betas
  )
  write_csv(comparacion_metodos, file.path(dir_salida, "comparacion_metodos_beta.csv"))
  saveRDS(conjunto, file.path(dir_salida, "resultado_calibracion_conjunta.rds"))
}

write_csv(metricas_globales, file.path(dir_salida, "metricas_globales.csv"))

# ------------------------------------------------------------------------------
# 8) Identificabilidad: perfil de perdida en un tramo + multistart
# ------------------------------------------------------------------------------
if (opciones$correr_perfil_identificabilidad) {
  cat("\nGenerando perfil de perdida beta-theta (tramo DISPO)...\n")
  perfil <- generar_perfil_identificabilidad(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    tramo_objetivo = 3, betas_referencia = baseline$betas,
    rango_beta = seq(valores$beta_min, valores$beta_max, length.out = 12),
    rango_theta = seq(valores$theta_min, valores$theta_max, length.out = 12)
  )
  write_csv(perfil, file.path(dir_salida, "perfil_identificabilidad_dispo.csv"))

  grafico_perfil <- ggplot(perfil, aes(x = beta, y = theta, fill = error_global)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", direction = -1) +
    labs(title = "Perfil de perdida beta-theta (tramo DISPO)",
         subtitle = "Una cresta diagonal indica no-identificabilidad entre beta y theta",
         x = "beta (tramo DISPO)", y = "theta", fill = "MSRE-log") +
    theme_minimal()
  ggsave(file.path(dir_salida, "perfil_identificabilidad_dispo.png"), grafico_perfil, width = 7, height = 5, dpi = 120)
}

if (opciones$correr_multistart) {
  cat("Corriendo multi-start (theta_inicial = 0, 0.2, 0.5, 1)...\n")
  multistart <- calibrar_multistart(
    thetas_iniciales = c(0, 0.2, 0.5, 1),
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    theta_min = valores$theta_min, theta_max = valores$theta_max, beta_min = valores$beta_min, beta_max = valores$beta_max,
    max_iteraciones = valores$max_iteraciones_externas, tolerancia = valores$tolerancia_convergencia
  )
  write_csv(multistart %>% select(-betas), file.path(dir_salida, "multistart_theta.csv"))
  cat(sprintf("  rango de theta_final entre starts: [%.4f, %.4f] (rango de error: [%.4f, %.4f])\n",
              min(multistart$theta_final), max(multistart$theta_final), min(multistart$mejor_error), max(multistart$mejor_error)))
}

# ------------------------------------------------------------------------------
# 9) Sensibilidad OAT -- OPCIONAL, costosa
# ------------------------------------------------------------------------------
if (opciones$correr_sensibilidad) {
  cat("\nCorriendo analisis de sensibilidad OAT...\n")
  parametros_sensibilidad <- list(
    dias_infecciosos = c(7, 14),
    ifr_base         = c(0.004, 0.009),
    rezago_semanas   = c(2, 3),
    ve_muerte        = c(0.6, 0.8),
    ve_contagio      = c(0.44, 0.68)
  )
  sensibilidad <- map_dfr(names(parametros_sensibilidad), function(p) {
    analisis_sensibilidad(
      cfg_valores = valores, parametro = p, valores_alternativos = parametros_sensibilidad[[p]],
      tramos_def = tramos_def, datos_semanal = datos_semanal, n_poblacion = n_poblacion,
      semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
      max_iteraciones = 3, tolerancia = valores$tolerancia_convergencia
    )
  })
  write_csv(sensibilidad, file.path(dir_salida, "sensibilidad_oat.csv"))
  print(sensibilidad %>% select(parametro, valor, mejor_error, theta))
}

# ------------------------------------------------------------------------------
# 10) Bootstrap de theta -- OPCIONAL
# ------------------------------------------------------------------------------
if (opciones$correr_bootstrap) {
  cat("\nCorriendo bootstrap de theta...\n")
  boot_theta <- bootstrap_theta(
    estado_inicial = estado_inicial, tramos_def = tramos_def, biologicos = biologicos,
    sigma = valores$sigma, ve_contagio = valores$ve_contagio,
    semana_inicio_aspo = semana_inicio_aspo, semana_inicio_vacunacion = semana_inicio_vacunacion,
    datos_semanal = datos_semanal, n_poblacion = n_poblacion, rezago_semanas = valores$rezago_semanas,
    betas_fijos = baseline$betas, theta_min = valores$theta_min, theta_max = valores$theta_max,
    n_rep = opciones$n_rep_bootstrap
  )
  write_csv(tibble(theta_bootstrap = boot_theta), file.path(dir_salida, "bootstrap_theta.csv"))
  cat(sprintf("  IC empirico 95%% (condicional a los beta del baseline): [%.4f, %.4f]\n",
              quantile(boot_theta, 0.025), quantile(boot_theta, 0.975)))
}

# ------------------------------------------------------------------------------
# 11) Graficos principales
# ------------------------------------------------------------------------------
grafico_compartimentos <- baseline$simulacion %>%
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

fechas_quiebre <- c(valores$fecha_aspo, valores$fecha_dispo, valores$fecha_vacunacion, valores$fecha_omicron)
comparacion_muertes <- alinear_muertes_con_rezago(baseline$emisiones, datos_semanal, n_poblacion, valores$rezago_semanas) %>%
  mutate(fecha = valores$fecha_inicio_serie + (semana - 1) * 7) %>%
  pivot_longer(cols = c(muertes_modelo, muertes_owid), names_to = "serie", values_to = "muertes") %>%
  mutate(serie = recode(serie, muertes_modelo = "Modelo SIVQRD", muertes_owid = "OWID (observado)"))

grafico_muertes <- ggplot(comparacion_muertes, aes(x = fecha, y = muertes, color = serie)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = fechas_quiebre, linetype = "dashed", color = "grey50") +
  scale_x_date(date_labels = "%b-%Y") +
  labs(title = "Muertes semanales: modelo vs. observado", x = NULL, y = "Muertes por semana", color = NULL) +
  theme_minimal()
ggsave(file.path(dir_salida, "grafico_muertes.png"), grafico_muertes, width = 8, height = 5, dpi = 120)

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
saveRDS(baseline, file.path(dir_salida, "resultado_calibracion_baseline.rds"))
saveRDS(list(cfg = cfg, valores = valores, tramos_def = tramos_def, n_poblacion = n_poblacion,
             datos_semanal = datos_semanal, owid_metadata = owid$metadata),
        file.path(dir_salida, "contexto_corrida.rds"))

writeLines(capture.output(sessionInfo()), file.path(dir_salida, "sessionInfo.txt"))

cat(sprintf("\nListo. Resultados en '%s/'.\n", dir_salida))
