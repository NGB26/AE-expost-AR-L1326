# Punto de entrada de testthat. Correr desde la raiz del repo con:
#   Rscript -e 'testthat::test_dir("tests/testthat")'
# o, mas simple, con el wrapper:
#   Rscript tests/testthat.R

suppressMessages(library(testthat))
suppressMessages(library(dplyr))

directorio_raiz <- if (basename(getwd()) == "tests") ".." else "."
for (archivo in c("config.R", "data.R", "model.R", "rates.R", "simulate.R", "metrics.R", "calibrate.R")) {
  source(file.path(directorio_raiz, "R", archivo))
}

test_dir(file.path(directorio_raiz, "tests", "testthat"), reporter = "summary")
