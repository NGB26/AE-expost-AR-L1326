# AE expost AR-L1326
This project is to build the economic analysis of the project AR-L1326 (COVID-19 prototype). The AE will be an ex-post cost-benefit analysis  

## Modelo epidemiologico SIVQRD

`scripts/modelo_sivqrd_argentina.R` implementa en R el modelo SIVQRD descrito
en `Modelo_SIR_para_PCR.docx` (ecuaciones 23-29), calibrado con la serie
semanal de OWID. Ya no depende de `data/raw/SIRD_Argentina_PCR.xlsx` (se
mantiene en el repo solo como referencia historica): todos los parametros
biologicos y las fechas de corte de beta se declaran a mano, en bloques
`# EDITAR` claramente marcados dentro del propio script (secciones 2.1 y
2.2), salvo la poblacion (N), que se toma directamente de los datos de OWID.

`sigma` y `ve_contagio` quedan en `NA_real_` hasta completarlos con fuente
bibliografica. `theta` ya no se fija a mano: se calibra junto con beta en un
proceso iterativo (Bloque 3) que alterna entre ajustar los 5 beta (uno por
tramo) y ajustar theta (un unico valor para toda la serie), quedandose con
el mejor par encontrado.
