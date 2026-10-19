# AE expost AR-L1326
This project is to build the economic analysis of the project AR-L1326 (COVID-19 prototype). The AE will be an ex-post cost-benefit analysis  

## Modelo epidemiologico SIVQRD

`scripts/modelo_sivqrd_argentina.R` implementa en R el modelo SIVQRD descrito
en `Modelo_SIR_para_PCR.docx` (ecuaciones 23-29), calibrado con series
semanales de OWID y datos.gob.ar. Parte de `data/raw/SIRD_Argentina_PCR.xlsx`
para los parametros biologicos y las fechas de los tramos de politica
(hojas `Parametros` y `Tramos_Beta`); ese Excel sigue siendo el lugar donde
se editan esos supuestos.

Antes de correrlo hay que completar en el script, con fuente bibliografica,
los parametros que no estaban en el Excel original (`theta`, `sigma`,
`ve_contagio`, y opcionalmente `dias_infecciosos_vac`): estan marcados como
`NA_real_` con un comentario `# EDITAR`.
