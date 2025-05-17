clear all
set more off

import excel "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/IGPA_Trabajar.xlsx", ///
    sheet("Sheet 1") firstrow clear

* Identificar el primer año en que entra == 1
egen first_entry = min(year / (entra == 1)), by(empresa)
gen ever_in_ipsa = (entra == 1 | mantiene == 1)
egen suma_in_ipsa = total(ever_in_ipsa), by(empresa)
gen never_treated = (suma_in_ipsa == 0)


* Rellenar para la variable first_entry
bysort empresa (year): replace first_entry = first_entry[_n-1] if missing(first_entry)

save "data_expanded.dta", replace

* Tiempo del evento
gen event_time = year - first_entry
keep if inrange(event_time, -10, 10)

gen uno = 1
collapse (sum) empresas = uno, by(event_time)

twoway bar empresas event_time, ///
    barwidth(0.8) ///
    color(steelblue) ///
    title("Distribución de empresas según años desde primera entrada al IPSA") ///
    ytitle("Número de empresas") ///
    xtitle("Años desde la primera entrada") ///
    graphregion(color(white)) ///
    ylabel(, angle(0)) ///
    xtick(, grid)
	
	graph export grafico1.png, replace
	
*Primera aproximación de la cantidad de empresas en cada año de tratamiento (primera entrada al ipsa)
*Ahora bien, pueden haber empresas para las cuales tenemos las dummies de entrada/salida pero no tenemos su información financiera...
*A continuación, nos hacemos cargo de eso 

clear all
set more off
use "data_expanded.dta", clear

* Marcar si tienen datos válidos
gen tiene_datos = !missing(totaldebt) & !missing(totalassets)

* Contar cuántos años tiene cada empresa con datos completos en ambas variables
bysort empresa: egen n_validos = total(tiene_datos)

* Filtrar solo empresas con al menos 8 años de datos de totaldebt y totalassets. ESTO LO PODEMOS AJUSTAR... MÁS VARIABLES, MÁS AÑOS
keep if n_validos >= 8

gen event_time = year - first_entry

keep if inrange(event_time, -25, 25) 

save "data_event_time.dta", replace


gen uno = 1
collapse (sum) empresas = uno, by(event_time)

twoway bar empresas event_time, ///
    barwidth(0.8) ///
    color(dkorange) ///
    title("Cantidad de empresas con datos válidos en torno al evento") ///
    ytitle("Número de empresas") ///
    xtitle("Años desde la primera entrada al IPSA") ///
    graphregion(color(white)) ///
    ylabel(, angle(0)) ///
    xtick(, grid)
	
*Se ve un pequeño desbalanceo... Quizás se podría aplicar weights... VER CÓMO NOS PODRÍAMOS HACER CARGO 


*Ahora la idea es evaluar la distribución de los tratados y nunca tratados. 
clear all
set more off

* === Cargar base original con first_entry ===
use "data_expanded.dta", clear



* === Dejar solo una fila por empresa ===
bysort empresa (year): keep if _n == 1

* === Crear etiqueta de entrada ===
gen entrada_label = string(first_entry)
replace entrada_label = "Never Treated" if never_treated

* === Contar empresas por grupo ===
gen uno = 1
collapse (sum) empresas = uno, by(entrada_label)

* === Orden para el gráfico: Never Treated al final ===
gen orden = real(entrada_label)
replace orden = 9999 if entrada_label == "Never Treated"
sort orden

* === Gráfico 1: Tratadas + Never Treated ===
graph bar empresas, over(entrada_label, sort(orden) label(angle(45))) ///
    bar(1, color(cranberry)) ///
    title("Distribución de empresas por año de entrada al IPSA") ///
    ytitle("Número de empresas") ///
    graphregion(color(white)) ///
    legend(off)
graph export grafico_tratados_y_nt.png, replace

* === Gráfico 2: Solo tratadas ===
drop if entrada_label == "Never Treated"
graph bar empresas, over(entrada_label, sort(orden) label(angle(45))) ///
    bar(1, color(navy)) ///
    title("Distribución de empresas tratadas por año de entrada") ///
    ytitle("Número de empresas") ///
    graphregion(color(white)) ///
    legend(off)
graph export grafico_solo_tratados.png, replace




*Incluir en estimaciones los nunca tratados como control 

* === ESTIMACIÓN TWFE ===  
*Ajustar variables por uf. FALTA


use "data_event_time.dta", clear

* Agrupar empresas como identificador numérico
egen empresa_id = group(empresa)

* Variables dependientes
gen debt_ratio  = totaldebt / totalassets
gen capex_ratio = capex / totalassets

save "data_event_time.dta", replace


*Empresas que sean menores 10 años se acumular como si fueran -10, y las mayores a 10 como si fueran 10
gen time_trim = event_time
replace time_trim = -10 if event_time < -10
replace time_trim =  10 if event_time > 10

gen aux = time_trim + 10  // Rango será de 0 a 20

label define aux ///
  0 "-10" 1 "-9" 2 "-8" 3 "-7" 4 "-6" 5 "-5" 6 "-4" 7 "-3" 8 "-2" 9 "-1" ///
  10 "0" 11 "1" 12 "2" 13 "3" 14 "4" 15 "5" 16 "6" 17 "7" 18 "8" 19 "9" 20 "10"
label values aux aux

* Normalizar β_{-1} = 0, que es aux == 10
reghdfe debt_ratio ib10.aux, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
eststo debt

reghdfe capex_ratio ib10.aux, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
eststo capex

reghdfe ebitdeprec~n ib10.aux, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
eststo ebit

* === GRÁFICOS ===

*Debt
coefplot debt, drop(_cons) vertical ///
    rename(1.aux = "-10" 2.aux = "-9" 3.aux = "-8" 4.aux = "-7" 5.aux = "-6" ///
           6.aux = "-5" 7.aux = "-4" 8.aux = "-3" 9.aux = "-2" 10.aux = "-1" ///
           11.aux = "0" 12.aux = "1" 13.aux = "2" 14.aux = "3" 15.aux = "4" ///
           16.aux = "5" 17.aux = "6" 18.aux = "7" 19.aux = "8" 20.aux = "9" ///
           21.aux = "10") ///
    title("Efecto sobre ratio de deuda (β_{-1} = 0)") ///
    yline(0, lpattern(dash)) ///
    ciopts(recast(rcap)) ///
    graphregion(color(white))


* Capex
coefplot capex, drop(_cons) vertical ///
    rename(1.aux = "-10" 2.aux = "-9" 3.aux = "-8" 4.aux = "-7" 5.aux = "-6" ///
           6.aux = "-5" 7.aux = "-4" 8.aux = "-3" 9.aux = "-2" 10.aux = "-1" ///
           11.aux = "0" 12.aux = "1" 13.aux = "2" 14.aux = "3" 15.aux = "4" ///
           16.aux = "5" 17.aux = "6" 18.aux = "7" 19.aux = "8" 20.aux = "9" ///
           21.aux = "10") ///
    title("Efecto sobre Capex / Total Assets (β_{-1} = 0)") ///
    yline(0, lpattern(dash)) ///
    ciopts(recast(rcap)) ///
    graphregion(color(white))


* EBITDA + depreciación
coefplot ebit, drop(_cons) vertical ///
    rename(1.aux = "-10" 2.aux = "-9" 3.aux = "-8" 4.aux = "-7" 5.aux = "-6" ///
           6.aux = "-5" 7.aux = "-4" 8.aux = "-3" 9.aux = "-2" 10.aux = "-1" ///
           11.aux = "0" 12.aux = "1" 13.aux = "2" 14.aux = "3" 15.aux = "4" ///
           16.aux = "5" 17.aux = "6" 18.aux = "7" 19.aux = "8" 20.aux = "9" ///
           21.aux = "10") ///
    title("Efecto sobre EBITDA + Depreciación (β_{-1} = 0)") ///
    yline(0, lpattern(dash)) ///
    ciopts(recast(rcap)) ///
    graphregion(color(white))

*Otra forma de graficar... Está centrado en 0... no se por qué lo anterior no... pero no aparece bien el t=0... 
 

*Deuda
coefplot debt, drop(_cons) vertical ///
    baselevels label ///
    title("Efecto sobre ratio de deuda (β_{-1} = 0)") ///
    yline(0, lpattern(dash)) ///
    ciopts(recast(rcap)) ///
    xlabel(0(1)20, valuelabel angle(45)) ///
    graphregion(color(white))


*Capex
coefplot capex, drop(_cons) vertical ///
    baselevels label ///
    title("Efecto sobre Capex / Total Assets (β_{-1} = 0)") ///
    yline(0, lpattern(dash)) ///
    ciopts(recast(rcap)) ///
    xlabel(0(1)20, valuelabel angle(45)) ///
    graphregion(color(white))

*Ebitda
coefplot ebit, drop(_cons) vertical ///
    baselevels label ///
    title("Efecto sobre EBITDA + Depreciación (β_{-1} = 0)") ///
    yline(0, lpattern(dash)) ///
    ciopts(recast(rcap)) ///
    xlabel(0(1)20, valuelabel angle(45)) ///
    graphregion(color(white))

	

*Recordar que esta levemente desbalanceado... quizás agregar peso...


*A continuación realizamos TWFE pero recortando a -5 y 5 
preserve

* Filtrar solo los años entre -5 y 5 desde el evento
keep if inrange(event_time, -5, 5)

* Recalcular la variable auxiliar
gen aux_trim = event_time + 5   // ahora va de 0 a 10

* Etiquetar correctamente
label define aux_trim 0 "-5" 1 "-4" 2 "-3" 3 "-2" 4 "-1" 5 "0" 6 "1" 7 "2" 8 "3" 9 "4" 10 "5"
label values aux_trim aux_trim


* Estimar modelo TWFE normalizando β_{-1} = 0 → aux_trim==4
reghdfe debt_ratio ib4.aux_trim, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
eststo debt_model3


reghdfe capex_ratio ib9.aux_trim, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
eststo capex_model3

reghdfe ebitdeprec~n ib9.aux_trim, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
eststo ebit_model3

restore

*DEBT
coefplot debt_model3, drop(_cons) vertical ///
    baselevels label ///
    title("TWFE en ventana [-5, 5]: Ratio de deuda (β_{-1} = 0)") ///
    yline(0, lpattern(dash)) ///
    xlabel(0(1)10, valuelabel angle(45)) ///
    ciopts(recast(rcap)) ///
    graphregion(color(white))
	
*CAPEX	
coefplot capex_model3, drop(_cons) vertical ///
    baselevels label ///
    title("TWFE en ventana [-5, 5]: Capex / Total Assets (β_{-1} = 0)") ///
    yline(0, lpattern(dash)) ///
    xlabel(0(1)10, valuelabel angle(45)) ///
    ciopts(recast(rcap)) ///
    graphregion(color(white))

*EBIT
coefplot ebit_model3, drop(_cons) vertical ///
    baselevels label ///
    title("TWFE en ventana [-5, 5]: EBITDA + Depreciación (β_{-1} = 0)") ///
    yline(0, lpattern(dash)) ///
    xlabel(0(1)10, valuelabel angle(45)) ///
    ciopts(recast(rcap)) ///
    graphregion(color(white))


*Estimaciones que buscan capturar dinámicas...

	
*Callaway and Sant'Anna (2020)
*-----------------------------------------------
use "data_event_time.dta", clear

* Variable de tratamiento: año de entrada al IPSA
gen G = first_entry
replace G = 0 if missing(G)   // No tratados


*-----------------------------------------------
* Modelo 1: CSDID - Ratio de deuda
*-----------------------------------------------
eststo M1: csdid debt_ratio, ivar(empresa_id) time(year) gvar(G) notyet long2
estat event, window(-7,7)
csdid_plot, ///
    title("CSDID: Ratio de deuda (ventana [-7, 7])") ///
    yline(0, lpattern(dash)) ///
    graphregion(color(white))

*-----------------------------------------------
* Modelo 2: CSDID - Capex / Total Assets
*-----------------------------------------------
eststo M2: csdid capex_ratio, ivar(empresa_id) time(year) gvar(G) notyet long2
estat event, window(-7,7)
csdid_plot, ///
    title("CSDID: Capex / Total Assets (ventana [-7, 7])") ///
    yline(0, lpattern(dash)) ///
    graphregion(color(white))

*-----------------------------------------------
* Modelo 3: CSDID - EBITDA + Depreciación
*-----------------------------------------------
eststo M3: csdid ebitdeprec~n, ivar(empresa_id) time(year) gvar(G) notyet long2
estat event, window(-7,7)
csdid_plot, ///
    title("CSDID: EBITDA + Depreciación (ventana [-7, 7])") ///
    yline(0, lpattern(dash)) ///
    graphregion(color(white))

*Chequear si es que tiene que estar centrado en 0...
*Se puede hacer otra estimación Robusta a efectos heterogéneos
	





*Realizamos lo anterior pero para salidas
	
	

	


