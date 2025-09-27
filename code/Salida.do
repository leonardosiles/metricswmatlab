*===========================================
* Salida Final - Paso 1: Base limpia para análisis de salida
*===========================================

clear all
set more off

* --- Cargar base original ---
import excel "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/IGPA_Trabajar.xlsx", ///
    sheet("Sheet 1") firstrow clear

* --- Marcar si estuvo alguna vez en el IPSA (ya sea entrando o manteniéndose) ---
gen estuvo_en_ipsa = (entra == 1 | mantiene == 1)
egen total_ipsa = total(estuvo_en_ipsa), by(empresa)

* --- Marcar si alguna vez salió del IPSA ---
egen total_salidas = total(sale == 1), by(empresa)

* --- Definir grupos ---
gen treated = (total_salidas > 0 & total_ipsa > 0)          // Tratados: estuvieron y salieron
gen never_exited = (total_salidas == 0 & total_ipsa > 0)     // Controles: estuvieron pero nunca salieron

* --- Quitar empresas que nunca estuvieron en el IPSA ---
keep if total_ipsa > 0

* --- Identificar el primer año en que salen del IPSA ---
egen first_exit = min(year / (sale == 1)), by(empresa)

* --- Rellenar hacia abajo ---
bysort empresa (year): replace first_exit = first_exit[_n-1] if missing(first_exit)

* --- Calcular el año relativo al evento de salida ---
gen event_time_exit = year - first_exit

* --- Guardar base intermedia para análisis de salida ---
save "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_raw.dta", replace


*===========================================
* Salida Final - Paso 2: Distribución general alrededor de la salida
*===========================================

use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_raw.dta", clear

* --- Filtrar empresas que salieron y están dentro de ventana [-15,15] ---
keep if (treated == 1 & inrange(event_time_exit, -15, 15)) | never_exited == 1

* --- Crear variable de conteo ---
gen uno = 1

* --- Colapsar: contar número de empresas por event_time_exit ---
collapse (sum) empresas = uno, by(event_time_exit)

* --- Graficar distribución SOLO de las tratadas (excluir missing) ---
twoway bar empresas event_time_exit if !missing(event_time_exit), ///
    barwidth(0.8) ///
    color(grey) ///
    title("Distribución empresas según años primera salida IPSA") ///
    ytitle("Número de empresas") ///
    xtitle("Años desde la salida") ///
    graphregion(color(white)) ///
    ylabel(, angle(0)) ///
    xtick(, grid)

*===========================================
* Salida Final - Paso 3: Filtrar empresas con datos válidos (10 años)
*===========================================

use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_raw.dta", clear

* --- Mantener empresas que salieron dentro de la ventana [-15,15] y las de control ---
keep if (treated == 1 & inrange(event_time_exit, -15, 15)) | never_exited == 1

* --- Crear indicador de validez de datos financieros ---
gen tiene_datos = !missing(totalassets) & !missing(totaldebt)

* --- Contar años válidos por empresa ---
bysort empresa: egen n_validos = total(tiene_datos)

* --- Filtrar empresas con al menos 10 años válidos ---
keep if n_validos >= 10

*===========================================
* Salida Final - Paso 4: Distribución empresas tratadas con datos válidos
*===========================================

use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_validos10.dta", clear

* --- Dejar solo empresas tratadas ---
keep if treated == 1 & !missing(event_time_exit)

* --- Variable de conteo ---
gen uno = 1

* --- Contar empresas por año relativo ---
collapse (sum) empresas = uno, by(event_time_exit)

* --- Gráfico de distribución ---
twoway bar empresas event_time_exit, ///
    barwidth(0.8) ///
    color(grey) ///
    title("Distribución empresas válidas según años desde salida") ///
    ytitle("Número de empresas") ///
    xtitle("Años desde la salida") ///
    graphregion(color(white)) ///
    ylabel(, angle(0)) ///
    xtick(, grid)
	
*===========================================
* Salida Final - Paso 5: Distribución por cohorte de salida
*===========================================

use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_validos10.dta", clear

* --- Dejar una fila por empresa ---
bysort empresa (year): keep if _n == 1

* --- Etiqueta con año de primera salida ---
gen salida_label = string(first_exit)
replace salida_label = "Nunca salió" if treated == 0

* --- Variable para conteo ---
gen uno = 1

* --- Colapsar por cohorte ---
collapse (sum) empresas = uno, by(salida_label)

* --- Ordenar: primero los años, luego "Nunca salió" al final ---
gen orden = real(salida_label)
replace orden = 9999 if salida_label == "Nunca salió"
sort orden

* --- Gráfico de barras por cohorte ---
graph bar empresas, over(salida_label, sort(orden) label(angle(45))) ///
    bar(1, color(grey)) ///
    title("Distribución empresas por año de primera salida del IPSA") ///
    ytitle("Número de empresas") ///
    graphregion(color(white)) ///
    ylabel(, angle(0)) ///
    legend(off)

*===========================================
* Salida Final - Paso 6: Ajuste a UF y construcción de ratios
*===========================================

use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_validos10.dta", clear

* --- Ajustar variables a UF ---
foreach var in ///
    totalassets ppe cashshortterminvestments totalcapital totaldebt ///
    longtermdebt totalshareholdersequity netcashflow_operatingactivs ///
    cashdividendspaid_total netsalesorrevenues capex ///
    depreciationanddepletion depreciationdepletionamort ///
    marketvalue marketcapitalization netdebt accumdeprc_othppe ///
    commonshareholdersequity totalliabilitiesshareholde ///
    commonsharesoutstanding ebitdepreciation earningsbefinteresttaxes ///
    interestexpenseondebt interestcapitalized currentliabilities_total ///
    currentassets_total researchdevelopment interestexpense_total ///
	interestincome_total   {

    gen `var'_uf = `var' / uf
}

* --- Crear identificador numérico de empresa ---
egen empresa_id = group(empresa)



* Liquidez
gen current_rat = currentassets_total_uf / currentliabilities_total_uf
gen k_trab = currentassets_total_uf - currentliabilities_total_uf
gen cash_ratio = cashshortterminvestments_uf / currentliabilities_total_uf
gen flujoefectivo_ventas = cashflowsales
gen cashflow_debt = netcashflow_operatingactivs_uf / totaldebt_uf
gen deuda_flujo = totaldebt_uf / netcashflow_operatingactivs_uf
gen gasto_financiero_fco = interestexpenseondebt_uf / netcashflow_operatingactivs_uf

* Endeudamiento
gen debt_assets = totaldebt_uf / totalassets_uf
gen deuda_capital = totaldebt_uf / totalcapital_uf
gen debt_equity = totaldebt_uf / totalshareholdersequity_uf
gen debt_equity2 = totaldebt_uf / commonshareholdersequity
gen deuda_equity2 = totaldebtcommonequity  // si existe esta variable
gen netdebt_ebitda = netdebt_uf / ebitdepreciation_uf
gen debt_ratio2 = totaldebttotalassets     // si existe esta variable
gen dk = totaldebttotalcapitalstd          // si existe esta variable
gen deudalarga_ratio = longtermdebt_uf / totalassets_uf
gen leverage_financiero = earningsbefinteresttaxes_uf / (earningsbefinteresttaxes_uf - interestexpenseondebt_uf)
gen Int_cove_ratio = earningsbefinteresttaxes_uf / interestexpense_total_uf
gen Int_cove_ratio2 = earningsbefinteresttaxes_uf / interestexpenseondebt_uf
gen pasivoscp = currentliabilities_total_uf / totalliabilitiesshareholde_uf

* Rentabilidad
gen margen_ebit_sales = earningsbefinteresttaxes_uf / netsalesorrevenues_uf
gen ebitda_ventas = ebitdepreciation_uf / netsalesorrevenues_uf
gen ROA = earningsbefinteresttaxes_uf / totalassets_uf
gen ROA_ebitda = ebitdepreciation_uf / totalassets_uf
gen ROE = returnonequity_total
gen retorno_capital = returnonequity_total
gen margen_operacional = operatingprofitmargin
gen margen_pretax = pretaxmargin
gen margen_bruto = grossprofitmargin

* Actividad / estructura
gen capex_assets = capex_uf / totalassets_uf
gen capex_depr = capex_uf / depreciationanddepletion_uf
gen ppe_assets = ppe_uf / totalassets_uf

* Valorización bursátil
gen market_to_book = marketcapitalization_uf / totalshareholdersequity_uf
gen dvds = divpershr




* --- Guardar base final con UF y ratios listos ---
save "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_final_ratios.dta", replace


*===========================================
* Salida Final - Paso 7: TWFE acumulado + recorte puro (loop por ratio)
*===========================================

clear all
set more off

* --- Definir lista de ratios a graficar ---
local ratios ROE ROA ROA_ebitda retorno_capital ///
             debt_assets debt_equity debt_equity2 deuda_capital deuda_equity2 ///
             netdebt_ebitda debt_ratio2 dk deudalarga_ratio pasivoscp ///
             cash_ratio current_rat k_trab cashflow_debt deuda_flujo ///
             flujoefectivo_ventas gasto_financiero_fco ///
             margen_ebit_sales ebitda_ventas margen_operacional margen_pretax  margen_bruto ///
			 retorno_capital leverage_financiero Int_cove_ratio Int_cove_ratio2 ///
             capex_assets capex_depr ppe_assets ///
             market_to_book dvds
	
* --- Loop por ratio ---
foreach var of local ratios {

    di "==> Ejecutando TWFE para `var'..."

    * --- Cargar base con ratios y UF ---
    use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_final_ratios.dta", clear

    *==============================
    * FILTRO: Tratadas con datos NO PERDIDOS en los años [-2,-1,0,1,2] del evento de salida
    * Se aplica solo a treated == 1 (las que salieron del IPSA)
    * Las controles (treated == 0) se mantienen siempre
    *==============================
    gen usable = treated & !missing(`var') & inlist(event_time_exit, -2,-1, 0, 1,2)
    bysort empresa: egen usable_count = total(usable)
    keep if (treated & usable_count == 5) | !treated
    drop usable usable_count

*==============================
* Estrategia 1: Acumulado en extremos [-5,5] (SALIDA)
*==============================
gen time_trim = event_time_exit
replace time_trim = -5 if event_time_exit < -5
replace time_trim =  5 if event_time_exit > 5

gen aux_acum = .
replace aux_acum = time_trim + 5 if !missing(time_trim)   // -5..5 -> 0..10
replace aux_acum = 11 if treated == 0                     // control

label define aux_lbl 0 "-5" 1 "-4" 2 "-3" 3 "-2" 4 "-1" 5 "0" 6 "1" 7 "2" 8 "3" 9 "4" 10 "5" 11 "Control", replace
label values aux_acum aux_lbl

quietly reghdfe `var' ib4.aux_acum, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
eststo acum

* === Eje X -5..5 y nombre en memoria compatible con 'graph combine'
local coefsS
local lblS
forvalues k = 0/10 {
    local coefsS `coefsS' `k'.aux_acum
    local et = `k' - 5
    local lblS `lblS' `k'.aux_acum = "`et'"
}

coefplot acum, drop(_cons 11.aux_acum) keep(`coefsS') order(`coefsS') ///
    vertical baselevels coeflabels(`lblS') xlabel(, angle(0)) ///
    title("`var': Acumulado [-5,5] (salida)") ///
    yline(0, lpattern(dash) lcolor(black)) ciopts(recast(rcap)) ///
    name(g_acum, replace) graphregion(color(white))   // <— usa g_acum
	xscale(range(0 10)) plotregion(margin(zero))

* === Guardar en la misma carpeta y formato que ENTRADA (prefijo S_)
local out_base "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/Graf final"
cap mkdir "`out_base'/gph"
graph save   "`out_base'/gph/S_`var'_acum.gph", replace
graph export "`out_base'/S_`var'_acum.png", replace width(2400)




    *==============================
    * Estrategia 2: Recorte puro [-5,5]
    *==============================
    gen dentro_rango = inrange(event_time_exit, -5, 5)
    keep if dentro_rango == 1 | treated == 0

    gen aux_rec = .
    replace aux_rec = event_time_exit + 5 if dentro_rango == 1
    replace aux_rec = 11 if missing(aux_rec)

    levelsof aux_rec, local(vals)
    local lbls
    foreach v of local vals {
        if `v' == 11 {
            local lbls `lbls' `v' "Control"
        }
        else {
            local etiq = `v' - 5
            local lbls `lbls' `v' "`etiq'"
        }
    }
    label define aux_lbl `lbls', replace
    label values aux_rec aux_lbl

    quietly reghdfe `var' ib4.aux_rec, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
    eststo rec

    coefplot rec, drop(_cons) vertical ///
        baselevels label ///
        title("`var': Recorte puro [-5,5]") ///
        yline(0, lpattern(dash) lcolor(black)) ///
        ciopts(recast(rcap)) ///
        xlabel(, angle(45)) ///
        name(g_rec, replace) ///
        graphregion(color(white))

    *==============================
    * COMBINAR Y EXPORTAR GRÁFICO
    *==============================
    graph combine g_acum g_rec, col(2) ///
        title("Event Study: Salida del IPSA sobre `var'") ///
        ycommon graphregion(color(white)) ///
        iscale(*0.8)

    graph export "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/Graf final/S_`var'_5.png", replace
}


*Con betas


*============================================================*
* SALIDA (-5,5): acumulado y recortado                       *
* + Opción A: pendientes pre/post sobre betas estimadas      *
* + Opción B: slopes pre/post en la data (TWFE formal)       *
*============================================================*

clear all
set more off

cap which reghdfe
if _rc ssc install reghdfe, replace
cap which coefplot
if _rc ssc install coefplot, replace
cap which estout
if _rc ssc install estout, replace

*===========================================
* Paso 1: Base limpia para análisis de salida
*===========================================
import excel "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/IGPA_Trabajar.xlsx", ///
    sheet("Sheet 1") firstrow clear

gen estuvo_en_ipsa = (entra == 1 | mantiene == 1)
egen total_ipsa = total(estuvo_en_ipsa), by(empresa)
egen total_salidas = total(sale == 1), by(empresa)

gen treated      = (total_salidas > 0 & total_ipsa > 0)      // tratadas: estuvieron y salieron
gen never_exited = (total_salidas == 0 & total_ipsa > 0)     // controles: estuvieron y nunca salieron

keep if total_ipsa > 0

egen first_exit = min(year / (sale == 1)), by(empresa)
bysort empresa (year): replace first_exit = first_exit[_n-1] if missing(first_exit)

gen event_time_exit = year - first_exit

save "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_raw.dta", replace

*===========================================
* Paso 2: Distribución general alrededor de la salida
*===========================================
use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_raw.dta", clear

keep if (treated == 1 & inrange(event_time_exit, -15, 15)) | never_exited == 1
gen uno = 1
collapse (sum) empresas = uno, by(event_time_exit)

twoway bar empresas event_time_exit if !missing(event_time_exit), ///
    barwidth(0.8) color(gray) ///
    title("Distribución empresas según años primera salida IPSA") ///
    ytitle("Número de empresas") xtitle("Años desde la salida") ///
    graphregion(color(white)) ylabel(, angle(0)) xtick(, grid)

*===========================================
* Paso 3: Filtrar empresas con datos válidos (10 años)
*===========================================
use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_raw.dta", clear
keep if (treated == 1 & inrange(event_time_exit, -15, 15)) | never_exited == 1

gen tiene_datos = !missing(totalassets) & !missing(totaldebt)
bysort empresa: egen n_validos = total(tiene_datos)
keep if n_validos >= 10

save "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_validos10.dta", replace

*===========================================
* Paso 4: Distribución empresas tratadas válidas
*===========================================
use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_validos10.dta", clear
keep if treated == 1 & !missing(event_time_exit)
gen uno = 1
collapse (sum) empresas = uno, by(event_time_exit)

twoway bar empresas event_time_exit, ///
    barwidth(0.8) color(gray) ///
    title("Distribución empresas válidas según años desde salida") ///
    ytitle("Número de empresas") xtitle("Años desde la salida") ///
    graphregion(color(white)) ylabel(, angle(0)) xtick(, grid)

*===========================================
* Paso 5: Distribución por cohorte de salida
*===========================================
use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_validos10.dta", clear
bysort empresa (year): keep if _n == 1

gen salida_label = string(first_exit)
replace salida_label = "Nunca salió" if treated == 0

gen uno = 1
collapse (sum) empresas = uno, by(salida_label)

gen orden = real(salida_label)
replace orden = 9999 if salida_label == "Nunca salió"
sort orden

graph bar empresas, over(salida_label, sort(orden) label(angle(45))) ///
    bar(1, color(gray)) ///
    title("Distribución empresas por año de primera salida del IPSA") ///
    ytitle("Número de empresas") graphregion(color(white)) ///
    ylabel(, angle(0)) legend(off)

*===========================================
* Paso 6: Ajuste a UF y construcción de ratios
*===========================================
use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_validos10.dta", clear

foreach var in ///
    totalassets ppe cashshortterminvestments totalcapital totaldebt ///
    longtermdebt totalshareholdersequity netcashflow_operatingactivs ///
    cashdividendspaid_total netsalesorrevenues capex ///
    depreciationanddepletion depreciationdepletionamort ///
    marketvalue marketcapitalization netdebt accumdeprc_othppe ///
    commonshareholdersequity totalliabilitiesshareholde ///
    commonsharesoutstanding ebitdepreciation earningsbefinteresttaxes ///
    interestexpenseondebt interestcapitalized currentliabilities_total ///
    currentassets_total researchdevelopment interestexpense_total ///
    interestincome_total {
    gen `var'_uf = `var' / uf
}

egen empresa_id = group(empresa)

* Liquidez (con guardas de división)
gen current_rat = currentassets_total_uf / currentliabilities_total_uf
gen k_trab = currentassets_total_uf - currentliabilities_total_uf
gen cash_ratio = cashshortterminvestments_uf / currentliabilities_total_uf
gen flujoefectivo_ventas = cashflowsales
gen cashflow_debt = netcashflow_operatingactivs_uf / totaldebt_uf
replace cashflow_debt = . if totaldebt_uf==0 | missing(totaldebt_uf)

gen deuda_flujo = .
replace deuda_flujo = totaldebt_uf / netcashflow_operatingactivs_uf if ///
    netcashflow_operatingactivs_uf!=. & netcashflow_operatingactivs_uf!=0

gen gasto_financiero_fco = interestexpenseondebt_uf / netcashflow_operatingactivs_uf
replace gasto_financiero_fco = . if netcashflow_operatingactivs_uf==0 | missing(netcashflow_operatingactivs_uf)

* Endeudamiento
gen debt_assets = totaldebt_uf / totalassets_uf
gen deuda_capital = totaldebt_uf / totalcapital_uf
gen debt_equity = totaldebt_uf / totalshareholdersequity_uf
gen debt_equity2 = totaldebt_uf / commonshareholdersequity
gen deuda_equity2 = totaldebtcommonequity
gen netdebt_ebitda = netdebt_uf / ebitdepreciation_uf
gen debt_ratio2 = totaldebttotalassets
gen dk = totaldebttotalcapitalstd
gen deudalarga_ratio = longtermdebt_uf / totalassets_uf
gen leverage_financiero = earningsbefinteresttaxes_uf / (earningsbefinteresttaxes_uf - interestexpenseondebt_uf)
gen Int_cove_ratio = earningsbefinteresttaxes_uf / interestexpense_total_uf
gen Int_cove_ratio2 = earningsbefinteresttaxes_uf / interestexpenseondebt_uf
gen pasivoscp = currentliabilities_total_uf / totalliabilitiesshareholde_uf

* Rentabilidad
gen margen_ebit_sales = earningsbefinteresttaxes_uf / netsalesorrevenues_uf
gen ebitda_ventas = ebitdepreciation_uf / netsalesorrevenues_uf
gen ROA = earningsbefinteresttaxes_uf / totalassets_uf
gen ROA_ebitda = ebitdepreciation_uf / totalassets_uf
gen ROE = returnonequity_total
gen retorno_capital = returnonequity_total
gen margen_operacional = operatingprofitmargin
gen margen_pretax = pretaxmargin
gen margen_bruto = grossprofitmargin

* Actividad / estructura
gen capex_assets = capex_uf / totalassets_uf
gen capex_depr = capex_uf / depreciationanddepletion_uf
gen ppe_assets = ppe_uf / totalassets_uf

* Valorización bursátil
gen market_to_book = marketcapitalization_uf / totalshareholdersequity_uf
gen dvds = divpershr

save "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_final_ratios.dta", replace

*===========================================
* Paso 7: TWFE + Gráficos + Opción A + Opción B
*===========================================
clear all
set more off

local ratios ROE ROA ROA_ebitda retorno_capital ///
             debt_assets debt_equity debt_equity2 deuda_capital deuda_equity2 ///
             netdebt_ebitda debt_ratio2 dk deudalarga_ratio pasivoscp ///
             cash_ratio current_rat k_trab cashflow_debt deuda_flujo ///
             flujoefectivo_ventas gasto_financiero_fco ///
             margen_ebit_sales ebitda_ventas margen_operacional margen_pretax margen_bruto ///
             retorno_capital leverage_financiero Int_cove_ratio Int_cove_ratio2 ///
             capex_assets capex_depr ppe_assets ///
             market_to_book dvds

* Resultados Opción A (pendientes sobre betas)
tempfile slopesA_exit_out
postfile SLa str30 ratio str4 spec double b_pre se_pre b_post se_post using `slopesA_exit_out', replace

* Resultados Opción B (slopes formales en la data)
tempfile slopesB_exit_out
postfile SLb str30 ratio double beta_pre se_pre beta_post se_post delta se_delta diff_post_pre se_diff using `slopesB_exit_out', replace

foreach var of local ratios {

    di as txt "==> Event-study SALIDA para `var'..."

    use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_exit_event_final_ratios.dta", clear

    * Filtro de calidad (tratadas con datos en -2,-1,0,1,2; controles se mantienen)
    gen usable = (treated==1) & !missing(`var') & inlist(event_time_exit,-2,-1,0,1,2)
    bysort empresa: egen usable_count = total(usable)
    keep if (treated & usable_count==5) | !treated
    drop usable usable_count

 *==============================
* ESTRATEGIA 1: ACUMULADO [-5,5] (SALIDA)
*==============================
gen time_trim = event_time_exit
replace time_trim = -5 if event_time_exit < -5
replace time_trim =  5 if event_time_exit > 5

gen aux_acum = .
replace aux_acum = time_trim + 5 if !missing(time_trim)      // -5..5 -> 0..10
replace aux_acum = 11 if treated==0                          // control

label define aux_lbl 0 "-5" 1 "-4" 2 "-3" 3 "-2" 4 "-1" 5 "0" 6 "1" 7 "2" 8 "3" 9 "4" 10 "5" 11 "Control", replace
label values aux_acum aux_lbl

quietly reghdfe `var' ib4.aux_acum, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
eststo acum

* === EJE X LIMPIO: -5..5 (sin control 11)
local coefsS
local lblS
forvalues k = 0/10 {
    local coefsS `coefsS' `k'.aux_acum
    local et = `k' - 5
    local lblS  `lblS'  `k'.aux_acum = "`et'"
}

coefplot acum, drop(_cons 11.aux_acum) keep(`coefsS') order(`coefsS') ///
    vertical baselevels coeflabels(`lblS') xlabel(, angle(0)) ///
    title("`var': Acumulado [-5,5] (salida)") ///
    yline(0, lpattern(dash) lcolor(black)) ciopts(recast(rcap)) ///
    name(g_acum, replace) graphregion(color(white))
	

* ---- Opción A: pendientes sobre betas (acum)
matrix b = e(b)
matrix V = e(V)
local K : colnames b

tempfile tb_acum
postfile PP1 double beta se int level using `tb_acum', replace
foreach c of local K {
    if strpos("`c'","aux_acum") {
        local dot = strpos("`c'",".")
        if `dot'>1 {
            local lev = real(substr("`c'",1,`dot'-1))
            if `lev'!=11 {                         // excluir control
                local col = colnumb(b,"`c'")
                scalar bb = b[1,`col']
                scalar ss = sqrt(V[`col',`col'])
                post PP1 (bb) (ss) (`lev')
            }
        }
    }
}
postclose PP1

preserve
    capture use `tb_acum', clear
    if _rc==0 {
        gen event_time = level - 5   // 0..10 -> -5..5

        tempfile base_acum
        save `base_acum', replace

        * PRE = [-5,-2]  (excluye -1)
        use `base_acum', clear
        keep if inrange(event_time,-5,-2)
        scalar bpre  = .
        scalar sepre = .
        count
        if r(N)>=2 {
            gen w = 1/(se^2)
            replace w = . if se<=0
            regress beta event_time [aw=w]
            scalar bpre  = _b[event_time]
            scalar sepre = _se[event_time]
        }

        * POST = [0,5]
        use `base_acum', clear
        keep if inrange(event_time,0,5)
        scalar bpost  = .
        scalar sepost = .
        count
        if r(N)>=2 {
            gen w = 1/(se^2)
            replace w = . if se<=0
            regress beta event_time [aw=w]
            scalar bpost  = _b[event_time]
            scalar sepost = _se[event_time]
        }

        post SLa ("`var'") ("acum") (bpre) (sepre) (bpost) (sepost)
    }
restore

    *==============================
    * ESTRATEGIA 2: RECORTE [-5,5]
    *==============================
    gen dentro_rango = inrange(event_time_exit,-5,5)
    keep if dentro_rango==1 | treated==0

    gen aux_rec = .
    replace aux_rec = event_time_exit + 5 if dentro_rango==1     // 0..10
    replace aux_rec = 11 if missing(aux_rec)                      // control

    levelsof aux_rec, local(vals)
    local lbls
    foreach v of local vals {
        if `v'==11 local lbls `lbls' `v' "Control"
        else {
            local etiq = `v' - 5
            local lbls `lbls' `v' "`etiq'"
        }
    }
    label define aux_lbl `lbls', replace
    label values aux_rec aux_lbl

    quietly reghdfe `var' ib4.aux_rec, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
    eststo rec

    coefplot rec, drop(_cons) vertical baselevels label ///
        title("`var': Recorte puro [-5,5] (salida)") ///
        yline(0, lpattern(dash) lcolor(black)) ciopts(recast(rcap)) ///
        xlabel(, angle(45)) name(g_rec, replace) graphregion(color(white))

    * ---- Opción A: pendientes sobre betas (rec)
    matrix b = e(b)
    matrix V = e(V)
    local K : colnames b

    tempfile tb_rec
    postfile PP2 double beta se int level using `tb_rec', replace
    foreach c of local K {
        if strpos("`c'","aux_rec") {
            local dot = strpos("`c'",".")
            if `dot'>1 {
                local lev = real(substr("`c'",1,`dot'-1))
                if `lev'!=11 {
                    local col = colnumb(b,"`c'")
                    scalar bb = b[1,`col']
                    scalar ss = sqrt(V[`col',`col'])
                    post PP2 (bb) (ss) (`lev')
                }
            }
        }
    }
    postclose PP2

    preserve
        capture use `tb_rec', clear
        if _rc==0 {
            gen event_time = level - 5

            tempfile base_rec
            save `base_rec', replace

* PRE = [-5,-2]
use `base_rec', clear
keep if inrange(event_time,-5,-2)
scalar bpre  = .
scalar sepre = .
count
if r(N)>=2 {
    gen w = 1/(se^2)
    replace w = . if se<=0
    regress beta event_time [aw=w]
    scalar bpre  = _b[event_time]
    scalar sepre = _se[event_time]
}

* POST = [0,5]
use `base_rec', clear
keep if inrange(event_time,0,5)
scalar bpost  = .
scalar sepost = .
count
if r(N)>=2 {
    gen w = 1/(se^2)
    replace w = . if se<=0
    regress beta event_time [aw=w]
    scalar bpost  = _b[event_time]
    scalar sepost = _se[event_time]
}


            post SLa ("`var'") ("rec") (bpre) (sepre) (bpost) (sepost)
        }
    restore

    *==============================
    * COMBINAR Y EXPORTAR GRÁFICO
    *==============================
    graph combine g_acum g_rec, col(2) ///
        title("Event Study (SALIDA) sobre `var'") ///
        ycommon graphregion(color(white)) iscale(*0.8)
    graph export "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/Graf final/S_`var'_5.png", replace

    *============================================================
    * Opción B: slopes formales en la DATA (TWFE, cluster empresa)
    *============================================================
    * piecewise sobre event_time_exit; 0 para controles (se quedan en la muestra)
* PRE = [-5,-2]; POST = [0,5]
gen tpre_treat  = cond(treated==1 & !missing(event_time_exit) & inrange(event_time_exit,-5,-2), event_time_exit, 0)
gen tpost_treat = cond(treated==1 & !missing(event_time_exit) & inrange(event_time_exit,0,5),  event_time_exit, 0)

* Mantén el salto (nivel) en t>=0 si te interesa capturarlo
gen post_treat  = cond(treated==1 & !missing(event_time_exit) & event_time_exit>=0, 1, 0)


    quietly reghdfe `var' tpre_treat tpost_treat post_treat, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons

    scalar bpreB   = _b[tpre_treat]
    scalar sepreB  = _se[tpre_treat]
    scalar bpostB  = _b[tpost_treat]
    scalar sepostB = _se[tpost_treat]
    scalar bdelB   = _b[post_treat]
    scalar sedelB  = _se[post_treat]

    lincom tpost_treat - tpre_treat
    scalar bdiffB  = r(estimate)
    scalar sediffB = r(se)

    post SLb ("`var'") (bpreB) (sepreB) (bpostB) (sepostB) (bdelB) (sedelB) (bdiffB) (sediffB)

    drop tpre_treat tpost_treat post_treat
}

*==== Exportar resultados (SALIDA) ====
postclose SLa
use `slopesA_exit_out', clear
order ratio spec b_pre se_pre b_post se_post
export delimited using "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/Graf final/pendientes_sobre_betas_exit.csv", replace

postclose SLb
use `slopesB_exit_out', clear
order ratio beta_pre se_pre beta_post se_post delta se_delta diff_post_pre se_diff
export delimited using "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/Graf final/tendencias_pre_post_B_exit.csv", replace


