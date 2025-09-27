*******************************************************
* ENTRADA FINAL — Event-study [-5,5], acumulado y recorte
* Con TWFE, gráficos y dos opciones de "slopes"
* (A) pendientes sobre betas estimadas
* (B) pendientes en la data (TWFE formal)
*******************************************************

clear all
set more off

*------------------------------------------------------
* Paquetes
*------------------------------------------------------
cap which reghdfe
if _rc ssc install reghdfe, replace
cap which coefplot
if _rc ssc install coefplot, replace
cap which estout
if _rc ssc install estout, replace

*******************************************************
* 1) Distribución general de empresas en torno a ENTRADA
*******************************************************
import excel "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/IGPA_Trabajar.xlsx", ///
    sheet("Sheet 1") firstrow clear

* Identificar primer año de entrada al IPSA
egen first_entry = min(year / (entra == 1)), by(empresa)

* Marcar empresas que alguna vez estuvieron en IPSA
gen ever_in_ipsa = (entra == 1 | mantiene == 1)
egen suma_in_ipsa = total(ever_in_ipsa), by(empresa)
gen never_treated = (suma_in_ipsa == 0)

* Rellenar first_entry hacia abajo
bysort empresa (year): replace first_entry = first_entry[_n-1] if missing(first_entry)

* Calcular año relativo a la primera entrada
gen event_time = year - first_entry

* Guardar base restringida a [-15,15] + controles (missing event_time)
keep if inrange(event_time, -15, 15) | missing(event_time)
save "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_event_-15_15.dta", replace

* Distribución general (tratadas)
gen uno = 1
collapse (sum) empresas = uno, by(event_time)
twoway bar empresas event_time if !missing(event_time), ///
    barwidth(0.8) color(gray) ///
    title("Distribución empresas según año desde primera entrada") ///
    ytitle("Número de empresas") xtitle("Años desde la primera entrada al IPSA") ///
    graphregion(color(white)) ylabel(, angle(0)) xtick(, grid)

*******************************************************
* 2) Filtro de calidad: >=10 años válidos en assets y debt
*******************************************************
clear all
set more off

use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_event_-15_15.dta", clear

gen tiene_datos = !missing(totalassets) & !missing(totaldebt)
bysort empresa: egen n_validos = total(tiene_datos)
keep if n_validos >= 10
save "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_event_-15_15_validos10.dta", replace

* Distribución con válidas (tratadas)
gen uno = 1
collapse (sum) empresas = uno, by(event_time)
twoway bar empresas event_time if !missing(event_time), ///
    barwidth(0.8) color(gray) ///
    title("Distribución empresas válidas según año desde primera entrada") ///
    ytitle("Número de empresas") xtitle("Años desde la primera entrada al IPSA") ///
    graphregion(color(white)) ylabel(, angle(0)) xtick(, grid)

*******************************************************
* 3) Distribución por cohortes de entrada
*******************************************************
use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_event_-15_15_validos10.dta", clear
bysort empresa (year): keep if _n == 1

gen entrada_label = string(first_entry)
replace entrada_label = "Never Treated" if never_treated

gen uno = 1
collapse (sum) empresas = uno, by(entrada_label)

gen orden = real(entrada_label)
replace orden = 9999 if entrada_label == "Never Treated"
sort orden

graph bar empresas, over(entrada_label, sort(orden) label(angle(45))) ///
    bar(1, color(gray)) ///
    title("Distribución de empresas por año de entrada al IPSA") ///
    ytitle("Número de empresas") ylabel(0(5)45, angle(0)) ///
    graphregion(color(white)) legend(off)

*******************************************************
* 4) Ajuste a UF y construcción de RATIOS
*******************************************************
clear all
set more off

use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_event_-15_15_validos10.dta", clear

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

* Liquidez (con guardas mínimas)
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

* Endeudamiento (algunas pueden venir ya como ratios)
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

save "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_event_final_ratios.dta", replace

*******************************************************
* 5) LOOP TWFE + Gráficos + Opción A + Opción B
*******************************************************
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

* Resultados Opción A (pendientes sobre betas del event-study)
tempfile slopesA_out
postfile SLa str30 ratio str4 spec double b_pre se_pre b_post se_post using `slopesA_out', replace

* Resultados Opción B (slopes formales en la data)
tempfile slopesB_out
postfile SLb str30 ratio double beta_pre se_pre beta_post se_post delta se_delta diff_post_pre se_diff using `slopesB_out', replace

foreach var of local ratios {

    di as txt "==> Ejecutando TWFE ENTRADA para `var'..."

    use "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/base_event_final_ratios.dta", clear

    * Tratadas = tienen event_time (controles quedan con missing)
    gen treated = !missing(event_time)

    * Filtro de calidad por ratio: presencia en -2,-1,0,1,2 para tratadas
    gen usable = treated & !missing(`var') & inlist(event_time,-2,-1,0,1,2)
    bysort empresa: egen usable_count = total(usable)
    keep if (treated & usable_count==5) | !treated
    drop usable usable_count

    eststo drop _all

*==============================
* ESTRATEGIA 1: ACUMULADO [-5,5]
*==============================
gen aux_acum = cond(missing(event_time), ., ///
                cond(event_time < -5, 0, ///
                cond(event_time > 5, 10, ///
                event_time + 5)))
replace aux_acum = 99 if missing(aux_acum)      // controles

label define aux_lbl 0 "-5" 1 "-4" 2 "-3" 3 "-2" 4 "-1" 5 "0" 6 "1" 7 "2" 8 "3" 9 "4" 10 "5" 99 "control", replace
label values aux_acum aux_lbl

quietly reghdfe `var' ib4.aux_acum, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
eststo acum

* === Etiquetas limpias para el eje X: -5..5 (y se saca el control 99)
local lblE
forvalues k = 0/10 {
    local et = `k' - 5
    local lblE `lblE' `k'.aux_acum = "`et'"
}

coefplot acum, drop(_cons 99.aux_acum) vertical baselevels ///
    coeflabels(`lblE') xlabel(, angle(0)) ///
    title("`var': Acumulado [-5,5] (entrada)") ///
    yline(0, lpattern(dash) lcolor(black)) ///
    ciopts(recast(rcap)) ///
    name(g_acum, replace) graphregion(color(white)) ///
	xscale(range(0 10)) plotregion(margin(zero))



* ---- Opción A: pendientes sobre betas (ACUM)
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
            if `lev'!=99 {
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
        gen event_time = level - 5
        tempfile base_acum
        save `base_acum', replace

        * PRE = [-5,-2] (excluye -1)
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

    * === Guardar .gph del acumulado de ENTRADA (con eje -5..5)
    local gph_dir "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/Graf final/gph"
    cap mkdir "`gph_dir'"

    * (reconstruyo el macro por si se perdió)
    local lblE
    forvalues k = 0/10 {
        local et = `k' - 5
        local lblE `lblE' `k'.aux_acum = "`et'"
    }

    coefplot acum, drop(_cons 99.aux_acum) vertical baselevels ///
        coeflabels(`lblE') xlabel(, angle(0)) ///
        title("`var': Acumulado [-5,5] (entrada)") ///
        yline(0, lpattern(dash) lcolor(black)) ///
        ciopts(recast(rcap)) ///
        name(gE_acum, replace) graphregion(color(white)) ///
		xscale(range(0 10)) plotregion(margin(zero))


    graph save "`gph_dir'/E_`var'_acum.gph", replace
restore


    *==============================
    * ESTRATEGIA 2: RECORTE [-5,5]
    *==============================
    keep if inrange(event_time, -5, 5) | missing(event_time)

    gen aux_rec = .
    replace aux_rec = event_time + 5 if !missing(event_time)
    replace aux_rec = 99 if missing(event_time)      // controles

    levelsof aux_rec if aux_rec != 99 & !missing(aux_rec), local(vals)
    local lbls
    foreach v of local vals {
        local etiq = `v' - 5
        local lbls `lbls' `v' "`etiq'"
    }
    local lbls `lbls' 99 "control"
    label define aux_lbl `lbls', replace
    label values aux_rec aux_lbl

    quietly reghdfe `var' ib4.aux_rec, absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons
    eststo rec

    coefplot rec, drop(_cons) vertical baselevels label ///
        title("`var': Recorte puro [-5,5] (entrada)") ///
        yline(0, lpattern(dash) lcolor(black)) ///
        ciopts(recast(rcap)) xlabel(, angle(45)) ///
        name(g_rec, replace) graphregion(color(white))

    * ---- Opción A: pendientes sobre betas (REC)
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
                if `lev'!=99 {
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
        title("Event Study: Entrada al IPSA sobre `var'") ///
        ycommon graphregion(color(white)) iscale(*0.8)
    graph export "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/Graf final/E_`var'_5.png", replace

    *============================================================
    * Opción B: slopes formales en la DATA (TWFE, cluster empresa)
    * PRE = [-5,-2]; POST = [0,5]
    *============================================================
    gen tpre_treat  = cond(treated==1 & !missing(event_time) & inrange(event_time,-5,-2), event_time, 0)
    gen tpost_treat = cond(treated==1 & !missing(event_time) & inrange(event_time,0,5),  event_time, 0)
    gen post_treat  = cond(treated==1 & !missing(event_time) & event_time>=0, 1, 0)

    quietly reghdfe `var' tpre_treat tpost_treat post_treat, ///
        absorb(i.empresa_id i.year) vce(cluster empresa_id) nocons

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

*==== Exportar resultados ====
postclose SLa
use `slopesA_out', clear
order ratio spec b_pre se_pre b_post se_post
export delimited using "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/Graf final/pendientes_sobre_betas.csv", replace

postclose SLb
use `slopesB_out', clear
order ratio beta_pre se_pre beta_post se_post delta se_delta diff_post_pre se_diff
export delimited using "/Users/sebamejias23/Desktop/Invest. GL/Bases finales/Graf final/tendencias_pre_post_B.csv", replace



*Combinar ambos gráficos, de entrada y salida. Lo hice al final en un nuevo código
