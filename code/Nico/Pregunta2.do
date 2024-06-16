clear all

set seed 12345
set obs 1000

* Configuración de parámetros del modelo
scalar beta0 = 1
scalar beta1 = 2
scalar beta2 = 5
scalar alpha0 = -4
scalar alpha2 = 3

* Valores de alpha1 a probar
local alpha1_values "0.1 0.5 1 5 10"

* Preparación para almacenar resultados
matrix results_beta1_IV = J(`=wordcount("`alpha1_values'")', 1, .)
matrix results_beta0_IV = J(`=wordcount("`alpha1_values'")', 1, .)
matrix results_IV_std_error = J(`=wordcount("`alpha1_values'")', 1, .)

* Generación de datos comunes a todos los loops

gen W = rnormal(2, 1)
gen nu = runiform()
gen Z = (nu < 0.8)
gen e = rnormal(0, 1)
gen u = rnormal(0, 1)

* Loop sobre diferentes valores de alpha1
foreach alpha1 of local alpha1_values {
    * Generación de X y Y específicos para cada iteración
    gen X = alpha0 + `alpha1' * Z + alpha2 * W + u
    gen Y = beta0 + beta1 * X + beta2 * W + e

    * Estimación por Variables Instrumentales (IV)
    ivregress 2sls Y (X = Z W)

    * Almacenar resultados
    matrix results_beta1_IV[_n, 1] = _b[X]
    matrix results_beta0_IV[_n, 1] = _b[_cons]
    matrix results_IV_std_error[_n, 1] = _se[X]

    * Limpiar variables generadas para la siguiente iteración
    drop X Y
}








* Comprobando los datos:

cd "C:\Users\nicol\OneDrive\Documentos\GitHub\metricswmatlab\code\Nico"

import delimited "datos_alpha_0.1.csv", clear
rename v1 Y_i
rename v2 X_i
rename v3 Z_i
rename v4 W_i

local alpha1_values "0.1 0.5 1.0 5.0 10.0"
foreach alpha1 of local alpha1_values {
	
    import delimited "datos_alpha_`alpha1'.csv", clear 
	rename v1 Y_i
	rename v2 X_i
	rename v3 Z_i
	rename v4 W_i
	
	di "Para un alpha1 = `alpha1'"
	reg Y_i X_i					// Pregunta 1
	reg Y_i X_i, robust			// Pregunta 1
	
}
	

local alpha1_values "0.1 0.5 1.0 5.0 10.0"
foreach alpha1 of local alpha1_values {
	
    import delimited "datos_alpha_`alpha1'.csv", clear 
	rename v1 Y_i
	rename v2 X_i
	rename v3 Z_i
	rename v4 W_i
	
	di "Para un alpha1 = `alpha1'"
	ivregress 2sls Y_i (X_i = Z_i), first
}
	
