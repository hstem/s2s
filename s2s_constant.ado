*********************************************************************
* 			                                        				*
*   		The when, the why, and the how to impute: 			  	*
* A practitioners' guide to survey-to-survey imputation of poverty	*
* 			                                        				*
*********************************************************************  

*** Authors: Paul Corral, Leonardo Lucchetti, Andres Ham, Henry Stemmler, P.F. Lanjouw ***
*** THIS DO-FILE: Prepares constant term for model ****
*** VERSION: 01/16/2025 ***

cap program drop s2s_constant

program s2s_constant, rclass
	version 18.0
	syntax [if] [in], [Country(string) constant(string) y1(string) y2(string)]

	// Warning messages and invalid values
	if "`constant'" ~= "" & "`constant'" ~= "constant" & "`constant'" ~= "growth" {
		noi di in red `"Constant option has an invalid value. Please use "constant", "growth", or leave out."'
		exit 1
	}

	if "`constant'" == "" {
			noi di in yellow "Constant option not specified, regular constant used."
	}

	// GDP growth
	qui {
		if "`constant'" == "growth" {
			// Fetch GDP pc
			preserve
				pip, country(`country') year(`y1' `y2') clear
				sum gdp if year == `y1'
				local gdp`y1' = `r(mean)'
				sum gdp if year == `y2'
				local gdp`y2' = `r(mean)'
				local gdpD`y2'`y1' = `gdp`y2''/`gdp`y1''
			restore
			// generate variable with gdp growth
			gen   	gdp_cons = 1 if year==`y1'
			replace gdp_cons = `gdpD`y2'`y1'' if year==`y2'
			global	gdp_cons 	"gdp_cons"
			// nocons option
			global 	constant_opt "nocons"
		}		
		else {
			// empty options to use normal constant
			global constant_opt ""
			global gdp_cons 	""
		}
	} // qui
	
end
	