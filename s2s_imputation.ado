*********************************************************************
* 			                                        				*
*   		The when, the why, and the how to impute: 			  	*
* A practitioners' guide to survey-to-survey imputation of poverty	*
* 			                                        				*
*********************************************************************  

*** Authors: Paul Corral, Leonardo Lucchetti, Andres Ham, Henry Stemmler, P.F. Lanjouw ***
*** THIS DO-FILE: Runs the imputation ****
*** VERSION: 01/16/2025 ***

cap program drop s2s_imputation
program s2s_imputation, rclass
	version 18.0
	syntax [if] [in], [	Country(string) 			///
						weight(string) 				///
						DEPVar(string) 				/// 
						INDVars(string) 			/// 
						BASELINEyears(string)		/// 
						TARGETyears(string)  		///
						id(string) 					///
						AGGregate 					///
						DPath(string)				///
						OPath(string)				///
						constant(string)  			///
						DIAGnostics(string)  		///
						VIFvalue(integer 3)  		///
						region(string)  			///
						modelselection(string) 		///
						stepwisepr(string)  		///
						lassoopts(string)  			///
						imputation(string)  		///
						impnum(string)  			///
						knn(string) 				///
						het(string)] 				
		**************************************************************************************************************************
	* Preparation of options and start loop
	
	//No constant option does not work with hetmireg and lassopmm
	if "`constant'" == "growth" & ("`imputation'" == "heteroskedastic_emp"  | ///
								   "`imputation'" == "lassp_pmm"  			) {
		noi di in red `"Using GDP growth instead of regular constant term is not possible with `imputation'. For this mehtod, please use constant(constant) option, or leave the constant option out."'
		exit 1						   	
	}	 
	
	//Weights
	if "`weight'" == ""  {
		noi di in red `"No weight specified, equal weights assumed."'
		gen weight = 1
	}
	
	//Years
	local years = `:word count `baselineyears''

	// Loop through all years
	forval y = 1 / `years' {
		local y1 = `: word `y' of `baselineyears''
		local y2 = `: word `y' of `targetyears''

	**************************************************************************************************************************
		* Preparation of data
		use "`dpath'/peru_inc_`y1'.dta", clear
		// Generate poverty lines from baseline data
		pctile myp = `depvar' [aw=`weight'], nq(100)							// Percentiles
		forval z = 5(5)95 {
			local pline_`z' = myp[`z']
		}
		if `y1' == `y2'	replace year = `y1'1 if year == `y1'					// For same year imputation, have to distinguish years 
		
		// Add target period
		qui append using "`dpath'/peru_inc_`y2'.dta"
		forval z = 5(5)95 {
			gen _pline`z' = `pline_`z''
		}

		// Dependent variable to impute
		foreach y of local depvar {
			clonevar `y'_real = `y'
			replace  `y' = . if year == `y2'
		}


	**************************************************************************************************************************		
		* 1) Preparation of constant term
		s2s_constant, c(`country') constant(`constant') y1(`y1') y2(`y2')
		
	**************************************************************************************************************************		
		* Save
		if `y1' == `y2' local y1 = `y1'1										// For same year imputation, have to distinguish years 
		tempfile orig_`y1'_`y2' imp`y2' reg_lab
		save `orig_`y1'_`y2'', replace
			
	**************************************************************************************************************************
		* 2) Region-based or national model
		if "`region'" ~= "" & "`region'" ~= "national" {
			qui sum `region', d
			// Save labels
			local reglbl: value label region
			label save `reglbl' using `reg_lab'	
			qui levelsof region, local(region_list)
		}
		else if "`region'" == "" | "`region'" == "national" {
			local region_list = 1
		}
		
		local num_regs = `:word count `region_list''
		di `num_regs'
		
		clear
		save `imp`y2'', replace emptyok
		
		// open loop
		foreach i of local region_list {												// loop through regions (only 1 loop if national)
			use `orig_`y1'_`y2''
			di `num_regs'
			if (`num_regs' > 1) {
				keep if region == `i'													// keep only region
			}
			
	**************************************************************************************************************************
		* 3) Model selection
			s2s_modelselection, weight(`weight') depvar(`depvar') indvars(`indvars') modelselection(`modelselection') 				///
								stepwisepr(`stepwisepr') lassoopts(`lassoopts')

	**************************************************************************************************************************
		* 4) Model diagnostic
			s2s_modeldiagnostic, weight(`weight') depvar(`depvar') modelindvars(${modelindvars}) diagnostics(`diagnostics')			///		// indvars are now modelindvars (from model selection)
								vif(`vifvalue')
							 
	**************************************************************************************************************************
		* 5) Imputation
			s2s_method, weight(`weight') depvar(`depvar') modelindvars(${modelindvars}) y1(`y1') y2(`y2') ///
								id(`id') imputation(`imputation') impnum(`impnum') knn(`knn')
				
			// Loop through all regions	
			qui append using `imp`y2''
			save `imp`y2'', replace
				
		} // region-specifc 
		
	**************************************************************************************************************************
		* 6) Aggregation
			// Imputed values for y1 are just real values, set to missing
		foreach var of varlist _*_`depvar' {
			qui replace `var' = . if year == `y1'
		}
		
		if "`aggregate'" ~= "" {
			sp_groupfunction [aw=`weight'], poverty(*_*_`depvar' `depvar'_real) povertyline(_pline*) by(year)
			gen 	sample = "real" if variable == "`depvar'_real"
			replace sample = "imputed" if sample != "real"
			keep if measure == "fgt0"
			
			groupfunction, mean(value) by(measure reference sample year)
			reshape wide value, i(reference sample) j(year)
			reshape wide value*, i(reference) j(sample, string)
			drop measure value`y1'imputed
			rename (value`y2'imputed value`y1'real value`y2'real) (imputed_povrate_`y2' real_povrate_`y1' real_povrate_`y2')
		}
		
		else if "`aggregate'" == "" {
			cap confirm variable `region'
			if !_rc { 
				keep year id depvar r2 indvars ${modelindvars} `region' *`depvar'* `depvar' `weight' _b* _pline* 
				order year id `region' `depvar' `depvar'_real `weight' depvar r2 indvars 
			}
			else {
				keep year id depvar r2 indvars ${modelindvars} *`depvar'* `depvar' `weight' _b_* _pline* 
				order year id `depvar' `depvar'_real `weight' depvar r2 indvars 
			}
		}
		
	**************************************************************************************************************************
		* Cleaning and save		
		qui {
		gen model_selection = "`modelselection'"
		gen method 			= "`imputation'"
		gen constant		= "`constant'"
		replace constant	= "constant" if "`constant'" == ""
		gen sample			= "region"	 if "`region'" ~= ""
		replace sample		= "national" if "`region'" == "" | "`region'" == "national" 
		gen diagnostics		= "`diagnostics'"
		
		// Shorten some expressions 
		local imp_n		"`imputation'"
		if "`imputation'" == "heteroskedastic" 		local imp_n		"het" 
		if "`imputation'" == "heteroskedastic_emp"  local imp_n		"het_emp"
		
		local diag_n	"`diagnostics'"
		if "`diagnostics'" == "all"					local diag_n	"allmd" 
		if "`diagnostics'" == ""					local diag_n	"nomd" 
		
		local const_n	"`constant'"
		if "`constant'" == ""						local const_n	"constant"  
		
		local reg_n		"`region'"
		if "`region'" == ""							local reg_n		"national"
		}	
		
		if "`aggregate'" ~= "" 						local agg		"agg"
		else 										local agg		"full"
		
		// Final steps
		if "`y1'" == "`y2'1" & "`aggregate'" != "" drop real_povrate_`y1'								// same year imputation, drop real pov rate of duplicate year
		if substr("`y1'",1,4) == "`y2'" local y1 = substr("`y1'",1,4)			// same year imputation, return correct y1
		save "`opath'/`country'_`y1'_`y2'_`modelselection'_`imp_n'_`const_n'_`diag_n'_`reg_n'_`agg'.dta", replace		
							 	
	} // years
	
end
