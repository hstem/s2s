*********************************************************************
* 			                                        				*
*   		The when, the why, and the how to impute: 			  	*
* A practitioners' guide to survey-to-survey imputation of poverty	*
* 			                                        				*
*********************************************************************  

*** Authors: Paul Corral, Leonardo Lucchetti, Andres Ham, Henry Stemmler, P.F. Lanjouw ***
*** THIS DO-FILE: Imputation ****
*** VERSION: 01/16/2025 ***

cap program drop s2s_method

program s2s_method, rclass
	version 18.0
	syntax [if] [in], [weight(string) DEPVar(string) modelindvars(string) y1(string) y2(string) id(string) imputation(string) impnum(string) knn(string)]		

	// Warning messages and invalid values
	if "`imputation'" ~= "ols_norm" & "`imputation'" ~= "ols_bootstrap" & "`imputation'" ~= "pmm" ///
		& "`imputation'" ~= "heteroskedastic" & "`imputation'" ~= "heteroskedastic_emp" 	  ///
		& "`imputation'" ~= "lasso_pmm" & "`imputation'" ~= "lasso_nonnormal" {
		noi di in red `"Imputation method not specified or invalid value. Please use "ols_norm", "ols_bootstrap", "pmm", "heteroskedastic", "heteroskedastic_emp", "lasso_pmm"."'
		exit 1
	}	
	
	if "`impnum'" == "" & ("`imputation'" == "ols_norm" | "`imputation'" == "ols_bootstrap" | "`imputation'" == "pmm"	///
		| "`imputation'" == "heteroskedastic" | "`imputation'" == "heteroskedastic_emp" 								///
		| "`imputation'" == "lasso_pmm" | "`imputation'" == "lasso_nonnormal") {
		local impnum 20
		noi di in yellow `"No number of imputations defined, default is 20."'
	}
	
	if "`knn'" == "" & ("`imputation'" == "pmm"  | "`imputation'" == "lasso_pmm") {
		local knn 5
		noi di in yellow `"No number of k-nearest-neighbors defined, default is 5."'
	}	

	if "`id'" == "" & ("`imputation'" == "heteroskedastic_emp"  | "`imputation'" == "lasso_pmm") {
		noi di in red `"Unique id needs to be specified for `imputation'."'
		exit 1
	}				

	if "`het'" == "" & ("`imputation'" == "heteroskedastic") {
		noi di in yellow `"For heteroskedastic estimation, default is to use predicted linear fit to model variance."'
	}	
	
	// Imputation
	qui {
		mi set wide
		mi register imputed `depvar'

		// 1) OLS with normal errors
		if "`imputation'" == "ols_norm" {
			set sortseed 94131																			// sortseed and rseed needed to get same results
			mi impute reg `depvar' ${gdp_cons} `modelindvars' [aw=`weight'], ${constant_opt} add(`impnum') rseed(94131)
		}
		
		// 2) OLS with non-normal errors
		if "`imputation'" == "ols_bootstrap" {
			set sortseed 94131																			
			mi impute reg `depvar' ${gdp_cons} `modelindvars', ${constant_opt} add(`impnum') bootstrap rseed(94131)
			noi di in yellow `"No weights are allowed with bootstrap option"."'
		}	

		// 3) PMM
		if "`imputation'" == "pmm" {
			set sortseed 94131																			
			mi impute pmm `depvar' ${gdp_cons} `modelindvars' [aw=`weight'], ${constant_opt} add(`impnum') knn(`knn') rseed(94131)
		}	
		
		// 4) Heteroskedastic with normal errors
		if "`imputation'" == "heteroskedastic" {
			set sortseed 94131	
			reg `depvar' ${gdp_cons} `modelindvars' [aw=`weight'], ${constant_opt}
			predict yhat, xb
			if "`het'" == "" local het yhat
			hetregress `depvar' ${gdp_cons} `modelindvars' [aw=`weight'], ${constant_opt} mle het(`het')
			predict xb, xb
			predict er, sigma
			forval z=1/`impnum' {
				gen double y_`z' = exp(rnormal(xb, er))
			}
			rename (y_*) (_*_`depvar')
			foreach i of varlist _*_`depvar' {
				replace `i' = log(`i')
			}
		}
		
		// 5) Heteroskedastic with empirical errors
		if "`imputation'" == "heteroskedastic_emp" {
			gen new = year == `y2'
			set sortseed 94131																						
			hetmireg `depvar' ${gdp_cons} `modelindvars' [aw=`weight'], sims(`impnum') uniqid(`id') by(new) errdraw(empirical) r seed(94131)
			noi di in yellow `"NOCONS OPTION TO BE ADDED"."'
			rename (yhat_*) (_*_`depvar')		
		}

		// 6) LASSO - PMM
		if "`imputation'" == "lasso_pmm" {
			set sortseed 94131				
			lassopmm `depvar' ${gdp_cons} `modelindvars' [aw=`weight'], add(`impnum') uniqid(`id' year) knn(`knn') seed(94131)
		}	
		
		// 7) LASSO with non-normal
		if "`imputation'" == "lasso_nonnormal" {
			clear mata
			gen new = year == `y2'
			lasso linear `depvar' ${gdp_cons} `modelindvars' [iw=`weight'], ${constant_opt} selection(cv) rseed(94131)
			predict xb, xb
			gen double res = `depvar' - xb
			gen touse = !missing(res)
			putmata e      = res if touse==1
			putmata xb     = xb  if new==1

			//Simulate vector a la EBP
			local the_y
			forval z=1/`impnum' {
				qui:gen double y_`z' = .
				local the_y `the_y'  y_`z'
				local nv = `z'
			}		
			mata: st_view(la_y=.,.,tokens("`the_y'"),"new")
			mata: la_y[.,.] = exp(xb:+_f_sampleepsi(`nv', rows(xb),e))
	
			rename (y_*) (_*_`depvar')
			foreach i of varlist _*_`depvar' {
				replace `i' = log(`i')
			}
			local modelindvars = e(allvars_sel)									// update independent variables
		}	// end of last method
					
		// Create additional variables	
		gen depvar  = "`depvar'"
		gen r2      = "`e(r2_a)'"
		gen indvars = "`modelindvars'"
		foreach i of local modelindvars {
			gen _b_`i' = _b[`i']
		}	
		else
		
	 } // qui

end	

clear mata
mata
	function _f_sampleepsi(real scalar n, real scalar dim, real matrix eps){				  
		sige2 = J(dim,n,0)
		N = rows(eps)
		if (cols(eps)==1) for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),1]
		else              for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),i]
		//for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(rows(eps)*runiform(dim,1)),i]
		return(sige2)	
	}
end
