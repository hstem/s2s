*********************************************************************
* 			                                        				*
*   		The when, the why, and the how to impute: 			  	*
* A practitioners' guide to survey-to-survey imputation of poverty	*
* 			                                        				*
*********************************************************************  

*** Authors: Paul Corral, Leonardo Lucchetti, Andres Ham, Henry Stemmler, P.F. Lanjouw ***
*** THIS DO-FILE: Model selection ****
*** VERSION: 01/16/2025 ***

cap program drop s2s_modelselection

program s2s_modelselection, rclass
	version 18.0
	syntax [if] [in], [weight(string) DEPVar(string) INDVars(string) modelselection(string) stepwisepr(string) lassoopts(string)]		

	// Warning messages and invalid values
	if "`modelselection'" ~= "lasso" & "`modelselection'" ~= "stepwise" & "`modelselection'" ~= "stepwise_cv"  & "`modelselection'" ~= "none" {
		noi di in red `"Model selection method not specified or invalid value. Please use "lasso", "stepwise", "stepwise_cv" or "none"."'
		exit 1
	}	
	
	if "`lassoopts'" == "" & "`modelselection'" == "lasso" {
		local lassoopts alllambdas serule
		noi di in yellow `"No lasso options specified for model selection, default is "alllambdas serule"."'
	}

	if "`stepwisepr'" == "" & "`modelselection'" == "stepwise" {
		local stepwisepr 0.10
		noi di in yellow `"No stepwise probability cut-off defined, default is 0.1."'
	}
	
	qui {
		
		//No model selection 
		if "`modelselection'" == "none"  {
			reg `depvar' `indvars' [aw=`weight']
			local vars : colnames e(b)
			local clean_vars
			* Loop over each variable name and remove those with "o." prefix (no observations)
			foreach var of local vars {
				if substr("`var'", 1, 2) != "o." {
					local clean_vars `clean_vars' `var'
				}
			}
			local cons _cons
			global modelindvars: list clean_vars - cons
		}	
		

		//LASSO model selection 
		if "`modelselection'" == "lasso"  {
			lasso linear `depvar' `indvars' [iw=`weight'], `lassoopts'
			local vars : colnames e(b)
			local cons _cons
			global modelindvars: list vars - cons
		}	
		
		
		//Stepwise model selection (backward)
		if "`modelselection'" == "stepwise"  {
			stepwise, pr(`stepwisepr'): regress `depvar' `indvars' [aw=`weight']
			local vars : colnames e(b)
			local cons _cons
			global modelindvars: list vars - cons
		}
		
		//Stepwise model selection with CV (backward)
		if "`modelselection'" == "stepwise_cv"  {
			// Include 10-fold CV (SWIFT) - first find out p-value that minimizes mse using CV
			capture drop random 
			gen random = runiform()
			xtile fold = random, n(10)
			local i = 1
			local row = 1
			mat A = J(340,3,.)
			local pe = 0.005
			local pr = `pe' + 0.00001
			while `pe' <= 0.1 {
				while `i' <= 10 {
					stepwise, pr(`pr') pe(`pe'): regress `depvar' `indvars' if fold != `i' [aw=`weight']
					local n = `pe' * 1000
					mat A[`row', 1] = `pe' 											//significance level 
					mat A[`row', 2] = `i' 											//fold
					mat A[`row', 3] = (`e(rmse)')^2									// mse
					local i = `i' + 1
					local row = `row' + 1
				}
			local i = 1
			local pe = `pe' + 0.005 
			local pr = `pe' + 0.00001
			}
			svmat A
			bys A1: egen mean_pe = mean(A3)											// what is the lowest average mse
			egen min_pe = min(mean_pe)
			qui sum A1 if mean_pe == min_pe
			local min_pr `r(min)'
			if `min_pr' > 0.1 local min_pr = 0.1
			local min_pe = `min_pr' - 0.00001
			drop min_pe mean_pe
			
			// Now run selection model with determined probability cut-off
			stepwise, pr(`min_pr') pe(`min_pe'): regress `depvar' `indvars' [aw=`weight']
			local vars : colnames e(b)
			local cons _cons
			global modelindvars: list vars - cons
		}
	} // qui
end	