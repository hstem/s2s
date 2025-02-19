*********************************************************************
* 			                                        				*
*   		The when, the why, and the how to impute: 			  	*
* A practitioners' guide to survey-to-survey imputation of poverty	*
* 			                                        				*
*********************************************************************  

*** Authors: Paul Corral, Leonardo Lucchetti, Andres Ham, Henry Stemmler, P.F. Lanjouw ***
*** THIS DO-FILE: Model diagnostic ****
*** VERSION: 01/16/2025 ***

cap program drop s2s_modeldiagnostic

program s2s_modeldiagnostic, rclass
	version 18.0
	syntax [if] [in], [weight(string) DEPVar(string) modelindvars(string) DIAGNostics(string) VIFvalue(integer 3)]		
		
		
	// Warning messages and invalid values
	if "`diagnostics'" ~= "" & "`diagnostics'" ~= "none" & "`diagnostics'" ~= "all" & "`diagnostics'" ~= "vif"  & "`diagnostics'" ~= "influence" {
		noi di in red `"Model diagnostic option has an invalid value. Please use "none", "all", "vif", "influence", or leave out."'
		exit 1
	}
	if "`diagnostics'" == "" {
			noi di in yellow "Model diagnostic option not specified, no model diagnostic implemented."
	}
		
	// VIF cut-off (default = 3)
	if "`vifvalue'"=="" local vifvalue 3
	
	qui {
		// Model diagnostics 1: Variance inflation factor (Checking for Multicollinearity)
		if "`diagnostics'" == "vif" | "`diagnostics'" == "all" {
			local stop = 0											
			while `stop' == 0 {															// Loop until no vif is > vifvalue
				reg `depvar' `modelindvars' [pw=`weight']
				estat vif															// variance inflation factor
				local var `r(name_1)'
				local vif r(vif_1)
				if `vif' >= `vifvalue' {
					local modelindvars: list modelindvars - var
				}
				else if `vif' < `vifvalue' {
					local stop = 1
				}
			} // stop
		}	// vif
		
		// Model diagnostics 2: Influence analysis
		gen outliers = .
		if "`diagnostics'" == "influence" | "`diagnostics'" == "all" {
			reg `depvar' `modelindvars'
			predict cdist, cooksd															// calculates the Cook´s D influence statistic
			predict rstud, rstudent															// calculates the Studentized (jackknifed) residuals
			reg `depvar' `modelindvars' [aw=`weight']
			predict lev, leverage															// calculates the diagonal elements of the proj matrix
			predict r, resid																// calculates the residuals
			local myN=e(N) 																	// # observations
			local myK=e(rank) 																// rank of k
			// Identified influential / outliers observations
			replace outliers = 1 if abs(rstud)>2 & cdist>4/`myN' & lev>(2*`myK'+2)/`myN' & `depvar' != .
			drop cdist rstud lev r 
		} // influence
	} //qui
	global modelindvars `modelindvars'
end
			