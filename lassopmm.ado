*! version 0.2 July 31-2018
*! R. Andres Castaneda - acastaneda@worldbank.org
*! Paul Corral - pcorralrodas@worldbank.org  -- Corresponding author
*! Leonardo Lucchetti - llucchetti@worldbank.org
*! World Bank Group - Poverty and Equity Global Practice 
*! Equity policy lab

cap program drop lassopmm
program define lassopmm, eclass byable(recall)
	version 13, missing
#delimit;
	syntax varlist(min=2 numeric fv) [if] [in] [aw],
		uniqid(varlist)
		add(integer)
	[
		SORTy
		knn(numlist max=1 int)
		lambda(real -1)
		numlambda(integer 100)
		numfolds(integer 10)
		psu(varlist max=1 numeric)
		strata(varlist max=1 numeric)
		seed(integer 12345)
		addvars(varlist numeric) 
		NOIsily
		EPSIlon(real 0.001) 
		postlasso
		putin(varlist numeric)
	];
#delimit cr
set more off

qui{
*===============================================================================
// 1. HOUSE KEEPING
*===============================================================================
	//Check lassoregress
	cap which lassoregress
	if (_rc) net install elasticregress

	//Mlong?
	local mlong `_dta[_mi_style]'
	
	//confirm mi registered
	if ("`_dta[_mi_marker]'"==""){
		dis as error "You must mi set your data before using this command, see help file"	
		error 119
		exit
	}
	
	if ("`putin'"!="" & "`postlasso'"==""){
		dis as error "putin option can only be used with postlasso"
		error 199
		exit
	}
	
	//mark sample for use
	marksample touse
	
	// Mark for the other variables in the estimation
	foreach j of varlist `psu' `uniqid' `putin'{
		replace `touse' = 0 if missing(`j')
	}
		
	tokenize `varlist'
	// Local for dependent variable
	local depvar `1'
	
	// obtain the independent variables
	macro shift 
	local _my_x `*'
	
	tempvar target
	gen `target'= `touse'==0 & missing(`depvar')
	foreach x of local _my_x{
		replace `target' = 0 if `target'==1 & missing(`x')
	}
	
	preserve
		keep if `target'==1
		tempfile target1
		save `target1'
	restore
	
	local isitin `_dta[_mi_ivars]' 
	local isitin : list isitin & depvar
	
	if ("`isitin'"==""){
		dis as error "You must mi register `depvar' as imputed before using this command"	
		error 198
		exit
	}	
	
	if ("`sorty'"!="" & "`knn'"!=""){
		//Only display warning!
		dis as error "Note that KNN option only works when the default option of PMM is used."
	}
	
	//Weights for estimation
	tempvar wi _psu _group
	if missing(`"`weight'"') generate byte `wi' = 1
	else generate `wi' `exp'
	
	//Check unique id
	isid `uniqid' 
	sort `uniqid'
	
*===============================================================================
// 2. Get bootstrap permutation vectors
*===============================================================================
	set seed `seed'	
	tempfile mydata
	save `mydata'
	
*===============================================================================
// 2.5 Ensure that add vars only works with knn
*===============================================================================
/*
	if ("`addvars'"!="") {	
		if ("`sorty'"!=""){
			dis as error "The addvars() option only works when sorty is not specified"
			error 198
			exit 
		}		
	}
*/	
*===============================================================================
// 3. Run lasso and imputation on different bootstrapped data
*===============================================================================
	//Get numsim vectors of Y	
	forval i=1/`add'{
		if (`i'!=1) use `mydata', clear
		
		bsample if `touse'==1, cluster(`psu') strata(`strata')
		append using `target1', force
		
		if ("`noisily'"!="") noi:lassoregress `depvar' `_my_x' if `touse'==1 [aw=`wi'], numlambda(`numlambda') numfolds(`numfolds') lambda(`lambda') epsilon(`epsilon')
		else                 lassoregress `depvar' `_my_x' if `touse'==1 [aw=`wi'], numlambda(`numlambda') numfolds(`numfolds') lambda(`lambda') epsilon(`epsilon')
		tempname _beta BB
		tempvar touse1 touse2
		
		mat `_beta' = e(b)
		
		gen `touse1' = e(sample)
		gen `touse2'= `touse1'==0 & missing(`depvar')
		
		local a=1
		local chosen
		foreach x of local _my_x{
			if (`_beta'[1,`a']!=0){
				local chosen `chosen' `x'
				mat `BB' = nullmat(`BB'),`_beta'[1,`a']
			}
			local ++a
		}
		mat `BB' = `BB',`_beta'[1,`a']
		
			mata: b = st_matrix("`BB'")
			mata: st_view(x=., .,"`chosen'", "`touse1'")
			mata: st_view(w=., .,"`wi'", "`touse1'")

			mata: st_view(y=., .,"`depvar'", "`touse1'")
			if ("`addvars'"!="") mata:st_view(add=.,.,"`addvars'", "`touse1'")			

			mata: st_view(x1=.,.,"`chosen'", "`touse2'")

			if ("`postlasso'"==""){
				mata: yhat1 = quadcross((x ,J(rows(x), 1,1))',b') + (1e-12):*runningsum(J(rows(x),1,1))
				mata: yhat2 = quadcross((x1,J(rows(x1),1,1))',b')
			}
			else{
				if ("`putin'"!=""){
					local isin: list chosen & putin
					if ("`isin'"==""){ //If putin have not been chosen for the model, we add them by force...
						if (`i'==1 & "`noisily'"!="") noi:reg  `depvar' `putin' if `touse'==1 [aw=`wi'], r
						mata:st_view(xputin=.,.,"`putin'","`touse1'")
						mata: x=x,xputin
						mata:st_view(xputin1=.,.,"`putin'","`touse2'")
						mata: x1=x1,xputin1
					}
				}
				mata: b = quadcross(invsym(quadcross((x ,J(rows(x), 1,1)),w,(x ,J(rows(x), 1,1)))),quadcross((x ,J(rows(x), 1,1)),w,y))
				mata: yhat1 = quadcross((x,J(rows(x),1,1))',b) + (1e-12):*runningsum(J(rows(x),1,1))
				mata: yhat2 =  quadcross((x1,J(rows(x1),1,1))',b)
			}
	
		if (`i'==1){
			if ("`addvars'"!="") mata: y = y,add
			if ("`sorty'"=="")   mata: y1 = y[_Mpmm(yhat1, yhat2, `knn'),.]
			else                 noi:mata: y1 = _randomleo(y,yhat2)
		}
		else{
			if ("`addvars'"!="") mata: y = y,add
			if ("`mlong'"=="mlong"){
				if ("`sorty'"=="") mata: y1 = y1 \	y[_Mpmm(yhat1, yhat2, `knn'),.]
				else               noi:mata: y1 = y1 \	_randomleo(y,yhat2)
			}	
			else{	
				if ("`sorty'"=="") mata: y1 = y1 ,	y[_Mpmm(yhat1, yhat2, `knn'),.]
				else               noi:mata: y1 = y1 ,	_randomleo(y,yhat2)			
			}
		}
	} //End of sim loop
*===============================================================================
// 3. Export simulations to mi set data
*===============================================================================
	if ("`mlong'"=="mlong"){
		use `mydata', clear
			tempvar touse2 nn
			gen `touse2' = `touse'==0 & missing(`depvar')
			foreach x of local _my_x{
				replace `touse2' = 0 if missing(`x')
			}
			qui:count 
			local hu0 = r(N)
				preserve
					qui:keep if `touse2'==1
					tempfile _nnio
					qui:save `_nnio'
					qui:count
					local hu=r(N)
				restore
			local a=1
			while (`a'<=`add'){
				qui:append using `_nnio', gen(`nn')
				qui:replace _mi_m=`a' if `nn'==1
				local ++a
				drop `nn'			
			}
			qui: char _dta[_mi_M] `add'
			qui: char _dta[_mi_n] `hu'
			qui: char _dta[_mi_N] `hu0'	
			replace `touse2'=`touse2'==1 & _mi_m>0
			if ("`addvars'"!="") mata: st_store(.,st_varindex(tokens("`depvar' `addvars'")),"`touse2'",y1)
			else mata: st_store(.,st_varindex(tokens("`depvar'")),"`touse2'",y1)				
	}
	else{
		use `mydata', clear
			tempvar touse2 nn
			gen `touse2' = `touse'==0 & missing(`depvar')
			replace _mi_miss = 1 if `touse2'==1
			local a=1
			local toimp
			while (`a'<=`add'){
				gen double _`a'_`depvar' = `depvar' if `touse'==1
				local toimp `toimp' _`a'_`depvar'
				foreach var of local addvars{
					gen double _`a'_`var' = `var' if `touse'==1
					local toimp `toimp' _`a'_`var'
				}

				local ++a				
			}
			qui: char _dta[_mi_M] `add'
			if ("`addvars'"!="") mata: st_store(.,st_varindex(tokens("`toimp'")),"`touse2'",y1)				
			else mata: st_store(.,st_varindex(tokens("`toimp'")),"`touse2'",y1)				

	}	
}
end

*===============================================================================
// Annex: MATA functions
*===============================================================================
mata
	//Function will return an index selection vector for PMM Y
	function _Mpmm(yh1, yh2, knn){
		//Search distance to yh2
		ry2 = rows(yh2)
		ry1 = rows(yh1)
		if (knn>1){
			for(i=1; i<=ry2; i++){
				myy = order(abs(yh1:-yh2[i]),1)[|1\knn|]
				if (i==1) mynn = _f_sampleepsi(1,1,myy)
				else mynn = mynn \ _f_sampleepsi(1,1,myy)
			}
		}
		else{
			for(i=1; i<=ry2; i++){
				if (i==1) mynn = order(abs(yh1:-yh2[i]),1)[1]	
				else mynn = mynn \ order(abs(yh1:-yh2[i]),1)[1]	
			}	
		}
		return(mynn)
	}	
	//Function, selects a random set of observed y, of length rows unobserved. 
	// sorts the set and assigns the value to the sorted xb from unobserved
	function _randomleo(yo, yh2){
		ry = rows(yh2)
		totcol = (cols(yo)+1)
		
		if (cols(yo)>1){		
			//random sample of y
			tosort     = sort((runningsum(J(ry,1,1)),yh2),2)
			tosort2 = sort(_f_sampleepsi22(1, ry, yo),1)
			tosort2 = sort((tosort[.,1],tosort2),1)
			cols(tosort2)
			return(tosort2[.,2..totcol])
		}
		else{
			tosort      = sort((runningsum(J(ry,1,1)),yh2),2)
			tosort[.,2] = sort(_f_sampleepsi(1, ry, yo),1)
			_sort(tosort,1)
		
			return(tosort[.,2])
		}
	}	
	//n is the number of simulations, dim is the number of rows of the output, epsi is the source
	function _f_sampleepsi(real scalar n, real scalar dim, real matrix eps){				  
		sige2 = J(dim,n,0)
		N = rows(eps)
		if (cols(eps)==1) for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),1]
		else              for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),i]
		//for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(rows(eps)*runiform(dim,1)),i]
		return(sige2)	
	}
	function _f_sampleepsi22(real scalar n, real scalar dim, real matrix eps){				  
		
		N = rows(eps)
        sige2=eps[ceil(N*runiform(dim,1)),.]
		//for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(rows(eps)*runiform(dim,1)),i]
		return(sige2)	
	}
	
	//Function for bootstrap permutation vectors
	function _getmyboot(_psu,sim){
		info = panelsetup(_psu,1)
		_index = runningsum(J(rows(_psu),1,1))
		rr = rows(info)
		
		for(i=1; i<=rr; i++){
			m = info[i,1],1 \ info[i,2],1
			r1= rows(_index[|m|])
			if (i==1) _X = _f_sampleepsi(sim,r1,_index[|m|])
			else      _X = _X \ _f_sampleepsi(sim,r1,_index[|m|])
		}
		
		return(_X)
	}	
end



