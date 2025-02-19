{smcl}
{* *! version 1.0.0  26Dec2018}{...}
{cmd:help lassopmm}
{hline}

{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{cmd:lassopmm} {hline 1}} PMM enabled Lasso regression for multiple imputation, 
command makes use of Wilbur Townsend's lassoregress command {p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 17 2}{opt lassopmm}   {depvar} [{indepvars}] {ifin} {weight} {cmd:,} 
uniqid(varlist) add(integer) [{it:options}]

{synoptset 20 tabbed}{...}
{synoptline}
{synopthdr}
{synoptline}
{syntab:{title:Required}}

{synopt:{opt uniqid(varlist)}}varlist of unique identifiers - used to ensure replicability{p_end}
{synopt:{opt add(integer)}}number of simulations to be run for multiple imputation{p_end}

{syntab:{title:Optional}}

{synopt:{opt SORTy}}method which takes a random draw of the observed dependent variable, 
sorts it and matches to the sorted yhat results for the non-observed data. When not specified, the default, pmm is used.{p_end}
{synopt:{opt psu(varlist)}}specifies the variables identifying resampling clusters.  If psu() is specified, the sample drawn during each replication is a bootstrap sample of clusters, and exp must be less than or equal to N_c (the number of clusters identified by the psu() option).  If strata() is also specified, exp must be less than or equal to the number of within-strata clusters.{p_end}
{synopt:{opt strata(varlist)}}  specifies the variables identifying strata.  If strata() is specified,
        bootstrap samples are selected within each stratum, and exp must be less than or equal
        to _N within the defined strata.{p_end}
{synopt:{opt lambda(real)}}penalty placed on larger coefficients — by default found 
by cross-validation {p_end}
{synopt:{opt numfolds(integer)}}number of folds used when cross-validating lambda or 
alpha — default is 10. {p_end}
{synopt:{opt numlambda(integer)}}number of lambda tested when lambda is found by 
cross-validation{p_end}
{synopt:{opt knn(integer)}}number of closest observations to draw result from {p_end}
{synopt:{opt seed(integer)}}Seed to ensure replicability, if not specified it uses 
the current seed's state. {p_end}
{synopt:{opt addvars(varlist numeric)}}Optional and allows the PMM match to bring over other variables into the imputation.{p_end}
{synopt:{opt NOIsily}}Display results of lassoregress{p_end}

{p 2 2 1}
{cmd:aweight}s are allowed; see {help weight}{break} 
{cmd:mi set wide} is allowed; see {help mi set}{break}
{cmd:mi set mlong} is allowed; see {help mi set}{break} 

{title:Example}

{p 6 10 1}.sysuse auto, clear{p_end}

{p 6 10 1}.gen psu = 7 if _n<74{p_end}
{p 6 10 1}.replace psu = 6 if _n<60{p_end}
{p 6 10 1}.replace psu = 5 if _n<50{p_end}
{p 6 10 1}.replace psu = 4 if _n<40{p_end}
{p 6 10 1}.replace psu = 3 if _n<30{p_end}
{p 6 10 1}.replace psu = 2 if _n<20{p_end}
{p 6 10 1}.replace psu = 1 if _n<10{p_end}

{p 6 10 1}.gen _numobs = _n{p_end}
{p 6 10 1}.preserve{p_end}
{p 6 10 1}.	sample 10{p_end}
{p 6 10 1}.	sum price{p_end}
{p 6 10 1}.	replace price = .{p_end}

{p 6 10 1}.tempfile uno{p_end}
{p 6 10 1}.save `uno'{p_end}
{p 6 10 1}.restore{p_end}
{p 6 10 1}.append using `uno', gen(samples){p_end}

{p 6 10 1}.gen _numobs11 = _n{p_end}

{p 6 10 1}//Local with all candidate variables{p_end}
{p 6 10 1}.local _x mpg headroom trunk weight length turn displacement gear_ratio foreign{p_end}
{p 6 10 1}//Local with dependent variable{p_end}
{p 6 10 1}.local _y price {p_end}

{p 6 10 1}.mi set wide{p_end}

{p 6 10 1}.mi register imputed price {p_end}

{p 6 10 1}.lassopmm `_y' `_x' [aw=weight], knn(5) add(5) psu(psu) seed(12388) uniqid(_numobs11){p_end}
{p 6 10 1}.mi estimate: mean price if samples==1 [aw=weight]{p_end}

{p 18 10 1}({stata lassopmm_ex ex1 :click to run}){p_end}

{title:Authors:}

{pstd}
Paul Corral{break}
The World Bank - Poverty and Equity Global Practice {break}
Washington, DC{break}
Corresponding author{break} 
pcorralrodas@worldbank.org{p_end}

{pstd}
Raul Andres Castaneda{break}
The World Bank - Poverty and Equity Global Practice {break}
Washington, DC{break}
acastaneda@worldbank.org{p_end}

{pstd}
Leonardo Ramiro Lucchetti{break}
The World Bank - Poverty and Equity Global Practice {break}
Washington, DC{break}
llucchetti@worldbank.org{p_end}

{title:References}

{pstd}
Wilbur Townsend, 2017. "ELASTICREGRESS: Stata module to perform elastic net regression, lasso regression, ridge regression," Statistical Software Components S458397, Boston College Department of Economics, revised 16 Apr 2018.

{title:Disclaimer}

{pstd}
Any error or omission is the author's responsibility alone.

