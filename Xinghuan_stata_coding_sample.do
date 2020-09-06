********************************************************************
// code_sample.do | By Xinghuan Luo (xinghuanluo@uchicago.edu)

// This is a code sample extracted from my homework of econometrics class at harris. 

// The questions in this homework are based heavily on the paper Almond et al. 2005. 
// The goal of this assignment is to examine the research question: what is the causal 
// effect of maternal smoking during pregnancy on infant birthweight and other infant health
// outcomes? The data for the problem set is an extract of all births from the 1993 National 
// Natality Detail Files for Pennsylvania. 

// This code file has 4 parts: 
// - Part 1 Initialization and Importing data
// - Part 2 Data Checking
// - Part 3 Estimating Propensity Score 
// - Part 4 Analysis with Propensity Score

// - Basically I first check the missing pattern of data and clean it. Then, I find that the covariates are unbalanced 
//   between smoking and non-smoking group. Also because the dataset is not the result of random experiement, to control 
//   for selection on observables, I use propensity score matching method. Finally, after estimating p-score, I use it 
//   with three different ways, simply including p-score as covariate, weighting both group by p-score and diving whole
//   population by p-score. 

************************Part 1 Initialization and Importing data***********************
clear
cap log close 
set more off

if inlist("`c(username)'", "lenovo" ) {
	cd "D:\OneDrive - The University of Chicago\summer\RA Application\Github\Xinghuan_Stata_Coding_Sample"
}

use "dataset.dta"

// From the codebook, below variables have missing values. 
//I use Little's test (mcartest) to see if the below variables are missing completely at random. 
preserve 

local var_mi_99 cigar alcohol wgain 
local var_mi_other tobacco cigar6 alcohol drink herpes

foreach var of varlist `var_mi_99'{
	replace `var' =. if `var' == 99
}

foreach var of varlist `var_mi_other'{
	replace `var' =. if inlist(1, tobacco==9, cigar6==6, alcohol==9, drink5==5, herpes==8)
}

mcartest tobacco cigar cigar6 alcohol drink drink5 wgain herpes 

// Drop variable with missing values.
restore  
drop if 1 == inlist(1, tobacco==9, cigar==99, cigar6==6, alcohol==9, drink==99, drink5==5, wgain==99, herpes==8) 

************************Part 2 Data Checking***********************

// List a group of predetermiend variables as covariates
local predetermined dmage mrace3 dmeduc dmar dfage dfeduc orfath cntocpop stresfip ormoth nprevist adequacy alcohol drink drink5 preterm pre4000 phyper monpre rectype anemia cardiac lung diabetes herpes chyper disllb isllb10 birmon stresfip pldel3 nlbnl dlivord dtotord totord9 weekday dgestat csex dplural 

// To find if there is any selection bias, I do balance check of all variables between treatment and control group. 
// If the pregnant women smokes, then tobacco = 1(treatment). Otherwise, tobacco = 0(control). 
replace tobacco = 0 if tobacco == 2
balancetable tobacco _all using B_table.tex, longtable pval oneline replace ctitles("No Smoking" "Smoking" "(i)Difference")

// After I control the predermined variables, I simply estimate the impact of smoking on birth weight, one minute apgar score and five minute agpar score.
estimates clear 
local i = 1
foreach var of varlist dbrwt omaps fmaps {
    eststo model_`i': qui reg `var' tobacco `predetermined', robust
    esttab model_`i' using model_`i'.tex, se keep(_cons tobacco ) replace label
    local i = `i' + 1 
}

************************Part 3 Estimating Propensity Score***********************

// To better control selection on observables, I use propensity score matching to find the treatment effect. 
// I use logit specification to generate propensity score 
logit tobacco `predetermined',vce(r)
predict pscore_1

// I only include predetermined covariates who are significant in the last logit regression
local predetermined_sig stresfip rectype adequacy dmage mrace3 dmeduc dmar dfage dfeduc orfath ormoth nprevist alcohol drink5 preterm pre4000 phyper  isllb10 pldel3 dlivord dtotord totord9
logit tobacco `predetermined_sig',vce(r)
predict pscore_2

// To better compare the above two propensity scores, I calculate the correlation between those 2 propensity scores 
// and simulate their density function by kdensity command
corr(pscore_1 pscore_2)
twoway (kdensity pscore_1) || (kdensity pscore_2),  xtitle("Propensity Score") legend(label(1 "Including All") label(2 "Including Signif Only"))

************************Part 4 Analysis with Propensity Score***********************

// First, I simply include ps_score as covariates. 
reg dbrwt tobacco pscore_2

// Second, I use them to reweight the outcomes and estimate the average treatment effect. 
gen ate_weight     = (1/pscore_2)   if tobacco==1
replace ate_weight = 1/(1-pscore_2) if tobacco==0
reg dbrwt tobacco [pweight=ate_weight]

// I use them to reweight the outcomes and estimate the average treatment effect on treated. 
gen att_weight     = 1                   if tobacco==1
replace att_weight = pscore_2/(1-pscore_2) if tobacco==0
reg dbrwt tobacco [pweight=att_weight]

//Third, I make blockings on propensity score. 
//I divide the data  into 100 approximately equally spaced bins based on the estimated propensity score (pscore_2).
sum pscore_2 
local bandwidth = (r(max)-r(min))/100
gen bin_num =.
forval i = 1(1)100{
	qui replace bin_num =`i' if pscore_2 <=`i'*`bandwidth' & bin_num ==.
}
replace bin_num = 100 if bin_num ==.

// I calculate the average treatment effect(b_tobacco) of each bin and calculate 
// the total treatment effect by weighting the b_tobacco of each bin based on the size of each bin
local b_tobacco = 0
foreach i of num 1/100 {
	qui reg dbrwt tobacco if bin_num == `i'
	local b_tobacco = `b_tobacco'+ _b[tobacco]*`e(N)'/_N
}
di `b_tobacco'

// I redo the last part using an indicator for low weight birth (less than 2500 grams) as the outcome. 
// I use different way here. I calculate the size of each bin first. 
gen low_bw = dbrwt < 2500
sum pscore_2
gen psc_bins = autocode(pscore_2, 100, r(min), r(max))
egen bin_size = sum(1), by(psc_bins)

// Then, for each bin, I calculate the mean value of low weight birth for smoking and non-smoking group respectively.  
egen mean_low_1 = mean(low) if tobacco == 1, by(psc_bins)
egen mean_low_0 = mean(low) if tobacco == 0, by(psc_bins)

//After that, for each bin, to calculate the average outcome of smoking and non-smoking group, I divided mean_low_1 and mean_low_0 
//by the size of bin. Finally, I get the average ate_low by weighting the smoking effect of each bin based on their size. 
preserve
collapse(mean) mean_low_0 mean_low_1 bin_size, by(psc_bins)
egen N = sum(bin_size)
egen ate_low = sum((mean_low_1 - mean_low_0)*(bin_size/N))
display ate_low
restore
