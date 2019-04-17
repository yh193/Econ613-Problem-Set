* ECON 613 Assignment 5 
* Name: Yaqi HU
* Net ID: yh193

* Options and miscellaneity
clear
set more off, perm
set scrollbufsize 2000000

* Set working directory.
cd "C:\Users\yh193\Desktop"
pwd

* HW2-Ex1 Data Creation
set seed 10000
set obs 10000
gen X1 = runiform(1,3)
gen X2 = rgamma(3,2)
gen X3 = rbinomial(1,0.3)
gen eps = rnormal(2,1)
gen Y = 0.5 + 1.2*X1 + (-0.9)*X2 + 0.1*X3 + eps
egen Y_bar = mean(Y)
gen Y_dum = 0
replace Y_dum = 1 if Y > Y_bar

* HW2-Ex2 OLS
* get statistical description about the Y
su Y 
* get the correlation as required
corr Y X1
* conduct OLS
reg Y X1 X2 X3
* conduct OLS by bootstrap 49 and 499 times
reg Y X1 X2 X3, vce(bootstrap, rep(49))
reg Y X1 X2 X3, vce(bootstrap, rep(499))

* HW2-Ex3 Numerical Optimization
* HW2-Ex4 Discrete choice
* conduct probit estimation
probit Y_dum X1 X2 X3
* conduct logit estimation
logit Y_dum X1 X2 X3
* conduct linear estimation
reg Y_dum X1 X2 X3	

* HW2- Ex5 Marginal Effect
* for probit model
quietly probit Y_dum X1 X2 X3
margins, dydx(X1 X2 X3)  	
* for logit model					
quietly logit Y_dum X1 X2 X3
margins, dydx(X1 X2 X3) 


***************************************************************
* HW3-Ex1 
clear

* install package to use eststo commad
ssc install estout,replace 

import delimited "product.csv"

* HW3-Ex1 Data Description
* summarize the statistical description of all product variables
summarize p*

* HW3-Ex2 First Model
rename (p*) (p_1 p_2 p_3 p_4 p_5 p_6 p_7 p_8 p_9 p_10)
reshape long p_, i(v1) j(product)
gen choice_condilogit = 0
replace choice_condilogit = 1 if choice==product
rename p_ price
* estimating the conditional logit model
asclogit choice_condilogit price, case(v1) alternatives(product) base(1)

* HW3-Ex3 Second Model
gen choice_multilogit = 0
replace choice_multilogit = choice if choice==product
* estimating multinomial logit model
mlogit choice_multilogit income fs3_4 fs5 college whtcollar retired, base(1)

* HW3 -Ex4 Calculating Marginal Effect
quietly asclogit choice_condilogit price, case(v1) alternatives(product) base(1)
estat mfx
quietly mlogit choice_multilogit income fs3_4 fs5 college whtcollar retired, base(1)
margins, dydx(income)

* HW3- Ex5 Mixed logit model
* Full mixed
eststo mixed_1: asmixlogit choice_condilogit price, case(v1) casevars(income fs3_4 fs5 college whtcollar retired) alternatives(product) base(1)
drop if choice == 1
drop if product == 1
* deleting some data
eststo mixed_2:asmixlogit choice_condilogit price, case(v1) casevars(income fs3_4 fs5 college whtcollar retired) alternatives(product) base(2)
* Conduct likelihood ratio test
lrtest ( eststo mixed_1) ( eststo mixed_2), stats force

*******************************************************************************
* HW4-Ex1 Data
* import the data
clear
import delimited "Koop-Tobias.csv"

* Declare the panel identifiers - otherwise it won't work.
xtset personid timetrnd
bysort personid: gen t = _n
xtdes
* it is clearly to see the data is unbalanced

* HW4-Ex2 Random Effect
xtreg logwage educ potexper, re

* HW4-Ex3 Fixed Effect Model
* Between Estimator
xtreg logwage educ potexper, be
* Within Estimator
xtreg logwage educ potexper, fe
* Firt Difference Estimator
* get the lag of the variables
gen logwage_F=F.logwage
gen educ_F = F.educ
gen potexper_F=F.potexper
gen logwage_diff=logwage_F-logwage
gen educ_diff = educ_F-educ
gen potexper_diff=potexper_F-potexper

reg logwage_diff educ_diff potexper_diff		

****************************************************
