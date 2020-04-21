
******************************************************************************
*Beata Luczywek
*ECON 220C: Econometrics C
*Problem Set 1
*04/18/2020
*************************

* Purpose: Practice working with panel data, analyze the effect of 'shall' laws
* on violent crime.
clear all

* Set working directory.
local maindir "C:\Users\beata\Documents\_UCSD\Spring 2020\Metrics C\Pset 1\Q5"
local dtadir "`maindir'\dta"

cd "`maindir'"

* Set the log file.
cap log close
log using log_metricsC_pset1_q5.log, replace

pwd /* Print location of working directory to log. */

*Import the handgun file.
use "`dtadir'\handguns"

****************************************************************************** 

// Part I //
gen log_vio = log(vio)
gen log_mur = log(mur)
gen log_rob = log(rob)

sum vio log_vio mur log_mur rob log_rob

reg log_vio shall, r
estimates store m1, title(ln(vio))

reg log_mur shall, r
estimates store m2, title(ln(mur))

reg log_rob shall, r
estimates store m3, title(ln(rob))


label variable shall "shall law"

estout m1 m2 m3, cells(b(star fmt(3)) se(par fmt(2)))   ///
   legend label varlabels(_cons constant)       	    ///
   stats(r2, fmt(3) label(R-sqr))

******************************************************************************
   
// Part II //
reg log_vio shall incarc_rate density pop pm1029 avginc, r
estimates store m1, title(ln(vio))

reg log_mur shall incarc_rate density pop pm1029 avginc, r
estimates store m2, title(ln(mur))

reg log_rob shall incarc_rate density pop pm1029 avginc, r
estimates store m3, title(ln(rob))

estout m1 m2 m3, cells(b(star fmt(3)) se(par fmt(2)))   ///
   legend label varlabels(_cons constant)       	    ///
   stats(r2, fmt(3) label(R-sqr))						///
   drop(incarc_rate density pop pm1029 avginc) 

******************************************************************************

// Part IV //
tab state, gen(statedummy)
tab year, gen(yeardummy)


*Model 2
reg log_vio shall incarc_rate density pop pm1029 avginc statedummy*, r
testparm statedummy*

reg log_mur shall incarc_rate density pop pm1029 avginc statedummy*, r
testparm statedummy*

reg log_rob shall incarc_rate density pop pm1029 avginc statedummy*, r
testparm statedummy*


*Model 3
reg log_vio shall incarc_rate density pop pm1029 avginc year* statedummy*, r
testparm statedummy*
testparm year*

reg log_mur shall incarc_rate density pop pm1029 avginc year* statedummy*, r
testparm statedummy*
testparm year*

reg log_rob shall incarc_rate density pop pm1029 avginc year* statedummy*, r
testparm statedummy*
testparm year*


*Model 4
reg log_vio shall incarc_rate density pop pm1029 avginc statedummy* year*, cluster(state) r

reg log_mur shall incarc_rate density pop pm1029 avginc statedummy* year*, cluster(state) r

reg log_rob shall incarc_rate density pop pm1029 avginc statedummy* year*, cluster(state) r

******************************************************************************

*Close log
log close