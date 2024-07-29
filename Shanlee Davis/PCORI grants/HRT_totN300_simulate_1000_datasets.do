* Oral, oral, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.3583)
generate group = "Oral, oral, oral"
replace id = id + 0.1
	
save temp, replace
	
* Oral, oral, transdermal
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.38253)
generate group = "Oral, oral, transdermal"
replace id = id + 0.2

append using temp

save temp, replace

* Oral, oral, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2399)
generate group = "Oral, oral, OCP"
replace id = id + 0.3

append using temp

save temp, replace

* Oral, transdermal, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.3583)
generate group = "Oral, transdermal, oral"
replace id = id + 0.4

append using temp

save temp, replace

* Oral, transdermal, transdermal
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.40676)
generate group = "Oral, transdermal, transdermal"
replace id = id + 0.5

append using temp

save temp, replace

* Oral, transdermal, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "Oral, transdermal, OCP"
replace id = id + 0.6

append using temp

save temp, replace

* Oral, OCP, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2399)
generate group = "Oral, OCP, oral"
replace id = id + 0.7

append using temp

save temp, replace

* Oral, OCP, transdermal
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "Oral, OCP, transdermal"
replace id = id + 0.8

append using temp

save temp, replace

* Oral, OCP, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.1215)
generate group = "Oral, OCP, OCP"
replace id = id + 0.9

append using temp

save temp, replace

* Transdermal, oral, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.3825)
generate group = "Transdermal, oral, oral"
replace id = id + 0.10

append using temp

save temp, replace

* Transdermal, oral, transdermal
aldsim, totn(12) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.40676)
generate group = "Transdermal, oral, transdermal"
replace id = id + 0.11

append using temp

save temp, replace

* Transdermal, oral, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "Transdermal, oral, OCP"
replace id = id + 0.12

append using temp

save temp, replace

* Transdermal, transdermal, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.40676)
generate group = "Transdermal, transdermal, oral"
replace id = id + 0.13

append using temp

save temp, replace

* Transdermal, transdermal, transdermal
aldsim, totn(12) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.431)
generate group = "Transdermal, transdermal, transdermal"
replace id = id + 0.14

append using temp

save temp, replace

* Transdermal, transdermal, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.28836)
generate group = "Transdermal, transdermal, OCP"
replace id = id + 0.15

append using temp

save temp, replace

* Transdermal, OCP, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "Transdermal, OCP, oral"
replace id = id + 0.16

append using temp

save temp, replace

* Transdermal, OCP, transdermal
aldsim, totn(12) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.28836)
generate group = "Transdermal, OCP, transdermal"
replace id = id + 0.17

append using temp

save temp, replace

* Transdermal, OCP, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.14573)
generate group = "Transdermal, OCP, OCP"
replace id = id + 0.18

append using temp

save temp, replace

* OCP, oral, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2399)
generate group = "OCP, oral, oral"
replace id = id + 0.19

append using temp

save temp, replace

* OCP, oral, transdermal
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "OCP, oral, transdermal"
replace id = id + 0.20

append using temp

save temp, replace

* OCP, oral, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.1215)
generate group = "OCP, oral, OCP"
replace id = id + 0.21

append using temp

save temp, replace

* OCP, transdermal, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "OCP, transdermal, oral"
replace id = id + 0.22

append using temp

save temp, replace

* OCP, transdermal, transdermal
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.28836)
generate group = "OCP, transdermal, transdermal"
replace id = id + 0.23

append using temp

save temp, replace

* OCP, transdermal, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.14573)
generate group = "OCP, transdermal, OCP"
replace id = id + 0.24

append using temp

save temp, replace

* OCP, OCP, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.1215)
generate group = "OCP, OCP, oral"
replace id = id + 0.25

append using temp

save temp, replace

* OCP, OCP, transdermal
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.14573)
generate group = "OCP, OCP, transdermal"
replace id = id + 0.26

append using temp

save temp, replace

* OCP, OCP, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.0031)
generate group = "OCP, OCP, OCP"
replace id = id + 0.27

append using temp

save temp, replace

generate sim = 1

save temp, replace

forvalues i=2/1000{
	
	* Oral, oral, oral
	aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.3583)
	generate group = "Oral, oral, oral"
	replace id = id + 0.1
    generate sim = `i'
	
	save temp1, replace
	
	* Oral, oral, transdermal
	aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.38253)
	generate group = "Oral, oral, transdermal"
	replace id = id + 0.2
	
    generate sim = `i'
	append using temp1
	
	save temp1, replace
	
	* Oral, oral, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2399)
generate group = "Oral, oral, OCP"
replace id = id + 0.3

generate sim = `i'
append using temp1

save temp1, replace

* Oral, transdermal, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.3583)
generate group = "Oral, transdermal, oral"
replace id = id + 0.4

generate sim = `i'
append using temp1

save temp1, replace

* Oral, transdermal, transdermal
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.40676)
generate group = "Oral, transdermal, transdermal"
replace id = id + 0.5

generate sim = `i'
append using temp1

save temp1, replace

* Oral, transdermal, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "Oral, transdermal, OCP"
replace id = id + 0.6

generate sim = `i'
append using temp1

save temp1, replace

* Oral, OCP, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2399)
generate group = "Oral, OCP, oral"
replace id = id + 0.7

generate sim = `i'
append using temp1

save temp1, replace

* Oral, OCP, transdermal
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "Oral, OCP, transdermal"
replace id = id + 0.8

generate sim = `i'
append using temp1

save temp1, replace

* Oral, OCP, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.1215)
generate group = "Oral, OCP, OCP"
replace id = id + 0.9

generate sim = `i'
append using temp1

save temp1, replace

* Transdermal, oral, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.3825)
generate group = "Transdermal, oral, oral"
replace id = id + 0.10

generate sim = `i'
append using temp1

save temp1, replace

* Transdermal, oral, transdermal
aldsim, totn(12) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.40676)
generate group = "Transdermal, oral, transdermal"
replace id = id + 0.11

generate sim = `i'
append using temp1

save temp1, replace

* Transdermal, oral, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "Transdermal, oral, OCP"
replace id = id + 0.12

generate sim = `i'
append using temp1

save temp1, replace

* Transdermal, transdermal, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.40676)
generate group = "Transdermal, transdermal, oral"
replace id = id + 0.13

generate sim = `i'
append using temp1

save temp1, replace

* Transdermal, transdermal, transdermal
aldsim, totn(12) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.431)
generate group = "Transdermal, transdermal, transdermal"
replace id = id + 0.14

generate sim = `i'
append using temp1

save temp1, replace

* Transdermal, transdermal, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.28836)
generate group = "Transdermal, transdermal, OCP"
replace id = id + 0.15

generate sim = `i'
append using temp1

save temp1, replace

* Transdermal, OCP, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "Transdermal, OCP, oral"
replace id = id + 0.16

generate sim = `i'
append using temp1

save temp1, replace

* Transdermal, OCP, transdermal
aldsim, totn(12) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.28836)
generate group = "Transdermal, OCP, transdermal"
replace id = id + 0.17

generate sim = `i'
append using temp1

save temp1, replace

* Transdermal, OCP, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.14573)
generate group = "Transdermal, OCP, OCP"
replace id = id + 0.18

generate sim = `i'
append using temp1

save temp1, replace

* OCP, oral, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2399)
generate group = "OCP, oral, oral"
replace id = id + 0.19

generate sim = `i'
append using temp1

save temp1, replace

* OCP, oral, transdermal
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "OCP, oral, transdermal"
replace id = id + 0.20

generate sim = `i'
append using temp1

save temp1, replace

* OCP, oral, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.1215)
generate group = "OCP, oral, OCP"
replace id = id + 0.21

generate sim = `i'
append using temp1

save temp1, replace

* OCP, transdermal, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.2641)
generate group = "OCP, transdermal, oral"
replace id = id + 0.22

generate sim = `i'
append using temp1

save temp1, replace

* OCP, transdermal, transdermal
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.28836)
generate group = "OCP, transdermal, transdermal"
replace id = id + 0.23

generate sim = `i'
append using temp1

save temp1, replace

* OCP, transdermal, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.14573)
generate group = "OCP, transdermal, OCP"
replace id = id + 0.24

generate sim = `i'
append using temp1

save temp1, replace

* OCP, OCP, oral
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.1215)
generate group = "OCP, OCP, oral"
replace id = id + 0.25

generate sim = `i'
append using temp1

save temp1, replace

* OCP, OCP, transdermal
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.14573)
generate group = "OCP, OCP, transdermal"
replace id = id + 0.26

generate sim = `i'
append using temp1

save temp1, replace

* OCP, OCP, OCP
aldsim, totn(6) cn(6) ci(2) pn(13) pi(1) agelb(14.1) ageub(14.9) intsd(5) slpsd(0.05) slpintcor(-0.1) resid(0.559017) effsize(0.0031)
generate group = "OCP, OCP, OCP"
replace id = id + 0.27

generate sim = `i'
append using temp1

save temp1, replace
	
	append using temp

	save temp, replace
}

sort sim group, stable
outsheet id sim cohort group period age y using "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Davis/PCORI proposals/HRT/ald_sim_totN_300.csv" , comma replace


OLD CODE BELOW


aldsim, totn(24) cn(4) ci(1) pn(3) pi(1) agelb(8) ageub(14) intsd(5) slpsd(1) slpintcor(-0.1) resid(0.559017) effsize(0.5)

generate group = "T1D normal weight"
replace id = id + 0.2

save temp, replace

aldsim, totn(24) cn(4) ci(1) pn(3) pi(1) agelb(8) ageub(14) intsd(5) slpsd(1) slpintcor(-0.1) resid(0.559017) effsize(1.1)

generate group = "T1D overweight"
replace id = id + 0.1

append using temp

generate sim = 1

save temp, replace

forvalues i=2/1000{
    aldsim, totn(40) cn(4) ci(1) pn(3) pi(1) agelb(8) ageub(14) intsd(1) slpsd(1) slpintcor(-0.1) resid(0.559017) effsize(0.5)
	generate group = "lean"
	replace id = id + 0.2
	
	save temp1, replace
	
	aldsim, totn(60) cn(4) ci(1) pn(3) pi(1) agelb(8) ageub(14) intsd(1) slpsd(1) slpintcor(-0.1) resid(0.559017) effsize(1.1)
	generate group = "obese"
	replace id = id + 0.1
	
	append using temp1
	
	generate sim = `i'

	append using temp
	
	save temp, replace
}

sort sim, stable

outsheet id sim cohort group period age y using "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Grants/Puberty and kidney structure and function R01 (PANTHER)/JDRF proposal to add T1D group/ald_sim_24_per_group.csv" , comma replace
