aldsim, totn(100) cn(4) ci(1) pn(3) pi(1) agelb(8) ageub(14) intsd(5) slpsd(1) slpintcor(-0.1) resid(0.559017) effsize(1.0)

generate group = "cSFK"
replace id = id + 0.2

save temp, replace

aldsim, totn(40) cn(4) ci(1) pn(3) pi(1) agelb(8) ageub(14) intsd(5) slpsd(1) slpintcor(-0.1) resid(0.559017) effsize(0.5)

generate group = "Control"
replace id = id + 0.1

append using temp

generate sim = 1

save temp, replace

forvalues i=2/1000{
    aldsim, totn(100) cn(4) ci(1) pn(3) pi(1) agelb(8) ageub(14) intsd(1) slpsd(1) slpintcor(-0.1) resid(0.559017) effsize(1.0)
	generate group = "cSFK"
	replace id = id + 0.2
	
	save temp1, replace
	
	aldsim, totn(40) cn(4) ci(1) pn(3) pi(1) agelb(8) ageub(14) intsd(1) slpsd(1) slpintcor(-0.1) resid(0.559017) effsize(0.5)
	generate group = "Control"
	replace id = id + 0.1
	
	append using temp1
	
	generate sim = `i'

	append using temp
	
	save temp, replace
}

sort sim, stable

outsheet id sim cohort group period age y using "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Grants/LEOPARD power and stats/ald_sim_100_vs_40_per_group.csv" , comma replace
