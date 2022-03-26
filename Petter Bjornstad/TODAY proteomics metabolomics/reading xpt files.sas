libname data "E:\Petter Bjornstad\TODAY subaward\Clinical data\";

/* comorb */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\COMORB.xpt" data=data.comorb ; 
run;
proc print data=data.comorb;
run;

/* neuro */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\NEURO.xpt" data=data.neuro ; 
run;
proc print data=data.neuro;
run;
