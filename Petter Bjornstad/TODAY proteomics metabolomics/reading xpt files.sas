libname data "E:\Petter Bjornstad\TODAY subaward\Clinical data\";

/* comorb */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\COMORB.xpt" data=data.comorb ; 
run;
proc print data=data.comorb;
run;
proc export data=data.comorb outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\COMORB.csv" dbms=csv replace; run;

/* neuro */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\NEURO.xpt" data=data.neuro ; 
run;
proc print data=data.neuro;
run;
proc export data=data.neuro outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\NEURO.csv" dbms=csv replace; run;

/* cbl */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\CBL.xpt" data=data.cbl ; 
run;
proc export data=data.cbl outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\CBL.csv" dbms=csv replace; run;

/* addcbl */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\ADDCBL.xpt" data=data.addcbl ; 
run;
proc export data=data.addcbl outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\ADDCBL.csv" dbms=csv replace; run;
