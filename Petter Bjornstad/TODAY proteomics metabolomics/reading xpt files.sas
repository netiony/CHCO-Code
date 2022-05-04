libname alldata "E:\Petter Bjornstad\TODAY subaward\Clinical data\";

/* comorb */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\COMORB.xpt" data=alldata.comorb ; 
run;
proc export data=alldata.comorb outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\COMORB.csv" dbms=csv replace; run;

/**********/
/* TODAY  */
/**********/

libname tdata "E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY";

/* TODAY addcbl */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\ADDCBL.xpt" data=tdata.ADDCBL ; 
run;
proc export data=tdata.ADDCBL outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\ADDCBL.csv" dbms=csv replace; run;

/* TODAY CBL */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sites\TODAY Repository Datasets\CBL.xpt" data=tdata.CBL ; 
run;
proc export data=tdata.CBL outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\CBL.csv" dbms=csv replace; run;

/* TODAY PE */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sites\TODAY Repository Datasets\PE.xpt" data=tdata.PE ; 
run;
proc export data=tdata.PE outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\PE.csv" dbms=csv replace; run;

/* TODAY BASELINE */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sites\TODAY Repository Datasets\BASELINE.xpt" data=tdata.BASELINE ; 
run;
proc export data=tdata.BASELINE outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\BASELINE.csv" dbms=csv replace; run;

/* TODAY PAT */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sites\TODAY Repository Datasets\PAT.xpt" data=tdata.PAT ; 
run;
proc export data=tdata.PAT outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\PAT.csv" dbms=csv replace; run;

/**********/
/* TODAY2 */
/**********/

libname t2data "E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2";

/* TODAY2 neuro */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\NEURO.xpt" data=t2data.neuro ; 
run;
proc export data=t2data.neuro outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\NEURO.csv" dbms=csv replace; run;

/* TODAY2 cbl */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\CBL.xpt" data=t2data.cbl ; 
run;
proc export data=t2data.cbl outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\CBL.csv" dbms=csv replace; run;
