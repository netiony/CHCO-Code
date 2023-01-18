options nofmterr;

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
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\ADDCBL.xpt" data=tdata.ADDCBL ; 
run;
proc export data=tdata.ADDCBL outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\ADDCBL.csv" dbms=csv replace; run;

/* TODAY CBL */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\CBL.xpt" data=tdata.CBL ; 
run;
proc export data=tdata.CBL outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\CBL.csv" dbms=csv replace; run;

/* TODAY PE */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\PE.xpt" data=tdata.PE ; 
run;
proc export data=tdata.PE outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\PE.csv" dbms=csv replace; run;

/* TODAY BASELINE */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\BASELINE.xpt" data=tdata.BASELINE ; 
run;
proc export data=tdata.BASELINE outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\BASELINE.csv" dbms=csv replace; run;

/* TODAY PAT */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\PAT.xpt" data=tdata.PAT ; 
run;
proc export data=tdata.PAT outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\PAT.csv" dbms=csv replace; run;

/* AGEBASE */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\AGEBASE.xpt" data=tdata.AGEBASE ; 
run;
proc export data=tdata.AGEBASE outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\AGEBASE.csv" dbms=csv replace; run;

/* 3DPAR */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\3DPAR.xpt" data=tdata.THREEDPAR ; 
run;
proc export data=tdata.THREEDPAR outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\3DPAR.csv" dbms=csv replace; run;

/* ACCEL */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\ACCEL.xpt" data=tdata.ACCEL ; 
run;
proc export data=tdata.ACCEL outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\ACCEL.csv" dbms=csv replace; run;

/* BDI */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\BDI.xpt" data=tdata.BDI ; 
run;
proc export data=tdata.BDI outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\BDI.csv" dbms=csv replace; run;

/* BPE */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\BPE.xpt" data=tdata.BPE ; 
run;
proc export data=tdata.BPE outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\BPE.csv" dbms=csv replace; run;

/* CDI */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\CDI.xpt" data=tdata.CDI ; 
run;
proc export data=tdata.CDI outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\CDI.csv" dbms=csv replace; run;

/* DEXA */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\DEXA.xpt" data=tdata.DEXA ; 
run;
proc export data=tdata.DEXA outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\DEXA.csv" dbms=csv replace; run;

/* ECHO */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\ECHO.xpt" data=tdata.ECHO ; 
run;
proc export data=tdata.ECHO outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\ECHO.csv" dbms=csv replace; run;

/* EDEQ-C */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\EDEQ-C.xpt" data=tdata.EDEQC ; 
run;
proc export data=tdata.EDEQC outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\EDEQ-C.csv" dbms=csv replace; run;

/* FFQ */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\FFQ.xpt" data=tdata.FFQ ; 
run;
proc export data=tdata.FFQ outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\FFQ.csv" dbms=csv replace; run;

/* FUNDUS */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\FUNDUS.xpt" data=tdata.FUNDUS ; 
run;
proc export data=tdata.FUNDUS outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\FUNDUS.csv" dbms=csv replace; run;

/* HUI */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\HUI.xpt" data=tdata.HUI ; 
run;
proc export data=tdata.HUI outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\HUI.csv" dbms=csv replace; run;

/* NEURO */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\NEURO.xpt" data=tdata.NEURO ; 
run;
proc export data=tdata.NEURO outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\NEURO.csv" dbms=csv replace; run;

/* OCT */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\OCT.xpt" data=tdata.OCT ; 
run;
proc export data=tdata.OCT outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\OCT.csv" dbms=csv replace; run;

/* PEDSQLDC */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\PEDSQLDC.xpt" data=tdata.PEDSQLDC ; 
run;
proc export data=tdata.PEDSQLDC outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\PEDSQLDC.csv" dbms=csv replace; run;

/* PEDSQLDT */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\PEDSQLDT.xpt" data=tdata.PEDSQLDT ; 
run;
proc export data=tdata.PEDSQLDT outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\PEDSQLDT.csv" dbms=csv replace; run;

/* PEDSQLGC */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\PEDSQLGC.xpt" data=tdata.PEDSQLGC ; 
run;
proc export data=tdata.PEDSQLGC outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\PEDSQLGC.csv" dbms=csv replace; run;

/* PEDSQLGT */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\PEDSQLGT.xpt" data=tdata.PEDSQLGT ; 
run;
proc export data=tdata.PEDSQLGT outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\PEDSQLGT.csv" dbms=csv replace; run;

/* PRIMOUT */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\PRIMOUT.xpt" data=tdata.PRIMOUT ; 
run;
proc export data=tdata.PRIMOUT outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\PRIMOUT.csv" dbms=csv replace; run;

/* PWC */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\PWC.xpt" data=tdata.PWC ; 
run;
proc export data=tdata.PWC outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\PWC.csv" dbms=csv replace; run;

/* RUN */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\RUN.xpt" data=tdata.RUN ; 
run;
proc export data=tdata.RUN outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\RUN.csv" dbms=csv replace; run;

/* TLP */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\TLP.xpt" data=tdata.TLP ; 
run;
proc export data=tdata.TLP outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TLP.csv" dbms=csv replace; run;

/* TODQUEST */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\TODQUEST.xpt" data=tdata.TODQUEST ; 
run;
proc export data=tdata.TODQUEST outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODQUEST.csv" dbms=csv replace; run;

/* VISIT */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\VISIT.xpt" data=tdata.VISIT ; 
run;
proc export data=tdata.VISIT outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\VISIT.csv" dbms=csv replace; run;

/* YLS */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\TODAY Data for Sub-Study Investigators (Uncollapsed)\TODAY Repository Datasets\YLS.xpt" data=tdata.YLS ; 
run;
proc export data=tdata.YLS outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY\YLS.csv" dbms=csv replace; run;

/**********/
/* TODAY2 */
/**********/

libname t2data "E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2";

/* TODAY2 neuro */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\NEURO.xpt" data=t2data.neuro ; 
run;
proc export data=t2data.neuro outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\NEURO.csv" dbms=csv replace; run;

/* TODAY2 cbl */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\CBL.xpt" data=t2data.cbl ; 
run;
proc export data=t2data.cbl outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\CBL.csv" dbms=csv replace; run;

/* ADDCBL */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\ADDCBL.xpt" data=t2data.ADDCBL ; 
run;
proc export data=t2data.ADDCBL outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\ADDCBL.csv" dbms=csv replace; run;

/* ADDDXA */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\ADDDXA.xpt" data=t2data.ADDDXA ; 
run;
proc export data=t2data.ADDDXA outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\ADDDXA.csv" dbms=csv replace; run;

/* AME */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\AME.xpt" data=t2data.AME ; 
run;
proc export data=t2data.AME outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\AME.csv" dbms=csv replace; run;

/* BDI */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\BDI.xpt" data=t2data.BDI ; 
run;
proc export data=t2data.BDI outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\BDI.csv" dbms=csv replace; run;

/* BERLIN */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\BERLIN.xpt" data=t2data.BERLIN ; 
run;
proc export data=t2data.BERLIN outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\BERLIN.csv" dbms=csv replace; run;

/* COMORB */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\COMORB.xpt" data=t2data.COMORB ; 
run;
proc export data=t2data.COMORB outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\COMORB.csv" dbms=csv replace; run;

/* DDS */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\DDS.xpt" data=t2data.DDS ; 
run;
proc export data=t2data.DDS outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\DDS.csv" dbms=csv replace; run;

/* EATING */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\EATING.xpt" data=t2data.EATING ; 
run;
proc export data=t2data.EATING outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\EATING.csv" dbms=csv replace; run;

/* ECHO */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\ECHO.xpt" data=t2data.ECHO ; 
run;
proc export data=t2data.ECHO outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\ECHO.csv" dbms=csv replace; run;

/* EPWRTH */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\EPWRTH.xpt" data=t2data.EPWRTH ; 
run;
proc export data=t2data.EPWRTH outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\EPWRTH.csv" dbms=csv replace; run;

/* FUNDUS */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\FUNDUS.xpt" data=t2data.FUNDUS ; 
run;
proc export data=t2data.FUNDUS outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\FUNDUS.csv" dbms=csv replace; run;

/* LIFE */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\LIFE.xpt" data=t2data.LIFE ; 
run;
proc export data=t2data.LIFE outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\LIFE.csv" dbms=csv replace; run;

/* LIPO */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\LIPO.xpt" data=t2data.LIPO ; 
run;
proc export data=t2data.LIPO outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\LIPO.csv" dbms=csv replace; run;

/* MEQ */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\MEQ.xpt" data=t2data.MEQ ; 
run;
proc export data=t2data.MEQ outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\MEQ.csv" dbms=csv replace; run;

/* MNI */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\MNI.xpt" data=t2data.MNI ; 
run;
proc export data=t2data.MNI outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\MNI.csv" dbms=csv replace; run;

/* OCT */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\OCT.xpt" data=t2data.OCT ; 
run;
proc export data=t2data.OCT outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\OCT.csv" dbms=csv replace; run;

/* OFFSP */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\OFFSP.xpt" data=t2data.OFFSP ; 
run;
proc export data=t2data.OFFSP outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\OFFSP.csv" dbms=csv replace; run;

/* PEDSQLGA */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\PEDSQLGA.xpt" data=t2data.PEDSQLGA ; 
run;
proc export data=t2data.PEDSQLGA outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\PEDSQLGA.csv" dbms=csv replace; run;

/* PEDSQLGT */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\PEDSQLGT.xpt" data=t2data.PEDSQLGT ; 
run;
proc export data=t2data.PEDSQLGT outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\PEDSQLGT.csv" dbms=csv replace; run;

/* PEDSQLGY */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\PEDSQLGY.xpt" data=t2data.PEDSQLGY ; 
run;
proc export data=t2data.PEDSQLGY outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\PEDSQLGY.csv" dbms=csv replace; run;

/* PEMD */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\PEMD.xpt" data=t2data.PEMD ; 
run;
proc export data=t2data.PEMD outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\PEMD.csv" dbms=csv replace; run;

/* PHQ */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\PHQ.xpt" data=t2data.PHQ ; 
run;
proc export data=t2data.PHQ outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\PHQ.csv" dbms=csv replace; run;

/* PREG */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\PREG.xpt" data=t2data.PREG ; 
run;
proc export data=t2data.PREG outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\PREG.csv" dbms=csv replace; run;

/* PSG */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\PSG.xpt" data=t2data.PSG ; 
run;
proc export data=t2data.PSG outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\PSG.csv" dbms=csv replace; run;

/* PSQI */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\PSQI.xpt" data=t2data.PSQI ; 
run;
proc export data=t2data.PSQI outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\PSQI.csv" dbms=csv replace; run;

/* PWV */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\PWV.xpt" data=t2data.PWV ; 
run;
proc export data=t2data.PWV outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\PWV.csv" dbms=csv replace; run;

/* SPECKLE */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\SPECKLE.xpt" data=t2data.SPECKLE ; 
run;
proc export data=t2data.SPECKLE outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\SPECKLE.csv" dbms=csv replace; run;

/* TME */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\TME.xpt" data=t2data.TME ; 
run;
proc export data=t2data.TME outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TME.csv" dbms=csv replace; run;

/* VISIT */
proc cimport infile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\TODAY2 Data for Sub-Study Investigators (Uncollapsed)\Datasets\VISIT.xpt" data=t2data.VISIT ; 
run;
proc export data=t2data.VISIT outfile="E:\Petter Bjornstad\TODAY subaward\Clinical data\TODAY2\VISIT.csv" dbms=csv replace; run;
