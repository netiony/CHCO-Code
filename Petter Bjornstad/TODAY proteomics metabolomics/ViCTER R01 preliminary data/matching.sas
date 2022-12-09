libname matching 'E:\Petter Bjornstad\TODAY subaward\ViCTER matching';

/* NEED TO ADD TX GROUP */

 /**********************************************************************
 *   PRODUCT:   SAS
 *   VERSION:   9.4
 *   CREATOR:   External File Interface
 *   DATE:      08DEC22
 *   DESC:      Generated SAS Datastep Code
 *   TEMPLATE SOURCE:  (None Specified.)
 ***********************************************************************/
    data WORK.MATCHING    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile 'E:\Petter Bjornstad\TODAY subaward\ViCTER matching\for_matching.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat releaseid $10. ;
       informat AGEBASE best32. ;
       informat sex best32. ;
       informat tanner best32. ;
       format releaseid $10. ;
       format AGEBASE best12. ;
       format sex best12. ;
       format tanner best12. ;
    input
                releaseid  $
                AGEBASE
                sex
                tanner
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;


* OLD CODE;

/* merge two sheets and fix status variable */
proc sort data=alldata1; by record_id; run;
proc sort data=alldata2; by record_id; run;
data alldata;
merge alldata1 alldata2;
by record_id; 
run;
data alldata;
set alldata;
if record_id>=1000 then do;
  status=1;
  diabetes=1;
end;
else do;
  status=2;
  diabetes=0;
end;
run;


data alldata;
set alldata;
age_at_dx=intck('years',dob,pcosdx_date);
if ethnicity=. then ethnicity=2;
bmi_round=round(pcosdx_bmi);
run;

/* create new variable for white/non-white */
data alldata;
set alldata;
if race=5 and ethnicity=2 then white=1;
else white=0;
run;

proc univariate data=alldata;
var  bmi_round ;
histogram bmi_round;
output out=x pctlpre=pv pctlpts=33 66;
run;
proc print data=x; run;

/* need to categorize BMI and age to try and get more matches */
data alldata;
set alldata;
/* these are the original matching critieria but not getting enough */
/* use quartiles of BMI */
/*if age_at_dx in (12,13) then age_cat=1;*/
/*else if age_at_dx in (14,15) then age_cat=2;*/
/*else if age_at_dx in (16,17) then age_cat=3;*/
/*else if age_at_dx  in (18,19,20) then age_cat=4;*/
/*if bmi_round ne . and bmi_round<=33 then bmi_cat=1;*/
/*else if bmi_round>33 and bmi_round<=37 then bmi_cat=2;*/
/*else if bmi_round>37 and bmi_round<=41 then bmi_cat=3;*/
/*else bmi_cat=4;*/
if age_at_dx in (12,13,14) then age_cat=1;
else if age_at_dx in (15,16,17) then age_cat=2;
else if age_at_dx  in (18,19,20) then age_cat=3;
if bmi_round ne . and bmi_round<=34 then bmi_cat=1;
else if bmi_round>34 and bmi_round<=39 then bmi_cat=2;
else if bmi_round>39 then bmi_cat=3;
run;

/* create index variable to match on */
data alldata;
set alldata;
index=  white || age_cat || bmi_cat;
run;
proc print data=alldata;
where record_id in (144,9,297);
run;

/* create separate datasets for cases and controls */
data t2d;
set alldata;
where status= 1;
run;
data nont2d;
set alldata;
where status=2;
run;
proc sort data=nont2d; by index descending person_years; run;
proc print data=nont2d;
by index;
run;


/* check numbers available to match */
proc freq data=t2d;
tables index / out=outa;
run;
proc freq data=nont2d;
tables index / out=outb;
run;
data outa;
set outa;
counta=count;
keep index counta;
run;
data outb;
set outb;
countb=count;
keep index countb;
run;
proc sort data=outa; by index; run;
proc sort data=outb; by index; run;
data temp;
merge outa outb;
by index;
run;
proc sort data=temp ;
by index ;
run;
*ods rtf file='C:\temp\matching categories.rtf' style=journal;
proc print data=temp noobs label; 
label index='Category' counta='Count of patients with T2D' countb='Count of patients without T2D';
where counta ne .;
run;
*ods rtf close;

/* do the matching */
/* 3 controls per case */
proc freq data=t2d noprint;
tables index/list missing out=casecnt (keep=index count rename=(count=casecnt));
run;
proc print data=casecnt; run;
proc freq data=nont2d noprint;
tables index/list missing out=ctrlcnt (keep=index count rename=(count=ctrlcnt));
run;
proc print data=ctrlcnt; run;
DATA ALLCOUNT;
 MERGE CASECNT (IN=A) CTRLCNT (IN=B);
 BY INDEX;
 IF CASECNT > 0;
 IF A AND NOT B THEN CTRLCNT = 0;
 _NSIZE_ = MIN(CASECNT,CTRLCNT);
 IF _NSIZE_ GT 0;
run;
proc print data=allcount; run;
PROC SQL;
CREATE TABLE WORK.ELIGIBLE_CONTROLS AS
SELECT *
FROM nont2d
WHERE INDEX IN (SELECT INDEX FROM ALLCOUNT);
run;
PROC SORT DATA = WORK.ELIGIBLE_CONTROLS;
BY INDEX;
run;
PROC SQL;
CREATE TABLE WORK.ELIGIBLE_CASES AS
SELECT *
FROM t2d
WHERE INDEX IN (SELECT INDEX FROM ALLCOUNT);
run;
PROC SORT DATA = WORK.ELIGIBLE_CASES;
BY INDEX;
run;
PROC SURVEYSELECT DATA = WORK.ELIGIBLE_CASES
 SAMPSIZE = ALLCOUNT
 METHOD = SRS
 SEED=499812
 OUT=WORK.SELECTED_CASES;
 STRATA INDEX;
run;
proc print data=SELECTED_CASES; run;
data allcount; 
set allcount;
if CTRLCNT<3 then _NSIZE_=CTRLCNT;
else _NSIZE_=_NSIZE_*3;
run;
proc print data=allcount; run;
PROC SURVEYSELECT DATA = WORK.ELIGIBLE_CONTROLS
 SAMPSIZE = ALLCOUNT
 METHOD = SRS
 SEED=499812
 OUT=WORK.SELECTED_CONTROLS;
 STRATA INDEX;
run;
proc sort data=selected_controls; by index; run;
proc print data=selected_controls; run;

DATA CC ;
 SET WORK.SELECTED_CONTROLS (IN=A KEEP=record_id INDEX)
 WORK.SELECTED_CASES (IN=B KEEP=record_id INDEX);
IF A THEN CCID = 1; *CONTROLS;
ELSE IF B THEN CCID = 0; *CASES;
run;
PROC SORT DATA= CC;
BY INDEX CCID;
run;

*ods rtf file="C:\temp\output.rtf" style=journal;
proc print data=cc noobs; 
run;
*ods rtf close;
 
/* now need to merge back to alldata and keep the ones selected */
proc sort data=cc; by record_id; run;
proc sort data=alldata; by record_id; run;
data alldata;
merge alldata cc(in=in2);
by record_id;
if in2;
run;
proc sort data=alldata; by index descending person_years; run;

*ods csv file="C:\temp\output.csv" style=journal;
proc print data=alldata;
var index record_id person_years;
run;
*ods csv close;

/* now read in the quad variable */
 /**********************************************************************
 *   PRODUCT:   SAS
 *   VERSION:   9.4
 *   CREATOR:   External File Interface
 *   DATE:      29AUG19
 *   DESC:      Generated SAS Datastep Code
 *   TEMPLATE SOURCE:  (None Specified.)
 ***********************************************************************/
    data WORK.quad    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile 'H:\Endocrinology\Green\PCOS T2D incidence\Revised matching\checking person years.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat index nlnum32. ;
       informat record_id best32. ;
       informat person_years best32. ;
       informat floor_years best32. ;
       informat quad best32. ;
       informat VAR6 $1. ;
       format index nlnum12. ;
       format record_id best12. ;
       format person_years best12. ;
       format floor_years best12. ;
       format quad best12. ;
       format VAR6 $1. ;
    input
                index
                record_id
                person_years
                floor_years
                quad
                VAR6 $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;
data quad;
set quad;
drop var6 person_years floor_years index;
if quad=. then delete;
run;
proc sort data=alldata; by record_id; run;
proc sort data=quad; by record_id; run;
data alldata;
merge quad alldata;
by record_id; 
run;
data alldata;
set alldata;
if quad=. then delete;
run;
proc freq data=alldata; table quad; run;
proc print data=alldata; run;

/* read in Jill's data so we have the predictor information */
data jill;
set jill.working16Aug19;
drop CCID Status quad diabetes;
run;
proc contents; run;
proc sort data=jill; by record_id; run;
data alldata;
merge alldata(in=in1) jill;
by record_id;
if in1;
run;

/* read in additional data pull */
 /**********************************************************************
 *   PRODUCT:   SAS
 *   VERSION:   9.4
 *   CREATOR:   External File Interface
 *   DATE:      06SEP19
 *   DESC:      Generated SAS Datastep Code
 *   TEMPLATE SOURCE:  (None Specified.)
 ***********************************************************************/
    data WORK.more    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile 'H:\Endocrinology\Green\PCOS T2D incidence\Revised matching\more pulls for matching.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat record_id best32. ;
       informat insurance best32. ;
       informat pcosdx_a1c best32. ;
       informat pcosdx_pmh_osa best32. ;
       informat pcosdx_pmh_anxiety best32. ;
       informat pcosdx_pmh_depression best32. ;
       informat pcosdx_fh_pcos best32. ;
       informat pcosdx_fh_t2d best32. ;
       informat pcosdx_fh_obesity best32. ;
       informat pcosdx_fh_anxiety best32. ;
       informat pcosdx_fh_depression best32. ;
       informat pcosdx_fh_osa best32. ;
       informat pcosdx_fh_gdm best32. ;
       informat pcosdx_tg best32. ;
       informat pcosdx_hdl best32. ;
       informat pcosdx_alt best32. ;
       informat pcosdx_FT best32. ;
       format record_id best12. ;
       format insurance best12. ;
       format pcosdx_a1c best12. ;
       format pcosdx_pmh_osa best12. ;
       format pcosdx_pmh_anxiety best12. ;
       format pcosdx_pmh_depression best12. ;
       format pcosdx_fh_pcos best12. ;
       format pcosdx_fh_t2d best12. ;
       format pcosdx_fh_obesity best12. ;
       format pcosdx_fh_anxiety best12. ;
       format pcosdx_fh_depression best12. ;
       format pcosdx_fh_osa best12. ;
       format pcosdx_fh_gdm best12. ;
       format pcosdx_tg best12. ;
       format pcosdx_hdl best12. ;
       format pcosdx_alt best12. ;
       format pcosdx_FT best12. ;
    input
                record_id
                insurance
                pcosdx_a1c
                pcosdx_pmh_osa
                pcosdx_pmh_anxiety
                pcosdx_pmh_depression
                pcosdx_fh_pcos
                pcosdx_fh_t2d
                pcosdx_fh_obesity
                pcosdx_fh_anxiety
                pcosdx_fh_depression
                pcosdx_fh_osa
                pcosdx_fh_gdm
                pcosdx_tg
                pcosdx_hdl
                pcosdx_alt
                pcosdx_FT
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;
proc sort data=more; by record_id; run;
proc sort data=alldata; by record_id; run;
data alldata;
merge alldata(in=in1) more;
by record_id;
if in1;
run;
proc freq data=alldata; tables quad; run;
proc freq data=alldata; tables diabetes; run;

/* for continuous covariates */
%macro contvar (covar);
proc logistic data=alldata;
strata quad;
model diabetes(event="1") = &covar ;
run;
quit;
%mend;
%contvar (pcosdx_a1c);
%contvar (pcosdx_tg);
%contvar (pcosdx_hdl);
%contvar (pcosdx_alt);
%contvar (pcosdx_ft);

/* for categorical covariates */
%macro catvar (covar);
proc logistic data=alldata;
strata quad;
class &covar;
model diabetes(event="1") = &covar ;
run;
quit;
%mend;
%catvar (insurance);
%catvar (pcosdx_pmh_osa);
%catvar (pcosdx_pmh_anxiety);
%catvar (pcosdx_pmh_depression);
%catvar (pcosdx_fh_pcos);
%catvar (pcosdx_fh_t2d);
%catvar (pcosdx_fh_obesity);
%catvar (pcosdx_fh_anxiety);
%catvar (pcosdx_fh_depression);
%catvar (pcosdx_fh_osa);
%catvar (pcosdx_fh_gdm);

proc ttest data=alldata;
var person_years age_at_dx pcosdx_bmiz pcosdx_a1c pcosdx_tg pcosdx_hdl  pcosdx_alt pcosdx_ft ;
class diabetes;
run;

%macro doit(var);
proc freq data=alldata;
table &var;
by diabetes;
run;
%mend;
%doit(white);
%doit(pcosdx_pmh_osa);
%doit(pcosdx_fh_pcos);
%doit(pcosdx_fh_t2d);
%doit(pcosdx_fh_obesity);
%doit(pcosdx_fh_osa);
%doit(pcosdx_fh_gdm);

 proc print data=alldata; 
where diabetes=0 and (pcosdx_pmh_osa=. or pcosdx_fh_pcos=. or pcosdx_fh_t2d=. or pcosdx_fh_obesity =. or pcosdx_fh_osa=. or pcosdx_fh_gdm); 
run;

/* export for Jill */
proc export data=alldata data=alldata outfile='H:\Endocrinology\Green\PCOS T2D incidence\Revised matching\export for jill.csv' dbms=csv;
