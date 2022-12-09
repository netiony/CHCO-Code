libname matching 'E:\Petter Bjornstad\TODAY subaward\ViCTER matching';

 /**********************************************************************
 *   PRODUCT:   SAS
 *   VERSION:   9.4
 *   CREATOR:   External File Interface
 *   DATE:      08DEC22
 *   DESC:      Generated SAS Datastep Code
 *   TEMPLATE SOURCE:  (None Specified.)
 ***********************************************************************/
    data WORK.matching    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile 'E:\Petter Bjornstad\TODAY subaward\ViCTER matching\for_matching.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat releaseid $10. ;
       informat HTN0 best32. ;
       informat HTN best32. ;
       informat DAYSTOHTN best32. ;
       informat LDLDLP0 best32. ;
       informat LDLDLP best32. ;
       informat DAYSTOLDL best32. ;
       informat TGDLP0 best32. ;
       informat TGDLP best32. ;
       informat DAYSTOTG best32. ;
       informat ANYDLP0 best32. ;
       informat ANYDLP best32. ;
       informat DAYSTOANYDLP best32. ;
       informat MIC0 best32. ;
       informat MIC best32. ;
       informat DAYSTOMIC best32. ;
       informat MAC0 best32. ;
       informat MAC best32. ;
       informat DAYSTOMAC best32. ;
       informat NEPHRO0 best32. ;
       informat NEPHRO best32. ;
       informat DAYSTONEPHRO best32. ;
       informat HYP0 $2. ;
       informat HYP $2. ;
       informat DAYSTOHYP $4. ;
       informat RAPID0 $2. ;
       informat RAPID $2. ;
       informat DAYSTORAPID $4. ;
       informat DNE0 best32. ;
       informat DNE best32. ;
       informat DAYSTODNE best32. ;
       informat FILAM0 best32. ;
       informat FILAM best32. ;
       informat DAYSTOFILAM best32. ;
       informat NEURO0 best32. ;
       informat NEURO best32. ;
       informat DAYSTONEURO best32. ;
       informat RETINO $2. ;
       informat DAYSTORETINO $4. ;
       informat MVD0 best32. ;
       informat MVD best32. ;
       informat DAYSTOMVD best32. ;
       informat NUMMVD best32. ;
       informat GLYC best32. ;
       informat DAYSTOGLYC best32. ;
       informat AGEBASE best32. ;
       informat days best32. ;
       informat sex best32. ;
       informat age best32. ;
       informat dxtime best32. ;
       informat birthwt $9. ;
       informat born $2. ;
       informat bf $2. ;
       informat bftime $2. ;
       informat bfstop $2. ;
       informat formula $2. ;
       informat formst $2. ;
       informat milkst $2. ;
       informat foodst $2. ;
       informat callhcp $2. ;
       informat outpt $2. ;
       informat uc $2. ;
       informat er $2. ;
       informat overnt $2. ;
       informat losts $3. ;
       informat lostw $2. ;
       informat losth $2. ;
       informat school $2. ;
       informat grade $2. ;
       informat work $2. ;
       informat hrswk $2. ;
       informat health $2. ;
       informat medicare $2. ;
       informat private $2. ;
       informat diabmed $2. ;
       informat syringe $2. ;
       informat monitor $2. ;
       informat period $2. ;
       informat perage $2. ;
       informat period6 $2. ;
       informat houseedu $2. ;
       informat housetot $2. ;
       informat lives $2. ;
       informat houseinc $2. ;
       informat momhist $2. ;
       informat momage $2. ;
       informat momwght $3. ;
       informat momhght $2. ;
       informat mompren $2. ;
       informat carest $2. ;
       informat momgdm $2. ;
       informat momgdma $2. ;
       informat mgdmoth $2. ;
       informat momdbnow $2. ;
       informat momhfat $2. ;
       informat momhbp $2. ;
       informat dadhist $2. ;
       informat dadage $2. ;
       informat dadwght $3. ;
       informat dadhght $2. ;
       informat daddbnow $2. ;
       informat dadhfat $2. ;
       informat dadhbp $2. ;
       informat fulldiab $2. ;
       informat halfdiab $2. ;
       informat grandiab $2. ;
       informat race best32. ;
       format releaseid $10. ;
       format HTN0 best12. ;
       format HTN best12. ;
       format DAYSTOHTN best12. ;
       format LDLDLP0 best12. ;
       format LDLDLP best12. ;
       format DAYSTOLDL best12. ;
       format TGDLP0 best12. ;
       format TGDLP best12. ;
       format DAYSTOTG best12. ;
       format ANYDLP0 best12. ;
       format ANYDLP best12. ;
       format DAYSTOANYDLP best12. ;
       format MIC0 best12. ;
       format MIC best12. ;
       format DAYSTOMIC best12. ;
       format MAC0 best12. ;
       format MAC best12. ;
       format DAYSTOMAC best12. ;
       format NEPHRO0 best12. ;
       format NEPHRO best12. ;
       format DAYSTONEPHRO best12. ;
       format HYP0 $2. ;
       format HYP $2. ;
       format DAYSTOHYP $4. ;
       format RAPID0 $2. ;
       format RAPID $2. ;
       format DAYSTORAPID $4. ;
       format DNE0 best12. ;
       format DNE best12. ;
       format DAYSTODNE best12. ;
       format FILAM0 best12. ;
       format FILAM best12. ;
       format DAYSTOFILAM best12. ;
       format NEURO0 best12. ;
       format NEURO best12. ;
       format DAYSTONEURO best12. ;
       format RETINO $2. ;
       format DAYSTORETINO $4. ;
       format MVD0 best12. ;
       format MVD best12. ;
       format DAYSTOMVD best12. ;
       format NUMMVD best12. ;
       format GLYC best12. ;
       format DAYSTOGLYC best12. ;
       format AGEBASE best12. ;
       format days best12. ;
       format sex best12. ;
       format age best12. ;
       format dxtime best12. ;
       format birthwt $9. ;
       format born $2. ;
       format bf $2. ;
       format bftime $2. ;
       format bfstop $2. ;
       format formula $2. ;
       format formst $2. ;
       format milkst $2. ;
       format foodst $2. ;
       format callhcp $2. ;
       format outpt $2. ;
       format uc $2. ;
       format er $2. ;
       format overnt $2. ;
       format losts $3. ;
       format lostw $2. ;
       format losth $2. ;
       format school $2. ;
       format grade $2. ;
       format work $2. ;
       format hrswk $2. ;
       format health $2. ;
       format medicare $2. ;
       format private $2. ;
       format diabmed $2. ;
       format syringe $2. ;
       format monitor $2. ;
       format period $2. ;
       format perage $2. ;
       format period6 $2. ;
       format houseedu $2. ;
       format housetot $2. ;
       format lives $2. ;
       format houseinc $2. ;
       format momhist $2. ;
       format momage $2. ;
       format momwght $3. ;
       format momhght $2. ;
       format mompren $2. ;
       format carest $2. ;
       format momgdm $2. ;
       format momgdma $2. ;
       format mgdmoth $2. ;
       format momdbnow $2. ;
       format momhfat $2. ;
       format momhbp $2. ;
       format dadhist $2. ;
       format dadage $2. ;
       format dadwght $3. ;
       format dadhght $2. ;
       format daddbnow $2. ;
       format dadhfat $2. ;
       format dadhbp $2. ;
       format fulldiab $2. ;
       format halfdiab $2. ;
       format grandiab $2. ;
       format race best12. ;
    input
                releaseid  $
                HTN0
                HTN
                DAYSTOHTN
                LDLDLP0
                LDLDLP
                DAYSTOLDL
                TGDLP0
                TGDLP
                DAYSTOTG
                ANYDLP0
                ANYDLP
                DAYSTOANYDLP
                MIC0
                MIC
                DAYSTOMIC
                MAC0
                MAC
                DAYSTOMAC
                NEPHRO0
                NEPHRO
                DAYSTONEPHRO
                HYP0  $
                HYP  $
                DAYSTOHYP  $
                RAPID0  $
                RAPID  $
                DAYSTORAPID  $
                DNE0
                DNE
                DAYSTODNE
                FILAM0
                FILAM
                DAYSTOFILAM
                NEURO0
                NEURO
                DAYSTONEURO
                RETINO  $
                DAYSTORETINO  $
                MVD0
                MVD
                DAYSTOMVD
                NUMMVD
                GLYC
                DAYSTOGLYC
                AGEBASE
                days
                sex
                age
                dxtime
                birthwt  $
                born  $
                bf  $
                bftime  $
                bfstop  $
                formula  $
                formst  $
                milkst  $
                foodst  $
                callhcp  $
                outpt  $
                uc  $
                er  $
                overnt  $
                losts  $
                lostw  $
                losth  $
                school  $
                grade  $
                work  $
                hrswk  $
                health  $
                medicare  $
                private  $
                diabmed  $
                syringe  $
                monitor  $
                period  $
                perage  $
                period6  $
                houseedu  $
                housetot  $
                lives  $
                houseinc  $
                momhist  $
                momage  $
                momwght  $
                momhght  $
                mompren  $
                carest  $
                momgdm  $
                momgdma  $
                mgdmoth  $
                momdbnow  $
                momhfat  $
                momhbp  $
                dadhist  $
                dadage  $
                dadwght  $
                dadhght  $
                daddbnow  $
                dadhfat  $
                dadhbp  $
                fulldiab  $
                halfdiab  $
                grandiab  $
                race
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* need to add Tanner and tx *./


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
