/* USER SETUP */
%let proj = YOUR_PROJECT_FOLDER_PATH;  
/* Example (Windows): C:\Users\Reshika\seer_project
   Example (Mac): /Users/reshika/seer_project */

/* 1) IMPORT DATA */
proc import datafile="&proj./data_raw/SEER.csv"
    out=work.sdata
    dbms=csv
    replace;
    guessingrows=max;
run;
proc contents data=work.sdata order=varnum;
run;

/* 2) CREATE VARIABLES */

/* Creating 4 categories for EP receptors + event indicator */
data work.data1;
set work.sdata;

if Estrogen_Status = "Positive" and Progesterone_Status = "Positive" then EP = 1;
else if Estrogen_Status = "Positive" and Progesterone_Status = "Negative" then EP = 2;
else if Estrogen_Status = "Negative" and Progesterone_Status = "Positive" then EP = 3;
else if Estrogen_Status = "Negative" and Progesterone_Status = "Negative" then EP = 4;
else EP = .;

if Status = "Alive" then Outcome = 0;
else if Status = "Dead" then Outcome = 1;

run;

/* Formats for readability */
proc format;
value cat
  1 = "Both Positive (ER+/PR+)"
  2 = "ER+ / PR-"
  3 = "ER- / PR+"
  4 = "Both Negative (ER-/PR-)";
value stat
  0 = "Alive"
  1 = "Dead";
run;

/* Check EP (missing?) */
proc freq data=work.data1;
tables EP;
format EP cat.;
run;

proc contents data=work.data1 order=varnum;
run;

/* 3) DESCRIPTIVE ANALYSIS */
proc freq data=work.data1;
tables Race Marital_Status Grade A_Stage Estrogen_Status Progesterone_Status EP Status;
format EP cat.;
run;

proc means data=work.data1 mean std median maxdec=2;
var Age Tumor_Size Survival_Months Reginol_Node_Positive;
run;

/* Optional descriptive by Status */
proc sort data=work.data1;
by Status;
run;

proc freq data=work.data1;
by Status;
tables Race Marital_Status Grade A_Stage Estrogen_Status Progesterone_Status;
run;

proc freq data=work.data1;
tables Status;
run;

/*4) KAPLAN–MEIER CURVE */
ods graphics on;

/* A) 4-category EP */
proc lifetest data=work.data1 method=km plots=(s(CL atrisk));
    time Survival_Months*Outcome(0);
    strata EP;
    format EP cat.;
run;

/* B) Estrogen_Status (2-category) */
proc lifetest data=work.data1 method=km plots=(s(CL atrisk));
    time Survival_Months*Outcome(0);
    strata Estrogen_Status;
run;

/* C) Progesterone_Status (2-category) */
proc lifetest data=work.data1 method=km plots=(s(CL atrisk));
    time Survival_Months*Outcome(0);
    strata Progesterone_Status;
run;

/* Pairwise EP comparisons (log-rank with Bonferroni) */
proc lifetest data=work.data1;
    time Survival_Months*Outcome(0);
    strata EP / test=logrank adjust=bon;
    format EP cat.;
run;

/*  5) COX PROPORTIONAL HAZARDS */

/* A) Model including ER and PR separately */
proc phreg data=work.data1;
class Race(param=ref ref="White")
      A_Stage(param=ref ref="Regional")
      Estrogen_Status(param=ref ref="Positive")
      Progesterone_Status(param=ref ref="Positive");

model Survival_Months*Outcome(0)= Age Race A_Stage Tumor_Size Reginol_Node_Positive
                                 Estrogen_Status Progesterone_Status / rl;

estimate "HR for Estrogen (Negative vs Positive)" Estrogen_Status 1 / exp cl;
estimate "HR for Progesterone (Negative vs Positive)" Progesterone_Status 1 / exp cl;
run;

/* B) Estrogen only */
proc phreg data=work.data1;
class Race(param=ref ref="White")
      A_Stage(param=ref ref="Regional")
      Estrogen_Status(param=ref ref="Positive");

model Survival_Months*Outcome(0)= Age Race A_Stage Tumor_Size Reginol_Node_Positive
                                 Estrogen_Status / rl;

estimate "HR for Estrogen (Negative vs Positive)" Estrogen_Status 1 / exp cl;
run;

/* C) Progesterone only */
proc phreg data=work.data1;
class Race(param=ref ref="White")
      A_Stage(param=ref ref="Regional")
      Progesterone_Status(param=ref ref="Positive");

model Survival_Months*Outcome(0)= Age Race A_Stage Tumor_Size Reginol_Node_Positive
                                 Progesterone_Status / rl;

estimate "HR for Progesterone (Negative vs Positive)" Progesterone_Status 1 / exp cl;
run;

/* Check association between ER and PR (chi-square) */
proc freq data=work.data1;
tables Progesterone_Status*Estrogen_Status / chisq;
run;

/* D) Cox model using combined EP (recommended) */
proc phreg data=work.data1;
class Race(param=ref ref="White")
      A_Stage(param=ref ref="Regional")
      EP(param=ref ref="4");

model Survival_Months*Outcome(0)= Age Race A_Stage Tumor_Size Reginol_Node_Positive EP / rl;

estimate "HR: EP1 (ER+/PR+) vs EP4 (ER-/PR-)" EP 1 0 0 / exp cl;
estimate "HR: EP2 (ER+/PR-) vs EP4 (ER-/PR-)" EP 0 1 0 / exp cl;
estimate "HR: EP3 (ER-/PR+) vs EP4 (ER-/PR-)" EP 0 0 1 / exp cl;

format EP cat.;
run;

/*  6) PH ASSUMPTION TEST (ASSESS PH)*/
ods graphics on;

proc phreg data=work.data1 plots(overlay=row)=(survival);
class Race(param=ref ref="White")
      A_Stage(param=ref ref="Regional")
      EP(param=ref ref="4");

model Survival_Months*Outcome(0)= Age Race A_Stage Tumor_Size Reginol_Node_Positive EP;
assess ph / resample;

format EP cat.;
run;

/*  7)TIME-DEPENDENT EXPLORATION*/

/* Create indicator variables for EP1 and EP2 so the model runs */
data work.data1_td;
set work.data1;
EP1 = (EP=1);
EP2 = (EP=2);
run;

proc phreg data=work.data1_td;
class Race(param=ref ref="White")
      A_Stage(param=ref ref="Regional")
      EP(param=ref ref="4");

model Survival_Months*Outcome(0)= Age Race Tumor_Size Reginol_Node_Positive A_Stage EP
                                 A_Stage*log(Survival_Months)
                                 EP1*log(Survival_Months)
                                 EP2*log(Survival_Months);
format EP cat.;
run;

/* 8) CHECK PH USING LOG-LOG SURVIVAL PLOTS */

/* EP */
proc lifetest data=work.data1 plots=(s,lls);
time Survival_Months*Outcome(0);
strata EP;
format EP cat.;
run;

/* Race */
proc lifetest data=work.data1 plots=(s,lls);
time Survival_Months*Outcome(0);
strata Race;
run;

/* Stage */
proc lifetest data=work.data1 plots=(s,lls);
time Survival_Months*Outcome(0);
strata A_Stage;
run;

/* Create cutpoints based on means (from your report): Age~54, Tumor_Size~30.47, Nodes~4.16 */
proc format;
value agef       low-54 = "Below Mean"
                 54<-high = "Mean or above";
value tsizef     low-30.47 = "Below Mean"
                 30.47<-high = "Mean or above";
value regposf    low-4.16 = "Below Mean"
                 4.16<-high = "Mean or above";
run;

/* Age (categorized) */
proc lifetest data=work.data1 plots=(s,lls);
time Survival_Months*Outcome(0);
strata Age;
format Age agef.;
run;

/* Tumor size (categorized) */
proc lifetest data=work.data1 plots=(s,lls);
time Survival_Months*Outcome(0);
strata Tumor_Size;
format Tumor_Size tsizef.;
run;

/* Regional node positive (categorized) */
proc lifetest data=work.data1 plots=(s,lls);
time Survival_Months*Outcome(0);
strata Reginol_Node_Positive;
format Reginol_Node_Positive regposf.;
run;

/* 9) ADJUSTED SURVIVAL CURVES (PHREG BASELINE)*/
data work.one_level;
length Race $59 A_Stage $8;
Age = 53.97;
Tumor_Size = 30.47;
Reginol_Node_Positive = 4.16;
Race = "White";
A_Stage = "Regional";
EP = 4;
run;

ods graphics on;
proc phreg data=work.data1 plots(overlay=row)=(survival);
class Race(param=ref ref="White")
      A_Stage(param=ref ref="Regional")
      EP(param=ref ref="4");

title "Adjusted survival curves from PHREG (covariates held constant)";
model Survival_Months*Outcome(0)= Age Race A_Stage Tumor_Size Reginol_Node_Positive;
strata EP;

baseline covariates=work.one_level out=work.curvedata
         survival=survgraph loglogs=llsurvgraph;

format EP cat.;
run;

proc sgplot data=work.curvedata;
title "Log-log survival vs survival time after adjustment";
series x=Survival_Months y=llsurvgraph / group=EP;
run;
quit;

10) SCHOENFELD RESIDUALS (CORR WITH FAILURE TIME RANK)
proc phreg data=work.data1 plots(overlay=row)=(survival);
class Race(param=ref ref="White")
      A_Stage(param=ref ref="Regional")
      EP(param=ref ref="4");

model Survival_Months*Outcome(0)= Age Race A_Stage Tumor_Size Reginol_Node_Positive;

output out=work.resid_plots
       wtressch=wtschage wtschrace wtschstage wtschsize wtschrnode
       ressch=schage schrace schstage schsize schrnode;

format EP cat.;
run;
data work.failed_only;
set work.resid_plots;
where Outcome=1;
run;

/* Rank failure times */
proc rank data=work.failed_only out=work.ranked ties=mean;
var Survival_Months;
ranks timerank;
run;

/* Correlation tests */
proc corr data=work.ranked;
with timerank;
var schage schrace schstage schsize schrnode;
run;

proc corr data=work.ranked;
with timerank;
var wtschage wtschrace wtschstage wtschsize wtschrnode;
run;

data work.resid_plots2;
set work.ranked;
logtime = log(Survival_Months);
run;

proc sgscatter data=work.resid_plots2;
label w1="wtschage" w2="wtschrace" w3="wtschstage" w4="wtschsize" w5="wtschcrnode"
      sch1="schage" sch2="schrace" sch3="schstage" sch4="schsize" sch5="schrnode";

plot wtschage*logtime schage*timerank
     wtschrace*logtime schrace*timerank
     wtschstage*logtime schstage*timerank
     wtschsize*logtime schsize*timerank
     wtschcrnode*logtime schrnode*timerank;
run;
quit;

ods graphics off;

title;
footnote;
