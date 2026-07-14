/* Descriptive analysis section from survival_analysis.sas */
/* PROC FREQ over the categorical predictors + PROC MEANS over the numerics. */
/* Runs against work.data1 (built in autoexec.sas from the derivation logic in the source). */

/* Check EP (missing?) */
proc freq data=work.data1;
tables EP;
format EP cat.;
run;

/* 3) DESCRIPTIVE ANALYSIS */
proc freq data=work.data1;
tables Race Marital_Status Grade A_Stage Estrogen_Status Progesterone_Status EP Status;
format EP cat.;
run;

proc means data=work.data1 mean std median maxdec=2;
var Age Tumor_Size Survival_Months Reginol_Node_Positive;
run;

proc freq data=work.data1;
tables Status;
run;
