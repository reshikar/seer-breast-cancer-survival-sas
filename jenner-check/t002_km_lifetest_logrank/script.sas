/* Kaplan-Meier + log-rank section from survival_analysis.sas */
/* PROC LIFETEST product-limit curves stratified by EP receptor group, plus the */
/* Bonferroni-adjusted log-rank test of equality over strata. */
/* Runs against work.data1 (built in autoexec.sas from the derivation logic in the source). */

/* 4) KAPLAN-MEIER CURVE */
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

/* Pairwise EP comparisons (log-rank with Bonferroni) */
proc lifetest data=work.data1;
    time Survival_Months*Outcome(0);
    strata EP / test=logrank adjust=bon;
    format EP cat.;
run;
