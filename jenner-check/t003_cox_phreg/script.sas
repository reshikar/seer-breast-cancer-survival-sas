/* Cox proportional-hazards section from survival_analysis.sas */
/* PROC PHREG on the combined EP receptor variable with ESTIMATE statements for the */
/* pairwise hazard ratios against the ER-/PR- reference group. */
/* Runs against work.data1 (built in autoexec.sas from the derivation logic in the source). */

/* 5) COX PROPORTIONAL HAZARDS - combined EP (recommended) */
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

/* Check association between ER and PR (chi-square) */
proc freq data=work.data1;
tables Progesterone_Status*Estrogen_Status / chisq;
run;
