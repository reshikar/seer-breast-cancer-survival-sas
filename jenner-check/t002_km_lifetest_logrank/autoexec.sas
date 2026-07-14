options obs=100;   /* cap input rows for the captured run */

/* --- MOCK SEER INPUT (substitutes the external data_raw/SEER.csv) --- */
/* Column shape derived from survival_analysis.sas variable usage. */
data work.sdata;
    length Race $8 Marital_Status $10 Grade $2 A_Stage $8
           Estrogen_Status $8 Progesterone_Status $8 Status $5;
    input Race $ Marital_Status $ Grade $ A_Stage $
          Estrogen_Status $ Progesterone_Status $ Status $
          Age Tumor_Size Survival_Months Reginol_Node_Positive;
    datalines;
Other Married 1 Distant Positive Positive Alive 45 33 21 22
Other Married 4 Regional Positive Negative Dead 31 16 18 20
White Single 4 Regional Negative Positive Alive 58 80 39 1
White Widowed 3 Distant Negative Negative Dead 39 32 60 3
Black Married 3 Distant Positive Positive Dead 68 38 64 15
Other Married 4 Regional Positive Negative Alive 65 42 84 12
Other Single 1 Regional Negative Positive Dead 72 34 61 3
White Married 4 Distant Negative Negative Dead 59 86 30 12
White Divorced 1 Regional Positive Positive Alive 64 36 24 9
Other Single 3 Regional Positive Negative Alive 44 9 107 9
White Single 3 Regional Negative Positive Alive 71 68 54 21
Black Single 3 Regional Negative Negative Dead 45 76 43 19
Black Widowed 3 Regional Positive Positive Dead 38 70 40 2
White Single 2 Distant Positive Negative Dead 68 13 31 15
Other Divorced 1 Regional Negative Positive Dead 73 73 60 21
Black Married 3 Distant Negative Negative Alive 40 63 4 24
Black Single 1 Distant Positive Positive Alive 83 86 68 5
Black Single 1 Distant Positive Negative Alive 61 7 18 10
White Married 2 Regional Negative Positive Alive 35 67 12 18
White Single 4 Regional Negative Negative Dead 46 72 48 7
Other Single 3 Distant Positive Positive Alive 72 88 51 17
Black Married 2 Regional Positive Negative Alive 34 48 6 8
Other Single 1 Regional Negative Positive Dead 75 85 6 2
Black Married 2 Distant Negative Negative Dead 72 67 18 24
Other Widowed 2 Distant Positive Positive Dead 81 57 16 22
Black Divorced 4 Distant Positive Negative Alive 59 11 90 21
White Married 4 Distant Negative Positive Dead 81 18 21 18
Black Single 4 Regional Negative Negative Alive 47 64 35 3
Black Married 1 Regional Positive Positive Alive 35 35 25 16
White Widowed 1 Regional Positive Negative Dead 54 5 31 15
Black Widowed 4 Regional Negative Positive Alive 42 42 31 19
Other Married 3 Regional Negative Negative Alive 33 79 65 17
White Married 1 Regional Positive Positive Alive 34 81 12 8
Black Married 2 Regional Positive Negative Alive 69 15 57 19
Other Divorced 3 Regional Negative Positive Dead 72 45 20 5
Other Divorced 4 Distant Negative Negative Dead 78 14 3 19
White Married 2 Distant Positive Positive Alive 38 49 12 12
Black Single 4 Distant Positive Negative Dead 69 88 42 18
Black Married 2 Distant Negative Positive Alive 37 18 99 9
Black Single 3 Regional Negative Negative Alive 73 86 37 9
White Widowed 3 Regional Positive Positive Dead 30 47 61 9
Other Widowed 4 Regional Positive Negative Alive 37 14 92 18
Black Single 4 Regional Positive Positive Alive 32 44 50 2
White Single 1 Distant Negative Positive Alive 79 76 56 24
White Single 2 Distant Positive Negative Alive 31 27 98 14
Black Single 1 Distant Positive Negative Dead 85 9 38 15
Black Single 2 Regional Negative Positive Alive 72 29 55 3
Black Widowed 3 Regional Negative Positive Alive 37 38 26 9
White Widowed 3 Distant Positive Positive Dead 57 82 41 19
Black Married 4 Regional Positive Negative Alive 63 73 91 24
Black Widowed 1 Distant Positive Negative Alive 69 45 88 24
Other Divorced 4 Distant Negative Positive Dead 55 42 44 14
Other Single 3 Distant Negative Negative Dead 65 5 25 14
Black Widowed 4 Regional Negative Positive Alive 62 65 105 24
Other Married 3 Distant Positive Negative Dead 35 35 54 7
White Married 2 Distant Positive Negative Dead 69 14 37 21
Other Widowed 4 Distant Positive Negative Alive 45 23 87 25
Black Single 2 Distant Positive Positive Alive 33 76 35 4
White Widowed 3 Distant Negative Negative Alive 69 69 58 18
White Widowed 4 Distant Negative Negative Dead 78 36 51 25
Other Single 3 Distant Negative Negative Dead 34 41 20 11
White Single 2 Distant Positive Positive Dead 74 24 56 14
Black Widowed 4 Regional Negative Negative Alive 43 58 53 19
Other Widowed 4 Regional Positive Positive Alive 52 43 100 14
Black Single 3 Distant Positive Negative Dead 61 8 31 22
Other Single 4 Regional Negative Negative Alive 69 73 7 19
White Widowed 2 Distant Positive Positive Alive 41 11 37 7
Black Divorced 4 Distant Negative Negative Alive 78 58 36 16
Other Married 3 Regional Positive Positive Alive 71 13 103 2
White Single 1 Regional Positive Positive Alive 45 21 64 19
;
run;

/* 2) CREATE VARIABLES  (verbatim from survival_analysis.sas) */
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
