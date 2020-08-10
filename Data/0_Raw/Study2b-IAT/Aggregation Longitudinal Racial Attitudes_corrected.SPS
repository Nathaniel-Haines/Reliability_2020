***We thank Christoph Klauer for making us aware of several errors in the earlier IAT data aggregation syntax [Aggregation Longitudinal Racial Attitudes.sps] posted at https://osf.io/792qj/***
***(1) Means and standard deviations of RTs were based on RTs from trials with correct responses instead of all trials.*** 
***(2) The earlier syntax presupposes 60 valid trials per combined block, which cannot assumed after eliminating trials with RTs > 10,000 ms.***
***(3) The earlier syntax included a typo that led to an omission of the 600 ms error penalty.***
***These errors are corrected in the current syntax, which follows the D-score aggregation guidelines described in Table 4 in Greenwald, Nosek, & Banaji (2003, p. 214).*** 
***The correction of these errors leads to small deviations in the descriptive and inferential statistics reported in Gawronski, Morrison, Phills, & Galdi (2017, Study 2b).***
***The corrections do not affect any of the main findings and substantive conclusions.***
***An official corrigendum with the correct statistical results has been submitted to Personality and Social Psychology Bulletin***


***IAT TIME 1***

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 1\iat1.sav'.
SELECT IF ((Block = 6) or (Block = 7) or (Block = 11) or (Block = 12)).
EXECUTE. 
RECODE rt  (10000 thru Highest=SYSMIS).
EXECUTE .
RECODE correct ('True'=1) INTO cor.
EXECUTE.
RECODE correct ('False'=1) INTO err.
EXECUTE.
IF (Block = 6) rt_all31_t1 = rt.
EXECUTE.
IF ((Block = 6) and (correct = 'True')) rt_cor31_t1 = rt.
EXECUTE.
IF ((Block = 6) and (cor = 1)) cor_31_t1 = cor.
EXECUTE.
IF ((Block = 6) and (err = 1)) err_31_t1 = err.
EXECUTE.
IF (Block = 7) rt_all32_t1 = rt.
EXECUTE.
IF ((Block = 7) and (correct = 'True')) rt_cor32_t1 = rt.
EXECUTE.
IF ((Block = 7) and (cor = 1)) cor_32_t1 = cor.
EXECUTE.
IF ((Block = 7) and (err = 1)) err_32_t1 = err.
EXECUTE.
IF (Block = 11) rt_all51_t1 = rt.
EXECUTE.
IF ((Block = 11) and (correct = 'True')) rt_cor51_t1 = rt.
EXECUTE.
IF ((Block = 11) and (cor = 1)) cor_51_t1 = cor.
EXECUTE.
IF ((Block = 11) and (err = 1)) err_51_t1 = err.
EXECUTE.
IF (Block = 12) rt_all52_t1 = rt.
EXECUTE.
IF ((Block = 12) and (correct = 'True')) rt_cor52_t1 = rt.
EXECUTE.
IF ((Block = 12) and (cor = 1)) cor_52_t1 = cor.
EXECUTE.
IF ((Block = 12) and (err = 1)) err_52_t1 = err.
EXECUTE.
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 1\iat_aggregate1.sav"
 /break subj1
 /rt_all31_t1 = MEAN(rt_all31_t1)
 /rt_cor31_t1 = MEAN(rt_cor31_t1)
 /sd_all31_t1 = SD(rt_all31_t1)
 /cor_31_t1 = SUM(cor_31_t1)
 /err_31_t1 = SUM(err_31_t1)
 /rt_all32_t1 = MEAN(rt_all32_t1)
 /rt_cor32_t1 = MEAN(rt_cor32_t1)
 /sd_all32_t1 = SD(rt_all32_t1)
 /cor_32_t1 = SUM(cor_32_t1)
 /err_32_t1 = SUM(err_32_t1)
 /rt_all51_t1 = MEAN(rt_all51_t1)
 /rt_cor51_t1 = MEAN(rt_cor51_t1)
 /sd_all51_t1 = SD(rt_all51_t1)
 /cor_51_t1 = SUM(cor_51_t1)
 /err_51_t1 = SUM(err_51_t1)
 /rt_all52_t1 = MEAN(rt_all52_t1)
 /rt_cor52_t1 = MEAN(rt_cor52_t1)
 /sd_all52_t1 = SD(rt_all52_t1)
 /cor_52_t1 = SUM(cor_52_t1)
 /err_52_t1 = SUM(err_52_t1) . 

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 1\iat_aggregate1.sav'.

SORT CASES BY subj1 (A) .

RECODE cor_31_t1 cor_32_t1 cor_51_t1 cor_52_t1 err_31_t1 err_32_t1 err_51_t1 err_52_t1 (MISSING=0).
EXECUTE.

COMPUTE sd_Block1_t1 = SQRT((((cor_31_t1-1 + err_31_t1-1) * sd_all31_t1 * sd_all31_t1) + ((cor_51_t1-1 + err_51_t1-1) * sd_all51_t1 * sd_all51_t1)) / (cor_31_t1 + err_31_t1-1 + cor_51_t1 + err_51_t1-1 - 2)).
EXECUTE.
COMPUTE sd_Block2_t1 = SQRT((((cor_32_t1-1 + err_32_t1-1) * sd_all32_t1 * sd_all32_t1) + ((cor_52_t1-1 + err_51_t1-1) * sd_all52_t1 * sd_all52_t1)) / (cor_32_t1 + err_32_t1-1 + cor_52_t1 + err_52_t1-1 - 2)).
EXECUTE.

COMPUTE rt_31n_t1 = ((cor_31_t1*rt_cor31_t1) + ((err_31_t1)*(rt_cor31_t1+600))) / (cor_31_t1 + err_31_t1).
EXECUTE.
COMPUTE rt_32n_t1 = ((cor_32_t1*rt_cor32_t1) + ((err_32_t1)*(rt_cor32_t1+600))) / (cor_32_t1 + err_32_t1).
EXECUTE.

COMPUTE rt_51n_t1 = ((cor_51_t1*rt_cor51_t1) + ((err_51_t1)*(rt_cor51_t1+600))) / (cor_51_t1 + err_51_t1).
EXECUTE.
COMPUTE rt_52n_t1 = ((cor_52_t1*rt_cor52_t1) + ((err_52_t1)*(rt_cor52_t1+600))) / (cor_52_t1 + err_52_t1).
EXECUTE.

COMPUTE iat_mx1new_t1 = (rt_51n_t1 - rt_31n_t1) / sd_Block1_t1.
EXECUTE.
COMPUTE iat_mx2new_t1 = (rt_52n_t1 - rt_32n_t1) / sd_Block2_t1.
EXECUTE.

COMPUTE iat_new_t1 = MEAN(iat_mx1new_t1, iat_mx2new_t1).
EXECUTE.
VARIABLE LABELS iat_new_t1 'IAT race bias t1 (new algorithm)'.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 1\iat_aggregate1.SAV'
  /COMPRESSED.


***AGGREGATION TIME 1***

MATCH FILES 
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 1\Start1.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 1\iat_aggregate1.SAV'
 /BY subj1.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 1\Longitudinal Race IAT Time 1.SAV'
  /COMPRESSED.

COMPUTE expl_pwb_ex01_t1 = w01t1 - b01t1.
EXECUTE.
COMPUTE expl_pwb_ex02_t1 = w02t1 - b02t1.
EXECUTE.
COMPUTE expl_pwb_ex03_t1 = w03t1 - b03t1.
EXECUTE.
COMPUTE expl_pwb_ex04_t1 = w04t1 - b04t1.
EXECUTE.
COMPUTE expl_pwb_ex05_t1 = w05t1 - b05t1.
EXECUTE.
COMPUTE expl_pwb_ex06_t1 = w06t1 - b06t1.
EXECUTE.
COMPUTE expl_pwb_ex07_t1 = w07t1 - b07t1.
EXECUTE.
COMPUTE expl_pwb_ex08_t1 = w08t1 - b08t1.
EXECUTE.
COMPUTE expl_pwb_ex09_t1 = w09t1 - b09t1.
EXECUTE.
COMPUTE expl_pwb_ex10_t1 = w10t1 - b10t1.
EXECUTE.

COMPUTE expl_pwb_cat1_t1 = cat1wt1 - cat1bt1.
EXECUTE.
COMPUTE expl_pwb_cat2_t1 = cat2wt1 - cat2bt1.
EXECUTE.
COMPUTE expl_pwb_cat3_t1 = cat3wt1 - cat3bt1.
EXECUTE.
COMPUTE expl_pwb_cat4_t1 = cat4wt1 - cat4bt1.
EXECUTE.
COMPUTE expl_pwb_cat5_t1 = cat5wt1 - cat5bt1.
EXECUTE.

COMPUTE expl_pwb_ex_t1 = MEAN(expl_pwb_ex01_t1, expl_pwb_ex02_t1, expl_pwb_ex03_t1, expl_pwb_ex04_t1, expl_pwb_ex05_t1, expl_pwb_ex06_t1, expl_pwb_ex07_t1, expl_pwb_ex08_t1, expl_pwb_ex09_t1, expl_pwb_ex10_t1).
EXECUTE.

COMPUTE expl_pwb_cat_t1 = MEAN(expl_pwb_cat1_t1, expl_pwb_cat2_t1, expl_pwb_cat3_t1, expl_pwb_cat4_t1, expl_pwb_cat5_t1).
EXECUTE.

COMPUTE errorrateB3_t1 = 1 - (cor_31_t1 + cor_32_t1) / (cor_31_t1 + cor_32_t1 + err_31_t1 + err_32_t1).
EXECUTE .
COMPUTE errorrateB5_t1 = 1 - (cor_51_t1 + cor_52_t1) / (cor_51_t1 + cor_52_t1 + err_51_t1 + err_52_t1).
EXECUTE .

SORT CASES BY code(A).

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 1\Longitudinal Race IAT Time 1.SAV'
  /COMPRESSED.


***IAT TIME 2***

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 2\iat2.sav'.
SELECT IF ((Block = 6) or (Block = 7) or (Block = 11) or (Block = 12)).
EXECUTE. 
RECODE rt  (10000 thru Highest=SYSMIS).
EXECUTE .
RECODE correct ('True'=1) INTO cor.
EXECUTE.
RECODE correct ('Fals'=1) INTO err.
EXECUTE.
IF (Block = 6) rt_all31_t2 = rt.
EXECUTE.
IF ((Block = 6) and (correct = 'True')) rt_cor31_t2 = rt.
EXECUTE.
IF ((Block = 6) and (cor = 1)) cor_31_t2 = cor.
EXECUTE.
IF ((Block = 6) and (err = 1)) err_31_t2 = err.
EXECUTE.
IF (Block = 7) rt_all32_t2 = rt.
EXECUTE.
IF ((Block = 7) and (correct = 'True')) rt_cor32_t2 = rt.
EXECUTE.
IF ((Block = 7) and (cor = 1)) cor_32_t2 = cor.
EXECUTE.
IF ((Block = 7) and (err = 1)) err_32_t2 = err.
EXECUTE.
IF (Block = 11) rt_all51_t2 = rt.
EXECUTE.
IF ((Block = 11) and (correct = 'True')) rt_cor51_t2 = rt.
EXECUTE.
IF ((Block = 11) and (cor = 1)) cor_51_t2 = cor.
EXECUTE.
IF ((Block = 11) and (err = 1)) err_51_t2 = err.
EXECUTE.
IF (Block = 12) rt_all52_t2 = rt.
EXECUTE.
IF ((Block = 12) and (correct = 'True')) rt_cor52_t2 = rt.
EXECUTE.
IF ((Block = 12) and (cor = 1)) cor_52_t2 = cor.
EXECUTE.
IF ((Block = 12) and (err = 1)) err_52_t2 = err.
EXECUTE.
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 2\iat_aggregate2.sav"
 /break subj2
 /rt_all31_t2 = MEAN(rt_all31_t2)
 /rt_cor31_t2 = MEAN(rt_cor31_t2)
 /sd_all31_t2 = SD(rt_all31_t2)
 /cor_31_t2 = SUM(cor_31_t2)
 /err_31_t2 = SUM(err_31_t2)
 /rt_all32_t2 = MEAN(rt_all32_t2)
 /rt_cor32_t2 = MEAN(rt_cor32_t2)
 /sd_all32_t2 = SD(rt_all32_t2)
 /cor_32_t2 = SUM(cor_32_t2)
 /err_32_t2 = SUM(err_32_t2)
 /rt_all51_t2 = MEAN(rt_all51_t2)
 /rt_cor51_t2 = MEAN(rt_cor51_t2)
 /sd_all51_t2 = SD(rt_all51_t2)
 /cor_51_t2 = SUM(cor_51_t2)
 /err_51_t2 = SUM(err_51_t2)
 /rt_all52_t2 = MEAN(rt_all52_t2)
 /rt_cor52_t2 = MEAN(rt_cor52_t2)
 /sd_all52_t2 = SD(rt_all52_t2)
 /cor_52_t2 = SUM(cor_52_t2)
 /err_52_t2 = SUM(err_52_t2) . 

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 2\iat_aggregate2.sav'.

SORT CASES BY subj2 (A) .

RECODE cor_31_t2 cor_32_t2 cor_51_t2 cor_52_t2 err_31_t2 err_32_t2 err_51_t2 err_52_t2 (MISSING=0).
EXECUTE.

COMPUTE sd_Block1_t2 = SQRT((((cor_31_t2-1 + err_31_t2-1) * sd_all31_t2 * sd_all31_t2) + ((cor_51_t2-1 + err_51_t2-1) * sd_all51_t2 * sd_all51_t2)) / (cor_31_t2 + err_31_t2-1 + cor_51_t2 + err_51_t2-1 - 2)).
EXECUTE.
COMPUTE sd_Block2_t2 = SQRT((((cor_32_t2-1 + err_32_t2-1) * sd_all32_t2 * sd_all32_t2) + ((cor_52_t2-1 + err_51_t2-1) * sd_all52_t2 * sd_all52_t2)) / (cor_32_t2 + err_32_t2-1 + cor_52_t2 + err_52_t2-1 - 2)).
EXECUTE.

COMPUTE rt_31n_t2 = ((cor_31_t2*rt_cor31_t2) + ((err_31_t2)*(rt_cor31_t2+600))) / (cor_31_t2 + err_31_t2).
EXECUTE.
COMPUTE rt_32n_t2 = ((cor_32_t2*rt_cor32_t2) + ((err_32_t2)*(rt_cor32_t2+600))) / (cor_32_t2 + err_32_t2).
EXECUTE.

COMPUTE rt_51n_t2 = ((cor_51_t2*rt_cor51_t2) + ((err_51_t2)*(rt_cor51_t2+600))) / (cor_51_t2 + err_51_t2).
EXECUTE.
COMPUTE rt_52n_t2 = ((cor_52_t2*rt_cor52_t2) + ((err_52_t2)*(rt_cor52_t2+600))) / (cor_52_t2 + err_52_t2).
EXECUTE.

COMPUTE iat_mx1new_t2 = (rt_51n_t2 - rt_31n_t2) / sd_Block1_t2.
EXECUTE.
COMPUTE iat_mx2new_t2 = (rt_52n_t2 - rt_32n_t2) / sd_Block2_t2.
EXECUTE.

COMPUTE iat_new_t2 = MEAN(iat_mx1new_t2, iat_mx2new_t2).
EXECUTE.
VARIABLE LABELS iat_new_t2 'IAT race bias t2 (new algorithm)'.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 2\iat_aggregate2.SAV'
  /COMPRESSED.


***AGGREGATION TIME 2***

MATCH FILES 
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 2\Start2.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 2\iat_aggregate2.SAV'
 /BY subj2.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 2\Longitudinal Race IAT Time 2.SAV'
  /COMPRESSED.

COMPUTE expl_pwb_ex01_t2 = w01t2 - b01t2.
EXECUTE.
COMPUTE expl_pwb_ex02_t2 = w02t2 - b02t2.
EXECUTE.
COMPUTE expl_pwb_ex03_t2 = w03t2 - b03t2.
EXECUTE.
COMPUTE expl_pwb_ex04_t2 = w04t2 - b04t2.
EXECUTE.
COMPUTE expl_pwb_ex05_t2 = w05t2 - b05t2.
EXECUTE.
COMPUTE expl_pwb_ex06_t2 = w06t2 - b06t2.
EXECUTE.
COMPUTE expl_pwb_ex07_t2 = w07t2 - b07t2.
EXECUTE.
COMPUTE expl_pwb_ex08_t2 = w08t2 - b08t2.
EXECUTE.
COMPUTE expl_pwb_ex09_t2 = w09t2 - b09t2.
EXECUTE.
COMPUTE expl_pwb_ex10_t2 = w10t2 - b10t2.
EXECUTE.

COMPUTE expl_pwb_cat1_t2 = cat1wt2 - cat1bt2.
EXECUTE.
COMPUTE expl_pwb_cat2_t2 = cat2wt2 - cat2bt2.
EXECUTE.
COMPUTE expl_pwb_cat3_t2 = cat3wt2 - cat3bt2.
EXECUTE.
COMPUTE expl_pwb_cat4_t2 = cat4wt2 - cat4bt2.
EXECUTE.
COMPUTE expl_pwb_cat5_t2 = cat5wt2 - cat5bt2.
EXECUTE.

COMPUTE expl_pwb_ex_t2 = MEAN(expl_pwb_ex01_t2, expl_pwb_ex02_t2, expl_pwb_ex03_t2, expl_pwb_ex04_t2, expl_pwb_ex05_t2, expl_pwb_ex06_t2, expl_pwb_ex07_t2, expl_pwb_ex08_t2, expl_pwb_ex09_t2, expl_pwb_ex10_t2).
EXECUTE.

COMPUTE expl_pwb_cat_t2 = MEAN(expl_pwb_cat1_t2, expl_pwb_cat2_t2, expl_pwb_cat3_t2, expl_pwb_cat4_t2, expl_pwb_cat5_t2).
EXECUTE.

COMPUTE errorrateB3_t2 = 1 - (cor_31_t2 + cor_32_t2) / (cor_31_t2 + cor_32_t2 + err_31_t2 + err_32_t2).
EXECUTE .
COMPUTE errorrateB5_t2 = 1 - (cor_51_t2 + cor_52_t2) / (cor_51_t2 + cor_52_t2 + err_51_t2 + err_52_t2).
EXECUTE .

SORT CASES BY code(A).

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 2\Longitudinal Race IAT Time 2.SAV'
  /COMPRESSED.


***MERGING FILES***

GET FILE = "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 1\Longitudinal Race IAT Time 1.sav".

MATCH FILES /FILE=*
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Session 2\Longitudinal Race IAT Time 2.sav'
 /BY code.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 2b - Racial Attitudes - IAT\Data\Longitudinal Race IAT.sav'
  /COMPRESSED.


***ATTRITION ANALYSIS***

T-TEST GROUPS=crit(1 2)
  /MISSING=ANALYSIS
  /VARIABLES=iat_new_t1 expl_pwb_cat_t1 expl_pwb_ex_t1 
  /CRITERIA=CI(.95).

USE ALL.
COMPUTE filter_$=(crit = 2).
VARIABLE LABEL filter_$ 'crit = 3 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMAT filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.


***RELIABILITY ANALYSIS***

RELIABILITY
  /VARIABLES=iat_mx1new_t1 iat_mx2new_t1
  /FORMAT=NOLABELS
  /SCALE(ALPHA)=ALL/MODEL=ALPHA.

RELIABILITY
  /VARIABLES=iat_mx1new_t2 iat_mx2new_t2
  /FORMAT=NOLABELS
  /SCALE(ALPHA)=ALL/MODEL=ALPHA.

RELIABILITY
  /VARIABLES=expl_pwb_cat1_t1 expl_pwb_cat2_t1 expl_pwb_cat3_t1 expl_pwb_cat4_t1 expl_pwb_cat5_t1 
  /SCALE('ALL VARIABLES')  ALL/MODEL=ALPHA
  /SUMMARY=TOTAL .

RELIABILITY
  /VARIABLES=expl_pwb_cat1_t2 expl_pwb_cat2_t2 expl_pwb_cat3_t2 expl_pwb_cat4_t2 expl_pwb_cat5_t2 
  /SCALE('ALL VARIABLES')  ALL/MODEL=ALPHA
  /SUMMARY=TOTAL .

RELIABILITY
  /VARIABLES=expl_pwb_ex01_t1 expl_pwb_ex02_t1 expl_pwb_ex03_t1 expl_pwb_ex04_t1 expl_pwb_ex05_t1 expl_pwb_ex06_t1 expl_pwb_ex07_t1 expl_pwb_ex08_t1 expl_pwb_ex09_t1 expl_pwb_ex10_t1
  /SCALE('ALL VARIABLES')  ALL/MODEL=ALPHA
  /SUMMARY=TOTAL .

RELIABILITY
  /VARIABLES=expl_pwb_ex01_t2 expl_pwb_ex02_t2 expl_pwb_ex03_t2 expl_pwb_ex04_t2 expl_pwb_ex05_t2 expl_pwb_ex06_t2 expl_pwb_ex07_t2 expl_pwb_ex08_t2 expl_pwb_ex09_t2 expl_pwb_ex10_t2
  /SCALE('ALL VARIABLES')  ALL/MODEL=ALPHA
  /SUMMARY=TOTAL .


***LONGITUDINAL ANALYSIS***

T-TEST PAIRS=iat_new_t1 expl_pwb_cat_t1 expl_pwb_ex_t1 WITH iat_new_t2 expl_pwb_cat_t2 expl_pwb_ex_t2 (PAIRED)
  /CRITERIA=CI(.9500)
  /MISSING=ANALYSIS.

CORRELATIONS
  /VARIABLES=iat_new_t1 iat_new_t2 expl_pwb_cat_t1 expl_pwb_cat_t2 expl_pwb_ex_t1 expl_pwb_ex_t2
  /PRINT=TWOTAIL NOSIG
  /MISSING=PAIRWISE.

FILTER OFF.
USE ALL.
EXECUTE.
