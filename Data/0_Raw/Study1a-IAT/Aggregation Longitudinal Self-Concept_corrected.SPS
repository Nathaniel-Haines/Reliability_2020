***We thank Christoph Klauer for making us aware of several errors in the earlier IAT data aggregation syntax [Aggregation Longitudinal Self-Concept.sps] posted at https://osf.io/792qj/***
***(1) Means and standard deviations of RTs were based on RTs from trials with correct responses instead of all trials.*** 
***(2) The earlier syntax presupposes 80 valid trials per combined block, which cannot assumed after eliminating trials with RTs > 10,000 ms.***
***(3) The earlier syntax included a typo that led to an omission of the 600 ms error penalty.***
***These errors are corrected in the current syntax, which follows the D-score aggregation guidelines described in Table 4 in Greenwald, Nosek, & Banaji (2003, p. 214).*** 
***The correction of these errors leads to small deviations in the descriptive and inferential statistics reported in Gawronski, Morrison, Phills, & Galdi (2017, Study 1a).***
***The corrections do not affect any of the main findings and substantive conclusions.***
***An official corrigendum with the correct statistical results has been submitted to Personality and Social Psychology Bulletin***


***NEW ALGORITHM TIME 1***

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\iat1.sav'.
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
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\iat_aggregate1.sav"
 /break subjt1
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

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\iat_aggregate1.sav'.

SORT CASES BY subjt1 (A) .

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
VARIABLE LABELS iat_new_t1 'IAT extraversion-introversion t1 (new algorithm)'.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\iat_aggregate1.SAV'
  /COMPRESSED.


***AGGREGATION TIME 1***

MATCH FILES 
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Start1.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\iat_aggregate1.SAV'
 /BY subjt1.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Longitudinal Self-Concept Time 1.SAV'
  /COMPRESSED.

RECODE int1 int2 int3 int4 int5 (1=5) (2=4) (3=3) (4=2) (5=1) INTO  int1r int2r int3r int4r int5r.
EXECUTE.

COMPUTE expl_traits_t1 = MEAN (ext1, ext2, ext3, ext4, ext5, int1r, int2r, int3r, int4r, int5r).
EXECUTE.

COMPUTE errorrateB3_t1 = 1 - (cor_31_t1 + cor_32_t1) / (cor_31_t1 + cor_32_t1 + err_31_t1 + err_32_t1).
EXECUTE .
COMPUTE errorrateB5_t1 = 1 - (cor_51_t1 + cor_52_t1) / (cor_51_t1 + cor_52_t1 + err_51_t1 + err_52_t1).
EXECUTE .

SORT CASES BY subjectnum(A).

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Longitudinal Self-Concept Time 1.SAV'
  /COMPRESSED.


***NEW ALGORITHM TIME 2***

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\iat2.sav'.
SELECT IF ((Block = 6) or (Block = 7) or (Block = 11) or (Block = 12)).
EXECUTE. 
RECODE rt  (10000 thru Highest=SYSMIS).
EXECUTE .
RECODE correct ('True'=1) INTO cor.
EXECUTE.
RECODE correct ('False'=1) INTO err.
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
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\iat_aggregate2.sav"
 /break subjt2
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

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\iat_aggregate2.sav'.

SORT CASES BY subjt2 (A) .

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
VARIABLE LABELS iat_new_t2 'IAT extraversion-introversion t2 (new algorithm)'.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\iat_aggregate2.SAV'
  /COMPRESSED.


***AGGREGATION TIME 2***

MATCH FILES 
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Start2.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\iat_aggregate2.SAV'
 /BY subjt2.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Longitudinal Self-Concept Time 2.SAV'
  /COMPRESSED.

RECODE int1t2 int2t2 int3t2 int4t2 int5t2 (1=5) (2=4) (3=3) (4=2) (5=1) INTO  int1t2r int2t2r int3t2r int4t2r int5t2r.
EXECUTE.

COMPUTE expl_traits_t2 = MEAN (ext1t2, ext2t2, ext3t2, ext4t2, ext5t2, int1t2r, int2t2r, int3t2r, int4t2r, int5t2r).
EXECUTE.

COMPUTE errorrateB3_t2 = 1 - (cor_31_t2 + cor_32_t2) / (cor_31_t2 + cor_32_t2 + err_31_t2 + err_32_t2).
EXECUTE .
COMPUTE errorrateB5_t2 = 1 - (cor_51_t2 + cor_52_t2) / (cor_51_t2 + cor_52_t2 + err_51_t2 + err_52_t2).
EXECUTE .

SORT CASES BY subjectnum(A).

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Longitudinal Self-Concept Time 2.SAV'
  /COMPRESSED.


***MERGING FILES***

GET FILE = "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Longitudinal Self-Concept Time 1.sav".

MATCH FILES /FILE=*
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Longitudinal Self-Concept Time 2.sav'
 /BY subjectnum.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Longitudinal Self-Concept.sav'
  /COMPRESSED.


***ATTRITION ANALYSIS***

T-TEST GROUPS=crit(1 2)
  /MISSING=ANALYSIS
  /VARIABLES=iat_new_t1 expl_traits_t1
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
  /VARIABLES=ext1 ext2 ext3 ext4 ext5 int1r int2r int3r int4r int5r 
  /SCALE('ALL VARIABLES')  ALL/MODEL=ALPHA
  /SUMMARY=TOTAL .

RELIABILITY
  /VARIABLES=ext1t2 ext2t2 ext3t2 ext4t2 ext5t2 int1t2r int2t2r int3t2r int4t2r int5t2r 
  /SCALE('ALL VARIABLES')  ALL/MODEL=ALPHA
  /SUMMARY=TOTAL .


***LONGITUDINAL ANALYSIS***

T-TEST PAIRS=expl_traits_t1 iat_new_t1 WITH expl_traits_t2 iat_new_t2 (PAIRED)
  /CRITERIA=CI(.9500)
  /MISSING=ANALYSIS.

CORRELATIONS
  /VARIABLES=expl_traits_t1 expl_traits_t2 iat_new_t1 iat_new_t2
  /PRINT=TWOTAIL NOSIG
  /MISSING=PAIRWISE.

FILTER OFF.
USE ALL.
EXECUTE.