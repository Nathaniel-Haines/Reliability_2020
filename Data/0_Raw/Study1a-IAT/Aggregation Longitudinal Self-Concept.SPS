***NEW ALGORITHM TIME 1***

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\iat1.sav'.
RECODE rt  (10000 thru Highest=SYSMIS).
EXECUTE .
SELECT IF Block=6.
EXECUTE.
RECODE correct ('True'=1) INTO cor.
EXECUTE.
SELECT IF cor=1.
EXECUTE.
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Block31.sav"
 /break subjt1
 /rt_31_t1 = MEAN(rt)
 /sd_31_t1 = SD(rt)
 /cor_31_t1 = SUM(cor).
EXECUTE.

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\iat1.sav'.
RECODE rt  (10000 thru Highest=SYSMIS).
EXECUTE .
SELECT IF Block=7.
EXECUTE.
RECODE correct ('True'=1) INTO cor.
EXECUTE.
SELECT IF cor=1.
EXECUTE.
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Block32.sav"
 /break subjt1
 /rt_32_t1 = MEAN(rt)
 /sd_32_t1 = SD(rt)
 /cor_32_t1 = SUM(cor).
EXECUTE.

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\iat1.sav'.
RECODE rt  (10000 thru Highest=SYSMIS).
EXECUTE .
SELECT IF Block=11.
EXECUTE.
RECODE correct ('True'=1) INTO cor.
EXECUTE.
SELECT IF cor=1.
EXECUTE.
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Block51.sav"
 /break subjt1
 /rt_51_t1 = MEAN(rt)
 /sd_51_t1 = SD(rt)
 /cor_51_t1 = SUM(cor).
EXECUTE.

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\iat1.sav'.
RECODE rt  (10000 thru Highest=SYSMIS).
EXECUTE .
SELECT IF Block=12.
EXECUTE.
RECODE correct ('True'=1) INTO cor.
EXECUTE.
SELECT IF cor=1.
EXECUTE.
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Block52.sav"
 /break subjt1
 /rt_52_t1 = MEAN(rt)
 /sd_52_t1 = SD(rt)
 /cor_52_t1 = SUM(cor).
EXECUTE.

MATCH FILES 
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Block31.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Block32.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Block51.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Block52.SAV'
 /BY subjt1.
EXECUTE.

SORT CASES BY subjt1 (A) .

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\implicit_new.SAV'.
EXECUTE.

COMPUTE sd_Block1_t1 = SQRT( (((cor_31_t1-1) * sd_31_t1 * sd_31_t1) + ((cor_51_t1-1) * sd_51_t1 * sd_51_t1)) / (cor_31_t1 + cor_51_t1 - 2) ).
EXECUTE.
COMPUTE sd_Block2_t1 = SQRT( (((cor_32_t1-1) * sd_32_t1 * sd_32_t1) + ((cor_52_t1-1) * sd_52_t1 * sd_52_t1)) / (cor_32_t1 + cor_52_t1 - 2) ).
EXECUTE.

COMPUTE rt_31n_t1 = ((cor_31_t1*rt_31_t1) + ((40 - cor_31_t1)*(rt_31_t1+0))) / 40.
EXECUTE.
COMPUTE rt_32n_t1 = ((cor_32_t1*rt_32_t1) + ((40 - cor_32_t1)*(rt_32_t1+0))) / 40.
EXECUTE.

COMPUTE rt_51n_t1 = ((cor_51_t1*rt_51_t1) + ((40 - cor_51_t1)*(rt_51_t1+0))) / 40.
EXECUTE.
COMPUTE rt_52n_t1 = ((cor_52_t1*rt_52_t1) + ((40 - cor_52_t1)*(rt_52_t1+0))) / 40.
EXECUTE.

COMPUTE iat_mx1new_t1 = (rt_51n_t1 - rt_31n_t1) / sd_Block1_t1.
EXECUTE.
COMPUTE iat_mx2new_t1 = (rt_52n_t1 - rt_32n_t1) / sd_Block2_t1.
EXECUTE.

COMPUTE iat_new_t1 = MEAN(iat_mx1new_t1, iat_mx2new_t1).
EXECUTE.
VARIABLE LABELS iat_new_t1 'IAT extraversion-introversion t1 (new algorithm)'.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\implicit_new.SAV'
  /COMPRESSED.


***AGGREGATION TIME 1***

MATCH FILES 
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Start1.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\implicit_new.SAV'
 /BY subjt1.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Longitudinal Self-Concept Time 1.SAV'
  /COMPRESSED.

RECODE int1 int2 int3 int4 int5 (1=5) (2=4) (3=3) (4=2) (5=1) INTO  int1r int2r int3r int4r int5r.
EXECUTE.

COMPUTE expl_traits_t1 = MEAN (ext1, ext2, ext3, ext4, ext5, int1r, int2r, int3r, int4r, int5r).
EXECUTE.

COMPUTE errorrateB3_t1 = 1 - ((cor_31_t1 + cor_32_t1) / 80).
EXECUTE .
COMPUTE errorrateB5_t1 = 1 - ((cor_51_t1 + cor_52_t1) / 80).
EXECUTE .

SORT CASES BY subjectnum(A).

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 1\Longitudinal Self-Concept Time 1.SAV'
  /COMPRESSED.


***NEW ALGORITHM TIME 2***

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\iat2.sav'.
RECODE rt  (10000 thru Highest=SYSMIS).
EXECUTE .
SELECT IF Block=6.
EXECUTE.
RECODE correct ('True'=1) INTO cor.
EXECUTE.
SELECT IF cor=1.
EXECUTE.
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Block31.sav"
 /break subjt2
 /rt_31_t2 = MEAN(rt)
 /sd_31_t2 = SD(rt)
 /cor_31_t2 = SUM(cor).
EXECUTE.

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\iat2.sav'.
RECODE rt  (10000 thru Highest=SYSMIS).
EXECUTE .
SELECT IF Block=7.
EXECUTE.
RECODE correct ('True'=1) INTO cor.
EXECUTE.
SELECT IF cor=1.
EXECUTE.
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Block32.sav"
 /break subjt2
 /rt_32_t2 = MEAN(rt)
 /sd_32_t2 = SD(rt)
 /cor_32_t2 = SUM(cor).
EXECUTE.

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\iat2.sav'.
RECODE rt  (10000 thru Highest=SYSMIS).
EXECUTE .
SELECT IF Block=11.
EXECUTE.
RECODE correct ('True'=1) INTO cor.
EXECUTE.
SELECT IF cor=1.
EXECUTE.
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Block51.sav"
 /break subjt2
 /rt_51_t2 = MEAN(rt)
 /sd_51_t2 = SD(rt)
 /cor_51_t2 = SUM(cor).
EXECUTE.

GET FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\iat2.sav'.
RECODE rt  (10000 thru Highest=SYSMIS).
EXECUTE .
SELECT IF Block=12.
EXECUTE.
RECODE correct ('True'=1) INTO cor.
EXECUTE.
SELECT IF cor=1.
EXECUTE.
AGGREGATE OUTFILE "C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Block52.sav"
 /break subjt2
 /rt_52_t2 = MEAN(rt)
 /sd_52_t2 = SD(rt)
 /cor_52_t2 = SUM(cor).
EXECUTE.

MATCH FILES 
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Block31.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Block32.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Block51.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Block52.SAV'
 /BY subjt2.
EXECUTE.

SORT CASES BY subjt2 (A) .

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\implicit_new.SAV'.
EXECUTE.

COMPUTE sd_Block1_t2 = SQRT( (((cor_31_t2-1) * sd_31_t2 * sd_31_t2) + ((cor_51_t2-1) * sd_51_t2 * sd_51_t2)) / (cor_31_t2 + cor_51_t2 - 2) ).
EXECUTE.
COMPUTE sd_Block2_t2 = SQRT( (((cor_32_t2-1) * sd_32_t2 * sd_32_t2) + ((cor_52_t2-1) * sd_52_t2 * sd_52_t2)) / (cor_32_t2 + cor_52_t2 - 2) ).
EXECUTE.

COMPUTE rt_31n_t2 = ((cor_31_t2*rt_31_t2) + ((40 - cor_31_t2)*(rt_31_t2+0))) / 40.
EXECUTE.
COMPUTE rt_32n_t2 = ((cor_32_t2*rt_32_t2) + ((40 - cor_32_t2)*(rt_32_t2+0))) / 40.
EXECUTE.

COMPUTE rt_51n_t2 = ((cor_51_t2*rt_51_t2) + ((40 - cor_51_t2)*(rt_51_t2+0))) / 40.
EXECUTE.
COMPUTE rt_52n_t2 = ((cor_52_t2*rt_52_t2) + ((40 - cor_52_t2)*(rt_52_t2+0))) / 40.
EXECUTE.

COMPUTE iat_mx1new_t2 = (rt_51n_t2 - rt_31n_t2) / sd_Block1_t2.
EXECUTE.
COMPUTE iat_mx2new_t2 = (rt_52n_t2 - rt_32n_t2) / sd_Block2_t2.
EXECUTE.


COMPUTE iat_new_t2 = MEAN(iat_mx1new_t2, iat_mx2new_t2).
EXECUTE.
VARIABLE LABELS iat_new_t2 'IAT extraversion-introversion t2 (new algorithm)'.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\implicit_new.SAV'
  /COMPRESSED.


***AGGREGATION TIME 2***

MATCH FILES 
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Start2.SAV'
 /FILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\implicit_new.SAV'
 /BY subjt2.
EXECUTE.

SAVE OUTFILE='C:\Users\bg24397\Documents\Documents\Research\Implicit Social Cognition\Longitudinal\Study 1a - Self-Concept - IAT\Data\Session 2\Longitudinal Self-Concept Time 2.SAV'
  /COMPRESSED.

RECODE int1t2 int2t2 int3t2 int4t2 int5t2 (1=5) (2=4) (3=3) (4=2) (5=1) INTO  int1t2r int2t2r int3t2r int4t2r int5t2r.
EXECUTE.

COMPUTE expl_traits_t2 = MEAN (ext1t2, ext2t2, ext3t2, ext4t2, ext5t2, int1t2r, int2t2r, int3t2r, int4t2r, int5t2r).
EXECUTE.

COMPUTE errorrateB3_t2 = 1 - ((cor_31_t2 + cor_32_t2) / 80).
EXECUTE .
COMPUTE errorrateB5_t2 = 1 - ((cor_51_t2 + cor_52_t2) / 80).
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