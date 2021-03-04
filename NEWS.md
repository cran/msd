
Version 0.3.1
=============

Modified the msd function to estimate parameters even when at least one dichotomization produces no parameter estimates. A dichotomization will produce no parameter estimates if an extreme rating category is never used (or rarely used). Previously, msd would not run unless all dichotomizations produced valid parameter estimates. 

Version 0.3.0
=============

Added function ims, which estimates item measures when person measures and thresholds are already known. This function can be used to estimate the item measure of any new item added to a questionnaire (i.e., use person measures estimated from the old items to calibrate the new items).

Added the option to specify a vector of person minus item measure differences in msdprob, making it easier to generate probability curves (a graphical representation of the probability of rating an item with any rating category given estimated rating category thresholds).

The dependency on optim is removed for function pms and replaced with an iterative method for estimating person measures and standard errors.

Standard errors in the function msd are now calculated without using imputed values from the individual dichotomizations.


Version 0.2.0
=============

Added function simdata, which generates simulated rating scale data and
allows msd to be tested for accuracy.

Added the option to specify a wider range of possible rating categories
for functions pms and misfit than may be present in the data. This is
important because pms and misfit are useful for estimating person
measures and misfit statistics for individual persons, and it is
possible that the person did not use all rating categories.


Version 0.1.0
=============

Created NEWS.md file
