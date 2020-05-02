
<!-- README.md is generated from README.Rmd. Please edit that file -->

# msd

<!-- badges: start -->

<!-- badges: end -->

The goal of msd is to implement functions related to the Method of
Successive Dichotomizations (MSD), which is a computationally efficient
method of estimating parameters for a polytomous Rasch model that
estimates ordered rating category thresholds. MSD applies to ordinal
rating scale data and is an improvement over polytomous Rasch models
such as the Andrich rating scale model and the Partial Credit Model that
often estimate disordered thresholds.

## Installation

You can install the released version of msd from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("msd")
```

## Example

Suppose a health status questionnaire is administered to patients who
rate the difficulty of items (e.g., tasks they must perform) using
ordinal rating categories (e.g., integers from 0-5). We can use msd to
estimate item measures (estimates of task difficulty), person measures
(estimates of patient abilities) and rating category thresholds that
define the average rating scale used by the sample of patients.

``` r
library(msd)

# Using a randomly generated ratings matrix with ratings from 0 to 5
d <- as.numeric(sample(0:5, 1000, replace = TRUE))
dm <- matrix(d, nrow = 50, ncol = 20)
m <- msd(dm, misfit = TRUE)

# m contains the estimated item measures, person measures, rating 
# category thresholds and their standard errors as well as other 
# diagnostic statistics such infit, outfit and reliability measures
```

If item measures and thresholds are already known from a large study and
we want to estimate person measures for a different sample of patients,
use the pms function.

``` r
# Estimate person measures based on known item measures and thresholds.
d2 <- as.numeric(sample(0:5, 200, replace = TRUE))
dm2 <- matrix(d2, nrow = 10, ncol = 20)
im = m$item_measures
th = m$thresholds
pm <- pms(data = dm2, items = im, thresholds = th, misfit = TRUE)
```

The accuracy of msd can be tested using the simdata function, which
generates simulated rating scale data given item measures, person
measures and thresholds. Because msd sets the axis origin at the mean
item measure, “true” item measures should have a mean of zero.
Similarly, because a non-zero mean threshold changes the person measure
by that same amount, the “true” mean threshold should be set to
zero.

``` r
# Test the accuracy of msd. First, generate simulated rating scale data with 15% missing 
# data and the rating scale going from 0 to 5.
im = runif(100, -2, 2)
im = im - mean(im)
pm = runif(200, -2, 2)
th = sort(runif(5, -2, 2))
th = th - mean(th)
d = simdata(im, pm, th, missingProb = 0.15, minRating = 0)

# Compare msd estimated parameters to true values.  Linear regression should yield a slope
# very close to 1 and an intercept very close to 0.
m = msd(d)
lm(m$item_measures ~ im)
#> 
#> Call:
#> lm(formula = m$item_measures ~ im)
#> 
#> Coefficients:
#> (Intercept)           im  
#>  -2.306e-17    1.026e+00
lm(m$person_measures ~ pm)
#> 
#> Call:
#> lm(formula = m$person_measures ~ pm)
#> 
#> Coefficients:
#> (Intercept)           pm  
#>    -0.01133      1.01436
lm(m$thresholds ~ th)
#> 
#> Call:
#> lm(formula = m$thresholds ~ th)
#> 
#> Coefficients:
#> (Intercept)           th  
#>    -0.00224      1.00751
```
