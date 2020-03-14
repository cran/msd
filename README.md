
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
m <- matrix(d, nrow = 50, ncol = 20)
z1 <- msd(m, misfit = TRUE)

# z1 contains the estimated item measures, person measures, rating 
# category thresholds and their standard errors as well as other 
# diagnostic statistics such infit, outfit and reliability measures
```

If item measures and thresholds are already known from a large study and
we want to estimate person measures for a different sample of patients,
use the pms function.

``` r
# Estimate person measures based on known item measures and thresholds.
d2 <- as.numeric(sample(0:5, 200, replace = TRUE))
m2 <- matrix(d2, nrow = 10, ncol = 20)
im = z1$item_measures
th = z1$thresholds
pm1 <- pms(data = m2, items = im, thresholds = th, misfit = TRUE)
```
