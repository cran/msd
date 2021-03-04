# Implements the Method of Successive Dichotomizations.

# REFERENCE:
# Bradley C & Massof RW. Method of successive dichotomizations: an
# improved method for estimating measures of latent variables from
# rating scale data. PLoS ONE 2018; 13(10): e0206106. https://doi.
# org/10.1371/journal.pone.0206106 PMID: 30335832

# Included functions:
# msd <- Method of Successive Dichotomizations
# rasch <- dichotomous Rasch model (special case of msd)
# thresh <- estimates msd thresholds given person and item measures
# pms <- estimates person measures given item measures and thresholds
# ims <- estimates item measures given person measures and thresholds
# msdprob <- calculates rating category probabilities
# expdata <- generates expected ratings for each person/item combination
# misfit <- calculates misfit statistics (infit and outfit)
# simdata <- generates simulated rating scale data

# LIBRARIES
library(stats)

##########################################################################

# Method of Successive Dichotomizations
msd <- function(data, items = NULL, persons = NULL, misfit = FALSE){

  # Estimates item measures, person measures, rating category thresholds
  # and their standard errors using the Method of Successive
  # Dichotomizations. Option provided for anchoring certain items and
  # persons while estimating the rest. Option also provided for
  # estimating infit and outfit statistics.

  # REQUIRED INPUTS
  # 'data' = a numeric matrix of ordinal rating scale data whose entries
  #          are integers with missing data set to NA. Rows are persons
  #          and columns are items. The ordinal rating scale is assumed
  #          to go from the smallest to largest integer in integer steps.

  # OPTIONAL INPUTS
  # 'items' = a numeric vector of anchored item measures. Item measures
  #           to be estimated are set to NA. The length of 'items' must
  #           equal the number of columns in 'data'. Default is NULL.
  # 'persons' = a numeric vector of anchored person measures. Person
  #             measures to be estimated are set to NA. The length of
  #             'persons' must equal the number of rows in 'data'.
  #             Default is NULL.
  # 'misfit' = logical for calculating infit and outfit statistics.
  #            Default is FALSE.


  # CHECK FOR MISSING REQUIRED INPUTS ------------------------------------

  if (missing("data")){
    stop('Argument "data" is missing with no default')
  }

  # CHECK DATA STRUCTURE AND DATA TYPE FOR RATING SCALE DATA -------------

  if (!is.matrix(data) || typeof(data) != "double"){
    stop('Argument "data" must be a matrix of mode numeric whose entries
         are either integers or NA.')
  }

  # DEFAULT VALUES -------------------------------------------------------

  # Default NA numeric vectors for optional inputs "items" and "persons"
  dm = dim(data)
  if (is.null(items)){
    items = as.numeric(NA*seq(1, dm[2]))
  }
  if (is.null(persons)){
    persons = as.numeric(NA*seq(1, dm[1]))
  }

  # CHECK OTHER DATA STRUCTURES AND DATA TYPES ---------------------------

  if (!is.vector(items) || typeof(items) != "double"){
    stop('Argument "items" must be a vector of mode numeric with all
         non-fixed item measures set to NA')
  }
  if (!is.vector(persons) || typeof(persons) != "double"){
    stop('Argument "persons" must be a vector of mode numeric with all
         non-fixed person measures set to NA')
  }
  if (!is.logical(misfit)){
    stop('Argument "misfit" must be a logical TRUE or FALSE')
  }

  # CHECK VALUES AND DIMENSIONS ------------------------------------------

  # Check if non-NA entries of data are integers
  if (all(data%%1 == 0, na.rm = TRUE) == FALSE){
    stop('Entries in "data" must consist of integers or NA')}

  # Check if there are at least two different integers
  if (length(unique(c(data[!is.na(data)]))) < 2){
    stop('There must be at least two different integers in "data"')
  }

  # Compare dimensions of data with vectors items and persons
  if (length(items) != dim(data)[2]){
    stop('Argument "items" must be a vector whose length is the
         number of columns in input matrix "data"')
  }
  if (length(persons) != dim(data)[1]){
    stop('Argument "persons" must be a vector whose length is the
         number of rows in input matrix "data"')
  }

  # INITIALIZE -----------------------------------------------------------

  # Set mean item measures to zero if there are no fixed items or persons.
  mZero = 0
  if (all(is.na(items)) && all(is.na(persons))){
    mZero = 1
  }

  # Make lowest rating zero
  rdata = data - min(data, na.rm = TRUE)

  # Highest rating category = number of thresholds
  R_max = max(rdata, na.rm = TRUE)

  # Use dichotomous Rasch model if there are only two rating categories
  if (R_max == 1){

    # OUTPUT -------------------------------------------------------------
    out = rasch(rdata, items, persons, misfit)

  } else {

    # Initialize person and item measure matrices for dichotomizations
    PM = PE = matrix(NA, dim(rdata)[1], R_max)
    IM = IE = matrix(NA, R_max, dim(rdata)[2])

    # PERSON AND ITEM MEASURES -------------------------------------------

    # Loop through dichotomizations to estimate person and item measures
    MM = rdata
    for (dd in seq(1, R_max)){

      # Dichotomize ratings matrix
      M2 = MM; M2[M2 < dd] = 0; M2[M2 >= dd] = 1

      # DICHOTOMOUS RASCH MODEL -----------------------------------------

      # Dichotomous Rasch model
      RD = rasch(M2, items, persons)

      # Person and item measures
      PM[,dd] = RD$person_measures
      IM[dd,] = RD$item_measures

      # Standard errors
      PE[,dd] = RD$person_std_errors
      IE[dd,] = RD$item_std_errors
    }

    # Impute item measures if estimate is NA for some dichotomization
    IM_dif_m = colMeans(IM, na.rm = TRUE)
    IM_dif_m2 = matrix(IM_dif_m, nrow = dim(IM)[1], ncol = dim(IM)[2],
                       byrow = TRUE)
    IM_imp = IM
    IM_imp[is.na(IM_imp)] = IM_dif_m2[is.na(IM_imp)]

    # Impute person measures if estimate is NA for some dichotomization
    PM_m = colMeans(PM, na.rm = TRUE)
    PM_m2 = matrix(PM_m, nrow = dim(PM)[1], ncol = dim(PM)[2], byrow = TRUE)
    PM_dif = PM - PM_m2
    PM_dif_m = rowMeans(PM_dif, na.rm = TRUE)
    PM_dif_m2 = matrix(PM_dif_m, nrow = dim(PM)[1], ncol = dim(PM)[2])
    PM_mdif = PM_m2 + PM_dif_m2
    PM_imp = PM
    PM_imp[is.na(PM_imp)] = PM_mdif[is.na(PM_imp)]

    # Final person and item measures
    IM_F = colMeans(IM_imp, na.rm = TRUE)
    PM_F = rowMeans(PM_imp, na.rm = TRUE)
    if (mZero == 1){
      IM_F = IM_F - mean(colMeans(IM, na.rm = TRUE), na.rm = TRUE)
      PM_F = PM_F - mean(colMeans(IM, na.rm = TRUE), na.rm = TRUE)
    }

    # Change non-finite numbers to NA
    IM_F[!is.finite(IM_F)] = NA
    PM_F[!is.finite(PM_F)] = NA

    # Fix items and persons
    IM_F[!is.na(items)] = items[!is.na(items)]
    PM_F[!is.na(persons)] = persons[!is.na(persons)]

    # Indexes of estimated person and item measures
    IM_ind = which(!is.na(IM_F))
    PM_ind = which(!is.na(PM_F))

    # THRESHOLDS ---------------------------------------------------------

    # Matrix of responses for estimated person and item measures
    MR = rdata[PM_ind,IM_ind]

    # Estimate MSD thresholds
    MSD_T = thresh(MR, IM_F[IM_ind], PM_F[PM_ind])

    # STANDARD ERRORS ----------------------------------------------------

    # Person and item standard errors
    PE_num = rowSums(!is.na(PE))
    PSE = sqrt(rowSums(PE^2, na.rm = TRUE)/PE_num)/sqrt(PE_num)
    IE_num = colSums(!is.na(IE))
    ISE = sqrt(colSums(IE^2, na.rm = TRUE)/IE_num)/sqrt(IE_num)

    # Change to NA if person or item measure is NA
    PSE[is.na(PM_F)] = NA
    ISE[is.na(IM_F)] = NA

    # Set standard error to NA if items or persons are fixed
    PSE[!is.na(persons)] = NA
    ISE[!is.na(items)] = NA

    # RELIABILITY --------------------------------------------------------

    # Reliability measures for items and persons
    RL_items = 1 - (mean(ISE, na.rm = TRUE))^2/
      (sd(IM_F, na.rm = TRUE))^2
    RL_persons = 1 - (mean(PSE, na.rm = TRUE))^2/
      (sd(PM_F, na.rm = TRUE))^2

    # Change to NA if reliability is Inf
    RL_items[!is.finite(RL_items)] = NA
    RL_persons[!is.finite(RL_persons)] = NA

    # Change to NA if reliability is less than zero or above 1
    RL_items[RL_items < 0] = NA
    RL_items[RL_items > 1] = NA
    RL_persons[RL_persons < 0] = NA
    RL_persons[RL_persons > 1] = NA

    # MISFIT STATISTICS --------------------------------------------------

    if (misfit == TRUE){
      MSF = misfit(rdata, IM_F, PM_F, MSD_T$thresholds)
    }

    # OUTPUT -------------------------------------------------------------
    out = c()
    out$item_measures = IM_F
    out$person_measures = PM_F
    out$thresholds = MSD_T$thresholds
    out$item_std_errors = ISE
    out$person_std_errors = PSE
    out$threshold_std_errors = MSD_T$threshold_std_errors
    out$item_reliability = RL_items
    out$person_reliability = RL_persons

    if (misfit == TRUE){
      out$infit_items = MSF$infit_items
      out$infit_persons = MSF$infit_persons
      out$outfit_items = MSF$outfit_items
      out$outfit_persons = MSF$outfit_persons
    }
  }

  return(out)
}



# Dichotomous Rasch model
rasch <- function(data, items = NULL, persons = NULL, misfit = FALSE){
  # Estimates item measures, person measures and their standard errors
  # using the dichotomous Rasch model. A special case of the function
  # msd when the rating scale consists of only two rating categories: 0
  # and 1. Option provided for anchoring certain items and persons while
  # estimating the rest. Option also provided for estimating infit and
  # outfit statistics.

  # REQUIRED INPUTS
  # 'data' = a numeric matrix of 0's and 1's with missing data set to NA.
  #          Rows are persons and columns are items.

  # OPTIONAL INPUTS
  # 'items' = a numeric vector of anchored item measures. Item measures
  #           to be estimated are set to NA. Default is NULL.
  # 'persons' = a numeric vector of anchored personm measures. Person
  #             measures to be estimated are set to NA. Default is NULL.
  # 'misfit' = logical for calculating infit and outfit statistics.
  #            Default is FALSE.


  # CHECK FOR MISSING REQUIRED INPUTS ------------------------------------

  if (missing("data")){
    stop('Argument "data" is missing with no default')
  }

  # CHECK DATA STRUCTURE AND DATA TYPE FOR RATING SCALE DATA -------------

  if (!is.matrix(data) || typeof(data) != "double"){
    stop('Argument "data" must be a matrix of mode numeric whose entries
         are 1, 0 or NA.')
  }

  # DEFAULT VALUES -------------------------------------------------------

  # Default NA numeric vectors for optional inputs "items" and "persons"
  dm = dim(data)
  if (is.null(items)){
    items = as.numeric(NA*seq(1,dm[2]))
  }
  if (is.null(persons)){
    persons = as.numeric(NA*seq(1,dm[1]))
  }

  # CHECK OTHER DATA STRUCTURES AND DATA TYPES ---------------------------

  if (!is.vector(items) || typeof(items) != "double"){
    stop('Argument "items" must be a vector of mode numeric with all
         non-fixed item measures set to NA')
  }
  if (!is.vector(persons) || typeof(persons) != "double"){
    stop('Argument "persons" must be a vector of mode numeric with all
         non-fixed person measures set to NA')
  }
  if (!is.logical(misfit)){
    stop('Argument "misfit" must be a logical TRUE or FALSE')
  }

  # CHECK VALUES AND DIMENSIONS ------------------------------------------

  # Check if 'data' consists of only 0's, 1's and NA
  if (!identical(unique(c(data[!is.na(data)])), c(0,1)) &&
      !identical(unique(c(data[!is.na(data)])), c(1,0))){
    stop('Entries in "data" must be either 0, 1 or NA with at
         least one 0 and at least on 1')
  }

  # Compare dimensions of data with vectors items and persons
  if (length(items) != dim(data)[2]){
    stop('Argument "items" must be a vector whose length is the
         number of columns in input matrix "data"')
  }
  if (length(persons) != dim(data)[1]){
    stop('Argument "persons" must be a vector whose length is the
         number of rows in input matrix "data"')
  }


  # INITIALIZE -----------------------------------------------------------

  # Set mean item measure to zero unless items or persons are anchored
  mZero = 0
  if (all(is.na(items)) && all(is.na(persons))){
    mZero = 1
  }

  # Proportion correct/incorrect
  PC_person = rowSums(data, na.rm=TRUE)/rowSums(!is.na(data))
  PI_item = colSums(data, na.rm=TRUE)/colSums(!is.na(data))

  # Convert into logits
  PCP_logit = log(PC_person/(1-PC_person))
  PII_logit = log((1-PI_item)/PI_item)

  # Change Inf to NA
  PCP_logit[!is.finite(PCP_logit)] = NA
  PII_logit[!is.finite(PII_logit)] = NA

  # Set mean item measure to zero if there are no anchored items or persons
  if (mZero == 1){
    PII_logit = PII_logit - mean(PII_logit, na.rm = TRUE)
  }

  # Anchored item and person measures
  PII_logit[!is.na(items)] = items[!is.na(items)]
  PCP_logit[!is.na(persons)] = persons[!is.na(persons)]

  # Iterate until residual (res) is close enough to zero
  res = Inf; iter = 1;
  while(iter==1){
    # Expected values
    PM = matrix(rep(PCP_logit, dim(data)[2]), ncol=dim(data)[2])
    PI = t(matrix(rep(PII_logit, dim(data)[1]), ncol=dim(data)[1]))
    EV = exp(PM-PI)/(1+exp(PM-PI))

    # Variance of expected values
    V = EV*(1-EV)
    V_col = -1*colSums(V, na.rm=TRUE)
    V_row = -1*rowSums(V, na.rm=TRUE)

    # Residuals
    R = data-EV
    R_col = -1*colSums(R, na.rm=TRUE)
    R_row = rowSums(R, na.rm=TRUE)

    # Updated person and item measures
    PCP_logit = PCP_logit - R_row/V_row
    PII_logit = PII_logit - R_col/V_col
    if (mZero == 1){
      # Set mean item measure to zero
      PII_logit = PII_logit - mean(PII_logit, na.rm = TRUE)
    }

    # Fix item and person measures
    PII_logit[!is.na(items)] = items[!is.na(items)]
    PCP_logit[!is.na(persons)] = persons[!is.na(persons)]

    # Stop iteration
    res1 = sum(R_col^2) + sum(R_row^2)
    if (abs(res1-res) <= 0.001){
      iter = 0
    } else {res = res1}
  }

  # Standard errors
  SE_items = sqrt(-1/V_col)
  SE_persons = sqrt(-1/V_row)

  # Change non-finite numbers to NA
  PII_logit[!is.finite(PII_logit) | !is.finite(SE_items)] = NA
  PCP_logit[!is.finite(PCP_logit) | !is.finite(SE_persons)] = NA
  SE_items[!is.finite(PII_logit) | !is.finite(SE_items)] = NA
  SE_persons[!is.finite(PCP_logit) | !is.finite(SE_persons)] = NA

  # Change to NA if standard error is too large (here: max 10 logits)
  maxSE = 10
  PII_logit[SE_items > maxSE] = NA
  PCP_logit[SE_persons > maxSE] = NA
  SE_items[SE_items > maxSE] = NA
  SE_persons[SE_persons > maxSE] = NA

  # Fix item and person measures
  PII_logit[!is.na(items)] = items[!is.na(items)]
  PCP_logit[!is.na(persons)] = persons[!is.na(persons)]
  SE_items[!is.na(items)] = NA
  SE_persons[!is.na(persons)] = NA

  # Reliability measures for items and persons
  RL_items = 1 - (mean(SE_items, na.rm=TRUE))^2/
    (sd(PII_logit, na.rm=TRUE))^2
  RL_persons = 1 - (mean(SE_persons, na.rm=TRUE))^2/
    (sd(PCP_logit, na.rm=TRUE))^2

  # Change to NA if reliability is infinite
  RL_items[!is.finite(RL_items)] = NA
  RL_persons[!is.finite(RL_persons)] = NA

  # Change to NA if reliability is less than zero or above 1
  RL_items[RL_items < 0] = NA
  RL_items[RL_items > 1] = NA
  RL_persons[RL_persons < 0] = NA
  RL_persons[RL_persons > 1] = NA

  # OUTPUT -------------------------------------------------------------
  out = c()
  out$item_measures = PII_logit
  out$person_measures = PCP_logit
  out$item_std_errors = SE_items
  out$person_std_errors = SE_persons
  out$item_reliability = RL_items
  out$person_reliability = RL_persons

  # MISFIT STATISTICS --------------------------------------------------
  if (misfit == TRUE){
    MSF = misfit(data, PII_logit, PCP_logit, 0)

    out$infit_items = MSF$infit_items
    out$infit_persons = MSF$infit_persons
    out$outfit_items = MSF$outfit_items
    out$outfit_persons = MSF$outfit_persons
  }

  return(out)
}



# Estimate MSD thresholds
thresh <- function(data, items, persons){
  # Estimates rating category thresholds for msd given rating scale data,
  # item measures and person measures.

  # REQUIRED INPUTS
  # 'data' = a numeric matrix of ordinal rating scale data whose entries
  #          are integers with missing data set to NA. Rows are persons
  #          and columns are items. The ordinal rating scale is assumed
  #          to go from the smallest integer to the largest integer in
  #          'data' in integer steps.
  # 'items' = a numeric vector of item measures with missing values set
  #           to NA. The length of 'items' must equal the number of
  #           columns in 'data'. There must be at least one non-NA
  #           entry in 'items'.
  # 'persons' = a numeric vector of person measures with missing values
  #             set to NA. The length of 'persons' must equal the
  #             number of rows in 'data'. There must be at least one
  #             non-NA entry in 'persons'.

  # CHECK FOR MISSING REQUIRED INPUTS ------------------------------------

  if (missing("data")){
    stop('Argument "data" is missing with no default')
  }
  if (missing("items")){
    stop('Argument "items" is missing with no default')
  }
  if (missing("persons")){
    stop('Argument "persons" is missing with no default')
  }

  # CHECK DATA STRUCTURE AND DATA TYPES ----------------------------------

  if (!is.matrix(data) || typeof(data) != "double"){
    stop('Argument "data" must be a matrix of mode numeric whose entries
         are 1, 0 or NA.')
  }
  if (!is.vector(items) || typeof(items) != "double"){
    stop('Argument "items" must be a vector of mode numeric with all
         unknown item measures set to NA')
  }
  if (!is.vector(persons) || typeof(persons) != "double"){
    stop('Argument "persons" must be a vector of mode numeric with all
         unknown item measures set to NA')
  }

  # CHECK VALUES AND DIMENSIONS ------------------------------------------

  # Check if non-NA entries of data are integers
  if (all(data%%1 == 0, na.rm = TRUE) == FALSE){
    stop('Entries in "data" must be either 1, 0 or NA')}

  # Check if there are at least two different integers
  if (length(unique(c(data[!is.na(data)]))) < 2){
    stop('There must be at least one 1 and one 0 in "data"')
  }

  # Compare dimensions of data with vectors items and persons
  if (length(items) != dim(data)[2]){
    stop('Argument "items" must be a vector whose length is the
         number of columns in input matrix "data"')
  }
  if (length(persons) != dim(data)[1]){
    stop('Argument "persons" must be a vector whose length is the
         number of rows in input matrix "data"')
  }

  # Check if items and persons are all NA
  if (all(is.na(items))){
    stop('Argument "items" cannot be all NA')
  }
  if (all(is.na(persons))){
    stop('Argument "persons" cannot be all NA')
  }


  # INITIALIZE -----------------------------------------------------------

  # Make lowest rating zero
  rdata = data - min(data, na.rm = TRUE)

  # Gamma matrix for estimated person and item measures
  IM_rep = t(matrix(rep(items, dim(rdata)[1]), ncol=dim(rdata)[1]))
  PM_rep = matrix(rep(persons, dim(rdata)[2]), ncol=dim(rdata)[2])
  GM = PM_rep - IM_rep

  # Loop through every dichotomization to estimate thresholds
  R_max = max(rdata, na.rm = TRUE)
  TH = NA*seq(1, R_max); TH_se = NA*seq(1,R_max)
  for(dd in seq(1, R_max)){
    # Dichotomize ratings matrix
    MR2 = rdata; MR2[MR2<dd] = 0; MR2[MR2>=dd] = 1

    # Get gamma values for 0's and 1's
    G0 = GM[which(MR2==0)]
    G1 = GM[which(MR2==1)]

    # Remove NA
    G0 = G0[!is.na(G0)]
    G1 = G1[!is.na(G1)]

    # Function to be minimized (negative of log likelihood)
    f <- function(x) sum(log(1/exp(G1-x)+1)) + sum(log(exp(G0-x)+1))

    # Estimate threshold
    T_opt = optim(0, f, method = "Brent" , lower = -100,
                  upper = 100, hessian = TRUE)
    TH[dd] = T_opt$par
    TH_se[dd] = 1/sqrt(T_opt$hessian)
  }

  # Output
  out = c()
  out$thresholds = TH
  out$threshold_std_errors = TH_se

  return(out)
}



# Estimate MSD person measures given population thresholds
pms <- function(data, items, thresholds, misfit = FALSE, minRating = NULL){
  # Estimates person measures assuming all persons use the same rating
  # category thresholds.

  # REQUIRED INPUTS
  # 'data' = a numeric matrix of ordinal rating scale data whose entries
  #          are integers with missing data set to NA. Rows are persons
  #          and columns are items. The ordinal rating scale is assumed
  #          to go from the smallest to largest integer in integer steps
  #          unless minRating is specified in which case the ordinal
  #          rating scale goes from minRating to the largest integer in
  #          integer steps.
  # 'items' = a numeric vector of item measures with missing values set
  #           to NA.
  # 'thresholds' = a numeric vector of ordered rating category thresholds
  #                with no NA.

  # OPTIONAL INPUTS
  # 'misfit' = logical for calculating infit and outfit statistics.
  #            Default is FALSE.
  # 'minRating' = integer representing the smallest ordinal rating
  #               category. Default is NULL.


  # CHECK FOR MISSING REQUIRED INPUTS ------------------------------------

  if (missing("data")){
    stop('Argument "data" is missing with no default')
  }
  if (missing("items")){
    stop('Argument "items" is missing with no default')
  }
  if (missing("thresholds")){
    stop('Argument "thresholds" is missing with no default')
  }

  # CHECK DATA STRUCTURES AND DATA TYPES ---------------------------------

  if (!is.matrix(data) || typeof(data) != "double"){
    stop('Argument "data" must be a matrix of mode numeric whose entries
         are either integers or NA.')
  }
  if (!is.vector(items) || typeof(items) != "double"){
    stop('Argument "items" must be a vector of mode numeric with all
         unknown item measures set to NA')
  }
  if (!is.vector(thresholds) || typeof(thresholds) != "double"){
    stop('Argument "thresholds" must be an ordered vector of mode numeric
         with no NA')
  }
  if (!is.logical(misfit)){
    stop('Argument "misfit" must be a logical TRUE or FALSE')
  }
  if (!is.null(minRating) && length(minRating) != 1){
    stop('Argument "minRating" must be an integer')
  }
  if (!is.null(minRating) && minRating%%1 != 0){
    stop('Argument "minRating" must be an integer')
  }
  if (length(minRating) == 1 && typeof(minRating) == "integer"){
    minRating = as.numeric(minRating)
  }


  # CHECK VALUES AND DIMENSIONS ------------------------------------------

  # Check if non-NA entries of data are integers
  if (all(data%%1 == 0, na.rm = TRUE) == FALSE){
    stop('Entries in "data" must consist of integers or NA')}

  # Check if there are at least two different integers
  if (length(unique(c(data[!is.na(data)]))) < 2){
    stop('There must be at least two different integers in "data"')
  }

  # Compare dimensions of data with vector items
  if (length(items) != dim(data)[2]){
    stop('Argument "items" must be a vector whose length is the
         number of columns in input matrix "data"')
  }

  # Check if there are any NA in thresholds and if thresholds is ordered
  if (any(is.na(thresholds)) ||
      identical(order(thresholds),
                seq(1, length(thresholds))) == FALSE){
    stop('Argument "thresholds" must be an ordered vector of real numbers
         with no NA')
  }

  # Check if the length of thresholds equals the maximum observed rating
  # in 'data' minus 0 if 'minRating' is not specified
  if (is.null(minRating) &&
      max(data, na.rm = TRUE) - min(data, na.rm = TRUE)
      != length(thresholds)){
    stop('The length of "thresholds" must equal the maximum minus
         minimum rating in "data" unless "minRating" is specified')
  }
  # Check if the maximum observed rating is too large given 'minRating'
  if (!is.null(minRating) &&
      max(data, na.rm = TRUE) - minRating > length(thresholds)){
    stop('The maximum rating in "data" cannot be larger than "minRating"
         plus the length of "thresholds"')
  }
  # Check if 'minRating' is too large
  if (!is.null(minRating) && minRating > min(data, na.rm = TRUE)){
    stop('"minRating" cannot be larger than the smallest integer in
         "data"')
  }


  # INITIALIZE -----------------------------------------------------------

  # Define lowest rating to be zero
  if (is.null(minRating)){
    rdata = data - min(data, na.rm = TRUE)
  } else {
    rdata = data - minRating
  }

  # Range of person minus item measures for lookup table
  rngI = max(items, na.rm = TRUE) - min(items, na.rm = TRUE)
  maxDiff = 10
  maxLogit = ceiling(rngI) + maxDiff

  # Create MSD probability lookup table
  inc = 0.001
  gamma = seq(-maxLogit, maxLogit, inc); lenG = length(gamma)
  GM = matrix(NA, nrow = length(thresholds) + 1, ncol = lenG)
  for (gg in seq(1, dim(GM)[2])){
    GM[,gg] = msdprob(gamma[gg], thresholds)
  }

  # Estimate each person measure independently
  PM = NA*seq(1, dim(rdata)[1]); PE = NA*seq(1, dim(rdata)[1])
  len0 = 10
  for (pp in seq(1, dim(rdata)[1])){
    # Responses of current person
    R_0 = rdata[pp, !is.na(rdata[pp,])]

    # Item measures for current person
    IP_0 = items[!is.na(rdata[pp,])]

    # Remove NA
    R = R_0[!is.na(IP_0)]
    IP = IP_0[!is.na(IP_0)]

    # Check if there is anything to estimate
    if ((length(R) == 0) || (length(IP) == 0)){
      PM[pp] = NA
      PE[pp] = NA

    } else {
      # Maximum likelihood estimation of person measures
      pMeas = 0; minP = min(IP) - maxDiff; maxP = max(IP) + maxDiff
      LLmax = -Inf; SDerr = NA; loop = 1
      while (loop == 1){

        # Current range of candidate person measures
        CR = seq(pMeas + minP, pMeas + maxP, length.out = len0)

        # Calculate log likelihood at all candidate person measures
        for (cc in seq(1, length(CR))){

          # Log likelihood
          LL = sum(log(GM[R + 1 + nrow(GM)*
                  round(maxLogit/inc + (CR[cc] - IP)/inc)]))

          # Maximum log likelihood
          if (LL > LLmax){
            LLmax = LL
            pMeas = CR[cc]
            indC = cc
          }
        }

        # Check if estimate is at edge of lookup table
        if ((pMeas == (min(IP) - maxDiff)) ||
            (pMeas == (max(IP) + maxDiff))){
          pMeas = NA
          loop = 0

        } else {
          # Threshold for ending search
          inc0 = CR[2] - CR[1]
          if (inc0 < inc){

              # Calculate standard error by fitting polynomial
              sV = seq(-3*inc, 3*inc, inc); LLse = rep(NA, 7)
              for (ss in seq(1, length(sV))){
                LLse[ss] = sum(log(GM[R + 1 + nrow(GM)*
                              round(maxLogit/inc +
                                  ((pMeas + sV[ss]) - IP)/inc)]))
              }

            # Polynomial fit
            PF = lm(LLse ~ poly(sV, 2, raw = TRUE))

            # Standard error
            if (-PF$coefficients[[3]] > 0){
              SDerr = 1/sqrt(-PF$coefficients[[3]])
            } else {
              SDerr = NA
            }
            loop = 0

          } else {
            # Update range of candidate person measures
            minP = -inc0
            maxP = inc0
          }
        }
      }

      # Store results
      PM[pp] = pMeas
      PE[pp] = SDerr
    }
  }

  # Change non-finite numbers to NA
  PM[!is.finite(PM) | !is.finite(PE)] = NA
  PE[!is.finite(PM) | !is.finite(PE)] = NA

  # MISFIT STATISTICS --------------------------------------------------

  if (misfit == TRUE){
    if (is.null(minRating)){
      MSF = misfit(rdata, items, PM, thresholds)
    } else {
      MSF = misfit(rdata + minRating, items, PM, thresholds, minRating)
    }
  }

  # OUTPUT -------------------------------------------------------------
  out = c()
  out$person_measures = PM
  out$person_std_errors = PE

  if (misfit == TRUE){
    out$infit_persons = MSF$infit_persons
    out$outfit_persons = MSF$outfit_persons
  }

  return(out)
}



# Estimates MSD item measures given person measures and thresholds
ims <- function(data, persons, thresholds, misfit = FALSE,
                minRating = NULL){
  # Estimates item measures assuming all persons use the same rating
  # category thresholds.

  # REQUIRED INPUTS
  # 'data' = a numeric matrix of ordinal rating scale data whose entries
  #          are integers with missing data set to NA. Rows are persons
  #          and columns are items. The ordinal rating scale is assumed
  #          to go from the smallest to largest integer in integer steps
  #          unless minRating is specified in which case the ordinal
  #          rating scale goes from minRating to the largest integer in
  #          integer steps.
  # 'persons' = a numeric vector of person measures with missing values
  #             set to NA.
  # 'thresholds' = a numeric vector of ordered rating category thresholds
  #                with no NA.

  # OPTIONAL INPUTS
  # 'misfit' = logical for calculating infit and outfit statistics.
  #            Default is FALSE.
  # 'minRating' = integer representing the smallest ordinal rating
  #               category. Default is NULL.


  # CHECK FOR MISSING REQUIRED INPUTS ------------------------------------

  if (missing("data")){
    stop('Argument "data" is missing with no default')
  }
  if (missing("persons")){
    stop('Argument "persons" is missing with no default')
  }
  if (missing("thresholds")){
    stop('Argument "thresholds" is missing with no default')
  }

  # CHECK DATA STRUCTURES AND DATA TYPES ---------------------------------

  if (!is.matrix(data) || typeof(data) != "double"){
    stop('Argument "data" must be a matrix of mode numeric whose entries
         are either integers or NA.')
  }
  if (!is.vector(persons) || typeof(persons) != "double"){
    stop('Argument "persons" must be a vector of mode numeric with all
         unknown person measures set to NA')
  }
  if (!is.vector(thresholds) || typeof(thresholds) != "double"){
    stop('Argument "thresholds" must be an ordered vector of mode numeric
         with no NA')
  }
  if (!is.logical(misfit)){
    stop('Argument "misfit" must be a logical TRUE or FALSE')
  }
  if (!is.null(minRating) && length(minRating) != 1){
    stop('Argument "minRating" must be an integer')
  }
  if (!is.null(minRating) && minRating%%1 != 0){
    stop('Argument "minRating" must be an integer')
  }
  if (length(minRating) == 1 && typeof(minRating) == "integer"){
    minRating = as.numeric(minRating)
  }


  # CHECK VALUES AND DIMENSIONS ------------------------------------------

  # Check if non-NA entries of data are integers
  if (all(data%%1 == 0, na.rm = TRUE) == FALSE){
    stop('Entries in "data" must consist of integers or NA')}

  # Check if there are at least two different integers
  if (length(unique(c(data[!is.na(data)]))) < 2){
    stop('There must be at least two different integers in "data"')
  }

  # Compare dimensions of data with vectors items
  if (length(persons) != dim(data)[1]){
    stop('Argument "persons" must be a vector whose length is the
         number of rows in input matrix "data"')
  }

  # Check if there are any NA in thresholds and if thresholds is ordered
  if (any(is.na(thresholds)) ||
      identical(order(thresholds),
                seq(1, length(thresholds))) == FALSE){
    stop('Argument "thresholds" must be an ordered vector of real numbers
         with no NA')
  }

  # Check if the length of thresholds equals the maximum observed rating
  # in 'data' minus 0 if 'minRating' is not specified
  if (is.null(minRating) &&
      max(data, na.rm = TRUE) - min(data, na.rm = TRUE)
      != length(thresholds)){
    stop('The length of "thresholds" must equal the maximum minus
         minimum rating in "data" unless "minRating" is specified')
  }
  # Check if the maximum observed rating is too large given 'minRating'
  if (!is.null(minRating) &&
      max(data, na.rm = TRUE) - minRating > length(thresholds)){
    stop('The maximum rating in "data" cannot be larger than "minRating"
         plus the length of "thresholds"')
  }
  # Check if 'minRating' is too large
  if (!is.null(minRating) && minRating > min(data, na.rm = TRUE)){
    stop('"minRating" cannot be larger than the smallest integer in
         "data"')
  }


  # INITIALIZE -----------------------------------------------------------

  # Define lowest rating to be zero
  if (is.null(minRating)){
    rdata = data - min(data, na.rm = TRUE)
  } else {
    rdata = data - minRating
  }

  # Range of person minus item measures for lookup table
  rngP = max(persons, na.rm = TRUE) - min(persons, na.rm = TRUE)
  maxDiff = 10
  maxLogit = ceiling(rngP) + maxDiff

  # Create MSD probability lookup table
  inc = 0.001
  gamma = seq(-maxLogit, maxLogit, inc); lenG = length(gamma)
  GM = matrix(NA, nrow = length(thresholds) + 1, ncol = lenG)
  for (gg in seq(1, dim(GM)[2])){
    GM[,gg] = msdprob(gamma[gg], thresholds)
  }

  # Estimate each item measure independently
  IM = NA*seq(1, dim(rdata)[2]); IE = NA*seq(1, dim(rdata)[2])
  len0 = 10
  for (ii in seq(1, dim(rdata)[2])){
    # Responses to current item
    R_0 = rdata[!is.na(rdata[,ii]), ii]

    # Person measures for current item
    PI_0 = persons[!is.na(rdata[,ii])]

    # Remove NA
    R = R_0[!is.na(PI_0)]
    PI = PI_0[!is.na(PI_0)]

    # Check if there is anything to estimate
    if ((length(R) == 0) || (length(PI) == 0)){
      IM[ii] = NA
      IE[ii] = NA

    } else {
      # Maximum likelihood estimation of item measures
      iMeas = 0; minI = min(PI) - maxDiff; maxI = max(PI) + maxDiff
      LLmax = -Inf; SDerr = NA; loop = 1
      while (loop == 1){

        # Current range of candidate item measures
        CR = seq(iMeas + minI, iMeas + maxI, length.out = len0)

        # Calculate log likelihood at all candidate item measures
        for (cc in seq(1, length(CR))){

          # Log likelihood
          LL = sum(log(GM[R + 1 + nrow(GM)*
                  round(maxLogit/inc + (PI - CR[cc])/inc)]))

          # Maximum log likelihood
          if (LL > LLmax){
            LLmax = LL
            iMeas = CR[cc]
            indC = cc
          }
        }

        # Check if estimate is at edge of lookup table
        if ((iMeas == (min(PI) - maxDiff)) ||
            (iMeas == (max(PI) + maxDiff))){
          iMeas = NA
          loop = 0

        } else {
          # Threshold for ending search
          inc0 = CR[2] - CR[1]
          if (inc0 < inc){

            # Calculate standard error by fitting polynomial
            sV = seq(-3*inc, 3*inc, inc); LLse = rep(NA, 7)
            for (ss in seq(1, length(sV))){
              LLse[ss] = sum(log(GM[R + 1 + nrow(GM)*
                          round(maxLogit/inc +
                              (PI - (iMeas + sV[ss]))/inc)]))
            }

            # Polynomial fit
            PF = lm(LLse ~ poly(sV, 2, raw = TRUE))

            # Standard error
            if (-PF$coefficients[[3]] > 0){
              SDerr = 1/sqrt(-PF$coefficients[[3]])
            } else {
              SDerr = NA
            }
            loop = 0

          } else {
            # Update range of candidate item measures
            minI = -inc0
            maxI = inc0
          }
        }
      }

      # Store results
      IM[ii] = iMeas
      IE[ii] = SDerr
    }
  }

  # Change non-finite numbers to NA
  IM[!is.finite(IM) | !is.finite(IE)] = NA
  IE[!is.finite(IM) | !is.finite(IE)] = NA

  # MISFIT STATISTICS --------------------------------------------------

  if (misfit == TRUE){
    if (is.null(minRating)){
      MSF = misfit(rdata, IM, persons, thresholds)
    } else {
      MSF = misfit(rdata + minRating, IM, persons, thresholds, minRating)
    }
  }

  # OUTPUT -------------------------------------------------------------
  out = c()
  out$item_measures = IM
  out$item_std_errors = IE

  if (misfit == TRUE){
    out$infit_items = MSF$infit_items
    out$outfit_items = MSF$outfit_items
  }

  return(out)
}



# MSD probabilities
msdprob <- function(x, thresholds){
  # Estimates the probability of observing each rating category.

  # REQUIRED INPUTS
  # 'x' = a numeric vector of real numbers representing a set of person
  #       minus item measures with no NA.
  # 'thresholds' = a numeric vector of ordered rating category thresholds
  #                with no NA.

  # CHECK FOR MISSING REQUIRED INPUTS ------------------------------------

  if (missing("x")){
    stop('Argument "x" is missing with no default')
  }
  if (missing("thresholds")){
    stop('Argument "thresholds" is missing with no default')
  }

  # CHECK DATA STRUCTURES AND DATA TYPES ---------------------------------

  if (!is.vector(x) || typeof(x) != "double"){
    stop('Argument "x" must be a vector of mode numeric with no NA')
  }
  if (!is.vector(thresholds) || typeof(thresholds) != "double"){
    stop('Argument "thresholds" must be an ordered vector of
         mode numeric with no NA')
  }

  # CHECK VALUES AND DIMENSIONS ------------------------------------------

  # Check if there are any NA in x
  if (any(is.na(x))){
    stop('Argument "x" must be a vector of mode numeric with no NA')
  }
  # Check if there are any NA in thresholds and if thresholds is ordered
  if (any(is.na(thresholds)) ||
      identical(order(thresholds),
                seq(1,length(thresholds))) == FALSE){
    stop('Argument "thresholds" must be an ordered vector of real numbers
         with no NA')
  }

  # INITIALIZE -----------------------------------------------------------

  PR = matrix(NA, length(thresholds)+1, length(x))
  for (gg in seq(0, length(x))){
    for (rr in seq(0, length(thresholds))){
      if (rr==0){                          # lowest rating
        PR[rr+1, gg] = 1/(1+exp(x[gg] - thresholds[rr+1]))
      } else if (rr==length(thresholds)){  # highest rating
        PR[rr+1, gg] = 1 - 1/(1+exp(x[gg] - thresholds[rr]))
      } else {
        PR[rr+1, gg] = 1/(1+exp(x[gg] - thresholds[rr+1])) -
          1/(1+exp(x[gg] - thresholds[rr]))
      }
    }
  }

  # Output
  return(PR)
}



# Expected ratings
expdata <- function(items, persons, thresholds, minRating){
  # Expected ratings matrix given item measures, person measures and
  # ordered rating category thresholds.

  # REQUIRED INPUTS
  # 'items' = a numeric vector of item measures with missing values set
  #           to NA.
  # 'persons' = a numeric vector of person measures with missing values
  #             set to NA.
  # 'thresholds' = a numeric vector of ordered rating category thresholds
  #                with no NA.
  # 'minRating' = integer representing the smallest ordinal rating
  #               category.


  # CHECK FOR MISSING REQUIRED INPUTS ------------------------------------

  if (missing("items")){
    stop('Argument "items" is missing with no default')
  }
  if (missing("persons")){
    stop('Argument "persons" is missing with no default')
  }
  if (missing("thresholds")){
    stop('Argument "thresholds" is missing with no default')
  }
  if (missing("minRating")){
    stop('Argument "minRating" is missing with no default')
  }

  # CHECK DATA STRUCTURES AND DATA TYPES ---------------------------------

  if (!is.vector(items) || typeof(items) != "double"){
    stop('Input "items" must be a vector of mode numeric with all
         unknown item measures set to NA')
  }
  if (!is.vector(persons) || typeof(persons) != "double"){
    stop('Input "persons" must be a vector of mode numeric with all
         unknown person measures set to NA')
  }
  if (!is.vector(thresholds) || typeof(thresholds) != "double"){
    stop('Input "thresholds" must be an ordered vector of mode numeric
         with no NA')
  }
  if (length(minRating) == 1 && typeof(minRating) == "integer"){
    minRating = as.numeric(minRating)
  }
  if (length(minRating) != 1 || typeof(minRating) != "double" ||
      minRating%%1 != 0){
    stop('Argument "minRating" must be an integer')
  }

  # CHECK VALUES AND DIMENSIONS ------------------------------------------

  # Check if there are any NA in thresholds and if thresholds is ordered
  if (any(is.na(thresholds)) ||
      identical(order(thresholds),seq(1,length(thresholds))) == FALSE){
    stop('Input "thresholds" must be an ordered vector of real numbers
         with no NA')
  }

  # INITIALIZE -----------------------------------------------------------

  D = matrix(NA, nrow=length(persons), ncol=length(items))
  RT = as.numeric(seq(minRating, minRating+length(thresholds)))

  for (pp in seq(1,length(persons))){
    for (ii in seq(1,length(items))){
      MR = msdprob(persons[pp]-items[ii],thresholds)
      D[pp,ii] = sum(RT*MR)
    }
  }

  # Output
  return(D)
}



# Misfit statistics (infit and outfit)
misfit <- function(data, items, persons, thresholds, minRating = NULL){
  # Calculates infit and outfit for items and persons.

  # REQUIRED INPUTS
  # 'data' = a numeric matrix of ordinal rating scale data whose entries
  #          are integers with missing data set to NA. Rows are persons
  #          and columns are items. The ordinal rating scale is assumed
  #          to go from the smallest to largest integer in integer steps
  #          unless minRating is specified in which case the ordinal
  #          rating scale goes from minRating to the largest integer in
  #          integer steps.
  # 'items' = a numeric vector of item measures with missing values set
  #           to NA.
  # 'persons' = a numeric vector of person measures with missing values
  #             set to NA.
  # 'thresholds' = a numeric vector of ordered rating category thresholds
  #                with no NA.

  # OPTIONAL INPUTS
  # 'minRating' = integer representing the smallest ordinal rating
  #               category. Default is NULL.

  # CHECK FOR MISSING REQUIRED INPUTS ------------------------------------

  if (missing("data")){
    stop('Argument "data" is missing with no default')
  }
  if (missing("items")){
    stop('Argument "items" is missing with no default')
  }
  if (missing("persons")){
    stop('Argument "persons" is missing with no default')
  }
  if (missing("thresholds")){
    stop('Argument "thresholds" is missing with no default')
  }

  # CHECK DATA STRUCTURES AND DATA TYPES ---------------------------------

  if (!is.matrix(data) || typeof(data) != "double"){
    stop('Input "data" must be a matrix of mode numeric or integers
         with missing values set to NA')
  }
  if (!is.vector(items) || typeof(items) != "double"){
    stop('Input "items" must be a vector of mode numeric with all
         non-fixed item measures set to NA')
  }
  if (!is.vector(persons) || typeof(persons) != "double"){
    stop('Input "persons" must be a vector of mode numeric with all
         non-fixed person measures set to NA')
  }
  if (!is.vector(thresholds) || typeof(thresholds) != "double"){
    stop('Input "thresholds" must be a vector of mode numeric with no NA')
  }
  if (!is.null(minRating) && length(minRating) != 1){
    stop('Argument "minRating" must be an integer')
  }
  if (!is.null(minRating) && minRating%%1 != 0){
    stop('Argument "minRating" must be an integer')
  }
  if (length(minRating) == 1 && typeof(minRating) == "integer"){
    minRating = as.numeric(minRating)
  }

  # CHECK VALUES AND DIMENSIONS ------------------------------------------

  # Check if non-NA entries of data are integers
  if (all(data%%1 == 0, na.rm = TRUE) == FALSE){
    stop('Input "data" must consist of integers and NA only')}

  # Compare dimensions of data with vectors items and persons
  if (length(items) != dim(data)[2]){
    stop('Input "items" must be a vector whose length is the
         number of columns in input matrix "data"')
  }
  if (length(persons) != dim(data)[1]){
    stop('Input "persons" must be a vector whose length is the
         number of rows in input matrix "data"')
  }

  # Check if there are any NA in thresholds and if thresholds is ordered
  if (any(is.na(thresholds)) ||
      identical(order(thresholds),
                seq(1,length(thresholds))) == FALSE){
    stop('Input "thresholds" must be an ordered vector of real numbers
         with no NA')
  }

  # Check if the length of thresholds equals the maximum observed rating
  # in 'data' minus 0 if 'minRating' is not specified
  if (is.null(minRating) &&
      max(data, na.rm = TRUE) - min(data, na.rm = TRUE)
      != length(thresholds)){
    stop('The length of "thresholds" must equal the maximum minus
         minimum rating in "data" unless "minRating" is specified')
  }
  # Check if the maximum observed rating is too large given 'minRating'
  if (!is.null(minRating) &&
      max(data, na.rm = TRUE) - minRating > length(thresholds)){
    stop('The maximum rating in "data" cannot be larger than "minRating"
         plus the length of "thresholds"')
  }
  # Check if 'minRating' is too large
  if (!is.null(minRating) && minRating > min(data, na.rm = TRUE)){
    stop('"minRating" cannot be larger than the smallest integer in
         "data"')
  }

  # INITIALIZE -----------------------------------------------------------

  # Define lowest rating to be zero
  if (is.null(minRating)){
    rdata = data - min(data, na.rm = TRUE)
  } else {
    rdata = data - minRating
  }

  # Define highest rating category = length of thresholds
  R_max = length(thresholds)

  RV = seq(0, R_max)
  IF_pN = rep(0, length(persons)); IF_pD = IF_pN; OF_p = IF_pD
  IF_iN = rep(0, length(items)); IF_iD = IF_iN; OF_i = IF_iD

  # Loop through all persons and items
  for (pp in seq(1,length(persons))){
    for (ii in seq(1,length(items))){
      if (!is.na(rdata[pp,ii]) && !is.na(persons[pp]) &&
          !is.na(items[ii])){
        # Expected frequencies of responses using MSD
        ER = msdprob(persons[pp]-items[ii], thresholds)

        # Check if person or item is NA
        if (!anyNA(ER)){
          # Add to infit
          ExR = sum(ER*RV, na.rm = TRUE)
          ExR2 = sum(ER*RV^2, na.rm = TRUE)
          IF_n = (rdata[pp,ii] - ExR)^2
          IF_d = ExR2 - ExR^2

          # Infit for persons and items
          IF_pN[pp] = IF_pN[pp] + IF_n
          IF_pD[pp] = IF_pD[pp] + IF_d
          IF_iN[ii] = IF_iN[ii] + IF_n
          IF_iD[ii] = IF_iD[ii] + IF_d

          # Add to output
          OF_pi = IF_n/IF_d

          # Output
          OF_p[pp] = OF_p[pp] + OF_pi
          OF_i[ii] = OF_i[ii] + OF_pi
        }
      }
    }
  }

  # Final infit and outfit
  Infit_I = IF_iN/IF_iD
  Outfit_I = OF_i/colSums(!is.na(rdata))
  Infit_P = IF_pN/IF_pD
  Outfit_P = OF_p/rowSums(!is.na(rdata))

  # Change non-finite numbers to NA
  Infit_I[!is.finite(Infit_I)] = NA
  Outfit_I[!is.finite(Outfit_I)] = NA
  Infit_P[!is.finite(Infit_P)] = NA
  Outfit_P[!is.finite(Outfit_P)] = NA

  # Force NA if item or person measure is NA
  Infit_I[is.na(items)] = NA
  Outfit_I[is.na(items)] = NA
  Infit_P[is.na(persons)] = NA
  Outfit_P[is.na(persons)] = NA

  # Output
  out = c()
  out$infit_items = Infit_I
  out$infit_persons = Infit_P
  out$outfit_items = Outfit_I
  out$outfit_persons = Outfit_P

  return(out)
}



# Simulated rating scale data
simdata <- function(items, persons, thresholds, missingProb = 0,
                    minRating = 0){
  # Simulated ratings matrix given item measures, person measures and
  # ordered rating category thresholds.

  # REQUIRED INPUTS
  # 'items' = a numeric vector of item measures with no NA.
  # 'persons' = a numeric vector of person measures with no NA.
  # 'thresholds' = a numeric vector of ordered rating category thresholds
  #                with no NA.

  # OPTIONAL INPUTS
  # 'missingProb' = a real number between 0 and 1 specifying the
  #                 probability of missing data (i.e., NA). Default is 0.
  # 'minRating' = integer representing the smallest ordinal rating
  #               category. Default is 0.

  # CHECK FOR MISSING REQUIRED INPUTS ------------------------------------

  if (missing("items")){
    stop('Argument "items" is missing with no default')
  }
  if (missing("persons")){
    stop('Argument "persons" is missing with no default')
  }
  if (missing("thresholds")){
    stop('Argument "thresholds" is missing with no default')
  }

  # CHECK DATA STRUCTURES AND DATA TYPES ---------------------------------

  if (!is.vector(items) || typeof(items) != "double"){
    stop('Argument "items" must be a vector of mode numeric with no NA')
  }
  if (!is.vector(persons) || typeof(persons) != "double"){
    stop('Argument "persons" must be a vector of mode numeric with no NA')
  }
  if (!is.vector(thresholds) || typeof(thresholds) != "double"){
    stop('Argument "thresholds" must be an ordered vector of
         mode numeric with no NA')
  }
  if (length(missingProb) == 1 && typeof(missingProb) == "integer"){
    missingProb = as.numeric(missingProb)
  }
  if (length(missingProb) != 1 || typeof(missingProb) != "double"){
    stop('Argument "missingProb" must be a real number between 0 and 1')
  }
  if (length(minRating) != 1 || typeof(minRating) != "double" ||
      minRating%%1 != 0){
    stop('Argument "minRating" must be an integer')
  }
  if (length(minRating) == 1 && typeof(minRating) == "integer"){
    minRating = as.numeric(minRating)
  }


  # CHECK VALUES AND DIMENSIONS ------------------------------------------

  # Check if there are any NA in items
  if (any(is.na(items))){
    stop('Argument "items" must be a vector of mode numeric with no NA')
  }
  # Check if there are any NA in persons
  if (any(is.na(items))){
    stop('Argument "persons" must be a vector of mode numeric with no NA')
  }
  # Check if there are any NA in thresholds and if thresholds is ordered
  if (any(is.na(thresholds)) ||
      identical(order(thresholds),
                seq(1, length(thresholds))) == FALSE){
    stop('Argument "thresholds" must be an ordered vector of real numbers
         with no NA')
  }
  # Check if optional argument 'missing' is less than 0 or greater than 1
  if (missingProb < 0 || missingProb > 1){
    stop('Argument "missingData" must be a real number between 0 and 1')
  }

  # INITIALIZE -----------------------------------------------------------

  # Output matrix with defaults set to NA
  M = matrix(NA, nrow = length(persons), ncol = length(items))

  # Loop through items and persons to simulate rating
  for (ii in seq(1, length(items))){
    for (pp in seq(1, length(persons))){

      # Check if NA
      if (!is.na(items[ii]) && !is.na(persons[pp])){

        # Probability of observing each rating category
        prb = msdprob(persons[pp] - items[ii], thresholds)

        # Randomly select rating category (lowest rating = 0)
        prb2 = c(0, cumsum(prb))
        M[pp, ii] = findInterval(runif(1), prb2) - 1.0

        # Probability of missing data (probability of NA)
        if (missingProb > 0){
          if (runif(1) < missingProb){
            M[pp, ii] = NA
          }
        }
      }
    }
  }

  # Set minimum rating category
  M = M + minRating

  # Output
  return(M)
}
