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

  # Use dichotomous Rasch model if only two rating categories
  if (R_max == 1){

    # OUTPUT -------------------------------------------------------------
    out = rasch(rdata, items, persons, misfit)

  } else {
    # Find all rows and columns with only one rating and replace with NA
    X_row = c(); X_col = c()
    UniqRow = apply(rdata, 1, function(x)length(unique(x)))
    UniqCol = apply(rdata, 2, function(x)length(unique(x)))
    X_row = c(X_row, which(UniqRow==1))
    X_col = c(X_col, which(UniqCol==1))
    MM=rdata; MM[X_row,]=NA; MM[,X_col]=NA

    # Initialize person and item measure matrices for dichotomizations
    PM = PE = matrix(NA, dim(rdata)[1], R_max)
    IM = IE = matrix(NA, R_max, dim(rdata)[2])

    # PERSON AND ITEM MEASURES -------------------------------------------

    # Loop through dichotomizations to estimate person and item measures
    for (dd in seq(1, R_max)){
      # Dichotomize ratings matrix
      M2 = MM; M2[M2<dd] = 0; M2[M2>=dd] = 1

      # Find rows and columns with all 1's or all 0's after removing NA
      X2_row = X_row; X2_col = X_col; new_RC = 1
      while(new_RC==1){
        # Current length of rows and columns to be removed
        lnR = length(X2_row); lnC = length(X2_col)

        # Find indexes of all 1's or all 0's rows and columns ignoring NA
        UniqRow = apply(M2, 1, function(x)length(unique(x[!is.na(x)])))
        UniqCol = apply(M2, 2, function(x)length(unique(x[!is.na(x)])))
        X2_row = unique(c(X2_row, which(UniqRow==1)))
        X2_col = unique(c(X2_col, which(UniqCol==1)))

        # Set all rows and columns to delete as NA
        M2[X2_row,]=NA; M2[,X2_col]=NA

        # When to stop iterating
        if (length(X2_row)==lnR && length(X2_col)==lnC){
          new_RC = 0
        }
      }

      # Remove all undesired rows and columns
      fixI=items; fixP=persons
      if (length(X2_row)>0){
        M2 = M2[-X2_row,]
        fixP = fixP[-X2_row]
      }
      if (length(X2_col)>0){
        M2 = M2[,-X2_col]
        fixI = fixI[-X2_col]
      }

      # CHECK IF MATRIX CONFORMS TO GUTTMAN SCALE -----------------------

      # Check if one of the dimensions is zero
      if (dim(M2)[1]==0 || dim(M2)[2]==0){
        stop('MSD cannot estimate parameters for this input matrix "data"
        because at least one dichotomization of the matrix conforms
             to a deterministic Guttman scale')
      }

      # DICHOTOMOUS RASCH MODEL -----------------------------------------

      # Dichotomous Rasch model
      RD = rasch(M2, fixI, fixP)

      # Indexes of estimated person and item measures
      Y2_row = setdiff(seq(1,dim(rdata)[1]), X2_row)
      Y2_col = setdiff(seq(1,dim(rdata)[2]), X2_col)

      # Person and item measures
      PM[Y2_row,dd] = RD$person_measures
      IM[dd,Y2_col] = RD$item_measures

      # Standard errors
      PE[Y2_row,dd] = RD$person_std_errors
      IE[dd,Y2_col] = RD$item_std_errors
    }

    # Impute person measures if estimate is NA for some dichotomization
    PM_m = colMeans(PM, na.rm=TRUE)
    PM_m2 = matrix(PM_m, nrow=dim(PM)[1], ncol=dim(PM)[2], byrow=TRUE)
    PM_dif = PM - PM_m2
    PM_dif_m = rowMeans(PM_dif, na.rm=TRUE)
    PM_dif_m2 = matrix(PM_dif_m, nrow=dim(PM)[1], ncol=dim(PM)[2])
    PM_mdif = PM_m2 + PM_dif_m2
    PM_imp = PM;
    PM_imp[is.na(PM_imp)] = PM_mdif[is.na(PM_imp)]

    # Final person and item measures
    IM_F = colMeans(IM, na.rm=TRUE)
    PM_F = rowMeans(PM_imp, na.rm=TRUE)
    if (mZero == 1){
      IM_F = IM_F - mean(colMeans(IM, na.rm=TRUE), na.rm=TRUE)
      PM_F = PM_F - mean(colMeans(IM, na.rm=TRUE), na.rm=TRUE)
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
    PSE = sqrt(rowSums(PE^2, na.rm=TRUE)/R_max)/sqrt(R_max)
    ISE = sqrt(colSums(IE^2, na.rm=TRUE)/R_max)/sqrt(R_max)

    # Change to NA if person or item measure is NA
    PSE[is.na(PM_F)] = NA
    ISE[is.na(IM_F)] = NA

    # RELIABILITY --------------------------------------------------------

    # Reliability measures for items and persons
    RL_items = 1 - (mean(ISE, na.rm=TRUE))^2/(sd(IM_F, na.rm=TRUE))^2
    RL_persons = 1 - (mean(PSE, na.rm=TRUE))^2/(sd(PM_F, na.rm=TRUE))^2

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

  # Determine if mean item measure should be set to zero
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
  if (mZero == 1){
    # Set mean item measure to zero
    PII_logit = PII_logit - mean(PII_logit, na.rm = TRUE)
  }

  # Fix item and person measures
  PII_logit[!is.na(items)] = items[!is.na(items)]
  PCP_logit[!is.na(persons)] = persons[!is.na(persons)]

  # Iterate until residuals (res) is close enough to zero
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
  PII_logit[!is.finite(PII_logit)] = NA
  PCP_logit[!is.finite(PCP_logit)] = NA
  SE_items[!is.finite(SE_items)] = NA
  SE_persons[!is.finite(SE_persons)] = NA

  # Change to NA if person or item measure is NA
  SE_items[is.na(PII_logit)] = NA
  SE_persons[is.na(PCP_logit)] = NA

  # Reliability measures for items and persons
  RL_items = 1 - (mean(SE_items, na.rm=TRUE))^2/
    (sd(PII_logit, na.rm=TRUE))^2
  RL_persons = 1 - (mean(SE_persons, na.rm=TRUE))^2/
    (sd(PCP_logit, na.rm=TRUE))^2

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

  # Compare dimensions of data with vectors items
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

  # Maximum absolute value person measure
  maxLogit = 10

  # Define lowest rating to be zero
  if (is.null(minRating)){
    rdata = data - min(data, na.rm = TRUE)
  } else {
    rdata = data - minRating
  }

  # Create MSD probability lookup table
  inc = 0.001; gamma = seq(-2.2*maxLogit,2.2*maxLogit,inc)
  GM = matrix(NA, length(thresholds)+1, length(gamma))
  for (gg in seq(1,dim(GM)[2])){
    GM[,gg] = msdprob(gamma[gg],thresholds)
  }

  # Estimate each person measure independently
  PM = NA*seq(1,dim(rdata)[1]); PE = NA*seq(1,dim(rdata)[1])
  for (pp in seq(1,dim(rdata)[1])){
    # Responses of current person
    R_0 = rdata[pp,!is.na(rdata[pp,])]

    # Item measures for current person
    IP_0 = items[!is.na(rdata[pp,])]

    # Remove NA
    R = R_0[!is.na(IP_0)]
    IP = IP_0[!is.na(IP_0)]

    # Minimize negative of likelihood function with x as person measure
    f <- function(x) -sum(log(GM[R+1
                          + nrow(GM)*(round((2.2*maxLogit)*(1/inc)+1) +
                          round((x-IP)*(1/inc))-1)]))

    # Estimate person measure
    P_opt = optim(0, f, method = "Brent" , lower = -1.1*maxLogit,
                  upper = 1.1*maxLogit, hessian = TRUE)
    PM[pp] = P_opt$par
    PE[pp] = 1/sqrt(P_opt$hessian)
  }

  # Set any estimates beyond cut-off point to NA
  PE[PM <= -maxLogit] = NA
  PE[PM >= maxLogit] = NA
  PM[PM <= -maxLogit] = NA
  PM[PM >= maxLogit] = NA

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



# MSD probabilities
msdprob <- function(x, thresholds){
  # Estimates the probability of observing each rating category.

  # REQUIRED INPUTS
  # 'x' = a real number representing a person minus item measure.
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

  if (length(x) == 1 && typeof(x) == "integer"){
    x = as.numeric(x)
  }
  if (length(x) != 1 || typeof(x) != "double"){
    stop('Argument "x" must be a real number')
  }
  if (!is.vector(thresholds) || typeof(thresholds) != "double"){
    stop('Argument "thresholds" must be an ordered vector of
         mode numeric with no NA')
  }

  # CHECK VALUES AND DIMENSIONS ------------------------------------------

  # Check if there are any NA in thresholds and if thresholds is ordered
  if (any(is.na(thresholds)) ||
      identical(order(thresholds),
      seq(1,length(thresholds))) == FALSE){
    stop('Argument "thresholds" must be an ordered vector of real numbers
         with no NA')
  }

  # INITIALIZE -----------------------------------------------------------

  PR = c()
  for (rr in seq(0,length(thresholds))){
    if (rr==0){                          # lowest rating
      PR[rr+1] = 1/(1+exp(x - thresholds[rr+1]))
    } else if (rr==length(thresholds)){  # highest rating
      PR[rr+1] = 1 - 1/(1+exp(x - thresholds[rr]))
    } else {
      PR[rr+1] = 1/(1+exp(x - thresholds[rr+1])) -
        1/(1+exp(x - thresholds[rr]))
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
      if (!is.na(rdata[pp,ii])){
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
