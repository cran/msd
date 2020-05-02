\name{msd}
\alias{msd}
\title{
Method of Successive Dichotomizations
}
\description{
Estimates item measures, person measures, rating category thresholds and their standard errors using the method of successive dichotomizations.  Option provided for anchoring certain items and persons while estimating the rest.  Option also provided for estimating infit and outfit statistics.
}
\usage{
msd(data, items = NULL, persons = NULL, misfit = FALSE)
}
\arguments{
  \item{data}{a numeric matrix of ordinal rating scale data whose entries are integers with missing data set to NA.  Rows are persons and columns are items.  The ordinal rating scale is assumed to go from the smallest integer to the largest integer in \code{data} in integer steps.}
  \item{items}{a numeric vector of anchored item measures.  Item measures to be estimated are set to NA.  Default is NULL (see Details).}
  \item{persons}{a numeric vector of anchored person measures.  Person measures to be estimated are set to NA.  Default is NULL (see Details).}
  \item{misfit}{logical for calculating infit and outfit statistics.  Default is FALSE.}
}
\details{
\code{items} and \code{persons} are optional numeric vectors that specify item and person measures that should be "anchored" and not estimated.  The length of \code{items} must equal the number of columns in \code{data} and the length of \code{persons} must equal the number of rows in \code{data}.  Only entries set to NA in \code{items} and \code{persons} are estimated.  Default for both \code{items} and \code{persons} is NULL, which is equivalent to a vector of NA so that all items and persons are estimated.
}
\value{
A list whose elements are:
  \item{item_measures}{a vector of item measures for each item}
  \item{person_measures}{a vector of person measures for each person}
  \item{thresholds}{a vector of average rating category thresholds used by the persons when rating the items}
  \item{item_std_errors}{a vector of standard errors for the items}
  \item{person_std_errors}{a vector of standard errors for the persons}
  \item{threshold_std_errors}{a vector of standard errors for the thresholds}
  \item{item_reliability}{reliability of the item measures}
  \item{person_reliability}{reliability of the person measures}
  \item{infit_items}{if \code{misfit} = TRUE, a vector of infit statistics for the items}
  \item{outfit_items}{if \code{misfit} = TRUE, a vector of outfit statistics for the items}
  \item{infit_persons}{if \code{misfit} = TRUE, a vector of infit statistics for the persons}
  \item{outfit_persons}{if \code{misfit} = TRUE, a vector of outfit statistics for the persons}
}
\references{
Bradley, C. and Massof, R. W. (2018) Method of successive dichotomizations: An improved method for estimating measures of latent variables from rating scale data. \emph{PLoS One}, \bold{13(10)} doi:10.1371/journal.pone.0206106
}
\author{
Chris Bradley (cbradley05@gmail.com)
}
\note{
The axis origin is set by convention at the mean item measure.  All item measures and person measures that cannot be estimated will return as NA (e.g., if a person responds with a single rating category to all items, that person's person measure cannot be estimated).

The accuracy of \code{msd} can be tested using the \code{simdata} function (see Examples).
}
\seealso{
\code{\link{simdata}}
}
\examples{
# Simple example using a randomly generated ratings matrix
d <- as.numeric(sample(0:5, 1000, replace = TRUE))
dm <- matrix(d, nrow = 50, ncol = 20)
m1 <- msd(dm, misfit = TRUE)

# Anchor first 10 item measures and first 10 person measures
im <- m1$item_measures
im[11:length(im)] <- NA
pm <- m1$person_measures
pm[11:length(pm)] <- NA
m2 <- msd(dm, items = im, persons = pm, misfit = TRUE)

# To test the accuracy of msd using simdata, set the mean item measure to zero
# (axis origin in msd is the mean item measure) and the mean threshold to
# zero (any non-zero mean threshold is reflected in the person measures).
im = runif(100, -2, 2)
im = im - mean(im)
pm = runif(200, -2, 2)
th = sort(runif(5, -2, 2))
th = th - mean(th)
d = simdata(im, pm, th, missingProb = 0.15, minRating = 0)
m = msd(d)

# Compare msd parameters to true values.  Linear regression should
# yield a slope very close to 1 and an intercept very close to 0.
lm(m$item_measures ~ im)
lm(m$person_measures ~ pm)
lm(m$thresholds ~ th)
}
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
