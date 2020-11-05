\name{rasch}
\alias{rasch}
\title{
Dichotomous Rasch Model
}
\description{
Estimates item measures, person measures and their standard errors using the dichotomous Rasch model.  A special case of the function \code{msd} when the rating scale consists of only two rating categories: 0 and 1.  Option provided for anchoring certain items and persons while estimating the rest.  Option also provided for estimating infit and outfit statistics.
}
\usage{
rasch(data, items = NULL, persons = NULL, misfit = FALSE)
}
\arguments{
  \item{data}{a numeric matrix of 0's and 1's with missing data set to NA.  Rows are persons and columns are items.}
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
  \item{item_std_errors}{a vector of standard errors for the items}
  \item{person_std_errors}{a vector of standard errors for the persons}
  \item{item_reliability}{reliability value for the items}
  \item{person_reliability}{reliability value for the persons}
  \item{infit_items}{if \code{misfit} = TRUE, a vector of infit statistics for the items}
  \item{outfit_items}{if \code{misfit} = TRUE, a vector of outfit statistics for the items}
  \item{infit_persons}{if \code{misfit} = TRUE, a vector of infit statistics for the persons}
  \item{outfit_persons}{if \code{misfit} = TRUE, a vector of outfit statistics for the persons}
}
\author{
Chris Bradley (cbradley05@gmail.com)
}
\note{
The axis origin is set by convention at the mean item measure.  All item measures and person measures that cannot be estimated will return as NA (e.g., if a person responds with a single rating category to all items, that person's person measure cannot be estimated).

\code{rasch} is the basis for the "successive dichotomizations" in \code{msd} and is repeatedly called by \code{msd} when there are three or more rating categories.

The accuracy of \code{rasch} can be tested using the \code{simdata} function (see Examples).
}
\seealso{
\code{\link{msd}}
\code{\link{simdata}}
}
\examples{
# Simple example using a randomly generated ratings matrix
d <- as.numeric(sample(0:1, 200, replace = TRUE))
dm <- matrix(d, nrow = 20, ncol = 10)
m1 <- rasch(dm, misfit = TRUE)

# Anchor first 5 item measures and first 10 person measures
im <- m1$item_measures
im[6:length(im)] <- NA
pm <- m1$person_measures
pm[11:length(pm)] <- NA
m2 <- rasch(dm, items = im, persons = pm)

# To test the accuracy of rasch using simdata, set the true mean item measure to
# zero (axis origin in rasch is the mean item measure).  Note that the threshold for
# dichotomous data is at 0.
im <- runif(100, -2, 2)
im <- im - mean(im)
pm <- runif(100, -2, 2)
th <- 0
d <- simdata(im, pm, th, missingProb = 0.15, minRating = 0)
m <- rasch(d)

# Compare rasch parameters to true values.  Linear regression should
# yield a slope very close to 1 and an intercept very close to 0.
lm(m$item_measures ~ im)
lm(m$person_measures ~ pm)
}
