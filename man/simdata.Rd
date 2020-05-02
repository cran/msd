\name{simdata}
\alias{simdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulated Rating Scale Data
}
\description{
Generates simulated rating scale data given item measures, person measures and
rating category thresholds.
}
\usage{
simdata(items, persons, thresholds, missingProb = 0, minRating = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{items}{a numeric vector of item measures with no NA.}
  \item{persons}{a numeric vector of person measures with no NA.}
  \item{thresholds}{a numeric vector of ordered rating category thresholds with no NA.}
  \item{missingProb}{a number between 0 and 1 specifying the probability of
  missing data.}
  \item{minRating}{integer representing the smallest ordinal rating category.
  Default is 0 (see Details).}
}
\details{
It is assumed that the set of ordinal rating categories consists of all integers from the lowest rating category specified by \code{minRating} to the highest rating category,
which is \code{minRating + length(thresholds)}.
}
\value{
A numeric matrix of simulated rating scale data.
}
\author{
Chris Bradley (cbradley05@gmail.com)
}
\note{
\code{simdata} can be used to test the accuracy of \code{msd} (see Examples).
}
\seealso{
\code{\link{msd}}
}
\examples{
# Using randomly generated item measures, person measures and thresholds with
# 15 percent missing data and ordinal rating categories from 0 to 5
im = runif(100, -2, 2)
pm = runif(200, -2, 2)
th = sort(runif(5, -2, 2))
d1 = simdata(im, pm, th, missingProb = 0.15, minRating = 0)

# To test the accuracy of msd using simdata, set the mean item measure to zero
# (axis origin in msd is the mean item measure) and the mean threshold to
# zero (any non-zero mean threshold is reflected in the person measures).
im = runif(100, -2, 2)
im = im - mean(im)
pm = runif(200, -2, 2)
th = sort(runif(5, -2, 2))
th = th - mean(th)
d2 = simdata(im, pm, th, missingProb = 0.15, minRating = 0)
m = msd(d2)

# Compare msd parameters to true values.  Linear regression should
# yield a slope very close to 1 and an intercept very close to 0.
lm(m$item_measures ~ im)
lm(m$person_measures ~ pm)
lm(m$thresholds ~ th)
}
\keyword{ datagen }% use one of  RShowDoc("KEYWORDS")
