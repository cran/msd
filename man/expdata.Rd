\name{expdata}
\alias{expdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expected Ratings Matrix
}
\description{
Expected ratings matrix given item measures, person measures and ordered rating category thresholds.
}
\usage{
expdata(items, persons, thresholds, minRating)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{items}{a numeric vector of item measures with missing values set to NA.}
  \item{persons}{a numeric vector of person measures with missing values set to NA.}
  \item{thresholds}{a numeric vector of ordered rating category thresholds with no NA.}
  \item{minRating}{integer representing the smallest ordinal rating category (see Details).}
}
\details{
It is assumed that the set of ordinal rating categories consists of all integers from the lowest rating category specified by \code{minRating} to the highest rating category, which is \code{minRating + length(thresholds)}.
}
\value{
A numeric matrix of expected ratings.
}
\author{
Chris Bradley (cbradley05@gmail.com)
}
\note{
Expected ratings are literally the expected value of the ordinal rating categories when treated as integers.  Expected ratings that cannot be calculated return as NA (e.g., if either the person or item measure is NA).  Intended use is for chi-squared tests or calculating infit and outfit statistics.
}
\seealso{
\code{\link{misfit}}
}
\examples{
# Using randomly generated values with minimum rating set to zero
im <- runif(20, -0.5, 0.5)
pm <- runif(50, -0.5, 0.5)
th <- c(-1.5, -0.8, 0.1, 0.6, 1.5)
m <- expdata(items = im, persons = pm, thresholds = th, minRating = 0)
}
