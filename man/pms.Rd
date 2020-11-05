\name{pms}
\alias{pms}
\title{
Person Measures
}
\description{
Estimates person measures assuming item measures are known and all persons use the same set of rating category thresholds.
}
\usage{
pms(data, items, thresholds, misfit = FALSE, minRating = NULL)
}
\arguments{
  \item{data}{a numeric matrix of ordinal rating scale data whose entries are integers with missing data set to NA.  Rows are persons and columns are items.  The ordinal rating scale is assumed to go from the smallest to largest integer in integer steps unless \code{minRating} is specified (see Details).}
  \item{items}{a numeric vector of item measures with missing values set to NA.  The length of \code{items} must equal the number of columns in \code{data}.}
  \item{thresholds}{a numeric vector of ordered rating category thresholds with no NA.}
  \item{misfit}{logical for calculating infit and outfit statistics. Default is FALSE.}
  \item{minRating}{integer representing the smallest ordinal rating category.
  Default is NULL (see Details).}
}
\details{
\code{minRating} must be specified if either the smallest or largest possible rating category is not in \code{data} (i.e., no person used one of the extreme rating categories).  If \code{minRating} is specified, the ordinal rating scale is assumed to go from \code{minRating} to \code{minRating + length(thresholds)} in integer steps.
}
\value{
A list whose elements are:
  \item{person_measures}{a vector of person measures for each person}
  \item{person_std_errors}{a vector of standard errors for the persons}
  \item{infit_persons}{if \code{misfit} = TRUE, a vector of infit statistics for the persons}
  \item{outfit_persons}{if \code{misfit} = TRUE, a vector of outfit statistics for the persons}
}
\author{
Chris Bradley (cbradley05@gmail.com)
}
\note{
Person measures estimated with \code{pms} differ from those estimated with \code{msd} because \code{pms} assumes all persons use the same rating category thresholds while \code{msd} does not.  Intended use of \code{pms} is with an anchored set of items and thresholds.  Person measures that cannot be estimated will return as NA (e.g., if a person responds to all items with only the highest rating category, or with only the lowest rating category, that person's person measure cannot be estimated).
}
\seealso{
\code{\link{msd}}
}
\examples{
# Simple example with randomly generated values and lowest rating category = 0
d <- as.numeric(sample(0:4, 1000, replace = TRUE))
dm <- matrix(d, nrow = 50, ncol = 20)
im <- runif(20, -2, 2)
th <- sort(runif(4, -2, 2))
pm <- pms(data = dm, items = im, thresholds = th, misfit = TRUE, minRating = 0)
}
