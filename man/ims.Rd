\name{ims}
\alias{ims}
\title{
Item Measures
}
\description{
Estimates item measures assuming person measures are known and all persons use the same set of rating category thresholds.
}
\usage{
ims(data, persons, thresholds, misfit = FALSE, minRating = NULL)
}
\arguments{
  \item{data}{a numeric matrix of ordinal rating scale data whose entries are integers with missing data set to NA.  Rows are persons and columns are items.  The ordinal rating scale is assumed to go from the smallest to largest integer in integer steps unless \code{minRating} is specified (see Details).}
  \item{persons}{a numeric vector of person measures with missing values set to NA.  The length of \code{persons} must equal the number of rows in \code{data}.}
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
  \item{item_measures}{a vector of person measures for each person}
  \item{item_std_errors}{a vector of standard errors for the persons}
  \item{infit_items}{if \code{misfit} = TRUE, a vector of infit statistics for the items}
  \item{outfit_items}{if \code{misfit} = TRUE, a vector of outfit statistics for the items}
}
\author{
Chris Bradley (cbradley05@gmail.com)
}
\note{
Item measures estimated with \code{ims} differ from those estimated with \code{msd} because \code{ims} assumes all persons use the same rating category thresholds while \code{msd} does not.  Intended use of \code{ims} is with an anchored set of persons and thresholds.  Item measures that cannot be estimated will return as NA (e.g., if all responses to an item consist of only the highest rating category, or of only the lowest rating category, that item's item measure cannot be estimated).
}
\seealso{
\code{\link{msd}}
}
\examples{
# Simple example with randomly generated values and lowest rating category = 0.
d <- as.numeric(sample(0:4, 500, replace = TRUE))
dm <- matrix(d, nrow = 50, ncol = 10)
pm <- runif(50, -2, 2)
th <- sort(runif(4, -2, 2))
im <- ims(data = dm, persons = pm, thresholds = th, misfit = TRUE, minRating = 0)
}
