\name{misfit}
\alias{misfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Infit and Outfit Statistics
}
\description{
Calculates infit and outfit for items and persons.
}
\usage{
misfit(data, items, persons, thresholds, minRating = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a numeric matrix of ordinal rating scale data whose entries are integers with missing data set to NA.  Rows are persons and columns are items.  The ordinal rating scale is assumed to go from the smallest to largest integer in integer steps unless \code{minRating} is specified.}
  \item{items}{a numeric vector of item measures with missing values set to NA.}
  \item{persons}{a numeric vector of person measures with missing values set to NA.}
  \item{thresholds}{a numeric vector of ordered rating category thresholds with no NA.}
  \item{minRating}{integer representing the smallest ordinal rating category.
  Default is NULL (see Details).}
}
\details{
\code{minRating} must be specified if either the smallest or largest possible rating category is not in \code{data} (no person used one of the extreme rating categories).  If \code{minRating} is specified, the ordinal rating scale is assumed to go from \code{minRating} to \code{minRating + length(thresholds)}.
}
\value{
A list whose elements are:
  \item{infit_items}{a vector of infit statistics for the items}
  \item{outfit_items}{a vector of outfit statistics for the items}
  \item{infit_persons}{a vector of infit statistics for the persons}
  \item{outfit_persons}{a vector of outfit statistics for the persons}
}
\author{
Chris Bradley (cbradley05@gmail.com)
}
\examples{
# Using randomly generated values
d <- as.numeric(sample(0:5, 1000, replace = TRUE))
dm <- matrix(d, nrow = 50, ncol = 20)
im <- runif(20, -2, 2)
pm <- runif(50, -2, 2)
th <- sort(runif(5, -2, 2))
m <- misfit(data = dm, items = im, persons = pm, thresholds = th)

# If the lowest or highest rating category is not in \code{data}, specify \code{minRating}
dm[dm == 0] <- NA
m2 <- misfit(data = dm, items = im, persons = pm, thresholds = th, minRating = 0)
}
