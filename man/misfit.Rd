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
misfit(data, items, persons, thresholds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a numeric matrix of ordinal rating scale data whose entries are integers with missing data set to NA.  Rows are persons and columns are items.  The ordinal rating scale is assumed to go from the smallest to largest integer in integer steps.}
  \item{items}{a numeric vector of item measures with missing values set to NA.}
  \item{persons}{a numeric vector of person measures with missing values set to NA.}
  \item{thresholds}{a numeric vector of ordered rating category thresholds with no NA.}
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
m <- matrix(d, nrow = 50, ncol = 20)
im <- runif(20, -0.5, 0.5)
pm <- runif(50, -0.5, 0.5)
th <- c(-1.5, -0.8, 0.1, 0.6, 1.5)
m <- misfit(data = m, items = im, persons = pm, thresholds = th)
}
