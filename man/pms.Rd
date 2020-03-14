\name{pms}
\alias{pms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Person Measures
}
\description{
Estimates person measures assuming all persons use the same rating category thresholds.
}
\usage{
pms(data, items, thresholds, misfit = FALSE)
}
\arguments{
  \item{data}{a numeric matrix of ordinal rating scale data whose entries are integers with missing data set to NA.  Rows are persons and columns are items.  The ordinal rating scale is assumed to go from the smallest to largest integer in integer steps.}
  \item{items}{a numeric vector of item measures with missing values set to NA.  The length of \code{items} must equal the number of columns in \code{data}.}
  \item{thresholds}{a numeric vector of ordered rating category thresholds with no NA.}
  \item{misfit}{logical for calculating infit and outfit statistics. Default is FALSE.}
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
Person measures estimated with \code{pms} differ from those estimated with \code{msd} because \code{pms} assumes all persons use the same rating category thresholds while \code{msd} does not.  Intended use of \code{pms} is with an anchored set of items and thresholds. Person measures that cannot be estimated will return as NA.
}
\seealso{
\code{\link{msd}}
}
\examples{
# Using randomly generated values
d <- as.numeric(sample(0:5, 1000, replace = TRUE))
m <- matrix(d, nrow = 50, ncol = 20)
im <- runif(20, -0.5, 0.5)
th <- c(-1.5, -0.8, 0.1, 0.6, 1.5)
p <- pms(data = m, items = im, thresholds = th, misfit = TRUE)
}
