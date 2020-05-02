\name{thresh}
\alias{thresh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rating Category Thresholds
}
\description{
Estimates rating category thresholds for \code{msd} given rating scale data, item measures and person measures.
}
\usage{
thresh(data, items, persons)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a numeric matrix of ordinal rating scale data whose entries are integers with missing data set to NA.  Rows are persons and columns are items.  The ordinal rating scale is assumed to go from the smallest integer to the largest integer in \code{data} in integer steps.}
  \item{items}{a numeric vector of item measures with missing values set to NA (see Details).}
  \item{persons}{a numeric vector of person measures with missing values set to NA (see Details).}
}
\details{
The length of \code{items} must equal the number of columns in \code{data} and the length of \code{persons} must equal the number of rows in \code{data}. Neither \code{items} nor \code{persons} can consist of only NA.
}
\value{
A list whose elements are:
  \item{thresholds}{a vector of average rating category thresholds used by the persons when rating the items}
  \item{threshold_std_errors}{a vector of standard errors for the thresholds}
}
\author{
Chris Bradley (cbradley05@gmail.com)
}
\note{
\code{thresh} is a special case of \code{msd} when item measures and person measures are known.
}
\seealso{
\code{\link{msd}}
}
\examples{
# Using randomly generated values
d <- as.numeric(sample(0:5, 1000, replace = TRUE))
m <- matrix(d, nrow = 50, ncol = 20)
im <- runif(20, -2, 2)
pm <- runif(50, -2, 2)
th1 <- thresh(m, items = im, persons = pm)

# Anchor first 10 item measures and first 10 person measures
im[11:length(im)] <- NA
pm[11:length(pm)] <- NA
th2 <- thresh(m, items = im, persons = pm)
}
