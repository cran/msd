\name{msdprob}
\alias{msdprob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rating Category Probabilities
}
\description{
Estimates the probability of observing each rating category.
}
\usage{
msdprob(x, thresholds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a real number representing a person minus item measure.}
  \item{thresholds}{a numeric vector of ordered rating category thresholds with no NA.}
}
\details{
It is assumed that \code{thresholds} partitions the real line into \code{length(thresholds)+1} ordered rating categories.
}
\value{
A vector of probabilities ordered from lowest to highest rating category.
}
\author{
Chris Bradley (cbradley05@gmail.com)
}
\examples{
p <- msdprob(1.4, thresholds = c(-1.1, -0.3, 0.5, 1.7, 2.2))
}
