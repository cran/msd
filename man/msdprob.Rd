\name{msdprob}
\alias{msdprob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rating Category Probabilities
}
\description{
Estimates the probability of observing each rating category given a set of ordered rating category thresholds.
}
\usage{
msdprob(x, thresholds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a real number or a vector of real numbers with no NA representing a set of person minus item measures.}
  \item{thresholds}{a numeric vector of ordered rating category thresholds with no NA.}
}
\details{
It is assumed that \code{thresholds} partitions the real line into \code{length(thresholds)+1} ordered intervals that represent the rating categories.
}
\value{
A matrix of probabilities where each of the \code{length(thresholds)+1} rows represents a different rating category (lowest rating category is the top row) and each of the \code{length(x)} columns represents a different person minus item measure.
}
\author{
Chris Bradley (cbradley05@gmail.com)
}
\note{
\code{msdprob} can be used to create probability curves, which represent
the probability of rating an item with each rating category as a function
of the person measure minus item measure (see Examples).
}
\examples{
# Simple example
p <- msdprob(c(1.4, -2.2), thresholds = c(-1.1, -0.3, 0.5, 1.7, 2.2))

# Plot probability curves — each curve represents the probability of
# rating an item with a given rating category as a function of the
# person measure minus item measure.
x <- seq(-6, 6, 0.1)
p <- msdprob(x, thresholds = c(-3.2, -1.4, 0.5, 1.7, 3.5))
plot(0, 0, xlim = c(-6, 6), ylim = c(0, 1), type = "n",
    xlab = "Person minus item measure", ylab = "Probability")
for (i in seq(1, dim(p)[1])){
  lines(x, p[i,], type = "l", lwd = "2" , col = rainbow(6)[i])
}
}
