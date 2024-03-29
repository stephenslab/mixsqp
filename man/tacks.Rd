% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tacks.R
\docType{data}
\name{tacks}
\alias{tacks}
\title{Beckett & Diaconis tack rolling example.}
\format{
\code{tacks} is a list with the following elements:

\describe{
  \item{L}{9 x 299 likelihood matrix.}

  \item{w}{Numeric vector of length 9 specifying the weights
    associated with the rows of \code{L}.}

  \item{x}{Solution provided by the \code{KWDual}
    solver.}
}
}
\description{
This data set contains the likelihood matrix and
weights for the Beckett-Diaconis tacks example, in which the data
are modeled using a binomial mixture. These data were generated by
running the "Bmix1" demo from the REBayes package, and saving the
arguments passed to \code{KWDual}, as well as the (normalized)
solution returned by the \code{KWDual} call.
}
\examples{

# The optimal solution for the tack example is extremely sparse.
data(tacks)
plot(tacks$x,type = "l",col = "royalblue")

}
\keyword{data}
