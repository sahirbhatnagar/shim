% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{coef.cv.shim}
\alias{coef.cv.shim}
\title{Make predictions from a "cv.shim" object}
\usage{
\method{coef}{cv.shim}(object, s = c("lambda.1se", "lambda.min"), ...)
}
\arguments{
\item{object}{object of class cv.shim from cv.shim function}

\item{s}{Value(s) of the penalty parameter lambda at which predictions are
required. Default is the value \code{s="lambda.1se"} stored on the cv.shim
object. Alternatively \code{s="lambda.min"} can be used.}
}
\description{
Make predictions from a "cv.shim" object
}

