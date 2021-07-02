\name{GetMCZ}
\alias{GetMCZ}
\title{
Convert disjointed structrual zeros to a dataframe, using the same setting from original structrual zero data. 
}
\description{
This is a utility function to convert the disjointed structrual zero matrix to a dataframe. This function will be implemented as a RCPP internal function later on.
}
\usage{
GetMCZ(dest, from, mcz, cols = 1:NCOL(from))
}
\arguments{
  \item{dest}{
the output data matrix for disjointed structrual zeros.
}
  \item{from}{the original input dataframe.}
  \item{mcz}{the original input dataframe for structrual zeros.}
  \item{cols}{optinal. Always use default for now.}
}

\value{
The returned dataframe object for disjointed structrual zeros.}
\references{
 Manrique-Vallier, D. and Reiter, J.P. (2013), "Bayesian Estimation of Discrete Multivariate Latent Structure Models with Structural Zeros", JCGS. 
  \cr\cr
  Si, Y. and Reiter, J.P. (2013), "Nonparametric Bayesian multiple imputation for incomplete categorical variables in large-scale assessment surveys", Journal of Educational and Behavioral Statistics, 38, 499 - 521 
  \cr\cr
  Manrique-Vallier, D. and Reiter, J.P. (2014), "Bayesian Multiple Imputation for Large-Scale Categorical Data with Structural Zeros", Survey Methodology.
}
