\name{UpdateX}
\alias{UpdateX}

\title{
  Allow user to update the model with data matrix of same kind.
}
\description{
Allow user to replace initial matrix with a new data matrix of same size and same number of factors. This is not intended for general use and is only useful for very specific circumstance. 
}
\usage{
UpdateX(model, X)
}

\arguments{
  \item{model}{The Rcpp model object created by the CreateModel function.}
  \item{X}{a data frame with the dataset with missing values. All variables must be unordered factors.}
}
