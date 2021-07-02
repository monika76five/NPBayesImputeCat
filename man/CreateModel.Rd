\name{CreateModel}
\alias{CreateModel}
\title{
Create and initialize the Lcm model object
}
\description{
CreateModel creates and initializes an Lcm \link{Lcm} object for non-parametric multiple imputation of discrete multivariate categorical data with or without structural zeros.
}
\usage{
CreateModel(X, MCZ, K, Nmax, aalpha, balpha,seed)
}
\arguments{
  \item{X}{a data frame with the dataset with missing values. All variables must be unordered factors.}
  \item{MCZ}{a dataframe with the definition of the structural zeros. Placeholder components are represented with NAs. Variables in MCZ must be factors with the same levels as X. Rows do not need to define disjoint regions of the contingency table. See Manrique-Vallier and Reiter (2014) for details of the definition of structural zeros.  MCZ should be set to NULL when there are no structure zeros.}
  \item{K}{the maximum number of mixture components.}
  \item{Nmax}{An upper truncation limit for the augmented sample size. This parameter will be ignored(set to 0) when there is no structural zeros.}
  \item{aalpha}{the hyper parameter 'a' for alpha in stick-breaking prior distribution.}
  \item{balpha}{the hyper parameter 'b' for alpha in stick-breaking prior distribution.}
  \item{seed}{the random seed for sampling. When setting to NULL(default), the random seed will be set randomly.}
}
\details{
This function should be the first function one should call to use the `NPBayesImputeCat`
library. The returned model is a `Rcpp lcm` object. See ? `Rcpp Lcm-class`
for more details on the fields available and their arguments. }

\value{

CreateModel returns an Rcpp_lcm object. The returned model object will be referenced in all subsequent calls.
}

\references{
  Manrique-Vallier, D. and Reiter, J.P. (2013), "Bayesian Estimation of Discrete Multivariate Latent Structure Models with Structural Zeros", JCGS. 
  \cr\cr
  Si, Y. and Reiter, J.P. (2013), "Nonparametric Bayesian multiple imputation for incomplete categorical variables in large-scale assessment surveys", Journal of Educational and Behavioral Statistics, 38, 499 - 521 
  \cr\cr
  Manrique-Vallier, D. and Reiter, J.P. (2014), "Bayesian Multiple Imputation for Large-Scale Categorical Data with Structural Zeros", Survey Methodology.
}
\examples{
require(NPBayesImputeCat)
#Please use NYexample data set for a more realistic example
data('NYMockexample')

#create the model
model <- CreateModel(X,MCZ,10,10000,0.25,0.25,8888)

#run 1 burnins, 2 mcmc iterations and thin every 2 iterations
model$Run(1,2,2,FALSE)

#retrieve parameters from the final iteration
result <- model$snapshot

#convert ImputedX matrix to dataframe, using proper factors/names etc.
ImputedX <- GetDataFrame(result$ImputedX,X)
#View(ImputedX)
}
