\name{NPBayesImputeCat-package}
\alias{NPBayesImputeCat-package}
\alias{NPBayesImputeCat}
\docType{package}
\title{
Bayesian Multiple Imputation for Large-Scale Categorical Data with Structural Zeros
}
\description{
This package implements a fully Bayesian, joint modeling approach to multiple imputation
for categorical data based on latent class models with structural zeros. The idea is to model the implied
contingency table of the categorical variables as a mixture of independent multinomial distributions,
estimating the mixture distributions nonparametrically with Dirichlet process prior
distributions. Mixtures of multinomials can describe arbitrarily complex dependencies and
are computationally expedient, so that they are effective general purpose multiple imputation
engines. In contrast to other approaches based on loglinear models or chained equations,
the mixture models avoid the need to specify (potentially many) models, which can be a
very time-consuming task with no guarantee of a theoretically coherent set of models. The
package is designed to include for structural zeros, i.e., certain combinations of variables are not possible
a priori.
}
\details{
\tabular{ll}{
Package: \tab NPBayesImputeCat\cr
Type: \tab Package\cr
Version: \tab 0.4\cr
Date: \tab 2021-06-30\cr
License: \tab GPL(>=3) \cr
}
}

\author{
Quanli Wang, Daniel Manrique-Vallier, Jerome P. Reiter and Jingchen Hu

Maintainer: Quanli Wang<quanli@stat.duke.edu>
}
\references{
  Manrique-Vallier, D. and Reiter, J.P. (2013), "Bayesian Estimation of Discrete Multivariate Latent Structure Models with Structural Zeros", JCGS. 
  \cr\cr
  Si, Y. and Reiter, J.P. (2013), "Nonparametric Bayesian multiple imputation for incomplete categorical variables in large-scale assessment surveys", Journal of Educational and Behavioral Statistics, 38, 499 - 521 
  \cr\cr
  Manrique-Vallier, D. and Reiter, J.P. (2014), "Bayesian Multiple Imputation for Large-Scale Categorical Data with Structural Zeros", Survey Methodology.
}
\keyword{ package }
\examples{
require(NPBayesImputeCat)
#Please use NYexample data set for a more realistic example
data('NYMockexample')

#create the model
model <- CreateModel(X,MCZ,10,10000,0.25,0.25,8888)

#run 1 burnins, 2 mcmc iterations and thin every 2 iterations
model$Run(1,2,2,TRUE)

#retrieve parameters from the final iteration
result <- model$snapshot

#convert ImputedX matrix to dataframe, using proper factors/names etc.
ImputedX <- GetDataFrame(result$ImputedX,X)
#View(ImputedX)

#Most exhauststic examples can be found in the demo below
#demo(example_short)
#demo(example)

}
