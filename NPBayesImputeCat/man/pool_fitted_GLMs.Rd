\name{pool_fitted_GLMs}
\alias{pool_fitted_GLMs}
\title{
Pool estimates of fitted GLM models in imputed or synthetic datasets
}
\description{
Pool estimates of fitted GLM models in imputed or synthetic datasets
}
\usage{
pool_fitted_GLMs(GLMResults, method = 
                      c("imputation", "synthesis_full", "synthesis_partial"))
}
\arguments{
  \item{GLMResults}{output from the fit_GLMs function}
  \item{method}{choose between "imputation", "synthesis_full", "synthesis_partial"}
}
\value{Results: a list of GLM results after combining rules}
