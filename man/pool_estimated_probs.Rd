\name{pool_estimated_probs}
\alias{pool_estimated_probs}
\title{
Pool probability estimates from imputed or synthetic datasets
}
\description{
Pool probability estimates from imputed or synthetic datasets
}
\usage{
pool_estimated_probs(ComputeProbsResults, method = 
                      c("imputation", "synthesis_full", "synthesis_partial"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ComputeProbsResults}{output from the compute_probs function}
  \item{method}{choose between "imputation", "synthesis_full", "synthesis_partial"}
}
\value{Results: a list of marginal and joint probability results after combining rules}
