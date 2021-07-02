\name{compute_probs}
\alias{compute_probs}
\title{
Estimating marginal and joint probabilities in imputed or synthetic datasets
}
\description{
Estimating marginal and joint probabilities in imputed or synthetic datasets
}
\usage{
compute_probs(InputData, varlist)
}
\arguments{
  \item{InputData}{a list of imputed or synthetic datasets}
  \item{varlist}{a list of variable names (or combination of names) to evaluate (marginal or joint) probabilities for}
}
\value{
  Results: a list of marginal and joint probability results after combining rules
}
