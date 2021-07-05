\name{marginal_compare_all_syn}
\alias{marginal_compare_all_syn}
\title{
Plot estimated marginal probabilities from observed data vs synthetic datasets
}
\description{
Plot estimated marginal probabilities from observed data vs synthetic datasets
}
\usage{
marginal_compare_all_syn(obsdata, syndata, vars)
}
\arguments{
  \item{obsdata}{the observed data}
  \item{syndata}{the list of m imputed datasets}
  \item{vars}{the variable of interest}
}
\value{ \item{Plot }{the barplot}
  \item{Comparison}{a table of marginal probabilies from observed data vs imputed datasets}
}

