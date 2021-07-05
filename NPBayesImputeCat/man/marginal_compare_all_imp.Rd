\name{marginal_compare_all_imp}
\alias{marginal_compare_all_imp}
\title{
Plot estimated marginal probabilities from observed data vs imputed datasets
}
\description{
Plot estimated marginal probabilities from observed data vs imputed datasets
}
\usage{
marginal_compare_all_imp(obsdata, impdata, vars)
}
\arguments{
  \item{obsdata}{he observed data}
  \item{impdata}{the list of m imputed datasets}
  \item{vars}{the variable of interest}
}

\value{
  \item{Plot }{the barplot}
  \item{Comparison}{a table of marginal probabilies from observed data vs imputed datasets}
}
