\name{kstar_MCMCdiag}
\alias{kstar_MCMCdiag}
\title{
Perform MCMC diagnostics for kstar
}
\description{
A helper function to perform MCMC diagnostics for kstar
}
\usage{
kstar_MCMCdiag(kstar, nrun, burn, thin)
}
\arguments{
  \item{kstar}{the vector output of kstar from running the DPMPM model}
  \item{nrun}{number of mcmc iterations used in running the DPMPM model}
  \item{burn}{number of burn-in iterations used in running the DPMPM model}
  \item{thin}{number of thinning used in running the DPMPM mode}
}

\value{

 \item{Traceplot }{the traceplot of kstar post burn-in and thinning}
 \item{Autocorrplot }{the autocorrelation plot of kstar post burn-in and thinning}

}
