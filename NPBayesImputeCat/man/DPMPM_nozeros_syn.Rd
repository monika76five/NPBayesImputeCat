\name{DPMPM_nozeros_syn}
\alias{DPMPM_nozeros_syn}
\title{
Use DPMPM models to synthesize data where there are no structural zeros
}
\description{Use DPMPM models to synthesize data where there are no structural zeros
}
\usage{
DPMPM_nozeros_syn(X, dj, nrun, burn, thin, K, aalpha, balpha, m, vars, seed, silent)
}
\arguments{
  \item{X}{data frame for the original data}
  \item{dj}{a vector recording the number of categories of the variables}
  \item{nrun}{number of mcmc iterations}
  \item{burn}{number of burn-in iterations}
  \item{thin}{thining parameter for outputing iterations}
  \item{K}{number of latent classes}
  \item{aalpha}{the hyperparameters in stick-breaking prior distribution for alpha}
  \item{balpha}{the hyperparameters in stick-breaking prior distribution for alpha}
  \item{m}{number of synthetic datasets}
  \item{vars}{the names of variables to be synthesized}
  \item{seed}{choice of random seed}
  \item{silent}{Default to TRUE. Set this parameter to FALSE if more iteration info are to be printed}
}

\value{
  \item{syndata }{m  synthetic datasets}
  \item{origdata }{original data}
  \item{alpha }{saved posterior draws of alpha, which can be used to check MCMC convergence}
  \item{kstar }{saved number of occupied mixture components, which can be used to track whether K is large enough}

}
