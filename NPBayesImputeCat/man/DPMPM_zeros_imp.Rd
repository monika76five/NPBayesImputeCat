\name{DPMPM_zeros_imp}
\alias{DPMPM_zeros_imp}
\title{
Use DPMPM models to impute missing data where there are no structural zeros
}
\description{
Use DPMPM models to impute missing data where there are no structural zeros
}
\usage{
DPMPM_zeros_imp(X, MCZ, Nmax, nrun, burn, thin, K, aalpha, balpha, m, seed, silent)
}
\arguments{
  \item{X}{data frame for the data containing missing values}
  \item{MCZ}{data frame containing the structural zeros definition}
  \item{Nmax}{an upper truncation limit for the augmented sample size}
  \item{nrun}{number of mcmc iterations}
  \item{burn}{number of burn-in iterations}
  \item{thin}{thining parameter for outputing iterations}
  \item{K}{number of latent classes}
  \item{aalpha}{the hyperparameters in stick-breaking prior distribution for alpha}
  \item{balpha}{the hyperparameters in stick-breaking prior distribution for alpha}
  \item{m}{number of imputations}
  \item{seed}{choice of random seed}
  \item{silent}{Default to TRUE. Set this parameter to FALSE if more iteration info are to be printed}
}
\value{
  \item{impdata }{m imputed datasets}
  \item{origdata }{original data containing missing values}
  \item{alpha }{save posterior draws of alpha, which can be used to check MCMC convergence}
  \item{kstar }{saved number of occupied mixture components, which can be used to track whether K is large enough}
  \item{Nmax }{saved posterior draws of the augmented sample size, which can be used to check MCMC convergence}

}
