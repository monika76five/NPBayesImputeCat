\name{Lcm}
\Rdversion{1.1}
\docType{class}
\alias{Lcm}

\title{Class \code{"Rcpp_Lcm"}}
\description{
This class implements the MCMC sampler for non-parametric imputation of discrete multivariate data described in Manrique-Vallier and Reiter (2014). It provides methods for updating and monitoring the sampler. 
}
\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}

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
}

\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{CurrentIteration}:}{the total number of iterations that have been run so far.}
    \item{\code{EnableTracer}:}{to check tracer status or to enable/disable the tracer.}
    \item{\code{MCZ}:}{the disjointed structural zero matrix.}
    \item{\code{snapshot}:}{retrieve a list with the current state of all the parameters in the sampler, including the imputed sample. A call the the "snapshot" method returns a list with the following components:\describe{
            \item{alpha:}{the concentration parameter of the stick breaking prior.}
            \item{k_star:}{the effective number number of latent classes (mixture components)}
            \item{Nmis: }{the size of the augmented sample.}
            \item{nu: }{a vector with the mixture weights}
            \item{z: }{a matrix with the current latent class assignment of each member of the sample}
            \item{ImputedX: }{the current raw imputed dataset. Use \link{GetDataFrame} to convert the raw data to a data frame of factors as defined in the input data set.}
            \item{psi:}{The conditional multinomial probabilties. A Lmax * K * J array, where Lmax is the maximum number of levels of all discrete factors in the dataset, J is the number of factors in the dataset, and K is the number of latent classes. Since variables might have different numbers of levels, unused entries in the first dimension are filled with NAs to complete Lmax.}
            }
  }
    \item{\code{traceable}:}{list of model parameters that can be traced by the tracer.}
    \item{\code{traced}:}{list of model parameters that are traced.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{SetTrace(paralist,num_of_iterations)}:}{set parameters to be traced.\describe{
            \item{paralist:}{a list of parameters to be traced.}
            \item{num_of_iterations:}{the maximum number of traced iterations.}
            }
    }
    \item{\code{Run(burnin, iter, thinning,silent)}:}{run MCMC iterations.
      \describe{
            \item{burnin:}{number of burn in iterations.}
            \item{iter:}{number of MCMC iterations.}
            \item{thinning:}{thinning parameter.}
            \item{silent:}{boolean indication if more iteration should be printed.}
            }
    }
    \item{\code{Resume()}:}{resume from an interrupted call to run method.}
    \item{\code{Parameters(paralist)}:}{retrieve a selected list of model parameters from last MCMC iteration.\describe{
            \item{paralist:}{a list of parameters to be traced.}
            }
    }
    \item{\code{GetTrace()}:}{retrieve all traced iterations. Returns a list with all the parameters set using the method SetTrace(). See description of snapshot{reference} method for a description of the parameters.  }
  }
}
\details{
Rcpp_lcm objects should be created with \link{CreateModel}. Please see the examples in the demo folder for more detailed explanation on model fitting and parameter tracing.
}
\references{
 Manrique-Vallier, D. and Reiter, J.P. (2013), "Bayesian Estimation of Discrete Multivariate Latent Structure Models with Structural Zeros", JCGS. 
  \cr\cr
  Si, Y. and Reiter, J.P. (2013), "Nonparametric Bayesian multiple imputation for incomplete categorical variables in large-scale assessment surveys", Journal of Educational and Behavioral Statistics, 38, 499 - 521 
  \cr\cr
  Manrique-Vallier, D. and Reiter, J.P. (2014), "Bayesian Multiple Imputation for Large-Scale Categorical Data with Structural Zeros", Survey Methodology.
}
