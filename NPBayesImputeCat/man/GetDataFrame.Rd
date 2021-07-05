\name{GetDataFrame}
\alias{GetDataFrame}
\title{
Convert imputed data to a dataframe, using the same setting from original input data.
}
\description{
This is a utility function to convert the imputed data matrix to a dataframe. This function will be implemented as a RCPP internal function later on.
}
\usage{
GetDataFrame(dest, from, cols = 1:NCOL(from))
}
\arguments{
  \item{dest}{the imputed output data matrix. }
  \item{from}{the original input dataframe. }
  \item{cols}{optinal. Always use default for now.}
}

\value{
The returned dataframe object for imputed data. 
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