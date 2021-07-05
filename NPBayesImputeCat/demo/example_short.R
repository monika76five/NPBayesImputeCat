#this R script demos the basic usage of model fitting and data #imputaiton

#load package and the example dataset
require(NPBayesImputeCat)
data('NYexample')

#create the model in one step
#CreateModel(X,MCZ,K,Nmiss_Max, alphaa,aplhab,seed)
model <- CreateModel(X,MCZ,50,200000,0.25,0.25,8888)
#to run a model without without structural zeros, set MCZ to NULL, and Nmax to 0
#model <- CreateModel(X,NULL,100,0,0.25,0.25,8888)


#run 100 burns, 1000 mcmc iterations and thinning every 10 #iterations
model$Run(100,1000,100)

#retrieve parameters from the final iteration
result <- model$snapshot

#convert ImputedX matrix to dataframe, using proper factors/names etc.
IX <- GetDataFrame(result$ImputedX,X)
View(IX)
