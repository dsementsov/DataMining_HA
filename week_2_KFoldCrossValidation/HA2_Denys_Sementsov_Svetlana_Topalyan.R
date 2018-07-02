###########################################################
#                      Data Mining                        #
#                   Home Assignment 2                     #
#           Denys Sementsov, Svetlana Topalyan            #
#        stu110401, Mn1001084     Mn1110072               #
###########################################################



# ===========================================1===========================================
set.seed(90) 
genSample = function(n)                             
# function that will take n as input and output x/y data frame of the form 
# y = 1 + x^3 + u, where x,u ~N(0,1)
{
  # generating random variables
  x = rnorm(n, mean = 0, sd = 1)    
  u = rnorm(n, mean = 0, sd = 1)
  
  # populating y vector with 0s to match length of x
  random_y = seq(0,200-1)
  
  # changing y to correspond the function
  for (i in 0:n)                      
  {
    y[i] = 1 + x[i]^3 + u[i]
  }
  
  # sewing x and y together and returning
  result = data.frame(x,y)

  return (result)            
}

# generating data with 200 observations     
data = genSample(200)     

# train and test samples split 80/20
trainSample = sample(1:nrow(data), 0.8*nrow(data))           
testSample = -trainSample

# train data is 80% of our data
trainData = data[trainSample,]                              

# plot train data
par(mfrow = c(2,2))                 
plot(trainData)                                              

# plot the true regression over training data
lm.fit = lm (y ~ 1 + I(x^3), data = data)
lines(sort(data$x), predict(lm.fit, data = data)[order(data$x)], 
      col = "red",type = 'l')

# ===========================================2===========================================


#' Function that calculates MSE given data
#' 
#'  @param input1 actual data
#'  @param imput2 predicted test data, generated over trained data
#'  @return mse value
#'  
#'  @seealso data should be atached, passing data$x will generate wrong results
calc_mse = function(actual, predicted)
{
  mse = mean((actual - predicted)^(2))
  return (mse)
}


#' @description Calculate k-fold cross-validation mse
#' @param input1 formula for the model you fit
#' @param input2 data model should be fitted for
#' @param #K of folds
#' 
#' @details Will produce a warning if data cannot be split evenly in K folds
my.kfold = function(model, mData, K)
{
  CV = NA
  if (nrow(mData) %% K != 0) warning("Data cannot be splitted into ", K," equal-sized parts")   # warning message
  
  # dividing data into k folds
  # gives us massive of repetitive 1 to K indexes we can use to split our data
  buckets = rep(1:K, nrow(mData)/K)     
  indexes = seq(1,200)
  
  # spliting data using buckets as indexes gives us K folds without replacement
  folds = split(indexes, buckets) 
  
  # initializing mse
  mse = c(0,K)                         
  
  # k-fold loop
  for (i in 1:K)
  {
    # Training a model using formula and K-1 training set
    fit = lm(model, data = mData[-folds[[i]],])
    
    # Predicting data 
    fit.prediction = predict(fit, newdata = mData[folds[[i]],])
    
    #calculating mse
    mse[i] = calc_mse(mData$y[folds[[i]]], fit.prediction)
    
  }  # predictions in data should always be named y for function to work
  CV = 1/K * sum(mse)
  return (CV)
}



# ===========================================3===========================================

# creating table for simulation, nrow will be polynomial order
# simTrainMean = matrix(1:1000, ncol = 100)
# simTestMean = matrix(1:1000, ncol = 100)
# simKFold = matrix(1:1000, ncol = 100)

# preparing data frame for simulation data
simulation = data.frame("Polynomial Degree" = rep(1:10,100), "Train MSE" = rep(0,1000), "Test MSE" = rep(0,1000), "KFold CV" = rep(0,1000))

# initializing temp vectors for choosing right polynomials
tmsepoly = rep(0,100)
vmsepoly = rep(0,100)
kfoldcvpoly = rep(0,100)

# number of replication of simulation study
for (r in 1:100)      
{
  # 1. Generate sample
  simData = genSample(200)

  # 2. Divide into train and test samples (subsets) 80/20
  simTrainSample = sample(1:nrow(simData), 0.8*nrow(simData))      
  simTestSample = -simTrainSample
  
  # 3. Fit polinomials 1-10 for each simulation
  # p for the degree of polynomial
  for (p in 1:10)     
  {
    # fitting polynomial
    fit = lm(y ~ poly(x, degree = p), data = simData, subset = simTrainSample)
    
    # 4. Generating results for the study
    # 4.1 Column1 Degree of the polynomial
    simulation[(r-1)*10 + p,1] = p
    
    # 4.2 Column2 Mse calculated for train data on train data (to show that it is a bad practice)
    simulation[(r-1)*10 + p,2] = calc_mse(actual = simData$y[simTrainSample], 
                                        predicted  = predict(fit, newdata = simData[simTrainSample,]))
    
    # 4.3 Column3 Mse calculated for test data
    simulation[(r-1)*10 + p,3] = calc_mse(actual = simData$y[simTestSample], 
                                          predicted = predict(fit, newdata = simData[simTestSample,]))
    
    # 4.4 Column4 K-fold Cross-validation mse
    simulation[(r-1)*10 + p,4] = my.kfold(y ~ poly(x, degree = p), simData, K = 10)
  }
  
  # 5. Choosing the polynomial according to different mses
  indexStart = (r-1)*10 + 1
  indexEnd = indexStart + 9
  tmsepoly[r] = simulation$Polynomial.Degree[indexStart - 1 + which.min(simulation$Train.MSE[indexStart:indexEnd])]  # choosing right polynomial
  vmsepoly[r] = simulation$Polynomial.Degree[indexStart - 1 + which.min(simulation$Test.MSE[indexStart:indexEnd])]
  kfoldcvpoly[r] = simulation$Polynomial.Degree[indexStart -1 + which.min(simulation$KFold.CV[indexStart:indexEnd])]
}

# Output table generation
outputTable = data.frame("Polynomial degree" = seq(1,10), "Mean, Train" = rep(0,10), "SD, Train" = rep(0,10),
                         "Mean, Val" = rep(0,10), "SD, Val" = rep(0,10),
                         "Mean, CV" = rep(0,10), "SD, CV" = rep(0,10))


# Feeding table with data
for (i in 1:10)
{
  outputTable$Mean..Train[i] = mean(simulation$Train.MSE[which(simulation$Polynomial.Degree == i)])
  outputTable$SD..Train[i] = sd(simulation$Train.MSE[which(simulation$Polynomial.Degree == i)])
  outputTable$Mean..Val[i] = mean(simulation$Test.MSE[which(simulation$Polynomial.Degree == i)])
  outputTable$SD..Val[i] = sd(simulation$Test.MSE[which(simulation$Polynomial.Degree == i)])
  outputTable$Mean..CV[i] = mean(simulation$KFold.CV[which(simulation$Polynomial.Degree == i)])
  outputTable$SD..CV[i] = sd(simulation$KFold.CV[which(simulation$Polynomial.Degree == i)])
}

options(digits = 2, scipen = 999)       # exactly 2 digits after comma, no scientific notation
outputTable



# ===========================================4===========================================
# how many times which polynomial was choosen

hist(tmsepoly, xlab = "Ploy Degree", ylab = "Times Chosen",  main = "Training error", breaks = c(7,8,9,10,11))          # training error
hist(vmsepoly, xlab = "Ploy Degree", ylab = "Times Chosen",  main = "Validation Set")          # validation set
hist(kfoldcvpoly, xlab = "Ploy Degree", ylab = "Times Chosen",  main = "10-Fold CV")           # 10-Fold CV

#nicier plots
#install.packages("ggplot2")
#library(ggplot2)
#qplot(tmsepoly)
#qplot(tmsepoly, 
#      geom = "histogram",
#      main = "Training Error",
#      xlab = "Ploy Degree", ylab = "Times Chosen")
#qplot(vmsepoly, 
#      geom = "histogram",
#      main = "Validation Set",
#      xlab = "Ploy Degree", ylab = "Times Chosen")
#qplot(kfoldcvpoly, 
#      geom = "histogram",
#      main = "10-Fold CV",
#      binwidth = 0.5,
#     xlab = "Ploy Degree", ylab = "Times Chosen")

outputTable

# ==========================================5===========================================

#  First and foremost you should never use train data for MSE and model validation. It will most likely never

#  show correct or meaningful validation for your model. As we have seen from the output table, train mse tend to 

#  lower even after the model doesn't show good results on test subset (Mean mse column suggest us to use 9-10 order

#  polynomials, however testing the model on the test data show significant decline in performance).

#  According to validation set and k-fold cross-validation methods for evaluating this simulated data we should use

#  3rd or 4th order polynomials with 3rd showing slight better results from test validation method.

#  Those polynomials they perform better with lower MSE and standard deviations. Even though the frequency of chosing 

#  higher degrees polynomials according to the 10-Fold CV is quite high, we observe enormous spike in mean and sd on average.

#  That underlines the importance of human factor in data analysis.
 

