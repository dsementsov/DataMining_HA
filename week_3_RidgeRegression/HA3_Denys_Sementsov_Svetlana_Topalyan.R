###########################################################
#                      Data Mining                        #
#                   Home Assignment 2                     #
#           Denys Sementsov, Svetlana Topalyan            #
#        stu110401, Mn1001084     Mn1110072               #
##########################################################






#### 1. Generating samples of size n = 100 of all random variables ####
options(scipen=999) 
set.seed(333)

# Generating random variables
n = 100

u = rnorm(n = 100, mean = 0, sd = 1)
x1 = rnorm(n = 100, mean = 0, sd = 1)
x2 = rnorm(n = 100, mean = 0, sd = 1)
x3 = rnorm(n = 100, mean = 0, sd = 1)

# Our true model 
y = as.matrix(1 + 0.5*x1^2 + u)

# refactoring data
data = data.frame(y, x1, x2, x3, u)

# creating our model
formula = model.matrix(y ~ I(x1^2) + I(x2^2) + I(x3^2) + u, data = data)

# getting only x's observations from out formula
x = formula[,2:4]

#### 2. Estimating model usig ridge regression ####

# We can do ridge regression using two approaches.

#### Frist idea is to take my.OLS function from tutorial 1 and edit it  ####
# as it would take lambda as part of equasion for Beta
# B = (X'X + lambda*I)**-1 * X'y, solving the problem numerically
# For this tutorial we know that all variables are drawn from the same distribution with mean 0, have standard deviation 1
# and will not scale dependent nor independent variables


#' @description  Fitting ridge regression 
#' 
#' @param Input1: Y - Dependent variables, vector kx1
#' @param Input2: X - Regressors, matrix k x n (model.matrix)
#' @param Input3: Lambda value (Float, Default = 0)
#' @param Input4: Estimate or center Beta?
#' 
#' @return Coefficients of ridge regression
#' 
my_ridge_regression = function (x, y, lambda = 0, shrink_beta = TRUE)
  {
  
  x = as.matrix(x)
  y = as.matrix(y)
  
  # As always, verifying input for the function
  if (nrow(y) != nrow(x)) warning("Input error! Make sure that number of rows of Y = number of rows of X")

  # converting data into matricies for better multiplication functionality
  mat_y = y
  mat_x = cbind(rep(1, nrow(x)), x)
  
  
  # preparing identity matrix times lambda, we will have to add lambda to a product of t(X) * X - that will be n+1 x n+1 matrix for matrix with intercept
  lambda_matrix = lambda * diag(ncol(mat_x))
  if (!shrink_beta)
  {
    # To avoid shrinkage effect on intercept
    lambda_matrix[1] = 0
  }
  
  
  # estimating betas using LS formula with shrinkage coefficient
  estimated_betas = solve( (t(mat_x)  %*%  mat_x) + lambda_matrix, t(mat_x) %*% mat_y )

  # finding RSS + Shrinkage estimate
  error = mat_y - mat_x %*% estimated_betas
  residuals_squared = t(error) %*% error
  error_rate = residuals_squared + lambda * sum(estimated_betas * estimated_betas)
  
  return(list("Coefficients" = as.matrix(estimated_betas), "Error" = as.vector(error_rate))) 
}


#### Checking true values of ridge regression ####

library(glmnet)
lambdas = 10^seq(10, -2, -0.1)
ridge.test = glmnet(x, y, alpha = 0, lambda = lambdas, standardize = FALSE)
best.lambda.index = which.min(ridge.test$lambda)
ridge.test$lambda[best.lambda.index]
coef(ridge.test)[,best.lambda.index]

#################################################

my_ridge_regression(x = x, y = y, lambda = ridge.test$lambda[best.lambda.index])

# my ridge regression returns slightly different coefficients than that of glmnet function - this might be of internall fitting procedure of predict
# the difference in coefficients is due to construction

#### The second idea is to use optim() method to minimize function of lambda instead of computing it by hand ####

f = function(par, lambda, y, x)
{
  betas = c(par[2], par[3], par[4])
  intercept = par[1]
  sum( ( y  - intercept -  x %*% betas)^2 )  + lambda * sum ( betas^2 )
}

nlm(p = c(0,0,0,0), f, lambda = ridge.test$lambda[best.lambda.index], y = y, x = x)
optim(par = c(0,0,0,0), f, lambda = ridge.test$lambda[best.lambda.index], y = y, x = x)

# both optim and nlm give us coefficient close to numerical solution of my_ridge_regression than to glmnet and further shrinking them towards zero

#### 3. Picking right shrinkage intensivity ####

# Straightforward iteration over lambda values and performing K-fold cross-validation to choose lambda with lowest corresponding CV.MSE

# K-fold function
#' @description This function will perform cross-validation based on your data and lambda values
#' 
#' @param Input: (x,y, K = 10, lambda = 0, type and shrink parameters for)
#' 
#' @return Cross-validated MSE for given inputs
#' 
my_kfold_ridge = function(x, y, K = 10, lambda = 0, fn = NA, type = "optim", shrink = TRUE) 
{
  x = as.matrix(x)
  y = as.matrix(y)
  
  CV = NA
  # Dimension check
  if (nrow(x) %% K != 0) warning("Data cannot be splitted into ", K," equal-sized parts")
  
  # Creating K equal folds indexies
  folds = cut(seq(1, nrow(data)), breaks = K,labels = FALSE)
  
  # Initialization of MSE vector
  MSE = c()
  
  # main loop
  for(i in 1:K)
  {
    # separating current validation and train data
    train_index = !folds == i
    validation_index = folds == i
    
    ### fitting with optim
    if (type == "optim")
    {
      coef =  optim(par = c(0,0,0,0), fn, lambda = lambda, y = y[train_index,], x = x[train_index,])$par
      intercept = coef[1]
      coef = coef[-1]
    }
    else if (type == "ridge")
    {
    # fitting with my_ridge
      coef = my_ridge_regression( x = x[train_index,], y = y[train_index,], lambda = lambda, shrink_beta = shrink)$Coefficients
      intercept = coef[1]
      coef = coef[-1]
    }
    else if (type == "nlm")
    {
      coef = nlm(p = c(0,0,0,0), fn, lambda = lambda, y = y[train_index,], x = x[train_index,])$estimate
      intercept = coef[1]
      coef = coef[-1]
    }
    
    
    
    # Predicting data of validation set
    predicted = x[validation_index,] %*% coef + intercept
    
    # actual 
    actual = y[validation_index,]
    
    # calculating mse for given fold
    MSE[i] = mean ((actual - predicted)^2)
  }
    # cv error
    CV = mean(MSE)
    return(CV)
}
  

#' @description This function will help choose optimal alpha for ridge regression
#'
#' @param Input1: x, y : Data for fitting the model
#' @param Input2: fn: function to be minimized
#' @param Input3: lambda: lambda range
#' @param Input4: K: amount of folds for crossvalidation
#' @param ...: type and shrink parameters for my_ridge regression
#' 
#' @return List with Lambdas and corresponding MSEs, $lambda.min is available
#' 
my_best_lambda = function(x, y, fn, lambda, K = 10, type = "optim", shrink = TRUE)
{
  # initializing rezult table
  results = list("lambda" = lambda, "cv.mse" = rep(0, length(lambda)), "lambda.best" = 0)
  
  ### main loop ###
  for (i in 1:length(results$lambda))
  {
    # Storing results of cross-validation for current lambda
    results$cv.mse[i] = my_kfold_ridge(x, y, K = K, lambda = results$lambda[i], f = fn, type, shrink)
  }
  
  # best lambda - the one that returns smallest mse
  results$lambda.best = results$lambda[which.min(results$cv.mse)]
  
  return (results)
}



#### Results ####

# According to my_ridge
results.ridge = my_best_lambda(x, y, f, lambda = lambdas, type = "ridge")
# best lambda
results.ridge$lambda.best
# min mse achieved
results.ridge$cv.mse[which.min(results.ridge$cv.mse)]
# corresponding coeffients
my_ridge_regression(x, y, results.ridge$lambda.best)

# According to optim
results.optim = my_best_lambda(x, y, f, lambda = lambdas, type = "optim")
# best lambda according to optim
results.optim$lambda.best
# min mse achieved with optim
results.optim$cv.mse[which.min(results.optim$cv.mse)]
# corresponding coeffients
optim(par = c(0,0,0,0), f, lambda = results.optim$lambda.best, y = y, x = x)$par

# According to nlm
results.nlm = my_best_lambda(x, y, f, lambda = lambdas, type = "nlm")
# best lambda according to nlm
results.nlm$lambda.best
# min mse achieved with nlm
results.nlm$cv.mse[which.min(results.nlm$cv.mse)]
# corresponding coeffients
nlm(p = c(0,0,0,0), f, lambda = results.nlm$lambda.best, y = y, x = x)$estimate

# According to ridhe regression
cv = cv.glmnet(x, y, alpha = 0, lambda = lambdas, nfolds = 10, standatize = FALSE)
# Best lambda
cv$lambda.min
# Corresponding coefficients
predict(ridge.test, s = cv$lambda.min, alpha = 0, type = "coefficients")


par(mfrow = c(2,2))
#### Plotting results ####
plot(log(results.ridge$lambda), results.ridge$cv.mse, col = "blue", type = "l", lwd = 2, main = "My Ridge", xlab = c("Lambda", results.ridge$lambda.best))
points(log(results.ridge$lambda.best), results.ridge$cv.mse[which.min(results.ridge$cv.mse)], col = "red", pch = "X")

plot(log(results.optim$lambda), results.optim$cv.mse, col = "yellow", type = "l", lwd = 2, main = "Optim", xlab = c("Lambda", results.optim$lambda.best))
points(log(results.optim$lambda.best), results.optim$cv.mse[which.min(results.optim$cv.mse)], col = "red", pch = "X")

plot(log(results.nlm$lambda), results.nlm$cv.mse, col = "magenta", type = "l", lwd = 2, main = "nlm", xlab = c("Lambda", results.nlm$lambda.best))
points(log(results.nlm$lambda.best), results.nlm$cv.mse[which.min(results.nlm$cv.mse)], col = "red", pch = "X")

plot(cv, main = "Ridge CV", xlab = c("Lambda", cv$lambda.min))


#' Despite the fact that none of the methods above were able to (as expected from ridge regression) 0 out x2 and x3 in our model (as y clearly does not depend on either)
#' all methods correctly shrinked them towards zero, with optim and nlm apply more penalty to the coefficients. Optim also has 
#' lowest mean(mse) out of three applied methods and is closer to the one r provides us with. However it doesnt spot lowest minumum optimum lambda for the
#' module (maybe for the best as we know that there is no correlation between x2, x3 and y) All methods seem to correctly estimate intecept close to 0 + sd(u) but
#' there is some deviation in this. Since we did not center data explicitly, my_ridge also brings intercept closer to 1.
#' nlm apears to return smaller mse and better overall coefficients for bigger range of lambda

my_ridge_regression(x, y, lambda = 10^15)
optim(par = c(0,0,0,0), f, lambda = 10^15, y = y, x = x)$par
nlm(p = c(0,0,0,0), f, lambda = 10^15, y = y, x = x)$estimate
predict(ridge.test, s = 10^15, alpha = 0, type = "coefficients")

#' Also it appears that when lambde gets higher optim is the only method besides glmnet that returns ridge-like coefficients with
#' unpenalized intercept and close to but not exactly zero coefficients. 
#' nlm is clearly follows different approach.

# if we decide to not shrink intercept in my_ridge approach, however, we get more sensible results
results.ridge.intercept = my_best_lambda(scale(x), scale(y), lambda = lambdas, K = 10, type = "ridge", shrink = FALSE)
results.ridge.intercept$lambda.best
my_ridge_regression(x, y, results.ridge.intercept$lambda.best)
# That seems to caught the behaviour of the model the better and return reasonable mse.

# Plot again
plot(log(results.ridge.intercept$lambda), results.ridge.intercept$cv.mse, col = "blue", type = "l", lwd = 2, main = "My Ridge", xlab = c("Lambda", results.ridge.intercept$lambda.best))
points(log(results.ridge.intercept$lambda.best), results.ridge.intercept$cv.mse[which.min(results.ridge.intercept$cv.mse)], col = "red", pch = "X")

plot(log(results.optim$lambda), results.optim$cv.mse, col = "yellow", type = "l", lwd = 2, main = "Optim", xlab = c("Lambda", results.optim$lambda.best))
points(log(results.optim$lambda.best), results.optim$cv.mse[which.min(results.optim$cv.mse)], col = "red", pch = "X")

plot(log(results.nlm$lambda), results.nlm$cv.mse, col = "magenta", type = "l", lwd = 2, main = "nlm", xlab = c("Lambda", results.nlm$lambda.best))
points(log(results.nlm$lambda.best), results.nlm$cv.mse[which.min(results.nlm$cv.mse)], col = "red", pch = "X")

plot(cv, main = "Ridge CV", xlab = c("Lambda", cv$lambda.min))



# all lines together
par(mfrow=c(1,1))
plot(cv)
lines(log(results.ridge.intercept$lambda), results.ridge.intercept$cv.mse, col = "blue", type = "l", lwd = 2, main = "My Ridge", xlab = c("Lambda", results.ridge.intercept$lambda.best), add = TRUE)
points(log(results.ridge.intercept$lambda.best), results.ridge.intercept$cv.mse[which.min(results.ridge.intercept$cv.mse)], col = "red", pch = "X")

lines(log(results.optim$lambda), results.optim$cv.mse, col = "yellow", type = "l", lwd = 2, main = "Optim", xlab = c("Lambda", results.optim$lambda.best), add = TRUE)
points(log(results.optim$lambda.best), results.optim$cv.mse[which.min(results.optim$cv.mse)], col = "red", pch = "X")

lines(log(results.nlm$lambda), results.nlm$cv.mse, col = "magenta", type = "l", lwd = 2, main = "nlm", xlab = c("Lambda", results.nlm$lambda.best), add = TRUE)
points(log(results.nlm$lambda.best), results.nlm$cv.mse[which.min(results.nlm$cv.mse)], col = "red", pch = "X")


