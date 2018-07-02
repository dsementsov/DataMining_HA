##################################################
#                  Data Mining                   #
#                Home Assignment 1               #
#                  Denys Sementsov               #
#                stu110401 MNr 1001084           #
#             stu110401@mail.uni-kiel.de         #
#             denys.sementsov@outlook.com        #
##################################################

library(stats)
library(utils)
library(methods)

# Part 1: Loading the data
# setting working directory
# please set the right working directory (location of a house_data.csv file)
setwd("~/Documents/Studies/Quantitative Economics/Data Minig/Home Assigment/1")

# importing data using read.csv function, we have a header, so we set it to TRUE, 
# choose right separator and leave everything else blank as we don't have decimals or omitted variables
# in the file. 

# The separator on my computer for some reason is ';'
# Please, choose appropriate separator if it differs on windows machine

house = read.csv('house_data.csv', header = TRUE, sep = ';')   # reading file
## the variable house is supposed to have 149 observations of 6 variables

options(scipen=999)         # turns off scientific notation

####################
# Part 2: Linear Fit

# fitting linear model of the form:
# price=β0 +β1 ·bedrooms+β2 ·sqft_living+β3 ·floors+β4 ·grade+β3 ·condition+ε

lm.fit = lm(formula = price ~ . , data = house)   # fits linear model with all features
summary(lm.fit)                                   # displays summary for the model

# out of summary we see that Rsquared is ~0.45, so we can conclude that straight linear model with chosen
# parameters may not be a good fit. Also, even though overall p-value approaches zero, F statistics isnt particularly high for 
# the model, and individual p values on bedrooms, floors and condition doesnt show significance.

lm.fitNoIntercept = lm(formula = price ~ bedrooms + sqft_living + floors + grade + condition + 0, 
            data = house)   # fits the same model with no intercept (for the future reference)
summary(lm.fitNoIntercept)

# fitting the model without intercept giving us way higher Rsquared, indicating better fit on the first glance
# It maybe the evidence of that is our data is supposed to go through origin 
# (however in practice, 0 price is not that common on the market)
# Some of the p-values are still too high, and Rsquared sometimes can be not the most reliable source of a fit=> 
# we should look for the better model

####################
# Part 3: OLS function
my.OLS = function(x, y, intercept = TRUE)   
  # x - matrix - independent vars; 
  # y - vector - dependent vars; 
  # intercept - bool var - whether we want intercept. ( default = TRUE, with intercept )
{
  # Verifying input:
  # this will get us out of the function if there is a mistake in the input 
  # without generating an error when multiplying matrices
  if ( nrow(y) != nrow(x) )  
  {
    return('Input error! Make sure that number of rows of Y = number of rows of X')
  }
  # add vector of 1 in the beginning of a X matrix for the intercept
  if (intercept)
  {
      x = cbind(matrix (seq(1, by = 0,  length.out = nrow(x)), ncol = 1), x) 
      colnames(x)[1] = 'intercept'
  }
  
  #converting from data_table to matrix
  matY = as.matrix(y)
  matX = as.matrix(x)
  
  # estimating Betas through OLS formula
  # function solve() gives us an inverse, t() gives us matrix transposed, %*% multiplies matrices
  estimatedBetas = solve((t(matX) %*% matX)) %*% (t(matX) %*% matY)
  # finding RSS
  error = matY - matX %*% estimatedBetas
  resSumSq = t(error) %*% error
  # output in a list for better readability
  return(list( "Estimated Betas" = estimatedBetas, "RSS" = resSumSq))
}

# Part 3.1: Checking my.OLS function against lm

#setting up data
X = house[,-1]
Y = matrix(house[,1])

# running my.OLS function
OLSoutput = my.OLS(X,Y)
OLSoutputNoIntercept = my.OLS(X,Y, FALSE)

# checking data against lm.fit
# Checking Betas

#with intercept
OLSoutput['Estimated Betas']
coef(lm.fit)

#without intercept
OLSoutputNoIntercept['Estimated Betas']
coef(lm.fitNoIntercept)

# Checking RSS
# function that will return us RSS given lm for better readability
my.RSS = function(lmF)
{
  return(sum(residuals(lmF)^2))
}
# with intercept
my.RSS(lm.fit)
OLSoutput['RSS']
# without intercept
my.RSS(lm.fitNoIntercept)
OLSoutputNoIntercept['RSS']

####################
# Part 3: Forward Selection
# forward selection

# function for choosing next regressor given initial function

my.forwardNext = function( checkFN )
{
  regressors = names(checkFN$coefficients)[-1]
  Rsq = 0                    # we will check adjusted Rsquared and chose regressors that will return max Rsq
  for (i in 2:ncol(house))   # iterate through data sample from 2 as house[1] is a price vector
  {
       if (names(house[i]) %in% regressors) # skip regressors that are already in the model
       {
         next
       }
       # unnecessary complicated approach to construct a formula
       regr = paste(regressors, collapse = '+')  
       regr = paste(regr, names(house[i]), sep = '+')   # we iterate through house[i] that are not already
                                                        # in the model to chose next regressor
       #print(regr)
       nextFit = lm(as.formula(paste('price ~ ', regr, sep = '')), data = house) # iterating through all lm
       #print(summary(nextFit)$adj.r.squared)
       if (summary(nextFit)$adj.r.squared >= Rsq)    # save the model only if Rsq is higher than that 
                                                     # of a previous lm
       {
         Rsq = summary(nextFit)$adj.r.squared
         bestFit = nextFit
         chosenRegr = names(house[i])
       }
  }
  return(bestFit)
}


forwardSelection = lm(price ~ 1, data = house)  # initial model with just intercept
allForwardModels = list()                       # list that will contain results of every step with adjRsq / Rsq

for (i in 1:(ncol(house)- 1))     # loop through all regressors
{
  forwardSelection = my.forwardNext(forwardSelection) # passing last step best function as a baseline for next
  allForwardModels[paste(names(forwardSelection$coefficients), collapse = "+")] = paste(summary(forwardSelection)$adj.r.squared, summary(forwardSelection)$r.squared, sep = '   /    ')
}

###########################
# Part 4: Backward Selection 

# function for choosing best regressor for improving adjRsq given previous model
my.backwardNext = function( checkFN )
{
  regressors = names(checkFN$coefficients)[-1]  # all but intercept
  Rsq = 0                    # we will check Rsquared and chose regressors that will return max Rsq
  if (length(regressors) == 0)       # breaking the loop when we pass model with not enough regressors (e.g only intercept)
  {
    return()
  }
  for (i in 1:length(regressors))   # iterate through data sample from 2 as house[1] is a price vector
  {
    regr = paste(regressors, collapse = '+')  
    regr = paste(regr, regressors[i], sep = '-')   # now we take all except one regressor
    #print(regr)
    nextFit = lm(as.formula(paste('price ~ ', regr, sep = ' ')), data = house) # iterating through all lm
    #print(summary(nextFit)$adj.r.squared)
    if (summary(nextFit)$adj.r.squared >= Rsq) 
    {
      Rsq = summary(nextFit)$adj.r.squared
      bestFit = nextFit
      chosenRegr = names(house[i])
    }
  }
  return(bestFit)
}

# initial model is lm with all variables
backwardSelection = lm(price ~ ., data = house)
allBackwardModels = list()
for (i in 1:(ncol(house)- 1))     # looping through all regressors
{
  backwardSelection = my.backwardNext(backwardSelection) # passing last step best function as a baseline for next
  allBackwardModels[paste(names(backwardSelection$coefficients), collapse = "+")] = paste(summary(backwardSelection)$adj.r.squared, summary(backwardSelection)$r.squared, sep = '   /    ')
}

# Part 5: Concluding output (which models you choose based on forward stepwise selection and 
# backward stepwise selection)


allForwardModels
# Based on forward selection and judging by adjusted Rsquared, 
# we would choose (Intercept)+sqft_living+grade model with adjRsq = 0.4374391, that is actually
# does not return highest Rsq value
bestForward = lm(price ~ sqft_living + grade, data = house)
summary(bestForward)
# checking results
stepForward = step(lm(price ~ . , data = house), direction = 'forward')
summary(stepForward)
# here, we are advised to add all regressors to the model, increasing Rsq on 0.001 and trading
# 0.12 adjusted Rsquared, which is not worth it for given data.
allBackwardModels
# surprisingly, backward step selection yields same result as a forward step selection.
bestBackward = lm(price ~ . - condition - floors - bedrooms, data = house)
summary(bestBackward)
# checking results
stepBackward = step(lm(price ~ . , data = house), direction = 'backward')
summary(stepBackward)
# backward selection given by R returns same result

# both methods chose regressors with lowest p-value
# given gathered data, we can try to find the best fit for the data, combining results of backward/forward
# selections and forcing the regression to go through the origin
conclusion = lm(price ~ 0 + sqft_living + grade, data = house)
summary(conclusion)
# with high Rsq and F-statistics


###############################
#Also, we couldve wrap our loops in the function to work much like Rs step function
###############################

my.Step = function(input, Y, steps = 1000, direction = "forward", maxRsq = 0) 
{
  # check number of steps
  if (steps > ncol(input)-1)
  {
    steps = ncol(input)-1
  }
  # check the direction of step
  if (direction == "forward")
  {
    # initial function
    output = lm(Y ~ 1, data = input) 
    for (i in 1:steps)    
    {
      # new fit with new parameters
      output = my.forwardNext(output) 
      # check whether this fit is better
      if (maxRsq <= summary(output)$adj.r.squared)
      {
        best = output
        maxRsq = summary(output)$adj.r.squared
      }
    }
  }
  # check the direction of step
  if (direction == "backward")
  {
    # initial function
    output = lm(Y ~ ., data = input)
    for (i in 1:steps)    
    {
      output = my.backwardNext(output)
      # check whether this fit is better
      if (maxRsq <= summary(output)$adj.r.squared)
      {
        best = output
        maxRsq = summary(output)$adj.r.squared
      }
    }
    
  }
  return(best)
}

bBest = my.Step(house, house$price, direction = "backward")
summary(bBest)
fBest = my.Step(house, house$price, direction = "forward")
summary(fBest)
