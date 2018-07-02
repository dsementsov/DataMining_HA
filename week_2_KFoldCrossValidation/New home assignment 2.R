# New home assignment
# Resamling

# 1. Function for generating random samples of a particular sample size for the model
# y = 1 + x^3 + u, where x is stN and u is stN.

gen_sim_data = function(sample_size) {
  u = rnorm(sample_size)
  x = rnorm(sample_size)
  y = 1 + x ^ 3 + u
  data.frame(x, y)
}

# Simulate a single dataset, split it into a train and 
# validation set. Here, the validation set is 20% of the data.
set.seed(90)
sim_data = gen_sim_data(sample_size = 200)
sim_idx  = sample(1:nrow(sim_data), 160)
sim_trn  = sim_data[sim_idx, ]
sim_val  = sim_data[-sim_idx, ]

# Plot this training data, as well as the true regression function.

plot(y ~ x, data = sim_trn, col = "dodgerblue", pch = 20)
grid()
curve(1 + x ^ 3, add = TRUE, col = "black", lwd = 2)


# 2. Write functions we need in 3.
# write a small function to calculate MSE
calc_mse = function(actual, predicted) {
  mean((actual - predicted) ^ 2)
}

# write a fuction to calculate a K-fold cross validation statistic (only for K>=2)
my.kfold = function(model,data,K){
  # check if the sample size is a multiple of K
  n = max(dim(data))
  folds = c()
  if(floor(n/K)==ceiling(n/K)){
    nf = n/K # number of obs in a fold
    for(i in 1:K){
      folds = c(folds,rep(i,times=nf))
    }
    MSE = c()
    for(i in 1:K){
      train = (folds!=i)
      valid = (folds==i)
      fit = lm(model,data=data[train,])
      y.pred = predict(fit,newdata=data[valid,])
      MSE[i] = calc_mse(data$y[valid], y.pred)
    }
    KF = mean(MSE)
  }else{
    cat("The data cannot be split evenly into ",K, "folds! \n")
    KF = NA
  }
  return(KF)
}

# 3. Simulation study 
# Compare:
# 3.1. use the whole dataset to calculate MSE and choose the model
# 3.2. use validation approach
# 3.3. use K-10 approach

num_sims = 100
num_degrees = 10
K = 10

whole_set = matrix(0, ncol = num_degrees, nrow = num_sims)
val_mse = matrix(0, ncol = num_degrees, nrow = num_sims)
cv_mse = matrix(0, ncol = num_degrees, nrow = num_sims)

set.seed(90)
choice_whole = c()
choice_val = c()
choice_cv = c()

for (i in 1:num_sims) {
  # simulate data
  sim_data = gen_sim_data(sample_size = 200)
  # set aside validation set
  sim_idx = sample(1:nrow(sim_data), 160)
  sim_trn = sim_data[sim_idx, ]
  sim_val = sim_data[-sim_idx, ]
  # fit models and store MSEs
  for (j in 1:num_degrees) {
    #fit model
    fit = glm(y ~ poly(x, degree = j), data = sim_trn)
    # calculate error
    whole_set[i, j] = calc_mse(actual = sim_trn$y, predicted = predict(fit, sim_trn))
    val_mse[i, j] = calc_mse(actual = sim_val$y, predicted = predict(fit, sim_val))
    cv_mse[i, j] = my.kfold(y ~ poly(x, degree = j),sim_data,K)
  }
  choice_whole[i] = which.min(whole_set[i, ])
  choice_val[i] = which.min(val_mse[i, ])
  choice_cv[i] = which.min(cv_mse[i, ])
  
}

whole.mean = colMeans(whole_set)
val.mean = colMeans(val_mse)
cv.mean = colMeans(cv_mse)

whole.sd = apply(whole_set,2,sd)
val.sd = apply(val_mse,2,sd)
cv.sd = apply(cv_mse,2,sd)

options(scipen=999)
results = data.frame(whole.mean, whole.sd, val.mean, val.sd,  cv.mean, cv.sd)

# 4. Tables and plots.
whole_freq <- table(choice_whole)/100
val_freq <- table(choice_val)/100
cv_freq <- table(choice_cv)/100

whole_freq
val_freq
cv_freq

install.packages("ggplot2")
library(ggplot2)
choice = data.frame(choice_whole, choice_val, choice_cv)
names_for_graph = c('training error', 'validation set error', 'CV error')
for (i in 1:3){
  print(qplot(choice[i],
        geom="histogram",
        binwidth=1,  
        main= substitute(paste("Model Chosen:", a), list(a = names_for_graph[i])),
        xlab="Polynomial Degree",
        ylab = 'Times Chosen',
        fill=I("blue"), 
        col=I("red"),
        breaks=c(0:10),
        xlim = c(0,10)
  ))
}


# 5. Explain your results.
# differences: whole is just wrong
# validation compared to cv: less variance, better selections