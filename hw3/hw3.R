library(tidyverse)
library(splines)
library(plotrix)

DGP <- function(){
  V =  rexp(200,1/2)
  X = 1 + V
  epsilon = rnorm(200,0,sqrt(0.5))
  y = 2 + 2 * log(X) + epsilon
  list(X, y)
}

# Selecting Optimal Number of Knots ----

spline <- function(X,y,x, numknots){
  Z = bs(X, degree=3, df = 4 + numknots, intercept = T)
  beta = solve(t(Z) %*% Z) %*% t(Z) %*% y
  z_x = predict(Z,x)
  list(Z, beta, z_x)
}

cv_spline <- function(X,y,numknots){
  tot_err = 0
  for (i in 1:(length(X)-1)) {
    bspl = spline(X[-i],y[-i],X[i],numknots)
    Z = bspl[[1]]
    beta = bspl[[2]]
    z_x = bspl[[3]]
    e_hat = (y[i] - z_x %*% beta)
    h = t(Z[i,]) %*% solve(t(Z) %*% Z) %*% Z[i,]
    err = (e_hat^2) / (1 + h)^2
    tot_err = tot_err + err
  }
  list(numknots,(tot_err/length(X)))
}

data = DGP()
X = data[[1]]
y = data[[2]]
num_knots <- c()
cv_err <- c()
for (knots in 1:30){
  result <- cv_spline(X,y,knots)
  num_knots <- c(num_knots,result[[1]])
  cv_err <- c(cv_err,result[[2]])
}

output = tibble(num_knots = num_knots, error = cv_err)
# View(output %>% arrange(error,descending=F))
opt_numknots = (output %>% arrange(error,descending=F))$num_knots[1]

# Constructing CI ----

opt_Z = bs(X, degree=3, df = 4 + opt_numknots, intercept = T)
opt_beta = solve(t(opt_Z) %*% opt_Z) %*% t(opt_Z) %*% y

pred_spline <- function(Z,x){
  zx = predict(Z,x)
  pred_y = zx %*% opt_beta
  pred_y
}

var_spline <- function(Z,x){
  a = predict(Z,x)
  Q <- t(Z) %*% Z / 200
  omega <- 0
  for (i in 1:length(y)) {
    e_hat <- y[i]- t(Z[i,]) %*% opt_beta
    omega <- omega + ((e_hat^2)[[1]] * (Z[i,] %*% t(Z[i,])))
  }
  omega <- omega / 200
  (a %*% Q %*% omega %*% Q %*% t(a))[[1]]
}

# Creating and Plotting Spline CIs with true values ----

pred_ci <- function(Z,x){
  center = pred_spline(Z,x)
  upper = center + 1100 * sqrt(var_spline(opt_Z,X) / length(X))
  lower = center - 1100 * sqrt(var_spline(opt_Z,X) / length(X))
  list(upper,center,lower)
}

x_grid = seq(1,5,0.2)
ci_data = as.data.frame(do.call(cbind,pred_ci(opt_Z,x_grid)))
colnames(ci_data) <- c('Upper','Center','Lower')

png(filename = '/Users/matthewzhao/Documents/Coding/ml4econs/hw3/spline_ci.png',
    width = 12, height = 10, units = 'cm',res=300)
plotCI(x_grid, ci_data$Center,
       ui=ci_data$Upper, 
       li=ci_data$Lower,
       col="red",xlab="X",ylab="m(X)")
points(X,y,col="blue")
dev.off()

# Plotting LL CIs with true values ----

llci = read.csv('hw3/ll_ci.csv')
lltruedata = read.csv('hw3/ll_nn_true_data.csv')

png(filename = '/Users/matthewzhao/Documents/Coding/ml4econs/hw3/ll_ci.png',
    width = 12, height = 10, units = 'cm',res=300)
plotCI(x_grid, llci$yhat,
       ui=llci$upper, 
       li=llci$lower,
       col="red",xlab="X",ylab="m(X)")
points(lltruedata$TrueX,lltruedata$TrueY,col="blue")
dev.off()

# Plotting NN CIs with true values ----

nn_vals = read.csv('hw3/nn_pred.csv')
nntruedata = read.csv('hw3/ll_nn_true_data.csv')

png(filename = '/Users/matthewzhao/Documents/Coding/ml4econs/hw3/nn_pred.png',
    width = 12, height = 10, units = 'cm',res=300)
plotCI(x_grid, nn_vals$yhat,
       ui=nn_vals$yhat, 
       li=nn_vals$yhat,
       col="red",xlab="X",ylab="m(X)")
points(nntruedata$TrueX,nntruedata$TrueY,col="blue")
dev.off()
