# Last updated: Feb 29, 2023

#### jive steinIV ----
##################################################
### JIVE
jive.est <- function(y,X,Z,SE=FALSE,n.bt=100){
  n  <- length(y);
  k  <- dim(X)[2];
  out <- list();
  B   <- jive.internal(y,X,Z);
  ### SE
  if(SE==TRUE){
    k       <- length(B);
    var.hat <- matrix(0,k,k);
    for(b in 1:n.bt){
      bt <- sample(1:n,replace=TRUE);
      y.bt <- y[bt];
      X.bt <- as.matrix(X[bt,]);
      Z.bt <- as.matrix(Z[bt,]);
      B.bt <- jive.internal(y.bt,X.bt,Z.bt);
      var.hat <- var.hat + ((B.bt-B)%*%t(B.bt-B))/(n.bt-k); 
    }#n.bt     
    ### Out.
    out$est <- B;
    out$se  <- sqrt(diag(var.hat));
    out$var <- var.hat;
  }else{
    out$est <- B;
  }#fi
  return(out)
}#jive.est

##################################################
### JIVE (Internal).
jive.internal <- function(y,X,Z){
  n <- length(y);
  k <- dim(X)[2];
  Xj <- matrix(0,n,k);
  iZZ<- solve(t(Z)%*%Z);
  Ga.hat <- (iZZ)%*%(t(Z)%*%X);
  h <- vector("numeric",n);
  for(i in 1:n){
    h[i]   <- t(Z[i,])%*%iZZ%*%Z[i,];    
    Xj[i,] <- (t(Z[i,])%*%Ga.hat - h[i]*X[i,])/(1-h[i]);    
  }#n
  XXj <- t(Xj)%*%X;
  iXXj<- solve(XXj);   
  #iXXj<- ginv(XXj);   ### Use Generalized inverse.
  B   <- (iXXj)%*%(t(Xj)%*%y);
  return(B)
}#jive.internal


#### rjive ----

library(fixest)

rjive.est <- function(y,w,d,z,gamma,df){
  # Inputs
  # The "n" long outcome vector
  y_vector = y
  # The "n x p" treatment variable matrix
  ## d is the endogenous regressor
  ## w are the excluded instruments/controls that appear in both the second and first stages
  x = cbind(d,w)
  x_matrix = x
  number_of_observations = dim(x_matrix)[1]
  number_of_variables = dim(x_matrix)[2] 
  # The "n x K" matrix of instruments
  z_matrix = z
  number_of_instruments = dim(z_matrix)[2]
  
  
  # Step 1: Ridge Regression
  # First, we use our lambda scalar to create
  # a penalty matrix for the Ridge Regularization
  # The scalar penalty parameter used in the Ridge Regularization
  indiv_cnty_fe <- colnames(df)[str_detect(colnames(df), "cntyfe_")]
  indiv_time_fe <- colnames(df)[str_detect(colnames(df),'timefe_')]
  model = feols(
    lprice ~ 1,
    fixef = c(indiv_time_fe,indiv_cnty_fe),
    data = df
  )
  x_res = model$residuals
  C = sd(x_res)
  gamma_sqrt = C * sqrt(number_of_instruments)
  
  LAMBDA = gamma_sqrt * diag(number_of_instruments)
  
  # Next, we create a dummy matrix
  # and use it to initialize a list of matrices.
  # This matrix will end up being replaced,
  # so its values and shap doesn't matter.
  PI_minusi_list = list()
  # Next, we are going to iterate through all of our observations.
  for (i in 1:number_of_observations){
    # In each iteration, we first create the
    # first term that will go into our Pi(-i) matrix.
    temp_first_term_matrix = (t(z_matrix) %*% z_matrix) - as.matrix(z_matrix[i,]) %*% as.matrix(t(z_matrix[i,])) + (t(LAMBDA) %*% LAMBDA)
    # Next, we invert the first term.
    temp_first_term_matrix_inverse = solve(temp_first_term_matrix)
    # Next, we construct the second term that will go into our matrix.
    temp_second_term_matrix = (t(z_matrix) %*% x_matrix) - as.matrix(z_matrix[i,]) %*% as.matrix(t(x_matrix[i,]))
    # Next, we matrix multiply the two terms together to get our matrix.
    temp_pi_matrix_minus_i_hat = temp_first_term_matrix_inverse %*% temp_second_term_matrix
    # Finally, we add our matrix to our list of matrices.
    PI_minusi_list[[i]] = temp_pi_matrix_minus_i_hat
  }
  
  # Step 2: JIVE
  # Now that we have our list of Pi(-i) matrices, we can perform JIVE.
  # Empty matrix to sum the inverted p x p matrix (LHS part of delta_tilde formula)
  LHS_dtilde = matrix(0,number_of_variables,number_of_variables)
  # Empty vector to sum the p x 1 matrix (RHS part of delta_tilde formula)
  RHS_dtilde = matrix(0,number_of_variables,1)

  # Next, we begin iterating through our observations.
  for (i in 1:number_of_observations){
    # First, we create the temporary matrix that corresponds
    # to the first term in our equation for observation i.
    d_tilda_first_matrix_temp = t(PI_minusi_list[[i]]) %*% as.matrix(z_matrix[i,]) %*% as.matrix(t(x_matrix[i,]))
    # Next, we add the temp matrix to our summation matrix
    LHS_dtilde = LHS_dtilde + d_tilda_first_matrix_temp
    # Next, we do the same process for the temporary matrix for
    # the second term corresponding to the current i value.
    d_tilda_second_matrix_temp = t(PI_minusi_list[[i]]) %*% as.matrix(z_matrix[i,]) %*% as.matrix(y_vector[i])
    # Then, we add our temporary matrix to our summation matrix
    RHS_dtilde = RHS_dtilde + d_tilda_second_matrix_temp
  }
  # Now that we have the summed matrices, we can get delta_tilde
  delta_tilde = ginv(LHS_dtilde) %*% RHS_dtilde
  
  # Step 3: Variance Estimator
  
  
  
  return(c(delta_tilde,delta_var))
}

#### post lasso ----
# see ipynb + section in analysis using hdm

