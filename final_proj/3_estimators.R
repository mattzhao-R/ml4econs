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

RJIVE_function = function(OUTCOME_VECTOR, TREATMENT_VARIABLE_MATRIX, INSTRUMENT_VARIABLE_MATRIX){
  # Inputs
  # Here I define the four inputs into RJIVE.
  # The "n" long outcome vector
  y_vector = OUTCOME_VECTOR
  # The "n x G" treatment variable matrix
  x_matrix = as.matrix(TREATMENT_VARIABLE_MATRIX)
  number_of_observations = dim(x_matrix)[1]
  number_of_variables = dim(x_matrix)[2]
  # The "n x K" matrix of instruments
  z_matrix = as.matrix(INSTRUMENT_VARIABLE_MATRIX)
  number_of_instruments = dim(z_matrix)[2]
  # The scalar penalty parameter used in the Ridge Regularization
  # First, we regress our X matrix on any columns of our original
  # matrix that are included in both X and Z.
  regression_of_x_on_controls = lm(lprice ~ ., as.data.frame(x_matrix))
  constant_of_proportionality = sd(regression_of_x_on_controls$residuals)
  lambda_scalar_sqrt = constant_of_proportionality * sqrt(number_of_instruments)
  
  
  # Step 1: Ridge Regression
  # First, we use our lambda scalar to create
  # a penalty matrix for the Ridge Regularization
  lambda_matrix = diag(number_of_instruments)
  lambda_matrix = lambda_scalar_sqrt * lambda_matrix
  
  # Next, we create a dummy matrix
  # and use it to initialize a list of matrices.
  # This matrix will end up being replaced,
  # so its values and shap doesn't matter.
  empty_matrix = matrix(1:6, ncol=2)
  pi_matrix_minus_i_hat_list = list(empty_matrix)
  # Next, we are going to iterate through all of our observations.
  for (i in 1:number_of_observations){
    # In each iteration, we first create the
    # first term that will go into our Pi(-i) matrix.
    temp_first_term_matrix = (t(z_matrix) %*% z_matrix) 
    - as.matrix(z_matrix[i,]) %*% as.matrix(t(z_matrix[i,]))
    + (t(lambda_matrix) %*% lambda_matrix)
    # Next, we invert the first term.
    temp_first_term_matrix_inverse = solve(temp_first_term_matrix)
    # Next, we construct the second term that will go into our matrix.
    temp_second_term_matrix = (t(z_matrix) %*% x_matrix)
    - as.matrix(z_matrix[i,]) %*% as.matrix(t(x_matrix[i,]))
    # Next, we matrix multiply the two terms together to get our matrix.
    temp_pi_matrix_minus_i_hat = temp_first_term_matrix_inverse %*% temp_second_term_matrix
    # Finally, we add our matrix to our list of matrices.
    pi_matrix_minus_i_hat_list[[i]] = temp_pi_matrix_minus_i_hat
  }
  
  
  # Step 2: JIVE
  # Now that we have our list of Pi(-i) matrices, we can preform JIVE.
  # First, we initialize two zero matrices to represent
  # our first and second terms in formula for our delta vector.
  # As we iterate, we will add our temporary matrices to these two matrices.
  d_tilda_first_term_sum_matrix = matrix(0, nrow = number_of_variables, ncol = number_of_variables)
  d_tilda_second_term_sum_matrix = matrix(0, nrow = number_of_variables, ncol = 1)
  # Next, we begin iterating through our observations.
  for (i in 1:number_of_observations){
    # First, we create the temporary matrix that corresponds
    # to the first term in our equation for observation i.
    d_tilda_first_matrix_temp = t(pi_matrix_minus_i_hat_list[[i]]) %*% as.matrix(z_matrix[i,]) %*% as.matrix(t(x_matrix[i,]))
    # Next, we add the temp matrix to our summation matrix
    d_tilda_first_term_sum_matrix = d_tilda_first_term_sum_matrix + d_tilda_first_matrix_temp
    # Next, we do the same process for the temporary matrix for
    # the second term corresponding to the current i value.
    d_tilda_second_matrix_temp = t(pi_matrix_minus_i_hat_list[[i]]) %*% as.matrix(z_matrix[i,]) %*% as.matrix(y_vector[i])
    # Then, we add our temporary matrix to our summation matrix
    d_tilda_second_term_sum_matrix = d_tilda_second_term_sum_matrix + d_tilda_second_matrix_temp
  }
  # Next, we invert the summation matrix for the first term.
  d_tilda_first_term_sum_matrix_inverse = solve(d_tilda_first_term_sum_matrix)
  # Finally, we matrix multiply the two terms together to get our delta tilda vector.
  d_tilda_vector = d_tilda_first_term_sum_matrix_inverse %*% d_tilda_second_term_sum_matrix
  
  
  # Step 3: Standard Error Estimator
  y_hat = x_matrix %*% d_tilda_vector
  y_minus_y_hat_squared = rep(0, number_of_observations)
  for (i in 1:number_of_observations){
    y_minus_y_hat_squared[i] = (y[i] - y_hat[i])^2 
  }
  sample_regression_sd = sqrt(sum(y_minus_y_hat_squared) / (number_of_observations - 1))
  sample_regression_se = sample_regression_sd / sqrt(number_of_observations)
  
  
  # Step 4: Variance Estimator
  p_lambda_matrix_middle_term = t(z_matrix) %*% z_matrix + t(lambda_matrix) %*% lambda_matrix
  p_lambda_matrix_middle_term_inverse = solve(p_lambda_matrix_middle_term)
  p_lambda_matrix = z_matrix %*% p_lambda_matrix_middle_term_inverse %*% t(z_matrix)
  
  x_tilda_matrix = x_matrix
  for (i in 1:number_of_observations){
    x_tilda_matrix[i,] = (x_matrix[i,]) / (1 - p_lambda_matrix[i,i])
  }
  h_tilda_matrix = t(x_matrix) %*% p_lambda_matrix %*% x_tilda_matrix
  for (i in 1:number_of_observations){
    h_tilda_matrix = h_tilda_matrix - (as.matrix(x_matrix[i,]) %*% p_lambda_matrix[i,i] %*% t(x_tilda_matrix[i,]))
  }
  h_tilda_inverse_matrix = solve(h_tilda_matrix)
  
  x_bar_matrix = p_lambda_matrix %*% x_matrix
  
  sigma_tilda_first_term_entries = rep(0, number_of_variables^2)
  sigma_tilda_first_term_matrix = matrix(data = sigma_tilda_first_term_entries, nrow = number_of_variables)
  for (i in 1:number_of_observations){
    sigma_tilda_first_term_first_term_temp = x_bar_matrix[i,] %*% t(x_bar_matrix[i,])
    sigma_tilda_first_term_second_term_temp = x_matrix[i,] %*% (p_lambda_matrix[i,i] * t(x_bar_matrix[i,]))
    sigma_tilda_first_term_third_term_temp = x_bar_matrix[i,] %*% (p_lambda_matrix[i,i] * t(x_matrix[i,]))
    sigma_tilda_first_term_temp = (sigma_tilda_first_term_first_term_temp - sigma_tilda_first_term_second_term_temp - sigma_tilda_first_term_third_term_temp)
    xi_tilda_first_term_temp = (1 - p_lambda_matrix[i,i])^-2 * (y_vector[i] - (t(x_matrix[i,]) %*% d_tilda_vector))^2
    sigma_tilda_first_term_matrix = sigma_tilda_first_term_matrix + sigma_tilda_first_term_temp * xi_tilda_first_term_temp[1, 1]
  }
  
  z_tilda_matrix_second_term = t(z_matrix) %*% z_matrix + t(lambda_matrix) %*% lambda_matrix
  z_tilda_matrix_second_term = solve(z_tilda_matrix_second_term)
  z_tilda_matrix = z_matrix %*% z_tilda_matrix_second_term
  
  sigma_tilda_second_term_entries = rep(0, number_of_variables^2)
  sigma_tilda_second_term_matrix = matrix(data = sigma_tilda_second_term_entries, nrow = number_of_variables)
  for (k in 1:number_of_instruments){
    for (l in 1:number_of_instruments){
      sigma_tilda_second_term_first_term_entries_temp = rep(0, number_of_variables)
      sigma_tilda_second_term_first_term_matrix_temp = matrix(data = sigma_tilda_second_term_first_term_entries_temp, nrow = number_of_variables)
      for (i in 1:number_of_observations){
        sigma_tilda_second_term_first_term_matrix_temp = sigma_tilda_second_term_first_term_matrix_temp + as.numeric(z_tilda_matrix[i, k]) * as.numeric(z_tilda_matrix[i, l]) * x_matrix[i,] * as.numeric((1 - p_lambda_matrix[i,i])^-1 * (y_vector[i] - (t(x_matrix[i,]) %*% d_tilda_vector)))
      }
      sigma_tilda_second_term_second_term_entries_temp = rep(0, number_of_variables)
      sigma_tilda_second_term_second_term_matrix_temp = matrix(data = sigma_tilda_second_term_second_term_entries_temp, nrow = number_of_variables)
      for (j in 1:number_of_observations){
        sigma_tilda_second_term_second_term_matrix_temp = sigma_tilda_second_term_second_term_matrix_temp + as.numeric(z_matrix[j, k]) * as.numeric(z_matrix[j, l]) * x_matrix[j,] * as.numeric((1 - p_lambda_matrix[j,j])^-1 * (y_vector[j] - (t(x_matrix[j,]) %*% d_tilda_vector)))
      }
      sigma_tilda_second_term_matrix = sigma_tilda_second_term_matrix + sigma_tilda_second_term_first_term_matrix_temp %*% t(sigma_tilda_second_term_second_term_matrix_temp)
    }
  }
  
  sigma_tilda_matrix = sigma_tilda_first_term_matrix + sigma_tilda_second_term_matrix
  
  variance_tilda = h_tilda_inverse_matrix %*% sigma_tilda_matrix %*% h_tilda_inverse_matrix
  
  return_list = list(c(1:3))
  return_list[[1]] = d_tilda_vector
  return_list[[2]] = sample_regression_se
  return_list[[3]] = variance_tilda
  
  return(return_list)
}

#### post lasso ----
# see ipynb + section in analysis using hdm

