# Last updated: Feb 27, 2023

#' Jackknife IV with Bootstrapped SEs
#'
#' This function returns point estimates and the variance-covariance matrix.
#' Written for data.table arguments.
#' @param Yvar Dependent variable
#' @param Xvar Endogenous regressor
#' @param Inc Included instruments
#' @param Exc Excluded instruments
#' @param data Of type data.table
#' @param intercept Logical, if including 1 as model intercept
#' @param tolerance Tolerance level for matrix inversion
#' @param n_boot Number of bootstrapped samples
#' @keywords jackknife linear bootstrap iv
#' @export
#' @examples
#' my.jackknife_iv_boot()

my.jackknife_iv_boot <- function(Yvar, Xvar, Inc, Exc, data, intercept=TRUE,
                                 tolerance=1e-16, n_boot=100) {
  # Select data and add intercept
  Y <- data[,..Yvar]
  if (length(Inc) > 0) { # If there are included instruments besides intercept
    X <- data[,c(..Inc, ..Xvar)] # Must be char vector
    Z <- data[,c(..Inc, ..Exc)] # Must be char vectors
  } else {
    X <- data[,..Xvar] # Must be char vector
    Z <- data[,..Exc] # Must be char vectors
  }
  if (intercept) {
    X <- X[,.intercept:=1]
    setcolorder(X, c(".intercept"))
    Z <- Z[,.intercept:=1]
    setcolorder(Z, c(".intercept"))
  }
  
  # Construct jackknife instrument
  ZtZinv <- solve(t(Z) %*% as.matrix(Z), tol=tolerance)
  h_i <- diag(as.matrix(Z) %*% ZtZinv %*% t(Z))
  Pi_hat <- t(X) %*% as.matrix(Z) %*% ZtZinv
  X_jack <- (1-h_i)^(-1)*(as.matrix(Z) %*% t(Pi_hat)  - h_i * X)
  
  # Coefficients
  XtXinv_jack <- solve(t(X_jack) %*% as.matrix(X), tol=tolerance)
  XtY_jack <- t(X_jack) %*% as.matrix(Y)
  beta <- XtXinv_jack %*% XtY_jack
  colnames(beta) <- c("coefs")
  
  # Predicted values, residuals, and squared residuals
  Yhat <- as.matrix(X) %*% as.matrix(unlist(beta))
  U <- as.matrix(as.matrix(Y) - Yhat)
  U2 <- c(U) * c(U)
  
  # Bootstrap coefficients matrix
  Boot_betas <- matrix(0, nrow=length(beta), ncol=n_boot)
  
  # Compute beta for each bootstrapped sample
  for (b in 1:n_boot) {
    # Bootstrapped data sample
    index_b <- sample(1:nrow(X), size=nrow(X), replace=TRUE)
    Y_b <- Y[index_b]
    X_b <- X[index_b]
    Z_b <- Z[index_b]
    
    # Construct jackknife instrument
    ZtZinv <- solve(t(Z_b) %*% as.matrix(Z_b), tol=tolerance)
    h_i <- diag(as.matrix(Z_b) %*% ZtZinv %*% t(Z_b))
    Pi_hat <- t(X_b) %*% as.matrix(Z_b) %*% ZtZinv
    X_jack <- (1-h_i)^(-1)*(as.matrix(Z_b) %*% t(Pi_hat)  - h_i * X_b)
    
    # Coefficients
    XtXinv_jack <- solve(t(X_jack) %*% as.matrix(X_b), tol=tolerance)
    XtY_jack <- t(X_jack) %*% as.matrix(Y_b)
    beta_b <- XtXinv_jack %*% XtY_jack
    
    # Store bootstrapped coefficient
    Boot_betas[,b] <- beta_b
  }
  
  # Mean of all bootstrapped coefficients
  beta_b_mean <- rowMeans(Boot_betas)
  
  # Subtract mean from coefficients, square, take mean, and then sqrt
  Boot_betas <- Boot_betas - beta_b_mean
  Boot_betas <- Boot_betas^2
  SEs <- sqrt(rowMeans(Boot_betas))
  
  # Return results
  results <- list(beta = t(beta), SEs = SEs)
  return(results)
}