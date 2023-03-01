# Last updated: Feb 28, 2023

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


#### post lasso ----
# see ipynb + section in analysis using hdm

