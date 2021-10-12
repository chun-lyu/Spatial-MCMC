sim_Panel_SAR = function(N,T,beta,rho,alpha){
  library(nimble)
  set.seed(1)
  alpha = rnorm(N,alpha,1);
  epsilon = rnorm(N*T,0,0.2)
  X_1 = rnorm(N*T);X_2 = rnorm(N*T);X = cbind(X_1,X_2)
  W = matrix(sample(c(0,0,0,0,0,0,0,0,1),N*N,replace=TRUE),N,N)
  W = t(W)%*%W
  diag(W) = 0
  W = as.matrix(t(apply(W,1,function(x)x/sum(x))))
  alpha = kronecker(rep(1,T),diag(N))%*%as.matrix(alpha,N,1)
  S_1 = solve(diag(N*T) - rho * kronecker(diag(T),W))
  Y = S_1%*%(X%*%beta+alpha+epsilon)
  return(list(Y=Y,X=cbind(X_1 = X_1,X_2=X_2),W=W,N=N,T=T))
}
