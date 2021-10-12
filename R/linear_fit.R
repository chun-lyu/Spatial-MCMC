#' @title linear model fitting function
#'
#' @description basic version of probit model fitting function.
#'
#' @param y y is the response;
#' @param X X is data frame ccontains data;
#' @param n n is the observation number;
#' @param ndraw ndraw is the draw times for MCMC;
#' @param nomit noimt is the draw times omit for burn in;
#'
#' @return a dataframe with (ndraw - nomit) row, each row is a MCMC draw;
#' @export
#' @examples #uncomment code to run test
#' @examples #dt = sim_linear(N=100,beta=c(0.2,-0.8));
#' @examples #md = linear_fit(y=dt$Y,X=dt$X,n = dt$N,ndraw=10000)
#' @examples #round(apply(md,2,mean),3);round(apply(md,2,mean)/apply(md,2,sd),3)
#' @examples #result
#' @examples #		       mean	   	sd     	true value
#' @examples #beta_1  0.203   1.806		         0.2
#' @examples #beta_2 -0.805  -7.702            -0.8

linear_fit = function(y,X,n,ndraw=2000,nomit=1000){
  # library(tmvtnorm)
  # library(MCMCpack)
  Sigma_0 = 0.001*diag(ncol(X))
  beta <- rep(0,ncol(X))
  xi_sqr <- 100
  B = matrix(0,nrow=(ndraw-nomit),ncol=ncol(X))
  pb <- txtProgressBar(style=3)

  for (i in 1:ndraw){
    ##draw beta
    Sigma_hat = solve(t(X)%*%X+Sigma_0)
    beta_hat = Sigma_hat %*%t(X)%*%y
    beta = t(rmvnorm(n = 1, mean = beta_hat, sigma = Sigma_hat))

    ##draw xi
    v_hat = 2+n
    lambda_hat = 0.01+crossprod(y - X%*%beta)
    xi_sqr = rinvgamma(n = 1, shape=v_hat/2, scale=lambda_hat/2)

    if (i > nomit) B[(i-nomit),]=c(beta)
    setTxtProgressBar(pb, i/ndraw)
  }
  return(B)
  close(pb)
}



sim_linear = function(N,beta){
  set.seed(1)
  X = matrix(rnorm(N*length(beta)),N,length(beta))
  Y = X%*%beta
  Y = Y+rnorm(N,0,0.1)
  return(list(Y=Y,X=X,N=N))
}
