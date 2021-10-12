#' @title fitting probit model
#'
#' @description A function for estimate probit model
#'
#' @param y y is binary response;
#' @param X X is data frame ccontains data;
#' @param n n is the observation number;
#' @param ndraw ndraw is the draw times for MCMC;
#' @param nomit noimt is the draw times omit for burn in;
#'
#' @return a dataframe with (ndraw - nomit) row, each row is a MCMC draw;
#' @export
#' @examples #uncomment code to run test
#' @examples #library(tmvtnorm)
#' @examples #dt=sim_probit(N=100,beta=c(0.5,-0.7));
#' @examples #md=probit(y=dt$Y,X=dt$X,n=dt$N,ndraw=3000)
#' @examples #round(apply(md,2,mean),3);round(apply(md,2,mean)/apply(md,2,sd),3)
#' @examples #  result
#' @examples #  [1]  0.394 -0.685
#' @examples #  [1]  2.430 -4.183


probit = function(y,X,n,ndraw=2000,nomit=1000){
  set.seed(1)
  # library(tmvtnorm)
  # library(MCMCpack)
  # library(Matrix)
  Sigma_hat = solve(t(X)%*%X+0.001*diag(ncol(X)))
  H = as(diag(n),"dgCMatrix")
  beta <- rep(0,ncol(X))
  xi_sqr <- 1
  B = matrix(0,nrow=(ndraw-nomit),ncol=ncol(X))
  lower = ifelse(y > 0, 0, -Inf)
  upper = ifelse(y > 0, Inf, 0)
  pb <- txtProgressBar(style=3)     #开始进度条

  for (i in 1:ndraw){
    ##draw y_star
    mu_1 = X%*% beta
    y_star = t(rtmvnorm(n = 1, mean = as.numeric(mu_1), H = H, lower = lower,
                        upper = upper,algorithm="gibbs",burn.in.samples=50))
    if(any(is.na(y_star)))stop("NA in y_star")

    ##draw beta
    beta_hat = Sigma_hat %*%t(X)%*%y_star
    beta = t(rmvnorm(n = 1, mean = beta_hat, sigma = Sigma_hat*xi_sqr))

    ##draw xi_sqr
    xi_sqr = rinvgamma(n = 1, shape=1+n/2, scale=0.5+crossprod(y_star-X%*%beta)/2)

    if (i > nomit) B[(i-nomit),]=c(beta)
    setTxtProgressBar(pb, i/ndraw)
  }
  return(B)
  close(pb)	#结束进度条
}


sim_probit = function(N,beta){
  set.seed(1)
  X = matrix(rnorm(N*length(beta)),N,length(beta))
  Y = X%*%beta
  Y = Y+rnorm(N,0,0.2)
  for (i in 1:N) Y[i] = ifelse(pnorm(Y[i])>runif(1),1,0)
  return(list(Y=Y,X=X,N=N))
}
