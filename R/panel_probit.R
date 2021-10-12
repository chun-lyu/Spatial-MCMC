#' @title fitting panel probit model
#'
#' @description A function for estimate panel probit model
#'
#' @param y y is binary response;
#' @param X X is data frame ccontains data;
#' @param n n is the observation number;
#' @param T T is the time period of data;
#' @param ndraw ndraw is the draw times for MCMC;
#' @param nomit noimt is the draw times omit for burn in;
#'
#' @return a dataframe with (ndraw - nomit) row, each row is a MCMC draw;
#' @export
#' @examples #uncomment code to run test
#' @examples #dt=sim_panel_probit(N=100,T=10,beta=c(0.7,-0.5,0.3),alpha=0.5);
#' @examples #md=panel_probit(y=dt$Y,X=dt$X,n=dt$N,T=dt$T,ndraw=1000)
#' @examples #round(apply(md,2,mean),3);round(apply(md,2,mean)/apply(md,2,sd),3)
#' @examples #result
#' @examples # 	    beta_1   beta_2  beta_3     mu  xi_sqr
#' @examples #[1]    0.724   -0.573   0.308  0.593   1.142
#' @examples #[1]   11.003   -9.525   5.246  4.824   4.195

panel_probit = function(y,X,n,T,ndraw=2000,nomit=100){
  # library(tmvtnorm)
  # library(MCMCpack)
  # library(Matrix)
  # library(MCMCglmm)
  set.seed(1)
  Sigma_hat = solve(t(X)%*%X+0.001*diag(ncol(X)))
  beta <- rep(0,ncol(X))
  mu <- 0.5
  xi_sqr <- 100
  alpha <- y_star <- rep(0,n)
  lower = ifelse(y > 0, 0, -Inf)
  upper = ifelse(y > 0, Inf, 0)
  B = matrix(0,nrow=(ndraw-nomit),ncol=(ncol(X)+2))
  pb <- txtProgressBar(style=3)

  for (i in 1:ndraw){
    ##draw y_star
    mu_1 = X%*%beta+kronecker(matrix(rep(1,T),T,1),alpha)
    for(k in 1:(n*T))y_star[k]<-rtnorm(1,mean=mu_1[k],sd=1,lower=lower[k],upper=upper[k])
    if(any(is.na(y_star)))stop("NA in y_star")

    ##draw beta
    beta_mu = t(X)%*%(y_star-kronecker(matrix(rep(1,T),T,1),diag(n))%*%alpha)
    beta_hat = Sigma_hat %*%beta_mu
    beta = t(rmvnorm(n = 1, mean = beta_hat, sigma = Sigma_hat))

    ##draw alpha
    xi_sqr_hat = 1/(T+1/xi_sqr)
    alpha_1 = mu/xi_sqr
    alpha_y=apply(matrix((y_star-X%*%beta),n,T),1,sum)
    for (j in 1:n){
      alpha_hat = xi_sqr_hat*(alpha_y[j]+alpha_1)
      alpha[j] = rnorm(1,mean=alpha_hat,sd=sqrt(xi_sqr_hat))
    }

    ##draw mu and xi
    mu_sigma = 1/(n/xi_sqr+0.01)
    mu_hat = mu_sigma*(sum(alpha)/xi_sqr+0.01)
    mu = rnorm(1,mean=mu_hat,sd=sqrt(mu_sigma))

    v_hat = 2+n
    lambda_hat = 1+crossprod(alpha-mu)
    xi_sqr = rinvgamma(n = 1, shape=v_hat/2, scale=lambda_hat/2)

    if (i > nomit)B[(i-nomit),]=c(beta,mu,xi_sqr)
    setTxtProgressBar(pb, i/ndraw)
  }
  close(pb)
  return(B)
}


sim_panel_probit = function(N,T,beta,rho,alpha){
  library(nimble)
  set.seed(1)
  alpha = rnorm(N,alpha,1);
  epsilon = rnorm(N*T,0,0.3)
  X = matrix(rnorm(N*T*length(beta)),N*T,length(beta))
  alpha = kronecker(rep(1,T),diag(N))%*%as.matrix(alpha,N,1)
  Y = X%*%beta+alpha+epsilon
  for (i in 1:(N*T)) Y[i] = ifelse(pnorm(Y[i])>runif(1),1,0)
  return(list(Y=as.numeric(Y),X=X,N=N,T=T))
}

