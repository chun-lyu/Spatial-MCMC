#' @title fitting panel spatial autoregressive model
#'
#' @description A function for estimate panel spatial autoregressive model
#'
#' @param y y is the response;
#' @param X X is data frame ccontains data;
#' @param W W is the spatial weight matrix;
#' @param T T is the time period of data;
#' @param n n is the observation number;
#' @param ndraw ndraw is the draw times for MCMC;
#' @param nomit noimt is the draw times omit for burn in;
#'
#' @return a dataframe with (ndraw - nomit) row, each row is a MCMC draw;
#' @export
#' @examples #uncomment code to run test
#' @examples #dt=sim_Panel_SAR(N=100,T=5,beta=c(0.5,-0.8),rho=-0.5,alpha=0.3);
#' @examples #md=Panel_SAR(y=dt$Y,X=apply(dt$X,2,scale),W=dt$W,n=dt$N,T=dt$T,ndraw=1000)
#' @examples #round(apply(md,2,mean),3);round(apply(md,2,mean)/apply(md,2,sd),3)
#' @examples #result
#' @examples #	    beta_1  beta_2     rho      mu   xi_sqr
#' @examples #[1]   0.537   -0.825  -0.454   0.384    0.635
#' @examples #[1]  13.673  -20.772  -1.853   3.528    5.342

Panel_SAR = function(y,X,W,n,T,ndraw=2000,nomit=100){

  # library(tmvtnorm)
  # library(MCMCpack)
  # library(Matrix)
  rho_prob = function(rho,W,y_star,X,beta,alpha,xi_sqr){
    T = nrow(X)/nrow(W)
    I_n = diag(nrow(X))
    I_n_W = I_n- kronecker(diag(rep(rho,T)), W)
    e = I_n_W%*%y_star - X%*%beta - kronecker(matrix(rep(1,T),T,1),alpha)
    prob = det(I_n_W)*exp(-0.5*crossprod(e)/xi_sqr)#
    return(as.numeric(prob))
  }

  set.seed(1)
  c = 0.2
  Sigma_hat = solve(t(X)%*%X+0.001*diag(ncol(X)))
  beta <- rep(0,ncol(X))
  acc_rt <- 0
  mu <- 0.5
  xi_sqr <- 100
  alpha = rep(0,n)
  rho = 0.5;

  I_n = diag(n*T)
#  W = as(W,"dgCMatrix")
  S_1 = I_n - rho * kronecker(diag(T),W)
  B = matrix(0,nrow=(ndraw-nomit),ncol=(ncol(X)+3))
  pb <- txtProgressBar(style=3)

  for (i in 1:ndraw){
    ##draw rho
    rho_new = rho + rnorm(1)*c
    rho_new_prb = rho_prob(rho_new,W,y,X,beta,alpha,xi_sqr)
    rho_old_prb = rho_prob(rho,W,y,X,beta,alpha,xi_sqr)
    #		cat("rho_new_prb = ",rho_new_prb,"rho_old_prb = ",rho_old_prb,"\n")
    prob = min(1,rho_new_prb/rho_old_prb)##rho_new_prb>rho_old_prb
    if (prob>runif(1) & rho_new > -1 & rho_new < 1){rho = rho_new;acc_rt = acc_rt+1;}
    if ((acc_rt/i)>.6){c=1.1*c}else if((acc_rt/i)<.4){c=c/1.1}else{c=c}
    S_1 = I_n - rho * kronecker(diag(T),W)

    ##draw beta
    Sigma_1 = t(X)%*%(S_1%*%y-kronecker(matrix(rep(1,T),T,1),diag(n))%*%alpha)
    beta_hat = Sigma_hat %*%Sigma_1
    beta = t(rmvnorm(n = 1, mean = beta_hat, sigma = Sigma_hat*xi_sqr))

    ##draw alpha
    xi_sqr_hat = 1/(T+1/xi_sqr)
    alpha_1 = mu/xi_sqr
    alpha_y=apply(matrix((S_1%*%y-X%*%beta),n,T),1,sum)
    for (j in 1:n){
      alpha_hat = xi_sqr_hat*(alpha_y[j]+alpha_1)
      alpha[j] = rnorm(1,mean=alpha_hat,sd=sqrt(xi_sqr_hat))
    }

    ##draw xi
    mu_sigma = 1/(n/xi_sqr+0.01)
    mu_hat = mu_sigma*(sum(alpha)/xi_sqr+0.01)
    mu = rnorm(1,mean=mu_hat,sd=sqrt(mu_sigma))

    v_hat = 2+n
    lambda_hat = 1+crossprod(alpha-mu)
    xi_sqr = rinvgamma(n = 1, shape=v_hat/2, scale=lambda_hat/2)

    if (i > nomit)B[(i-nomit),]=c(beta,rho,mu,xi_sqr)
    setTxtProgressBar(pb, i/ndraw)
  }
  close(pb)
  return(B)
}




