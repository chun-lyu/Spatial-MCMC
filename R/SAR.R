#' @title fitting spatial autoregressive model
#'
#' @description A function for estimate spatial autoregressive model
#'
#' @param y y is the response;
#' @param X X is data frame ccontains data;
#' @param W W is the spatial weight matrix;
#' @param n n is the observation number;
#' @param ndraw ndraw is the draw times for MCMC;
#' @param nomit noimt is the draw times omit for burn in;
#'
#' @return a dataframe with (ndraw - nomit) row, each row is a MCMC draw;
#' @export
#' @examples #uncomment code to run test
#' @examples #dt = sim_SAR(N=100,beta=c(0.5,-0.8),rho=-0.5);
#' @examples #md = SAR(y=dt$Y,X=apply(dt$X,2,scale),W=dt$W,n=dt$N,ndraw=2000)
#' @examples #round(apply(md,2,mean),3);round(apply(md,2,mean)/apply(md,2,sd),3)
#' @examples #Result
#' @examples #[1]   0.479   -0.825  -0.402
#' @examples #[1]  22.391  -58.843  -1.800

SAR = function(y,X,W,n,ndraw=2000,nomit=100){
  # library(tmvtnorm)
  # library(MCMCpack)
  # library(Matrix)
  rho_prob = function(rho,W,y,X,beta,xi_sqr){
    T = nrow(X)/nrow(W)
  #  I_n = as(diag(nrow(X)),"dgCMatrix")
    I_n = diag(nrow(X))
    I_n_W = I_n- rho * W
    e = I_n_W%*%y - X%*%beta
    prob = det(I_n_W)*exp(-0.5*crossprod(e)/xi_sqr)#
    return(as.numeric(prob))
  }

  set.seed(1)
  c = 0.5
  Sigma_hat = solve(t(X)%*%X+0.001*diag(ncol(X)))
  beta = rep(0,ncol(X))
  xi_sqr = 1
  acc_rt = 0
  rho = 0.5;

  I_n = diag(n)
#  I_n = as(diag(n),"dgCMatrix")
 # W = as(W,"dgCMatrix")
  B = matrix(0,nrow=(ndraw-nomit),ncol=(ncol(X)+1))
  pb <- txtProgressBar(style=3)

  for (i in 1:ndraw){
    ##draw rho
    rho_new = rho + rnorm(1)*c
    rho_new_prb = rho_prob(rho_new,W,y,X,beta,xi_sqr)
    rho_old_prb = rho_prob(rho,W,y,X,beta,xi_sqr)
    prob = min(1,rho_new_prb/rho_old_prb)
    if (prob>runif(1) & rho_new > -1 & rho_new < 1){
      rho = rho_new;
      acc_rt = acc_rt+1;
      }
    if ((acc_rt/i)>.6){
      c=1.1*c
    }else if((acc_rt/i)<.4){
        c=c/1.1
    }else{
        c=c
    }
    S = I_n - rho * W

    ##draw beta
    beta_hat = Sigma_hat%*%t(X)%*%S%*%y
    beta = t(rmvnorm(n = 1, mean = beta_hat, sigma = Sigma_hat*xi_sqr))

    ##draw xi_sqr
    xi_sqr = rinvgamma(n = 1, shape=1+n/2, scale=0.5+crossprod(S%*%y-X%*%beta)/2)

    if (i > nomit)B[(i-nomit),]=c(beta,rho)
    setTxtProgressBar(pb, i/ndraw)
  }
  close(pb)
  return(B)
}




sim_SAR = function(N,beta,rho){
  set.seed(1)
  epsilon = rnorm(N,0,0.1)
  X_1 = rnorm(N*T);
  X_2 = rnorm(N*T);X = cbind(X_1,X_2)
  W = matrix(sample(c(0,0,0,0,0,0,1),N*N,replace=TRUE),N,N)
  W = t(W)%*%W
  diag(W) = 0
  W = as.matrix(t(apply(W,1,function(x)x/sum(x))))
  S_1 = solve(diag(N) - rho * W)
  Y = S_1%*%(X%*%beta+epsilon)
  return(list(Y=Y,X=cbind(X_1 = X_1,X_2=X_2),W=W,N=N))
}

