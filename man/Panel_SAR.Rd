% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/panel_SAR.R
\name{Panel_SAR}
\alias{Panel_SAR}
\title{fitting panel spatial autoregressive model}
\usage{
Panel_SAR(y, X, W, n, T, ndraw = 2000, nomit = 100)
}
\arguments{
\item{y}{y is the response;}

\item{X}{X is data frame ccontains data;}

\item{W}{W is the spatial weight matrix;}

\item{n}{n is the observation number;}

\item{T}{T is the time period of data;}

\item{ndraw}{ndraw is the draw times for MCMC;}

\item{nomit}{noimt is the draw times omit for burn in;}
}
\value{
a dataframe with (ndraw - nomit) row, each row is a MCMC draw;
}
\description{
A function for estimate panel spatial autoregressive model
}
\examples{
#uncomment code to run test
#dt=sim_Panel_SAR(N=100,T=5,beta=c(0.5,-0.8),rho=-0.5,alpha=0.3);
#md=Panel_SAR(y=dt$Y,X=apply(dt$X,2,scale),W=dt$W,n=dt$N,T=dt$T,ndraw=1000)
#round(apply(md,2,mean),3);round(apply(md,2,mean)/apply(md,2,sd),3)
#result
#	    beta_1  beta_2     rho      mu   xi_sqr
#[1]   0.537   -0.825  -0.454   0.384    0.635
#[1]  13.673  -20.772  -1.853   3.528    5.342
}
