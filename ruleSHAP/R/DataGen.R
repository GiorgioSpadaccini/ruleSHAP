#' Simulate data using the Friedman generating function
#'
#' This function creates synthetic data with the Friedman data-generating scheme:
#' \deqn{x_1,\ldots,x_{5} \sim \mathcal{U}(0,1),}
#' \deqn{\epsilon \sim \mathcal{N}(0,\sigma^2),}
#' \deqn{y=10\sin(\pi x_1x_2)+20(x_3-0.5)^2+10x_4+5x_5+\epsilon.}
#' This function may also binarize the outcome in the following manner:
#' \deqn{\eta=10\sin(\pi x_1x_2)+20(x_3-0.5)^2+10x_4+5x_5,}
#' \deqn{p=\frac{exp(\eta-\overline{\eta})}{exp(\eta-\overline{\eta})+1},}
#' \deqn{y \sim \mathcal{B}\textbf{e}(p),}
#' with \eqn{\overline{\eta} \approx 14.4} computed as the approximate mean of
#' \eqn{\eta} under this data-generation scheme, in the case of independent features.
#' Features may be generated to be correlated, using a gaussian copula.
#'
#'
#' @param n number of observations that should be generated.
#' @param p number of predictors that should be generated. Should be at least 5. For any
#' number higher than \eqn{p=5}, the remaining features are all noise features.
#' @param rho correlation induced across the features.
#' @param sd In the case of continuous outcome, the standard deviation of the normally-distributed
#' irreducible noise. For binary outcome, this parameter is ignored.
#' @param binary A logical parameter determining whether the outcome should binarized (\code{TRUE})
#' or should be left continuous (\code{FALSE}).
#'
#' @return A dataframe containing the synthetic dataset.
#'
#' @export
#'
gendata.friedman1=function(n,p=100,rho=NULL,sd=1,binary=F){
  if(is.null(rho)){
    rho=diag(p)
  }

  #use gaussian copula to have Corr[data$x]=rho
  #adjust correlations to account for copula
  corr_u=2*sin(pi/6*rho)
  u=mvtnorm::rmvnorm(n=n,mean=rep(0,p),sigma = corr_u)
  x=pnorm(u)


  y=10*sin(x[,1]*x[,2]*pi)+20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]

  #Convert outcome to binary, if requested
  if(binary){
    avg_eta=10*rmutil::int2(function(x,z){sin(x*z*pi)},a=c(0,0),b=c(1,1))+
      5/3+5+2.5
    y=rbinom(n=n,size=1,prob=1/(exp(-y+avg_eta)+1))
  }
  else{
    y=y+rnorm(n,sd=sd)
  }

  data <- data.frame(x = x, y = y)
  return(data=data)
}
