test_that("RuleSHAP model fits correctly", {
  #Define input arguments
  data=gendata.friedman1(n=100,p=5)
  formula=y~.
  ntree=3
  burn.in=2
  nmc=2

  #Fit the model
  fitted_model=expect_no_error(ruleSHAP(formula,data,ntree=3,burn.in=2,nmc=2,verbose = F))
})












#' Fit a horseshoe model
#'
#' This function fits a plain horseshoe Bayesian regression for continuous and binary data with
#' a MCMC algorithm based on Gibbs sampling. This function is also compatible with a grouped
#' local shrinkage for linear terms and rule terms separately.
#'
#' @param X a numeric \eqn{n \times p} matrix whose columns are predictors and whose rows are observations.
#' This is used as training data to fit the model.
#' @param y a numeric vector of outcomes used as training data to fit the model. For binary outcome, this should
#' be converted into zeroes and ones.
#' @param p.lin number of columns of X which are linear terms and should thus get their
#' own global shrinkage separately from the last remaining columns of X. By default, all
#' columns are linear terms and therefore all columns share the same shrinkage.
#' @param family A character string denoting the type of (Generalized) linear model
#' that is being fitted
#' @param method.tau.rules character string denoting the prior distribution of the
#' global shrinkage assigned to rules: the choice is between a half-Cauchy prior
#' (\code{"halfCauchy"}),a truncated Cauchy prior (\code{"truncatedCauchy"}) or
#' a fixed value (\code{"fixed"})
#' @param tau.rules a numeric entry denoting which value the global shrinkage for rules should take,
#' in the case where \code{method.tau.rules} is set to \code{"fixed"}. If other values of \code{method.tau.lin}
#' are chosen, this parameter is used as initial value to start the Gibbs sampling from.
#' @param method.tau.lin character string denoting the prior distribution of the
#' global shrinkage assigned to linear terms: the choice is between a half-Cauchy prior
#' (\code{"halfCauchy"}),a truncated Cauchy prior (\code{"truncatedCauchy"}) or
#' a fixed value (\code{"fixed"})
#' @param tau.rules a numeric entry denoting which value the global shrinkage for linear terms should take,
#' in the case where \code{method.tau.lin} is set to \code{"fixed"}. If other values of \code{method.tau.lin}
#' are chosen, this parameter is used as initial value to start the Gibbs sampling from.
#' @param method.sigma character string denoting which prior distribution the error variance \eqn{\sigma^2}
#' should have: the choice is between Jeffre's prior (\code{"Jeffreys"}) or a fixed
#' value (\code{"fixed"}).
#' @param Sigma2 a numeric entry denoting which value the error variance \eqn{\sigma^2} should take,
#' in the case where \code{method.sigma} is set to \code{"fixed"}. If other values of \code{method.sigma}
#' are chosen, this parameter is ignored.
#' @param prior_lambda a vector of \eqn{p} positive numeric entries representing local shrinkage relaxation
#' (\code{prior_lambda > 1}) or contraction (\code{prior_lambda < 1}). The \eqn{j}-th entry of this vector
#' denotes by how much the local shrinkage should be relaxed for the \eqn{j}-th predictor.
#' @param burn.in Number of burn-in MCMC iterations that are discarded before
#' the draws start to be recorded.
#' @param nmc Number of MCMC iterations after the initial burn-in period
#' @param thin Thinning parameter: if saving all MCMC draws is computationally
#' infeasible, thinning makes the model only record every \code{thin}-th iteration.
#' @param ping If the model is selected to be verbose, this denotes every how many
#' MCMC draws an update message should be printed
#' @param Beta_init a numeric vector of lenth \eqn{p} representing the initial coefficients to start
#' the Gibbs sampler with.
#' @param Lambda_init a numeric vector of lenth \eqn{p} representing the initial values of
#' local shrinkage to start the Gibbs sampler with.
#'
#' @return BetaSamples A matrix where every column is a MCMC sample of the coefficients.
#' @return BetaHat A vector containing the posterior mean of the coefficients.
#' @return TauLinSamples A vector where every entry is a MCMC sample of the global shrinkage of the linear terms.
#' @return TauLinHat The posterior mean of the global shrinkage of the linear terms.
#' @return TauRulesSamples A vector where every entry is a MCMC sample of the global shrinkage of the rules.
#' @return TauRulesHat The posterior mean of the global shrinkage of the rules.
#' @return LambdaSamples A matrix where every column is a MCMC sample of the local shrinkages, one entry per term.
#' @return LambdaHat A vector containing the posterior mean of the local shrinkages.
#' @export
hs=function(X, y, family=c('gaussian',"binomial"), p.lin=ncol(X),
            method.tau.rules = c("halfCauchy", "truncatedCauchy", "fixed"),
            tau.rules = 1, method.tau.lin = c("halfCauchy", "truncatedCauchy", "fixed"),
            tau.lin = 1, method.sigma = c("Jeffreys","fixed"), Sigma2 = 1,
            burn.in = 1000, nmc = 5000, thin = 1, prior_lambda = NULL,
            Beta_init=rep(0,ncol(X)),Lambda_init=rep(1,ncol(X)),ping=1e3){

  #Retrieve arguments and basic quantities
  family=match.arg(family)
  method.tau.lin = match.arg(method.tau.lin)
  method.tau.rules = match.arg(method.tau.rules)
  method.sigma = match.arg(method.sigma)
  niter <- burn.in + nmc
  effsamp = (niter - burn.in)/thin

  n = nrow(X)
  p <- ncol(X)
  Beta <- Beta_init
  lambda <- Lambda_init
  Xt=t(X)
  omega=sqrt.omega=rep(1,n)
  betaout <- matrix(0, p, effsamp)
  lambdaout <- matrix(0, p, effsamp)
  taulinout <- tauruleout <- rep(0, effsamp)

  algoBhatt = p>=n

  #If prior for local shrinkage is not specified, it's just uniform
  if(is.null(prior_lambda)){
    prior_lambda = rep(1, times=p)
  }

  #If logistic, re-adjust y and sigma as per Polya Inverse Gamma augmentation
  if(family=='binomial'){
    y <- y - 0.5
    method.sigma='fixed'
    Sigma2=1
  }

  #If family is gaussian and algorithm is Rue, Q_star must only be computed once
  if(family=='gaussian' & algoBhatt == F){
    Q_star=crossprod(X)
  }


  message("Markov chain monte carlo is running")
  for (i in 1:niter) {
    #Update lambda directly, without slice sampling
    #As in horseshoe, impose values not too large lambda from previous iteration
    #Determine threshold for eta
    #tempps = (Beta/(tau*prior_lambda))^2/2/Sigma2
    #thresh = -log(1-1e-4)/tempps
    #Now sample, imposing that 1/lambda^2 larger than threshold
    #nu.inv=stats::rexp(p, 1+pmax(1/lambda^2,thresh))
    #lambda=sqrt(1/stats::rexp(p, nu.inv + tempps))

    #Combine the vector of tau.rules for rules and tau.lin for linear terms
    tau_vec=c(rep(tau.lin,p.lin),rep(tau.rules,p-p.lin))

    #As in horserule
    tempps = (Beta/(tau_vec*prior_lambda))^2/2/Sigma2
    nu.inv=stats::rexp(p, 1+1/lambda^2)
    lambda=sqrt(1/stats::rexp(p, nu.inv + tempps))

    #Update lambda (as in horseshoe)
    #eta = 1/(lambda^2)
    #upsi = stats::runif(p, 0, 1/(1 + eta))
    #tempps = (Beta/(tau*prior_lambda))^2/2/Sigma2
    #ub = (1 - upsi)/upsi
    #Fub = 1 - exp(-tempps * ub)
    #Fub[Fub < (1e-04)] = 1e-04
    #up = stats::runif(p, 0, Fub)
    #eta <- -log(1 - up)/tempps
    #lambda = 1/sqrt(eta)




    #Update tau.lin
    if (method.tau.lin == "halfCauchy") {
      #directly, without slice sampling
      #As in horseshoe, impose values not too large lambda from previous iteration
      #Determine threshold for xi
      #tempt = sum(Beta^2/(2 * (lambda*prior_lambda)^2 * Sigma2))
      #thresh = stats::qgamma(1e-8, (p + 1)/2, scale = 1/tempt)
      #Now sample, imposing that 1/lambda^2 larger than threshold
      #xi.inv = rexp(1, 1 + max(1/tau^2,thresh))
      #tau = sqrt(1/rgamma(1, shape = (p + 1)/2,rate = xi.inv + tempt))

      #As in horserule, no threshold
      xi = 1/rexp(1, 1 + 1/tau.lin^2)
      tau.lin = sqrt(1/rgamma(1, shape = (p.lin + 1)/2,
                              rate = 1/xi + sum(Beta[1:p.lin]^2/(2 * (lambda[1:p.lin]*prior_lambda[1:p.lin])^2 * Sigma2))))

      #as in horseshoe
      #tempt <- sum(Beta^2/(2 * (lambda*prior_lambda)^2 * Sigma2))
      #et = 1/tau^2
      #utau = stats::runif(1, 0, 1/(1 + et))
      #ubt = (1 - utau)/utau
      #Fubt = stats::pgamma(ubt, (p + 1)/2, scale = 1/tempt)
      #Fubt = max(Fubt, 1e-08)
      #ut = stats::runif(1, 0, Fubt)
      #et = stats::qgamma(ut, (p + 1)/2, scale = 1/tempt)
      #tau = 1/sqrt(et)
    }
    if (method.tau.lin == "truncatedCauchy") {
      #as in horseshoe
      tempt <- sum(Beta^2/(2 * (lambda[1:p.lin]*prior_lambda[1:p.lin])^2 * Sigma2))
      et = 1/tau.lin^2
      utau = stats::runif(1, 0, 1/(1 + et))
      ubt_1 = 1
      ubt_2 = min((1 - utau)/utau, p.lin^2)
      Fubt_1 = stats::pgamma(ubt_1, (p.lin + 1)/2, scale = 1/tempt)
      Fubt_2 = stats::pgamma(ubt_2, (p.lin + 1)/2, scale = 1/tempt)
      ut = stats::runif(1, Fubt_1, Fubt_2)
      et = stats::qgamma(ut, (p.lin + 1)/2, scale = 1/tempt)
      tau.lin = 1/sqrt(et)
    }



    #Update tau.rules
    if (method.tau.rules == "halfCauchy") {
      #directly, without slice sampling
      #As in horseshoe, impose values not too large lambda from previous iteration
      #Determine threshold for xi
      #tempt = sum(Beta^2/(2 * (lambda*prior_lambda)^2 * Sigma2))
      #thresh = stats::qgamma(1e-8, (p + 1)/2, scale = 1/tempt)
      #Now sample, imposing that 1/lambda^2 larger than threshold
      #xi.inv = rexp(1, 1 + max(1/tau^2,thresh))
      #tau = sqrt(1/rgamma(1, shape = (p + 1)/2,rate = xi.inv + tempt))

      #As in horserule, no threshold
      xi = 1/rexp(1, 1 + 1/tau.rules^2)
      tau.rules = sqrt(1/rgamma(1, shape = (p-p.lin + 1)/2,
                                rate = 1/xi + sum(Beta[-(1:p.lin)]^2/(2 * (lambda[-(1:p.lin)]*prior_lambda[-(1:p.lin)])^2 * Sigma2))))

      #as in horseshoe
      #tempt <- sum(Beta^2/(2 * (lambda*prior_lambda)^2 * Sigma2))
      #et = 1/tau^2
      #utau = stats::runif(1, 0, 1/(1 + et))
      #ubt = (1 - utau)/utau
      #Fubt = stats::pgamma(ubt, (p + 1)/2, scale = 1/tempt)
      #Fubt = max(Fubt, 1e-08)
      #ut = stats::runif(1, 0, Fubt)
      #et = stats::qgamma(ut, (p + 1)/2, scale = 1/tempt)
      #tau = 1/sqrt(et)
    }
    if (method.tau.rules == "truncatedCauchy") {
      #as in horseshoe
      tempt <- sum(Beta^2/(2 * (lambda[-(1:p.lin)]*prior_lambda[-(1:p.lin)])^2 * Sigma2))
      et = 1/tau.rules^2
      utau = stats::runif(1, 0, 1/(1 + et))
      ubt_1 = 1
      ubt_2 = min((1 - utau)/utau, (p-p.lin)^2)
      Fubt_1 = stats::pgamma(ubt_1, (p-p.lin + 1)/2, scale = 1/tempt)
      Fubt_2 = stats::pgamma(ubt_2, (p-p.lin + 1)/2, scale = 1/tempt)
      ut = stats::runif(1, Fubt_1, Fubt_2)
      et = stats::qgamma(ut, (p-p.lin + 1)/2, scale = 1/tempt)
      tau.rules = 1/sqrt(et)
    }


    #Update sigma
    if (method.sigma == "Jeffreys") {
      if (algoBhatt) {
        E_1 = max(t(y - X %*% Beta) %*% (y - X %*% Beta),
                  (1e-10))
        E_2 = max(sum(Beta^2/((tau_vec * prior_lambda * lambda))^2), 1e-10)
      }
      else {
        E_1 = max(t(y - X %*% Beta) %*% (y - X %*% Beta),
                  1e-08)
        E_2 = max(sum(Beta^2/((tau_vec * prior_lambda * lambda))^2), 1e-08)
      }
      Sigma2 = 1/stats::rgamma(1, (n + p)/2, scale = 2/(E_1 +
                                                          E_2))
    }





    #Sample omega if bernoulli family (otherwise it stays 1 so no damage done)
    #In doing so, update z
    if(family=='binomial'){
      omega=pgdraw::pgdraw(1, X%*%Beta)
      sqrt.omega <- sqrt(omega)
      Xt=t(sqrt.omega * X)
      z=y/omega
    } else{
      z=y/sqrt(Sigma2)
    }

    #Update beta
    lambda_star = tau_vec * lambda * prior_lambda

    if (algoBhatt) {
      U = lambda_star^2 * Xt #this is D * Phi^t
      u = lambda_star*stats::rnorm(p)
      v = sqrt.omega * (X %*% u) + stats::rnorm(n)
      #Use left mult. by sqrt(Omega) to re-write the following in a more efficient:
      #v_star = solve((sqrt.omega*X %*% U + I_n), z*sqrt.omega - v)
      v_star = solve((X %*% U + diag(1/sqrt.omega)), z - v/sqrt.omega)
      Beta = sqrt(Sigma2) * c(u + U %*% v_star)
    }
    else{
      if(family=='binomial'){
        Q_star=tcrossprod(Xt)
      }
      L = chol((Q_star + diag(1/lambda_star^2,p,p))/Sigma2)
      v = solve(t(L), Xt %*% (z*sqrt.omega/sqrt(Sigma2)))
      mu = solve(L, v)
      u = solve(L, stats::rnorm(p))
      Beta = c(mu + u)
    }





    #Print update if requested
    if (ping > 0 & i%%ping == 0) {
      message("iteration = ", i)
    }
    if (i > burn.in && i%%thin == 0) {
      betaout[, (i - burn.in)/thin] <- Beta
      lambdaout[, (i - burn.in)/thin] <- lambda
      taulinout[(i - burn.in)/thin] <- tau.lin
      tauruleout[(i - burn.in)/thin] <- tau.rules
    }
  }
  pLambda <- apply(lambdaout, 1, mean)
  result = list(BetaSamples = betaout,
                BetaHat= rowMeans(betaout),
                LambdaSamples = lambdaout,
                LambdaHat = rowMeans(lambdaout),
                TauLinSamples = taulinout,
                TauLinHat= mean(taulinout),
                TauRulesSamples = tauruleout,
                TauRulesHat= mean(tauruleout))
  return(result)
}
