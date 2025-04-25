#' Fit a RuleSHAP model
#'
#' This function fits a RuleSHAP model for continuous and binary data.
#' It is based upon fitting the following tweaked horseshoe prior:
#' \deqn{y|X,a,b,\sigma^2 \sim \mathcal{N}(\beta_0+\sum_{j=1}^pb_jx_j+\sum_{k=1}^qa_kr_k,\sigma^2I),}
#' \deqn{b_j|\lambda^{(l)}_j,\tau,\sigma^2 \sim \mathcal{N}\Big(0,\big(\lambda^{(l)}_j\tau\tau_0\sigma\big)^2\Big),}
#' \deqn{a_k|\lambda^{(r)}_k,\tau,\sigma^2 \sim \mathcal{N}\Big(0,\big(\lambda^{(r)}_k\tau\tau_1\sigma\big)^2),}
#' \deqn{\lambda^{(l)}_j \sim \mathcal{C}^+(0,A^{(l)}),}
#' \deqn{\lambda^{(r)}_j \sim \mathcal{C}^+(0,A^{(r)}_j),}
#' \deqn{A^{(r)}_j=\frac{\big(2\cdot \min(1-\overline{r}_j,\overline{r}_j)\big)^\mu}{(m_j)^\eta},}
#' \deqn{\tau_0,\tau_1,\tau \sim \mathcal{C}^+(0,1),}
#' \deqn{\sigma^2 \sim \sigma^{-2}d\sigma^2,}
#' where \eqn{m_j} is the number of predictors involved in the rule \eqn{r_j}.
#'
#' @param formula an object of class "\link{formula}" which represents the
#' desired predicted outcome (Left-hand side) and the desired predictors
#' (Right-hand side). If intercept should not be included, do not actively
#' remove it from the formula. Instead, use the parameter \code{intercept}
#' @param data a data frame containing the data to be used to fit the model
#' @param intercept a logical entry determining whether an intercept should be
#' used (\code{TRUE}) or not used (\code{FALSE}). The former is recommended.
#' @param family A character string denoting the type of (Generalized) linear model
#' that is being fitted.
#' @param mu,eta hyparameters defining how strong the prior shrinkage of complex
#' rules should be. A higher value of \code{mu} induces stronger shrinkage on
#' rules with extreme support, a higher value of \code{eta} induces a stronger
#' shrinkage on the rules involving many predictors.
#' @param lin.relax the level of shrinkage of linear terms and intercept.
#' This is equal to the hyperparameter \eqn{A^{(l)}}.
#' @param method.tau.lin the type of prior distribution for the global shrinkage
#' \eqn{\tau_0} which is shared by the coefficients of the linear terms. Set to
#' \code{"halfCauchy"} by default, but may be set to \code{"fixed"} to have it
#' fixed to a value equal to \eqn{A^{(l)}} (see parameter \code{lin.relax}).
#' @param ntree number of trees that should be used when generating the rules
#' @param maxdepth maximum depth of the trees used to generate the rules. In particular,
#' this induces a maximum on the order of interactions between predictors.
#' @param disaggregate a logical value determining whether rules should be extracted
#' by also considering all possible subrules (\code{TRUE}) or if rules should only
#' be extracted in a top-down manner as in RuleFit and HorseRule (\code{FALSE}).
#' @param rounding.digits Number of decimals that the splitting points
#' are rounded to in the rule generation. This rounding simplifies readability
#' of rules and prevents rules from being very similar but not identical.
#' @param burn.in Number of burn-in MCMC iterations that are discarded before
#' the draws start to be recorded.
#' @param nmc Number of MCMC iterations after the initial burn-in period
#' @param thin Thinning parameter: if saving all MCMC draws is computationally
#' infeasible, thinning makes the model only record every \code{thin}-th iteration.
#' @param ping If the model is selected to be verbose, this denotes every how many
#' MCMC draws an update message should be printed
#' @param verbose logical value denoting whether the model should send printed
#' updates on how far the fitting procedure is
#' @param rules Either \code{NULL} or a character vector containing the rules that
#' should be used to fit the RuleSHAP model, if one wants to specify them in advance.
#' If left to its default \code{NULL}, rules will be generated using a
#' Parametric Random Forest.
#' @param wins.frac A numeric entry representing the fraction by which the dataset
#' is winsorized. It also applies to rules by inducing a minimum node size.
#'
#' @return rules a character vector containing the rules that were generated and used
#' for the model fitting
#' @return BetaSamples A matrix where every column is a MCMC sample of the coefficients.
#' The first row represents the intercept. The following \eqn{p} represent
#' the linear terms. The remaining rows are for the rules, in the same order as in
#' output \code{rules}
#' @return BetaHat A vector containing the posterior mean of the coefficients.
#' The first entry represents the intercept. The following \eqn{p} represent
#' the linear terms. The remaining entries are for the rules, in the same order as in
#' output \code{rules}
#' @export
ruleSHAP=function(formula,data,intercept=T,family=c('gaussian','binomial'),
                  mu=0.5,eta=2,lin.relax=1,method.tau.lin="halfCauchy",
                  ntree=500,maxdepth=3,disaggregate=T,rounding.digits=2,
                  burn.in=2e3,nmc=2e4,thin=1,ping=1e3,
                  verbose=T,rules=NULL,wins.frac=0.025){


  #Retrieve basic quantities
  family=match.arg(family)
  y_name=as.character(formula[[2]])
  X_names=labels(terms(formula, data=data))


  #Winsorize data and store winsorization points
  wins.points=matrix(nrow=2,ncol=ncol(data))
  which_numeric=which(sapply(data,is.numeric) & (names(data) %in% X_names))
  for(i in which_numeric){
    wins.points[,i]=quantile(data[,i],probs=c(wins.frac,1-wins.frac))
    data[data[,i]<wins.points[1,i],i]=wins.points[1,i]
    data[data[,i]>wins.points[2,i],i]=wins.points[2,i]
  }


  #Subset to relevant variables only
  y=data[, y_name]
  X=data[, X_names]
  data=data[,c(y_name,X_names)]


  #Retrieve linear terms (in case some variables are factors) and re-scale them
  X_lin=model.matrix(formula,data=data)
  #Remove the intercept if requested
  if(!intercept){X_lin=X_lin[,-1]}
  #Standardize X_lin (but store info before doing so)
  muX_lin=apply(X_lin,MARGIN=2,FUN=mean)
  sdX_lin=apply(X_lin,MARGIN=2,FUN=sd)
  #Keep into account intercept in scaling
  if(intercept){
    X_lin[,-1]=scale(X_lin[,-1])
  } else{
    X_lin=scale(X_lin)
  }

  #Now actually for categorical entries do contrast-coding standardized as if balanced
  #Obtain number of categories and dummies used
  nlevs=unlist(lapply(X,nlevels))
  if(intercept){nlevs=c(0,nlevs)}
  ndummies=pmax(1,nlevs-1)
  #Compute which columns each variable uses
  end_i=cumsum(pmax(ndummies,1))
  start_i=c(1,end_i[-length(end_i)]+1)
  #For each categorical predictor, overwrite previous stdzation
  for(i in 1:length(nlevs)){
    if(nlevs[i]==0){ #Skip continuous columns (including intercept)
      next
    }
    else{ #Re-scale otherwise (use X_lin on RHS instead of X_mat)
      X_lin[,start_i[i]:end_i[i]]=(nlevs[i]*X_lin[,start_i[i]:end_i[i]]-1)/sqrt(nlevs[i]-1)
      #Update muX and sdX_lin acccordingly
      muX_lin[i]=1/nlevs[i]
      sdX_lin[i]=sqrt(nlevs[i]-1)/nlevs[i]
    }
  }

  #Get p (notice that p potentially includes the intercept) and n
  p=ncol(X_lin)
  n=nrow(data)


  #If gaussian, standardize Y
  if(family=='gaussian'){
    MuY=mean(y)
    SDy=sd(y)
    if(intercept){
      y=(y-MuY)/SDy
    }
  }


  #Fit parametric random forest on residuals if rules not pre-specified
  if(is.null(rules)){
    #Use the linear terms to compute residuals
    if(verbose){print('Generating rules')}
    HS_fit=hs(y=y, X=X_lin,family=family,
              method.tau.lin=method.tau.lin, tau.lin=1,
              method.tau.rules = 'fixed', tau.rules = 1,
              method.tau = 'fixed', tau = 1,
              burn=burn.in, nmc=nmc,
              method.sigma="Jeffreys")

    #Round numeric predictors and create the dataset for rule generation
    which_numeric=sapply(X,is.numeric)
    rulegen_data=X
    rulegen_data[,which_numeric]=round(X[,which_numeric],digits=rounding.digits)
    if(family=='gaussian'){
      rulegen_data$res=y-c(X_lin%*%HS_fit$BetaHat)
    } else if(family=='binomial'){
      #For now I'll define the working residuals
      phat=1/(1+exp(-c(X_lin%*%HS_fit$BetaHat)))
      rulegen_data$res=(y-phat)/(phat*(1-phat))
    }


    rules=PRF(x=rulegen_data[,X_names],y=rulegen_data$res,
              ntree=ntree,rounding.digits=rounding.digits,
              maxdepth=maxdepth,disaggregate=disaggregate)

    #Filter out duplicate rules
    rules=filter_rules(rules,X,filtering_tol=0,
                       priority=1/rowSums(RuleMats(rules,X)$RulePredMat))
  }


  #Extend X_lin with the rules and standardize.
  #Since intercept is controlled externally, you can't just take the linear terms
  #from get_modmat, cause those change according to whether formula contains
  #the intercept or not. So, to be safe, attach rules to X_lin manually
  X_rules=pre:::get_modmat(formula = formula, data = data,
                           rules = rules, type = "both",
                           x_names = X_names, winsfrac = 0,
                           normalize = F, y_names = y_name)$x
  X_rules=X_rules[,(ncol(X_rules)-length(rules)+1):ncol(X_rules)]

  #Store info of X before rescaling
  muX=c(muX_lin,apply(X_rules,MARGIN=2,FUN=mean))
  X_mat=cbind(X_lin,scale(X_rules,scale=F))




  ####### Start the fitting process

  #Compute structured penalization for rules
  penaliz=(2*pmin(muX[-(1:p)],1-muX[-(1:p)]))^mu/
    sqrt(2*pmax(muX[-(1:p)],1-muX[-(1:p)]))/
    (rowSums(RuleMats(rules,X)$RulePredMat)^eta)
  prior_lambda=c(rep(lin.relax,p),penaliz)

  #Use CV Ridge regression to guesstimate a good MCMC starting point
  Ridge_fit=glmnet::cv.glmnet(x=t(t(X_rules)*penaliz),alpha=0,
                      y=rulegen_data$res,standardize=F,
                      family='gaussian')
  Beta_init=c(HS_fit$BetaHat,Ridge_fit$glmnet.fit$beta[,Ridge_fit$index[2]]*penaliz)
  Lambda_init=c(HS_fit$LambdaHat,rep(1,ncol(X_rules)))
  Sigma2_init=Ridge_fit$cvm[Ridge_fit$index[2]]
  TauRules_init=c(sqrt(Sigma2_init/n/Ridge_fit$lambda[Ridge_fit$index[2]]))
  TauLin_init=c(HS_fit$TauLinHat)


  #Fit model (it's okay if you specify method.sigma and Sigma2 for logistic,
  #hs() ignores it anyhow)
  if(verbose){print("Fitting model")}

  HS_fit=hs(y=y, X=X_mat, family=family, method.tau.lin = method.tau.lin,
            method.tau.rules = 'halfCauchy', p.lin = p,
            method.sigma = "Jeffreys", burn.in = burn.in, nmc = nmc,
            thin = thin, prior_lambda = prior_lambda, ping=ping,
            Beta_init=Beta_init,Sigma2=Sigma2_init,Lambda_init=Lambda_init,
            tau.lin=TauLin_init,tau.rules=TauRules_init,verbose=verbose)

  #For gaussian case, Y was standardized.
  #Keep that into account, also for intercept
  if(family=='gaussian'){
    HS_fit$BetaSamples=HS_fit$BetaSamples*SDy
    if(intercept){
      HS_fit$BetaSamples[1,]=HS_fit$BetaSamples[1,]+MuY
    }
  }

  #Define output
  output=HS_fit
  output$rules=rules
  output$wins.points=wins.points

  #X was also centered and partially scaled.
  #Keep that into account for intercept
  if(intercept){
    output$BetaSamples[2:p,]=output$BetaSamples[2:p,]/sdX_lin[-1]
    output$BetaSamples[1,]=
      output$BetaSamples[1,]-muX[-1]%*%output$BetaSamples[-1,]
  } else{
    output$BetaSamples[1:p,]=output$BetaSamples[1:p,]/sdX_lin
  }


  #Compute estimated coefficients
  output$BetaHat=rowMeans(output$BetaSamples)

  #Return fit model
  return(output)
}












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
#' @param tau.lin a numeric entry denoting which value the global shrinkage for linear terms should take,
#' in the case where \code{method.tau.lin} is set to \code{"fixed"}. If other values of \code{method.tau.lin}
#' are chosen, this parameter is used as initial value to start the Gibbs sampling from.
#' @param method.tau character string denoting the prior distribution of the
#' overall global shrinkage assigned to all terms: the choice is between a half-Cauchy prior
#' (\code{"halfCauchy"}),a truncated Cauchy prior (\code{"truncatedCauchy"}) or
#' a fixed value (\code{"fixed"})
#' @param tau a numeric entry denoting which value the global shrinkage for linear terms should take,
#' in the case where \code{method.tau} is set to \code{"fixed"}. If other values of \code{method.tau}
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
#' @param verbose Logical value determining whether updates should be given on the progress
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
hs=function(y, X, family=c('gaussian',"binomial"), p.lin=ncol(X),
            method.tau.rules = c("halfCauchy", "truncatedCauchy", "fixed"),
            tau.rules = 1, method.tau.lin = c("halfCauchy", "truncatedCauchy", "fixed"),
            tau.lin = 1, method.tau = c("halfCauchy", "truncatedCauchy", "fixed"),
            tau = 1, method.sigma = c("Jeffreys","fixed"), Sigma2 = 1,
            burn.in = 1000, nmc = 5000, thin = 1, prior_lambda = NULL,
            Beta_init=rep(0,ncol(X)),Lambda_init=rep(1,ncol(X)),ping=1e3,
            verbose=F){

  #Retrieve arguments and basic quantities
  family=match.arg(family)
  method.tau = match.arg(method.tau)
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
  taulinout <- tauruleout <- tauout <- sigma2out <- rep(0, effsamp)

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


  if(verbose){message("Markov chain monte carlo is running")}
  for (i in 1:niter) {
    #Update lambda

    #Combine the vector of tau.rules for rules and tau.lin for linear terms
    tau_vec=c(rep(tau.lin,p.lin),rep(tau.rules,p-p.lin))*tau

    tempps = (Beta/(tau_vec*prior_lambda))^2/2/Sigma2
    nu.inv=stats::rexp(p, 1+1/lambda^2)
    lambda=sqrt(1/stats::rexp(p, nu.inv + tempps))


    #Update tau.lin
    if (method.tau.lin == "halfCauchy") {
      xi.inv = rexp(1, 1 + 1/tau.lin^2)
      tau.lin = sqrt(1/rgamma(1, shape = (p.lin + 1)/2,
                              rate = xi.inv + sum(Beta[1:p.lin]^2/(2 * (lambda[1:p.lin]*prior_lambda[1:p.lin]*tau)^2 * Sigma2))))
    }
    if (method.tau.lin == "truncatedCauchy") {
      #as in horseshoe
      tempt <- sum(Beta[1:p.lin]^2/(2 * (lambda[1:p.lin]*prior_lambda[1:p.lin]*tau)^2 * Sigma2))
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
      xi.inv = rexp(1, 1 + 1/tau.rules^2)
      tau.rules = sqrt(1/rgamma(1, shape = (p-p.lin + 1)/2,
                                rate = xi.inv + sum(Beta[-(1:p.lin)]^2/(2 * (lambda[-(1:p.lin)]*prior_lambda[-(1:p.lin)]*tau)^2 * Sigma2))))
    }
    if (method.tau.rules == "truncatedCauchy") {
      #as in horseshoe
      tempt <- sum(Beta[-(1:p.lin)]^2/(2 * (lambda[-(1:p.lin)]*prior_lambda[-(1:p.lin)]*tau)^2 * Sigma2))
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



    #Update tau
    if (method.tau == "halfCauchy") {
      xi.inv = 1/rexp(1, 1 + 1/tau^2)
      tau = sqrt(1/rgamma(1, shape = (p + 1)/2,
                          rate = xi.inv + sum(Beta^2/(2 * (lambda*prior_lambda*c(rep(tau.lin,p.lin),rep(tau.rules,p-p.lin)))^2 * Sigma2))))
    }
    if (method.tau == "truncatedCauchy") {
      #as in horseshoe
      tempt <- sum(Beta^2/(2 * (lambda*prior_lambda*c(rep(tau.lin,p.lin),rep(tau.rules,p-p.lin)))^2 * Sigma2))
      et = 1/tau^2
      utau = stats::runif(1, 0, 1/(1 + et))
      ubt_1 = 1
      ubt_2 = min((1 - utau)/utau, p^2)
      Fubt_1 = stats::pgamma(ubt_1, (p + 1)/2, scale = 1/tempt)
      Fubt_2 = stats::pgamma(ubt_2, (p + 1)/2, scale = 1/tempt)
      ut = stats::runif(1, Fubt_1, Fubt_2)
      et = stats::qgamma(ut, (p + 1)/2, scale = 1/tempt)
      tau.rules = 1/sqrt(et)
    }


    #Update sigma as in horseshoe
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
    if (verbose & i%%ping == 0) {
      message("iteration = ", i)
    }
    if (i > burn.in && i%%thin == 0) {
      betaout[, (i - burn.in)/thin] <- Beta
      lambdaout[, (i - burn.in)/thin] <- lambda
      taulinout[(i - burn.in)/thin] <- tau.lin
      tauruleout[(i - burn.in)/thin] <- tau.rules
      tauout[(i - burn.in)/thin] <- tau
      sigma2out[(i - burn.in)/thin] <- Sigma2
    }
  }

  result = list(LambdaSamples = lambdaout,
                TauSamples = taulinout,
                TauLinSamples = taulinout,
                TauRulesSamples = tauruleout,
                Sigma2Samples = sigma2out,
                BetaSamples = betaout,
                BetaHat= rowMeans(betaout),
                TauHat= mean(tauout),
                TauLinHat= mean(taulinout),
                TauRulesHat= mean(tauruleout),
                LambdaHat = rowMeans(lambdaout),
                Sigma2Hat = mean(sigma2out))
  return(result)
}
