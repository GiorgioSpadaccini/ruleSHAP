n=1e3
p=30
formula = y ~ .
nrep=2

## Load packages
library("pre")
library('dplyr')
library('ruleSHAP')
library("partykit")
library("randomForest")
library("pre")
library("ggplot2")
library("ggthemes")
library('treeshap')
library('xgboost')
library('horseshoe')
library('horseshoenlm')
library('ggforce')
library("randomForestSRC")




set.seed(1)


RejVec=rep(0,2)
HowManyNoise=c(nrep*(p-5),0)

for(i in 1:nrep){
  #Generate data
  mydata=gendata.friedman1(n=n,p=p,sd=10,rho=0.3+diag(0.7,nrow=p,ncol=p))
  
  
  
  #Use RF-VImp
  RF_fit <- rfsrc(formula, data = mydata)
  reg.smp.o <- subsample(RF_fit)
  CIs=extract.subsample(reg.smp.o, alpha = .05, target = 0,
                        m.target = 'y', standardize = T, raw = TRUE)$ci.jk.Z[c(1,5),-(1:5)]
  RejVec[1]=RejVec[1]+sum(apply(CIs,2,function(x){sign(x[1])==sign(x[2])}))
  
  
  
  #Use two-step procedure
  param <- list(objective = "reg:squarederror",
                max_depth = 5,
                eta = 0.01)
  xgb_model <- xgboost::xgboost(as.matrix(mydata[,1:p]), params = param, label = mydata$y,
                                nrounds = 500, verbose = 0)
  
  #Base it off of shapleys
  unified <- treeshap::unify(xgb_model, mydata)
  treeshap2 <- treeshap::treeshap(unified, mydata, verbose = 0)
  shap_plot=treeshap::plot_feature_importance(treeshap2, max_vars=8)
  lm_model=lm(formula,mydata[,c(as.character(shap_plot$data$variable),'y')])
  pvals=summary(lm_model)$coefficients[,4]
  RejVec[2]=RejVec[2]+sum(pvals[!(names(pvals) %in% c('(Intercept)',paste0('x.',1:5)))]<.05)
  HowManyNoise[2]=HowManyNoise[2]+sum(!(names(pvals) %in% c('(Intercept)',paste0('x.',1:5))))
}


names(RejVec)=names(HowManyNoise)=c('RF-VImp','Two-step')
print('Absolute frequencies:')
print(RejVec)
print('Total noise predictors tested:')
print(HowManyNoise)
print('Estimated type I error rate:')
print(RejVec/HowManyNoise)
