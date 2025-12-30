lin_relax=100
burn.in=2e3
nmc=2e4
maxdepth=3
rounding.digits=2
mu=0.5
eta=2
n_vec=c(100,300,500,1e3,3e3,5e3)
nrep=5
formula = y ~ .
ntree=500
p=30


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


#Shuffle samples and split into trains and test
set.seed(1)
friedman_trains=gendata.friedman1(n=nrep*max(n_vec),p=p,sd=10,rho=0.3+diag(0.7,nrow=p,ncol=p), binary=T)
friedman_test=gendata.friedman1(n=1e4,p=p,sd=10,rho=0.3+diag(0.7,nrow=p,ncol=p), binary=T)


#Preallocate data frame with predicted outcomes
preds <- data.frame(pre = rep(NA, times = nrow(friedman_test)),
                    true = friedman_test[ , all.vars(formula)[1L]])
preds$lasso=preds$ols=preds$rf=preds$hr1=preds$hr2=preds$ruleshap=preds$pre


#Preallocate matrix to store results in
AUCMat=matrix(ncol= nrep*length(n_vec), nrow=ncol(preds)-1)
LinRegMat=LassoMat=RuleFitMat=HorseRule1Mat=HorseRule2Mat=RuleSHAPMat=
  matrix(ncol= nrep*length(n_vec), nrow=ncol(model.matrix(formula,data=friedman_trains))-1)
rownames(AUCMat)=names(preds)[-2]


#Fit models
set.seed(42)
for(i in 1:nrep){
  for(n_index in 1:length(n_vec)){
    #Set repetition-specific quantities
    set.seed(length(n_vec)*(i-1)+n_index)
    n=n_vec[n_index]
    data = friedman_trains[(i-1)*max(n_vec)+(1:n),]
    
    #Verbose
    print(paste('Fitting on sample size n =',n))
    
    #Fit OLS,LASSO, PRE, RF
    ols <- glm(formula=formula, data=data, family='binomial')
    prel <- cv.glmnet(x=as.matrix(data[,1:p]),y=data$y,family='binomial')
    data_rf=data
    data_rf$y=as.factor(data_rf$y)
    preb <- pre(formula, data = data_rf, family='binomial')
    rf <- randomForest(formula, data = data_rf)
    
    #Fit HorseRule in two ways
    #I'll just assume intercept is included in the formula given in the input
    hr_data=as.data.frame(model.matrix(formula,data))[,-1]
    hr_data[,as.character(formula[[2]])]=data_rf[,as.character(formula[[2]])]
    hr_test=as.data.frame(model.matrix(formula,friedman_test))[,-1]
    
    hr1 <- horserule:::HorseRuleFit(model.formula = formula, data=hr_data,
                                    restricted=0,linterms = 1:(ncol(hr_data)-1),
                                    niter = nmc+burn.in, burnin = burn.in, thin = 1,
                                    ntree = ntree,intercept=T, Xtest=hr_test,
                                    linp=1, ensemble='both',mix=0.3)
    
    hr2 <- horserule:::HorseRuleFit(model.formula = formula, data=hr_data,
                                    restricted=0,linterms = 1:(ncol(hr_data)-1),
                                    niter = nmc+burn.in, burnin = burn.in, thin = 1,
                                    ntree = ntree,intercept=T, Xtest=hr_test,
                                    linp=lin_relax, ensemble='RF')
    
    
    ###########Store predictions
    preds$ols <- predict(ols, newdata = friedman_test,type='response')
    preds$pre <- predict(preb, newdata = friedman_test,type='response')
    preds$lasso <- predict(prel, newx = as.matrix(friedman_test[,1:p]),type='response')
    preds$rf <- predict(rf, newdata = friedman_test,type='prob')[,2]
    preds$hr1 <- hr1$pred
    preds$hr2 <- hr2$pred
    
    
    
    #############Store coefficients
    LinRegMat[,(n_index-1)*nrep+i]=coef(ols)[-1]
    LassoMat[,(n_index-1)*nrep+i]=as.vector(coef(prel)[-1])
    
    #PRE re-orders the coefficients, so you need to use the correct order
    RuleFitMat[,(n_index-1)*nrep+i]=coef(preb)$coefficient[order(as.numeric(rownames(coef(preb))))[1+(1:nrow(RuleFitMat))]]
    
    #From the code of predict.HorseRulemodel we can see that stored coefficients
    #are not rescaled to refer to the standardization process, so that needs to be included
    HorseRule1Mat[,(n_index-1)*nrep+i]=hr1$postmean[1+(1:nrow(HorseRule1Mat))]/
      hr1$modelstuff$sdl
    
    HorseRule2Mat[,(n_index-1)*nrep+i]=hr2$postmean[1+(1:nrow(HorseRule2Mat))]/
      hr2$modelstuff$sdl
    
    
    
    #Now do the same for RuleSHAP
    RS_fit=ruleSHAP(formula=formula, data=data, family='binomial',verbose=T,nmc=nmc,burn.in = burn.in)
    
    #Store predictions
    Ytr=data[, as.character(formula[[2]])]
    Xtr=data[, colnames(data) != as.character(formula[[2]])]
    my_modmat=pre:::get_modmat(formula = formula, data = friedman_test,
                               rules = RS_fit$rules,
                               type = "both", x_names = colnames(Xtr),
                               winsfrac = 0, normalize = F, y_names = as.character(formula[[2]]))
    preds$ruleshap <- c(1/(1+exp(-cbind(1,my_modmat$x) %*% RS_fit$BetaHat)))
    
    #Store coefficients
    RuleSHAPMat[,(n_index-1)*nrep+i]=RS_fit$BetaHat[1+(1:nrow(RuleSHAPMat))]
    
    #Compute and store AUC values
    AUCMat[,(n_index-1)*nrep+i]=apply(preds[,-2],2,function(x){pROC:::roc(preds$true,x)$auc})
    
    #Print AUC
    print(AUCMat[,(n_index-1)*nrep+i])
  }
}


#Save everything
saveRDS(AUCMat,'output/AUCMat_p30_logit.Rda')
saveRDS(LinRegMat,'output/LinRegMat_p30_logit.Rda')
saveRDS(LassoMat,'output/LassoMat_p30_logit.Rda')
saveRDS(RuleFitMat,'output/RuleFitMat_p30_logit.Rda')
saveRDS(HorseRule1Mat,'output/HorseRule1Mat_p30_logit.Rda')
saveRDS(HorseRule2Mat,'output/HorseRule2Mat_p30_logit.Rda')
saveRDS(RuleSHAPMat,'output/RuleSHAPMat_p30_logit.Rda')









#Plot result:

#Plot coefficients:
#Create dataframe
nmethods=6
p_show=6
models=c('OLS','LASSO','RuleSHAP','RuleFit','HR1','HR2','BART','RF')
df=data.frame(coef=c(readRDS('output/RuleFitMat_p30_logit.Rda')[1:p_show,],
                     readRDS('output/RuleSHAPMat_p30_logit.Rda')[1:p_show,],
                     readRDS('output/HorseRule1Mat_p30_logit.Rda')[1:p_show,],
                     readRDS('output/HorseRule2Mat_p30_logit.Rda')[1:p_show,],
                     readRDS('output/LinRegMat_p30_logit.Rda')[1:p_show,],
                     readRDS('output/LassoMat_p30_logit.Rda')[1:p_show,]),
              model=rep(c('RuleFit','RuleSHAP','HR1','HR2','OLS','LASSO'),
                        each=nrep*length(n_vec)*p_show),
              rep=rep(1:nrep,each=nmethods),
              n=rep(n_vec,each=nrep*nmethods),predictor=paste0('x.',1:p_show))

#Define target (use lm() for first two entries)
int_coef=mean(coef(lm(formula=y~.,data=gendata.friedman1(n=1e6,p=6,sd=0)))[2:3])
lines_df=data.frame(predictor=paste0('x.',1:6),y=c(rep(int_coef,2),0,10,5,0),model=NA,
                    Reference='Reference')


#Plot it
ggplot(df, aes(x=predictor, y=coef, fill=factor(model,levels=models))) +
  facet_grid(factor(paste('n =',n),levels=paste('n =',n_vec))~.)+
  geom_boxplot(width=0.5, position = position_dodge(width = .75))+
  stat_summary(fun = median, geom = "point", aes(group = factor(model,levels=models),shape = factor(model,levels=models)),
               size = 3, position = position_dodge(width = .75))+
  theme(legend.position = "top", 
        legend.box = "horizontal",
        legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  geom_errorbar(data=lines_df,aes(x=predictor,y=y,ymin=y,ymax=y,color=Reference), group='model',
                linetype = 'dashed')+
  xlab('Predictor')+ylab('Coefficient')+
  scale_fill_manual(values=c(colorspace::rainbow_hcl(length(models))),
                    breaks=models,name='Model')+
  scale_shape_manual(values=1:(length(models)),
                     breaks=models, name='Model')+ylim(c(-10,15))





#Plot predictive performance
#Load dataframes and join them
nmethods=8
models=c('OLS','LASSO','RuleSHAP','RuleFit','HR1','HR2','BART','RF')
df=rbind(data.frame(AUC=c(readRDS('output/AUCMat_p10_logit.Rda')),p=10,
                    model=c('RuleFit','RuleSHAP','HR2','HR1','RF','OLS','LASSO'),
                    rep=rep(1:nrep,each=nmethods),
                    n=rep(n_vec,each=nrep*nmethods)),
         data.frame(AUC=c(readRDS('output/AUCMat_p30_logit.Rda')),p=30,
                    model=c('RuleFit','RuleSHAP','HR2','HR1','RF','OLS','LASSO'),
                    rep=rep(1:nrep,each=nmethods),
                    n=rep(n_vec,each=nrep*nmethods)))
#Plot it
df=df[complete.cases(df),]
ggplot(df, aes(x=factor(model,levels=models), y=AUC,
               fill=factor(model,levels=models))) +
  facet_grid(paste('p =',p) ~ 
               factor(paste('n =',n),levels=paste('n =',n_vec)))+
  geom_boxplot(width=0.5)+
  stat_summary(fun = median, geom = "point", aes(group = factor(model,levels=models),
                                                 shape = factor(model,levels=models)),
               size = 3)+
  theme(legend.position = "top", 
        legend.box = "horizontal",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(nrow = 1))+
  scale_fill_manual(values=c(colorspace::rainbow_hcl(length(models))),
                    breaks=models,name='Model')+
  scale_shape_manual(values=1:(length(models)),
                     breaks=models, name='Model')+
  xlab('Model')+ylab('Area Under the ROC Curve')
