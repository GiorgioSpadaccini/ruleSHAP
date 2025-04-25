lin_relax=100
burn.in=2e3
nmc=2e4
maxdepth=3
rounding.digits=2
mu=0.5
eta=2
n_vec=c(100,300,500,1e3,3e3,5e3)
nrep=5
formula = sbp ~ .
ntree=500


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
library('glmnet')

## Load data
fulldata=readRDS('helius.Rda')

#Compute number of predictors
p=ncol(fulldata)-2

#Select outcome
heliusdata=cbind(fulldata[,1:p],sbp=fulldata$sbp)

#Shuffle samples and split into trains and test
set.seed(1)
helius_data=heliusdata[sample(nrow(heliusdata)),]
helius_trains=helius_data[1:(max(n_vec[-6])*nrep),]
helius_test=helius_data[-(1:(max(n_vec[-6])*nrep)),]


#Preallocate data frame with predicted outcomes
preds <- data.frame(pre = rep(NA, times = nrow(helius_test)),
                    true = helius_test[ , all.vars(formula)[1L]])
preds$tree=preds$lasso=preds$ols=preds$rf=preds$hr1=preds$hr2=preds$ruleshap=preds$pre


#Preallocate matrix to store results in
MSEMat=matrix(ncol= nrep*length(n_vec), nrow=ncol(preds)-1)
LinRegMat=LassoMat=RuleFitMat=HorseRule1Mat=HorseRule2Mat=RuleSHAPMat=
  matrix(ncol= nrep*length(n_vec), nrow=ncol(model.matrix(formula,data=helius_trains))-1)
rownames(MSEMat)=names(preds)[-2]


#Fit models
set.seed(42)
for(i in 1:nrep){
  for(n_index in 1:length(n_vec)){
    #Set repetition-specific quantities
    set.seed(length(n_vec)*(i-1)+n_index)
    n=n_vec[n_index]
    data = helius_trains[(i-1)*max(n_vec[-6])+(1:n),]
    modmat=model.matrix(formula,data=data)[,-1]
    
    #Verbose
    print(paste('Fitting on sample size n =',n))
    
    #Fit OLS,LASSO, PRE, RF and CTree
    ols <- lm(formula=formula, data=data)
    preb <- pre(formula, data = data)
    prel <- cv.glmnet(x=modmat,y=data$sbp)
    rf <- randomForest(formula, data = data)
    tree <- ctree(formula, data = data)
    
    #Fit HorseRule
    #I'll just assume intercept is included in the formula given in the input
    hr_data=as.data.frame(model.matrix(formula,data))[,-1]
    hr_data[,as.character(formula[[2]])]=data[,as.character(formula[[2]])]
    hr_test=as.data.frame(model.matrix(formula,helius_test))[,-1]
    
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
    preds$ols <- predict(ols, newdata = helius_test)
    preds$pre <- predict(preb, newdata = helius_test)
    preds$lasso <- predict(prel, newx=model.matrix(formula,data=helius_test)[,-1])
    preds$rf <- predict(rf, newdata = helius_test)
    preds$tree <- predict(tree, newdata = helius_test)
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
    RS_fit=ruleSHAP(formula=formula, data=data, family='gaussian',verbose=T,
                    nmc=nmc,burn.in = burn.in)
    
    #Store predictions
    Ytr=data[, as.character(formula[[2]])]
    Xtr=data[, colnames(data) != as.character(formula[[2]])]
    my_modmat=pre:::get_modmat(formula = formula, data = helius_test,
                               rules = RS_fit$rules,
                               type = "both", x_names = colnames(Xtr),
                               winsfrac = 0, normalize = F, y_names = as.character(formula[[2]]))
    preds$ruleshap <- c(cbind(1,my_modmat$x) %*% RS_fit$BetaHat)
    
    #Store coefficients
    RuleSHAPMat[,(n_index-1)*nrep+i]=RS_fit$BetaHat[1+(1:nrow(RuleSHAPMat))]
    
    #Compute and store MSE values
    MSEMat[,(n_index-1)*nrep+i]=colMeans((preds[,-2] - preds$true)^2)
    
    #Print R squared values
    print(1-MSEMat[,(n_index-1)*nrep+i]/var(helius_test$sbp))
  }
}


#Save everything
saveRDS(MSEMat,'output/MSEMat_sbp.Rda')
saveRDS(LinRegMat,'output/LinRegMat_sbp.Rda')
saveRDS(LassoMat,'output/LassoMat_sbp.Rda')
saveRDS(RuleFitMat,'output/RuleFitMat_sbp.Rda')
saveRDS(HorseRule1Mat,'output/HorseRule1Mat_sbp.Rda')
saveRDS(HorseRule2Mat,'output/HorseRule2Mat_sbp.Rda')
saveRDS(RuleSHAPMat,'output/RuleSHAPMat_sbp.Rda')










#Plot results:

#Plot coefficients:
#Create dataframe
nmethods=6
p_show=6
models=c('OLS','LASSO','RuleSHAP','RuleFit','HR1','HR2','BART','RF','cTree')
df=data.frame(coef=c(readRDS('output/RuleFitMat_sbp.Rda')[1:p_show,],
                     readRDS('output/RuleSHAPMat_sbp.Rda')[1:p_show,],
                     readRDS('output/HorseRule1Mat_sbp.Rda')[1:p_show,],
                     readRDS('output/HorseRule2Mat_sbp.Rda')[1:p_show,],
                     readRDS('output/LinRegMat_sbp.Rda')[1:p_show,],
                     readRDS('output/LassoMat_sbp.Rda')[1:p_show,]),
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
                     breaks=models, name='Model')





#Plot predictive performance
#Load dataframe
nmethods=8
models=c('OLS','LASSO','RuleSHAP','RuleFit','HR1','HR2','BART','RF','cTree')
df=data.frame(MSE=c(readRDS('output/MSEMat_sbp.Rda')),
              model=c('RuleFit','RuleSHAP','HR2','HR1','RF','OLS','LASSO','cTree'),
              rep=rep(1:nrep,each=nmethods),
              n=rep(n_vec,each=nrep*nmethods))
#Plot it
ggplot(df, aes(x=factor(model,levels=models), y=MSE,
               fill=factor(model,levels=models))) +
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
  xlab('Model')+ylab('MSE')