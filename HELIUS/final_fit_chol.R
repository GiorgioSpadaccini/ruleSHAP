burn.in=2e3
nmc=2e4
maxdepth=3
rounding.digits=2
mu=0.5
eta=2
n=1e4
formula = chol ~ .
ntree=500
family='gaussian'
verbose=T
intercept=T
disaggregate=T
method.tau.lin="halfCauchy"
thin=1
ping=1e3
rules=NULL
wins.frac=0.025


## Load packages and data
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
fulldata=readRDS('helius.Rda')

#Compute number of predictors
p=ncol(fulldata)-2

#Select outcome
heliusdata=cbind(fulldata[,1:p],chol=fulldata$chol)

#Shuffle samples and split into trains and test
set.seed(1)
base_ind=sample(nrow(heliusdata))[1:10000]
helius_train=heliusdata[base_ind,]


#Fit model's first 10k draws (it's too slow to do it all at once)
set.seed(42)

data=helius_train

#Retrieve basic quantities
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

#Retrieve basic quantities
y_name=as.character(formula[[2]])
X_names=labels(terms(formula, data=data))
y=data[, y_name]
X=data[, X_names]
data=data[,c(X_names,y_name)]


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




#Use the linear terms to compute residuals
if(verbose){print('Generating rules')}
RS_fit=hs(y=y, X=X_lin,family=family,
          method.tau.lin=method.tau.lin, tau.lin=1,
          method.tau.rules = 'fixed', tau.rules = 1,
          method.tau = 'fixed', tau = 1,
          burn=burn.in, nmc=nmc,
          method.sigma="Jeffreys")

#Round numeric predictors and create the dataset for rule generation
rulegen_data=X
which_numeric=sapply(X,is.numeric)
rulegen_data[,which_numeric]=round(X[,which_numeric],digits=rounding.digits)
if(family=='gaussian'){
  rulegen_data$res=y-c(X_lin%*%RS_fit$BetaHat)
} else if(family=='binomial'){
  #For now I'll define the working residuals
  phat=1/(1+exp(-c(X_lin%*%RS_fit$BetaHat)))
  rulegen_data$res=(y-phat)/(phat*(1-phat))
}


#Fit parametric random forest on residuals and filter out duplicate rules
if(is.null(rules)){
  rules=PRF(x=rulegen_data[,X_names],y=rulegen_data$res,
            ntree=ntree,rounding.digits=rounding.digits,
            maxdepth=maxdepth,disaggregate=disaggregate)
}
rules=filter_rules(rules,X,filtering_tol=0,
                   priority=1/rowSums(RuleMats(rules,X)$RulePredMat))

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
prior_lambda=c(rep(1,p),penaliz)

#Use CV Ridge regression to guesstimate a good MCMC starting point
Ridge_fit=cv.glmnet(x=t(t(X_rules)*penaliz),alpha=0,
                    y=rulegen_data$res,standardize=F,
                    family='gaussian')
Beta_init=c(RS_fit$BetaHat,Ridge_fit$glmnet.fit$beta[,Ridge_fit$index[2]]*penaliz)
Lambda_init=c(RS_fit$LambdaHat,rep(1,ncol(X_rules)))
Sigma2_init=Ridge_fit$cvm[Ridge_fit$index[2]]
TauRules_init=c(sqrt(Sigma2_init/n/Ridge_fit$lambda[Ridge_fit$index[2]]))
TauLin_init=c(RS_fit$TauLinHat)


#Fit model
if(verbose){print("Fitting model")}













#Fit first 10k observations and save them

RS_fit_pt1=hs(y=y, X=X_mat, family=family, method.tau.lin = method.tau.lin,
              method.tau.rules = 'halfCauchy', p.lin = p,
              method.sigma = "Jeffreys", burn = burn.in, nmc = nmc/2,
              thin = thin, prior_lambda = prior_lambda, ping=ping,
              Beta_init=Beta_init,Sigma2=Sigma2_init,Lambda_init=Lambda_init,
              tau.lin=TauLin_init,tau.rules=TauRules_init)

saveRDS(RS_fit_pt1,'output/RS_fit_pt1_chol.Rda')    
saveRDS(rules,'output/rules_chol.Rda')
saveRDS(base_ind,'output/base_ind_chol.Rda')












#Now fit second half and finish off the rest      

#Retrieve where you left off
RS_fit_pt1=readRDS('output/RS_fit_pt1_chol.Rda')

#Now fit the second part
set.seed(1)
RS_fit_pt2=hs(y=y, X=X_mat, family=family, method.tau.lin = method.tau.lin,
              method.tau.rules = 'halfCauchy', p.lin = p,
              method.sigma = "Jeffreys", burn = 0, nmc = nmc/2,
              thin = thin, prior_lambda = prior_lambda, ping=ping,
              Beta_init=RS_fit_pt1$BetaSamples[,nmc/2],
              Sigma2=RS_fit_pt1$Sigma2Samples[nmc/2],
              Lambda_init=RS_fit_pt1$LambdaSamples[,nmc/2],
              tau.lin=RS_fit_pt1$TauLinSamples[nmc/2],
              tau.rules=RS_fit_pt1$TauRulesSamples[nmc/2],
              tau=RS_fit_pt1$TauRulesSamples[nmc/2])

saveRDS(RS_fit_pt2,'output/RS_fit_pt2_chol.Rda')


#Combine the two parts together
RS_fit=list()
RS_fit$LambdaSamples=cbind(RS_fit_pt1$LambdaSamples,RS_fit_pt2$LambdaSamples)
RS_fit$TauLinSamples=c(RS_fit_pt1$TauLinSamples,RS_fit_pt2$TauLinSamples)
RS_fit$TauRulesSamples=c(RS_fit_pt1$TauRulesSamples,RS_fit_pt2$TauRulesSamples)
RS_fit$BetaSamples=cbind(RS_fit_pt1$BetaSamples,RS_fit_pt2$BetaSamples)
RS_fit$Sigma2Samples=c(RS_fit_pt1$Sigma2Samples,RS_fit_pt2$Sigma2Samples)
RS_fit$BetaHat=rowMeans(cbind(RS_fit_pt1$BetaHat,RS_fit_pt2$BetaHat))
RS_fit$TauLinHat=mean(c(RS_fit_pt1$LambdaSamples,RS_fit_pt2$LambdaSamples))
RS_fit$TauRulesHat=mean(c(RS_fit_pt1$TauRulesHat,RS_fit_pt2$TauRulesHat))
RS_fit$LambdaHat=rowMeans(cbind(RS_fit_pt1$LambdaHat,RS_fit_pt2$LambdaHat))
RS_fit$Sigma2Hat=mean(c(RS_fit_pt1$Sigma2Hat,RS_fit_pt2$Sigma2Hat))


#For gaussian case, Y was standardized.
#Keep that into account, also for intercept
if(family=='gaussian'){
  RS_fit$BetaSamples=RS_fit$BetaSamples*SDy
  if(intercept){
    RS_fit$BetaSamples[1,]=RS_fit$BetaSamples[1,]+MuY
  }
}

#Define output
RS_fit$rules=rules
RS_fit$wins.points=wins.points

#X was also centered and partially scaled.
#Keep that into account for intercept
if(intercept){
  RS_fit$BetaSamples[2:p,]=RS_fit$BetaSamples[2:p,]/sdX_lin[-1]
  RS_fit$BetaSamples[1,]=
    RS_fit$BetaSamples[1,]-muX[-1]%*%RS_fit$BetaSamples[-1,]
} else{
  RS_fit$BetaSamples[1:p,]=RS_fit$BetaSamples[1:p,]/sdX_lin
}


#Compute estimated coefficients
RS_fit$BetaHat=rowMeans(RS_fit$BetaSamples)

#Save final model
saveRDS(RS_fit,'output/RS_fit_chol.Rda')


#Compute Shapley values of e.g. first 3k observations
p=ncol(fulldata)-2
shapleys_dfs=compute_SHAP(RS_fit,data[1:3e3,1:p],interactions=T)
saveRDS(shapleys_dfs,'output/final_shapleys_chol.Rda')






#Plot heatmap of interactions
library(dplyr)
library(ggplot2)
ruleshap_df=readRDS('output/final_shapleys_chol.Rda')[[2]]
ruleshap_df$sig=sign(ruleshap_df$CIinf*ruleshap_df$CIsup)==1
ruleshap_inters=ruleshap_df %>% group_by(predictor1,predictor2) %>% summarize(sigval=mean(abs(value))*(sum(sig)>0))
ggplot(data = ruleshap_inters, mapping = aes(x = factor(predictor1,levels=X_names),
                                             y = factor(predictor2,levels=X_names),
                                             fill = sigval)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_continuous(low = "white", high = "#56b1f7")+
  theme(legend.position = "top", 
        legend.box = "horizontal",
        legend.direction = "horizontal")+
  guides(fill = guide_colorbar(title='Interaction intensity (significant only)',
                               title.position = 'top',
                               barwidth = 11.5))+
  xlab('Predictor')+ylab('Predictor')+
  geom_hline(yintercept=8.5, linetype='dashed')+
  geom_vline(xintercept = 8.5, linetype='dashed')


#Explore interaction between sex and age
ruleshap_df=readRDS('output/final_shapleys_chol.Rda')[[1]]
sub_df=ruleshap_df[ruleshap_df$predictor %in% c('age','genderf'),]
sub_df=cbind(sub_df[1:(nrow(sub_df)/2),],
             sub_df[-(1:(nrow(sub_df)/2)),])
names(sub_df)[6:10]=paste0(names(sub_df)[1:5],'.2')
sub_df$x.2=c('Male','Female')[(sub_df$x.2+3)/2]
sub_df$val=sub_df$value+sub_df$value.2

lines_df=sub_df %>% group_by(x,x.2) %>% summarize(val=mean(val)) 

ggplot(data=sub_df,
       aes(x = x,y=val,color = factor(x.2), group=factor(x.2)))+
  geom_point(aes(shape=factor(x.2)),alpha=.5)+
  geom_line(data=lines_df, linewidth=1)+
  #theme(legend.position = "top", 
  #      legend.box = "horizontal",
  #      legend.direction = "horizontal")+
  guides(color = guide_legend(nrow = 1))+
  xlab('Age')+ylab('Joint Age + Sex effect')+
  scale_color_manual(values=c(colorspace::rainbow_hcl(2)[2:1]),
                     breaks=c('Male','Female'),name='Sex')+
  scale_shape_manual(values=1:2,
                     breaks=c('Male','Female'), name='Sex')



#Plot overall Shapley values
shapleys_df=readRDS('output/final_shapleys_chol.Rda')[[1]]
shapleys_df$sig=sign(shapleys_df$CIinf*shapleys_df$CIsup)==1

ggplot(shapleys_df, aes(x = x, y = value)) +
  geom_linerange(aes(x = x, ymin = CIinf, ymax = CIsup, color=sig),size=1, show.legend=F)+
  geom_point(size=0.5) +
  geom_line()+
  facet_wrap(~factor(predictor,levels = X_names),ncol=6,scales='free_x')+
  theme_bw()+
  theme(strip.text.x = element_text(size = 10))+
  scale_color_manual(name = "Significance", 
                     values = c("TRUE" = "lightblue", "FALSE" = "#cfcfcf"))+
  labs(x = "Feature value", y = "Shapley value")
