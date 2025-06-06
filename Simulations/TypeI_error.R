lin_relax=100
burn.in=2e3
nmc=2e4
maxdepth=3
rounding.digits=2
mu=0.5
eta=2
p_vec=c(10,30,50)
nrep=100
formula = y ~ .
ntree=500
n=1e3



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
friedman_trains=gendata.friedman1(n=nrep*n,p=max(p_vec),sd=10,rho=0.3+diag(.7,nrow=max(p_vec),ncol=max(p_vec)))
friedman_test=gendata.friedman1(n=1e4,p=max(p_vec),sd=10,rho=0.3+diag(.7,nrow=max(p_vec),ncol=max(p_vec)))

ruleshap_dfs=bart_dfs=hr1_dfs=hr2_dfs=list()


#Fit models
set.seed(42)
for(i in 1:nrep){
  for(p_index in 1:length(p_vec)){
    
    #Set repetition-specific quantities
    set.seed(length(p_vec)*(i-1)+p_index)
    p=p_vec[p_index]
    data = friedman_trains[(i-1)*n+(1:n),c(1:p,max(p_vec)+1)]
    
    #Verbose
    print(paste('Fitting on sample with p =',p))
    
    
    #Start off with BART
    bartFit <- dbarts::bart(data[,1:p], data$y, keeptrees = TRUE, verbose = FALSE)
    bart_dfs[[(p_index-1)*nrep+i]]=bart_df=bartSHAP(bartFit,data[,1:p])
    
    #Now Horserule
    hr_data=as.data.frame(model.matrix(formula,data))[,-1]
    hr_data[,as.character(formula[[2]])]=data[,as.character(formula[[2]])]
    
    hr1 <- horserule:::HorseRuleFit(model.formula = formula, data=hr_data,
                                    restricted=0,linterms = 1:(ncol(hr_data)-1),
                                    niter = nmc+burn.in, burnin = burn.in, thin = 1,
                                    ntree = ntree,intercept=T,
                                    linp=1, ensemble='both',mix=0.3)
    
    hr2 <- horserule:::HorseRuleFit(model.formula = formula, data=hr_data,
                                    restricted=0,linterms = 1:(ncol(hr_data)-1),
                                    niter = nmc+burn.in, burnin = burn.in, thin = 1,
                                    ntree = ntree,intercept=T,
                                    linp=lin_relax, ensemble='RF')
    
    
    #Shapleys for HR1
    hr1_new=list(BetaSamples=t(hr1$postdraws$beta[(1:burn.in),])/c(1,hr1$modelstuff$sdl,hr1$modelstuff$sdr)*hr1$modelstuff$sdy,
                 rules=gsub("]"," ",gsub("X[,","x.",hr1$rules,fixed=T),fixed=T),wins.points=matrix(ncol=p,nrow=2))
    shapleys=compute_SHAP(hr1_new,data[,1:p])
    hr1_dfs[[(p_index-1)*nrep+i]]=hr1_df=shapleys$marginal
    
    #Shapleys for HR2
    hr2_new=list(BetaSamples=t(hr2$postdraws$beta[(1:burn.in),])/c(1,hr2$modelstuff$sdl,hr2$modelstuff$sdr)*hr2$modelstuff$sdy,
                 rules=gsub("]"," ",gsub("X[,","x.",hr2$rules,fixed=T),fixed=T),wins.points=matrix(ncol=p,nrow=2))
    shapleys=compute_SHAP(hr2_new,data[,1:p])
    hr2_dfs[[(p_index-1)*nrep+i]]=hr2_df=shapleys$marginal
    
    
    #Now RuleSHAP
    RS_fit=ruleSHAP(formula=formula, data=data,family='gaussian',verbose=T,burn.in = burn.in,nmc = nmc)
    ruleshap_dfs[[(p_index-1)*nrep+i]]=ruleshap_df=compute_SHAP(RS_fit,data[,1:p])$marginal
    
    
    
    #### Print counts as a preview
    bart_df$sig=sign(bart_df$CIinf*bart_df$CIsup)==1
    ruleshap_df$sig=sign(ruleshap_df$CIinf*ruleshap_df$CIsup)==1
    hr1_df$sig=sign(hr1_df$CIinf*hr1_df$CIsup)==1
    hr2_df$sig=sign(hr2_df$CIinf*hr2_df$CIsup)==1
    
    bart_counts=bart_df %>% group_by(predictor) %>% summarize(npoints=sum(sig))
    ruleshap_counts=ruleshap_df %>% group_by(predictor) %>% summarize(npoints=sum(sig))
    hr1_counts=hr1_df %>% group_by(predictor) %>% summarize(npoints=sum(sig))
    hr2_counts=hr2_df %>% group_by(predictor) %>% summarize(npoints=sum(sig))
    
    print('BART')
    print(head(bart_counts[order(bart_counts$npoints,decreasing = T),],15))
    print('HR1')
    print(head(hr1_counts[order(hr1_counts$npoints,decreasing = T),],15))
    print('HR2')
    print(head(hr2_counts[order(hr2_counts$npoints,decreasing = T),],15))
    print('ruleSHAP')
    print(head(ruleshap_counts[order(ruleshap_counts$npoints,decreasing = T),],15))
  }
}







#Show results

SigCounts=NoiseCounts=matrix(0,ncol=3,nrow=4)
rownames(SigCounts)=rownames(NoiseCounts)=c('BART','HR1','HR2','RuleSHAP')
p_vec=rep(p_vec,each=nrep)

for(i in 1:length(p_vec)){
  p=p_vec[i]
  
  bart_df=bart_dfs[[i]]
  ruleshap_df=ruleshap_dfs[[i]]
  hr1_df=hr1_dfs[[i]]
  hr2_df=hr2_dfs[[i]]
  
  if(is.null(bart_df)){next}
  
  bart_df$sig=sign(bart_df$CIinf*bart_df$CIsup)==1
  ruleshap_df$sig=sign(ruleshap_df$CIinf*ruleshap_df$CIsup)==1
  hr1_df$sig=sign(hr1_df$CIinf*hr1_df$CIsup)==1
  hr2_df$sig=sign(hr2_df$CIinf*hr2_df$CIsup)==1
  
  bart_counts=bart_df %>% group_by(predictor) %>% summarize(npoints=sum(sig))
  ruleshap_counts=ruleshap_df %>% group_by(predictor) %>% summarize(npoints=sum(sig))
  hr1_counts=hr1_df %>% group_by(predictor) %>% summarize(npoints=sum(sig))
  hr2_counts=hr2_df %>% group_by(predictor) %>% summarize(npoints=sum(sig))
  
  bart_counts$signal=bart_counts$predictor %in% paste0('x.',1:5)
  ruleshap_counts$signal=ruleshap_counts$predictor %in% paste0('x.',1:5)
  hr1_counts$signal=hr1_counts$predictor %in% paste0('x.',1:5)
  hr2_counts$signal=hr2_counts$predictor %in% paste0('x.',1:5)
  
  bart_counts=bart_counts %>% group_by(signal) %>% summarize(nvars=sum(npoints>0))
  ruleshap_counts=ruleshap_counts %>% group_by(signal) %>% summarize(nvars=sum(npoints>0))
  hr1_counts=hr1_counts %>% group_by(signal) %>% summarize(nvars=sum(npoints>0))
  hr2_counts=hr2_counts %>% group_by(signal) %>% summarize(nvars=sum(npoints>0))
  
  SigCounts[,ceiling(i/100)]=SigCounts[,ceiling(i/100)]+c(bart_counts$nvars[bart_counts$signal],
                                                          hr1_counts$nvars[hr1_counts$signal],
                                                          hr2_counts$nvars[hr2_counts$signal],
                                                          ruleshap_counts$nvars[ruleshap_counts$signal])
  NoiseCounts[,ceiling(i/100)]=NoiseCounts[,ceiling(i/100)]+c(bart_counts$nvars[!bart_counts$signal],
                                                              hr1_counts$nvars[!hr1_counts$signal],
                                                              hr2_counts$nvars[!hr2_counts$signal],
                                                              ruleshap_counts$nvars[!ruleshap_counts$signal])
}

#Power:
SigCounts/(5*nrep)
#Type I error:
NoiseCounts/rep(100*c(5,25,45),each=4)
#FDR:
NoiseCounts/(SigCounts+NoiseCounts)






#Example of BART vs. ruleSHAP
library(ggplot2)
shapleys_df=rbind(cbind(bart_dfs[[13]],
                        model='BART'),
                  cbind(ruleshap_dfs[[13]],
                        model='RuleSHAP'))
shapleys_df$sig=sign(shapleys_df$CIinf*shapleys_df$CIsup)==1
X_names=paste0('x.',1:10)

ggplot(shapleys_df[shapleys_df$predictor %in% X_names,], aes(x = x, y = value)) +
  geom_linerange(aes(x = x, ymin = CIinf, ymax = CIsup, color=sig),size=1, show.legend=F)+
  geom_point(size=0.5) +
  geom_line()+
  facet_grid(model~factor(predictor,levels=X_names), scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 10))+
  scale_color_manual(name = "Significance", 
                     values = c("TRUE" = "lightblue", "FALSE" = "#cfcfcf"))+
  labs(x = "Feature value", y = "Shapley value")+
  scale_x_continuous(breaks=c(0,0.50,1), labels=c(0,0.5,1))
