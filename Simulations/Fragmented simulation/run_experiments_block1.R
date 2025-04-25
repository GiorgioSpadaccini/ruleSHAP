lin_relax=100
burn.in=2e3
nmc=2e4
maxdepth=3
rounding.digits=2
mu=0.5
eta=2
p_vec=c(10,30,50)
nrep=100
block_size=5
i_block=1       #Change this block index for different runs
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

shapleys_dfs=shapleys_bart_dfs=hr1_shapleys=hr2_shapleys=list()


#Fit models
set.seed(42)
for(j in 1:block_size){
  i=(i_block-1)*block_size+j
  for(p_index in 1:length(p_vec)){
    
    #Set repetition-specific quantities
    set.seed(length(p_vec)*(i-1)+p_index)
    p=p_vec[p_index]
    data = friedman_trains[(i-1)*n+(1:n),c(1:p,max(p_vec)+1)]
    
    #Verbose
    print(paste('Fitting on sample with p =',p))
    
    
    #Start off with BART
    bartFit <- dbarts::bart(data[,1:p], data$y, keeptrees = TRUE, verbose = FALSE)
    shapleys_bart_dfs[[(p_index-1)*nrep+i]]=bart_df=bartSHAP(bartFit,data[,1:p])
    
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
    hr1_shapleys[[(p_index-1)*nrep+i]]=hr1_df=shapleys$marginal
    
    #Shapleys for HR2
    hr2_new=list(BetaSamples=t(hr2$postdraws$beta[(1:burn.in),])/c(1,hr2$modelstuff$sdl,hr2$modelstuff$sdr)*hr2$modelstuff$sdy,
                 rules=gsub("]"," ",gsub("X[,","x.",hr2$rules,fixed=T),fixed=T),wins.points=matrix(ncol=p,nrow=2))
    shapleys=compute_SHAP(hr2_new,data[,1:p])
    hr2_shapleys[[(p_index-1)*nrep+i]]=hr2_df=shapleys$marginal
    
    
    #Now RuleSHAP
    RS_fit=ruleSHAP(formula=formula, data=data,family='gaussian',verbose=T,burn.in = burn.in,nmc = nmc)
    shapleys_dfs[[(p_index-1)*nrep+i]]=ruleshap_df=compute_SHAP(RS_fit,data[,1:p])$marginal
    
    
    
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
    
    #Save everything
    saveRDS(shapleys_dfs,paste0('output/shapleys_dfs_',i_block,'.Rda'))
    saveRDS(shapleys_bart_dfs,paste0('output/shapleys_bart_dfs_',i_block,'.Rda'))
    saveRDS(hr1_shapleys,paste0('output/hr1_shapleys_',i_block,'.Rda'))
    saveRDS(hr2_shapleys,paste0('output/hr2_shapleys_',i_block,'.Rda'))
  }
}