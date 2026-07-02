burn.in=2e3
nmc=2e4
maxdepth=3
rounding.digits=2
mu=0.5
eta=2
n_vec=c(500,1e3,3e3)
nrep=100
formula = y ~ .
ntree=500
n_max=max(n_vec)
p=10



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
friedman_trains=gendata.friedman1(n=nrep*n_max,p=p,sd=10,rho=0.3+diag(.7,nrow=p,ncol=p))
friedman_test=gendata.friedman1(n=1e4,p=p,sd=10,rho=0.3+diag(.7,nrow=p,ncol=p))

ruleshap_dfs=bart_dfs=hr_dfs=list()


#Fit models
set.seed(42)
for(i in 1:nrep){
  for(n_index in 1:length(n_vec)){
    
    #Set repetition-specific quantities
    set.seed(length(n_vec)*(i-1)+n_index)
    n=n_vec[n_index]
    data = friedman_trains[(i-1)*n_max+(1:n),]
    
    #Verbose
    print(paste('Fitting on sample size n =',n))
    
    
    #Start off with BART
    bartFit <- dbarts::bart(data[,1:p], data$y, keeptrees = TRUE, verbose = FALSE)
    bart_dfs[[(n_index-1)*nrep+i]]=bart_df=bartSHAP(bartFit,data[,1:p])
    
    #Now Horserule
    hr_data=as.data.frame(model.matrix(formula,data))[,-1]
    hr_data[,as.character(formula[[2]])]=data[,as.character(formula[[2]])]
    
    hr <- horserule:::HorseRuleFit(model.formula = formula, data=hr_data,
                                    restricted=0,linterms = 1:(ncol(hr_data)-1),
                                    niter = nmc+burn.in, burnin = burn.in, thin = 1,
                                    ntree = ntree,intercept=T,
                                    linp=1, ensemble='both',mix=0.3)
    
    
    #Shapleys for HR
    hr_new=list(BetaSamples=t(hr$postdraws$beta[(1:burn.in),])/c(1,hr$modelstuff$sdl,hr$modelstuff$sdr)*hr$modelstuff$sdy,
                 rules=gsub("]"," ",gsub("X[,","x.",hr$rules,fixed=T),fixed=T),wins.points=matrix(ncol=p,nrow=2))
    shapleys=compute_SHAP(hr_new,data[,1:p])
    hr_dfs[[(n_index-1)*nrep+i]]=hr_df=shapleys$marginal
    
    
    #Now RuleSHAP
    RS_fit=ruleSHAP(formula=formula, data=data,family='gaussian',verbose=T,burn.in = burn.in,nmc = nmc)
    ruleshap_dfs[[(n_index-1)*nrep+i]]=ruleshap_df=compute_SHAP(RS_fit,data[,1:p])$marginal
    
    #Save results
    saveRDS(bart_dfs,'output/bart_dfs.Rda')
    saveRDS(hr_dfs,'output/hr_dfs.Rda')
    saveRDS(ruleshap_dfs,'output/ruleshap_dfs.Rda')
  }
}


#Plot results
lowsignal_quantile=0.10
SmallSigCounts=LargeSigCounts=NoiseCounts=list(matrix(0,ncol=100,nrow=3),
                                               matrix(0,ncol=100,nrow=3),
                                               matrix(0,ncol=100,nrow=3))
rownames(LargeSigCounts[[1]])=rownames(LargeSigCounts[[2]])=
  rownames(LargeSigCounts[[3]])=rownames(SmallSigCounts[[1]])=
  rownames(SmallSigCounts[[2]])=
  rownames(SmallSigCounts[[3]])=rownames(NoiseCounts[[1]])=
  rownames(NoiseCounts[[2]])=rownames(NoiseCounts[[3]])=
  c('BART','HorseRule','RuleSHAP')

names(SmallSigCounts)=names(LargeSigCounts)=
  names(NoiseCounts)=paste0('n=',c(500,1000,3000))
nrep=100
n_vec=rep(1:3,each=nrep)
ns=c(500,1000,3000)

for(i in 1:length(n_vec)){
  
  n_index=n_vec[i]
  n=ns[n_index]
  
  bart_df=bart_dfs[[i]]
  ruleshap_df=ruleshap_dfs[[i]]
  hr_df=hr_dfs[[i]]
  
  #Assess local significance
  bart_df$sig=sign(bart_df$CIinf*bart_df$CIsup)==1
  ruleshap_df$sig=sign(ruleshap_df$CIinf*ruleshap_df$CIsup)==1
  hr_df$sig=sign(hr_df$CIinf*hr_df$CIsup)==1
  
  #Asses size of signal
  target=bart_df[,c('x','value','predictor')]
  target$value=0
  target$value[1:n]=
    (10*sapply(target$x[1:n],\(x){integrate(\(z){sin(x*z*pi)},0,1)$value})+
       10*sin(target$x[1:n]*target$x[n+(1:n)]*pi))/2
  target$value[1:n]=target$value[1:n]-mean(target$value[1:n])
  target$value[n+(1:n)]=
    (10*sapply(target$x[n+(1:n)],\(x){integrate(\(z){sin(x*z*pi)},0,1)$value})+
       10*sin(target$x[1:n]*target$x[n+(1:n)]*pi))/2
  target$value[n+(1:n)]=target$value[n+(1:n)]-mean(target$value[n+(1:n)])
  target$value[2*n+(1:n)]=20*(target$x[2*n+(1:n)]-0.5)^2
  target$value[2*n+(1:n)]=target$value[2*n+(1:n)]-mean(target$value[2*n+(1:n)])
  target$value[3*n+(1:n)]=10*(target$x[3*n+(1:n)]-0.5)
  target$value[4*n+(1:n)]=5*(target$x[4*n+(1:n)]-0.5)
  
  target = target %>% 
    group_by(predictor) %>%
    mutate(thresh = quantile(abs(value),lowsignal_quantile)) %>%
    ungroup() %>%
    as.data.frame()
  
  small_signal=which(abs(target$value[1:(5*n)])<target$thresh[1:(5*n)])
  
  #Create variable determining signal/noise
  bart_df$signal=2*(bart_df$predictor %in% paste0('x.',1:5))
  ruleshap_df$signal=2*(ruleshap_df$predictor %in% paste0('x.',1:5))
  hr_df$signal=2*(hr_df$predictor %in% paste0('x.',1:5))
  
  bart_df$signal[small_signal]=ruleshap_df$signal[small_signal]=
    hr_df$signal[small_signal]=1
  
  bart_counts=bart_df %>% group_by(signal) %>% summarize(npoints=sum(sig))
  ruleshap_counts=ruleshap_df %>% group_by(signal) %>% summarize(npoints=sum(sig))
  hr_counts=hr_df %>% group_by(signal) %>% summarize(npoints=sum(sig))
  
  
  SmallSigCounts[[n_index]][1,i%%100+1]=
    sum(bart_counts$npoints*(bart_counts$signal==1))
  SmallSigCounts[[n_index]][2,i%%100+1]=
    sum(hr_counts$npoints*(hr_counts$signal==1))
  SmallSigCounts[[n_index]][3,i%%100+1]=
    sum(ruleshap_counts$npoints*(ruleshap_counts$signal==1))
  
  LargeSigCounts[[n_index]][1,i%%100+1]=
    sum(bart_counts$npoints*(bart_counts$signal==2))
  LargeSigCounts[[n_index]][2,i%%100+1]=
    sum(hr_counts$npoints*(hr_counts$signal==2))
  LargeSigCounts[[n_index]][3,i%%100+1]=
    sum(ruleshap_counts$npoints*(ruleshap_counts$signal==2))
  
  
  NoiseCounts[[n_index]][1,i%%100+1]=
    sum(bart_counts$npoints*(bart_counts$signal==0))
  NoiseCounts[[n_index]][2,i%%100+1]=
    sum(hr_counts$npoints*(hr_counts$signal==0))
  NoiseCounts[[n_index]][3,i%%100+1]=
    sum(ruleshap_counts$npoints*(ruleshap_counts$signal==0))
}


#Ridgeline density plots for local rejection rates
library(ggplot2)
library(ggridges)
lines_df=data.frame(npoints=c(t(NoiseCounts[[1]])/2500,
                              t(NoiseCounts[[2]])/5000,
                              t(NoiseCounts[[3]])/15000,
                              t(SmallSigCounts[[1]])/(lowsignal_quantile*2500),
                              t(SmallSigCounts[[2]])/(lowsignal_quantile*5000),
                              t(SmallSigCounts[[3]])/(lowsignal_quantile*15000),
                              t(LargeSigCounts[[1]])/((1-lowsignal_quantile)*2500),
                              t(LargeSigCounts[[2]])/((1-lowsignal_quantile)*5000),
                              t(LargeSigCounts[[3]])/((1-lowsignal_quantile)*15000)),
                    n=paste('n =',rep(c(500,1000,3000),each=100*3)),
                    model=rep(c('BART','HorseRule','RuleSHAP'),each=100),
                    order=rep(c(3,1,2),each=100),
                    signal=rep(c('Noise','Small signal','Large signal'),each=100*3*3))
lines_df$signal=factor(lines_df$signal, levels=c('Noise','Small signal','Large signal'))
lines_df = lines_df %>%
  group_by(n,model,signal) %>%
  mutate(avg_npoints=mean(npoints,na.rm = T)) %>%
  ungroup() %>%
  as.data.frame()


#For the sake of plotting visualization, set any NoiseCount larger than 6% to 6%
NoiseCounts[[1]]=pmin(NoiseCounts[[1]],0.06*2500)
NoiseCounts[[2]]=pmin(NoiseCounts[[2]],0.06*5000)
NoiseCounts[[3]]=pmin(NoiseCounts[[3]],0.06*15000)
#And now redefine npoints (this way avg_npoints is correct)
lines_df$npoints=c(t(NoiseCounts[[1]])/2500,
                   t(NoiseCounts[[2]])/5000,
                   t(NoiseCounts[[3]])/15000,
                   t(SmallSigCounts[[1]])/(lowsignal_quantile*2500),
                   t(SmallSigCounts[[2]])/(lowsignal_quantile*5000),
                   t(SmallSigCounts[[3]])/(lowsignal_quantile*15000),
                   t(LargeSigCounts[[1]])/((1-lowsignal_quantile)*2500),
                   t(LargeSigCounts[[2]])/((1-lowsignal_quantile)*5000),
                   t(LargeSigCounts[[3]])/((1-lowsignal_quantile)*15000))

models=c('OLS','LASSO','RuleSHAP','RuleFit','HorseRule','BART','RF','cTree')


ggplot(lines_df, aes(x = npoints,
                                             fill = factor(model,levels=models),
                                             color = factor(model,levels=models),
                                             y=factor(model,levels=models[c(6,5,3)]))) +
  #geom_segment(aes(x = avg_npoints, xend = avg_npoints,
  #                 y = factor(model,levels=models[c(6,5,3)]),
  #                 yend = factor(model,levels=models[c(6,5,3)])+0.05,
  #                 color = factor(model,levels=models[c(6,5,3)])))+
  geom_point(aes(x = avg_npoints,
                 y = factor(model,levels=models[c(6,5,3)]),
                 color = factor(model,levels=models[c(6,5,3)]),
                 #shape = factor(model,levels=models)
  ))+
  geom_density_ridges(alpha=0.45,
                      scale=2,rel_min_height=1e-2) +
  facet_grid(factor(n,levels=paste('n =',c(500,1000,3000)))~signal,scales='free')+
  theme_bw()+
  theme(strip.text.x = element_text(size = 10))+
  xlim(c(0,NA))+
  labs(x = "Fraction of observations with significant feature effect",
       y = "Frequency across replicates (density)")+
  scale_fill_manual(values=c(colorspace::rainbow_hcl(length(models))),
                    breaks=models,name='Model')+
  scale_color_manual(values=c(colorspace::rainbow_hcl(length(models))),
                     breaks=models,name='Model')+
  scale_shape_manual(values=1:length(models),
                     breaks=models, name='Model')+
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 1),
    shape = guide_legend(order = 1)
  )+
  theme(legend.position = "right", 
        legend.box = "vertical",
        legend.direction = "vertical",
        panel.spacing = unit(0.5, "lines"))
