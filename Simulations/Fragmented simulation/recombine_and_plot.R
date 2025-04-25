#Recombine separate files

p_vec=c(10,30,50)
nrep=100
block_size=5



all_bart_dfs=all_ruleshap_dfs=all_hr1_dfs=all_hr2_dfs=list()
for(i in 1:(nrep/block_size)){
  all_ruleshap_dfs=c(all_ruleshap_dfs,list(readRDS(paste0('output/shapleys_dfs_',i,'.Rda'))))
  all_bart_dfs=c(all_bart_dfs,list(readRDS(paste0('output/shapleys_bart_dfs_',i,'.Rda'))))
  all_hr1_dfs=c(all_hr1_dfs,list(readRDS(paste0('output/hr1_shapleys_',i,'.Rda'))))
  all_hr2_dfs=c(all_hr2_dfs,list(readRDS(paste0('output/hr2_shapleys_',i,'.Rda'))))
}



bart_dfs=ruleshap_dfs=hr1_dfs=hr2_dfs=list()
file_index=rep(rep(1:20,each=block_size),times=length(p_vec))


for(i in 1:300){
  bart_dfs[[i]]=all_bart_dfs[[file_index[i]]][[i]]
  ruleshap_dfs[[i]]=all_ruleshap_dfs[[file_index[i]]][[i]]
  hr1_dfs[[i]]=all_hr1_dfs[[file_index[i]]][[i]]
  hr2_dfs[[i]]=all_hr2_dfs[[file_index[i]]][[i]]
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
SigCounts/500
#Type I error:
NoiseCounts/rep(100*c(5,25,45),each=4)
#FDR:
NoiseCounts/(SigCounts+NoiseCounts)
