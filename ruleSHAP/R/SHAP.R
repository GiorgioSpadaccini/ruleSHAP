#' Compute Shapley values for BART model
#'
#' This function computes the marginal Shapley values for a BART model as fitted with dbarts.
#'
#' @param model An object of class \code{dbartsSampler} as obtained by fitting a BART with the \code{dbarts} function
#' @param data A data frame containing the data used to compute the shapley values. Should not contain the outcome \eqn{y}
#' @param n_test While using the whole dataset to compute the shapley values, one might only
#' want to compute the Shapley values of the first few observations. Such number can be regulated with \code{n_test}
#' @param ping An update on the progress of the fit will be given every \code{ping} iterations.
#' To receive no updates at all, set this to a higher value than the number of MCMC samples of the model.
#'
#' @return A dataframe containing the shapley values of each of the first \code{n_test} points of the dataset,
#' per predictor, with confidence intervals (\code{CIinf} being the lower bound and \code{CIsup} being the upper bound)
#' @export
bartSHAP=function(model,data,n_test=nrow(data),ping=1e3){
  #Basic info
  p=ncol(data)
  n=nrow(data)

  #Get dataframe
  bartdf=dbarts::extract(model,type='trees')

  #Preallocate matrix storing the shapleys
  SHAPmat=matrix(ncol=max(bartdf$sample),nrow=n_test*p)

  #Get shapleys
  for(j in 1:max(bartdf$sample)){
    if(j %% ping ==0){
      print(paste('Working on MCMC sample',j))
    }
    #Get subdataframe
    df2un=bartdf[bartdf$sample == j,-1]

    #Convert it into unified format for treeSHAP

    #Adjust tree numbering to start from 0
    df2un$tree=df2un$tree-1

    #Add prediction and splitting variables and decision type and splitting feature
    df2un$Prediction=NA
    df2un$Prediction[df2un$var == -1] = df2un$value[df2un$var == -1]
    df2un$Split=NA
    df2un$Split[df2un$var != -1] = df2un$value[df2un$var != -1]
    df2un$Decision.type=factor(x=rep('<=',nrow(df2un)),levels=c('<=','<'))
    df2un$Decision.type[df2un$var == -1] = NA
    df2un$Feature=NA
    df2un$Feature[df2un$var != -1] = paste0('x.',df2un$var[df2un$var != -1])

    #Add Node
    #Count nodes per tree
    NodeXTree=aggregate(x=rep(1,nrow(df2un)),by=list(df2un$tree),FUN = sum)
    NodeXTree=NodeXTree$x[order(NodeXTree$Group.1)]
    #Repeat sequence of 1:numnodes for each tree, but start it from 0
    df2un$Node=c(unlist(apply(matrix(NodeXTree,nrow=1),2,function(n){0:(n-1)})))

    #Add Yes, No
    df2un$Yes=df2un$No=NA
    #Create a vector to keep track of which children have already been identified
    #Also keep track of which entries of the dataframe we're working on
    children=df2un$var==-1
    indices=1:nrow(df2un)
    while(length(children) > max(df2un$tree)+1){ #Stop when each tree has no children nodes left
      i=1
      while(i <= length(children)-2){
        #Check if parent is followed by its to children (they need to be in same tree)
        if((!children[i]) & children[i+1] & children[i+2] &
           df2un$tree[indices[i]] == df2un$tree[indices[i+1]] &
           df2un$tree[indices[i]] == df2un$tree[indices[i+2]]){

          #Update Yes and No
          df2un$Yes[indices[i]]=indices[i+1]
          df2un$No[indices[i]]=indices[i+2]

          #Remove the children, they've been dealt with. Turn parent into child
          indices=indices[-c(i+1,i+2)]
          children=children[-c(i+1,i+2)]
          children[i]=T
        }
        i=i+1
      }
    }

    #Assimilate Missing as Yes
    df2un$Missing=df2un$Yes

    #Put everything together in the unified object
    names(df2un)[1:2]=c('Tree','Cover')

    #Create unified object
    bart_unif <- list(model=
                        df2un[,c('Tree','Node','Feature','Decision.type','Split','Yes','No','Missing','Prediction','Cover')],
                      data=data, feature_names=names(data))
    attr(bart_unif,'missing_support')=T
    attr(bart_unif,'model')='bart'
    class(bart_unif) <- "model_unified"


    #Compute shapleys
    shapleys=treeshap:::treeshap(bart_unif,  data[1:n_test,], verbose = 0)

    #Store shapleys in the matrix
    SHAPmat[,j]=c(as.matrix(shapleys$shaps))
  }

  #Get CIs and point estimate, plot, same as usual
  vals=rowMeans(SHAPmat)
  CIs=apply(SHAPmat,MARGIN=1,FUN=function(x){quantile(x,probs=c(0.025,0.975))})


  #For some reason, BART rescales estimated outcomes for extract()
  #Retrieve by what they have been re-scaled and update Shapley values
  #accordingly, to put them back on the original scale

  #Use a linear model to retrieve the re-scaling. e.g. 5 observations are enough
  yhat=colMeans(predict(bartFit,data[1:5,]))
  y_shap=rowSums(matrix(vals,nrow=1e3))[1:5]
  rescale_fac=coef(lm(yhat ~ ., data=data.frame(yhat=yhat,y_shap=y_shap)))[2]
  vals=vals*rescale_fac
  CIs=CIs*rescale_fac


  #Define and return output
  shapleys_df=data.frame(x=c(unlist(lapply(data[1:n_test,],as.numeric))),
                         value=c(vals),
                         CIinf=c(CIs[1,]),
                         CIsup=c(CIs[2,]),
                         predictor=rep(names(data),each=n_test))

  return(shapleys_df)
}




#' Compute Shapley values for RuleSHAP model
#'
#' This function computes the marginal (interaction) Shapley values for a RuleSHAP model as fitted with the
#' \code{ruleSHAP} function. The computation includes confidence intervals.
#'
#' @param model An object of class \code{dbartsSampler} as obtained by fitting a BART with the \code{dbarts} function
#' @param data A data frame containing the data used to compute the Shapley values. Should not contain the outcome \eqn{y}
#' @param intercept Logical value denoting whether the fitted modeluses an intercept
#' @param alpha The significance level of the confidence intervals.
#' @param interactions A logical value denoting whether marginal interaction Shapley values should be computed,
#' on top of the overall values.
#' @param combine A list containing all combinations of Shapley values that should be additively
#' combined to explore the joint effect of two variables. For instance, if the joint effect
#' of variables \code{x.1} and \code{x.2} should be added to the dataframe, the list will
#' contain a vector \code{c('x.1','x.2')}.
#' @param block_size The computation of these Shapley values involves computations on large matrices which,
#' for computational reasons, are split into small submatrices. This integer denotes the size of a submatrix
#' in terms of how many rows it should contain.
#'
#' @return A dataframe containing the shapley values of each point of the dataset,
#' per predictor, with confidence intervals (\code{CIinf} being the lower bound and \code{CIsup} being the upper bound)
#' @export
compute_SHAP=function(model,data,intercept=T,alpha=0.05,interactions=F,
                      combine=NULL,block_size=5e3){
  #Winsorize data (store unwinsorized version too)
  data_unwins=data
  for(i in 1:ncol(data)){
    if(!any(is.na(model$wins.points[,i]))){
      data[data[,i]<model$wins.points[1,i],i]=model$wins.points[1,i]
      data[data[,i]>model$wins.points[2,i],i]=model$wins.points[2,i]
    }
  }


  #Compute basic quantities
  p=ncol(data)
  n_test=nrow(data)
  rule_objects=RuleMats(model$rules,data_unwins)
  WeightMatrices=ShapleyMats(data=data,Rs=rule_objects$Rs,
                             id_mat=rule_objects$RulePredMat,
                             interactions=interactions)
  if(intercept){CoeffMat=model$BetaSamples[-1,]}else{CoeffMat=model$BetaSamples}

  #If combinations are requested, add them to WeightMatrix of marginal
  #To do that, do left-matrix-multiplication by a sparse matrix.
  if(!is.null(combine)){
    #Compute the matrix that groups up weights
    combi_lengths=sapply(combine,length)
    combi_i=rep(1:(n_test*length(combine)),times=rep(combi_lengths,each=n_test))
    combi_j=sapply(combine,function(combination){
      rep((which(names(data) %in% combination)-1)*n_test,times=n_test)+
        rep(1:n_test,each=length(combination))
    })
    combi_mat=sparseMatrix(i=combi_i,j=combi_j,x=1,
                           dims=c(length(combine)*n_test,nrow(WeightMatrices$marginal)))

    #Left multiply and append to WeightMatrices
    WeightMatrices$marginal=rbind(WeightMatrices$marginal,
                                  combi_mat%*%WeightMatrices$marginal)
  }


  #First do main effects
  CIs=matrix(nrow=2,ncol=nrow(WeightMatrices$marginal))
  vals=rep(NA,nrow(WeightMatrices$marginal))
  for(i in 1:ceiling(nrow(WeightMatrices$marginal)/block_size)){
    block_i=(block_size*(i-1)+1):min(block_size*i,nrow(WeightMatrices$marginal))
    SHAP_vecs=WeightMatrices$marginal[block_i,,drop=F]%*%CoeffMat
    CIs[,block_i]=apply(SHAP_vecs,1,function(x){quantile(x,probs=c(alpha/2,1-alpha/2))})
    vals[block_i]=rowMeans(as.matrix(SHAP_vecs))
  }

  shapleys_df=data.frame(x=c(unlist(lapply(data_unwins,as.numeric)),
                             rep(NA,nrow(data)*length(combine))),
                         value=c(vals),
                         CIinf=c(CIs[1,]),
                         CIsup=c(CIs[2,]),
                         predictor=rep(c(names(data),
                                         unlist(sapply(combine,paste,collapse=' & '))),
                                       each=nrow(data)))

  #Now do interactions, if requested
  shapleys_inter_df=NULL
  if(interactions){
    CIs=matrix(nrow=2,ncol=nrow(WeightMatrices$interactions))
    vals=rep(NA,nrow(WeightMatrices$interactions))
    for(i in 1:ceiling(nrow(WeightMatrices$interactions)/block_size)){
      block_i=(block_size*(i-1)+1):min(block_size*i,nrow(WeightMatrices$interactions))
      SHAP_vecs=WeightMatrices$interactions[block_i,,drop=F]%*%CoeffMat
      CIs[,block_i]=apply(SHAP_vecs,1,function(x){quantile(x,probs=c(alpha/2,1-alpha/2))})
      vals[block_i]=rowMeans(as.matrix(SHAP_vecs))
    }

    shapleys_inter_df=data.frame(x1=c(unlist(lapply(data_unwins,as.numeric))),
                                 x2=c(unlist(lapply(data_unwins,function(x){rep(as.numeric(x),times=p)}))),
                                 value=c(vals),
                                 CIinf=c(CIs[1,]),
                                 CIsup=c(CIs[2,]),
                                 predictor1=rep(names(data),each=n_test),
                                 predictor2=rep(names(data),each=n_test*p))
  }

  return(list(marginal=shapleys_df,interactions=shapleys_inter_df))
}











#' Compute a matrix with Shapley values of each rule
#'
#' This function computes the Shapley values of each rule of a RuleSHAP model and
#' arranges them in a matrix that can be used to compute Shapley values of the
#' model as a whole.
#'
#' @param data a dataframe containing the data used to estimate the expectations in Shapley values
#' @param data_test a dataframe containing the points to compute the Shapley values of. By default,
#' this coincides with the datapoints used to estimate Shapley values.
#' @param Rs a list of as many matrices as there are rules to compute the Shapley values for.
#' The \eqn{j}-th element is a matrix corresponding to the \eqn{j}-th rule. Each of its columns
#' corresponds to the 0-1 encoding of a subrule of the \eqn{j}-th rule, as observed in the data
#' provided with input parameter \code{data}. Can be computed with the \link{RuleMats} function.
#' @param Rs_test same as Rs, but computed for the (possibly different) observations provided
#' from the \code{data_test} parameter.
#' @param RulePredMat A matrix with as many rows as there are rules and as many columns as there
#' are predictors. The \eqn{j}-th column has ones on the entries corresponding to rules that
#' involve the \eqn{j}-th predictor, while all remaining entries are zeroes. Can be computed
#' with the \link{RuleMats} function.
#' @param interactions A logical parameter determining whether interaction Shapley values should
#' also be computed
#'
#' @return marginal A matrix with \eqn{n\cdot p} rows and as many columns as
#' there are terms (both linear and rules). It is obtained by vertically
#' stacking matrices of \eqn{n} rows. Each submatrix focuses on the shapley values
#' of a different predictor: the \eqn{(i,k)}-th entry of the \eqn{j}-th of
#' such submatrices represents the contribution of the \eqn{k}-th term to
#' the Shapley value of the \eqn{i}-th datapoint for the \eqn{j}-th predictor.
#' The first \eqn{p} terms are the linear terms, and the remaining columns
#' refer to the rules.
#' @return interaction A matrix with \eqn{n\cdot p^2} rows and as many columns as
#' there are terms (both linear and rules). It is obtained by vertically
#' stacking matrices of \eqn{n \cdot p} rows. Each submatrix is in turn split
#' into \eqn{p} subsubmatrices which focuses on the interaction shapley values
#' of a different pair of predictors: the \eqn{(i,k)}-th entry of the \eqn{j}-th
#' subsubmatrix of the \eqn{j'}-th submatrix represents the contribution of the
#' \eqn{k}-th term to the Shapley value of the \eqn{i}-th datapoint for the
#' interaction between the \eqn{j}-th and the \eqn{j'}-th predictor.
#' The first \eqn{p} terms are the linear terms, and the remaining columns
#' refer to the rules.
#'
#' @export
ShapleyMats=function(data,data_test=data,Rs,Rs_test=Rs,id_mat, interactions=F){
  #Obtain basic quantities
  n=nrow(data)
  n_test=nrow(data_test)
  X=model.matrix(~ ., data)[,-1] #remove intercept, model.matrix adds it and doing ~ .-1 does not fix it
  X_test=model.matrix(~ ., data_test)[,-1]
  p=ncol(data)
  pX=ncol(X)
  q=nrow(id_mat)
  P=pX+q

  #compute how many columns each predictor takes up
  ndummies=unlist(lapply(data,nlevels))-1
  ndummies[ndummies < 0]=1

  #To avoid the computational cost of repeatedly updating sparse Matrices, i.e.
  #their dataframes, store everything as preallocated vectors i,j,x (SparseDataFrame).
  #in marginal/main effects, each rule produces n_test*d non-zero coefficients
  #in interactions case, each rule produces n_test*d*(d-1) non-zero coefficients
  #for each of these, use a pointer to keep track of the first entry available
  #for writing (all i,j,x of the same type share same pointer)
  d_vec=sapply(Rs,ncol)
  i_marg=j_marg=x_marg=integer(sum(d_vec)*n_test)
  i_inter=j_inter=x_inter=integer(sum(d_vec*(d_vec-1))*n_test)
  first_wrt_marg=first_wrt_inter=1

  #The first p columns are for the linear terms:
  i_marg[1:(n_test*pX)]=(rep(1:p,times=ndummies*n_test)-1)*n_test+1:n_test
  j_marg[1:(n_test*pX)]=rep(1:pX,each=n_test)
  x_marg[1:(n_test*pX)]=c(X_test)-rep(colMeans(X),each=n_test)
  first_wrt_marg=n_test*pX+1


  #The remaining ones are for the rules:
  #create a progress bar
  pb = utils::txtProgressBar(min = 0, max = length(Rs), initial = 0,style=3)

  #Go through all the rule contributions one by one
  for(i in 1:length(Rs)){
    d=d_vec[i]
    #If d=1, then the contribution is like for linear terms (and interactions=0)
    if(d==1){
      #Involved predictor
      inv_pred=which(id_mat[i,]==1)

      #Update matrix
      i_marg[seq(first_wrt_marg,length.out=n_test)]=
        (inv_pred-1)*n_test+1:n_test
      j_marg[seq(first_wrt_marg,length.out=n_test)]=pX+i
      x_marg[seq(first_wrt_marg,length.out=n_test)]=
        c(Rs_test[[i]])-mean(Rs[[i]])

      #Update pointer
      first_wrt_marg=first_wrt_marg+n_test

      #Update progress bar
      utils::setTxtProgressBar(pb,i)
      next
    }

    #Otherwise, id d>1, we compute the Shapleys properly
    SSt=tcrossprod(!Rs_test[[i]],!Rs[[i]])
    SSt_checks=which(SSt==0,arr.ind = T)
    SSt_points=data.frame(x=SSt_checks[,1],y=SSt_checks[,2])

    #Now go through all involved predictors and update the contributions
    inv_pred=which(id_mat[i,]==1)
    for(j in 1:d){
      #Start with updating the marginal shapleys

      #check the extra condition R_i(x_i)R_i(y_i)=0. These points are contributing
      contributions=SSt_points[(Rs_test[[i]][SSt_points$x,j]*Rs[[i]][SSt_points$y,j])==0,]

      #For these points, compute the weights. Arrange them in a matrix
      Weights=matrix(0,nrow=n_test,ncol=n)
      qx=rowSums(Rs_test[[i]][,-j,drop=F])
      qy=rowSums(Rs[[i]][,-j,drop=F])
      Weights[as.matrix(contributions)]=
        (Rs_test[[i]][contributions$x,j]-Rs[[i]][contributions$y,j])/
        (n*(d-qx[contributions$x])*
           choose(2*d-qx[contributions$x]-qy[contributions$y]-1,d-qx[contributions$x]))

      sum_contributions=rowSums(Weights)

      #Add the weights to the matrix
      #Update matrix
      i_marg[seq(first_wrt_marg,length.out=n_test)]=
        (inv_pred[j]-1)*n_test+1:n_test
      j_marg[seq(first_wrt_marg,length.out=n_test)]=pX+i
      x_marg[seq(first_wrt_marg,length.out=n_test)]=
        sum_contributions
      #Update pointer
      first_wrt_marg=first_wrt_marg+n_test


      #Now update the interaction shapleys, if so requested
      if(interactions){
        #(to avoid double counting, only do predictors before j)
        for(k in seq(1,length.out=j-1)){
          #Compute interaction-specific q(x)
          qx_int=qx-Rs_test[[i]][,k]
          qy_int=qy-Rs[[i]][,k]

          #Compute interaction specific weights (you can use contributions df, since
          #R_j(x_j)=R_j(y_j) implies no contribution to interactions)
          Weights=matrix(0,nrow=n_test,ncol=n)
          Weights[as.matrix(contributions)]=
            (Rs_test[[i]][contributions$x,j]*Rs_test[[i]][contributions$x,k]
             +Rs[[i]][contributions$y,j]*Rs[[i]][contributions$y,k]
             -Rs[[i]][contributions$y,j]*Rs_test[[i]][contributions$x,k]
             -Rs_test[[i]][contributions$x,j]*Rs[[i]][contributions$y,k])/
            (2*n*(d-1-qx_int[contributions$x])*
               choose(2*d-qx_int[contributions$x]-qy_int[contributions$y]-3,d-1-qx_int[contributions$x]))

          sum_contributions=rowSums(Weights)

          #Add the weights to the interaction matrix (both as j,k and k,j
          #since we wrote the combination only as k<j)
          #Update matrix
          i_inter[seq(first_wrt_inter,length.out=2*n_test)]=
            c((inv_pred[j]-1)*n_test*p+(inv_pred[k]-1)*n_test+1:n_test,
              (inv_pred[k]-1)*n_test*p+(inv_pred[j]-1)*n_test+1:n_test)
          j_inter[seq(first_wrt_inter,length.out=2*n_test)]=pX+i
          x_inter[seq(first_wrt_inter,length.out=2*n_test)]=
            sum_contributions #this whole vector will be repeated twice
          #Update pointer
          first_wrt_inter=first_wrt_inter+2*n_test
        }
      }

    }

    #Update progress bar
    utils::setTxtProgressBar(pb,i)
  }

  #close progress bar
  close(pb)

  #Using i,j,x computed above, build the sparse matrices
  #(if contributions$x does not contain all points, then i,j,x are shorter
  #on the off chance that this happens, shorten the vectors to prevent error)

  Shapleymat=Matrix::sparseMatrix(i=i_marg[seq(1,length.out=first_wrt_marg-1)],
                                  j=j_marg[seq(1,length.out=first_wrt_marg-1)],
                                  x=x_marg[seq(1,length.out=first_wrt_marg-1)],
                                  dims=c(n_test*p,P))

  InterShapleymat=Matrix::sparseMatrix(i=i_inter[seq(1,length.out=first_wrt_inter-1)],
                                       j=j_inter[seq(1,length.out=first_wrt_inter-1)],
                                       x=x_inter[seq(1,length.out=first_wrt_inter-1)],
                                       dims=c(n_test*p^2,P))

  #The interaction Shapley matrix is missing the main effects. Compute them by
  #deducting the interactions from the marginal shapleys
  MainShapleymat=InterShapleymat
  dim(MainShapleymat)=c(n_test,p^2*P)
  #Prepare matrix to perform block sum
  blocksum=Matrix::sparseMatrix(i=1:(p*P),j=1:(p*P),x=1) %x% matrix(1,nrow=p,ncol=1)
  #Sum up all interactions per predictor, per point
  MainShapleymat=MainShapleymat %*%blocksum
  #Put them back in the same format as marginal shapley values
  dim(MainShapleymat)=c(n_test*p,P)
  #Subtract the sums of interactions from the marginal shapleys. That's the main effect
  MainShapleymat = Shapleymat-MainShapleymat

  #Return the matrix
  return(list(marginal=Shapleymat,main=MainShapleymat,
              interactions=InterShapleymat))
}
