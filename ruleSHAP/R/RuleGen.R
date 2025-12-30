#Other

#' @export
`%notin%` = Negate(`%in%`)




#' Compute matrices that identify the rules
#'
#' This function computes the matrices that are needed to compute marginal
#' Shapley values of a RuleSHAP model. These matrices check which points satisfy
#' each of the subrules of each rule. They also check which predictors are involved
#' in the definition of each rule.
#'
#' @param rules a character vector containing the rules to compute the matrices for
#' @param x_df a dataframe containing the observations to compute the matrices for
#'
#' @return RulePredMat A matrix with as many rows as there are rules and as many
#' columns as there are predictors. The \eqn{j}-th column has ones on the entries
#' corresponding to rules that involve the \eqn{j}-th predictor, while all
#' remaining entries are zeroes. Can be computed with the \link{RuleMats} function.
#' @return Rs a list of as many matrices as there are rules. The \eqn{j}-th element
#' is a matrix corresponding to the \eqn{j}-th rule. Each of its columns corresponds
#' to the 0-1 encoding of a subrule of the \eqn{j}-th rule, as observed in the data
#' provided with input parameter \code{x}.
#'
#' @export
RuleMats=function(rules,x_df){
  rules_separated=stringr::str_split_fixed(rules," & ",n=1+max(stringr::str_count(rules,'&')))
  vars_in_term <- gsub( " .*$", "",rules_separated)

  #Codify rules_separated into coding language
  Rulesmat=matrix(paste0("x_df$",rules_separated),ncol=ncol(rules_separated))
  Rulesmat[Rulesmat=="x_df$"]=""

  #Some subrules might involve the same predictor. Those, we'd still like together
  for(i in 1:nrow(vars_in_term)){
    for(x_name in unique(vars_in_term[i,])){
      #unique(vars_in_term[i,]) also includes the name "". Skip that case
      if(x_name==''){
        next
      }

      #Find all spots sharing the same predictor x_name
      indices=which(x_name==vars_in_term[i,])

      #Join back these subrules in rules_separated, place them in first occurrance
      rules_separated[i,indices[1]]=paste(rules_separated[i,indices],collapse=' & ')
      #In the remaining spots, we need to delete the subrule, it was already merged
      rules_separated[i,indices[-1]]=''

      #Do the same with Rulesmat
      Rulesmat[i,indices[1]]=paste(Rulesmat[i,indices],collapse=' & ')
      Rulesmat[i,indices[-1]]=''

      #Do the same with vars_in_term (nothing to collapse, only delete copies)
      vars_in_term[i,indices[-1]]=''
    }
  }
  #Define a matrix where M_{i,j} tells if x_j is involved in the i-th rule
  x_names <- names(x_df)
  id_mat <- t(apply(vars_in_term,MARGIN=1,FUN=function(x){x_names %in% x}))
  colnames(id_mat)=x_names

  #Rulesmat defined the rules in coding. Create a function that runs it
  parseval=function(text){return(eval(parse(text=text)))}

  #Use this function to create a list of matrices. Each matrix is \{R_i(x^{(j)}_i)\}_{i,j}
  Rs=list()
  for(i in 1:nrow(Rulesmat)){
    #Build matrix R for the i-th rule
    R=matrix(unlist(apply(Rulesmat[i,,drop=F],MARGIN=2,FUN = parseval)),
             nrow=nrow(x_df))
    colnames(R)=vars_in_term[i,vars_in_term[i,]!='']

    #Using vars_in_mat, reorder the matrices by predictor
    Rs[[i]]=R[,order(match(colnames(R),x_names)),drop=F]
  }

  return(list(RulePredMat=id_mat,Rs=Rs))
}




#' Generate rules with a Parametric Random Forest
#'
#' This function fits a Parametric Random Forest to generate dychotomous rules that
#' may be used to fit a RuleSHAP model. This is designed for a continuous outcome.
#'
#' @param x_df a dataframe containing the observations to use when computing the rules.
#' The dataframe should not contain the outcome \eqn{y}.
#' @param y a numerical vector representing the outcome \eqn{y}.
#' @param ntree number of trees that should be used when generating the rules
#' @param maxdepth maximum depth of the trees used to generate the rules. In particular,
#' this induces a maximum on the order of interactions between predictors.
#' @param disaggregate a logical value determining whether rules should be extracted
#' by also considering all possible subrules (\code{TRUE}) or if rules should only
#' be extracted in a top-down manner as in RuleFit and HorseRule (\code{FALSE}).
#' @param rounding.digits Number of decimals that the splitting points
#' are rounded to in the rule generation. This rounding simplifies readability
#' of rules and prevents rules from being very similar but not identical.
#' @param type Assumption underlying generation of synthetic copies of the outcome.
#' Possible values \code{"normal"} (generates errors from a normal distribution with 
#' $\sigma =$ random forest's OOB MSE) or \code{"empirical"} (errors are reshuffled 
#' OOB errors). 
#' @return A character vector containing all the rules that were generated. The rules
#' are not filtered for duplicates. To filter rules, check the function \code{filter_rules}.
#'
#' @export
PRF=function(x_df, y, ntree=500, maxdepth = 3, rounding.digits=3,
             disaggregate=TRUE, type = "normal"){

  #Fit initial forest to obtain OOB predictions (mu) and errors
  first_fit=randomForest(x=x_df, y=y)
  mu <- predict(first_fit)
  if (type == "normal"){
    MSE=mean((y-mu)^2) # OOB MSE
  } else if (type == "empirical"){
    errors=y-mu
  }

  #Fit small forests to get rules from
  n=nrow(x_df)
  rules=NULL
  for(i in 1:ntree){
    #Define bootstrap subsample to fit the forest on
    indices=sample(1:n,n,replace=TRUE)
    X=x_df[indices,]
    Y=mu[indices] + if (type == "normal") {
      rnorm(n, sd=sqrt(MSE))
    } else if (type == "empirical") {
      errors[sample(1:n)] ## samples errors from full sample
      ## alternatively, to sample errors from current bootstrap sample only:
      ## errors[indices][sample(1:n)]
    }

    inner_fit=randomForest::randomForest(X, y=Y, ntree=1, maxnodes=2^maxdepth,
                           replace=F,sampsize = n,nodesize=max(10,floor(0.025*n)))

    unified_RF = RF2unified(inner_fit, X,rounding.digits=rounding.digits)
    rules=c(rules,
            unified2rules(unified_RF,maxdepth = maxdepth, disaggregate = disaggregate))
  }

  return(rules)
}





#' Filter out duplicate rules
#'
#' This function filters out rules that are duplicates in the sense that their support
#' coincides on a given dataset. Rules can be given a priority when choosing which one to delete.
#'
#' @param rules A character vector containing the rules.
#' @param x_df a dataframe containing the observations to use when computing the rules.
#' The dataframe should not contain the outcome \eqn{y}.
#' @param filtering_tol By how many points can rules differ in terms of support and
#' still be considered duplicates. By default, \code{filtering_tol=0} implies that
#' only rules that coincide on each point of the dataset.
#' @param priority numerical vector of priorities: one value per rule. A higher value
#' means a more desirable rule. If left to NULL, the priority is given based on the order
#' of appearance in the vector.
#' @param verbose Logical value denoting whether a progress bar should be printed to show
#' the progress of the filtering.
#'
#' @return A character vector containing all the rules that remained after filtering.
#'
#' @export
filter_rules=function(rules,x_df,filtering_tol=0,priority=NULL,verbose=F){
  n=nrow(x_df)

  #Prepare for rule filtering: determine by how many points rules can be collinear
  if(filtering_tol < 1){ #then this is a proportion of points
    n_tol=n*filtering_tol
  }
  else if(filtering_tol %% 1 == 0){
    n_tol=filtering_tol
  }
  else{ #then it's in the wrong format
    stop('filtering_tol should be a proportion or a positive integer')
  }

  #Decide priority of rules through their ordering
  if(!is.null(priority)){
    rules=rules[order(priority, decreasing=T)]
  }

  #Filter out rules. Use the with() environment to keep variable names the same
  parseval=function(text){return(eval(parse(text=text)))}
  sampleXrule=sapply(paste0('with(as.data.frame(x_df),',rules,')'),parseval)

  #Do a pre-filtering: two rules will coincide or complement if sums do so too:
  #(adapted to accommodate tolerance)
  supports=colSums(sampleXrule)
  same_supports=sapply(supports,function(x){abs(supports-x)<=n_tol})
  compl_supports=sapply(supports,function(x){abs(supports-n+x)==n_tol})
  which_same=which(lower.tri(same_supports)&same_supports,arr.ind = T)
  which_compl=which(lower.tri(compl_supports)&compl_supports,arr.ind = T)

  #Set up progress bar
  pb_counter=0
  if(verbose){
    pb = utils::txtProgressBar(min = 0, max = nrow(which_same)+nrow(which_compl),
                        initial = 0,style=3)
  }

  #Find rules with same support, up to tolerance (loop over it to save memory)
  same_checks=c()
  for(i in seq(1,length.out=nrow(which_same))){
    if(sum(sampleXrule[,which_same[i,1]]==sampleXrule[,which_same[i,2]])>=n-n_tol){
      same_checks=c(same_checks,i)
    }

    if(verbose){
      pb_counter=pb_counter+1
      utils::setTxtProgressBar(pb,pb_counter)
    }
  }
  #Find rules with complementary support, up to tolerance (loop over it to save memory)
  compl_checks=c()
  for(i in seq(1,length.out=nrow(which_compl))){
    if(sum(sampleXrule[,which_compl[i,1]]==sampleXrule[,which_compl[i,2]])>=n-n_tol){
      compl_checks=c(compl_checks,i)
    }

    if(verbose){
      pb_counter=pb_counter+1
      utils::setTxtProgressBar(pb,pb_counter)
    }
  }


  #Since we have lower triangular, first entry is the later rule. Remove that one
  same_rules=which_same[same_checks,1]
  compl_rules=which_same[compl_checks,1]

  #Find rules that are almost constant
  const_rules=which(supports<=n_tol | supports>=n-n_tol)

  #Remove all rules, if there are any to remove
  if(length(c(compl_rules,same_rules,const_rules)) != 0){
    rules=rules[-c(compl_rules,same_rules,const_rules)]
  }

  return(rules)
}








#' Convert a unified dataframe to a vector of rules
#'
#' This function takes a tree ensemble as input and extracts all the rules that define it.
#'
#' @param unified A dataframe containing all essential information about the tree ensemble.
#' Can be produced with the \code{*.unify} functions from the \code{treeshap} package,
#' or from the \code{RF2unified} and \code{XGB2unified} functions in this package.
#' @param maxdepth maximum depth of the trees used to generate the rules. In particular,
#' this induces a maximum on the order of interactions between predictors.
#' @param disaggregate a logical value determining whether rules should be extracted
#' by also considering all possible subrules (\code{TRUE}) or if rules should only
#' be extracted in a top-down manner as in RuleFit and HorseRule (\code{FALSE}).
#'
#' @return A character vector containing all the rules extracted from the tree ensemble.
#'
#' @export
unified2rules=function(unified,maxdepth=3,disaggregate=F){
  Rs_df=unified

  #Add a column that gives an individual splitting rule
  Rs_df$rule=paste(Rs_df$Feature,Rs_df$Decision.type,Rs_df$Split)

  #To reconstruct rules, check which entries are the first rule of their tree
  deep_nodes=c(Rs_df$Yes,Rs_df$No)
  deep_nodes=deep_nodes[!is.na(deep_nodes)]
  lvl1=(1:nrow(Rs_df))[-deep_nodes]

  #Define function that gives the opposite of an inequality symbol, to negate the rule
  swapineq=function(x){
    x=gsub('%in%','sym0',x)
    x=gsub('%notin%','%in%',x)
    x=gsub('sym0','%notin%',x)
    x=gsub('<=','sym1',x)
    x=gsub('>=','sym2',x)
    x=gsub('<','sym3',x)
    x=gsub('>','<=',x)
    x=gsub('sym1','>',x)
    x=gsub('sym2','<',x)
    return(gsub('sym3','>=',x))
  }

  #Build recursive function to list all rules in a tree, given its root node
  rule_restructuring=function(i,rule,rules,counter){
    #i is the current index explored in the dataframe Rs_df
    #rule is the (whole) rule that brought the function to the current index
    #rules is a vector containing all the rules that have been extracted so far
    #counter counts the depth of the nesting of the functions

    #Also save linear combinations of nodes, to reduce how many shapleys to compute
    #Let's decide that right child is parent minus left child

    #If threshold is reached for depth, stop
    if(counter >= maxdepth){
      return(rules)
    }

    #Check if rule has a left and/or right subnode
    left_check=is.na(Rs_df$Yes[i])
    right_check=is.na(Rs_df$No[i])

    #Any inner node has left and right child nodes. Any terminal node has no subnodes
    if(left_check != right_check){
      stop('Inner node at row ',i,' of the dataframe does not have both left and right child.
           Dataframe might be wrongly formatted.')
    }

    #If it doesn't have children, then it's a terminal node, so no more recursion
    #If it does have child nodes, then do recursion and save the rules.

    if(!left_check){
      #Go down the left branch
      new_rule=paste(rule,'&',Rs_df$rule[i])
      rules=c(new_rule,
              rule_restructuring(Rs_df$Yes[i],new_rule,rules,counter+1))

      #Go down the right branch
      #If we're at depth one, the Left and Right are complementary. Don't save Right
      new_rule=paste(rule,'&',swapineq(Rs_df$rule[i]))
      if(rule==''){
        rules=rule_restructuring(Rs_df$No[i],new_rule,rules,counter+1)
      }
      else{
        rules=c(new_rule,
                rule_restructuring(Rs_df$No[i],new_rule,rules,counter+1))
      }

      #if disaggregate=T, on top of attaching left child or right child, add the
      #option to add nothing
      if(disaggregate){
        rules=c(rules,rule_restructuring(Rs_df$Yes[i],rule,NULL,counter+1),
                rule_restructuring(Rs_df$No[i],rule,NULL,counter+1))
      }
    }
    return(rules)
  }

  #For each root node, call the recursive function
  rules=NULL
  for(i in lvl1){
    #Call recursive function, starting already with the root rule
    restructured=rule_restructuring(i,'',NULL,0)
    rules=c(rules,substr(restructured,4,nchar(restructured)))
  }

  return(rules)
}







#' Convert Random Forest model to data frame format
#'
#' This function converts a Random Forest fit into a data frame format compatible with rule manipulation
#'
#' @param rf A Random Forest model as fitted from the \code{randomForest} package.
#' @param data a data frame containing the observations to use when computing the rules.
#' The data frame should not contain the outcome \eqn{y}.
#' @param rounding.digits Number of decimals that the splitting points
#' are rounded to in the rule generation. This rounding simplifies readability
#' of rules and prevents rules from being very similar but not identical.
#'
#' @return A data frame expressing the tree ensemble as a sequence of decision rules.
#'
#' @export
RF2unified=function(rf,data,rounding.digits=3){
  #Retrieve relevant info. For terminal nodes, most info is NA.
  var_names=names(rf$forest$xlevels)
  bestvar=rf$forest$bestvar
  bestvar[bestvar==0]=NA
  xbestsplit=rf$forest$xbestsplit
  xbestsplit[is.na(bestvar)]=NA
  leftDaughter=rf$forest$leftDaughter
  leftDaughter[is.na(bestvar)]=NA
  rightDaughter=rf$forest$rightDaughter
  rightDaughter[is.na(bestvar)]=NA
  nodepred=rf$forest$nodepred
  nodepred[!is.na(bestvar)]=NA
  ntree=rf$forest$ntree

  #Trees might have different lengths, meaning that not all values of matrix
  #should be read (every column only should have indices 1:numnodes, for that tree)
  #We create vector index_vec that specifies which entries are relevant
  times_vec=rep(0,2*ntree)
  times_vec[c(T,F)]=rf$forest$ndbigtree
  times_vec[c(F,T)]=nrow(rf$forest$leftDaughter)-rf$forest$ndbigtree
  index_vec=rep(rep(c(T,F),times=ntree),times=times_vec)

  #Define node_base, which summed up to the Daughters gives position in dataframe
  used_nodes=cumsum(rf$forest$ndbigtree)
  node_base=rep(c(0,used_nodes[-ntree]),times=rf$forest$ndbigtree)

  #Treat categorical predictors: first find which variables are non-numeric
  which_numeric=sapply(data,class)[var_names]=='numeric'

  #For categorical variables, the decision type is %in%
  destype_vec=ifelse(which_numeric,'<=','%in%')

  #Apply rounding. The categorical splits are represented with integers, so they remain unchanged
  #RF et simila take the middle point and put <=. So I take floor and keep the <=
  xbestsplit=floor(xbestsplit*(10^rounding.digits))/(10^rounding.digits)

  #For categorical variables, split corresponds to relevant categories (binary, see ?getTree)
  #Start from all splits and just edit the categorical entries
  split_vec=c(xbestsplit)[index_vec]

  #Find which splits need to be changed, which predictor it is and which levels it has
  #but only change them in the case that at least one of the splits is categorical
  if(sum(!(which_numeric[bestvar[index_vec]]),na.rm = T)!=0){
    cat_splits=which(!(which_numeric[bestvar[index_vec]]))
    cat_vars=names(cat_splits)
    all_levels=sapply(data,levels)

    #Apply changes
    for(i in 1:length(cat_splits)){
      levs=all_levels[[cat_vars[i]]]
      split_vec[cat_splits[i]]=
        paste0("c('",
               paste(levs[as.numeric(intToBits(split_vec[cat_splits[i]]))[1:length(levs)]==1],
                     collapse = "','"),"')")
    }
  }

  #Return data frame
  return(data.frame(Tree=rep(1:ntree,times=rf$forest$ndbigtree),#Node=1:rf$forest$ndbigtree,
                    Feature=var_names[bestvar[index_vec]],
                    Decision.type=destype_vec[bestvar[index_vec]],
                    Split=split_vec,
                    Yes=node_base+c(leftDaughter)[index_vec],
                    No=node_base+c(rightDaughter)[index_vec],
                    Prediction=c(nodepred)[index_vec]))
}







#' Convert XGBoost model to data frame format
#'
#' This function converts a Random Forest fit into a data frame format compatible with rule manipulation
#'
#' @param model Ax XGBoost model as fitted from the \code{xgboost} package.
#' @param data a data frame containing the observations to use when computing the rules.
#' The data frame should not contain the outcome \eqn{y}.
#' @param rounding.digits Number of decimals that the splitting points
#' are rounded to in the rule generation. This rounding simplifies readability
#' of rules and prevents rules from being very similar but not identical.
#'
#' @return A data frame expressing the tree ensemble as a sequence of decision rules.
#'
#' @export
#Function that converts XGB to unified
XGB2unified=function(model,data,rounding.digits=3){
  #Extract rom XGB
  unified=xgboost::xgb.model.dt.tree(model = model)
  unified$Feature[is.na(unified$Split)]=NA

  #Get the node base, i.e. how many nodes used before this tree (nodes start from 0, add +1)
  tree_sizes=(unified %>% group_by(Tree) %>% summarize(size=max(Node)+1))$size
  used_nodes=cumsum(tree_sizes)
  node_base=rep(c(0,used_nodes[-length(used_nodes)]),times=tree_sizes)+1

  #Use the node nr to add to node base
  unified[!is.na(unified$Split),c('Yes','No')]=
    unified[!is.na(unified$Split),c('Yes','No')] %>%
    lapply(X=,FUN=\(x){gsub(pattern=".*-",replacement="",x)}) %>%
    as.data.frame()
  unified[,c('Yes','No')]=as.data.frame(lapply(unified[,c('Yes','No')],as.numeric))
  unified[,c('Yes','No')]=unified[,c('Yes','No')]+node_base
  unified$Decision.type='<='
  unified$Decision.type[is.na(unified$Split)]=NA
  unified[,5]=floor(unified[,5]*(10^rounding.digits))/(10^rounding.digits)
  unified=unified[,c('Tree','Node','Feature','Decision.type','Split','Yes','No')]

  return(unified)
}

