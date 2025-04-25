# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'





### Using this as an archive

mywRandomForest_rules=function(x,y,r,weights,mtry=ncol(x)/3,ntree=500,maxdepth=3){
  #Define the custom loss, which will depend on the weights
  custom_loss <- function(preds, dtrain, weights) {
    labels <- getinfo(dtrain, "label")
    grad <- weights*(preds - labels)  # Gradient (first derivative)
    hess <- weights  # Hessian (second derivative)
    return(list(grad = grad, hess = hess))
  }

  #Set parameters to turn XGB into a RF
  params=list('eta'=1, 'gamma'=0, 'max_depth'=maxdepth,# 'min_child_weight'=0,
              'colsample_bynode'=mtry/ncol(x), 'lambda'=0, 'alpha'=0,
              'tree_method'='exact', 'subsample'=1)

  #Fit initial forest to determine p_hat
  first_fit=randomForest(x=x, y=as.factor(y))
  p_hat=predict(first_fit,type='prob')[,2]

  #Fit trees one at the time with xgb
  n=nrow(x)
  rules=NULL
  for(i in 1:ntree){
    #Define bootstrap subsample to fit the forest on
    indices=sample(1:n,n,replace=T)
    X_tree=x[indices,]
    weights_tree=weights[indices]
    i_tree=rbinom(n,size=1,prob=p_hat[indices])+1
    y_tree=c(t(r))[2*(indices-1)+i_tree]

    #Fit single tree on it
    inner_fit=xgboost(data = as.matrix(X_tree),
                      label = y_tree,
                      objective = function(preds, dtrain){custom_loss(preds, dtrain,weights_tree)},
                      nrounds = 1,verbose=0,params=params)

    #Extract rules (if it fails, the tree is contant, so just ignore it)
    tryCatch(expr = {
      unified_XGB = XGB2unified(inner_fit)
      rules=c(rules,
              unified2rules(unified_XGB,maxdepth = maxdepth, disaggregate = T))
    },
    error = function(e){})
  }

  return(rules)
}


wRandomForest_rules=function(x,y,weights,mtry=ncol(x)/3,ntree=500,maxdepth=3){
  #Define the custom loss, which will depend on the weights
  custom_loss <- function(preds, dtrain, weights) {
    labels <- getinfo(dtrain, "label")
    grad <- weights*(preds - labels)  # Gradient (first derivative)
    hess <- weights  # Hessian (second derivative)
    return(list(grad = grad, hess = hess))
  }

  #Set parameters to turn XGB into a RF
  params=list('eta'=1, 'gamma'=0, 'max_depth'=maxdepth, 'min_child_weight'=5,
              'colsample_bynode'=mtry/ncol(x), 'lambda'=0, 'alpha'=0,
              'tree_method'='exact','subsample'=1)


  #Fit trees one at the time with xgb
  n=nrow(x)
  rules=NULL
  for(i in 1:ntree){
    #Define bootstrap subsample to fit the forest on
    indices=sample(1:n,n,replace=T)
    X_tree=x[indices,]
    y_tree=y[indices]

    #Fit single tree on it
    inner_fit=xgboost(data = as.matrix(X_tree),
                      label = y_tree,
                      objective = function(preds, dtrain){custom_loss(preds, dtrain,weights[indices])},
                      nrounds = 1, verbose=0,params=params)

    #Extract rules (flat trees give error, ignore those)
    tryCatch(expr = {
      unified_XGB = XGB2unified(inner_fit)
      rules=c(rules,
              unified2rules(unified_XGB,maxdepth = maxdepth, disaggregate = T))
    },
    error = function(e){})
  }

  return(rules)
}


myRandomForest_rules_logit=function(x,y,r, ntree=500, maxdepth = 3){

  #Fit initial forest to determine p_hat
  first_fit=randomForest(x=x, y=as.factor(y))
  p_hat=predict(first_fit,type='prob')[,2]

  #Fit trees one at the time with RF
  n=nrow(x)
  rules=NULL
  for(i in 1:ntree){
    #Define bootstrap subsample to fit the forest on
    indices=sample(1:n,n,replace=T)
    X_tree=x[indices,]
    i_tree=rbinom(n,size=1,prob=p_hat[indices])+1
    y_tree=c(t(r))[2*(indices-1)+y[indices]+1]#i_tree]

    inner_fit=randomForest(X_tree, y=y_tree, ntree=1, maxnodes=2^maxdepth,
                           replace=F,sampsize = n)

    unified_RF = RF2unified(inner_fit, as.data.frame(X_tree))
    rules=c(rules,
            unified2rules(unified_RF,maxdepth = maxdepth, disaggregate = T))
  }

  return(rules)
}

