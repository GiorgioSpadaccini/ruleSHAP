test_that("Shapleys computed correctly", {
  #Define input arguments
  data=gendata.friedman1(n=100,p=5)
  formula=y~.
  ntree=3
  burn.in=2
  nmc=2

  #Fit the model
  fitted_model=ruleSHAP(formula,data,ntree=3,burn.in=2,nmc=2,verbose = F)

  #Compute Shapleys
  shaps=expect_no_error(compute_SHAP(fitted_model,data[,1:5],intercept=T,
                     alpha=0.05,interactions=F,block_size=1e2))
})
