library(Matrix)
jaccard <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  
  ## only non-zero values of common
  Aim = A[im]
  
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
  
  return( J )
}

getMeans=function(v) {
  samples=v@sim$samples[colnames(v@sim$samples) != "lp__"]
  means=apply(samples,2,mean)
  idx_wo_lp <- which(v@model_pars != "lp__")
  skeleton <- rstan:::create_skeleton(v@model_pars[idx_wo_lp], v@par_dims[idx_wo_lp])
  rstan:::rstan_relist(means, skeleton)
}

unscale=function(y,ys) {
  atty=attributes(ys)
  y*atty[["scaled:scale"]]+atty[["scaled:center"]]
}

require(dplyr)

cor_se=function(x,y) {
  if (length(x) <= 2) return(1)
  ct=cor.test(x,y)
  ct$estimate / ct$statistic # se
}

get_score=function(pred, test) {
  test$my_pred=pred
  df=test %>% group_by(COMBINATION_ID) %>% summarise(pearson=cor(SYNERGY_SCORE, my_pred), corv=cor_se(SYNERGY_SCORE,my_pred)^2, n=length(SYNERGY_SCORE)) %>% as.data.frame
  df$w=sqrt(df$n-1) * (df$n > 2)
  score=sum( df$w * df$pearson ) / sum(df$w)
  attr(score,"df")=df
  attr(score,"se")=sqrt( sum( df$w^2 * df$corv ) / sum(df$w)^2 )
  score
}