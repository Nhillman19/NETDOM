# helper functions for statFun() function
getCmat = function(full.mod.columns, columns.to.test){

  Cmat = matrix(0, nrow = length(columns.to.test),
                ncol = length(full.mod.columns),
                dimnames = list(columns.to.test,
                                full.mod.columns))
  diag(Cmat[,columns.to.test]) = 1

  return(Cmat)

}

# H0: C\beta = 0
## constraint matrix C has dimension q x p, where each row corresponds to one constraint (i.e. one constraint for each beta we are testing to be = 0)
mv_wald_stat = function(columns.to.test, vcov.v, beta.v){
  Cmat = getCmat(full.mod.columns = colnames(vcov.v), columns.to.test = columns.to.test)
  to.invert = Cmat %*% vcov.v %*% t(Cmat)
  middle.term = solve(to.invert)
  Cbeta = Cmat %*% beta.v
  stat = t(Cbeta) %*% middle.term %*% Cbeta
  return(stat)
}

get_vcov = function(X.v, model.matrix.full, MTM.inv, beta.v){
  resid.v = (X.v - model.matrix.full%*%beta.v)
  vcov.out = MTM.inv*sum(resid.v^2)/(nrow(model.matrix.full) - length(beta.v)) # n = number of subjects, p = number of parameters in full model
  return(vcov.out)
}

get_odc <- function(x1,x2,t_trim){
  t <- seq(t_trim,1,length = length(x1))
  cons <- sqrt(length(x1)*length(x2)/(length(x1) + length(x2)))
  x2_cdf <- ecdf(x2)
  x2_comp_x1_quant <- cons*(x2_cdf(quantile(x1,probs = t)) - t)
  int_val <- trapz(t,x2_comp_x1_quant)
  return(int_val)
  }
