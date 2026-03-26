# this is a custom statFun NEST function for computing delta adjusted rsq
#' statFun.gam.deltaRsq() will be called by NEST_audrey() if NEST argument statFun=="gam.deltaRsq"
#'@param X n x p matrix of image measurements (n = number of subjects, p = number of image locations)
#'@param dat matrix including phenotype and any covariates. make sure all columns in dat are names
#'@param gam.formula user-specified formula for gam fit at each vertex (in each gam, outcome will be a different column of X and the right side of the formula should include variables included in dat matrix)
#'@param lm.formula user-specified formula for lm fit at each vertex (purpose is to get the sign of the delta adjusted R-squared)
#'@param y.in.gam character specifying how the phenotype variable can be identified in gam. if phenotype is included as a linear term in the gam formula, this could just be the variable name. if there's a smooth term for the phenotype, it may be something like s(y)
#'@param y.in.lm character specifying how the phenotype variable can be identified in lm. again, might just be the name of the variable, but could be different (e.g., if testing an interaction term)
#'@param y.permute specify names of columns in dat matrix that should be permuted when getNull = TRUE
#'@param n.cores for parallelization, number of cores to use (default is 1 / no parallelization)
#'@param seed optional to set seed
#'@param n.perm number of permutations to conduct for inference. default is 999 (i.e. minimum p-value 1/1000)
#'@param getNull whether to obtain null distribution vs. just get observed map of statistics. default will be TRUE inside NEST function, but the statFun function will then be recursively called to get null distribution and getNull will then switch to FALSE
#'@export
statFun.gam.deltaRsq = function(X, dat, gam.full.formula, gam.null.formula, lm.formula, y.in.gam, y.in.lm, y.permute = NULL, n.cores = 1, seed = NULL, n.perm = 999, getNull = TRUE){
  
  xyz = cbind(X.v = X[,1], # just use the first image location for building model.matrix (will end up refitting at every location though)
              dat)
  
  # pre-specify the permutations (i.e., do the same permutations of people across vertices)
  if (getNull == TRUE){
    if (is.null(y.permute)){
      message("need to specify which columns of `dat` should be permuted")
      return(NULL)
    }
    if (!is.null(seed)){set.seed(seed)}
    perm.ind = lapply(1:n.perm, FUN = function(k){
      sample(1:nrow(X), replace = FALSE)
    })
  }
  
  # Calculate delta adjusted R-squared (ΔR^2adj) for observed data:
  stat.obs = unlist(mclapply(1:ncol(X), FUN = function(v){
    xyz[,"X.v"] = X[,v]
    
    # Full model (with smooth term)
    gam.fullmodel = gam(formula = gam.full.formula, data = data.frame(xyz))
    gam.fullmodel.results = summary(gam.fullmodel)
    
    # Null model (without smooth term)

    #nullmodel.formula = update(gam.full.formula, . ~ . - s(age, k = 3, fx = TRUE))  # Remove the smooth term
    gam.nullmodel = gam(formula = gam.null.formula, data = data.frame(xyz))
    gam.nullmodel.results = summary(gam.nullmodel)
    
    # Calculate delta adjusted R-squared (ΔR^2adj)
    adjRsq = gam.fullmodel.results$r.sq - gam.nullmodel.results$r.sq
    
    # Determine the sign using the linear model
    lm.model = lm(lm.formula, data = data.frame(xyz))
    # lm.model.t = summary(lm.model)$coefficients[2, 3]  # t-value for the first covariate (assumed to be the smooth term)
    
    # if(lm.model.t < 0){
    #   adjRsq = adjRsq * -1
    # }

    sign.t = sign(coef(lm.model)[y.in.lm])
    adjRsq = sign.t * adjRsq

    return(adjRsq)
  }, mc.cores = n.cores))
  
  # If getNull is TRUE, calculate the null distribution of ΔR^2adj
  if (getNull == TRUE){
    stat.null = lapply(1:n.perm, FUN = function(k){
      
      dat.k = dat
      dat.k[,y.permute] = dat.k[perm.ind[[k]], y.permute]
      
      statFun.gam.deltaRsq(X = X,
                           dat = dat.k,
                           gam.full.formula = gam.full.formula,
                           gam.null.formula = gam.null.formula,
                           lm.formula = lm.formula,
                           y.in.gam = y.in.gam,
                           y.in.lm = y.in.lm, n.cores = n.cores, seed = seed,
                           n.perm = NULL,  # no permutations to do with getNull = FALSE
                           getNull = FALSE)
      
    })
    return(list(T.obs = stat.obs,
                T.null = stat.null))
  } else{
    return(stat.obs)
  }
}