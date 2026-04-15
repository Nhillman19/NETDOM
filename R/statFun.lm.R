#' statFun.lm() will be called by NETDOM() if NETDOM argument statFun=="lm"
#'@param X n x p matrix (n = number of subjects, p = number of image locations)
#'@param y vector length n with phenotype measurements for each subject
#'@param Z vector (length n) or matrix (number of rows = n) with covariates to be adjusted for in lm. default Z = 1 assumes no covariates (i.e. Z=1 becomes placeholder for the intercept)
#'@param type type of statistic from lm() output to use. default is coefficient (will be the coefficient associated with phenotype variable y)
#'@param n.cores for parallelization, number of cores to use (default is 1 / no parallelization)
#'@param seed optional to set seed
#'@param FL whether to use Freedman-Lane approach to permutation. default is FALSE. this option may be preferable if phenotype/covariates are not independent
#'@param n.perm number of permutations to conduct for inference. default is 999 (i.e. minimum p-value 1/1000)
#'@param getNull whether to obtain null distribution vs. just get observed map of statistics. default will be TRUE inside NEST function, but the statFun function will then be recursively called to get null distribution and getNull will then switch to FALSE
#'@export
statFun.lm = function(X, y, Z = 1, type = "coef", n.cores = 1, seed = NULL, FL = FALSE, n.perm = 999, getNull = TRUE){

  # determine what column from summary(lm()) output will need to be extracted to match 'type' (if none of the options below are specified, will default to coef)
  type.ind = ifelse(type=="coef", 1, # summary(lm(X[,v] ~ y + Z))$coefficients["y",1] = coefficient
                    ifelse(type == "tvalue", 3, # summary(lm(X[,v] ~ y + Z))$coefficients["y",3] = t-value
                           ifelse(type == "pvalue", 4, # summary(lm(X[,v] ~ y + Z))$coefficients["y",4] = p-value
                                  1))) # revert to coef if unspecified

  # pre-specify the permutations (i.e. do the same permutations of people across vertices)
  if (getNull == TRUE){
    if (!is.null(seed)){ set.seed(seed)}
    perm.ind = lapply(1:n.perm, FUN = function(k){
      sample(1:nrow(X), replace = F)
    })
  }

  if (FL == FALSE){ # no freedman-lane (assume independence between y and Z)
    stat.obs = unlist(mclapply(1:ncol(X), FUN = function(v){
      # observed statistic:
      if(identical(Z,1)){
        fullmod = lm(X[,v] ~ y)
      } else {
        fullmod = lm(X[,v] ~ y + Z)
      }
      
      obs.stat.v = summary(fullmod)$coefficients["y",type.ind]

      return(obs.stat.v)

    }, mc.cores = n.cores))

    if (getNull == TRUE){ # recursive function used below -- specifying getNull = FALSE so that we don't get a null distribution for the null iterations
      null.stat = lapply(1:n.perm, FUN = function(k){
        statFun.lm(X = X, y = y[perm.ind[[k]]], Z = Z, type = type, n.cores = n.cores, seed = seed, FL = FALSE, getNull = FALSE)$T.obs
      })

      stat.null = do.call("rbind",null.stat)
      return(list(T.obs = stat.obs,
                  T.null = stat.null))
    } else{
      return(list(T.obs = stat.obs))
    }

  }  else{ # do freedman-lane
    if (getNull == TRUE){
      FL.out = mclapply(1:ncol(X), FUN = function(v){
        # observed statistic:
        if(identical(Z,1)){
          fullmod = lm(X[,v] ~ y)
        } else {
          fullmod = lm(X[,v] ~ y + Z)
        }
        obs.stat.v = summary(fullmod)$coefficients["y",type.ind]
        # under H0:
        if(identical(Z,1)){
          reducedmod.v = lm(X[,v] ~ 1)
        } else {
          reducedmod.v = lm(X[,v] ~ Z)
        }

        reducedmod.resid.v = matrix(resid(reducedmod.v), ncol = 1)

        null.stat.v = sapply(1:n.perm, FUN = function(k){
          permX.v.k = reducedmod.resid.v[perm.ind[[k]],] + reducedmod.v$fitted.values # permute residuals + add back fitted values from reduced model
          # Should this be permX.v.k? Because permX.k is not defined and could cause an error.
          if(identical(Z,1)){
            fullmod.permX.v.k = lm(permX.v.k ~ y)
          } else {
            fullmod.permX.v.k = lm(permX.v.k ~ y + Z)
          }
           # refit full model using permuted data
          #fullmod.permX.v.k = lm(permX.k ~ y + Z) # refit full model using permuted data
          null.stat.v.k = summary(fullmod.permX.v.k)$coefficients["y",type.ind]
          return(null.stat.v.k)
        })

        return(list(obs = obs.stat.v, null = null.stat.v))

      }, mc.cores = n.cores)

      # observed T(v):
      stat.obs = sapply(FL.out, FUN = function(x){x$obs})

      # null T(v) (i.e., maps of null T(v) which will be used to get null distribution of enrichment scores):
      stat.null = do.call("cbind",lapply(FL.out, FUN = function(x){x$null}))

      return(list(T.obs = stat.obs,
                  T.null = stat.null))

    } else{
      message("specify FL = FALSE if getNull = FALSE")
      break
    }

  }
}
