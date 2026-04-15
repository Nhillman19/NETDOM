#' statFun.lm.fast() will be called by NETDOM() if NETDOM argument statFun=="lm.fast"
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
statFun.lm.fast <- function(X, y, Z = 1, type = "coef", n.cores = 1, seed = NULL, FL = FALSE, n.perm = 999, getNull = TRUE){

# determine what column from summary(lm()) output will need to be extracted to match 'type' (if none of the options below are specified, will default to coef)
type_ind <- ifelse(type=="coef", 1, # summary(lm(X[,v] ~ y + Z))$coefficients["y",1] = coefficient
                ifelse(type == "tvalue", 3, # summary(lm(X[,v] ~ y + Z))$coefficients["y",3] = t-value
                        ifelse(type == "pvalue", 4, # summary(lm(X[,v] ~ y + Z))$coefficients["y",4] = p-value
                                1))) # revert to coef if unspecified

# pre-specify the permutations (i.e. do the same permutations of people across vertices)
if (getNull == TRUE){
  if (!is.null(seed)){ set.seed(seed)}
  perm_ind <- lapply(1:n.perm, FUN = function(k){
    sample(1:nrow(X), replace = FALSE)
  })
}

  if (FL == FALSE){ # no freedman-lane (assume independence between y and Z)
    if (is.matrix(Z)) {
        mat <- cbind(1,Z,y)
        mat_red <- cbind(1,Z)
      } else if (identical(Z,1)) {
        mat <- cbind(1,y)
        mat_red <- mat[,1]
    }
  MTM_inv <- solve(t(mat) %*% mat)
  deg_f <- nrow(mat) - ncol(mat)
  bs <- (MTM_inv %*% t(mat) %*% X)
  resids <- X - (mat%*%bs)
  sig2 <- colSums(resids^2)/deg_f
  n_col_mat <- ncol(mat)
  ts <- bs["y",]/sqrt(sig2 * MTM_inv[n_col_mat,n_col_mat])
  ps <- pt(ts,df = deg_f)
  ps_both <- ifelse(ps > 1 - ps,1-ps,ps)*2
  if (type_ind == 1) {
    stat_obs <- bs
  } else if (type_ind == 3){
    stat_obs <- ts
  } else if (type_ind == 4){
    stat_obs <- ps_both
  }
  # Recursive function used below 
  if (getNull == TRUE) {
    stat_null <- do.call("rbind",mclapply(1:n.perm, FUN = function(k){
    statFun.lm.fast(X = X, y = y[perm_ind[[k]]], Z = Z, type = type, n.cores = n.cores, seed = seed, FL = FALSE, getNull = FALSE)$T.obs
  },mc.cores = n.cores))
    return(list(T.obs = stat_obs,T.null = stat_null))
    } else {
      return(list(T.obs = stat_obs))
    }
  } else { # do freedman-lane
    if (getNull == TRUE) {
      if (is.matrix(Z)) {
        mat <- cbind(1,Z,y)
        mat_red <- cbind(1,Z)
      } else if (identical(Z,1)) {
        mat <- cbind(1,y)
        mat_red <- mat[,1]
      }
      # Get residuals and fitted values from reduced model
      MTM_inv_red <- solve(t(mat_red)%*%mat_red)
      bs_red <- MTM_inv_red%*%t(mat_red)%*%X
      resids_red <- X - mat_red%*%bs_red
      # Get list of observed + null statistics (each element contains a vector) 
      stat <- mclapply(1:(n.perm+1), FUN = function(k){ 
        if (k == (n.perm + 1)) {
          Xperm <- X
        } else {
          # Permute residuals + add fitted values
          Xperm <- resids_red[perm_ind[[k]],] + mat_red%*%bs_red 
        }
        MTM_inv <- solve(t(mat)%*%mat)
        deg_f <- nrow(mat) - ncol(mat)
        bs <- (MTM_inv%*%t(mat)%*%Xperm)
        resids <- Xperm - (mat%*%bs)
        sig2 <- colSums(resids^2)/deg_f
        n_col_mat <- ncol(mat)
        ts <- bs["y",]/sqrt(sig2*MTM_inv[n_col_mat,n_col_mat])
        ps <- pt(ts,df = deg_f)
        ps_both <- ifelse(ps > 1 - ps,1-ps,ps)*2
        if(type_ind == 1){
            stat_ret <- bs
        } else if(type_ind == 3){
            stat_ret <- ts
        } else if(type_ind == 4){
            stat_ret <- ps_both
        }
        return(stat_ret)
      },mc.cores = n.cores)

      stat_obs <- stat[[(n.perm + 1)]]
      stat_null <- do.call("rbind",stat[1:n.perm])
      return(list(T.obs = stat_obs,T.null = stat_null))
    } else {
      message("specify FL = FALSE if getNull = FALSE")
      break
    }
  }
}
