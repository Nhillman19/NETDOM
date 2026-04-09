#' NETDOM function
#'@importFrom parallel mclapply
#'@param statFun A string specifying the statistical method to use. Options include:
#'   \describe{
#'     \item{"lm"}{Linear regression model (`statFun.lm`).}
#'     \item{"gam.mvwald"}{Generalized additive model with multivariate Wald tests (`statFun.gam.mvwald`).}
#'     \item{"gam.deltaRsq"}{Generalized additive model for delta adjusted R-squared (`statFun.gam.deltaRsq`).}
#'     \item{"custom"}{Custom statistical function provided via `statFun.custom`.}
#'   }
#'@param args A list of arguments required by the selected `statFun`. Specific arguments depend on the chosen method:
#'   \describe{
#'     \item{"lm"}{Requires `X` (design matrix), `y` (response variable). Optional: `Z` (covariates), `type`, `FL`, `getNull`, `n.perm`.}
#'     \item{"lm.fast"}{Requires `X` (design matrix), `y` (response variable). Optional: `Z` (covariates), `type`, `FL`, `getNull`, `n.perm`.}
#'     \item{"gam.mvwald"}{Requires `X`, `dat`, `gam.formula`, `lm.formula`, `y.in.gam`, `y.in.lm`.}
#'     \item{"gam.deltaRsq"}{Requires `X`, `dat`, `gam.full.formula`, `gam.null.formula`, `lm.formula`, `y.in.gam`, `y.in.lm`.}
#'     \item{"custom"}{Depends on the custom function implementation and should match the arguments specified in `statFun.custom`.}
#'   }
#'@param net.maps A list of vectors representing network maps. Each vector contains 0s (outside network/ROI) and 1s (inside network/ROI). The length of each vector must match the number of columns in `X`.
#'@param direction A string indicating the direction of testing. Options are `"right"` for testing positive associations and `"left"` for testing negative associations.
#'@param one.sided Logical. If `TRUE`, performs one-sided testing. Default is `TRUE`.
#'@param n.cores Integer. The number of cores to use for parallel processing. Default is `1`.
#'@param seed Integer. Random seed for reproducibility. Default is `NULL`.
#'@param what.to.return A character vector specifying the outputs to return. Options include:
#'   \describe{
#'     \item{"pval"}{P-values for each network (default).}
#'     \item{"ES"}{Enrichment scores (observed and null distributions).}
#'     \item{"ES.obs"}{Observed enrichment scores.}
#'     \item{"T.obs"}{Observed test statistics.}
#'     \item{"T.null"}{Null test statistics.}
#'     \item{"everything"}{All of the above outputs.}
#'   }
#'@param statFun.custom A user-defined function for custom statistical analysis. Required if `statFun` is set to `"custom"`. This function should accept arguments passed via `args`.
#'@return A list containing outputs specified in `what.to.return`. Default is p-values for each network.
#'@export
#' 
NETDOM = function(statFun, args, net.maps, direction, one.sided = TRUE, n.cores = 1, seed = NULL, what.to.return = c("pval"), statFun.custom=NULL){

  # check input requirements:
  if (!is.list(net.maps)){
    message("net.maps argument should be a list, of vectors (one for each network)")
    return(NULL)
  }

  if (ncol(args$X) != length(net.maps[[1]])){
    message("each network should be specified as a vector with length equal to number of columns of X.")
    return(NULL)
  }

  if (!identical(as.numeric(sort(unique(unlist(net.maps)))), as.numeric(c(0,1)))){ # check that network maps are all 0s and 1s
    message("net.maps should include only 0's and 1's (1 = in network/ROI; 0 = outside network/ROI)")
    return(NULL)
  }

  # map brain-phenotype associations
  if (statFun == "lm" | statFun == "lm.fast"){
    required.args = c("X", "y")
    optional.args = c("Z", "type", "FL", "getNull","n.perm")

    args = checkArgs(args = args, required.args = required.args, optional.args = optional.args)
    if (!isFALSE(args)) {
      if (statFun == "lm") {
        statFun.out <- statFun.lm(X = args$X,
                               y = args$y,
                               Z = args$Z,
                               type = args$type,
                               n.cores = n.cores,
                               seed = seed,
                               FL = args$FL,
                               n.perm = args$n.perm,
                               getNull = args$getNull)
      } else {
        statFun.out <-  statFun.lm.fast(X = args$X,
                               y = args$y,
                               Z = args$Z,
                               type = args$type,
                               n.cores = n.cores,
                               seed = seed,
                               FL = args$FL,
                               n.perm = args$n.perm,
                               getNull = args$getNull)
      }
    } else {
      message("fix args!")
      return(NULL)
    }
  } else if (statFun == "gam.mvwald"){
    # check args:
    required.args = c("X","dat","gam.formula","lm.formula","y.in.gam","y.in.lm")
    optional.args = c("n.perm")

    args = checkArgs(args = args, required.args = required.args, optional.args = optional.args)

    if (!isFALSE(args)){
      statFun.out = statFun.gam.mvwald(X = args$X,
                                      dat = args$dat,
                                      gam.formula = args$gam.formula,
                                      lm.formula = args$lm.formula,
                                      y.in.gam = args$y.in.gam,
                                      y.in.lm = args$y.in.lm,
                                      y.permute = args$y.permute,
                                      n.cores = n.cores, seed = seed,
                                      n.perm = args$n.perm,
                                      getNull = TRUE # if doing this inside NEST function, assume testing is being done (if just want to get map, could just use the statFun function directly)
                                      )
    }else{
      message("fix args!")
      return(NULL)
    }

  }

  else if (statFun == "gam.deltaRsq"){
    if (!isFALSE(args)){
      statFun.out = statFun.gam.deltaRsq(X = args$X,
                                      dat = args$dat,
                                      gam.full.formula = args$gam.full.formula,
                                      gam.null.formula = args$gam.null.formula,
                                      lm.formula = args$lm.formula,
                                      y.in.gam = args$y.in.gam,
                                      y.in.lm = args$y.in.lm,
                                      y.permute = args$y.permute,
                                      n.cores = n.cores, seed = seed,
                                      n.perm = args$n.perm,
                                      getNull = TRUE)
    }else{
      message("fix args!")
      return(NULL)
    }
  }

  # Get p-values for NETDOM for each network map

  pval.NETDOM.list = lapply(net.maps,FUN = function(net.map){
    t_trim <- seq(0,.95,by = .05)
    if(direction == "right") {
      in_obs_ord <- statFun.out$T.obs[net.map == 1] 
      out_obs_ord <- statFun.out$T.obs[net.map == 0] 
    } else if(direction == "left") {
      in_obs_ord <- -1*statFun.out$T.obs[net.map == 1] 
      out_obs_ord <- -1*statFun.out$T.obs[net.map == 0] 
    }

    out_cdf <- ecdf(out_obs_ord)
    t <- seq(0,1,length = length(out_obs_ord))
    x2_cdf <- ecdf(out_obs_ord)
    x2_comp_x1_quant <- x2_cdf(quantile(in_obs_ord,probs = t)) - t
    obs.odc <- unlist(lapply(t_trim,FUN = function(x) {
      indx <- (1+floor(x*length(t))):length(t)
      trapz(t[indx],x2_comp_x1_quant[indx])
    }))

    
    l_out <- sum(net.map == 0)
    t <- seq(0,1,length = l_out)
    if(direction == "left"){
      mult <- -1
    } else {
      mult <- 1
    }
    null.odc <- mclapply(1:args$n.perm,FUN = function(k) {
      stat_vec <- statFun.out$T.null[[k]]
      in_stat_ord <- mult*stat_vec[net.map == 1]
      out_stat_ord <- mult*stat_vec[net.map == 0]
      out_cdf <- ecdf(out_stat_ord)
      out_comp_in_quant <- out_cdf(quantile(in_stat_ord,probs = t)) - t
      odc_val_cum <- cumtrapz(t,out_comp_in_quant)
      odc_val <- odc_val_cum[l_out] - odc_val_cum[(1+floor(t_trim*l_out))]
      return(odc_val)
    },mc.cores = n.cores)
  
  null.odc.mat <- do.call("rbind",null.odc)
  null_odc_out <- apply(null.odc.mat,MARGIN = 1,FUN = function(r){
    pval_vec <- c()
    for(c in 1:ncol(null.odc.mat)){
      pval_vec[c] <- pvalFun(obs = r[c],null.dist = as.numeric(null.odc.mat[,c]))
    }
    min_pval <- min(pval_vec)
    gamma_opt_idx <- which.min(pval_vec)
    gamma_opt <- t_trim[gamma_opt_idx]
    return(c(gamma_opt,min_pval))
  })

  null_gamma_opts <- null_odc_out[1,]
  null_pvals_odc <- null_odc_out[2,]
  
  pvals_obs_odc <- c()
    for(c in 1:ncol(null.odc.mat)){
      pvals_obs_odc[c] <- pvalFun(obs = obs.odc[c],null.dist = as.numeric(null.odc.mat[,c]))
    }
    min_pval_obs_odc <- min(pvals_obs_odc)
    min_pval_idx <- which.min(pvals_obs_odc)
    gamma_opt <- t_trim[min_pval_idx]
    final_pval_odc <- pvalFun(-1*min_pval_obs_odc,-1*null_pvals_odc)
    return(list("Gamma_opt" = gamma_opt, "Null_gamma_opts" = null_gamma_opts, "Diff_pval_odc" = final_pval_odc))
  })

  # Get p-values for distance from 0

  coef.zero.list = mclapply(net.maps, FUN = function(net.map){
    T.in.mean = mean(statFun.out$T.obs[net.map == 1])

    T.in.null = unlist(lapply(1:args$n.perm, FUN = function(k) {  
          mean(statFun.out$T.null[[k]][net.map == 1])
    }))
    return(list(T.obs.mean = T.in.mean,
                T.null.mean = T.in.null))
}, mc.cores = n.cores)

  pval.zeroCoef.list = mclapply(1:length(net.maps), FUN = function(net){
    if (direction == "right") {
      pvalFun(obs = coef.zero.list[[net]]$T.obs.mean, null.dist = coef.zero.list[[net]]$T.null.mean)
    } else {
      pvalFun(obs = -1*coef.zero.list[[net]]$T.obs.mean, null.dist = -1*coef.zero.list[[net]]$T.null.mean)
    }
  }, mc.cores = n.cores) 

  out = list()

  if ("everything" %in% what.to.return){
    out$pval.NETDOC.list = pval.NETDOC.list
    out$pval.zeroCoef.list = pval.zeroCoef.list
    out$pval.RIGEA.list = pval.RIGEA.list
    out$ES = ES.list
    out$T.obs = statFun.out$T.obs
    out$T.null = statFun.out$T.null
    return(out)
  } else{
    if ("pval" %in% what.to.return){ # default is to just return p-values
      out$pval.NETDOM.list = pval.NETDOM.list
      out$pval.zeroCoef.list = pval.zeroCoef.list
    }

    if ("T.obs" %in% what.to.return){
      out$T.obs = statFun.out$T.obs
    }

    if ("T.null" %in% what.to.return){
      out$T.null = statFun.out$T.null
    }
    return(out)
  }
}
