############################################################################################
## Tetraploid likelihood, LOD and recombination frequency functions
## Assumes random bivalent pairing. 
## Peter Bourke, Wageningen UR Plant Breeding. August 2016
############################################################################################

#' Calculate recombination frequency, LOD and log-likelihood from frequency tables in a random pairing tetraploid
#' @description This group of functions is called by \code{\link{linkage}}.
#' @param x A frequency table of the different classes of dosages in the progeny. The column names start with \code{"n_"}. Followed by the dosage of the first marker and then of the second.
#' @param ncores Number of cores to use for parallel processing (deprecated).
#' @return
#' A list with the following items:
#' \item{r_mat}{A matrix with recombination frequencies for the different phases}
#' \item{LOD_mat}{A matrix with LOD scores for the different phases}
#' \item{logL_mat}{A matrix with log likelihood ratios for the different phases}
#' \item{phasing_strategy}{A character string specifying the phasing strategy. \code{"MLL"} for maximum likelihood en \code{"MINR"} for minimum recombination frequency.}
#' \item{possible_phases}{The phases between markers that are possible. Same order and length as column names of output matrices.}
#' @name r4_functions
NULL

#' @rdname r4_functions
#' @noRd
r4_1.0_1.0<-function(x, ncores=1){
  # I have disabled the disomic scripts as they are not appropriate here (these are functions for polysomic behaviour)

  r_c <- (x[,"n_01"] + x[,"n_10"])/(x[,"n_00"] + x[,"n_01"] + x[,"n_10"] + x[,"n_11"]) # recombination frequencies coupling
  r_r <- (2*x[,"n_00"] - x[,"n_01"] - x[,"n_10"] + 2*x[,"n_11"])/(x[,"n_00"] + x[,"n_01"] + x[,"n_10"] + x[,"n_11"]) # recombination frequencies repulsion   
  # r_rdis <- (x[,"n_00"] + x[,"n_11"])/(x[,"n_00"] + x[,"n_01"] + x[,"n_10"] + x[,"n_11"]) 
  
  LOD_c <- (x[,"n_00"] + x[,"n_11"])*log10(pmax(1 - r_c,1e-6)) + (x[,"n_01"] + x[,"n_10"])*log10(pmax(r_c,1e-6)) - (x[,"n_00"] + x[,"n_11"])*log10(0.5) - (x[,"n_01"] + x[,"n_10"])*log10(0.5)
  LOD_r <- (x[,"n_00"] + x[,"n_11"])*log10(pmax(1 - r_c,1e-6)) + (x[,"n_01"] + x[,"n_10"])*log10(pmax(r_c,1e-6)) - (x[,"n_00"] + x[,"n_11"])*log10(0.5) - (x[,"n_01"] + x[,"n_10"])*log10(0.5)
  # LOD_rdis <- (x[,"n_00"] + x[,"n_01"] + x[,"n_10"] + x[,"n_11"])*log10(2) + (x[,"n_00"] + x[,"n_11"])*log10(pmax(r_rdis,1e-6)) + (x[,"n_01"] + x[,"n_10"])*log10(pmax(1-r_rdis,1e-6))
  
  # logL_c <- 2^(-x[,"n_00"] - x[,"n_01"] - x[,"n_10"] - x[,"n_11"])*(1 - r_c)^(x[,"n_00"] + x[,"n_11"])*r_c^(x[,"n_01"] + x[,"n_10"])
  # logL_r <- 6^(-x[,"n_00"] - x[,"n_01"] - x[,"n_10"] - x[,"n_11"])*(2 - r_r)^(x[,"n_01"] + x[,"n_10"])*(1 + r_r)^(x[,"n_00"] + x[,"n_11"])
  
  logL_c <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_10"] - x[,"n_11"])*log(2) + 
    (x[,"n_01"] + x[,"n_10"])*log(r_c) + (x[,"n_00"] + x[,"n_11"])*log(1 - r_c)
  
  logL_r <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_10"] - x[,"n_11"])*(log(2) + log(3)) + 
    (x[,"n_01"] + x[,"n_10"])*log(2 - r_r) + (x[,"n_00"] + x[,"n_11"])*log(1 + r_r)
  
  # logL_rdis <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_10"] - x[,"n_11"])*log(2) + 
  #   (x[,"n_01"] + x[,"n_10"])*log(1 - r_rdis) + (x[,"n_00"] + x[,"n_11"])*log(r_rdis)
   
  return(list(r_mat=cbind(r_c, r_r),# r_rdis),
              LOD_mat=cbind(LOD_c, LOD_r),# LOD_rdis),
              logL_mat=cbind(logL_c, logL_r),# logL_rdis),
              phasing_strategy="MLL", 
              possible_phases=c("coupling",
                                "repulsion"
                                )))#,"repulsion_disomic")))
  
}

#' @rdname r4_functions
#' @noRd
r4_1.0_2.0<-function(x, ncores=1){
  
  #attach(x)
  
  common_sum <- x[,"n_00"] + x[,"n_02"] + x[,"n_10"] + x[,"n_12"]
  
  r_c <- (x[,"n_02"] + x[,"n_10"])/common_sum
  r_r <- (x[,"n_00"] + x[,"n_12"])/common_sum
  
  ############################################################################
  ## These are not the true likelihood functions, as I have removed the parts that are
  ## independent of r, which cancel in the likelihood ratio anyway...
  Lc <- function(r) {
    L <- 2^(-common_sum)*(1 - r)^(x[,"n_00"] + x[,"n_12"])*r^(x[,"n_02"] + x[,"n_10"])
    return(L)}
  
  Lr <- function(r) {
    L <- 2^(-common_sum)*(1 - r)^(x[,"n_02"] + x[,"n_10"])*r^(x[,"n_00"] + x[,"n_12"])
    return(L)}
  
  LOD_c <- log10(Lc(r_c)/Lc(0.5))
  LOD_r <- log10(Lr(r_r)/Lr(0.5))
  
  ## Record the logL also:
  logL_c <- 2^(-x[,"n_00"] - x[,"n_02"] - x[,"n_10"] - x[,"n_12"])*3^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"])*(1 - r_c)^(x[,"n_00"] + x[,"n_12"])*r_c^(x[,"n_02"] + x[,"n_10"])
  logL_r <- 2^(-x[,"n_00"] - x[,"n_02"] - x[,"n_10"] - x[,"n_12"])*3^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"])*(1 - r_r)^(x[,"n_02"] + x[,"n_10"])*r_r^(x[,"n_00"] + x[,"n_12"])
  
  #detach(x)
  return(list(r_mat=cbind(r_c, r_r),
              LOD_mat=cbind(LOD_c, LOD_r),
              logL_mat=cbind(logL_c, logL_r),
              phasing_strategy="MLL", 
              possible_phases=c("coupling",
                                "repulsion")))
}

#' @rdname r4_functions
#' @noRd
r4_1.0_1.1<-function(x, ncores=1){
  r_c <- (x[,"n_02"] + x[,"n_10"])/(x[,"n_00"]+x[,"n_02"]+x[,"n_10"]+x[,"n_12"])
  r_r <- (2*x[,"n_00"]-x[,"n_02"]-x[,"n_10"]+2*x[,"n_12"])/(x[,"n_00"]+x[,"n_02"]+x[,"n_10"]+x[,"n_12"])
  
  ############################################################################
  LOD_c <- log10(((1 - r_c)^(x[,"n_00"] + x[,"n_12"])*r_c^(x[,"n_02"] + x[,"n_10"]))/((1 - 0.5)^(x[,"n_00"] + x[,"n_12"])*0.5^(x[,"n_02"] + x[,"n_10"])))
  LOD_r <- log10(((2 - r_r)^(x[,"n_02"] + x[,"n_10"])*(1 + r_r)^(x[,"n_00"] + x[,"n_12"]))/((2 - 0.5)^(x[,"n_02"] + x[,"n_10"])*(1 + 0.5)^(x[,"n_00"] + x[,"n_12"])))
  
  ## Record the logL also:
  logL_c <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"])*log(4) + (x[,"n_00"] + x[,"n_12"])*log(1 - r_c) + (x[,"n_02"] + x[,"n_10"])*log(r_c)
  logL_r <- (-x[,"n_00"] - x[,"n_02"] - x[,"n_10"] - x[,"n_12"])*log(3) + (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"])*log(4) + (x[,"n_02"] + x[,"n_10"])*log(2 - r_r) + (x[,"n_00"] + x[,"n_12"])*log(1 + r_r)
  
  return(list(r_mat=cbind(r_c, r_r),
              LOD_mat=cbind(LOD_c, LOD_r),
              logL_mat=cbind(logL_c, logL_r),
              phasing_strategy="MLL", 
              possible_phases=c("coupling",
                                "repulsion")))
}

#' @rdname r4_functions
#' @noRd
r4_1.0_1.3<-function(x, ncores=1){
  r_c <- (x[,"n_03"] + x[,"n_11"])/(x[,"n_01"]+x[,"n_03"]+x[,"n_11"]+x[,"n_13"])
  r_r <- (2*x[,"n_01"]-x[,"n_03"]-x[,"n_11"]+2*x[,"n_13"])/(x[,"n_01"]+x[,"n_03"]+x[,"n_11"]+x[,"n_13"])
  
  ############################################################################
  LOD_c <- log10(((1 - r_c)^(x[,"n_01"] + x[,"n_13"])*r_c^(x[,"n_03"] + x[,"n_11"]))/((1 - 0.5)^(x[,"n_01"] + x[,"n_13"])*0.5^(x[,"n_03"] + x[,"n_11"])))
  LOD_r <- log10(((2 - r_r)^(x[,"n_03"] + x[,"n_11"])*(1 + r_r)^(x[,"n_01"] + x[,"n_13"]))/((2 - 0.5)^(x[,"n_03"] + x[,"n_11"])*(1 + 0.5)^(x[,"n_01"] + x[,"n_13"])))
  
  ## Record the logL also:
  logL_c <- (-x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"])*log(4) + (x[,"n_01"] + x[,"n_13"])*log(1 - r_c) + (x[,"n_03"] + x[,"n_11"])*log(r_c)
  logL_r <- (-x[,"n_01"] - x[,"n_03"] - x[,"n_11"] - x[,"n_13"])*log(3) + (-x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"])*log(4) + (x[,"n_03"] + x[,"n_11"])*log(2 - r_r) + (x[,"n_01"] + x[,"n_13"])*log(1 + r_r)
  
  return(list(r_mat=cbind(r_c, r_r),
              LOD_mat=cbind(LOD_c, LOD_r),
              logL_mat=cbind(logL_c, logL_r),
              phasing_strategy="MLL", 
              possible_phases=c("coupling",
                                "repulsion")))
}

#' @rdname r4_functions
#' @noRd
r4_1.0_2.1<-function(x, ncores=1){

  logLc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13)*log(12) + 
      (n00 + n13)*log(1 - r) + (n01 + n12)*log(3 - r) + (n03 + n10)*log(r) + (n02 + n11)*log(2 + r)
    return(L)}
  
  inter_logLc <- function(n00,n01,n02,n03,n10,n11,n12,n13){
    optimize(logLc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  logLr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13)*log(12) + 
      (n03 + n10)*log(1 - r) + (n02 + n11)*log(3 - r) + (n00 + n13)*log(r) + (n01 + n12)*log(2 + r)
    return(L)}
  
  inter_logLr <- function(n00,n01,n02,n03,n10,n11,n12,n13){
    optimize(logLr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_c <- parallel::mcmapply(inter_logLc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],
                            mc.cores = ncores) # keep one for now.. 
  r_r <- parallel::mcmapply(inter_logLr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],
                            mc.cores = ncores)
  
  LOD_c <- (x[,"n_00"] + x[,"n_13"])*log10(1 - r_c) + (x[,"n_01"] + x[,"n_12"])*log10(3 - r_c) + (x[,"n_03"] + x[,"n_10"])*log10(r_c) + (x[,"n_02"] + x[,"n_11"])*log10(2 + r_c) - (x[,"n_00"] + x[,"n_13"])*log10(1 - 0.5) - (x[,"n_01"] + x[,"n_12"])*log10(3 - 0.5) - (x[,"n_03"] + x[,"n_10"])*log10(0.5) - (x[,"n_02"] + x[,"n_11"])*log10(2 + 0.5)
  LOD_r <- (x[,"n_03"] + x[,"n_10"])*log10(1 - r_r) + (x[,"n_02"] + x[,"n_11"])*log10(3 - r_r) + (x[,"n_00"] + x[,"n_13"])*log10(r_r) + (x[,"n_01"] + x[,"n_12"])*log10(2 + r_r) - (x[,"n_03"] + x[,"n_10"])*log10(1 - 0.5) - (x[,"n_02"] + x[,"n_11"])*log10(3 - 0.5) - (x[,"n_00"] + x[,"n_13"])*log10(0.5) - (x[,"n_01"] + x[,"n_12"])*log10(2 + 0.5)
  
  return(list(r_mat=cbind(r_c, r_r),
              LOD_mat=cbind(LOD_c, LOD_r),
              logL_mat=NULL,
              phasing_strategy="MINR", 
              possible_phases=c("coupling",
                                "repulsion")))
}

#' @rdname r4_functions
#' @noRd
r4_1.0_1.2<-function(x, ncores=1){
  
  ############################
  ## COUPLING 
  ############################
  logLc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13){ 
    L <- (-n00 - n02 - n03 - n10 - n11 - n13)*log(12) + 
      (n00 + n13)*log(1 - r) + (n01 + n12)*log(1/3 - r/4) + 
      (n03 + n10)*log(r) + (n02 + n11)*log(1 + 3*r)
    return(L)}
  
  inter_logLc <- function(n00,n01,n02,n03,n10,n11,n12,n13){
    optimize(logLc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,
             maximum=T,lower = 0, upper = 0.5)$maximum}

  ############################
  ## REPULSION
  ############################
  logLr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13){ 
    L <- (-2*n00 - n01 - n02 - 2*n03 - 2*n10 - n11 - n12 - 2*n13)*log(3) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13)*log(4) + 
      (n03 + n10)*log(2 - r) + (n02 + n11)*log(3 - r) + (n00 + n13)*log(1 + r) + 
      (n01 + n12)*log(2 + r)
    return(L)}
  
  inter_logLr <- function(n00,n01,n02,n03,n10,n11,n12,n13){
    optimize(logLr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_c <- parallel::mcmapply(inter_logLc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],
                            mc.cores = ncores)
  
  r_r <- parallel::mcmapply(inter_logLr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],
                            mc.cores = ncores)
  
  LOD_c <- (x[,"n_00"] + x[,"n_13"])*log10(1 - r_c) + (x[,"n_01"] + x[,"n_12"])*log10(1/3 - r_c/4) + (x[,"n_03"] + x[,"n_10"])*log10(r_c) + (x[,"n_02"] + x[,"n_11"])*log10(1 + 3*r_c) - (x[,"n_00"] + x[,"n_13"])*log10(1 - 0.5) - (x[,"n_01"] + x[,"n_12"])*log10(1/3 - 0.5/4) - (x[,"n_03"] + x[,"n_10"])*log10(0.5) - (x[,"n_02"] + x[,"n_11"])*log10(1 + 3*0.5)
  LOD_r <- (x[,"n_03"] + x[,"n_10"])*log10(2 - r_r) + (x[,"n_02"] + x[,"n_11"])*log10(3 - r_r) + (x[,"n_00"] + x[,"n_13"])*log10(1 + r_r) + (x[,"n_01"] + x[,"n_12"])*log10(2 + r_r) - (x[,"n_03"] + x[,"n_10"])*log10(2 - 0.5) - (x[,"n_02"] + x[,"n_11"])*log10(3 - 0.5) - (x[,"n_00"] + x[,"n_13"])*log10(1 + 0.5) - (x[,"n_01"] + x[,"n_12"])*log10(2 + 0.5)

  logL_c <- (-x[,"n_00"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_13"])*log(12) + 
    (x[,"n_00"] + x[,"n_13"])*log(1 - r_c) + (x[,"n_01"] + x[,"n_12"])*log(1/3 - r_c/4) + 
    (x[,"n_03"] + x[,"n_10"])*log(r_c) + (x[,"n_02"] + x[,"n_11"])*log(1 + 3*r_c)
  
  logL_r <- (-2*x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"])*log(3) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"])*log(4) + 
    (x[,"n_03"] + x[,"n_10"])*log(2 - r_r) + (x[,"n_02"] + x[,"n_11"])*log(3 - r_r) + (x[,"n_00"] + x[,"n_13"])*log(1 + r_r) + 
    (x[,"n_01"] + x[,"n_12"])*log(2 + r_r)
  
  return(list(r_mat=cbind(r_c, r_r),
              LOD_mat=cbind(LOD_c, LOD_r),
              logL_mat=cbind(logL_c, logL_r),
              phasing_strategy="MLL", 
              possible_phases=c("coupling",
                                "repulsion")))
}

#' @rdname r4_functions
#' @noRd
r4_1.0_2.2<-function(x, ncores=1){
  
  logLc <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14){ 
    L <- (-2*n00 - 2*n02 - n03 - 2*n04 - 2*n10 - n11 - 2*n12 - 2*n14)*log(2) + 
      (-n00 - n03 - n04 - n10 - n11 - n14)*log(9) + (n00 + n14)*log(1 - r) + (n01 + n13)*log(1/6 - r/9) + 
      (n04 + n10)*log(r) + (n03 + n11)*log(1 + 2*r)
    return(L)}
  
  inter_logLc <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14){
    optimize(logLc,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  logLr <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14){ 
    L <- (-2*n00 - n01 - 2*n02 - 2*n04 - 2*n10 - 2*n12 - n13 - 2*n14)*log(2) + 
      (-n00 - n01 - n04 - n10 - n13 - n14)*log(9) + (n04 + n10)*log(1 - r) + (n03 + n11)*log(1/6 - r/9) + 
      (n00 + n14)*log(r) + (n01 + n13)*log(1 + 2*r)
    return(L)}
  
  inter_logLr <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14){
    optimize(logLr,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_c <- parallel::mcmapply(inter_logLc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                            mc.cores = ncores)

  r_r <- parallel::mcmapply(inter_logLr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                            mc.cores = ncores)
  
  LOD_c <- (x[,"n_00"] + x[,"n_14"])*log10(1 - r_c) + (x[,"n_01"] + x[,"n_13"])*log10(1/6 - r_c/9) + (x[,"n_04"] + x[,"n_10"])*log10(r_c) + (x[,"n_03"] + x[,"n_11"])*log10(1 + 2*r_c) - 
    (x[,"n_00"] + x[,"n_14"])*log10(1 - 0.5) - (x[,"n_01"] + x[,"n_13"])*log10(1/6 - 0.5/9) - (x[,"n_04"] + x[,"n_10"])*log10(0.5) - (x[,"n_03"] + x[,"n_11"])*log10(1 + 2*0.5)
  
  LOD_r <- (x[,"n_04"] + x[,"n_10"])*log10(1 - r_r) + (x[,"n_03"] + x[,"n_11"])*log10(1/6 - r_r/9) + (x[,"n_00"] + x[,"n_14"])*log10(r_r) + (x[,"n_01"] + x[,"n_13"])*log10(1 + 2*r_r) - 
    (x[,"n_04"] + x[,"n_10"])*log10(1 - 0.5) - (x[,"n_03"] + x[,"n_11"])*log10(1/6 - 0.5/9) - (x[,"n_00"] + x[,"n_14"])*log10(0.5) - (x[,"n_01"] + x[,"n_13"])*log10(1 + 2*0.5)
  
  return(list(r_mat=cbind(r_c, r_r),
              LOD_mat=cbind(LOD_c, LOD_r),
              logL_mat=NULL,
              phasing_strategy="MINR", 
              possible_phases=c("coupling",
                                "repulsion")))
}

#' @rdname r4_functions
#' @noRd
r4_2.0_1.1<-function(x, ncores=1){
  
  common_sum <- x[,"n_00"] + x[,"n_02"] + x[,"n_20"] + x[,"n_22"] 
  r_c <- (x[,"n_02"] + x[,"n_20"])/common_sum
  r_r <- (x[,"n_00"] + x[,"n_22"])/common_sum
  Lc <- function(r) {
    L <- (1 - r)^(x[,"n_00"] + x[,"n_22"])*r^(x[,"n_02"] + x[,"n_20"])
    return(L)}
  
  Lr <- function(r) {
    L <- (1 - r)^(x[,"n_02"] + x[,"n_20"])*r^(x[,"n_00"] + x[,"n_22"])
    return(L)}
  
  LOD_c <- log10(Lc(r_c)/Lc(0.5))
  LOD_r <- log10(Lr(r_r)/Lr(0.5))

  return(list(r_mat=cbind(r_c, r_r), 
              LOD_mat=cbind(LOD_c, LOD_r),
              logL_mat=NULL,
              phasing_strategy="MINR",
              possible_phases=c("coupling",
                                "repulsion")
  ))
}

#' @rdname r4_functions
#' @noRd
r4_2.0_2.1<-function(x, ncores=1){
  logLc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){ 
    L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - n10 - n11 - n12 - n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(3) + 
      (n00 + n23)*log((-1 + r)^2) + (n02 + n21)*log(-((-2 + r)*r)) + (n10 + n13)*log(-((-1 + r)*r)) + 
      (n03 + n20)*log(r^2) + (n01 + n22)*log(1 - r^2) + (n11 + n12)*log(2 - r + r^2)
    return(L)}
  
  inter_logLc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){
    optimize(logLc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  logLm <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){ 
    L <- (-3*n00 - 3*n01 - 3*n02 - 3*n03 - 2*n10 - 2*n11 - 2*n12 - 2*n13 - 3*n20 - 3*n21 - 3*n22 - 3*n23)*log(2) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(3) + 
      (n00 + n03 + n20 + n23)*log(-((-1 + r)*r)) + (n11 + n12)*log(3 + r - r^2) + (n10 + n13)*log(1 - r + r^2) + 
      (n01 + n02 + n21 + n22)*log(2 - r + r^2)
    return(L)}
  
  inter_logLm <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){
    optimize(logLm,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  logLr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){ 
    L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - n10 - n11 - n12 - n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(3) + 
      (n03 + n20)*log((-1 + r)^2) + (n01 + n22)*log(-((-2 + r)*r)) + (n10 + n13)*log(-((-1 + r)*r)) + 
      (n00 + n23)*log(r^2) + (n02 + n21)*log(1 - r^2) + (n11 + n12)*log(2 - r + r^2)
    return(L)}
  
  inter_logLr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){
    optimize(logLr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_c <- parallel::mcmapply(inter_logLc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                            mc.cores = ncores)

  r_m <- parallel::mcmapply(inter_logLm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                            mc.cores = ncores)

  r_r <- parallel::mcmapply(inter_logLr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                            mc.cores = ncores)
  
  LOD_c <- (x[,"n_00"] + x[,"n_23"])*log10((-1 + r_c)^2) + (x[,"n_02"] + x[,"n_21"])*log10(-((-2 + r_c)*r_c)) + (x[,"n_10"] + x[,"n_13"])*log10(-((-1 + r_c)*r_c)) + 
    (x[,"n_03"] + x[,"n_20"])*log10(r_c^2) + (x[,"n_01"] + x[,"n_22"])*log10(1 - r_c^2) + (x[,"n_11"] + x[,"n_12"])*log10(2 - r_c + r_c^2) - 
    (x[,"n_00"] + x[,"n_23"])*log10((-1 + 0.5)^2) - (x[,"n_02"] + x[,"n_21"])*log10(-((-2 + 0.5)*0.5)) - (x[,"n_10"] + x[,"n_13"])*log10(-((-1 + 0.5)*0.5)) - 
    (x[,"n_03"] + x[,"n_20"])*log10(0.5^2) - (x[,"n_01"] + x[,"n_22"])*log10(1 - 0.5^2) - (x[,"n_11"] + x[,"n_12"])*log10(2 - 0.5 + 0.5^2)
  
  LOD_m <- (x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*log10(-((-1 + r_m)*r_m)) + (x[,"n_11"] + x[,"n_12"])*log10(3 + r_m - r_m^2) + 
    (x[,"n_10"] + x[,"n_13"])*log10(1 - r_m + r_m^2) + (x[,"n_01"] + x[,"n_02"] + x[,"n_21"] + x[,"n_22"])*log10(2 - r_m + r_m^2) -
    (x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*log10(-((-1 + 0.5)*0.5)) - (x[,"n_11"] + x[,"n_12"])*log10(3 + 0.5 - 0.5^2) - 
    (x[,"n_10"] + x[,"n_13"])*log10(1 - 0.5 + 0.5^2) - (x[,"n_01"] + x[,"n_02"] + x[,"n_21"] + x[,"n_22"])*log10(2 - 0.5 + 0.5^2)
  
  LOD_r <- (x[,"n_03"] + x[,"n_20"])*log10((-1 + r_r)^2) + (x[,"n_01"] + x[,"n_22"])*log10(-((-2 + r_r)*r_r)) + (x[,"n_10"] + x[,"n_13"])*log10(-((-1 + r_r)*r_r)) + 
    (x[,"n_00"] + x[,"n_23"])*log10(r_r^2) + (x[,"n_02"] + x[,"n_21"])*log10(1 - r_r^2) + (x[,"n_11"] + x[,"n_12"])*log10(2 - r_r + r_r^2) - 
    (x[,"n_03"] + x[,"n_20"])*log10((-1 + 0.5)^2) - (x[,"n_01"] + x[,"n_22"])*log10(-((-2 + 0.5)*0.5)) - (x[,"n_10"] + x[,"n_13"])*log10(-((-1 + 0.5)*0.5)) - 
    (x[,"n_00"] + x[,"n_23"])*log10(0.5^2) - (x[,"n_02"] + x[,"n_21"])*log10(1 - 0.5^2) - (x[,"n_11"] + x[,"n_12"])*log10(2 - 0.5 + 0.5^2)
  
  ## We determine the logL for phasing:
  logL_c <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - 2*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(3) + 
    (x[,"n_00"] + x[,"n_23"])*log((-1 + r_c)^2) + (x[,"n_02"] + x[,"n_21"])*log(-((-2 + r_c)*r_c)) + (x[,"n_10"] + x[,"n_13"])*log(-((-1 + r_c)*r_c)) + 
    (x[,"n_03"] + x[,"n_20"])*log(r_c^2) + (x[,"n_01"] + x[,"n_22"])*log(1 - r_c^2) + (x[,"n_11"] + x[,"n_12"])*log(2 - r_c + r_c^2)
  
  logL_m <- (-3*x[,"n_00"] - 3*x[,"n_01"] - 3*x[,"n_02"] - 3*x[,"n_03"] - 2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - 3*x[,"n_20"] - 3*x[,"n_21"] - 3*x[,"n_22"] - 3*x[,"n_23"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(3) + 
    (x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*log(-((-1 + r_m)*r_m)) + (x[,"n_11"] + x[,"n_12"])*log(3 + r_m - r_m^2) + (x[,"n_10"] + x[,"n_13"])*log(1 - r_m + r_m^2) + 
    (x[,"n_01"] + x[,"n_02"] + x[,"n_21"] + x[,"n_22"])*log(2 - r_m + r_m^2)
  
  logL_r <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - 2*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(3) + 
    (x[,"n_03"] + x[,"n_20"])*log((-1 + r_r)^2) + (x[,"n_01"] + x[,"n_22"])*log(-((-2 + r_r)*r_r)) + (x[,"n_10"] + x[,"n_13"])*log(-((-1 + r_r)*r_r)) + 
    (x[,"n_00"] + x[,"n_23"])*log(r_r^2) + (x[,"n_02"] + x[,"n_21"])*log(1 - r_r^2) + (x[,"n_11"] + x[,"n_12"])*log(2 - r_r + r_r^2)
  
  return(list(r_mat=cbind(r_c, r_m, r_r),
              LOD_mat=cbind(LOD_c, LOD_m, LOD_r),
              logL_mat=cbind(logL_c, logL_m, logL_r),
              phasing_strategy="MLL", 
              possible_phases=c("coupling",
                                "mixed",
                                "repulsion")))

  

}

#' @rdname r4_functions
#' @noRd
r4_2.2_2.2 <- function(x,ncores=1){
  ## Function re-written Aug 2016 because of bug in older version. Should be correct now.
  
  logL_cc <- function(r,n00,n02,n04,n11,n13,n20,n24,n31,n33,n40,n42,n44,n01,n03,n10,n12,n14,n21,n23,n30,n32,n34,n41,n43,n22) {
    L <- (-2*n00 - n02 - 2*n04 + n11 + n13 - n20 - n24 + n31 + n33 - 2*n40 - n42 - 2*n44)*log(2) + (-2*n00 - 2*n01 - n02 - 2*n03 - 2*n04 - 2*n10 - 2*n11 - 2*n12 - 2*n13 - 2*n14 - n20 - 2*n21 - 2*n23 - n24 - 2*n30 - 2*n31 - 2*n32 - 2*n33 - 2*n34 - 2*n40 - 2*n41 - n42 - 2*n43 - 2*n44)*log(3) + 4*(n00 + n44)*log(pmax(1e-6,1 - r)) + 4*(n04 + n40)*log(pmax(1e-6,r)) + (n01 + n10 + n34 + n43)*(3*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n02 + n20 + n24 + n42)*(2*log(pmax(1e-6,1 - r)) + 2*log(pmax(1e-6,r))) + (n03 + n14 + n30 + n41)*(log(pmax(1e-6,1 - r)) + 3*log(pmax(1e-6,r))) + (n13 + n31)*(2*log(pmax(1e-6,r)) + log(2 - 3*r + 2*r^2)) + (n11 + n33)*(2*log(pmax(1e-6,1 - r)) + log(1 - r + 2*r^2)) + (n12 + n21 + n23 + n32)*(log(pmax(1e-6,r)) + log(5 - 11*r + 12*r^2 - 6*r^3)) + n22*log(1/2 - (10*r)/9 + (19*r^2)/9 - 2*r^3 + r^4)
    return(L)}
  interlogL_cc <- function(n00,n02,n04,n11,n13,n20,n24,n31,n33,n40,n42,n44,n01,n03,n10,n12,n14,n21,n23,n30,n32,n34,n41,n43,n22) {
    optimize(logL_cc,c(0,0.5), n00,n02,n04,n11,n13,n20,n24,n31,n33,n40,n42,n44,n01,n03,n10,n12,n14,n21,n23,n30,n32,n34,n41,n43,n22, maximum=TRUE, lower=0, upper=0.5)$maximum}
  
  
  r_cc <- parallel::mcmapply(interlogL_cc,x[,"n_00"],x[,"n_02"],x[,"n_04"],x[,"n_11"],x[,"n_13"],x[,"n_20"],x[,"n_24"],x[,"n_31"],x[,"n_33"],x[,"n_40"],x[,"n_42"],x[,"n_44"],x[,"n_01"],x[,"n_03"],x[,"n_10"],x[,"n_12"],x[,"n_14"],x[,"n_21"],x[,"n_23"],x[,"n_30"],x[,"n_32"],x[,"n_34"],x[,"n_41"],x[,"n_43"],x[,"n_22"], mc.cores = ncores)
  
  
  LOD_cc <- 2*(x[,"n_13"] + x[,"n_31"])*log10(2) + 2*(x[,"n_11"] + x[,"n_33"])*log10(2) + 4*(x[,"n_04"] + x[,"n_40"])*log10(2) + 4*(x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*log10(2) + 4*(x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*log10(2) + 4*(x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*log10(2) + 4*(x[,"n_00"] + x[,"n_44"])*log10(2) - (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(-3*log10(2) + log10(7)) + x[,"n_22"]*(4*log10(2) + 2*log10(3) - log10(41)) + 4*(x[,"n_00"] + x[,"n_44"])*log10(pmax(1e-6,1 - r_cc)) + 4*(x[,"n_04"] + x[,"n_40"])*log10(pmax(1e-6,r_cc)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(3*log10(pmax(1e-6,1 - r_cc)) + log10(pmax(1e-6,r_cc))) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(2*log10(pmax(1e-6,1 - r_cc)) + 2*log10(pmax(1e-6,r_cc))) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(log10(pmax(1e-6,1 - r_cc)) + 3*log10(pmax(1e-6,r_cc))) + (x[,"n_13"] + x[,"n_31"])*(2*log10(pmax(1e-6,r_cc)) + log10(2 - 3*r_cc + 2*r_cc^2)) + (x[,"n_11"] + x[,"n_33"])*(2*log10(pmax(1e-6,1 - r_cc)) + log10(1 - r_cc + 2*r_cc^2)) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(log10(pmax(1e-6,r_cc)) + log10(5 - 11*r_cc + 12*r_cc^2 - 6*r_cc^3)) + x[,"n_22"]*log10(1/2 - (10*r_cc)/9 + (19*r_cc^2)/9 - 2*r_cc^3 + r_cc^4)
  
  
  logL_cc <- (-2*x[,"n_00"] - x[,"n_02"] - 2*x[,"n_04"] + x[,"n_11"] + x[,"n_13"] - x[,"n_20"] - x[,"n_24"] + x[,"n_31"] + x[,"n_33"] - 2*x[,"n_40"] - x[,"n_42"] - 2*x[,"n_44"])*log(2) + (-2*x[,"n_00"] - 2*x[,"n_01"] - x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_14"] - x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_23"] - x[,"n_24"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"] - 2*x[,"n_34"] - 2*x[,"n_40"] - 2*x[,"n_41"] - x[,"n_42"] - 2*x[,"n_43"] - 2*x[,"n_44"])*log(3) + 4*(x[,"n_00"] + x[,"n_44"])*log(pmax(1e-6,1 - r_cc)) + 4*(x[,"n_04"] + x[,"n_40"])*log(pmax(1e-6,r_cc)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(3*log(pmax(1e-6,1 - r_cc)) + log(pmax(1e-6,r_cc))) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(2*log(pmax(1e-6,1 - r_cc)) + 2*log(pmax(1e-6,r_cc))) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(log(pmax(1e-6,1 - r_cc)) + 3*log(pmax(1e-6,r_cc))) + (x[,"n_13"] + x[,"n_31"])*(2*log(pmax(1e-6,r_cc)) + log(2 - 3*r_cc + 2*r_cc^2)) + (x[,"n_11"] + x[,"n_33"])*(2*log(pmax(1e-6,1 - r_cc)) + log(1 - r_cc + 2*r_cc^2)) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(log(pmax(1e-6,r_cc)) + log(5 - 11*r_cc + 12*r_cc^2 - 6*r_cc^3)) + x[,"n_22"]*log(1/2 - (10*r_cc)/9 + (19*r_cc^2)/9 - 2*r_cc^3 + r_cc^4)
  
  
  logL_cm <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22) {
    L <- (-3*n00 - 2*n01 - 3*n02 - 2*n03 - 3*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 3*n20 - 2*n21 - 2*n23 - 3*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34 - 3*n40 - 2*n41 - 3*n42 - 2*n43 - 3*n44)*log(2) + 2*(-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n23 - n24 - n30 - n31 - n32 - n33 - n34 - n40 - n41 - n42 - n43 - n44)*log(3) + (n00 + n44)*(3*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n04 + n40)*(log(pmax(1e-6,1 - r)) + 3*log(pmax(1e-6,r))) + (n03 + n14 + n30 + n41)*(2*log(pmax(1e-6,r)) + log(2 - 3*r + 2*r^2)) + (n01 + n10 + n34 + n43)*(2*log(pmax(1e-6,1 - r)) + log(1 - r + 2*r^2)) + (n02 + n20 + n24 + n42)*(log(pmax(1e-6,r)) + log(5 - 11*r + 12*r^2 - 6*r^3)) + (n13 + n31)*(log(pmax(1e-6,r)) + log(3 - 5*r + 7*r^2 - 4*r^3)) + (n11 + n33)*log(1 + 2*r - 8*r^2 + 9*r^3 - 4*r^4) + n22*log(2/9 + r/4 - (3*r^2)/4 + r^3 - r^4/2) + (n12 + n21 + n23 + n32)*log(5 - 7*r + 19*r^2 - 24*r^3 + 12*r^4)
    return(L)}
  interlogL_cm <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22) {
    optimize(logL_cm,c(0,0.5), n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22, maximum=TRUE, lower=0, upper=0.5)$maximum}
  
  
  r_cm <- parallel::mcmapply(interlogL_cm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],x[,"n_40"],x[,"n_41"],x[,"n_42"],x[,"n_43"],x[,"n_44"],x[,"n_22"], mc.cores = ncores)
  
  
  LOD_cm <- -2*(x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(2) + 4*(x[,"n_04"] + x[,"n_40"])*log10(2) + 2*(x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*log10(2) + 2*(x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*log10(2) + 4*(x[,"n_00"] + x[,"n_44"])*log10(2) + (x[,"n_11"] + x[,"n_33"])*(3*log10(2) - log10(7)) - (x[,"n_13"] + x[,"n_31"])*(-3*log10(2) + log10(7)) - (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(-3*log10(2) + log10(7)) + x[,"n_22"]*(5*log10(2) + 2*log10(3) - log10(73)) + (x[,"n_00"] + x[,"n_44"])*(3*log10(pmax(1e-6,1 - r_cm)) + log10(pmax(1e-6,r_cm))) + (x[,"n_04"] + x[,"n_40"])*(log10(pmax(1e-6,1 - r_cm)) + 3*log10(pmax(1e-6,r_cm))) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(2*log10(pmax(1e-6,r_cm)) + log10(2 - 3*r_cm + 2*r_cm^2)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(2*log10(pmax(1e-6,1 - r_cm)) + log10(1 - r_cm + 2*r_cm^2)) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(log10(pmax(1e-6,r_cm)) + log10(5 - 11*r_cm + 12*r_cm^2 - 6*r_cm^3)) + (x[,"n_13"] + x[,"n_31"])*(log10(pmax(1e-6,r_cm)) + log10(3 - 5*r_cm + 7*r_cm^2 - 4*r_cm^3)) + (x[,"n_11"] + x[,"n_33"])*log10(1 + 2*r_cm - 8*r_cm^2 + 9*r_cm^3 - 4*r_cm^4) + x[,"n_22"]*log10(2/9 + r_cm/4 - (3*r_cm^2)/4 + r_cm^3 - r_cm^4/2) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(5 - 7*r_cm + 19*r_cm^2 - 24*r_cm^3 + 12*r_cm^4)
  
  
  logL_cm <- (-3*x[,"n_00"] - 2*x[,"n_01"] - 3*x[,"n_02"] - 2*x[,"n_03"] - 3*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 3*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_23"] - 3*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"] - 3*x[,"n_40"] - 2*x[,"n_41"] - 3*x[,"n_42"] - 2*x[,"n_43"] - 3*x[,"n_44"])*log(2) + 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"] - x[,"n_40"] - x[,"n_41"] - x[,"n_42"] - x[,"n_43"] - x[,"n_44"])*log(3) + (x[,"n_00"] + x[,"n_44"])*(3*log(pmax(1e-6,1 - r_cm)) + log(pmax(1e-6,r_cm))) + (x[,"n_04"] + x[,"n_40"])*(log(pmax(1e-6,1 - r_cm)) + 3*log(pmax(1e-6,r_cm))) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(2*log(pmax(1e-6,r_cm)) + log(2 - 3*r_cm + 2*r_cm^2)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(2*log(pmax(1e-6,1 - r_cm)) + log(1 - r_cm + 2*r_cm^2)) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(log(pmax(1e-6,r_cm)) + log(5 - 11*r_cm + 12*r_cm^2 - 6*r_cm^3)) + (x[,"n_13"] + x[,"n_31"])*(log(pmax(1e-6,r_cm)) + log(3 - 5*r_cm + 7*r_cm^2 - 4*r_cm^3)) + (x[,"n_11"] + x[,"n_33"])*log(1 + 2*r_cm - 8*r_cm^2 + 9*r_cm^3 - 4*r_cm^4) + x[,"n_22"]*log(2/9 + r_cm/4 - (3*r_cm^2)/4 + r_cm^3 - r_cm^4/2) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log(5 - 7*r_cm + 19*r_cm^2 - 24*r_cm^3 + 12*r_cm^4)
  
  
  logL_cr <- function(r,n00,n01,n02,n03,n04,n10,n14,n20,n24,n30,n34,n40,n41,n42,n43,n44,n11,n12,n13,n21,n23,n31,n32,n33,n22) {
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - n10 - n14 - 2*n20 - 2*n24 - n30 - n34 - 2*n40 - n41 - 2*n42 - n43 - 2*n44)*log(2) + 2*(-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n23 - n24 - n30 - n31 - n32 - n33 - n34 - n40 - n41 - n42 - n43 - n44)*log(3) + (n00 + n04 + n40 + n44)*(2*log(pmax(1e-6,1 - r)) + 2*log(pmax(1e-6,r))) + (n12 + n21 + n23 + n32)*(log(pmax(1e-6,r)) + log(5 - 11*r + 12*r^2 - 6*r^3)) + (n01 + n03 + n10 + n14 + n30 + n34 + n41 + n43)*log(r - 3*r^2 + 4*r^3 - 2*r^4) + n22*log(4/9 - (8*r)/9 + (17*r^2)/9 - 2*r^3 + r^4) + (n11 + n13 + n31 + n33)*log(1 - 3*r + 7*r^2 - 8*r^3 + 4*r^4) + (n02 + n20 + n24 + n42)*log(1 - 4*r + 10*r^2 - 12*r^3 + 6*r^4)
    return(L)}
  interlogL_cr <- function(n00,n01,n02,n03,n04,n10,n14,n20,n24,n30,n34,n40,n41,n42,n43,n44,n11,n12,n13,n21,n23,n31,n32,n33,n22) {
    optimize(logL_cr,c(0,0.5), n00,n01,n02,n03,n04,n10,n14,n20,n24,n30,n34,n40,n41,n42,n43,n44,n11,n12,n13,n21,n23,n31,n32,n33,n22, maximum=TRUE, lower=0, upper=0.5)$maximum}
  
  
  r_cr <- parallel::mcmapply(interlogL_cr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_14"],x[,"n_20"],x[,"n_24"],x[,"n_30"],x[,"n_34"],x[,"n_40"],x[,"n_41"],x[,"n_42"],x[,"n_43"],x[,"n_44"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_21"],x[,"n_23"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_22"], mc.cores = ncores)
  
  
  LOD_cr <- (x[,"n_11"] + x[,"n_13"] + x[,"n_31"] + x[,"n_33"])*log10(2) + 3*(x[,"n_01"] + x[,"n_03"] + x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"] + x[,"n_41"] + x[,"n_43"])*log10(2) + 4*(x[,"n_00"] + x[,"n_04"] + x[,"n_40"] + x[,"n_44"])*log10(2) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(3*log10(2) - log10(3)) - (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(-3*log10(2) + log10(7)) + x[,"n_22"]*(4*log10(2) + 2*log10(3) - log10(41)) + (x[,"n_00"] + x[,"n_04"] + x[,"n_40"] + x[,"n_44"])*(2*log10(pmax(1e-6,1 - r_cr)) + 2*log10(pmax(1e-6,r_cr))) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(log10(pmax(1e-6,r_cr)) + log10(5 - 11*r_cr + 12*r_cr^2 - 6*r_cr^3)) + (x[,"n_01"] + x[,"n_03"] + x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"] + x[,"n_41"] + x[,"n_43"])*log10(r_cr - 3*r_cr^2 + 4*r_cr^3 - 2*r_cr^4) + x[,"n_22"]*log10(4/9 - (8*r_cr)/9 + (17*r_cr^2)/9 - 2*r_cr^3 + r_cr^4) + (x[,"n_11"] + x[,"n_13"] + x[,"n_31"] + x[,"n_33"])*log10(1 - 3*r_cr + 7*r_cr^2 - 8*r_cr^3 + 4*r_cr^4) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*log10(1 - 4*r_cr + 10*r_cr^2 - 12*r_cr^3 + 6*r_cr^4)
  
  
  logL_cr <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - x[,"n_10"] - x[,"n_14"] - 2*x[,"n_20"] - 2*x[,"n_24"] - x[,"n_30"] - x[,"n_34"] - 2*x[,"n_40"] - x[,"n_41"] - 2*x[,"n_42"] - x[,"n_43"] - 2*x[,"n_44"])*log(2) + 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"] - x[,"n_40"] - x[,"n_41"] - x[,"n_42"] - x[,"n_43"] - x[,"n_44"])*log(3) + (x[,"n_00"] + x[,"n_04"] + x[,"n_40"] + x[,"n_44"])*(2*log(pmax(1e-6,1 - r_cr)) + 2*log(pmax(1e-6,r_cr))) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(log(pmax(1e-6,r_cr)) + log(5 - 11*r_cr + 12*r_cr^2 - 6*r_cr^3)) + (x[,"n_01"] + x[,"n_03"] + x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"] + x[,"n_41"] + x[,"n_43"])*log(r_cr - 3*r_cr^2 + 4*r_cr^3 - 2*r_cr^4) + x[,"n_22"]*log(4/9 - (8*r_cr)/9 + (17*r_cr^2)/9 - 2*r_cr^3 + r_cr^4) + (x[,"n_11"] + x[,"n_13"] + x[,"n_31"] + x[,"n_33"])*log(1 - 3*r_cr + 7*r_cr^2 - 8*r_cr^3 + 4*r_cr^4) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*log(1 - 4*r_cr + 10*r_cr^2 - 12*r_cr^3 + 6*r_cr^4)
  
  
  logL_mc <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22) {
    L <- (-3*n00 - 2*n01 - 3*n02 - 2*n03 - 3*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 3*n20 - 2*n21 - 2*n23 - 3*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34 - 3*n40 - 2*n41 - 3*n42 - 2*n43 - 3*n44)*log(2) + 2*(-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n23 - n24 - n30 - n31 - n32 - n33 - n34 - n40 - n41 - n42 - n43 - n44)*log(3) + (n00 + n44)*(3*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n04 + n40)*(log(pmax(1e-6,1 - r)) + 3*log(pmax(1e-6,r))) + (n03 + n14 + n30 + n41)*(2*log(pmax(1e-6,r)) + log(2 - 3*r + 2*r^2)) + (n01 + n10 + n34 + n43)*(2*log(pmax(1e-6,1 - r)) + log(1 - r + 2*r^2)) + (n02 + n20 + n24 + n42)*(log(pmax(1e-6,r)) + log(5 - 11*r + 12*r^2 - 6*r^3)) + (n13 + n31)*(log(pmax(1e-6,r)) + log(3 - 5*r + 7*r^2 - 4*r^3)) + (n11 + n33)*log(1 + 2*r - 8*r^2 + 9*r^3 - 4*r^4) + n22*log(2/9 + r/4 - (3*r^2)/4 + r^3 - r^4/2) + (n12 + n21 + n23 + n32)*log(5 - 7*r + 19*r^2 - 24*r^3 + 12*r^4)
    return(L)}
  interlogL_mc <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22) {
    optimize(logL_mc,c(0,0.5), n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22, maximum=TRUE, lower=0, upper=0.5)$maximum}
  
  
  r_mc <- parallel::mcmapply(interlogL_mc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],x[,"n_40"],x[,"n_41"],x[,"n_42"],x[,"n_43"],x[,"n_44"],x[,"n_22"], mc.cores = ncores)
  
  
  LOD_mc <- -2*(x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(2) + 4*(x[,"n_04"] + x[,"n_40"])*log10(2) + 2*(x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*log10(2) + 2*(x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*log10(2) + 4*(x[,"n_00"] + x[,"n_44"])*log10(2) + (x[,"n_11"] + x[,"n_33"])*(3*log10(2) - log10(7)) - (x[,"n_13"] + x[,"n_31"])*(-3*log10(2) + log10(7)) - (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(-3*log10(2) + log10(7)) + x[,"n_22"]*(5*log10(2) + 2*log10(3) - log10(73)) + (x[,"n_00"] + x[,"n_44"])*(3*log10(pmax(1e-6,1 - r_mc)) + log10(pmax(1e-6,r_mc))) + (x[,"n_04"] + x[,"n_40"])*(log10(pmax(1e-6,1 - r_mc)) + 3*log10(pmax(1e-6,r_mc))) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(2*log10(pmax(1e-6,r_mc)) + log10(2 - 3*r_mc + 2*r_mc^2)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(2*log10(pmax(1e-6,1 - r_mc)) + log10(1 - r_mc + 2*r_mc^2)) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(log10(pmax(1e-6,r_mc)) + log10(5 - 11*r_mc + 12*r_mc^2 - 6*r_mc^3)) + (x[,"n_13"] + x[,"n_31"])*(log10(pmax(1e-6,r_mc)) + log10(3 - 5*r_mc + 7*r_mc^2 - 4*r_mc^3)) + (x[,"n_11"] + x[,"n_33"])*log10(1 + 2*r_mc - 8*r_mc^2 + 9*r_mc^3 - 4*r_mc^4) + x[,"n_22"]*log10(2/9 + r_mc/4 - (3*r_mc^2)/4 + r_mc^3 - r_mc^4/2) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(5 - 7*r_mc + 19*r_mc^2 - 24*r_mc^3 + 12*r_mc^4)
  
  
  logL_mc <- (-3*x[,"n_00"] - 2*x[,"n_01"] - 3*x[,"n_02"] - 2*x[,"n_03"] - 3*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 3*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_23"] - 3*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"] - 3*x[,"n_40"] - 2*x[,"n_41"] - 3*x[,"n_42"] - 2*x[,"n_43"] - 3*x[,"n_44"])*log(2) + 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"] - x[,"n_40"] - x[,"n_41"] - x[,"n_42"] - x[,"n_43"] - x[,"n_44"])*log(3) + (x[,"n_00"] + x[,"n_44"])*(3*log(pmax(1e-6,1 - r_mc)) + log(pmax(1e-6,r_mc))) + (x[,"n_04"] + x[,"n_40"])*(log(pmax(1e-6,1 - r_mc)) + 3*log(pmax(1e-6,r_mc))) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(2*log(pmax(1e-6,r_mc)) + log(2 - 3*r_mc + 2*r_mc^2)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(2*log(pmax(1e-6,1 - r_mc)) + log(1 - r_mc + 2*r_mc^2)) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(log(pmax(1e-6,r_mc)) + log(5 - 11*r_mc + 12*r_mc^2 - 6*r_mc^3)) + (x[,"n_13"] + x[,"n_31"])*(log(pmax(1e-6,r_mc)) + log(3 - 5*r_mc + 7*r_mc^2 - 4*r_mc^3)) + (x[,"n_11"] + x[,"n_33"])*log(1 + 2*r_mc - 8*r_mc^2 + 9*r_mc^3 - 4*r_mc^4) + x[,"n_22"]*log(2/9 + r_mc/4 - (3*r_mc^2)/4 + r_mc^3 - r_mc^4/2) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log(5 - 7*r_mc + 19*r_mc^2 - 24*r_mc^3 + 12*r_mc^4)
  
  
  logL_mm <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44) {
    L <- (-4*n00 - 2*n01 - 3*n02 - 2*n03 - 4*n04 - 2*n10 - n11 - n12 - n13 - 2*n14 - 3*n20 - n21 - 2*n22 - n23 - 3*n24 - 2*n30 - n31 - n32 - n33 - 2*n34 - 4*n40 - 2*n41 - 3*n42 - 2*n43 - 4*n44)*log(2) + 2*(-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34 - n40 - n41 - n42 - n43 - n44)*log(3) + (n00 + n04 + n40 + n44)*(2*log(pmax(1e-6,1 - r)) + 2*log(pmax(1e-6,r))) + (n01 + n03 + n10 + n14 + n30 + n34 + n41 + n43)*(log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r)) + log(1 - r + r^2)) + (n12 + n21 + n23 + n32)*log(2 + r - 4*r^2 + 6*r^3 - 3*r^4) + (n11 + n13 + n31 + n33)*log(1 - r + 3*r^2 - 4*r^3 + 2*r^4) + (n02 + n20 + n24 + n42)*log(2 - 4*r + 7*r^2 - 6*r^3 + 3*r^4) + n22*log(8 + 9*r^2 - 18*r^3 + 9*r^4)
    return(L)}
  interlogL_mm <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44) {
    optimize(logL_mm,c(0,0.5), n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44, maximum=TRUE, lower=0, upper=0.5)$maximum}
  
  
  r_mm <- parallel::mcmapply(interlogL_mm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],x[,"n_40"],x[,"n_41"],x[,"n_42"],x[,"n_43"],x[,"n_44"], mc.cores = ncores)
  
  
  LOD_mm <- 4*(x[,"n_00"] + x[,"n_04"] + x[,"n_40"] + x[,"n_44"])*log10(2) - (x[,"n_01"] + x[,"n_03"] + x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"] + x[,"n_41"] + x[,"n_43"])*(-4*log10(2) + log10(3)) + (x[,"n_11"] + x[,"n_13"] + x[,"n_31"] + x[,"n_33"])*(3*log10(2) - log10(7)) - (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(-4*log10(2) + log10(3) + log10(11)) - (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(-4*log10(2) + log10(19)) - x[,"n_22"]*(-4*log10(2) + log10(137)) + (x[,"n_00"] + x[,"n_04"] + x[,"n_40"] + x[,"n_44"])*(2*log10(pmax(1e-6,1 - r_mm)) + 2*log10(pmax(1e-6,r_mm))) + (x[,"n_01"] + x[,"n_03"] + x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"] + x[,"n_41"] + x[,"n_43"])*(log10(pmax(1e-6,1 - r_mm)) + log10(pmax(1e-6,r_mm)) + log10(1 - r_mm + r_mm^2)) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(2 + r_mm - 4*r_mm^2 + 6*r_mm^3 - 3*r_mm^4) + (x[,"n_11"] + x[,"n_13"] + x[,"n_31"] + x[,"n_33"])*log10(1 - r_mm + 3*r_mm^2 - 4*r_mm^3 + 2*r_mm^4) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*log10(2 - 4*r_mm + 7*r_mm^2 - 6*r_mm^3 + 3*r_mm^4) + x[,"n_22"]*log10(8 + 9*r_mm^2 - 18*r_mm^3 + 9*r_mm^4)
  
  
  logL_mm <- (-4*x[,"n_00"] - 2*x[,"n_01"] - 3*x[,"n_02"] - 2*x[,"n_03"] - 4*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 3*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 3*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"] - 4*x[,"n_40"] - 2*x[,"n_41"] - 3*x[,"n_42"] - 2*x[,"n_43"] - 4*x[,"n_44"])*log(2) + 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"] - x[,"n_40"] - x[,"n_41"] - x[,"n_42"] - x[,"n_43"] - x[,"n_44"])*log(3) + (x[,"n_00"] + x[,"n_04"] + x[,"n_40"] + x[,"n_44"])*(2*log(pmax(1e-6,1 - r_mm)) + 2*log(pmax(1e-6,r_mm))) + (x[,"n_01"] + x[,"n_03"] + x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"] + x[,"n_41"] + x[,"n_43"])*(log(pmax(1e-6,1 - r_mm)) + log(pmax(1e-6,r_mm)) + log(1 - r_mm + r_mm^2)) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log(2 + r_mm - 4*r_mm^2 + 6*r_mm^3 - 3*r_mm^4) + (x[,"n_11"] + x[,"n_13"] + x[,"n_31"] + x[,"n_33"])*log(1 - r_mm + 3*r_mm^2 - 4*r_mm^3 + 2*r_mm^4) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*log(2 - 4*r_mm + 7*r_mm^2 - 6*r_mm^3 + 3*r_mm^4) + x[,"n_22"]*log(8 + 9*r_mm^2 - 18*r_mm^3 + 9*r_mm^4)
  
  
  logL_mr <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22) {
    L <- (-3*n00 - 2*n01 - 3*n02 - 2*n03 - 3*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 3*n20 - 2*n21 - 2*n23 - 3*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34 - 3*n40 - 2*n41 - 3*n42 - 2*n43 - 3*n44)*log(2) + 2*(-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n23 - n24 - n30 - n31 - n32 - n33 - n34 - n40 - n41 - n42 - n43 - n44)*log(3) + (n04 + n40)*(3*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n00 + n44)*(log(pmax(1e-6,1 - r)) + 3*log(pmax(1e-6,r))) + (n01 + n10 + n34 + n43)*(2*log(pmax(1e-6,r)) + log(2 - 3*r + 2*r^2)) + (n03 + n14 + n30 + n41)*(2*log(pmax(1e-6,1 - r)) + log(1 - r + 2*r^2)) + (n02 + n20 + n24 + n42)*(log(pmax(1e-6,r)) + log(5 - 11*r + 12*r^2 - 6*r^3)) + (n11 + n33)*(log(pmax(1e-6,r)) + log(3 - 5*r + 7*r^2 - 4*r^3)) + (n13 + n31)*log(1 + 2*r - 8*r^2 + 9*r^3 - 4*r^4) + n22*log(2/9 + r/4 - (3*r^2)/4 + r^3 - r^4/2) + (n12 + n21 + n23 + n32)*log(5 - 7*r + 19*r^2 - 24*r^3 + 12*r^4)
    return(L)}
  interlogL_mr <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22) {
    optimize(logL_mr,c(0,0.5), n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22, maximum=TRUE, lower=0, upper=0.5)$maximum}
  
  
  r_mr <- parallel::mcmapply(interlogL_mr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],x[,"n_40"],x[,"n_41"],x[,"n_42"],x[,"n_43"],x[,"n_44"],x[,"n_22"], mc.cores = ncores)
  
  
  LOD_mr <- -2*(x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(2) + 4*(x[,"n_04"] + x[,"n_40"])*log10(2) + 2*(x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*log10(2) + 2*(x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*log10(2) + 4*(x[,"n_00"] + x[,"n_44"])*log10(2) + (x[,"n_13"] + x[,"n_31"])*(3*log10(2) - log10(7)) - (x[,"n_11"] + x[,"n_33"])*(-3*log10(2) + log10(7)) - (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(-3*log10(2) + log10(7)) + x[,"n_22"]*(5*log10(2) + 2*log10(3) - log10(73)) + (x[,"n_04"] + x[,"n_40"])*(3*log10(pmax(1e-6,1 - r_mr)) + log10(pmax(1e-6,r_mr))) + (x[,"n_00"] + x[,"n_44"])*(log10(pmax(1e-6,1 - r_mr)) + 3*log10(pmax(1e-6,r_mr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(2*log10(pmax(1e-6,r_mr)) + log10(2 - 3*r_mr + 2*r_mr^2)) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(2*log10(pmax(1e-6,1 - r_mr)) + log10(1 - r_mr + 2*r_mr^2)) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(log10(pmax(1e-6,r_mr)) + log10(5 - 11*r_mr + 12*r_mr^2 - 6*r_mr^3)) + (x[,"n_11"] + x[,"n_33"])*(log10(pmax(1e-6,r_mr)) + log10(3 - 5*r_mr + 7*r_mr^2 - 4*r_mr^3)) + (x[,"n_13"] + x[,"n_31"])*log10(1 + 2*r_mr - 8*r_mr^2 + 9*r_mr^3 - 4*r_mr^4) + x[,"n_22"]*log10(2/9 + r_mr/4 - (3*r_mr^2)/4 + r_mr^3 - r_mr^4/2) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(5 - 7*r_mr + 19*r_mr^2 - 24*r_mr^3 + 12*r_mr^4)
  
  
  logL_mr <- (-3*x[,"n_00"] - 2*x[,"n_01"] - 3*x[,"n_02"] - 2*x[,"n_03"] - 3*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 3*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_23"] - 3*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"] - 3*x[,"n_40"] - 2*x[,"n_41"] - 3*x[,"n_42"] - 2*x[,"n_43"] - 3*x[,"n_44"])*log(2) + 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"] - x[,"n_40"] - x[,"n_41"] - x[,"n_42"] - x[,"n_43"] - x[,"n_44"])*log(3) + (x[,"n_04"] + x[,"n_40"])*(3*log(pmax(1e-6,1 - r_mr)) + log(pmax(1e-6,r_mr))) + (x[,"n_00"] + x[,"n_44"])*(log(pmax(1e-6,1 - r_mr)) + 3*log(pmax(1e-6,r_mr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(2*log(pmax(1e-6,r_mr)) + log(2 - 3*r_mr + 2*r_mr^2)) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(2*log(pmax(1e-6,1 - r_mr)) + log(1 - r_mr + 2*r_mr^2)) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(log(pmax(1e-6,r_mr)) + log(5 - 11*r_mr + 12*r_mr^2 - 6*r_mr^3)) + (x[,"n_11"] + x[,"n_33"])*(log(pmax(1e-6,r_mr)) + log(3 - 5*r_mr + 7*r_mr^2 - 4*r_mr^3)) + (x[,"n_13"] + x[,"n_31"])*log(1 + 2*r_mr - 8*r_mr^2 + 9*r_mr^3 - 4*r_mr^4) + x[,"n_22"]*log(2/9 + r_mr/4 - (3*r_mr^2)/4 + r_mr^3 - r_mr^4/2) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log(5 - 7*r_mr + 19*r_mr^2 - 24*r_mr^3 + 12*r_mr^4)
  
  
  logL_rc <- function(r,n00,n01,n02,n03,n04,n10,n14,n20,n24,n30,n34,n40,n41,n42,n43,n44,n11,n12,n13,n21,n23,n31,n32,n33,n22) {
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - n10 - n14 - 2*n20 - 2*n24 - n30 - n34 - 2*n40 - n41 - 2*n42 - n43 - 2*n44)*log(2) + 2*(-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n23 - n24 - n30 - n31 - n32 - n33 - n34 - n40 - n41 - n42 - n43 - n44)*log(3) + (n00 + n04 + n40 + n44)*(2*log(pmax(1e-6,1 - r)) + 2*log(pmax(1e-6,r))) + (n12 + n21 + n23 + n32)*(log(pmax(1e-6,r)) + log(5 - 11*r + 12*r^2 - 6*r^3)) + (n01 + n03 + n10 + n14 + n30 + n34 + n41 + n43)*log(r - 3*r^2 + 4*r^3 - 2*r^4) + n22*log(4/9 - (8*r)/9 + (17*r^2)/9 - 2*r^3 + r^4) + (n11 + n13 + n31 + n33)*log(1 - 3*r + 7*r^2 - 8*r^3 + 4*r^4) + (n02 + n20 + n24 + n42)*log(1 - 4*r + 10*r^2 - 12*r^3 + 6*r^4)
    return(L)}
  interlogL_rc <- function(n00,n01,n02,n03,n04,n10,n14,n20,n24,n30,n34,n40,n41,n42,n43,n44,n11,n12,n13,n21,n23,n31,n32,n33,n22) {
    optimize(logL_rc,c(0,0.5), n00,n01,n02,n03,n04,n10,n14,n20,n24,n30,n34,n40,n41,n42,n43,n44,n11,n12,n13,n21,n23,n31,n32,n33,n22, maximum=TRUE, lower=0, upper=0.5)$maximum}
  
  
  r_rc <- parallel::mcmapply(interlogL_rc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_14"],x[,"n_20"],x[,"n_24"],x[,"n_30"],x[,"n_34"],x[,"n_40"],x[,"n_41"],x[,"n_42"],x[,"n_43"],x[,"n_44"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_21"],x[,"n_23"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_22"], mc.cores = ncores)
  
  
  LOD_rc <- (x[,"n_11"] + x[,"n_13"] + x[,"n_31"] + x[,"n_33"])*log10(2) + 3*(x[,"n_01"] + x[,"n_03"] + x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"] + x[,"n_41"] + x[,"n_43"])*log10(2) + 4*(x[,"n_00"] + x[,"n_04"] + x[,"n_40"] + x[,"n_44"])*log10(2) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(3*log10(2) - log10(3)) - (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(-3*log10(2) + log10(7)) + x[,"n_22"]*(4*log10(2) + 2*log10(3) - log10(41)) + (x[,"n_00"] + x[,"n_04"] + x[,"n_40"] + x[,"n_44"])*(2*log10(pmax(1e-6,1 - r_rc)) + 2*log10(pmax(1e-6,r_rc))) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(log10(pmax(1e-6,r_rc)) + log10(5 - 11*r_rc + 12*r_rc^2 - 6*r_rc^3)) + (x[,"n_01"] + x[,"n_03"] + x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"] + x[,"n_41"] + x[,"n_43"])*log10(r_rc - 3*r_rc^2 + 4*r_rc^3 - 2*r_rc^4) + x[,"n_22"]*log10(4/9 - (8*r_rc)/9 + (17*r_rc^2)/9 - 2*r_rc^3 + r_rc^4) + (x[,"n_11"] + x[,"n_13"] + x[,"n_31"] + x[,"n_33"])*log10(1 - 3*r_rc + 7*r_rc^2 - 8*r_rc^3 + 4*r_rc^4) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*log10(1 - 4*r_rc + 10*r_rc^2 - 12*r_rc^3 + 6*r_rc^4)
  
  
  logL_rc <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - x[,"n_10"] - x[,"n_14"] - 2*x[,"n_20"] - 2*x[,"n_24"] - x[,"n_30"] - x[,"n_34"] - 2*x[,"n_40"] - x[,"n_41"] - 2*x[,"n_42"] - x[,"n_43"] - 2*x[,"n_44"])*log(2) + 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"] - x[,"n_40"] - x[,"n_41"] - x[,"n_42"] - x[,"n_43"] - x[,"n_44"])*log(3) + (x[,"n_00"] + x[,"n_04"] + x[,"n_40"] + x[,"n_44"])*(2*log(pmax(1e-6,1 - r_rc)) + 2*log(pmax(1e-6,r_rc))) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(log(pmax(1e-6,r_rc)) + log(5 - 11*r_rc + 12*r_rc^2 - 6*r_rc^3)) + (x[,"n_01"] + x[,"n_03"] + x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"] + x[,"n_41"] + x[,"n_43"])*log(r_rc - 3*r_rc^2 + 4*r_rc^3 - 2*r_rc^4) + x[,"n_22"]*log(4/9 - (8*r_rc)/9 + (17*r_rc^2)/9 - 2*r_rc^3 + r_rc^4) + (x[,"n_11"] + x[,"n_13"] + x[,"n_31"] + x[,"n_33"])*log(1 - 3*r_rc + 7*r_rc^2 - 8*r_rc^3 + 4*r_rc^4) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*log(1 - 4*r_rc + 10*r_rc^2 - 12*r_rc^3 + 6*r_rc^4)
  
  
  logL_rm <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22) {
    L <- (-3*n00 - 2*n01 - 3*n02 - 2*n03 - 3*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 3*n20 - 2*n21 - 2*n23 - 3*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34 - 3*n40 - 2*n41 - 3*n42 - 2*n43 - 3*n44)*log(2) + 2*(-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n23 - n24 - n30 - n31 - n32 - n33 - n34 - n40 - n41 - n42 - n43 - n44)*log(3) + (n04 + n40)*(3*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n00 + n44)*(log(pmax(1e-6,1 - r)) + 3*log(pmax(1e-6,r))) + (n01 + n10 + n34 + n43)*(2*log(pmax(1e-6,r)) + log(2 - 3*r + 2*r^2)) + (n03 + n14 + n30 + n41)*(2*log(pmax(1e-6,1 - r)) + log(1 - r + 2*r^2)) + (n02 + n20 + n24 + n42)*(log(pmax(1e-6,r)) + log(5 - 11*r + 12*r^2 - 6*r^3)) + (n11 + n33)*(log(pmax(1e-6,r)) + log(3 - 5*r + 7*r^2 - 4*r^3)) + (n13 + n31)*log(1 + 2*r - 8*r^2 + 9*r^3 - 4*r^4) + n22*log(2/9 + r/4 - (3*r^2)/4 + r^3 - r^4/2) + (n12 + n21 + n23 + n32)*log(5 - 7*r + 19*r^2 - 24*r^3 + 12*r^4)
    return(L)}
  interlogL_rm <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22) {
    optimize(logL_rm,c(0,0.5), n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n23,n24,n30,n31,n32,n33,n34,n40,n41,n42,n43,n44,n22, maximum=TRUE, lower=0, upper=0.5)$maximum}
  
  
  r_rm <- parallel::mcmapply(interlogL_rm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],x[,"n_40"],x[,"n_41"],x[,"n_42"],x[,"n_43"],x[,"n_44"],x[,"n_22"], mc.cores = ncores)
  
  
  LOD_rm <- -2*(x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(2) + 4*(x[,"n_04"] + x[,"n_40"])*log10(2) + 2*(x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*log10(2) + 2*(x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*log10(2) + 4*(x[,"n_00"] + x[,"n_44"])*log10(2) + (x[,"n_13"] + x[,"n_31"])*(3*log10(2) - log10(7)) - (x[,"n_11"] + x[,"n_33"])*(-3*log10(2) + log10(7)) - (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(-3*log10(2) + log10(7)) + x[,"n_22"]*(5*log10(2) + 2*log10(3) - log10(73)) + (x[,"n_04"] + x[,"n_40"])*(3*log10(pmax(1e-6,1 - r_rm)) + log10(pmax(1e-6,r_rm))) + (x[,"n_00"] + x[,"n_44"])*(log10(pmax(1e-6,1 - r_rm)) + 3*log10(pmax(1e-6,r_rm))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(2*log10(pmax(1e-6,r_rm)) + log10(2 - 3*r_rm + 2*r_rm^2)) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(2*log10(pmax(1e-6,1 - r_rm)) + log10(1 - r_rm + 2*r_rm^2)) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(log10(pmax(1e-6,r_rm)) + log10(5 - 11*r_rm + 12*r_rm^2 - 6*r_rm^3)) + (x[,"n_11"] + x[,"n_33"])*(log10(pmax(1e-6,r_rm)) + log10(3 - 5*r_rm + 7*r_rm^2 - 4*r_rm^3)) + (x[,"n_13"] + x[,"n_31"])*log10(1 + 2*r_rm - 8*r_rm^2 + 9*r_rm^3 - 4*r_rm^4) + x[,"n_22"]*log10(2/9 + r_rm/4 - (3*r_rm^2)/4 + r_rm^3 - r_rm^4/2) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(5 - 7*r_rm + 19*r_rm^2 - 24*r_rm^3 + 12*r_rm^4)
  
  
  logL_rm <- (-3*x[,"n_00"] - 2*x[,"n_01"] - 3*x[,"n_02"] - 2*x[,"n_03"] - 3*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 3*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_23"] - 3*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"] - 3*x[,"n_40"] - 2*x[,"n_41"] - 3*x[,"n_42"] - 2*x[,"n_43"] - 3*x[,"n_44"])*log(2) + 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"] - x[,"n_40"] - x[,"n_41"] - x[,"n_42"] - x[,"n_43"] - x[,"n_44"])*log(3) + (x[,"n_04"] + x[,"n_40"])*(3*log(pmax(1e-6,1 - r_rm)) + log(pmax(1e-6,r_rm))) + (x[,"n_00"] + x[,"n_44"])*(log(pmax(1e-6,1 - r_rm)) + 3*log(pmax(1e-6,r_rm))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(2*log(pmax(1e-6,r_rm)) + log(2 - 3*r_rm + 2*r_rm^2)) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(2*log(pmax(1e-6,1 - r_rm)) + log(1 - r_rm + 2*r_rm^2)) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(log(pmax(1e-6,r_rm)) + log(5 - 11*r_rm + 12*r_rm^2 - 6*r_rm^3)) + (x[,"n_11"] + x[,"n_33"])*(log(pmax(1e-6,r_rm)) + log(3 - 5*r_rm + 7*r_rm^2 - 4*r_rm^3)) + (x[,"n_13"] + x[,"n_31"])*log(1 + 2*r_rm - 8*r_rm^2 + 9*r_rm^3 - 4*r_rm^4) + x[,"n_22"]*log(2/9 + r_rm/4 - (3*r_rm^2)/4 + r_rm^3 - r_rm^4/2) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log(5 - 7*r_rm + 19*r_rm^2 - 24*r_rm^3 + 12*r_rm^4)
  
  
  logL_rr <- function(r,n00,n02,n04,n11,n13,n20,n24,n31,n33,n40,n42,n44,n01,n03,n10,n12,n14,n21,n23,n30,n32,n34,n41,n43,n22) {
    L <- (-2*n00 - n02 - 2*n04 + n11 + n13 - n20 - n24 + n31 + n33 - 2*n40 - n42 - 2*n44)*log(2) + (-2*n00 - 2*n01 - n02 - 2*n03 - 2*n04 - 2*n10 - 2*n11 - 2*n12 - 2*n13 - 2*n14 - n20 - 2*n21 - 2*n23 - n24 - 2*n30 - 2*n31 - 2*n32 - 2*n33 - 2*n34 - 2*n40 - 2*n41 - n42 - 2*n43 - 2*n44)*log(3) + 4*(n04 + n40)*log(pmax(1e-6,1 - r)) + 4*(n00 + n44)*log(pmax(1e-6,r)) + (n03 + n14 + n30 + n41)*(3*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n02 + n20 + n24 + n42)*(2*log(pmax(1e-6,1 - r)) + 2*log(pmax(1e-6,r))) + (n01 + n10 + n34 + n43)*(log(pmax(1e-6,1 - r)) + 3*log(pmax(1e-6,r))) + (n11 + n33)*(2*log(pmax(1e-6,r)) + log(2 - 3*r + 2*r^2)) + (n13 + n31)*(2*log(pmax(1e-6,1 - r)) + log(1 - r + 2*r^2)) + (n12 + n21 + n23 + n32)*(log(pmax(1e-6,r)) + log(5 - 11*r + 12*r^2 - 6*r^3)) + n22*log(1/2 - (10*r)/9 + (19*r^2)/9 - 2*r^3 + r^4)
    return(L)}
  interlogL_rr <- function(n00,n02,n04,n11,n13,n20,n24,n31,n33,n40,n42,n44,n01,n03,n10,n12,n14,n21,n23,n30,n32,n34,n41,n43,n22) {
    optimize(logL_rr,c(0,0.5), n00,n02,n04,n11,n13,n20,n24,n31,n33,n40,n42,n44,n01,n03,n10,n12,n14,n21,n23,n30,n32,n34,n41,n43,n22, maximum=TRUE, lower=0, upper=0.5)$maximum}
  
  
  r_rr <- parallel::mcmapply(interlogL_rr,x[,"n_00"],x[,"n_02"],x[,"n_04"],x[,"n_11"],x[,"n_13"],x[,"n_20"],x[,"n_24"],x[,"n_31"],x[,"n_33"],x[,"n_40"],x[,"n_42"],x[,"n_44"],x[,"n_01"],x[,"n_03"],x[,"n_10"],x[,"n_12"],x[,"n_14"],x[,"n_21"],x[,"n_23"],x[,"n_30"],x[,"n_32"],x[,"n_34"],x[,"n_41"],x[,"n_43"],x[,"n_22"], mc.cores = ncores)
  
  
  LOD_rr <- 2*(x[,"n_13"] + x[,"n_31"])*log10(2) + 2*(x[,"n_11"] + x[,"n_33"])*log10(2) + 4*(x[,"n_04"] + x[,"n_40"])*log10(2) + 4*(x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*log10(2) + 4*(x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*log10(2) + 4*(x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*log10(2) + 4*(x[,"n_00"] + x[,"n_44"])*log10(2) - (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(-3*log10(2) + log10(7)) + x[,"n_22"]*(4*log10(2) + 2*log10(3) - log10(41)) + 4*(x[,"n_04"] + x[,"n_40"])*log10(pmax(1e-6,1 - r_rr)) + 4*(x[,"n_00"] + x[,"n_44"])*log10(pmax(1e-6,r_rr)) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(3*log10(pmax(1e-6,1 - r_rr)) + log10(pmax(1e-6,r_rr))) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(2*log10(pmax(1e-6,1 - r_rr)) + 2*log10(pmax(1e-6,r_rr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(log10(pmax(1e-6,1 - r_rr)) + 3*log10(pmax(1e-6,r_rr))) + (x[,"n_11"] + x[,"n_33"])*(2*log10(pmax(1e-6,r_rr)) + log10(2 - 3*r_rr + 2*r_rr^2)) + (x[,"n_13"] + x[,"n_31"])*(2*log10(pmax(1e-6,1 - r_rr)) + log10(1 - r_rr + 2*r_rr^2)) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(log10(pmax(1e-6,r_rr)) + log10(5 - 11*r_rr + 12*r_rr^2 - 6*r_rr^3)) + x[,"n_22"]*log10(1/2 - (10*r_rr)/9 + (19*r_rr^2)/9 - 2*r_rr^3 + r_rr^4)
  
  
  logL_rr <- (-2*x[,"n_00"] - x[,"n_02"] - 2*x[,"n_04"] + x[,"n_11"] + x[,"n_13"] - x[,"n_20"] - x[,"n_24"] + x[,"n_31"] + x[,"n_33"] - 2*x[,"n_40"] - x[,"n_42"] - 2*x[,"n_44"])*log(2) + (-2*x[,"n_00"] - 2*x[,"n_01"] - x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_14"] - x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_23"] - x[,"n_24"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"] - 2*x[,"n_34"] - 2*x[,"n_40"] - 2*x[,"n_41"] - x[,"n_42"] - 2*x[,"n_43"] - 2*x[,"n_44"])*log(3) + 4*(x[,"n_04"] + x[,"n_40"])*log(pmax(1e-6,1 - r_rr)) + 4*(x[,"n_00"] + x[,"n_44"])*log(pmax(1e-6,r_rr)) + (x[,"n_03"] + x[,"n_14"] + x[,"n_30"] + x[,"n_41"])*(3*log(pmax(1e-6,1 - r_rr)) + log(pmax(1e-6,r_rr))) + (x[,"n_02"] + x[,"n_20"] + x[,"n_24"] + x[,"n_42"])*(2*log(pmax(1e-6,1 - r_rr)) + 2*log(pmax(1e-6,r_rr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_34"] + x[,"n_43"])*(log(pmax(1e-6,1 - r_rr)) + 3*log(pmax(1e-6,r_rr))) + (x[,"n_11"] + x[,"n_33"])*(2*log(pmax(1e-6,r_rr)) + log(2 - 3*r_rr + 2*r_rr^2)) + (x[,"n_13"] + x[,"n_31"])*(2*log(pmax(1e-6,1 - r_rr)) + log(1 - r_rr + 2*r_rr^2)) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(log(pmax(1e-6,r_rr)) + log(5 - 11*r_rr + 12*r_rr^2 - 6*r_rr^3)) + x[,"n_22"]*log(1/2 - (10*r_rr)/9 + (19*r_rr^2)/9 - 2*r_rr^3 + r_rr^4)
  
  
  return(list(
    r_mat = cbind(r_cc,r_cm,r_cr,r_mc,r_mm,r_mr,r_rc,r_rm,r_rr,0.499),
    LOD_mat = cbind(LOD_cc,LOD_cm,LOD_cr,LOD_mc,LOD_mm,LOD_mr,LOD_rc,LOD_rm,LOD_rr,0),
    logL_mat = cbind(logL_cc,logL_cm,logL_cr,logL_mc,logL_mm,logL_mr,logL_rc,logL_rm,logL_rr,-1e6),
    phasing_strategy = "MLL",
    possible_phases = c("coupling coupling","coupling mixed","coupling repulsion","mixed coupling","mixed mixed","mixed repulsion","repulsion coupling","repulsion mixed","repulsion repulsion","unknown")
  )
  )
}

#' @rdname r4_functions
#' @noRd
r4_2.0_1.2<-function(x, ncores=1){
  logLc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){ 
    L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - n10 - n11 - n12 - n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + 
      (n11 + n12)*log(5) + (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(9) + 
      (n01 + n22)*log(4 - 3*r) + (n00 + n23)*log(1 - r) + (n03 + n20)*log(r) + (n02 + n21)*log(1 + 3*r)
    return(L)}
  
  inter_logLc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){
    optimize(logLc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION
  ############################
  logLr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){ 
    L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - n10 - n11 - n12 - n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + 
      (n11 + n12)*log(5) + (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(9) + 
      (n02 + n21)*log(4 - 3*r) + (n03 + n20)*log(1 - r) + (n00 + n23)*log(r) + (n01 + n22)*log(1 + 3*r)
    return(L)}
  
  inter_logLr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){
    optimize(logLr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_c <- parallel::mcmapply(inter_logLc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                                mc.cores = ncores)
  
  r_r <- parallel::mcmapply(inter_logLr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                                 mc.cores = ncores)
  
  LOD_c <- log10(((4 - 3*r_c)^(x[,"n_01"] + x[,"n_22"])*(1 - r_c)^(x[,"n_00"] + x[,"n_23"])*r_c^(x[,"n_03"] + x[,"n_20"])*(1 + 3*r_c)^(x[,"n_02"] + x[,"n_21"]))/((4 - 3*0.5)^(x[,"n_01"] + x[,"n_22"])*(1 - 0.5)^(x[,"n_00"] + x[,"n_23"])*0.5^(x[,"n_03"] + x[,"n_20"])*(1 + 3*0.5)^(x[,"n_02"] + x[,"n_21"])))
  LOD_r <- log10(((4 - 3*r_r)^(x[,"n_02"] + x[,"n_21"])*(1 - r_r)^(x[,"n_03"] + x[,"n_20"])*r_r^(x[,"n_00"] + x[,"n_23"])*(1 + 3*r_r)^(x[,"n_01"] + x[,"n_22"]))/((4 - 3*0.5)^(x[,"n_02"] + x[,"n_21"])*(1 - 0.5)^(x[,"n_03"] + x[,"n_20"])*0.5^(x[,"n_00"] + x[,"n_23"])*(1 + 3*0.5)^(x[,"n_01"] + x[,"n_22"])))
  
  return(list(r_mat=cbind(r_c, r_r),
              LOD_mat=cbind(LOD_c, LOD_r),
              logL_mat=NULL,
              phasing_strategy="MINR", 
              possible_phases=c("coupling",
                                "repulsion")))
  
}

#' @rdname r4_functions
#' @noRd
r4_2.0_1.3<-function(x, ncores=1){
  common_sum <- x[,"n_01"] + x[,"n_03"] + x[,"n_21"] + x[,"n_23"]
  
  r_c <- (x[,"n_03"] + x[,"n_21"])/common_sum
  r_r <- (x[,"n_01"] + x[,"n_23"])/common_sum
  
  ############################################################################
  ## These are not the true likelihood functions, as I have removed the parts that are
  ## independent of r, which cancel in the likelihood ratio anyway...
  Lc <- function(r) {
    L <- (1 - r)^(x[,"n_01"] + x[,"n_23"])*r^(x[,"n_03"] + x[,"n_21"])
    return(L)}
  
  Lr <- function(r) {
    L <- (1 - r)^(x[,"n_03"] + x[,"n_21"])*r^(x[,"n_01"] + x[,"n_23"])
    return(L)}
  
  LOD_c <- log10(Lc(r_c)/Lc(0.5))
  LOD_r <- log10(Lr(r_r)/Lr(0.5))
  
  return(list(r_mat=cbind(r_c, r_r),
              LOD_mat=cbind(LOD_c, LOD_r),
              logL_mat=NULL,
              phasing_strategy="MINR", 
              possible_phases=c("coupling",
                                "repulsion")))
}

#' @rdname r4_functions
#' @noRd
r4_2.0_2.2<-function(x, ncores=1){
  logLc <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - n10 - n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24)*log(9) + 
      (n00 + n24)*log((-1 + r)^2) + (n10 + n14)*log(-((-1 + r)*r)) + (n04 + n20)*log(r^2) + 
      (n03 + n21)*log(r*(1 + r)) + (n02 + n22)*log(1 + 6*r - 6*r^2) + (n11 + n13)*log(1 + r - r^2) + 
      (n01 + n23)*log(2 - 3*r + r^2) + n12*log(4 - 3*r + 3*r^2)
    return(L)}
  
  inter_logLc <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){
    optimize(logLc,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  
  ############################
  ## MIXED
  ############################
  logLm <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){ 
    L <- (-3*n00 - 2*n01 - 2*n02 - 2*n03 - 3*n04 - 2*n10 - n11 - n12 - n13 - 2*n14 - 3*n20 - 2*n21 - 2*n22 - 2*n23 - 3*n24)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24)*log(9) + 
      (n00 + n04 + n20 + n24)*log(-((-1 + r)*r)) + n12*log(5 + 3*r - 3*r^2) + (n01 + n03 + n21 + n23)*log(1 + r - r^2) + 
      (n10 + n14)*log(1 - r + r^2) + (n11 + n13)*log(3 - r + r^2) + (n02 + n22)*log(4 - 3*r + 3*r^2)
    return(L)}
  
  inter_logLm <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){
    optimize(logLm,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION
  ############################
  logLr <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - n10 - n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24)*log(9) + 
      (n04 + n20)*log((-1 + r)^2) + (n10 + n14)*log(-((-1 + r)*r)) + (n00 + n24)*log(r^2) + (n01 + n23)*log(r*(1 + r)) + 
      (n02 + n22)*log(1 + 6*r - 6*r^2) + (n11 + n13)*log(1 + r - r^2) + (n03 + n21)*log(2 - 3*r + r^2) + n12*log(4 - 3*r + 3*r^2)
    return(L)}
  
  inter_logLr <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){
    optimize(logLr,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_c <- parallel::mcmapply(inter_logLc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],
                                mc.cores = ncores)
  
  r_m <- parallel::mcmapply(inter_logLm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],
                                mc.cores = ncores)
  
  r_r <- parallel::mcmapply(inter_logLr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],
                                      mc.cores = ncores)
  
  LOD_c <- (x[,"n_00"] + x[,"n_24"])*log10((-1 + r_c)^2) + (x[,"n_10"] + x[,"n_14"])*log10(-((-1 + r_c)*r_c)) + (x[,"n_04"] + x[,"n_20"])*log10(r_c^2) + 
    (x[,"n_03"] + x[,"n_21"])*log10(r_c*(1 + r_c)) + (x[,"n_02"] + x[,"n_22"])*log10(1 + 6*r_c - 6*r_c^2) + (x[,"n_11"] + x[,"n_13"])*log10(1 + r_c - r_c^2) + 
    (x[,"n_01"] + x[,"n_23"])*log10(2 - 3*r_c + r_c^2) + x[,"n_12"]*log10(4 - 3*r_c + 3*r_c^2) - 
    (x[,"n_00"] + x[,"n_24"])*log10((-1 + 0.5)^2) - (x[,"n_10"] + x[,"n_14"])*log10(-((-1 + 0.5)*0.5)) - (x[,"n_04"] + x[,"n_20"])*log10(0.5^2) - 
    (x[,"n_03"] + x[,"n_21"])*log10(0.5*(1 + 0.5)) - (x[,"n_02"] + x[,"n_22"])*log10(1 + 6*0.5 - 6*0.5^2) - (x[,"n_11"] + x[,"n_13"])*log10(1 + 0.5 - 0.5^2) - 
    (x[,"n_01"] + x[,"n_23"])*log10(2 - 3*0.5 + 0.5^2) - x[,"n_12"]*log10(4 - 3*0.5 + 3*0.5^2)
  
  LOD_m <- (x[,"n_00"] + x[,"n_04"] + x[,"n_20"] + x[,"n_24"])*log10(-((-1 + r_m)*r_m)) + x[,"n_12"]*log10(5 + 3*r_m - 3*r_m^2) + (x[,"n_01"] + x[,"n_03"] + x[,"n_21"] + x[,"n_23"])*log10(1 + r_m - r_m^2) + 
    (x[,"n_10"] + x[,"n_14"])*log10(1 - r_m + r_m^2) + (x[,"n_11"] + x[,"n_13"])*log10(3 - r_m + r_m^2) + (x[,"n_02"] + x[,"n_22"])*log10(4 - 3*r_m + 3*r_m^2) - 
    (x[,"n_00"] + x[,"n_04"] + x[,"n_20"] + x[,"n_24"])*log10(-((-1 + 0.5)*0.5)) - x[,"n_12"]*log10(5 + 3*0.5 - 3*0.5^2) - (x[,"n_01"] + x[,"n_03"] + x[,"n_21"] + x[,"n_23"])*log10(1 + 0.5 - 0.5^2) - 
    (x[,"n_10"] + x[,"n_14"])*log10(1 - 0.5 + 0.5^2) - (x[,"n_11"] + x[,"n_13"])*log10(3 - 0.5 + 0.5^2) - (x[,"n_02"] + x[,"n_22"])*log10(4 - 3*0.5 + 3*0.5^2)
  
  LOD_r <- (x[,"n_04"] + x[,"n_20"])*log10((-1 + r_r)^2) + (x[,"n_10"] + x[,"n_14"])*log10(-((-1 + r_r)*r_r)) + (x[,"n_00"] + x[,"n_24"])*log10(r_r^2) + (x[,"n_01"] + x[,"n_23"])*log10(r_r*(1 + r_r)) + 
    (x[,"n_02"] + x[,"n_22"])*log10(1 + 6*r_r - 6*r_r^2) + (x[,"n_11"] + x[,"n_13"])*log10(1 + r_r - r_r^2) + (x[,"n_03"] + x[,"n_21"])*log10(2 - 3*r_r + r_r^2) + x[,"n_12"]*log10(4 - 3*r_r + 3*r_r^2) - 
    (x[,"n_04"] + x[,"n_20"])*log10((-1 + 0.5)^2) - (x[,"n_10"] + x[,"n_14"])*log10(-((-1 + 0.5)*0.5)) - (x[,"n_00"] + x[,"n_24"])*log10(0.5^2) - (x[,"n_01"] + x[,"n_23"])*log10(0.5*(1 + 0.5)) -
    (x[,"n_02"] + x[,"n_22"])*log10(1 + 6*0.5 - 6*0.5^2) - (x[,"n_11"] + x[,"n_13"])*log10(1 + 0.5 - 0.5^2) - (x[,"n_03"] + x[,"n_21"])*log10(2 - 3*0.5 + 0.5^2) - x[,"n_12"]*log10(4 - 3*0.5 + 3*0.5^2)
  
  logL_c <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - x[,"n_10"] - x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"])*log(9) + 
    (x[,"n_00"] + x[,"n_24"])*log((-1 + r_c)^2) + (x[,"n_10"] + x[,"n_14"])*log(-((-1 + r_c)*r_c)) + (x[,"n_04"] + x[,"n_20"])*log(r_c^2) + 
    (x[,"n_03"] + x[,"n_21"])*log(r_c*(1 + r_c)) + (x[,"n_02"] + x[,"n_22"])*log(1 + 6*r_c - 6*r_c^2) + (x[,"n_11"] + x[,"n_13"])*log(1 + r_c - r_c^2) + 
    (x[,"n_01"] + x[,"n_23"])*log(2 - 3*r_c + r_c^2) + x[,"n_12"]*log(4 - 3*r_c + 3*r_c^2)
  
  logL_m <- (-3*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - 3*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 3*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"] - 3*x[,"n_24"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"])*log(9) + 
    (x[,"n_00"] + x[,"n_04"] + x[,"n_20"] + x[,"n_24"])*log(-((-1 + r_m)*r_m)) + x[,"n_12"]*log(5 + 3*r_m - 3*r_m^2) + (x[,"n_01"] + x[,"n_03"] + x[,"n_21"] + x[,"n_23"])*log(1 + r_m - r_m^2) + 
    (x[,"n_10"] + x[,"n_14"])*log(1 - r_m + r_m^2) + (x[,"n_11"] + x[,"n_13"])*log(3 - r_m + r_m^2) + (x[,"n_02"] + x[,"n_22"])*log(4 - 3*r_m + 3*r_m^2)
  
  logL_r <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - x[,"n_10"] - x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"])*log(9) + 
    (x[,"n_04"] + x[,"n_20"])*log((-1 + r_r)^2) + (x[,"n_10"] + x[,"n_14"])*log(-((-1 + r_r)*r_r)) + (x[,"n_00"] + x[,"n_24"])*log(r_r^2) + (x[,"n_01"] + x[,"n_23"])*log(r_r*(1 + r_r)) + 
    (x[,"n_02"] + x[,"n_22"])*log(1 + 6*r_r - 6*r_r^2) + (x[,"n_11"] + x[,"n_13"])*log(1 + r_r - r_r^2) + (x[,"n_03"] + x[,"n_21"])*log(2 - 3*r_r + r_r^2) + x[,"n_12"]*log(4 - 3*r_r + 3*r_r^2)
  
  return(list(r_mat=cbind(r_c, r_m, r_r),
              LOD_mat=cbind(LOD_c, LOD_m, LOD_r),
              logL_mat=cbind(logL_c, logL_m, logL_r),
              phasing_strategy="MLL", 
              possible_phases=c("coupling",
                                "mixed",
                                "repulsion")))
  
}

#' @rdname r4_functions
#' @noRd
r4_1.1_1.2<-function(x, ncores=1){
  ############################
  ## COUPLING COUPLING
  ############################
  logLcc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){ 
    L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - n10 - 2*n11 - 2*n12 - n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(3) + 
      (n00 + n23)*log((-1 + r)^2) + (n02 + n21)*log(-((-3 + r)*r)) + (n10 + n13)*log(-((-1 + r)*r)) + 
      (n03 + n20)*log(r^2) + (n01 + n22)*log(2 - r - r^2) + (n11 + n12)*log(3 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLcc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){
    optimize(logLcc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  
  ############################
  ## COUPLING REPULSION
  ############################
  logLcr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){ 
    L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - 2*n10 - n11 - n12 - 2*n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(3) + 
      (n00 + n03 + n20 + n23)*log(-((-1 + r)*r)) + (n11 + n12)*log(1 + r - r^2) + (n02 + n21)*log(1 + r^2) + 
      (n01 + n22)*log(2 - 2*r + r^2) + (n10 + n13)*log(1 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLcr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){
    optimize(logLcr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION COUPLING
  ############################
  logLrc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){ 
    L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - n10 - 2*n11 - 2*n12 - n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(9) + 
      (n03 + n20)*log(-((-2 + r)*r)) + (n00 + n23)*log(-((-1 + r)*(1 + r))) + (n11 + n12)*log(7 + 2*r - 2*r^2) + 
      (n10 + n13)*log(1 - r + r^2) + (n01 + n02 + n21 + n22)*log(4 - r + r^2)
    return(L)}
  
  inter_logLrc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){
    optimize(logLrc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION REPULSION
  ############################
  logLrr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){ 
    L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - 2*n10 - n11 - n12 - 2*n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(9) + 
      (n03 + n20)*log((-2 + r)*(-1 + r)) + (n00 + n23)*log(r*(1 + r)) + (n10 + n13)*log(1 + 2*r - 2*r^2) +
      (n02 + n21)*log(5 - 2*r - r^2) + (n01 + n22)*log(2 + 4*r - r^2) + (n11 + n12)*log(4 - r + r^2)
    return(L)}
  
  inter_logLrr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23){
    optimize(logLrr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_cc <- parallel::mcmapply(inter_logLcc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                             mc.cores = ncores)
  
  r_cr <-  parallel::mcmapply(inter_logLcr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                              mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                              mc.cores = ncores)
  
  r_rr <- parallel::mcmapply(inter_logLrr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                             mc.cores = ncores)
  
  LOD_cc <- (x[,"n_00"] + x[,"n_23"])*log10((-1 + r_cc)^2) + (x[,"n_02"] + x[,"n_21"])*log10(-((-3 + r_cc)*r_cc)) + (x[,"n_10"] + x[,"n_13"])*log10(-((-1 + r_cc)*r_cc)) + 
    (x[,"n_03"] + x[,"n_20"])*log10(r_cc^2) + (x[,"n_01"] + x[,"n_22"])*log10(2 - r_cc - r_cc^2) + (x[,"n_11"] + x[,"n_12"])*log10(3 - 2*r_cc + 2*r_cc^2) - 
    (x[,"n_00"] + x[,"n_23"])*log10((-1 + 0.5)^2) - (x[,"n_02"] + x[,"n_21"])*log10(-((-3 + 0.5)*0.5)) - (x[,"n_10"] + x[,"n_13"])*log10(-((-1 + 0.5)*0.5)) - 
    (x[,"n_03"] + x[,"n_20"])*log10(0.5^2) - (x[,"n_01"] + x[,"n_22"])*log10(2 - 0.5 - 0.5^2) - (x[,"n_11"] + x[,"n_12"])*log10(3 - 2*0.5 + 2*0.5^2)
  
  LOD_cr <- (x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*log10(-((-1 + r_cr)*r_cr)) + (x[,"n_11"] + x[,"n_12"])*log10(1 + r_cr - r_cr^2) + (x[,"n_02"] + x[,"n_21"])*log10(1 + r_cr^2) + 
    (x[,"n_01"] + x[,"n_22"])*log10(2 - 2*r_cr + r_cr^2) + (x[,"n_10"] + x[,"n_13"])*log10(1 - 2*r_cr + 2*r_cr^2) - (x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*log10(-((-1 + 0.5)*0.5)) - 
    (x[,"n_11"] + x[,"n_12"])*log10(1 + 0.5 - 0.5^2) - (x[,"n_02"] + x[,"n_21"])*log10(1 + 0.5^2) - (x[,"n_01"] + x[,"n_22"])*log10(2 - 2*0.5 + 0.5^2) - (x[,"n_10"] + x[,"n_13"])*log10(1 - 2*0.5 + 2*0.5^2)
  
  LOD_rc <- (x[,"n_03"] + x[,"n_20"])*log10(-((-2 + r_rc)*r_rc)) + (x[,"n_00"] + x[,"n_23"])*log10(-((-1 + r_rc)*(1 + r_rc))) + (x[,"n_11"] + x[,"n_12"])*log10(7 + 2*r_rc - 2*r_rc^2) + 
    (x[,"n_10"] + x[,"n_13"])*log10(1 - r_rc + r_rc^2) + (x[,"n_01"] + x[,"n_02"] + x[,"n_21"] + x[,"n_22"])*log10(4 - r_rc + r_rc^2) - (x[,"n_03"] + x[,"n_20"])*log10(-((-2 + 0.5)*0.5)) -
    (x[,"n_00"] + x[,"n_23"])*log10(-((-1 + 0.5)*(1 + 0.5))) - (x[,"n_11"] + x[,"n_12"])*log10(7 + 2*0.5 - 2*0.5^2) - (x[,"n_10"] + x[,"n_13"])*log10(1 - 0.5 + 0.5^2) - (x[,"n_01"] + x[,"n_02"] + x[,"n_21"] + x[,"n_22"])*log10(4 - 0.5 + 0.5^2)
  
  LOD_rr <- (x[,"n_03"] + x[,"n_20"])*log10((-2 + r_rr)*(-1 + r_rr)) + (x[,"n_00"] + x[,"n_23"])*log10(r_rr*(1 + r_rr)) + (x[,"n_10"] + x[,"n_13"])*log10(1 + 2*r_rr - 2*r_rr^2) + 
    (x[,"n_02"] + x[,"n_21"])*log10(5 - 2*r_rr - r_rr^2) + (x[,"n_01"] + x[,"n_22"])*log10(2 + 4*r_rr - r_rr^2) + (x[,"n_11"] + x[,"n_12"])*log10(4 - r_rr + r_rr^2) - 
    (x[,"n_03"] + x[,"n_20"])*log10((-2 + 0.5)*(-1 + 0.5)) - (x[,"n_00"] + x[,"n_23"])*log10(0.5*(1 + 0.5)) - (x[,"n_10"] + x[,"n_13"])*log10(1 + 2*0.5 - 2*0.5^2) - 
    (x[,"n_02"] + x[,"n_21"])*log10(5 - 2*0.5 - 0.5^2) - (x[,"n_01"] + x[,"n_22"])*log10(2 + 4*0.5 - 0.5^2) - (x[,"n_11"] + x[,"n_12"])*log10(4 - 0.5 + 0.5^2)
  
  
  ## Record the log likelihoods also:
  logL_cc <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(3) + 
    (x[,"n_00"] + x[,"n_23"])*log((-1 + r_cc)^2) + (x[,"n_02"] + x[,"n_21"])*log(-((-3 + r_cc)*r_cc)) + (x[,"n_10"] + x[,"n_13"])*log(-((-1 + r_cc)*r_cc)) + 
    (x[,"n_03"] + x[,"n_20"])*log(r_cc^2) + (x[,"n_01"] + x[,"n_22"])*log(2 - r_cc - r_cc^2) + (x[,"n_11"] + x[,"n_12"])*log(3 - 2*r_cc + 2*r_cc^2)
  
  logL_cr <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(3) + 
    (x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*log(-((-1 + r_cr)*r_cr)) + (x[,"n_11"] + x[,"n_12"])*log(1 + r_cr - r_cr^2) + (x[,"n_02"] + x[,"n_21"])*log(1 + r_cr^2) + 
    (x[,"n_01"] + x[,"n_22"])*log(2 - 2*r_cr + r_cr^2) + (x[,"n_10"] + x[,"n_13"])*log(1 - 2*r_cr + 2*r_cr^2)
  
  logL_rc <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(9) + 
    (x[,"n_03"] + x[,"n_20"])*log(-((-2 + r_rc)*r_rc)) + (x[,"n_00"] + x[,"n_23"])*log(-((-1 + r_rc)*(1 + r_rc))) + (x[,"n_11"] + x[,"n_12"])*log(7 + 2*r_rc - 2*r_rc^2) + 
    (x[,"n_10"] + x[,"n_13"])*log(1 - r_rc + r_rc^2) + (x[,"n_01"] + x[,"n_02"] + x[,"n_21"] + x[,"n_22"])*log(4 - r_rc + r_rc^2)
  
  logL_rr <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(9) + 
    (x[,"n_03"] + x[,"n_20"])*log((-2 + r_rr)*(-1 + r_rr)) + (x[,"n_00"] + x[,"n_23"])*log(r_rr*(1 + r_rr)) + (x[,"n_10"] + x[,"n_13"])*log(1 + 2*r_rr - 2*r_rr^2) +
    (x[,"n_02"] + x[,"n_21"])*log(5 - 2*r_rr - r_rr^2) + (x[,"n_01"] + x[,"n_22"])*log(2 + 4*r_rr - r_rr^2) + (x[,"n_11"] + x[,"n_12"])*log(4 - r_rr + r_rr^2)
  
  return(list(
    r_mat=cbind(r_cc, r_cr, r_rc, r_rr),
    LOD_mat=cbind(LOD_cc, LOD_cr, LOD_rc, LOD_rr),
    logL_mat=cbind(logL_cc, logL_cr, logL_rc, logL_rr),
    phasing_strategy="MLL", 
    possible_phases=c("coupling coupling",
                      "coupling repulsion",
                      "repulsion coupling",
                      "repulsion repulsion")
  ))
}

#' @rdname r4_functions
#' @noRd
r4_1.1_1.3<-function(x, ncores=1){
  ############################
  ## COUPLING COUPLING
  ############################
  logLcc <- function(r,n01,n02,n03,n11,n12,n13,n21,n22,n23){ 
    L <- (-n01 - n02 - n03 - n11 - n12 - n13 - n21 - n22 - n23)*log(3) + 
      (-n01 - n02 - n03 - n11 - n13 - n21 - n22 - n23)*log(4) + 
      (n01 + n23)*log((-2 + r)*(-1 + r)) + (n03 + n21)*log(r*(1 + r)) + 
      (n02 + n11 + n13 + n22)*log(1 + 2*r - 2*r^2) + n12*log(1 - r + r^2)
    return(L)}
  
  inter_logLcc <- function(n01,n02,n03,n11,n12,n13,n21,n22,n23){
    optimize(logLcc,c(0,0.5),n01,n02,n03,n11,n12,n13,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcc <- function(r){
  #     L <- 3^(-x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*4^(-x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_11"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*((-2 + r)*(-1 + r))^(x[,"n_01"] + x[,"n_23"])*(r*(1 + r))^(x[,"n_03"] + x[,"n_21"])*(1 + 2*r - 2*r^2)^(x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*(1 - r + r^2)^x[,"n_12"]
  #     return(L)}
  
  ############################
  ## MIXED COUPLING REPULSION
  ############################
  logLcr <- function(r,n01,n02,n03,n11,n12,n13,n21,n22,n23){ 
    L <- (-n01 - n02 - n03 - n11 - n13 - n21 - n22 - n23)*log(4) + 
      (n01 + n03 + n12 + n21 + n23)*log(-((-1 + r)*r)) + 
      (n02 + n11 + n13 + n22)*log(1 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLcr <- function(n01,n02,n03,n11,n12,n13,n21,n22,n23){
    optimize(logLcr,c(0,0.5),n01,n02,n03,n11,n12,n13,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## MIXED REPULSION COUPLING
  ############################
  logLrc <- function(r,n01,n02,n03,n11,n12,n13,n21,n22,n23){ 
    L <- (-n01 - n02 - n03 - n11 - n13 - n21 - n22 - n23)*log(4) + 
      (-n01 - n02 - n03 - n11 - n12 - n13 - n21 - n22 - n23)*log(9) + 
      (n01 + n03 + n12 + n21 + n23)*log(-((-2 + r)*(1 + r))) + (n02 + n11 + n13 + n22)*log(5 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLrc <- function(n01,n02,n03,n11,n12,n13,n21,n22,n23){
    optimize(logLrc,c(0,0.5),n01,n02,n03,n11,n12,n13,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION
  ############################
  logLrr <- function(r,n01,n02,n03,n11,n12,n13,n21,n22,n23){ 
    L <- (-n01 - n02 - n03 - n11 - n12 - n13 - n21 - n22 - n23)*log(3) + 
      (-n01 - n02 - n03 - n11 - n13 - n21 - n22 - n23)*log(4) + (n03 + n21)*log((-2 + r)*(-1 + r)) + 
      (n01 + n23)*log(r*(1 + r)) + (n02 + n11 + n13 + n22)*log(1 + 2*r - 2*r^2) + n12*log(1 - r + r^2)
    return(L)}
  
  inter_logLrr <- function(n01,n02,n03,n11,n12,n13,n21,n22,n23){
    optimize(logLrr,c(0,0.5),n01,n02,n03,n11,n12,n13,n21,n22,n23,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_cc <- parallel::mcmapply(inter_logLcc,x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                             mc.cores = ncores)
  
  r_cr <-  parallel::mcmapply(inter_logLcr,x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                              mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc,x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                              mc.cores = ncores)
  
  r_rr <- parallel::mcmapply(inter_logLrr,x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_21"],x[,"n_22"],x[,"n_23"],
                             mc.cores = ncores)
  
  LOD_cc <- (x[,"n_01"] + x[,"n_23"])*log10((-2 + r_cc)*(-1 + r_cc)) + (x[,"n_03"] + x[,"n_21"])*log10(r_cc*(1 + r_cc)) + (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log10(1 + 2*r_cc - 2*r_cc^2) + 
    x[,"n_12"]*log10(1 - r_cc + r_cc^2) - (x[,"n_01"] + x[,"n_23"])*log10((-2 + 0.5)*(-1 + 0.5)) - (x[,"n_03"] + x[,"n_21"])*log10(0.5*(1 + 0.5)) - 
    (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log10(1 + 2*0.5 - 2*0.5^2) - x[,"n_12"]*log10(1 - 0.5 + 0.5^2)
  
  LOD_cr <- (x[,"n_01"] + x[,"n_03"] + x[,"n_12"] + x[,"n_21"] + x[,"n_23"])*log10(-((-1 + r_cr)*r_cr)) + (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log10(1 - 2*r_cr + 2*r_cr^2) - 
    (x[,"n_01"] + x[,"n_03"] + x[,"n_12"] + x[,"n_21"] + x[,"n_23"])*log10(-((-1 + 0.5)*0.5)) - (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log10(1 - 2*0.5 + 2*0.5^2)
  
  # LOD_rc <- (x[,"n_01"] + x[,"n_03"] + x[,"n_12"] + x[,"n_21"] + x[,"n_23"])*log10(-((-2 + r_rc)*(1 + r_rc))) + (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log10(5 - 2*r_rc + 2*r_rc^2) - 
  #   (x[,"n_01"] + x[,"n_03"] + x[,"n_12"] + x[,"n_21"] + x[,"n_23"])*log10(-((-2 + 0.5)*(1 + 0.5))) - (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log10(5 - 2*0.5 + 2*0.5^2)
  
  LOD_rr <- (x[,"n_03"] + x[,"n_21"])*log10((-2 + r_rr)*(-1 + r_rr)) + (x[,"n_01"] + x[,"n_23"])*log10(r_rr*(1 + r_rr)) + (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log10(1 + 2*r_rr - 2*r_rr^2) + 
    x[,"n_12"]*log10(1 - r_rr + r_rr^2) - (x[,"n_03"] + x[,"n_21"])*log10((-2 + 0.5)*(-1 + 0.5)) - (x[,"n_01"] + x[,"n_23"])*log10(0.5*(1 + 0.5)) - 
    (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log10(1 + 2*0.5 - 2*0.5^2) - x[,"n_12"]*log10(1 - 0.5 + 0.5^2)
  
  ## Record the logL also:
  logL_cc <- (-x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(3) + 
    (-x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_11"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(4) + 
    (x[,"n_01"] + x[,"n_23"])*log((-2 + r_cc)*(-1 + r_cc)) + (x[,"n_03"] + x[,"n_21"])*log(r_cc*(1 + r_cc)) + 
    (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log(1 + 2*r_cc - 2*r_cc^2) + x[,"n_12"]*log(1 - r_cc + r_cc^2)
  
  logL_cr <- (-x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_11"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(4) + 
    (x[,"n_01"] + x[,"n_03"] + x[,"n_12"] + x[,"n_21"] + x[,"n_23"])*log(-((-1 + r_cr)*r_cr)) + 
    (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log(1 - 2*r_cr + 2*r_cr^2)
  
  logL_rc <- (-x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_11"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(4) + 
    (-x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(9) + 
    (x[,"n_01"] + x[,"n_03"] + x[,"n_12"] + x[,"n_21"] + x[,"n_23"])*log(-((-2 + r_rc)*(1 + r_rc))) + (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log(5 - 2*r_rc + 2*r_rc^2)
  
  logL_rr <- (-x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(3) + 
    (-x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_11"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(4) + (x[,"n_03"] + x[,"n_21"])*log((-2 + r_rr)*(-1 + r_rr)) + 
    (x[,"n_01"] + x[,"n_23"])*log(r_rr*(1 + r_rr)) + (x[,"n_02"] + x[,"n_11"] + x[,"n_13"] + x[,"n_22"])*log(1 + 2*r_rr - 2*r_rr^2) + x[,"n_12"]*log(1 - r_rr + r_rr^2)
  
  r_rc<-0.499 #This is because the r estimate in this phase is very suspect
  LOD_rc<-0 #This is because the r estimate in this phase is very suspect
  
  logL_cr[r_cr>0.2]<-NA #logL estimates in this case are unreliable
  
  return(list(
    r_mat=cbind(r_cc, r_cr, r_rc, r_rr),
    LOD_mat=cbind(LOD_cc, LOD_cr, LOD_rc, LOD_rr),
    logL_mat=cbind(logL_cc, logL_cr, logL_rc, logL_rr),
    phasing_strategy="MLL", 
    possible_phases=c("coupling coupling",
                      "coupling repulsion",
                      "repulsion coupling",
                      "repulsion repulsion")
  ))
}

#' @rdname r4_functions
#' @noRd
r4_1.1_2.2<-function(x, ncores=1){
  ############################
  ## COUPLING COUPLING
  ############################
  logLcc <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){ 
    L <- (-2*n00 - n02 - 2*n04 - n10 - n12 - n14 - 2*n20 - n22 - 2*n24)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24)*log(9) + 
      (n01 + n23)*log(1 - r) + (n00 + n24)*log((-1 + r)^2) + (n03 + n21)*log(r) + (n10 + n14)*log(-((-1 + r)*r)) + 
      (n04 + n20)*log(r^2) + (n02 + n22)*log(2 + r - r^2) + n12*log(5 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLcc <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){
    optimize(logLcc,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcc <- function(r){
  #     L <- 2^(-2*x[,"n_00"] - x[,"n_02"] - 2*x[,"n_04"] - x[,"n_10"] - x[,"n_12"] - x[,"n_14"] - 2*x[,"n_20"] - x[,"n_22"] - 2*x[,"n_24"])*9^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"])*(1 - r)^(x[,"n_01"] + x[,"n_23"])*((-1 + r)^2)^(x[,"n_00"] + x[,"n_24"])*r^(x[,"n_03"] + x[,"n_21"])*(-((-1 + r)*r))^(x[,"n_10"] + x[,"n_14"])*(r^2)^(x[,"n_04"] + x[,"n_20"])*(2 + r - r^2)^(x[,"n_02"] + x[,"n_22"])*(5 - 2*r + 2*r^2)^x[,"n_12"]
  #     return(L)}
  
  ############################
  ## COUPLING REPULSION
  ############################
  logLcr <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - 2*n10 - 2*n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24)*log(9) + 
      (n00 + n04 + n20 + n24)*log(-((-1 + r)*r)) + n12*log(2 + r - r^2) + (n10 + n14)*log(1 - 2*r + 2*r^2) + 
      (n02 + n22)*log(5 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLcr <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){
    optimize(logLcr,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcr <- function(r){
  #     L <- 2^(-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"])*9^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"])*(-((-1 + r)*r))^(x[,"n_00"] + x[,"n_04"] + x[,"n_20"] + x[,"n_24"])*(2 + r - r^2)^x[,"n_12"]*(1 - 2*r + 2*r^2)^(x[,"n_10"] + x[,"n_14"])*(5 - 2*r + 2*r^2)^(x[,"n_02"] + x[,"n_22"])
  #     return(L)}
  
  ############################
  ## REPULSION COUPLING
  ############################
  logLrc <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - 2*n10 - 2*n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24)*log(9) + 
      (n00 + n04 + n20 + n24)*log(-((-1 + r)*r)) + n12*log(2 + r - r^2) + (n10 + n14)*log(1 - 2*r + 2*r^2) + 
      (n02 + n22)*log(5 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLrc <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){
    optimize(logLrc,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lrc <- function(r){
  #     L <- 2^(-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"])*9^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"])*(-((-1 + r)*r))^(x[,"n_00"] + x[,"n_04"] + x[,"n_20"] + x[,"n_24"])*(2 + r - r^2)^x[,"n_12"]*(1 - 2*r + 2*r^2)^(x[,"n_10"] + x[,"n_14"])*(5 - 2*r + 2*r^2)^(x[,"n_02"] + x[,"n_22"])
  #     return(L)}
  
  
  ############################
  ## REPULSION REPULSION
  ############################
  logLrr <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){ 
    L <- (-2*n00 - n02 - 2*n04 - n10 - n12 - n14 - 2*n20 - n22 - 2*n24)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24)*log(9) + 
      (n03 + n21)*log(1 - r) + (n04 + n20)*log((-1 + r)^2) + (n01 + n23)*log(r) + (n10 + n14)*log(-((-1 + r)*r)) + 
      (n00 + n24)*log(r^2) + (n02 + n22)*log(2 + r - r^2) + n12*log(5 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLrr <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24){
    optimize(logLrr,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lrr <- function(r){
  #     L <- 2^(-2*x[,"n_00"] - x[,"n_02"] - 2*x[,"n_04"] - x[,"n_10"] - x[,"n_12"] - x[,"n_14"] - 2*x[,"n_20"] - x[,"n_22"] - 2*x[,"n_24"])*9^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"])*(1 - r)^(x[,"n_03"] + x[,"n_21"])*((-1 + r)^2)^(x[,"n_04"] + x[,"n_20"])*r^(x[,"n_01"] + x[,"n_23"])*(-((-1 + r)*r))^(x[,"n_10"] + x[,"n_14"])*(r^2)^(x[,"n_00"] + x[,"n_24"])*(2 + r - r^2)^(x[,"n_02"] + x[,"n_22"])*(5 - 2*r + 2*r^2)^x[,"n_12"]
  #     return(L)}
  r_cc <- parallel::mcmapply(inter_logLcc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],
                             mc.cores = ncores)
  
  r_cr <-  parallel::mcmapply(inter_logLcr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],
                              mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],
                              mc.cores = ncores)
  
  r_rr <- parallel::mcmapply(inter_logLrr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],
                             mc.cores = ncores)
  
  LOD_cc <- log10(((1 - r_cc)^(x[,"n_01"] + x[,"n_23"])*((-1 + r_cc)^2)^(x[,"n_00"] + x[,"n_24"])*r_cc^(x[,"n_03"] + x[,"n_21"])*(-((-1 + r_cc)*r_cc))^(x[,"n_10"] + x[,"n_14"])*(r_cc^2)^(x[,"n_04"] + x[,"n_20"])*(2 + r_cc - r_cc^2)^(x[,"n_02"] + x[,"n_22"])*(5 - 2*r_cc + 2*r_cc^2)^x[,"n_12"])/((1 - 0.5)^(x[,"n_01"] + x[,"n_23"])*((-1 + 0.5)^2)^(x[,"n_00"] + x[,"n_24"])*0.5^(x[,"n_03"] + x[,"n_21"])*(-((-1 + 0.5)*0.5))^(x[,"n_10"] + x[,"n_14"])*(0.5^2)^(x[,"n_04"] + x[,"n_20"])*(2 + 0.5 - 0.5^2)^(x[,"n_02"] + x[,"n_22"])*(5 - 2*0.5 + 2*0.5^2)^x[,"n_12"]))
  LOD_cr <- log10(((-((-1 + r_cr)*r_cr))^(x[,"n_00"] + x[,"n_04"] + x[,"n_20"] + x[,"n_24"])*(2 + r_cr - r_cr^2)^x[,"n_12"]*(1 - 2*r_cr + 2*r_cr^2)^(x[,"n_10"] + x[,"n_14"])*(5 - 2*r_cr + 2*r_cr^2)^(x[,"n_02"] + x[,"n_22"]))/((-((-1 + 0.5)*0.5))^(x[,"n_00"] + x[,"n_04"] + x[,"n_20"] + x[,"n_24"])*(2 + 0.5 - 0.5^2)^x[,"n_12"]*(1 - 2*0.5 + 2*0.5^2)^(x[,"n_10"] + x[,"n_14"])*(5 - 2*0.5 + 2*0.5^2)^(x[,"n_02"] + x[,"n_22"])))
  LOD_rc <- log10(((-((-1 + r_rc)*r_rc))^(x[,"n_00"] + x[,"n_04"] + x[,"n_20"] + x[,"n_24"])*(2 + r_rc - r_rc^2)^x[,"n_12"]*(1 - 2*r_rc + 2*r_rc^2)^(x[,"n_10"] + x[,"n_14"])*(5 - 2*r_rc + 2*r_rc^2)^(x[,"n_02"] + x[,"n_22"]))/((-((-1 + 0.5)*0.5))^(x[,"n_00"] + x[,"n_04"] + x[,"n_20"] + x[,"n_24"])*(2 + 0.5 - 0.5^2)^x[,"n_12"]*(1 - 2*0.5 + 2*0.5^2)^(x[,"n_10"] + x[,"n_14"])*(5 - 2*0.5 + 2*0.5^2)^(x[,"n_02"] + x[,"n_22"])))
  LOD_rr <- log10(((1 - r_rr)^(x[,"n_03"] + x[,"n_21"])*((-1 + r_rr)^2)^(x[,"n_04"] + x[,"n_20"])*r_rr^(x[,"n_01"] + x[,"n_23"])*(-((-1 + r_rr)*r_rr))^(x[,"n_10"] + x[,"n_14"])*(r_rr^2)^(x[,"n_00"] + x[,"n_24"])*(2 + r_rr - r_rr^2)^(x[,"n_02"] + x[,"n_22"])*(5 - 2*r_rr + 2*r_rr^2)^x[,"n_12"])/((1 - 0.5)^(x[,"n_03"] + x[,"n_21"])*((-1 + 0.5)^2)^(x[,"n_04"] + x[,"n_20"])*0.5^(x[,"n_01"] + x[,"n_23"])*(-((-1 + 0.5)*0.5))^(x[,"n_10"] + x[,"n_14"])*(0.5^2)^(x[,"n_00"] + x[,"n_24"])*(2 + 0.5 - 0.5^2)^(x[,"n_02"] + x[,"n_22"])*(5 - 2*0.5 + 2*0.5^2)^x[,"n_12"]))
  
  ## Record the logL:
  logL_cc <- (-2*x[,"n_00"] - x[,"n_02"] - 2*x[,"n_04"] - x[,"n_10"] - x[,"n_12"] - x[,"n_14"] - 2*x[,"n_20"] - x[,"n_22"] - 2*x[,"n_24"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"])*log(9) + 
    (x[,"n_01"] + x[,"n_23"])*log(1 - r_cc) + (x[,"n_00"] + x[,"n_24"])*log((-1 + r_cc)^2) + (x[,"n_03"] + x[,"n_21"])*log(r_cc) + (x[,"n_10"] + x[,"n_14"])*log(-((-1 + r_cc)*r_cc)) + 
    (x[,"n_04"] + x[,"n_20"])*log(r_cc^2) + (x[,"n_02"] + x[,"n_22"])*log(2 + r_cc - r_cc^2) + x[,"n_12"]*log(5 - 2*r_cc + 2*r_cc^2)
  
  logL_cr <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"])*log(9) + 
    (x[,"n_00"] + x[,"n_04"] + x[,"n_20"] + x[,"n_24"])*log(-((-1 + r_cr)*r_cr)) + x[,"n_12"]*log(2 + r_cr - r_cr^2) + (x[,"n_10"] + x[,"n_14"])*log(1 - 2*r_cr + 2*r_cr^2) + 
    (x[,"n_02"] + x[,"n_22"])*log(5 - 2*r_cr + 2*r_cr^2)
  
  logL_rc <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"])*log(9) + 
    (x[,"n_00"] + x[,"n_04"] + x[,"n_20"] + x[,"n_24"])*log(-((-1 + r_rc)*r_rc)) + x[,"n_12"]*log(2 + r_rc - r_rc^2) + (x[,"n_10"] + x[,"n_14"])*log(1 - 2*r_rc + 2*r_rc^2) + 
    (x[,"n_02"] + x[,"n_22"])*log(5 - 2*r_rc + 2*r_rc^2)
  
  logL_rr <- (-2*x[,"n_00"] - x[,"n_02"] - 2*x[,"n_04"] - x[,"n_10"] - x[,"n_12"] - x[,"n_14"] - 2*x[,"n_20"] - x[,"n_22"] - 2*x[,"n_24"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"])*log(9) + 
    (x[,"n_03"] + x[,"n_21"])*log(1 - r_rr) + (x[,"n_04"] + x[,"n_20"])*log((-1 + r_rr)^2) + (x[,"n_01"] + x[,"n_23"])*log(r_rr) + (x[,"n_10"] + x[,"n_14"])*log(-((-1 + r_rr)*r_rr)) + 
    (x[,"n_00"] + x[,"n_24"])*log(r_rr^2) + (x[,"n_02"] + x[,"n_22"])*log(2 + r_rr - r_rr^2) + x[,"n_12"]*log(5 - 2*r_rr + 2*r_rr^2)
  
  return(list(
    r_mat=cbind(r_cc, r_cr, r_rc, r_rr),
    LOD_mat=cbind(LOD_cc, LOD_cr, LOD_rc, LOD_rr),
    logL_mat=cbind(logL_cc, logL_cr, logL_rc, logL_rr),
    phasing_strategy="MLL", 
    possible_phases=c("coupling coupling",
                      "coupling repulsion",
                      "repulsion coupling",
                      "repulsion repulsion")
  ))
  
}

#' @rdname r4_functions
#' @noRd
r4_1.2_2.1<-function(x, ncores=1){
  ############################
  ## COUPLING COUPLING
  ############################
  logLcc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(36) + 
      (n00 + n33)*log((-1 + r)^2) + (n02 + n13 + n20 + n31)*log(-((-3 + r)*r)) + (n03 + n30)*log(r^2) + 
      (n01 + n10 + n23 + n32)*log(2 - r - r^2) + (n11 + n22)*log(8 - 4*r + r^2) + (n12 + n21)*log(5 + 2*r + r^2)
    return(L)}
  
  inter_logLcc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## COUPLING REPULSION
  ############################
  logLcr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(36) + 
      (n00 + n03 + n30 + n33)*log(-((-1 + r)*r)) + (n11 + n12 + n21 + n22)*log(6 + r - r^2) + 
      (n01 + n13 + n20 + n32)*log(1 + r^2) + (n02 + n10 + n23 + n31)*log(2 - 2*r + r^2)
    return(L)}
  
  inter_logLcr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION COUPLING
  ############################
  logLrc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(36) + 
      (n00 + n03 + n30 + n33)*log(-((-1 + r)*r)) + (n11 + n12 + n21 + n22)*log(6 + r - r^2) + 
      (n02 + n10 + n23 + n31)*log(1 + r^2) + (n01 + n13 + n20 + n32)*log(2 - 2*r + r^2)
    return(L)}
  
  inter_logLrc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION REPULSION
  ############################
  logLrr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(36) + 
      (n03 + n30)*log((-1 + r)^2) + (n01 + n10 + n23 + n32)*log(-((-3 + r)*r)) + (n00 + n33)*log(r^2) + 
      (n02 + n13 + n20 + n31)*log(2 - r - r^2) + (n12 + n21)*log(8 - 4*r + r^2) + (n11 + n22)*log(5 + 2*r + r^2)
    return(L)}
  
  inter_logLrr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_cc <- parallel::mcmapply(inter_logLcc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_cr <- parallel::mcmapply(inter_logLcr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  r_rr <-  parallel::mcmapply(inter_logLrr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  LOD_cc <- (x[,"n_00"] + x[,"n_33"])*log10((-1 + r_cc)^2) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(-((-3 + r_cc)*r_cc)) + (x[,"n_03"] + x[,"n_30"])*log10(r_cc^2) + 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 - r_cc - r_cc^2) + (x[,"n_11"] + x[,"n_22"])*log10(8 - 4*r_cc + r_cc^2) + (x[,"n_12"] + x[,"n_21"])*log10(5 + 2*r_cc + r_cc^2) - 
    (x[,"n_00"] + x[,"n_33"])*log10((-1 + 0.5)^2) - (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(-((-3 + 0.5)*0.5)) - (x[,"n_03"] + x[,"n_30"])*log10(0.5^2) - 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 - 0.5 - 0.5^2) - (x[,"n_11"] + x[,"n_22"])*log10(8 - 4*0.5 + 0.5^2) - (x[,"n_12"] + x[,"n_21"])*log10(5 + 2*0.5 + 0.5^2)
  
  LOD_cr <- (x[,"n_00"] + x[,"n_03"] + x[,"n_30"] + x[,"n_33"])*log10(-((-1 + r_cr)*r_cr)) + (x[,"n_11"] + x[,"n_12"] + x[,"n_21"] + x[,"n_22"])*log10(6 + r_cr - r_cr^2) + (x[,"n_01"] + x[,"n_13"] + x[,"n_20"] + x[,"n_32"])*log10(1 + r_cr^2) + 
    (x[,"n_02"] + x[,"n_10"] + x[,"n_23"] + x[,"n_31"])*log10(2 - 2*r_cr + r_cr^2) - (x[,"n_00"] + x[,"n_03"] + x[,"n_30"] + x[,"n_33"])*log10(-((-1 + 0.5)*0.5)) - (x[,"n_11"] + x[,"n_12"] + x[,"n_21"] + x[,"n_22"])*log10(6 + 0.5 - 0.5^2) - 
    (x[,"n_01"] + x[,"n_13"] + x[,"n_20"] + x[,"n_32"])*log10(1 + 0.5^2) - (x[,"n_02"] + x[,"n_10"] + x[,"n_23"] + x[,"n_31"])*log10(2 - 2*0.5 + 0.5^2)
  
  LOD_rc <- (x[,"n_00"] + x[,"n_03"] + x[,"n_30"] + x[,"n_33"])*log10(-((-1 + r_rc)*r_rc)) + (x[,"n_11"] + x[,"n_12"] + x[,"n_21"] + x[,"n_22"])*log10(6 + r_rc - r_rc^2) + (x[,"n_02"] + x[,"n_10"] + x[,"n_23"] + x[,"n_31"])*log10(1 + r_rc^2) + 
    (x[,"n_01"] + x[,"n_13"] + x[,"n_20"] + x[,"n_32"])*log10(2 - 2*r_rc + r_rc^2) - (x[,"n_00"] + x[,"n_03"] + x[,"n_30"] + x[,"n_33"])*log10(-((-1 + 0.5)*0.5)) - (x[,"n_11"] + x[,"n_12"] + x[,"n_21"] + x[,"n_22"])*log10(6 + 0.5 - 0.5^2) -
    (x[,"n_02"] + x[,"n_10"] + x[,"n_23"] + x[,"n_31"])*log10(1 + 0.5^2) - (x[,"n_01"] + x[,"n_13"] + x[,"n_20"] + x[,"n_32"])*log10(2 - 2*0.5 + 0.5^2)
  
  LOD_rr <- (x[,"n_03"] + x[,"n_30"])*log10((-1 + r_rr)^2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(-((-3 + r_rr)*r_rr)) + (x[,"n_00"] + x[,"n_33"])*log10(r_rr^2) + 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(2 - r_rr - r_rr^2) + (x[,"n_12"] + x[,"n_21"])*log10(8 - 4*r_rr + r_rr^2) + (x[,"n_11"] + x[,"n_22"])*log10(5 + 2*r_rr + r_rr^2) - 
    (x[,"n_03"] + x[,"n_30"])*log10((-1 + 0.5)^2) - (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(-((-3 + 0.5)*0.5)) - (x[,"n_00"] + x[,"n_33"])*log10(0.5^2) - 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(2 - 0.5 - 0.5^2) - (x[,"n_12"] + x[,"n_21"])*log10(8 - 4*0.5 + 0.5^2) - (x[,"n_11"] + x[,"n_22"])*log10(5 + 2*0.5 + 0.5^2)
  
  ## Record the logL values:
  logL_cc <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(36) + 
    (x[,"n_00"] + x[,"n_33"])*log((-1 + r_cc)^2) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(-((-3 + r_cc)*r_cc)) + (x[,"n_03"] + x[,"n_30"])*log(r_cc^2) + 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(2 - r_cc - r_cc^2) + (x[,"n_11"] + x[,"n_22"])*log(8 - 4*r_cc + r_cc^2) + (x[,"n_12"] + x[,"n_21"])*log(5 + 2*r_cc + r_cc^2)
  
  logL_cr <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(36) + 
    (x[,"n_00"] + x[,"n_03"] + x[,"n_30"] + x[,"n_33"])*log(-((-1 + r_cr)*r_cr)) + (x[,"n_11"] + x[,"n_12"] + x[,"n_21"] + x[,"n_22"])*log(6 + r_cr - r_cr^2) + 
    (x[,"n_01"] + x[,"n_13"] + x[,"n_20"] + x[,"n_32"])*log(1 + r_cr^2) + (x[,"n_02"] + x[,"n_10"] + x[,"n_23"] + x[,"n_31"])*log(2 - 2*r_cr + r_cr^2)
  
  logL_rc <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(36) + 
    (x[,"n_00"] + x[,"n_03"] + x[,"n_30"] + x[,"n_33"])*log(-((-1 + r_rc)*r_rc)) + (x[,"n_11"] + x[,"n_12"] + x[,"n_21"] + x[,"n_22"])*log(6 + r_rc - r_rc^2) + 
    (x[,"n_02"] + x[,"n_10"] + x[,"n_23"] + x[,"n_31"])*log(1 + r_rc^2) + (x[,"n_01"] + x[,"n_13"] + x[,"n_20"] + x[,"n_32"])*log(2 - 2*r_rc + r_rc^2)
  
  logL_rr <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(36) + 
    (x[,"n_03"] + x[,"n_30"])*log((-1 + r_rr)^2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(-((-3 + r_rr)*r_rr)) + (x[,"n_00"] + x[,"n_33"])*log(r_rr^2) + 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(2 - r_rr - r_rr^2) + (x[,"n_12"] + x[,"n_21"])*log(8 - 4*r_rr + r_rr^2) + (x[,"n_11"] + x[,"n_22"])*log(5 + 2*r_rr + r_rr^2)
  
  return(list(
    r_mat=cbind(r_cc, r_cr, r_rc, r_rr),
    LOD_mat=cbind(LOD_cc, LOD_cr, LOD_rc, LOD_rr),
    logL_mat=cbind(logL_cc, logL_cr, logL_rc, logL_rr),
    phasing_strategy="MLL", 
    possible_phases=c("coupling coupling",
                      "coupling repulsion",
                      "repulsion coupling",
                      "repulsion repulsion")
  ))
  
}

#' @rdname r4_functions
#' @noRd
r4_1.3_2.1<-function(x, ncores=1){ # simplex triplex does not give nulliplex progeny?
  ############################
  ## COUPLING COUPLING
  ############################
  logLcc <- function(r,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-2*n10 - 2*n11 - 2*n12 - 2*n13 - 2*n20 - n21 - n22 - 2*n23 - 2*n30 - 2*n31 - 2*n32 - 2*n33)*log(2) + 
      (-n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(9) + 
      (n10 + n33)*log((-2 + r)*(-1 + r)) + (n13 + n30)*log(r*(1 + r)) + (n20 + n23)*log(1 + 2*r - 2*r^2) + 
      (n11 + n32)*log(5 - 2*r - r^2) + (n12 + n31)*log(2 + 4*r - r^2) + (n21 + n22)*log(4 - r + r^2)
    return(L)}
  
  inter_logLcc <- function(n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcc,c(0,0.5),n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## COUPLING REPULSION
  ############################
  logLcr <- function(r,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-2*n10 - 2*n11 - 2*n12 - 2*n13 - 2*n20 - n21 - n22 - 2*n23 - 2*n30 - 2*n31 - 2*n32 - 2*n33)*log(2) + 
      (-n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(3) + 
      (n10 + n13 + n30 + n33)*log(-((-1 + r)*r)) + (n21 + n22)*log(1 + r - r^2) + (n11 + n32)*log(1 + r^2) + 
      (n12 + n31)*log(2 - 2*r + r^2) + (n20 + n23)*log(1 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLcr <- function(n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcr,c(0,0.5),n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION COUPLING
  ############################
  logLrc <- function(r,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-2*n10 - 2*n11 - 2*n12 - 2*n13 - n20 - 2*n21 - 2*n22 - n23 - 2*n30 - 2*n31 - 2*n32 - 2*n33)*log(2) + 
      (-n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(9) + 
      (n10 + n33)*log(-((-2 + r)*r)) + (n13 + n30)*log(-((-1 + r)*(1 + r))) + (n21 + n22)*log(7 + 2*r - 2*r^2) + 
      (n20 + n23)*log(1 - r + r^2) + (n11 + n12 + n31 + n32)*log(4 - r + r^2)
    return(L)}
  
  inter_logLrc <- function(n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrc,c(0,0.5),n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION REPULSION
  ############################
  logLrr <- function(r,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-2*n10 - 2*n11 - 2*n12 - 2*n13 - n20 - 2*n21 - 2*n22 - n23 - 2*n30 - 2*n31 - 2*n32 - 2*n33)*log(2) + 
      (-n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(3) + 
      (n13 + n30)*log((-1 + r)^2) + (n11 + n32)*log(-((-3 + r)*r)) + (n20 + n23)*log(-((-1 + r)*r)) +
      (n10 + n33)*log(r^2) + (n12 + n31)*log(2 - r - r^2) + (n21 + n22)*log(3 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLrr <- function(n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrr,c(0,0.5),n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_cc <- parallel::mcmapply(inter_logLcc,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_cr <- parallel::mcmapply(inter_logLcr,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  r_rr <-  parallel::mcmapply(inter_logLrr,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  LOD_cc <- (x[,"n_10"] + x[,"n_33"])*log10((-2 + r_cc)*(-1 + r_cc)) + (x[,"n_13"] + x[,"n_30"])*log10(r_cc*(1 + r_cc)) + (x[,"n_20"] + x[,"n_23"])*log10(1 + 2*r_cc - 2*r_cc^2) + 
    (x[,"n_11"] + x[,"n_32"])*log10(5 - 2*r_cc - r_cc^2) + (x[,"n_12"] + x[,"n_31"])*log10(2 + 4*r_cc - r_cc^2) + (x[,"n_21"] + x[,"n_22"])*log10(4 - r_cc + r_cc^2) - 
    (x[,"n_10"] + x[,"n_33"])*log10((-2 + 0.5)*(-1 + 0.5)) - (x[,"n_13"] + x[,"n_30"])*log10(0.5*(1 + 0.5)) - (x[,"n_20"] + x[,"n_23"])*log10(1 + 2*0.5 - 2*0.5^2) - 
    (x[,"n_11"] + x[,"n_32"])*log10(5 - 2*0.5 - 0.5^2) - (x[,"n_12"] + x[,"n_31"])*log10(2 + 4*0.5 - 0.5^2) - (x[,"n_21"] + x[,"n_22"])*log10(4 - 0.5 + 0.5^2)
  
  LOD_cr <- (x[,"n_10"] + x[,"n_13"] + x[,"n_30"] + x[,"n_33"])*log10(-((-1 + r_cr)*r_cr)) + (x[,"n_21"] + x[,"n_22"])*log10(1 + r_cr - r_cr^2) + (x[,"n_11"] + x[,"n_32"])*log10(1 + r_cr^2) + 
    (x[,"n_12"] + x[,"n_31"])*log10(2 - 2*r_cr + r_cr^2) + (x[,"n_20"] + x[,"n_23"])*log10(1 - 2*r_cr + 2*r_cr^2) - (x[,"n_10"] + x[,"n_13"] + x[,"n_30"] + x[,"n_33"])*log10(-((-1 + 0.5)*0.5)) - 
    (x[,"n_21"] + x[,"n_22"])*log10(1 + 0.5 - 0.5^2) - (x[,"n_11"] + x[,"n_32"])*log10(1 + 0.5^2) - (x[,"n_12"] + x[,"n_31"])*log10(2 - 2*0.5 + 0.5^2) - (x[,"n_20"] + x[,"n_23"])*log10(1 - 2*0.5 + 2*0.5^2)
  
  LOD_rc <- (x[,"n_10"] + x[,"n_33"])*log10(-((-2 + r_rc)*r_rc)) + (x[,"n_13"] + x[,"n_30"])*log10(-((-1 + r_rc)*(1 + r_rc))) + (x[,"n_21"] + x[,"n_22"])*log10(7 + 2*r_rc - 2*r_rc^2) + 
    (x[,"n_20"] + x[,"n_23"])*log10(1 - r_rc + r_rc^2) + (x[,"n_11"] + x[,"n_12"] + x[,"n_31"] + x[,"n_32"])*log10(4 - r_rc + r_rc^2) - (x[,"n_10"] + x[,"n_33"])*log10(-((-2 + 0.5)*0.5)) - 
    (x[,"n_13"] + x[,"n_30"])*log10(-((-1 + 0.5)*(1 + 0.5))) - (x[,"n_21"] + x[,"n_22"])*log10(7 + 2*0.5 - 2*0.5^2) - (x[,"n_20"] + x[,"n_23"])*log10(1 - 0.5 + 0.5^2) - 
    (x[,"n_11"] + x[,"n_12"] + x[,"n_31"] + x[,"n_32"])*log10(4 - 0.5 + 0.5^2)
  
  LOD_rr <- (x[,"n_13"] + x[,"n_30"])*log10((-1 + r_rr)^2) + (x[,"n_11"] + x[,"n_32"])*log10(-((-3 + r_rr)*r_rr)) + (x[,"n_20"] + x[,"n_23"])*log10(-((-1 + r_rr)*r_rr)) + 
    (x[,"n_10"] + x[,"n_33"])*log10(r_rr^2) + (x[,"n_12"] + x[,"n_31"])*log10(2 - r_rr - r_rr^2) + (x[,"n_21"] + x[,"n_22"])*log10(3 - 2*r_rr + 2*r_rr^2) - (x[,"n_13"] + x[,"n_30"])*log10((-1 + 0.5)^2) -
    (x[,"n_11"] + x[,"n_32"])*log10(-((-3 + 0.5)*0.5)) - (x[,"n_20"] + x[,"n_23"])*log10(-((-1 + 0.5)*0.5)) - (x[,"n_10"] + x[,"n_33"])*log10(0.5^2) - 
    (x[,"n_12"] + x[,"n_31"])*log10(2 - 0.5 - 0.5^2) - (x[,"n_21"] + x[,"n_22"])*log10(3 - 2*0.5 + 2*0.5^2)
  
  ## Record the logL:
  logL_cc <- (-2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - 2*x[,"n_23"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(9) + 
    (x[,"n_10"] + x[,"n_33"])*log((-2 + r_cc)*(-1 + r_cc)) + (x[,"n_13"] + x[,"n_30"])*log(r_cc*(1 + r_cc)) + (x[,"n_20"] + x[,"n_23"])*log(1 + 2*r_cc - 2*r_cc^2) + 
    (x[,"n_11"] + x[,"n_32"])*log(5 - 2*r_cc - r_cc^2) + (x[,"n_12"] + x[,"n_31"])*log(2 + 4*r_cc - r_cc^2) + (x[,"n_21"] + x[,"n_22"])*log(4 - r_cc + r_cc^2)
  
  logL_cr <- (-2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - 2*x[,"n_23"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(3) + 
    (x[,"n_10"] + x[,"n_13"] + x[,"n_30"] + x[,"n_33"])*log(-((-1 + r_cr)*r_cr)) + (x[,"n_21"] + x[,"n_22"])*log(1 + r_cr - r_cr^2) + (x[,"n_11"] + x[,"n_32"])*log(1 + r_cr^2) + 
    (x[,"n_12"] + x[,"n_31"])*log(2 - 2*r_cr + r_cr^2) + (x[,"n_20"] + x[,"n_23"])*log(1 - 2*r_cr + 2*r_cr^2)
  
  logL_rc <- (-2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(9) + 
    (x[,"n_10"] + x[,"n_33"])*log(-((-2 + r_rc)*r_rc)) + (x[,"n_13"] + x[,"n_30"])*log(-((-1 + r_rc)*(1 + r_rc))) + (x[,"n_21"] + x[,"n_22"])*log(7 + 2*r_rc - 2*r_rc^2) + 
    (x[,"n_20"] + x[,"n_23"])*log(1 - r_rc + r_rc^2) + (x[,"n_11"] + x[,"n_12"] + x[,"n_31"] + x[,"n_32"])*log(4 - r_rc + r_rc^2)
  
  logL_rr <- (-2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(3) + 
    (x[,"n_13"] + x[,"n_30"])*log((-1 + r_rr)^2) + (x[,"n_11"] + x[,"n_32"])*log(-((-3 + r_rr)*r_rr)) + (x[,"n_20"] + x[,"n_23"])*log(-((-1 + r_rr)*r_rr)) +
    (x[,"n_10"] + x[,"n_33"])*log(r_rr^2) + (x[,"n_12"] + x[,"n_31"])*log(2 - r_rr - r_rr^2) + (x[,"n_21"] + x[,"n_22"])*log(3 - 2*r_rr + 2*r_rr^2)
  
  return(list(
    r_mat=cbind(r_cc, r_cr, r_rc, r_rr),
    LOD_mat=cbind(LOD_cc, LOD_cr, LOD_rc, LOD_rr),
    logL_mat=cbind(logL_cc, logL_cr, logL_rc, logL_rr),
    phasing_strategy="MLL", 
    possible_phases=c("coupling coupling",
                      "coupling repulsion",
                      "repulsion coupling",
                      "repulsion repulsion")
  ))
}

#' @rdname r4_functions
#' @noRd
r4_2.1_2.2<-function(x, ncores=1){
  ############################
  ## COUPLING COUPLING
  ############################
  logLcc <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34)*log(2) + 
      (-2*n00 - 2*n01 - 2*n02 - 2*n03 - 2*n04 - n10 - 2*n11 - 2*n12 - 2*n13 - n14 - n20 - 2*n21 - 2*n22 - 2*n23 - n24 - 2*n30 - 2*n31 - 2*n32 - 2*n33 - 2*n34)*log(3) + 
      (n02 + n32)*log(5) + (n00 + n34)*log(-(-1 + r)^3) + (n02 + n32)*log(-((-1 + r)*r)) + (n10 + n24)*log((-1 + r)^2*r) + (n03 + n31)*log(-((-2 + r)*r^2)) +
      (n14 + n20)*log(-((-1 + r)*r^2)) + (n04 + n30)*log(r^3) + (n01 + n33)*log((-1 + r)^2*(1 + r)) + (n13 + n21)*log(r*(5 - 5*r + 3*r^2)) + 
      (n12 + n22)*log(9 - 5*r + 5*r^2) + (n11 + n23)*log(3 - 4*r + 4*r^2 - 3*r^3)
    return(L)}
  
  inter_logLcc <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLcc,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## COUPLING REPULSION
  ############################
  logLcr <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n00 + n34)*log((-1 + r)^2*r) + (n04 + n30)*log(-((-1 + r)*r^2)) + (n02 + n32)*log(1 + r - r^2) + (n03 + n31)*log(r*(1 - r + r^2)) + 
      (n12 + n22)*log(8 - r + r^2) + (n14 + n20)*log(r*(2 - 4*r + 3*r^2)) + (n13 + n21)*log(2 - 2*r + 4*r^2 - 3*r^3) + 
      (n10 + n24)*log(1 - 3*r + 5*r^2 - 3*r^3) + (n01 + n33)*log(1 - 2*r + 2*r^2 - r^3) + (n11 + n23)*log(1 + 3*r - 5*r^2 + 3*r^3)
    return(L)}
  
  inter_logLcr <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLcr,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## MIXED COUPLING
  ############################
  logLmc <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-3*n00 - 2*n01 - 3*n02 - 2*n03 - 3*n04 - 3*n10 - 2*n11 - 3*n12 - 2*n13 - 3*n14 - 3*n20 - 2*n21 - 3*n22 - 2*n23 - 3*n24 - 3*n30 - 2*n31 - 3*n32 - 2*n33 - 3*n34)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n00 + n34)*log((-1 + r)^2*r) + (n04 + n30)*log(-((-1 + r)*r^2)) + (n12 + n22)*log(14 + 3*r - 3*r^2) + (n03 + n31)*log(r*(2 - 2*r + r^2)) + 
      (n14 + n20)*log(r*(3 - 4*r + 3*r^2)) + (n02 + n32)*log(4 - 3*r + 3*r^2) + (n10 + n24)*log(2 - 4*r + 5*r^2 - 3*r^3) + 
      (n13 + n21)*log(3 - r + 5*r^2 - 3*r^3) + (n01 + n33)*log(1 - r + r^2 - r^3) + (n11 + n23)*log(4 - 4*r^2 + 3*r^3)
    return(L)}
  
  inter_logLmc <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLmc,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## MIXED REPULSION 
  ############################
  logLmr <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-3*n00 - 2*n01 - 3*n02 - 2*n03 - 3*n04 - 3*n10 - 2*n11 - 3*n12 - 2*n13 - 3*n14 - 3*n20 - 2*n21 - 3*n22 - 2*n23 - 3*n24 - 3*n30 - 2*n31 - 3*n32 - 2*n33 - 3*n34)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n04 + n30)*log((-1 + r)^2*r) + (n00 + n34)*log(-((-1 + r)*r^2)) + (n12 + n22)*log(14 + 3*r - 3*r^2) + (n01 + n33)*log(r*(2 - 2*r + r^2)) + 
      (n10 + n24)*log(r*(3 - 4*r + 3*r^2)) + (n02 + n32)*log(4 - 3*r + 3*r^2) + (n14 + n20)*log(2 - 4*r + 5*r^2 - 3*r^3) + 
      (n11 + n23)*log(3 - r + 5*r^2 - 3*r^3) + (n03 + n31)*log(1 - r + r^2 - r^3) + (n13 + n21)*log(4 - 4*r^2 + 3*r^3)
    return(L)}
  
  inter_logLmr <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLmr,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION COUPLING
  ############################
  logLrc <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n04 + n30)*log((-1 + r)^2*r) + (n00 + n34)*log(-((-1 + r)*r^2)) + (n02 + n32)*log(1 + r - r^2) + (n01 + n33)*log(r*(1 - r + r^2)) + 
      (n12 + n22)*log(8 - r + r^2) + (n10 + n24)*log(r*(2 - 4*r + 3*r^2)) + (n11 + n23)*log(2 - 2*r + 4*r^2 - 3*r^3) + 
      (n14 + n20)*log(1 - 3*r + 5*r^2 - 3*r^3) + (n03 + n31)*log(1 - 2*r + 2*r^2 - r^3) + (n13 + n21)*log(1 + 3*r - 5*r^2 + 3*r^3)
    return(L)}
  
  inter_logLrc <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLrc,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION REPULSION
  ############################
  logLrr <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34)*log(2) + 
      (-2*n00 - 2*n01 - 2*n02 - 2*n03 - 2*n04 - n10 - 2*n11 - 2*n12 - 2*n13 - n14 - n20 - 2*n21 - 2*n22 - 2*n23 - n24 - 2*n30 - 2*n31 - 2*n32 - 2*n33 - 2*n34)*log(3) + 
      (n02 + n32)*log(5) + (n04 + n30)*log(-(-1 + r)^3) + (n02 + n32)*log(-((-1 + r)*r)) + (n14 + n20)*log((-1 + r)^2*r) + 
      (n01 + n33)*log(-((-2 + r)*r^2)) + (n10 + n24)*log(-((-1 + r)*r^2)) + (n00 + n34)*log(r^3) + (n03 + n31)*log((-1 + r)^2*(1 + r)) + 
      (n11 + n23)*log(r*(5 - 5*r + 3*r^2)) + (n12 + n22)*log(9 - 5*r + 5*r^2) + (n13 + n21)*log(3 - 4*r + 4*r^2 - 3*r^3)
    return(L)}
  
  inter_logLrr <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLrr,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_cc <- parallel::mcmapply(inter_logLcc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                             x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                             mc.cores = ncores)
  
  r_cr <- parallel::mcmapply(inter_logLcr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                             x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                             mc.cores = ncores)
  
  r_mc <- parallel::mcmapply(inter_logLmc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                             x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                             mc.cores = ncores)
  
  r_mr <-  parallel::mcmapply(inter_logLmr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                              x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                              mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                              x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                              mc.cores = ncores)
  
  r_rr <-  parallel::mcmapply(inter_logLrr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                              x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                              mc.cores = ncores)
  
  LOD_cc <- (x[,"n_00"] + x[,"n_34"])*log10(-(-1 + r_cc)^3) + (x[,"n_02"] + x[,"n_32"])*log10(-((-1 + r_cc)*r_cc)) + (x[,"n_10"] + x[,"n_24"])*log10((-1 + r_cc)^2*r_cc) + 
    (x[,"n_03"] + x[,"n_31"])*log10(-((-2 + r_cc)*r_cc^2)) + (x[,"n_14"] + x[,"n_20"])*log10(-((-1 + r_cc)*r_cc^2)) + (x[,"n_04"] + x[,"n_30"])*log10(r_cc^3) + 
    (x[,"n_01"] + x[,"n_33"])*log10((-1 + r_cc)^2*(1 + r_cc)) + (x[,"n_13"] + x[,"n_21"])*log10(r_cc*(5 - 5*r_cc + 3*r_cc^2)) + (x[,"n_12"] + x[,"n_22"])*log10(9 - 5*r_cc + 5*r_cc^2) + 
    (x[,"n_11"] + x[,"n_23"])*log10(3 - 4*r_cc + 4*r_cc^2 - 3*r_cc^3)  -  ( (x[,"n_00"] + x[,"n_34"])*log10(-(-1 + 0.5)^3) + (x[,"n_02"] + x[,"n_32"])*log10(-((-1 + 0.5)*0.5)) + (x[,"n_10"] + x[,"n_24"])*log10((-1 + 0.5)^2*0.5) + 
                                                                              (x[,"n_03"] + x[,"n_31"])*log10(-((-2 + 0.5)*0.5^2)) + (x[,"n_14"] + x[,"n_20"])*log10(-((-1 + 0.5)*0.5^2)) + (x[,"n_04"] + x[,"n_30"])*log10(0.5^3) + 
                                                                              (x[,"n_01"] + x[,"n_33"])*log10((-1 + 0.5)^2*(1 + 0.5)) + (x[,"n_13"] + x[,"n_21"])*log10(0.5*(5 - 5*0.5 + 3*0.5^2)) + (x[,"n_12"] + x[,"n_22"])*log10(9 - 5*0.5 + 5*0.5^2) + 
                                                                              (x[,"n_11"] + x[,"n_23"])*log10(3 - 4*0.5 + 4*0.5^2 - 3*0.5^3)  )
  
  LOD_cr <- (x[,"n_00"] + x[,"n_34"])*log10((-1 + r_cr)^2*r_cr) + (x[,"n_04"] + x[,"n_30"])*log10(-((-1 + r_cr)*r_cr^2)) + (x[,"n_02"] + x[,"n_32"])*log10(1 + r_cr - r_cr^2) + 
    (x[,"n_03"] + x[,"n_31"])*log10(r_cr*(1 - r_cr + r_cr^2)) + (x[,"n_12"] + x[,"n_22"])*log10(8 - r_cr + r_cr^2) + (x[,"n_14"] + x[,"n_20"])*log10(r_cr*(2 - 4*r_cr + 3*r_cr^2)) + 
    (x[,"n_13"] + x[,"n_21"])*log10(2 - 2*r_cr + 4*r_cr^2 - 3*r_cr^3) + (x[,"n_10"] + x[,"n_24"])*log10(1 - 3*r_cr + 5*r_cr^2 - 3*r_cr^3) + (x[,"n_01"] + x[,"n_33"])*log10(1 - 2*r_cr + 2*r_cr^2 - r_cr^3) + 
    (x[,"n_11"] + x[,"n_23"])*log10(1 + 3*r_cr - 5*r_cr^2 + 3*r_cr^3)  - ( (x[,"n_00"] + x[,"n_34"])*log10((-1 + 0.5)^2*0.5) + (x[,"n_04"] + x[,"n_30"])*log10(-((-1 + 0.5)*0.5^2)) + (x[,"n_02"] + x[,"n_32"])*log10(1 + 0.5 - 0.5^2) + 
                                                                             (x[,"n_03"] + x[,"n_31"])*log10(0.5*(1 - 0.5 + 0.5^2)) + (x[,"n_12"] + x[,"n_22"])*log10(8 - 0.5 + 0.5^2) + (x[,"n_14"] + x[,"n_20"])*log10(0.5*(2 - 4*0.5 + 3*0.5^2)) + 
                                                                             (x[,"n_13"] + x[,"n_21"])*log10(2 - 2*0.5 + 4*0.5^2 - 3*0.5^3) + (x[,"n_10"] + x[,"n_24"])*log10(1 - 3*0.5 + 5*0.5^2 - 3*0.5^3) + (x[,"n_01"] + x[,"n_33"])*log10(1 - 2*0.5 + 2*0.5^2 - 0.5^3) + 
                                                                             (x[,"n_11"] + x[,"n_23"])*log10(1 + 3*0.5 - 5*0.5^2 + 3*0.5^3) )
  
  LOD_mc <- (x[,"n_00"] + x[,"n_34"])*log10((-1 + r_mc)^2*r_mc) + (x[,"n_04"] + x[,"n_30"])*log10(-((-1 + r_mc)*r_mc^2)) + (x[,"n_12"] + x[,"n_22"])*log10(14 + 3*r_mc - 3*r_mc^2) + 
    (x[,"n_03"] + x[,"n_31"])*log10(r_mc*(2 - 2*r_mc + r_mc^2)) + (x[,"n_14"] + x[,"n_20"])*log10(r_mc*(3 - 4*r_mc + 3*r_mc^2)) + (x[,"n_02"] + x[,"n_32"])*log10(4 - 3*r_mc + 3*r_mc^2) + 
    (x[,"n_10"] + x[,"n_24"])*log10(2 - 4*r_mc + 5*r_mc^2 - 3*r_mc^3) + (x[,"n_13"] + x[,"n_21"])*log10(3 - r_mc + 5*r_mc^2 - 3*r_mc^3) + (x[,"n_01"] + x[,"n_33"])*log10(1 - r_mc + r_mc^2 - r_mc^3) + 
    (x[,"n_11"] + x[,"n_23"])*log10(4 - 4*r_mc^2 + 3*r_mc^3)  -  ( (x[,"n_00"] + x[,"n_34"])*log10((-1 + 0.5)^2*0.5) + (x[,"n_04"] + x[,"n_30"])*log10(-((-1 + 0.5)*0.5^2)) + (x[,"n_12"] + x[,"n_22"])*log10(14 + 3*0.5 - 3*0.5^2) + 
                                                                     (x[,"n_03"] + x[,"n_31"])*log10(0.5*(2 - 2*0.5 + 0.5^2)) + (x[,"n_14"] + x[,"n_20"])*log10(0.5*(3 - 4*0.5 + 3*0.5^2)) + (x[,"n_02"] + x[,"n_32"])*log10(4 - 3*0.5 + 3*0.5^2) + 
                                                                     (x[,"n_10"] + x[,"n_24"])*log10(2 - 4*0.5 + 5*0.5^2 - 3*0.5^3) + (x[,"n_13"] + x[,"n_21"])*log10(3 - 0.5 + 5*0.5^2 - 3*0.5^3) + (x[,"n_01"] + x[,"n_33"])*log10(1 - 0.5 + 0.5^2 - 0.5^3) + 
                                                                     (x[,"n_11"] + x[,"n_23"])*log10(4 - 4*0.5^2 + 3*0.5^3) )
  
  LOD_mr <- (x[,"n_04"] + x[,"n_30"])*log10((-1 + r_mr)^2*r_mr) + (x[,"n_00"] + x[,"n_34"])*log10(-((-1 + r_mr)*r_mr^2)) + (x[,"n_12"] + x[,"n_22"])*log10(14 + 3*r_mr - 3*r_mr^2) + 
    (x[,"n_01"] + x[,"n_33"])*log10(r_mr*(2 - 2*r_mr + r_mr^2)) + (x[,"n_10"] + x[,"n_24"])*log10(r_mr*(3 - 4*r_mr + 3*r_mr^2)) + (x[,"n_02"] + x[,"n_32"])*log10(4 - 3*r_mr + 3*r_mr^2) + 
    (x[,"n_14"] + x[,"n_20"])*log10(2 - 4*r_mr + 5*r_mr^2 - 3*r_mr^3) + (x[,"n_11"] + x[,"n_23"])*log10(3 - r_mr + 5*r_mr^2 - 3*r_mr^3) + (x[,"n_03"] + x[,"n_31"])*log10(1 - r_mr + r_mr^2 - r_mr^3) + 
    (x[,"n_13"] + x[,"n_21"])*log10(4 - 4*r_mr^2 + 3*r_mr^3)  - ( (x[,"n_04"] + x[,"n_30"])*log10((-1 + 0.5)^2*0.5) + (x[,"n_00"] + x[,"n_34"])*log10(-((-1 + 0.5)*0.5^2)) + (x[,"n_12"] + x[,"n_22"])*log10(14 + 3*0.5 - 3*0.5^2) + 
                                                                    (x[,"n_01"] + x[,"n_33"])*log10(0.5*(2 - 2*0.5 + 0.5^2)) + (x[,"n_10"] + x[,"n_24"])*log10(0.5*(3 - 4*0.5 + 3*0.5^2)) + (x[,"n_02"] + x[,"n_32"])*log10(4 - 3*0.5 + 3*0.5^2) + 
                                                                    (x[,"n_14"] + x[,"n_20"])*log10(2 - 4*0.5 + 5*0.5^2 - 3*0.5^3) + (x[,"n_11"] + x[,"n_23"])*log10(3 - 0.5 + 5*0.5^2 - 3*0.5^3) + (x[,"n_03"] + x[,"n_31"])*log10(1 - 0.5 + 0.5^2 - 0.5^3) + 
                                                                    (x[,"n_13"] + x[,"n_21"])*log10(4 - 4*0.5^2 + 3*0.5^3)  )
  
  LOD_rc <- (x[,"n_04"] + x[,"n_30"])*log10((-1 + r_rc)^2*r_rc) + (x[,"n_00"] + x[,"n_34"])*log10(-((-1 + r_rc)*r_rc^2)) + (x[,"n_02"] + x[,"n_32"])*log10(1 + r_rc - r_rc^2) + 
    (x[,"n_01"] + x[,"n_33"])*log10(r_rc*(1 - r_rc + r_rc^2)) + (x[,"n_12"] + x[,"n_22"])*log10(8 - r_rc + r_rc^2) + (x[,"n_10"] + x[,"n_24"])*log10(r_rc*(2 - 4*r_rc + 3*r_rc^2)) + 
    (x[,"n_11"] + x[,"n_23"])*log10(2 - 2*r_rc + 4*r_rc^2 - 3*r_rc^3) + (x[,"n_14"] + x[,"n_20"])*log10(1 - 3*r_rc + 5*r_rc^2 - 3*r_rc^3) + (x[,"n_03"] + x[,"n_31"])*log10(1 - 2*r_rc + 2*r_rc^2 - r_rc^3) + 
    (x[,"n_13"] + x[,"n_21"])*log10(1 + 3*r_rc - 5*r_rc^2 + 3*r_rc^3)  -  ( (x[,"n_04"] + x[,"n_30"])*log10((-1 + 0.5)^2*0.5) + (x[,"n_00"] + x[,"n_34"])*log10(-((-1 + 0.5)*0.5^2)) + (x[,"n_02"] + x[,"n_32"])*log10(1 + 0.5 - 0.5^2) + 
                                                                              (x[,"n_01"] + x[,"n_33"])*log10(0.5*(1 - 0.5 + 0.5^2)) + (x[,"n_12"] + x[,"n_22"])*log10(8 - 0.5 + 0.5^2) + (x[,"n_10"] + x[,"n_24"])*log10(0.5*(2 - 4*0.5 + 3*0.5^2)) + 
                                                                              (x[,"n_11"] + x[,"n_23"])*log10(2 - 2*0.5 + 4*0.5^2 - 3*0.5^3) + (x[,"n_14"] + x[,"n_20"])*log10(1 - 3*0.5 + 5*0.5^2 - 3*0.5^3) + (x[,"n_03"] + x[,"n_31"])*log10(1 - 2*0.5 + 2*0.5^2 - 0.5^3) + 
                                                                              (x[,"n_13"] + x[,"n_21"])*log10(1 + 3*0.5 - 5*0.5^2 + 3*0.5^3) )   
  
  LOD_rr <- (x[,"n_04"] + x[,"n_30"])*log10(-(-1 + r_rr)^3) + (x[,"n_02"] + x[,"n_32"])*log10(-((-1 + r_rr)*r_rr)) + (x[,"n_14"] + x[,"n_20"])*log10((-1 + r_rr)^2*r_rr) + 
    (x[,"n_01"] + x[,"n_33"])*log10(-((-2 + r_rr)*r_rr^2)) + (x[,"n_10"] + x[,"n_24"])*log10(-((-1 + r_rr)*r_rr^2)) + (x[,"n_00"] + x[,"n_34"])*log10(r_rr^3) + 
    (x[,"n_03"] + x[,"n_31"])*log10((-1 + r_rr)^2*(1 + r_rr)) + (x[,"n_11"] + x[,"n_23"])*log10(r_rr*(5 - 5*r_rr + 3*r_rr^2)) + (x[,"n_12"] + x[,"n_22"])*log10(9 - 5*r_rr + 5*r_rr^2) + 
    (x[,"n_13"] + x[,"n_21"])*log10(3 - 4*r_rr + 4*r_rr^2 - 3*r_rr^3)  -  (  (x[,"n_04"] + x[,"n_30"])*log10(-(-1 + 0.5)^3) + (x[,"n_02"] + x[,"n_32"])*log10(-((-1 + 0.5)*0.5)) + (x[,"n_14"] + x[,"n_20"])*log10((-1 + 0.5)^2*0.5) + 
                                                                               (x[,"n_01"] + x[,"n_33"])*log10(-((-2 + 0.5)*0.5^2)) + (x[,"n_10"] + x[,"n_24"])*log10(-((-1 + 0.5)*0.5^2)) + (x[,"n_00"] + x[,"n_34"])*log10(0.5^3) + 
                                                                               (x[,"n_03"] + x[,"n_31"])*log10((-1 + 0.5)^2*(1 + 0.5)) + (x[,"n_11"] + x[,"n_23"])*log10(0.5*(5 - 5*0.5 + 3*0.5^2)) + (x[,"n_12"] + x[,"n_22"])*log10(9 - 5*0.5 + 5*0.5^2) + 
                                                                               (x[,"n_13"] + x[,"n_21"])*log10(3 - 4*0.5 + 4*0.5^2 - 3*0.5^3)  )
  
  
  ## Record the logL:
  logL_cc <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*log(2) + 
    (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_04"] - x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"] - x[,"n_24"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"] - 2*x[,"n_34"])*log(3) + 
    (x[,"n_02"] + x[,"n_32"])*log(5) + (x[,"n_00"] + x[,"n_34"])*log(-(-1 + r_cc)^3) + (x[,"n_02"] + x[,"n_32"])*log(-((-1 + r_cc)*r_cc)) + (x[,"n_10"] + x[,"n_24"])*log((-1 + r_cc)^2*r_cc) + (x[,"n_03"] + x[,"n_31"])*log(-((-2 + r_cc)*r_cc^2)) +
    (x[,"n_14"] + x[,"n_20"])*log(-((-1 + r_cc)*r_cc^2)) + (x[,"n_04"] + x[,"n_30"])*log(r_cc^3) + (x[,"n_01"] + x[,"n_33"])*log((-1 + r_cc)^2*(1 + r_cc)) + (x[,"n_13"] + x[,"n_21"])*log(r_cc*(5 - 5*r_cc + 3*r_cc^2)) + 
    (x[,"n_12"] + x[,"n_22"])*log(9 - 5*r_cc + 5*r_cc^2) + (x[,"n_11"] + x[,"n_23"])*log(3 - 4*r_cc + 4*r_cc^2 - 3*r_cc^3)
  
  logL_cr <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_00"] + x[,"n_34"])*log((-1 + r_cr)^2*r_cr) + (x[,"n_04"] + x[,"n_30"])*log(-((-1 + r_cr)*r_cr^2)) + (x[,"n_02"] + x[,"n_32"])*log(1 + r_cr - r_cr^2) + (x[,"n_03"] + x[,"n_31"])*log(r_cr*(1 - r_cr + r_cr^2)) + 
    (x[,"n_12"] + x[,"n_22"])*log(8 - r_cr + r_cr^2) + (x[,"n_14"] + x[,"n_20"])*log(r_cr*(2 - 4*r_cr + 3*r_cr^2)) + (x[,"n_13"] + x[,"n_21"])*log(2 - 2*r_cr + 4*r_cr^2 - 3*r_cr^3) + 
    (x[,"n_10"] + x[,"n_24"])*log(1 - 3*r_cr + 5*r_cr^2 - 3*r_cr^3) + (x[,"n_01"] + x[,"n_33"])*log(1 - 2*r_cr + 2*r_cr^2 - r_cr^3) + (x[,"n_11"] + x[,"n_23"])*log(1 + 3*r_cr - 5*r_cr^2 + 3*r_cr^3)
  
  logL_mc <- (-3*x[,"n_00"] - 2*x[,"n_01"] - 3*x[,"n_02"] - 2*x[,"n_03"] - 3*x[,"n_04"] - 3*x[,"n_10"] - 2*x[,"n_11"] - 3*x[,"n_12"] - 2*x[,"n_13"] - 3*x[,"n_14"] - 3*x[,"n_20"] - 2*x[,"n_21"] - 3*x[,"n_22"] - 2*x[,"n_23"] - 3*x[,"n_24"] - 3*x[,"n_30"] - 2*x[,"n_31"] - 3*x[,"n_32"] - 2*x[,"n_33"] - 3*x[,"n_34"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_00"] + x[,"n_34"])*log((-1 + r_mc)^2*r_mc) + (x[,"n_04"] + x[,"n_30"])*log(-((-1 + r_mc)*r_mc^2)) + (x[,"n_12"] + x[,"n_22"])*log(14 + 3*r_mc - 3*r_mc^2) + (x[,"n_03"] + x[,"n_31"])*log(r_mc*(2 - 2*r_mc + r_mc^2)) + 
    (x[,"n_14"] + x[,"n_20"])*log(r_mc*(3 - 4*r_mc + 3*r_mc^2)) + (x[,"n_02"] + x[,"n_32"])*log(4 - 3*r_mc + 3*r_mc^2) + (x[,"n_10"] + x[,"n_24"])*log(2 - 4*r_mc + 5*r_mc^2 - 3*r_mc^3) + 
    (x[,"n_13"] + x[,"n_21"])*log(3 - r_mc + 5*r_mc^2 - 3*r_mc^3) + (x[,"n_01"] + x[,"n_33"])*log(1 - r_mc + r_mc^2 - r_mc^3) + (x[,"n_11"] + x[,"n_23"])*log(4 - 4*r_mc^2 + 3*r_mc^3)
  
  logL_mr <- (-3*x[,"n_00"] - 2*x[,"n_01"] - 3*x[,"n_02"] - 2*x[,"n_03"] - 3*x[,"n_04"] - 3*x[,"n_10"] - 2*x[,"n_11"] - 3*x[,"n_12"] - 2*x[,"n_13"] - 3*x[,"n_14"] - 3*x[,"n_20"] - 2*x[,"n_21"] - 3*x[,"n_22"] - 2*x[,"n_23"] - 3*x[,"n_24"] - 3*x[,"n_30"] - 2*x[,"n_31"] - 3*x[,"n_32"] - 2*x[,"n_33"] - 3*x[,"n_34"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_04"] + x[,"n_30"])*log((-1 + r_mr)^2*r_mr) + (x[,"n_00"] + x[,"n_34"])*log(-((-1 + r_mr)*r_mr^2)) + (x[,"n_12"] + x[,"n_22"])*log(14 + 3*r_mr - 3*r_mr^2) + (x[,"n_01"] + x[,"n_33"])*log(r_mr*(2 - 2*r_mr + r_mr^2)) + 
    (x[,"n_10"] + x[,"n_24"])*log(r_mr*(3 - 4*r_mr + 3*r_mr^2)) + (x[,"n_02"] + x[,"n_32"])*log(4 - 3*r_mr + 3*r_mr^2) + (x[,"n_14"] + x[,"n_20"])*log(2 - 4*r_mr + 5*r_mr^2 - 3*r_mr^3) + 
    (x[,"n_11"] + x[,"n_23"])*log(3 - r_mr + 5*r_mr^2 - 3*r_mr^3) + (x[,"n_03"] + x[,"n_31"])*log(1 - r_mr + r_mr^2 - r_mr^3) + (x[,"n_13"] + x[,"n_21"])*log(4 - 4*r_mr^2 + 3*r_mr^3)
  
  logL_rc <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_04"] + x[,"n_30"])*log((-1 + r_rc)^2*r_rc) + (x[,"n_00"] + x[,"n_34"])*log(-((-1 + r_rc)*r_rc^2)) + (x[,"n_02"] + x[,"n_32"])*log(1 + r_rc - r_rc^2) + (x[,"n_01"] + x[,"n_33"])*log(r_rc*(1 - r_rc + r_rc^2)) + 
    (x[,"n_12"] + x[,"n_22"])*log(8 - r_rc + r_rc^2) + (x[,"n_10"] + x[,"n_24"])*log(r_rc*(2 - 4*r_rc + 3*r_rc^2)) + (x[,"n_11"] + x[,"n_23"])*log(2 - 2*r_rc + 4*r_rc^2 - 3*r_rc^3) + 
    (x[,"n_14"] + x[,"n_20"])*log(1 - 3*r_rc + 5*r_rc^2 - 3*r_rc^3) + (x[,"n_03"] + x[,"n_31"])*log(1 - 2*r_rc + 2*r_rc^2 - r_rc^3) + (x[,"n_13"] + x[,"n_21"])*log(1 + 3*r_rc - 5*r_rc^2 + 3*r_rc^3)
  
  logL_rr <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*log(2) + 
    (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_04"] - x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"] - x[,"n_24"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"] - 2*x[,"n_34"])*log(3) + 
    (x[,"n_02"] + x[,"n_32"])*log(5) + (x[,"n_04"] + x[,"n_30"])*log(-(-1 + r_rr)^3) + (x[,"n_02"] + x[,"n_32"])*log(-((-1 + r_rr)*r_rr)) + (x[,"n_14"] + x[,"n_20"])*log((-1 + r_rr)^2*r_rr) + 
    (x[,"n_01"] + x[,"n_33"])*log(-((-2 + r_rr)*r_rr^2)) + (x[,"n_10"] + x[,"n_24"])*log(-((-1 + r_rr)*r_rr^2)) + (x[,"n_00"] + x[,"n_34"])*log(r_rr^3) + (x[,"n_03"] + x[,"n_31"])*log((-1 + r_rr)^2*(1 + r_rr)) + 
    (x[,"n_11"] + x[,"n_23"])*log(r_rr*(5 - 5*r_rr + 3*r_rr^2)) + (x[,"n_12"] + x[,"n_22"])*log(9 - 5*r_rr + 5*r_rr^2) + (x[,"n_13"] + x[,"n_21"])*log(3 - 4*r_rr + 4*r_rr^2 - 3*r_rr^3)
  
  return(list(r_mat=cbind(r_cc, r_cr, r_mc, r_mr, r_rc, r_rr),
              LOD_mat=cbind(LOD_cc, LOD_cr, LOD_mc, LOD_mr, LOD_rc, LOD_rr),
              logL_mat=cbind(logL_cc, logL_cr, logL_mc, logL_mr, logL_rc, logL_rr),
              phasing_strategy="MLL", 
              possible_phases=c("coupling coupling", 
                                "coupling repulsion",
                                "mixed coupling",
                                "mixed repulsion",
                                "repulsion coupling",
                                "repulsion repulsion")))
  
}

#' @rdname r4_functions
#' @noRd
r4_1.3_1.2<-function(x, ncores=1){
  ############################
  ## COUPLING COUPLING
  ############################
  logLcc <- function(r,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-2*n10 - 2*n11 - 2*n12 - 2*n13 - n20 - 2*n21 - 2*n22 - n23 - 2*n30 - 2*n31 - 2*n32 - 2*n33)*log(2) + 
      (-n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(3) + 
      (n10 + n33)*log((-1 + r)^2) + (n12 + n31)*log(-((-3 + r)*r)) + (n20 + n23)*log(-((-1 + r)*r)) + 
      (n13 + n30)*log(r^2) + (n11 + n32)*log(2 - r - r^2) + (n21 + n22)*log(3 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLcc <- function(n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcc,c(0,0.5),n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcc <- function(r){
  #     L <- 2^(-2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"])*3^(-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*((-1 + r)^2)^(x[,"n_10"] + x[,"n_33"])*(-((-3 + r)*r))^(x[,"n_12"] + x[,"n_31"])*(-((-1 + r)*r))^(x[,"n_20"] + x[,"n_23"])*(r^2)^(x[,"n_13"] + x[,"n_30"])*(2 - r - r^2)^(x[,"n_11"] + x[,"n_32"])*(3 - 2*r + 2*r^2)^(x[,"n_21"] + x[,"n_22"])
  #     return(L)}
  
  ############################
  ## COUPLING MIXED
  ############################
  logLcm <- function(r,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-2*n10 - 2*n11 - 2*n12 - 2*n13 - 2*n20 - n21 - n22 - 2*n23 - 2*n30 - 2*n31 - 2*n32 - 2*n33)*log(2) + 
      (-n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(3) + 
      (n10 + n13 + n30 + n33)*log(-((-1 + r)*r)) + (n21 + n22)*log(1 + r - r^2) + (n12 + n31)*log(1 + r^2) + 
      (n11 + n32)*log(2 - 2*r + r^2) + (n20 + n23)*log(1 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLcm <- function(n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcm,c(0,0.5),n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION COUPLING
  ############################
  logLrc <- function(r,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-2*n10 - 2*n11 - 2*n12 - 2*n13 - n20 - 2*n21 - 2*n22 - n23 - 2*n30 - 2*n31 - 2*n32 - 2*n33)*log(2) + 
      (-n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(9) + 
      (n13 + n30)*log(-((-2 + r)*r)) + (n10 + n33)*log(-((-1 + r)*(1 + r))) + (n21 + n22)*log(7 + 2*r - 2*r^2) + 
      (n20 + n23)*log(1 - r + r^2) + (n11 + n12 + n31 + n32)*log(4 - r + r^2)
    return(L)}
  
  inter_logLrc <- function(n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrc,c(0,0.5),n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION MIXED
  ############################
  logLrm <- function(r,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-2*n10 - 2*n11 - 2*n12 - 2*n13 - 2*n20 - n21 - n22 - 2*n23 - 2*n30 - 2*n31 - 2*n32 - 2*n33)*log(2) + 
      (-n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(9) + 
      (n13 + n30)*log((-2 + r)*(-1 + r)) + (n10 + n33)*log(r*(1 + r)) + (n20 + n23)*log(1 + 2*r - 2*r^2) + 
      (n12 + n31)*log(5 - 2*r - r^2) + (n11 + n32)*log(2 + 4*r - r^2) + (n21 + n22)*log(4 - r + r^2)
    return(L)}
  
  inter_logLrm <- function(n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrm,c(0,0.5),n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_cc <- parallel::mcmapply(inter_logLcc,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_rm <- parallel::mcmapply(inter_logLrm,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_cm <-  parallel::mcmapply(inter_logLcm,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  LOD_cc <- (x[,"n_10"] + x[,"n_33"])*log10((-1 + r_cc)^2) + (x[,"n_12"] + x[,"n_31"])*log10(-((-3 + r_cc)*r_cc)) + (x[,"n_20"] + x[,"n_23"])*log10(-((-1 + r_cc)*r_cc)) + 
    (x[,"n_13"] + x[,"n_30"])*log10(r_cc^2) + (x[,"n_11"] + x[,"n_32"])*log10(2 - r_cc - r_cc^2) + (x[,"n_21"] + x[,"n_22"])*log10(3 - 2*r_cc + 2*r_cc^2) - 
    (x[,"n_10"] + x[,"n_33"])*log10((-1 + 0.5)^2) - (x[,"n_12"] + x[,"n_31"])*log10(-((-3 + 0.5)*0.5)) - (x[,"n_20"] + x[,"n_23"])*log10(-((-1 + 0.5)*0.5)) - (x[,"n_13"] + x[,"n_30"])*log10(0.5^2) - 
    (x[,"n_11"] + x[,"n_32"])*log10(2 - 0.5 - 0.5^2) - (x[,"n_21"] + x[,"n_22"])*log10(3 - 2*0.5 + 2*0.5^2)
  
  LOD_cm <- (x[,"n_10"] + x[,"n_13"] + x[,"n_30"] + x[,"n_33"])*log10(-((-1 + r_cm)*r_cm)) + (x[,"n_21"] + x[,"n_22"])*log10(1 + r_cm - r_cm^2) + (x[,"n_12"] + x[,"n_31"])*log10(1 + r_cm^2) + 
    (x[,"n_11"] + x[,"n_32"])*log10(2 - 2*r_cm + r_cm^2) + (x[,"n_20"] + x[,"n_23"])*log10(1 - 2*r_cm + 2*r_cm^2) - (x[,"n_10"] + x[,"n_13"] + x[,"n_30"] + x[,"n_33"])*log10(-((-1 + 0.5)*0.5)) - 
    (x[,"n_21"] + x[,"n_22"])*log10(1 + 0.5 - 0.5^2) - (x[,"n_12"] + x[,"n_31"])*log10(1 + 0.5^2) - (x[,"n_11"] + x[,"n_32"])*log10(2 - 2*0.5 + 0.5^2) - (x[,"n_20"] + x[,"n_23"])*log10(1 - 2*0.5 + 2*0.5^2)
  
  LOD_rc <- (x[,"n_13"] + x[,"n_30"])*log10(-((-2 + r_rc)*r_rc)) + (x[,"n_10"] + x[,"n_33"])*log10(-((-1 + r_rc)*(1 + r_rc))) + (x[,"n_21"] + x[,"n_22"])*log10(7 + 2*r_rc - 2*r_rc^2) + 
    (x[,"n_20"] + x[,"n_23"])*log10(1 - r_rc + r_rc^2) + (x[,"n_11"] + x[,"n_12"] + x[,"n_31"] + x[,"n_32"])*log10(4 - r_rc + r_rc^2) - (x[,"n_13"] + x[,"n_30"])*log10(-((-2 + 0.5)*0.5)) - 
    (x[,"n_10"] + x[,"n_33"])*log10(-((-1 + 0.5)*(1 + 0.5))) - (x[,"n_21"] + x[,"n_22"])*log10(7 + 2*0.5 - 2*0.5^2) - (x[,"n_20"] + x[,"n_23"])*log10(1 - 0.5 + 0.5^2) - 
    (x[,"n_11"] + x[,"n_12"] + x[,"n_31"] + x[,"n_32"])*log10(4 - 0.5 + 0.5^2)
  
  LOD_rm <- (x[,"n_13"] + x[,"n_30"])*log10((-2 + r_rm)*(-1 + r_rm)) + (x[,"n_10"] + x[,"n_33"])*log10(r_rm*(1 + r_rm)) + (x[,"n_20"] + x[,"n_23"])*log10(1 + 2*r_rm - 2*r_rm^2) + 
    (x[,"n_12"] + x[,"n_31"])*log10(5 - 2*r_rm - r_rm^2) + (x[,"n_11"] + x[,"n_32"])*log10(2 + 4*r_rm - r_rm^2) + (x[,"n_21"] + x[,"n_22"])*log10(4 - r_rm + r_rm^2) - 
    (x[,"n_13"] + x[,"n_30"])*log10((-2 + 0.5)*(-1 + 0.5)) - (x[,"n_10"] + x[,"n_33"])*log10(0.5*(1 + 0.5)) - (x[,"n_20"] + x[,"n_23"])*log10(1 + 2*0.5 - 2*0.5^2) - 
    (x[,"n_12"] + x[,"n_31"])*log10(5 - 2*0.5 - 0.5^2) - (x[,"n_11"] + x[,"n_32"])*log10(2 + 4*0.5 - 0.5^2) - (x[,"n_21"] + x[,"n_22"])*log10(4 - 0.5 + 0.5^2)
  
  ## Determine the logL:
  logL_cc <- (-2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(3) + 
    (x[,"n_10"] + x[,"n_33"])*log((-1 + r_cc)^2) + (x[,"n_12"] + x[,"n_31"])*log(-((-3 + r_cc)*r_cc)) + (x[,"n_20"] + x[,"n_23"])*log(-((-1 + r_cc)*r_cc)) + 
    (x[,"n_13"] + x[,"n_30"])*log(r_cc^2) + (x[,"n_11"] + x[,"n_32"])*log(2 - r_cc - r_cc^2) + (x[,"n_21"] + x[,"n_22"])*log(3 - 2*r_cc + 2*r_cc^2)
  
  logL_cm <- (-2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - 2*x[,"n_23"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(3) + 
    (x[,"n_10"] + x[,"n_13"] + x[,"n_30"] + x[,"n_33"])*log(-((-1 + r_cm)*r_cm)) + (x[,"n_21"] + x[,"n_22"])*log(1 + r_cm - r_cm^2) + (x[,"n_12"] + x[,"n_31"])*log(1 + r_cm^2) + 
    (x[,"n_11"] + x[,"n_32"])*log(2 - 2*r_cm + r_cm^2) + (x[,"n_20"] + x[,"n_23"])*log(1 - 2*r_cm + 2*r_cm^2)
  
  logL_rc <- (-2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(9) + 
    (x[,"n_13"] + x[,"n_30"])*log(-((-2 + r_rc)*r_rc)) + (x[,"n_10"] + x[,"n_33"])*log(-((-1 + r_rc)*(1 + r_rc))) + (x[,"n_21"] + x[,"n_22"])*log(7 + 2*r_rc - 2*r_rc^2) + 
    (x[,"n_20"] + x[,"n_23"])*log(1 - r_rc + r_rc^2) + (x[,"n_11"] + x[,"n_12"] + x[,"n_31"] + x[,"n_32"])*log(4 - r_rc + r_rc^2)
  
  logL_rm <- (-2*x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - 2*x[,"n_23"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(9) + 
    (x[,"n_13"] + x[,"n_30"])*log((-2 + r_rm)*(-1 + r_rm)) + (x[,"n_10"] + x[,"n_33"])*log(r_rm*(1 + r_rm)) + (x[,"n_20"] + x[,"n_23"])*log(1 + 2*r_rm - 2*r_rm^2) + 
    (x[,"n_12"] + x[,"n_31"])*log(5 - 2*r_rm - r_rm^2) + (x[,"n_11"] + x[,"n_32"])*log(2 + 4*r_rm - r_rm^2) + (x[,"n_21"] + x[,"n_22"])*log(4 - r_rm + r_rm^2)
  
  return(list(r_mat=cbind(r_cc, r_cm, r_rc, r_rm),
              LOD_mat=cbind(LOD_cc, LOD_cm, LOD_rc, LOD_rm),
              logL_mat=cbind(logL_cc, logL_cm, logL_rc, logL_rm),
              phasing_strategy="MLL", 
              possible_phases=c("coupling coupling", 
                                "coupling mixed",
                                "repulsion coupling",
                                "repulsion mixed")))
}

#' @rdname r4_functions
#' @noRd
r4_1.2_2.2<-function(x, ncores=1){
  ############################
  ## COUPLING COUPLING
  ############################
  logLcc <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34)*log(2) + 
      (-2*n00 - 2*n01 - 2*n02 - 2*n03 - 2*n04 - n10 - 2*n11 - 2*n12 - 2*n13 - n14 - n20 - 2*n21 - 2*n22 - 2*n23 - n24 - 2*n30 - 2*n31 - 2*n32 - 2*n33 - 2*n34)*log(3) + 
      (n02 + n32)*log(5) + (n00 + n34)*log(-(-1 + r)^3) + (n02 + n32)*log(-((-1 + r)*r)) + (n10 + n24)*log((-1 + r)^2*r) + 
      (n03 + n31)*log(-((-2 + r)*r^2)) + (n14 + n20)*log(-((-1 + r)*r^2)) + (n04 + n30)*log(r^3) + (n01 + n33)*log((-1 + r)^2*(1 + r)) + 
      (n13 + n21)*log(r*(5 - 5*r + 3*r^2)) + (n12 + n22)*log(9 - 5*r + 5*r^2) + (n11 + n23)*log(3 - 4*r + 4*r^2 - 3*r^3)
    return(L)}
  
  inter_logLcc <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLcc,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## COUPLING MIXED
  ############################
  logLcm <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-3*n00 - 2*n01 - 3*n02 - 2*n03 - 3*n04 - 3*n10 - 2*n11 - 3*n12 - 2*n13 - 3*n14 - 3*n20 - 2*n21 - 3*n22 - 2*n23 - 3*n24 - 3*n30 - 2*n31 - 3*n32 - 2*n33 - 3*n34)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n00 + n34)*log((-1 + r)^2*r) + (n04 + n30)*log(-((-1 + r)*r^2)) + (n12 + n22)*log(14 + 3*r - 3*r^2) + (n03 + n31)*log(r*(2 - 2*r + r^2)) + 
      (n14 + n20)*log(r*(3 - 4*r + 3*r^2)) + (n02 + n32)*log(4 - 3*r + 3*r^2) + (n10 + n24)*log(2 - 4*r + 5*r^2 - 3*r^3) + (n13 + n21)*log(3 - r + 5*r^2 - 3*r^3) + 
      (n01 + n33)*log(1 - r + r^2 - r^3) + (n11 + n23)*log(4 - 4*r^2 + 3*r^3)
    return(L)}
  
  inter_logLcm <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLcm,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## COUPLING REPULSION
  ############################
  logLcr <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n04 + n30)*log((-1 + r)^2*r) + (n00 + n34)*log(-((-1 + r)*r^2)) + (n02 + n32)*log(1 + r - r^2) + (n01 + n33)*log(r*(1 - r + r^2)) + 
      (n12 + n22)*log(8 - r + r^2) + (n10 + n24)*log(r*(2 - 4*r + 3*r^2)) + (n11 + n23)*log(2 - 2*r + 4*r^2 - 3*r^3) + 
      (n14 + n20)*log(1 - 3*r + 5*r^2 - 3*r^3) + (n03 + n31)*log(1 - 2*r + 2*r^2 - r^3) + (n13 + n21)*log(1 + 3*r - 5*r^2 + 3*r^3)
    return(L)}
  
  inter_logLcr <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLcr,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION COUPLING
  ############################
  logLrc <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n00 + n34)*log((-1 + r)^2*r) + (n04 + n30)*log(-((-1 + r)*r^2)) + (n02 + n32)*log(1 + r - r^2) + (n03 + n31)*log(r*(1 - r + r^2)) + 
      (n12 + n22)*log(8 - r + r^2) + (n14 + n20)*log(r*(2 - 4*r + 3*r^2)) + (n13 + n21)*log(2 - 2*r + 4*r^2 - 3*r^3) + (n10 + n24)*log(1 - 3*r + 5*r^2 - 3*r^3) + 
      (n01 + n33)*log(1 - 2*r + 2*r^2 - r^3) + (n11 + n23)*log(1 + 3*r - 5*r^2 + 3*r^3)
    return(L)}
  
  inter_logLrc <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLrc,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION MIXED
  ############################
  logLrm <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-3*n00 - 2*n01 - 3*n02 - 2*n03 - 3*n04 - 3*n10 - 2*n11 - 3*n12 - 2*n13 - 3*n14 - 3*n20 - 2*n21 - 3*n22 - 2*n23 - 3*n24 - 3*n30 - 2*n31 - 3*n32 - 2*n33 - 3*n34)*log(2) + 
      (-n00 - n01 - n02 - n03 - n04 - n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n04 + n30)*log((-1 + r)^2*r) + (n00 + n34)*log(-((-1 + r)*r^2)) + (n12 + n22)*log(14 + 3*r - 3*r^2) + (n01 + n33)*log(r*(2 - 2*r + r^2)) + 
      (n10 + n24)*log(r*(3 - 4*r + 3*r^2)) + (n02 + n32)*log(4 - 3*r + 3*r^2) + (n14 + n20)*log(2 - 4*r + 5*r^2 - 3*r^3) + (n11 + n23)*log(3 - r + 5*r^2 - 3*r^3) + 
      (n03 + n31)*log(1 - r + r^2 - r^3) + (n13 + n21)*log(4 - 4*r^2 + 3*r^3)
    return(L)}
  
  inter_logLrm <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLrm,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  ############################
  ## REPULSION REPULSION
  ############################
  logLrr <- function(r,n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n00 - n01 - 2*n02 - n03 - 2*n04 - 2*n10 - n11 - 2*n12 - n13 - 2*n14 - 2*n20 - n21 - 2*n22 - n23 - 2*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34)*log(2) + 
      (-2*n00 - 2*n01 - 2*n02 - 2*n03 - 2*n04 - n10 - 2*n11 - 2*n12 - 2*n13 - n14 - n20 - 2*n21 - 2*n22 - 2*n23 - n24 - 2*n30 - 2*n31 - 2*n32 - 2*n33 - 2*n34)*log(3) + 
      (n02 + n32)*log(5) + (n04 + n30)*log(-(-1 + r)^3) + (n02 + n32)*log(-((-1 + r)*r)) + (n14 + n20)*log((-1 + r)^2*r) + (n01 + n33)*log(-((-2 + r)*r^2)) + 
      (n10 + n24)*log(-((-1 + r)*r^2)) + (n00 + n34)*log(r^3) + (n03 + n31)*log((-1 + r)^2*(1 + r)) + (n11 + n23)*log(r*(5 - 5*r + 3*r^2)) + 
      (n12 + n22)*log(9 - 5*r + 5*r^2) + (n13 + n21)*log(3 - 4*r + 4*r^2 - 3*r^3)
    return(L)}
  
  inter_logLrr <- function(n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLrr,c(0,0.5),n00,n01,n02,n03,n04,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_cc <- parallel::mcmapply(inter_logLcc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                             x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                             mc.cores = ncores)
  
  r_cm <- parallel::mcmapply(inter_logLcm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                             x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                             mc.cores = ncores)
  
  r_cr <- parallel::mcmapply(inter_logLcr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                             x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                             mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                              x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                              mc.cores = ncores)
  
  r_rm <-  parallel::mcmapply(inter_logLrm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                              x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                              mc.cores = ncores)
  
  r_rr <-  parallel::mcmapply(inter_logLrr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_04"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],
                              x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                              mc.cores = ncores)
  
  LOD_cc <- (x[,"n_00"] + x[,"n_34"])*log10(-(-1 + r_cc)^3) + (x[,"n_02"] + x[,"n_32"])*log10(-((-1 + r_cc)*r_cc)) + (x[,"n_10"] + x[,"n_24"])*log10((-1 + r_cc)^2*r_cc) + 
    (x[,"n_03"] + x[,"n_31"])*log10(-((-2 + r_cc)*r_cc^2)) + (x[,"n_14"] + x[,"n_20"])*log10(-((-1 + r_cc)*r_cc^2)) + (x[,"n_04"] + x[,"n_30"])*log10(r_cc^3) + 
    (x[,"n_01"] + x[,"n_33"])*log10((-1 + r_cc)^2*(1 + r_cc)) + (x[,"n_13"] + x[,"n_21"])*log10(r_cc*(5 - 5*r_cc + 3*r_cc^2)) + (x[,"n_12"] + x[,"n_22"])*log10(9 - 5*r_cc + 5*r_cc^2) + 
    (x[,"n_11"] + x[,"n_23"])*log10(3 - 4*r_cc + 4*r_cc^2 - 3*r_cc^3) - (x[,"n_00"] + x[,"n_34"])*log10(-(-1 + 0.5)^3) - (x[,"n_02"] + x[,"n_32"])*log10(-((-1 + 0.5)*0.5)) - 
    (x[,"n_10"] + x[,"n_24"])*log10((-1 + 0.5)^2*0.5) - (x[,"n_03"] + x[,"n_31"])*log10(-((-2 + 0.5)*0.5^2)) - (x[,"n_14"] + x[,"n_20"])*log10(-((-1 + 0.5)*0.5^2)) - 
    (x[,"n_04"] + x[,"n_30"])*log10(0.5^3) - (x[,"n_01"] + x[,"n_33"])*log10((-1 + 0.5)^2*(1 + 0.5)) - (x[,"n_13"] + x[,"n_21"])*log10(0.5*(5 - 5*0.5 + 3*0.5^2)) - 
    (x[,"n_12"] + x[,"n_22"])*log10(9 - 5*0.5 + 5*0.5^2) - (x[,"n_11"] + x[,"n_23"])*log10(3 - 4*0.5 + 4*0.5^2 - 3*0.5^3)
  
  LOD_cr <- (x[,"n_04"] + x[,"n_30"])*log10((-1 + r_cr)^2*r_cr) + (x[,"n_00"] + x[,"n_34"])*log10(-((-1 + r_cr)*r_cr^2)) + (x[,"n_02"] + x[,"n_32"])*log10(1 + r_cr - r_cr^2) + 
    (x[,"n_01"] + x[,"n_33"])*log10(r_cr*(1 - r_cr + r_cr^2)) + (x[,"n_12"] + x[,"n_22"])*log10(8 - r_cr + r_cr^2) + (x[,"n_10"] + x[,"n_24"])*log10(r_cr*(2 - 4*r_cr + 3*r_cr^2)) +
    (x[,"n_11"] + x[,"n_23"])*log10(2 - 2*r_cr + 4*r_cr^2 - 3*r_cr^3) + (x[,"n_14"] + x[,"n_20"])*log10(1 - 3*r_cr + 5*r_cr^2 - 3*r_cr^3) + (x[,"n_03"] + x[,"n_31"])*log10(1 - 2*r_cr + 2*r_cr^2 - r_cr^3) +
    (x[,"n_13"] + x[,"n_21"])*log10(1 + 3*r_cr - 5*r_cr^2 + 3*r_cr^3) - (x[,"n_04"] + x[,"n_30"])*log10((-1 + 0.5)^2*0.5) - (x[,"n_00"] + x[,"n_34"])*log10(-((-1 + 0.5)*0.5^2)) - 
    (x[,"n_02"] + x[,"n_32"])*log10(1 + 0.5 - 0.5^2) - (x[,"n_01"] + x[,"n_33"])*log10(0.5*(1 - 0.5 + 0.5^2)) - (x[,"n_12"] + x[,"n_22"])*log10(8 - 0.5 + 0.5^2) - 
    (x[,"n_10"] + x[,"n_24"])*log10(0.5*(2 - 4*0.5 + 3*0.5^2)) - (x[,"n_11"] + x[,"n_23"])*log10(2 - 2*0.5 + 4*0.5^2 - 3*0.5^3) - 
    (x[,"n_14"] + x[,"n_20"])*log10(1 - 3*0.5 + 5*0.5^2 - 3*0.5^3) - (x[,"n_03"] + x[,"n_31"])*log10(1 - 2*0.5 + 2*0.5^2 - 0.5^3) - (x[,"n_13"] + x[,"n_21"])*log10(1 + 3*0.5 - 5*0.5^2 + 3*0.5^3)
  
  LOD_cm <- (x[,"n_00"] + x[,"n_34"])*log10((-1 + r_cm)^2*r_cm) + (x[,"n_04"] + x[,"n_30"])*log10(-((-1 + r_cm)*r_cm^2)) + (x[,"n_12"] + x[,"n_22"])*log10(14 + 3*r_cm - 3*r_cm^2) + 
    (x[,"n_03"] + x[,"n_31"])*log10(r_cm*(2 - 2*r_cm + r_cm^2)) + (x[,"n_14"] + x[,"n_20"])*log10(r_cm*(3 - 4*r_cm + 3*r_cm^2)) + (x[,"n_02"] + x[,"n_32"])*log10(4 - 3*r_cm + 3*r_cm^2) + 
    (x[,"n_10"] + x[,"n_24"])*log10(2 - 4*r_cm + 5*r_cm^2 - 3*r_cm^3) + (x[,"n_13"] + x[,"n_21"])*log10(3 - r_cm + 5*r_cm^2 - 3*r_cm^3) + (x[,"n_01"] + x[,"n_33"])*log10(1 - r_cm + r_cm^2 - r_cm^3) +
    (x[,"n_11"] + x[,"n_23"])*log10(4 - 4*r_cm^2 + 3*r_cm^3) - (x[,"n_00"] + x[,"n_34"])*log10((-1 + 0.5)^2*0.5) - (x[,"n_04"] + x[,"n_30"])*log10(-((-1 + 0.5)*0.5^2)) - 
    (x[,"n_12"] + x[,"n_22"])*log10(14 + 3*0.5 - 3*0.5^2) - (x[,"n_03"] + x[,"n_31"])*log10(0.5*(2 - 2*0.5 + 0.5^2)) - (x[,"n_14"] + x[,"n_20"])*log10(0.5*(3 - 4*0.5 + 3*0.5^2)) - 
    (x[,"n_02"] + x[,"n_32"])*log10(4 - 3*0.5 + 3*0.5^2) - (x[,"n_10"] + x[,"n_24"])*log10(2 - 4*0.5 + 5*0.5^2 - 3*0.5^3) - (x[,"n_13"] + x[,"n_21"])*log10(3 - 0.5 + 5*0.5^2 - 3*0.5^3) - 
    (x[,"n_01"] + x[,"n_33"])*log10(1 - 0.5 + 0.5^2 - 0.5^3) - (x[,"n_11"] + x[,"n_23"])*log10(4 - 4*0.5^2 + 3*0.5^3)
  
  LOD_rc <- (x[,"n_00"] + x[,"n_34"])*log10((-1 + r_rc)^2*r_rc) + (x[,"n_04"] + x[,"n_30"])*log10(-((-1 + r_rc)*r_rc^2)) + (x[,"n_02"] + x[,"n_32"])*log10(1 + r_rc - r_rc^2) + 
    (x[,"n_03"] + x[,"n_31"])*log10(r_rc*(1 - r_rc + r_rc^2)) + (x[,"n_12"] + x[,"n_22"])*log10(8 - r_rc + r_rc^2) + (x[,"n_14"] + x[,"n_20"])*log10(r_rc*(2 - 4*r_rc + 3*r_rc^2)) + 
    (x[,"n_13"] + x[,"n_21"])*log10(2 - 2*r_rc + 4*r_rc^2 - 3*r_rc^3) + (x[,"n_10"] + x[,"n_24"])*log10(1 - 3*r_rc + 5*r_rc^2 - 3*r_rc^3) + (x[,"n_01"] + x[,"n_33"])*log10(1 - 2*r_rc + 2*r_rc^2 - r_rc^3) + 
    (x[,"n_11"] + x[,"n_23"])*log10(1 + 3*r_rc - 5*r_rc^2 + 3*r_rc^3) - (x[,"n_00"] + x[,"n_34"])*log10((-1 + 0.5)^2*0.5) - (x[,"n_04"] + x[,"n_30"])*log10(-((-1 + 0.5)*0.5^2)) - 
    (x[,"n_02"] + x[,"n_32"])*log10(1 + 0.5 - 0.5^2) - (x[,"n_03"] + x[,"n_31"])*log10(0.5*(1 - 0.5 + 0.5^2)) - (x[,"n_12"] + x[,"n_22"])*log10(8 - 0.5 + 0.5^2) - 
    (x[,"n_14"] + x[,"n_20"])*log10(0.5*(2 - 4*0.5 + 3*0.5^2)) - (x[,"n_13"] + x[,"n_21"])*log10(2 - 2*0.5 + 4*0.5^2 - 3*0.5^3) - 
    (x[,"n_10"] + x[,"n_24"])*log10(1 - 3*0.5 + 5*0.5^2 - 3*0.5^3) - (x[,"n_01"] + x[,"n_33"])*log10(1 - 2*0.5 + 2*0.5^2 - 0.5^3) - 
    (x[,"n_11"] + x[,"n_23"])*log10(1 + 3*0.5 - 5*0.5^2 + 3*0.5^3)
  
  LOD_rm <- (x[,"n_04"] + x[,"n_30"])*log10((-1 + r_rm)^2*r_rm) + (x[,"n_00"] + x[,"n_34"])*log10(-((-1 + r_rm)*r_rm^2)) + (x[,"n_12"] + x[,"n_22"])*log10(14 + 3*r_rm - 3*r_rm^2) + 
    (x[,"n_01"] + x[,"n_33"])*log10(r_rm*(2 - 2*r_rm + r_rm^2)) + (x[,"n_10"] + x[,"n_24"])*log10(r_rm*(3 - 4*r_rm + 3*r_rm^2)) + (x[,"n_02"] + x[,"n_32"])*log10(4 - 3*r_rm + 3*r_rm^2) + 
    (x[,"n_14"] + x[,"n_20"])*log10(2 - 4*r_rm + 5*r_rm^2 - 3*r_rm^3) + (x[,"n_11"] + x[,"n_23"])*log10(3 - r_rm + 5*r_rm^2 - 3*r_rm^3) + (x[,"n_03"] + x[,"n_31"])*log10(1 - r_rm + r_rm^2 - r_rm^3) + 
    (x[,"n_13"] + x[,"n_21"])*log10(4 - 4*r_rm^2 + 3*r_rm^3) - (x[,"n_04"] + x[,"n_30"])*log10((-1 + 0.5)^2*0.5) - (x[,"n_00"] + x[,"n_34"])*log10(-((-1 + 0.5)*0.5^2)) - 
    (x[,"n_12"] + x[,"n_22"])*log10(14 + 3*0.5 - 3*0.5^2) - (x[,"n_01"] + x[,"n_33"])*log10(0.5*(2 - 2*0.5 + 0.5^2)) - (x[,"n_10"] + x[,"n_24"])*log10(0.5*(3 - 4*0.5 + 3*0.5^2)) - 
    (x[,"n_02"] + x[,"n_32"])*log10(4 - 3*0.5 + 3*0.5^2) - (x[,"n_14"] + x[,"n_20"])*log10(2 - 4*0.5 + 5*0.5^2 - 3*0.5^3) - (x[,"n_11"] + x[,"n_23"])*log10(3 - 0.5 + 5*0.5^2 - 3*0.5^3) - 
    (x[,"n_03"] + x[,"n_31"])*log10(1 - 0.5 + 0.5^2 - 0.5^3) - (x[,"n_13"] + x[,"n_21"])*log10(4 - 4*0.5^2 + 3*0.5^3)
  
  LOD_rr <- (x[,"n_04"] + x[,"n_30"])*log10(-(-1 + r_rr)^3) + (x[,"n_02"] + x[,"n_32"])*log10(-((-1 + r_rr)*r_rr)) + (x[,"n_14"] + x[,"n_20"])*log10((-1 + r_rr)^2*r_rr) + 
    (x[,"n_01"] + x[,"n_33"])*log10(-((-2 + r_rr)*r_rr^2)) + (x[,"n_10"] + x[,"n_24"])*log10(-((-1 + r_rr)*r_rr^2)) + (x[,"n_00"] + x[,"n_34"])*log10(r_rr^3) + (x[,"n_03"] + x[,"n_31"])*log10((-1 + r_rr)^2*(1 + r_rr)) + 
    (x[,"n_11"] + x[,"n_23"])*log10(r_rr*(5 - 5*r_rr + 3*r_rr^2)) + (x[,"n_12"] + x[,"n_22"])*log10(9 - 5*r_rr + 5*r_rr^2) + (x[,"n_13"] + x[,"n_21"])*log10(3 - 4*r_rr + 4*r_rr^2 - 3*r_rr^3) - 
    (x[,"n_04"] + x[,"n_30"])*log10(-(-1 + 0.5)^3) - (x[,"n_02"] + x[,"n_32"])*log10(-((-1 + 0.5)*0.5)) - (x[,"n_14"] + x[,"n_20"])*log10((-1 + 0.5)^2*0.5) - 
    (x[,"n_01"] + x[,"n_33"])*log10(-((-2 + 0.5)*0.5^2)) - (x[,"n_10"] + x[,"n_24"])*log10(-((-1 + 0.5)*0.5^2)) - (x[,"n_00"] + x[,"n_34"])*log10(0.5^3) - 
    (x[,"n_03"] + x[,"n_31"])*log10((-1 + 0.5)^2*(1 + 0.5)) - (x[,"n_11"] + x[,"n_23"])*log10(0.5*(5 - 5*0.5 + 3*0.5^2)) - (x[,"n_12"] + x[,"n_22"])*log10(9 - 5*0.5 + 5*0.5^2) - 
    (x[,"n_13"] + x[,"n_21"])*log10(3 - 4*0.5 + 4*0.5^2 - 3*0.5^3)
  
  ## Record the logL also:
  logL_cc <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*log(2) + 
    (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_04"] - x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"] - x[,"n_24"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"] - 2*x[,"n_34"])*log(3) + 
    (x[,"n_02"] + x[,"n_32"])*log(5) + (x[,"n_00"] + x[,"n_34"])*log(-(-1 + r_cc)^3) + (x[,"n_02"] + x[,"n_32"])*log(-((-1 + r_cc)*r_cc)) + (x[,"n_10"] + x[,"n_24"])*log((-1 + r_cc)^2*r_cc) + 
    (x[,"n_03"] + x[,"n_31"])*log(-((-2 + r_cc)*r_cc^2)) + (x[,"n_14"] + x[,"n_20"])*log(-((-1 + r_cc)*r_cc^2)) + (x[,"n_04"] + x[,"n_30"])*log(r_cc^3) + (x[,"n_01"] + x[,"n_33"])*log((-1 + r_cc)^2*(1 + r_cc)) + 
    (x[,"n_13"] + x[,"n_21"])*log(r_cc*(5 - 5*r_cc + 3*r_cc^2)) + (x[,"n_12"] + x[,"n_22"])*log(9 - 5*r_cc + 5*r_cc^2) + (x[,"n_11"] + x[,"n_23"])*log(3 - 4*r_cc + 4*r_cc^2 - 3*r_cc^3)
  
  logL_cm <- (-3*x[,"n_00"] - 2*x[,"n_01"] - 3*x[,"n_02"] - 2*x[,"n_03"] - 3*x[,"n_04"] - 3*x[,"n_10"] - 2*x[,"n_11"] - 3*x[,"n_12"] - 2*x[,"n_13"] - 3*x[,"n_14"] - 3*x[,"n_20"] - 2*x[,"n_21"] - 3*x[,"n_22"] - 2*x[,"n_23"] - 3*x[,"n_24"] - 3*x[,"n_30"] - 2*x[,"n_31"] - 3*x[,"n_32"] - 2*x[,"n_33"] - 3*x[,"n_34"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_00"] + x[,"n_34"])*log((-1 + r_cm)^2*r_cm) + (x[,"n_04"] + x[,"n_30"])*log(-((-1 + r_cm)*r_cm^2)) + (x[,"n_12"] + x[,"n_22"])*log(14 + 3*r_cm - 3*r_cm^2) + (x[,"n_03"] + x[,"n_31"])*log(r_cm*(2 - 2*r_cm + r_cm^2)) + 
    (x[,"n_14"] + x[,"n_20"])*log(r_cm*(3 - 4*r_cm + 3*r_cm^2)) + (x[,"n_02"] + x[,"n_32"])*log(4 - 3*r_cm + 3*r_cm^2) + (x[,"n_10"] + x[,"n_24"])*log(2 - 4*r_cm + 5*r_cm^2 - 3*r_cm^3) + (x[,"n_13"] + x[,"n_21"])*log(3 - r_cm + 5*r_cm^2 - 3*r_cm^3) + 
    (x[,"n_01"] + x[,"n_33"])*log(1 - r_cm + r_cm^2 - r_cm^3) + (x[,"n_11"] + x[,"n_23"])*log(4 - 4*r_cm^2 + 3*r_cm^3)
  
  logL_cr <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_04"] + x[,"n_30"])*log((-1 + r_cr)^2*r_cr) + (x[,"n_00"] + x[,"n_34"])*log(-((-1 + r_cr)*r_cr^2)) + (x[,"n_02"] + x[,"n_32"])*log(1 + r_cr - r_cr^2) + (x[,"n_01"] + x[,"n_33"])*log(r_cr*(1 - r_cr + r_cr^2)) + 
    (x[,"n_12"] + x[,"n_22"])*log(8 - r_cr + r_cr^2) + (x[,"n_10"] + x[,"n_24"])*log(r_cr*(2 - 4*r_cr + 3*r_cr^2)) + (x[,"n_11"] + x[,"n_23"])*log(2 - 2*r_cr + 4*r_cr^2 - 3*r_cr^3) + 
    (x[,"n_14"] + x[,"n_20"])*log(1 - 3*r_cr + 5*r_cr^2 - 3*r_cr^3) + (x[,"n_03"] + x[,"n_31"])*log(1 - 2*r_cr + 2*r_cr^2 - r_cr^3) + (x[,"n_13"] + x[,"n_21"])*log(1 + 3*r_cr - 5*r_cr^2 + 3*r_cr^3)
  
  logL_rc <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_00"] + x[,"n_34"])*log((-1 + r_rc)^2*r_rc) + (x[,"n_04"] + x[,"n_30"])*log(-((-1 + r_rc)*r_rc^2)) + (x[,"n_02"] + x[,"n_32"])*log(1 + r_rc - r_rc^2) + (x[,"n_03"] + x[,"n_31"])*log(r_rc*(1 - r_rc + r_rc^2)) + 
    (x[,"n_12"] + x[,"n_22"])*log(8 - r_rc + r_rc^2) + (x[,"n_14"] + x[,"n_20"])*log(r_rc*(2 - 4*r_rc + 3*r_rc^2)) + (x[,"n_13"] + x[,"n_21"])*log(2 - 2*r_rc + 4*r_rc^2 - 3*r_rc^3) + (x[,"n_10"] + x[,"n_24"])*log(1 - 3*r_rc + 5*r_rc^2 - 3*r_rc^3) + 
    (x[,"n_01"] + x[,"n_33"])*log(1 - 2*r_rc + 2*r_rc^2 - r_rc^3) + (x[,"n_11"] + x[,"n_23"])*log(1 + 3*r_rc - 5*r_rc^2 + 3*r_rc^3)
  
  logL_rm <- (-3*x[,"n_00"] - 2*x[,"n_01"] - 3*x[,"n_02"] - 2*x[,"n_03"] - 3*x[,"n_04"] - 3*x[,"n_10"] - 2*x[,"n_11"] - 3*x[,"n_12"] - 2*x[,"n_13"] - 3*x[,"n_14"] - 3*x[,"n_20"] - 2*x[,"n_21"] - 3*x[,"n_22"] - 2*x[,"n_23"] - 3*x[,"n_24"] - 3*x[,"n_30"] - 2*x[,"n_31"] - 3*x[,"n_32"] - 2*x[,"n_33"] - 3*x[,"n_34"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_04"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_04"] + x[,"n_30"])*log((-1 + r_rm)^2*r_rm) + (x[,"n_00"] + x[,"n_34"])*log(-((-1 + r_rm)*r_rm^2)) + (x[,"n_12"] + x[,"n_22"])*log(14 + 3*r_rm - 3*r_rm^2) + (x[,"n_01"] + x[,"n_33"])*log(r_rm*(2 - 2*r_rm + r_rm^2)) + 
    (x[,"n_10"] + x[,"n_24"])*log(r_rm*(3 - 4*r_rm + 3*r_rm^2)) + (x[,"n_02"] + x[,"n_32"])*log(4 - 3*r_rm + 3*r_rm^2) + (x[,"n_14"] + x[,"n_20"])*log(2 - 4*r_rm + 5*r_rm^2 - 3*r_rm^3) + (x[,"n_11"] + x[,"n_23"])*log(3 - r_rm + 5*r_rm^2 - 3*r_rm^3) + 
    (x[,"n_03"] + x[,"n_31"])*log(1 - r_rm + r_rm^2 - r_rm^3) + (x[,"n_13"] + x[,"n_21"])*log(4 - 4*r_rm^2 + 3*r_rm^3)
  
  logL_rr <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_03"] - 2*x[,"n_04"] - 2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"] - x[,"n_23"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*log(2) + 
    (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_04"] - x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - 2*x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"] - x[,"n_24"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"] - 2*x[,"n_34"])*log(3) + 
    (x[,"n_02"] + x[,"n_32"])*log(5) + (x[,"n_04"] + x[,"n_30"])*log(-(-1 + r_rr)^3) + (x[,"n_02"] + x[,"n_32"])*log(-((-1 + r_rr)*r_rr)) + (x[,"n_14"] + x[,"n_20"])*log((-1 + r_rr)^2*r_rr) + (x[,"n_01"] + x[,"n_33"])*log(-((-2 + r_rr)*r_rr^2)) + 
    (x[,"n_10"] + x[,"n_24"])*log(-((-1 + r_rr)*r_rr^2)) + (x[,"n_00"] + x[,"n_34"])*log(r_rr^3) + (x[,"n_03"] + x[,"n_31"])*log((-1 + r_rr)^2*(1 + r_rr)) + (x[,"n_11"] + x[,"n_23"])*log(r_rr*(5 - 5*r_rr + 3*r_rr^2)) + 
    (x[,"n_12"] + x[,"n_22"])*log(9 - 5*r_rr + 5*r_rr^2) + (x[,"n_13"] + x[,"n_21"])*log(3 - 4*r_rr + 4*r_rr^2 - 3*r_rr^3)
  
  return(list(r_mat=cbind(r_cc, r_cm, r_cr, r_rc, r_rm, r_rr),
              LOD_mat=cbind(LOD_cc, LOD_cm, LOD_cr, LOD_rc, LOD_rm, LOD_rr),
              logL_mat=cbind(logL_cc, logL_cm, logL_cr, logL_rc, logL_rm, logL_rr),
              phasing_strategy="MLL", 
              possible_phases=c("coupling coupling", 
                                "coupling mixed",
                                "coupling repulsion",
                                "repulsion coupling",
                                "repulsion mixed",
                                "repulsion repulsion")))
}

#' @rdname r4_functions
#' @noRd
r4_1.3_2.2<-function(x, ncores=1){
  ############################
  ## COUPLING COUPLING
  ############################
  logLcc <- function(r,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n10 - n12 - 2*n14 - n20 - n22 - n24 - 2*n30 - n32 - 2*n34)*log(2) + 
      (-n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n11 + n33)*log(1 - r) + (n10 + n34)*log((-1 + r)^2) + (n13 + n31)*log(r) + 
      (n20 + n24)*log(-((-1 + r)*r)) + (n14 + n30)*log(r^2) + (n12 + n32)*log(2 + r - r^2) + n22*log(5 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLcc <- function(n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLcc,c(0,0.5),n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcc <- function(r){
  #     L <- 2^(-2*x[,"n_10"] - x[,"n_12"] - 2*x[,"n_14"] - x[,"n_20"] - x[,"n_22"] - x[,"n_24"] - 2*x[,"n_30"] - x[,"n_32"] - 2*x[,"n_34"])*9^(-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*(1 - r)^(x[,"n_11"] + x[,"n_33"])*((-1 + r)^2)^(x[,"n_10"] + x[,"n_34"])*r^(x[,"n_13"] + x[,"n_31"])*(-((-1 + r)*r))^(x[,"n_20"] + x[,"n_24"])*(r^2)^(x[,"n_14"] + x[,"n_30"])*(2 + r - r^2)^(x[,"n_12"] + x[,"n_32"])*(5 - 2*r + 2*r^2)^x[,"n_22"]
  #     return(L)}
  
  ############################
  ## COUPLING MIXED
  ############################
  logLcm <- function(r,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n10 - n11 - 2*n12 - n13 - 2*n14 - 2*n20 - 2*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34)*log(2) + 
      (-n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n10 + n14 + n30 + n34)*log(-((-1 + r)*r)) + n22*log(2 + r - r^2) + (n20 + n24)*log(1 - 2*r + 2*r^2) + 
      (n12 + n32)*log(5 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLcm <- function(n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLcm,c(0,0.5),n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcm <- function(r){
  #     L <- 2^(-2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*9^(-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*(-((-1 + r)*r))^(x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"])*(2 + r - r^2)^x[,"n_22"]*(1 - 2*r + 2*r^2)^(x[,"n_20"] + x[,"n_24"])*(5 - 2*r + 2*r^2)^(x[,"n_12"] + x[,"n_32"])
  #     return(L)}
  
  ############################
  ## REPULSION COUPLING
  ############################
  logLrc <- function(r,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n10 - n11 - 2*n12 - n13 - 2*n14 - 2*n20 - 2*n24 - 2*n30 - n31 - 2*n32 - n33 - 2*n34)*log(2) + 
      (-n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n10 + n14 + n30 + n34)*log(-((-1 + r)*r)) + n22*log(2 + r - r^2) + (n20 + n24)*log(1 - 2*r + 2*r^2) + 
      (n12 + n32)*log(5 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLrc <- function(n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLrc,c(0,0.5),n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lrc <- function(r){
  #     L <- 2^(-2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*9^(-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*(-((-1 + r)*r))^(x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"])*(2 + r - r^2)^x[,"n_22"]*(1 - 2*r + 2*r^2)^(x[,"n_20"] + x[,"n_24"])*(5 - 2*r + 2*r^2)^(x[,"n_12"] + x[,"n_32"])
  #     return(L)}
  
  ############################
  ## REPULSION MIXED
  ############################
  logLrm <- function(r,n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){ 
    L <- (-2*n10 - n12 - 2*n14 - n20 - n22 - n24 - 2*n30 - n32 - 2*n34)*log(2) + 
      (-n10 - n11 - n12 - n13 - n14 - n20 - n21 - n22 - n23 - n24 - n30 - n31 - n32 - n33 - n34)*log(9) + 
      (n13 + n31)*log(1 - r) + (n14 + n30)*log((-1 + r)^2) + (n11 + n33)*log(r) + 
      (n20 + n24)*log(-((-1 + r)*r)) + (n10 + n34)*log(r^2) + (n12 + n32)*log(2 + r - r^2) + 
      n22*log(5 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLrm <- function(n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34){
    optimize(logLrm,c(0,0.5),n10,n11,n12,n13,n14,n20,n21,n22,n23,n24,n30,n31,n32,n33,n34,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lrm <- function(r){
  #     L <- 2^(-2*x[,"n_10"] - x[,"n_12"] - 2*x[,"n_14"] - x[,"n_20"] - x[,"n_22"] - x[,"n_24"] - 2*x[,"n_30"] - x[,"n_32"] - 2*x[,"n_34"])*9^(-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*(1 - r)^(x[,"n_13"] + x[,"n_31"])*((-1 + r)^2)^(x[,"n_14"] + x[,"n_30"])*r^(x[,"n_11"] + x[,"n_33"])*(-((-1 + r)*r))^(x[,"n_20"] + x[,"n_24"])*(r^2)^(x[,"n_10"] + x[,"n_34"])*(2 + r - r^2)^(x[,"n_12"] + x[,"n_32"])*(5 - 2*r + 2*r^2)^x[,"n_22"]
  #     return(L)} 
  
  r_cc <- parallel::mcmapply(inter_logLcc,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                             mc.cores = ncores)
  
  r_cm <- parallel::mcmapply(inter_logLcm,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                             mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                              mc.cores = ncores)
  
  r_rm <-  parallel::mcmapply(inter_logLrm,x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_14"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_24"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],x[,"n_34"],
                              mc.cores = ncores)
  
  LOD_cc <- log10(((1 - r_cc)^(x[,"n_11"] + x[,"n_33"])*((-1 + r_cc)^2)^(x[,"n_10"] + x[,"n_34"])*r_cc^(x[,"n_13"] + x[,"n_31"])*(-((-1 + r_cc)*r_cc))^(x[,"n_20"] + x[,"n_24"])*(r_cc^2)^(x[,"n_14"] + x[,"n_30"])*(2 + r_cc - r_cc^2)^(x[,"n_12"] + x[,"n_32"])*(5 - 2*r_cc + 2*r_cc^2)^x[,"n_22"])/((1 - 0.5)^(x[,"n_11"] + x[,"n_33"])*((-1 + 0.5)^2)^(x[,"n_10"] + x[,"n_34"])*0.5^(x[,"n_13"] + x[,"n_31"])*(-((-1 + 0.5)*0.5))^(x[,"n_20"] + x[,"n_24"])*(0.5^2)^(x[,"n_14"] + x[,"n_30"])*(2 + 0.5 - 0.5^2)^(x[,"n_12"] + x[,"n_32"])*(5 - 2*0.5 + 2*0.5^2)^x[,"n_22"]))
  LOD_cm <- log10(((-((-1 + r_cm)*r_cm))^(x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"])*(2 + r_cm - r_cm^2)^x[,"n_22"]*(1 - 2*r_cm + 2*r_cm^2)^(x[,"n_20"] + x[,"n_24"])*(5 - 2*r_cm + 2*r_cm^2)^(x[,"n_12"] + x[,"n_32"]))/((-((-1 + 0.5)*0.5))^(x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"])*(2 + 0.5 - 0.5^2)^x[,"n_22"]*(1 - 2*0.5 + 2*0.5^2)^(x[,"n_20"] + x[,"n_24"])*(5 - 2*0.5 + 2*0.5^2)^(x[,"n_12"] + x[,"n_32"])))
  LOD_rc <- log10(((-((-1 + r_rc)*r_rc))^(x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"])*(2 + r_rc - r_rc^2)^x[,"n_22"]*(1 - 2*r_rc + 2*r_rc^2)^(x[,"n_20"] + x[,"n_24"])*(5 - 2*r_rc + 2*r_rc^2)^(x[,"n_12"] + x[,"n_32"]))/((-((-1 + 0.5)*0.5))^(x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"])*(2 + 0.5 - 0.5^2)^x[,"n_22"]*(1 - 2*0.5 + 2*0.5^2)^(x[,"n_20"] + x[,"n_24"])*(5 - 2*0.5 + 2*0.5^2)^(x[,"n_12"] + x[,"n_32"])))
  LOD_rm <- log10(((1 - r_rm)^(x[,"n_13"] + x[,"n_31"])*((-1 + r_rm)^2)^(x[,"n_14"] + x[,"n_30"])*r_rm^(x[,"n_11"] + x[,"n_33"])*(-((-1 + r_rm)*r_rm))^(x[,"n_20"] + x[,"n_24"])*(r_rm^2)^(x[,"n_10"] + x[,"n_34"])*(2 + r_rm - r_rm^2)^(x[,"n_12"] + x[,"n_32"])*(5 - 2*r_rm + 2*r_rm^2)^x[,"n_22"])/((1 - 0.5)^(x[,"n_13"] + x[,"n_31"])*((-1 + 0.5)^2)^(x[,"n_14"] + x[,"n_30"])*0.5^(x[,"n_11"] + x[,"n_33"])*(-((-1 + 0.5)*0.5))^(x[,"n_20"] + x[,"n_24"])*(0.5^2)^(x[,"n_10"] + x[,"n_34"])*(2 + 0.5 - 0.5^2)^(x[,"n_12"] + x[,"n_32"])*(5 - 2*0.5 + 2*0.5^2)^x[,"n_22"]))
  
  ## Record the logL:
  logL_cc <- (-2*x[,"n_10"] - x[,"n_12"] - 2*x[,"n_14"] - x[,"n_20"] - x[,"n_22"] - x[,"n_24"] - 2*x[,"n_30"] - x[,"n_32"] - 2*x[,"n_34"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_11"] + x[,"n_33"])*log(1 - r_cc) + (x[,"n_10"] + x[,"n_34"])*log((-1 + r_cc)^2) + (x[,"n_13"] + x[,"n_31"])*log(r_cc) + 
    (x[,"n_20"] + x[,"n_24"])*log(-((-1 + r_cc)*r_cc)) + (x[,"n_14"] + x[,"n_30"])*log(r_cc^2) + (x[,"n_12"] + x[,"n_32"])*log(2 + r_cc - r_cc^2) + x[,"n_22"]*log(5 - 2*r_cc + 2*r_cc^2)
  
  logL_cm <- (-2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"])*log(-((-1 + r_cm)*r_cm)) + x[,"n_22"]*log(2 + r_cm - r_cm^2) + (x[,"n_20"] + x[,"n_24"])*log(1 - 2*r_cm + 2*r_cm^2) + 
    (x[,"n_12"] + x[,"n_32"])*log(5 - 2*r_cm + 2*r_cm^2)
  
  logL_rc <- (-2*x[,"n_10"] - x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_14"] - 2*x[,"n_20"] - 2*x[,"n_24"] - 2*x[,"n_30"] - x[,"n_31"] - 2*x[,"n_32"] - x[,"n_33"] - 2*x[,"n_34"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_10"] + x[,"n_14"] + x[,"n_30"] + x[,"n_34"])*log(-((-1 + r_rc)*r_rc)) + x[,"n_22"]*log(2 + r_rc - r_rc^2) + (x[,"n_20"] + x[,"n_24"])*log(1 - 2*r_rc + 2*r_rc^2) + 
    (x[,"n_12"] + x[,"n_32"])*log(5 - 2*r_rc + 2*r_rc^2)
  
  logL_rm <- (-2*x[,"n_10"] - x[,"n_12"] - 2*x[,"n_14"] - x[,"n_20"] - x[,"n_22"] - x[,"n_24"] - 2*x[,"n_30"] - x[,"n_32"] - 2*x[,"n_34"])*log(2) + 
    (-x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_14"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_24"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"] - x[,"n_34"])*log(9) + 
    (x[,"n_13"] + x[,"n_31"])*log(1 - r_rm) + (x[,"n_14"] + x[,"n_30"])*log((-1 + r_rm)^2) + (x[,"n_11"] + x[,"n_33"])*log(r_rm) + 
    (x[,"n_20"] + x[,"n_24"])*log(-((-1 + r_rm)*r_rm)) + (x[,"n_10"] + x[,"n_34"])*log(r_rm^2) + (x[,"n_12"] + x[,"n_32"])*log(2 + r_rm - r_rm^2) + 
    x[,"n_22"]*log(5 - 2*r_rm + 2*r_rm^2)
  
  return(list(r_mat=cbind(r_cc, r_cm, r_rc, r_rm),
              LOD_mat=cbind(LOD_cc, LOD_cm, LOD_rc, LOD_rm),
              logL_mat=cbind(logL_cc, logL_cm, logL_rc, logL_rm),
              phasing_strategy="MLL", 
              possible_phases=c("coupling coupling", 
                                "coupling mixed",
                                "repulsion coupling",
                                "repulsion mixed")))
  
}

#' @rdname r4_functions
#' @noRd
r4_2.0_2.0<-function(x, ncores=1){
  #########################
  ## COUPLING
  #########################
  
  logLc<-function(r,n00,n01,n02,n10,n11,n12,n20,n21,n22){ #Coupling function
    L <- (-n00 - n02 + n11 - n20 - n22)*log(2) + 
      (-n00 - n01 - n02 - n10 - n11 - n12 - n20 - n21 - n22)*log(3) + 
      (n00 + n22)*log(((-1 + r)^2)) + (n01 + n10 + n12 + n21)*log((1 - r)*r) + 
      (n02 + n20)*log(r^2) + n11*log(1 - r + r*r)
    return(L)}
  
  inter_logLc <- function(n00,n01,n02,n10,n11,n12,n20,n21,n22){
    optimize(logLc,c(0,0.5),n00,n01,n02,n10,n11,n12,n20,n21,n22,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lc<-function(r){ #Coupling function
  #     L <- 2^(-x[,"n_00"] - x[,"n_02"] + x[,"n_11"] - x[,"n_20"] - x[,"n_22"])*3^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*((-1 + r)^2)^(x[,"n_00"] + x[,"n_22"])*(-((-1 + r)*r))^(x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*(r^2)^(x[,"n_02"] + x[,"n_20"])*(1 - r + r^2)^x[,"n_11"]
  #     return(L)}
  
  #########################
  ## MIXED
  #########################
  logLm<-function(r,n00,n01,n02,n10,n11,n12,n20,n21,n22){ 
    L <- (-2*n00 - n01 - 2*n02 - n10 - n12 - 2*n20 - n21 - 2*n22)*log(2) + 
      (-n00 - n01 - n02 - n10 - n11 - n12 - n20 - n21 - n22)*log(3) + 
      (n00 + n02 + n20 + n22)*log((1 - r)*r) + n11*log(1 + r - r*r) + 
      (n01 + n10 + n12 + n21)*log(1 - r + r*r)
    return(L)}
  
  inter_logLm <- function(n00,n01,n02,n10,n11,n12,n20,n21,n22){
    optimize(logLm,c(0,0.5),n00,n01,n02,n10,n11,n12,n20,n21,n22,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lm<-function(r){ #Mixed coupling and Repulsion 
  #     L <- 2^(-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*3^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*(-((-1 + r)*r))^(x[,"n_00"] + x[,"n_02"] + x[,"n_20"] + x[,"n_22"])*(1 + r - r^2)^x[,"n_11"]*(1 - r + r^2)^(x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])
  #     return(L)}
  
  #########################
  ## REPULSION
  #########################
  
  logLr<-function(r,n00,n01,n02,n10,n11,n12,n20,n21,n22){ #Repulsion function
    L <- (-n00 - n02 + n11 - n20 - n22)*log(2) + 
      (-n00 - n01 - n02 - n10 - n11 - n12 - n20 - n21 - n22)*log(3) + 
      (n02 + n20)*log((-1 + r)^2)+ (n01 + n10 + n12 + n21)*log((1 - r)*r) + 
      (n00 + n22)*log(r^2) + n11*log(1 - r + r*r)
    return(L)}
  
  inter_logLr <- function(n00,n01,n02,n10,n11,n12,n20,n21,n22){
    optimize(logLr,c(0,0.5),n00,n01,n02,n10,n11,n12,n20,n21,n22,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lr<-function(r){ #Repulsion function
  #     L <- 2^(-x[,"n_00"] - x[,"n_02"] + x[,"n_11"] - x[,"n_20"] - x[,"n_22"])*3^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*((-1 + r)^2)^(x[,"n_02"] + x[,"n_20"])*(-((-1 + r)*r))^(x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*(r^2)^(x[,"n_00"] + x[,"n_22"])*(1 - r + r^2)^x[,"n_11"]
  #     return(L)}
  
  r_c <- parallel::mcmapply(inter_logLc, x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],
                            mc.cores = ncores)
  
  r_m <-  parallel::mcmapply(inter_logLm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],
                             mc.cores = ncores)
  
  r_r <- parallel::mcmapply(inter_logLr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],
                            mc.cores = ncores)
  
  LOD_c <- (x[,"n_00"] + x[,"n_22"])*log10((-1 + r_c)^2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(-((-1 + r_c)*r_c)) + (x[,"n_02"] + x[,"n_20"])*log10(r_c^2) + x[,"n_11"]*log10(1 - r_c + r_c^2) - (x[,"n_00"] + x[,"n_22"])*log10((-1 + 0.5)^2) - (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(-((-1 + 0.5)*0.5)) - (x[,"n_02"] + x[,"n_20"])*log10(0.5^2) - x[,"n_11"]*log10(1 - 0.5 + 0.5^2)
  LOD_m <- log10((-((-1 + r_m)*r_m))^(x[,"n_00"] + x[,"n_02"] + x[,"n_20"] + x[,"n_22"])*(1 + r_m - r_m^2)^x[,"n_11"]*(1 - r_m + r_m^2)^(x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])) - log10((-((-1 + 0.5)*0.5))^(x[,"n_00"] + x[,"n_02"] + x[,"n_20"] + x[,"n_22"])*(1 + 0.5 - 0.5^2)^x[,"n_11"]*(1 - 0.5 + 0.5^2)^(x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"]))
  LOD_r <- (x[,"n_02"] + x[,"n_20"])*log10((-1 + r_r)^2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(-((-1 + r_r)*r_r)) + (x[,"n_00"] + x[,"n_22"])*log10(r_r^2) + x[,"n_11"]*log10(1 - r_r + r_r^2) - (x[,"n_02"] + x[,"n_20"])*log10((-1 + 0.5)^2) - (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(-((-1 + 0.5)*0.5)) -(x[,"n_00"] + x[,"n_22"])*log10(0.5^2) - x[,"n_11"]*log10(1 - 0.5 + 0.5^2)
  
  ## Having checked, phasing using max logL is more accurate than min r...
  logL_c <- (-x[,"n_00"] - x[,"n_02"] + x[,"n_11"] - x[,"n_20"] - x[,"n_22"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*log(3) + 
    (x[,"n_00"] + x[,"n_22"])*log(((-1 + r_c)^2)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log((1 - r_c)*r_c) + 
    (x[,"n_02"] + x[,"n_20"])*log(r_c^2) + x[,"n_11"]*log(1 - r_c + r_c*r_c)
  
  logL_m <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*log(3) + 
    (x[,"n_00"] + x[,"n_02"] + x[,"n_20"] + x[,"n_22"])*log((1 - r_m)*r_m) + x[,"n_11"]*log(1 + r_m - r_m*r_m) + 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log(1 - r_m + r_m*r_m)
  
  logL_r <- (-x[,"n_00"] - x[,"n_02"] + x[,"n_11"] - x[,"n_20"] - x[,"n_22"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*log(3) + 
    (x[,"n_02"] + x[,"n_20"])*log((-1 + r_r)^2)+ (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log((1 - r_r)*r_r) + 
    (x[,"n_00"] + x[,"n_22"])*log(r_r^2) + x[,"n_11"]*log(1 - r_r + r_r*r_r)
  
  return(list(r_mat=cbind(r_c, r_m, r_r),
              LOD_mat=cbind(LOD_c, LOD_m, LOD_r),
              logL_mat=cbind(logL_c, logL_m, logL_r),
              phasing_strategy="MLL", 
              possible_phases=c("coupling",
                                "mixed",
                                "repulsion")))
}

#' @rdname r4_functions
#' @noRd
r4_1.1_1.1<-function(x, ncores=1){
  #########################
  ## COUPLING COUPLING
  #########################
  
  logLcc<-function(r,n00,n01,n02,n10,n11,n12,n20,n21,n22){ 
    L <- (-2*n00 - n01 - 2*n02 - n10 - n12 - 2*n20 - n21 - 2*n22)*log(2) + (n00 + n22)*log((-1 + r)^2) + (n01 + n10 + n12 + n21)*log(-((-1 + r)*r)) + 
      (n02 + n20)*log(r^2) + n11*log(1/2 - r + r^2)
    return(L)}
  
  inter_logLcc <- function(n00,n01,n02,n10,n11,n12,n20,n21,n22){
    optimize(logLcc,c(0,0.5),n00,n01,n02,n10,n11,n12,n20,n21,n22,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcc<-function(r){
  #     L <- 2^(-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*((-1 + r)^2)^(x[,"n_00"] + x[,"n_22"])*(-((-1 + r)*r))^(x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*(r^2)^(x[,"n_02"] + x[,"n_20"])*(1/2 - r + r^2)^x[,"n_11"]
  #     return(L)}
  
  #########################
  ## COUPLING REPULSION or REPULSION COUPLING
  #########################
  
  logLcr<-function(r,n00,n01,n02,n10,n11,n12,n20,n21,n22){ 
    L <- (-2*n00 - n01 - 2*n02 - n10 - n11 - n12 - 2*n20 - n21 - 2*n22)*log(2) + (-n00 - n01 - n02 - n10 - n11 - n12 - n20 - n21 - n22)*log(3) + 
      (n02 + n20)*log(-((-2 + r)*r)) + (n00 + n22)*log(-((-1 + r)*(1 + r))) + n11*log(1 + 2*r - 2*r^2) + (n01 + n10 + n12 + n21)*log(1 - r + r^2)
    return(L)}
  
  inter_logLcr <- function(n00,n01,n02,n10,n11,n12,n20,n21,n22){
    optimize(logLcr,c(0,0.5),n00,n01,n02,n10,n11,n12,n20,n21,n22,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcr<-function(r){
  #     L <- 2^(-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*3^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*(-((-2 + r)*r))^(x[,"n_02"] + x[,"n_20"])*(-((-1 + r)*(1 + r)))^(x[,"n_00"] + x[,"n_22"])*(1 + 2*r - 2*r^2)^x[,"n_11"]*(1 - r + r^2)^(x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])
  #     return(L)}
  
  #########################
  ## REPULSION REPULSION
  #########################
  
  logLrr<-function(r,n00,n01,n02,n10,n11,n12,n20,n21,n22){ 
    L <- (-2*n00 - n01 - 2*n02 - n10 - n11 - n12 - 2*n20 - n21 - 2*n22)*log(2) + 
      (-n00 - n01 - n02 - n10 - n11 - n12 - n20 - n21 - n22)*log(9) + (n02 + n20)*log((-2 + r)^2) + (n01 + n10 + n12 + n21)*log(-((-2 + r)*(1 + r))) + 
      (n00 + n22)*log((1 + r)^2) + n11*log(5 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLrr <- function(n00,n01,n02,n10,n11,n12,n20,n21,n22){
    optimize(logLrr,c(0,0.5),n00,n01,n02,n10,n11,n12,n20,n21,n22,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lrr<-function(r){
  #     L <- 2^(-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*9^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*((-2 + r)^2)^(x[,"n_02"] + x[,"n_20"])*(-((-2 + r)*(1 + r)))^(x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*((1 + r)^2)^(x[,"n_00"] + x[,"n_22"])*(5 - 2*r + 2*r^2)^x[,"n_11"]
  #     return(L)}
  
  r_cc <- parallel::mcmapply(inter_logLcc, x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],
                             mc.cores = ncores)
  
  r_mix <-  parallel::mcmapply(inter_logLcr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],
                               mc.cores = ncores)
  
  r_rr <- parallel::mcmapply(inter_logLrr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],
                             mc.cores = ncores)
  
  LOD_cc <- (x[,"n_00"] + x[,"n_22"])*log10((-1 + r_cc)^2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(-((-1 + r_cc)*r_cc)) + (x[,"n_02"] + x[,"n_20"])*log10(r_cc^2) + 
    x[,"n_11"]*log10(1/2 - r_cc + r_cc^2) - (x[,"n_00"] + x[,"n_22"])*log10((-1 + 0.5)^2) - (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(-((-1 + 0.5)*0.5)) - 
    (x[,"n_02"] + x[,"n_20"])*log10(0.5^2) - x[,"n_11"]*log10(1/2 - 0.5 + 0.5^2)
  
  LOD_mix <- (x[,"n_02"] + x[,"n_20"])*log10(-((-2 + r_mix)*r_mix)) + (x[,"n_00"] + x[,"n_22"])*log10(-((-1 + r_mix)*(1 + r_mix))) + x[,"n_11"]*log10(1 + 2*r_mix - 2*r_mix^2) + 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(1 - r_mix + r_mix^2) -  (x[,"n_02"] + x[,"n_20"])*log10(-((-2 + 0.5)*0.5)) - (x[,"n_00"] + x[,"n_22"])*log10(-((-1 + 0.5)*(1 + 0.5))) -
    x[,"n_11"]*log10(1 + 2*0.5 - 2*0.5^2) - (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(1 - 0.5 + 0.5^2)
  
  LOD_rr <- (x[,"n_02"] + x[,"n_20"])*log10((-2 + r_rr)^2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(-((-2 + r_rr)*(1 + r_rr))) + (x[,"n_00"] + x[,"n_22"])*log10((1 + r_rr)^2) + 
    x[,"n_11"]*log10(5 - 2*r_rr + 2*r_rr^2) - (x[,"n_02"] + x[,"n_20"])*log10((-2 + 0.5)^2) - (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(-((-2 + 0.5)*(1 + 0.5))) - 
    (x[,"n_00"] + x[,"n_22"])*log10((1 + 0.5)^2) - x[,"n_11"]*log10(5 - 2*0.5 + 2*0.5^2)
  
  ## Record the logL values:
  logL_cc <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*log(2) + (x[,"n_00"] + x[,"n_22"])*log((-1 + r_cc)^2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log(-((-1 + r_cc)*r_cc)) + 
    (x[,"n_02"] + x[,"n_20"])*log(r_cc^2) + x[,"n_11"]*log(1/2 - r_cc + r_cc^2)
  
  logL_mix <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*log(2) + (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*log(3) + 
    (x[,"n_02"] + x[,"n_20"])*log(-((-2 + r_mix)*r_mix)) + (x[,"n_00"] + x[,"n_22"])*log(-((-1 + r_mix)*(1 + r_mix))) + x[,"n_11"]*log(1 + 2*r_mix - 2*r_mix^2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log(1 - r_mix + r_mix^2)
  
  logL_rr <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*log(2) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*log(9) + (x[,"n_02"] + x[,"n_20"])*log((-2 + r_rr)^2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log(-((-2 + r_rr)*(1 + r_rr))) + 
    (x[,"n_00"] + x[,"n_22"])*log((1 + r_rr)^2) + x[,"n_11"]*log(5 - 2*r_rr + 2*r_rr^2)
  
  return(list(r_mat=cbind(r_cc, r_mix, r_rr),
              LOD_mat=cbind(LOD_cc, LOD_mix, LOD_rr),
              logL_mat=cbind(logL_cc, logL_mix, logL_rr),
              phasing_strategy="MLL", 
              possible_phases=c("coupling coupling", # should be "coupling"?
                                "mixed",
                                "repulsion repulsion"))) # should be "repulsion"?
  
}

#' @rdname r4_functions
#' @noRd
r4_1.2_1.2<-function(x, ncores=1){
  #########################
  ## COUPLING COUPLING
  #########################
  
  logLcc<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n03 - n11 - n12 - n21 - n22 - n30 - n33)*log(3) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(4) + 
      (n00 + n33)*log(-(-1 + r)^3) + (n01 + n10 + n23 + n32)*log((-1 + r)^2*r) + (n02 + n13 + n20 + n31)*log(-((-1 + r)*r^2)) + 
      (n03 + n30)*log(r^3) + (n12 + n21)*log(r*(8 - 12*r + 9*r^2)) + (n11 + n22)*log(5 - 11*r + 15*r^2 - 9*r^3)
    return(L)}
  
  inter_logLcc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcc<-function(r){ 
  #     L <- 3^(-x[,"n_00"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_21"] - x[,"n_22"] - x[,"n_30"] - x[,"n_33"])*4^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*(-(-1 + r)^3)^(x[,"n_00"] + x[,"n_33"])*((-1 + r)^2*r)^(x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(-((-1 + r)*r^2))^(x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(r^3)^(x[,"n_03"] + x[,"n_30"])*(r*(8 - 12*r + 9*r^2))^(x[,"n_12"] + x[,"n_21"])*(5 - 11*r + 15*r^2 - 9*r^3)^(x[,"n_11"] + x[,"n_22"])
  #     return(L)}
  
  #########################
  ## COUPLING MIXED
  #########################
  logLcm<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(24) + 
      (n00 + n33)*log((-1 + r)^2*r) + (n03 + n30)*log(-((-1 + r)*r^2)) + (n02 + n13 + n20 + n31)*log(r*(3 - 4*r + 3*r^2)) + 
      (n12 + n21)*log(4 - 4*r + 13*r^2 - 9*r^3) + (n01 + n10 + n23 + n32)*log(2 - 4*r + 5*r^2 - 3*r^3) + 
      (n11 + n22)*log(4 + 5*r - 14*r^2 + 9*r^3)
    return(L)}
  
  inter_logLcm <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcm,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #########################
  ## COUPLING REPULSION
  #########################
  logLcr<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(12) + 
      (n03 + n30)*log((-1 + r)^2*r) + (n00 + n33)*log(-((-1 + r)*r^2)) + (n01 + n10 + n23 + n32)*log(r*(2 - 4*r + 3*r^2)) + 
      (n12 + n21)*log(r*(9 - 14*r + 9*r^2)) + (n11 + n22)*log(4 - 8*r + 13*r^2 - 9*r^3) + 
      (n02 + n13 + n20 + n31)*log(1 - 3*r + 5*r^2 - 3*r^3)
    return(L)}
  
  inter_logLcr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #########################
  ## REPULSION COUPLING
  #########################
  logLrc<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(36) + 
      (n03 + n30)*log(-((-2 + r)*r^2)) + (n00 + n33)*log((-1 + r)^2*(1 + r)) + (n02 + n13 + n20 + n31)*log(r*(4 - 5*r + 3*r^2)) + 
      (n12 + n21)*log(8 - 8*r + 14*r^2 - 9*r^3) + (n01 + n10 + n23 + n32)*log(2 - 3*r + 4*r^2 - 3*r^3) + 
      (n11 + n22)*log(5 + 7*r - 13*r^2 + 9*r^3)
    return(L)}
  
  inter_logLrc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #########################
  ## REPULSION MIXED
  #########################
  logLrm<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - 2*n10 - n11 - n12 - 2*n13 - 2*n20 - n21 - n22 - 2*n23 - 2*n30 - 2*n31 - 2*n32 - 2*n33)*log(3) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(8) + 
      (n03 + n30)*log((-2 + r)*(-1 + r)*r) + (n00 + n33)*log(-((-1 + r)*r*(1 + r))) + (n11 + n22)*log(4 - r + 4*r^2 - 3*r^3) + 
      (n02 + n13 + n20 + n31)*log(4 - 5*r + 6*r^2 - 3*r^3) + (n12 + n21)*log(4 + 2*r - 5*r^2 + 3*r^3) + 
      (n01 + n10 + n23 + n32)*log(2 + 2*r - 3*r^2 + 3*r^3)
    return(L)}
  
  inter_logLrm <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrm,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #########################
  ## REPULSION REPULSION
  #########################
  logLrr<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(36) + 
      (n03 + n30)*log(-((-2 + r)*(-1 + r)^2)) + (n00 + n33)*log(r^2*(1 + r)) + (n01 + n10 + n23 + n32)*log(r*(2 + 2*r - 3*r^2)) + 
      (n12 + n21)*log(10 - 13*r + 16*r^2 - 9*r^3) + (n02 + n13 + n20 + n31)*log(1 + 3*r - 7*r^2 + 3*r^3) + 
      (n11 + n22)*log(4 + 8*r - 11*r^2 + 9*r^3)
    return(L)}
  
  inter_logLrr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_cc <- parallel::mcmapply(inter_logLcc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_cm <- parallel::mcmapply(inter_logLcm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_cr <-  parallel::mcmapply(inter_logLcr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  r_rm <-  parallel::mcmapply(inter_logLrm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  r_rr <-  parallel::mcmapply(inter_logLrr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  LOD_cc <- (x[,"n_00"] + x[,"n_33"])*log10(-(-1 + r_cc)^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10((-1 + r_cc)^2*r_cc) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(-((-1 + r_cc)*r_cc^2)) + 
    (x[,"n_03"] + x[,"n_30"])*log10(r_cc^3) + (x[,"n_12"] + x[,"n_21"])*log10(r_cc*(8 - 12*r_cc + 9*r_cc^2)) + (x[,"n_11"] + x[,"n_22"])*log10(5 - 11*r_cc + 15*r_cc^2 - 9*r_cc^3) - 
    (x[,"n_00"] + x[,"n_33"])*log10(-(-1 + 0.5)^3) - (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10((-1 + 0.5)^2*0.5) - (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(-((-1 + 0.5)*0.5^2)) - 
    (x[,"n_03"] + x[,"n_30"])*log10(0.5^3) - (x[,"n_12"] + x[,"n_21"])*log10(0.5*(8 - 12*0.5 + 9*0.5^2)) - (x[,"n_11"] + x[,"n_22"])*log10(5 - 11*0.5 + 15*0.5^2 - 9*0.5^3)
  
  LOD_cm <- (x[,"n_00"] + x[,"n_33"])*log10((-1 + r_cm)^2*r_cm) + (x[,"n_03"] + x[,"n_30"])*log10(-((-1 + r_cm)*r_cm^2)) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(r_cm*(3 - 4*r_cm + 3*r_cm^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log10(4 - 4*r_cm + 13*r_cm^2 - 9*r_cm^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 - 4*r_cm + 5*r_cm^2 - 3*r_cm^3) + 
    (x[,"n_11"] + x[,"n_22"])*log10(4 + 5*r_cm - 14*r_cm^2 + 9*r_cm^3) - (x[,"n_00"] + x[,"n_33"])*log10((-1 + 0.5)^2*0.5) - (x[,"n_03"] + x[,"n_30"])*log10(-((-1 + 0.5)*0.5^2)) - 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(0.5*(3 - 4*0.5 + 3*0.5^2)) - (x[,"n_12"] + x[,"n_21"])*log10(4 - 4*0.5 + 13*0.5^2 - 9*0.5^3) - 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 - 4*0.5 + 5*0.5^2 - 3*0.5^3) - (x[,"n_11"] + x[,"n_22"])*log10(4 + 5*0.5 - 14*0.5^2 + 9*0.5^3)
  
  LOD_cr <- (x[,"n_03"] + x[,"n_30"])*log10((-1 + r_cr)^2*r_cr) + (x[,"n_00"] + x[,"n_33"])*log10(-((-1 + r_cr)*r_cr^2)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(r_cr*(2 - 4*r_cr + 3*r_cr^2)) +
    (x[,"n_12"] + x[,"n_21"])*log10(r_cr*(9 - 14*r_cr + 9*r_cr^2)) + (x[,"n_11"] + x[,"n_22"])*log10(4 - 8*r_cr + 13*r_cr^2 - 9*r_cr^3) + 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(1 - 3*r_cr + 5*r_cr^2 - 3*r_cr^3) - (x[,"n_03"] + x[,"n_30"])*log10((-1 + 0.5)^2*0.5) - (x[,"n_00"] + x[,"n_33"])*log10(-((-1 + 0.5)*0.5^2)) - 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(0.5*(2 - 4*0.5 + 3*0.5^2)) - (x[,"n_12"] + x[,"n_21"])*log10(0.5*(9 - 14*0.5 + 9*0.5^2)) - 
    (x[,"n_11"] + x[,"n_22"])*log10(4 - 8*0.5 + 13*0.5^2 - 9*0.5^3) - (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(1 - 3*0.5 + 5*0.5^2 - 3*0.5^3)
  
  LOD_rc <- (x[,"n_03"] + x[,"n_30"])*log10(-((-2 + r_rc)*r_rc^2)) + (x[,"n_00"] + x[,"n_33"])*log10((-1 + r_rc)^2*(1 + r_rc)) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(r_rc*(4 - 5*r_rc + 3*r_rc^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log10(8 - 8*r_rc + 14*r_rc^2 - 9*r_rc^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 - 3*r_rc + 4*r_rc^2 - 3*r_rc^3) + 
    (x[,"n_11"] + x[,"n_22"])*log10(5 + 7*r_rc - 13*r_rc^2 + 9*r_rc^3) - (x[,"n_03"] + x[,"n_30"])*log10(-((-2 + 0.5)*0.5^2)) - (x[,"n_00"] + x[,"n_33"])*log10((-1 + 0.5)^2*(1 + 0.5)) -
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(0.5*(4 - 5*0.5 + 3*0.5^2)) - (x[,"n_12"] + x[,"n_21"])*log10(8 - 8*0.5 + 14*0.5^2 - 9*0.5^3) - 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 - 3*0.5 + 4*0.5^2 - 3*0.5^3) - (x[,"n_11"] + x[,"n_22"])*log10(5 + 7*0.5 - 13*0.5^2 + 9*0.5^3)
  
  LOD_rm <- (x[,"n_03"] + x[,"n_30"])*log10((-2 + r_rm)*(-1 + r_rm)*r_rm) + (x[,"n_00"] + x[,"n_33"])*log10(-((-1 + r_rm)*r_rm*(1 + r_rm))) + (x[,"n_11"] + x[,"n_22"])*log10(4 - r_rm + 4*r_rm^2 - 3*r_rm^3) +
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(4 - 5*r_rm + 6*r_rm^2 - 3*r_rm^3) + (x[,"n_12"] + x[,"n_21"])*log10(4 + 2*r_rm - 5*r_rm^2 + 3*r_rm^3) + 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 + 2*r_rm - 3*r_rm^2 + 3*r_rm^3) - (x[,"n_03"] + x[,"n_30"])*log10((-2 + 0.5)*(-1 + 0.5)*0.5) - 
    (x[,"n_00"] + x[,"n_33"])*log10(-((-1 + 0.5)*0.5*(1 + 0.5))) - (x[,"n_11"] + x[,"n_22"])*log10(4 - 0.5 + 4*0.5^2 - 3*0.5^3) - 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(4 - 5*0.5 + 6*0.5^2 - 3*0.5^3) - (x[,"n_12"] + x[,"n_21"])*log10(4 + 2*0.5 - 5*0.5^2 + 3*0.5^3) - 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 + 2*0.5 - 3*0.5^2 + 3*0.5^3)
  
  LOD_rr <- (x[,"n_03"] + x[,"n_30"])*log10(-((-2 + r_rr)*(-1 + r_rr)^2)) + (x[,"n_00"] + x[,"n_33"])*log10(r_rr^2*(1 + r_rr)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(r_rr*(2 + 2*r_rr - 3*r_rr^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log10(10 - 13*r_rr + 16*r_rr^2 - 9*r_rr^3) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(1 + 3*r_rr - 7*r_rr^2 + 3*r_rr^3) + 
    (x[,"n_11"] + x[,"n_22"])*log10(4 + 8*r_rr - 11*r_rr^2 + 9*r_rr^3) - (x[,"n_03"] + x[,"n_30"])*log10(-((-2 + 0.5)*(-1 + 0.5)^2)) - (x[,"n_00"] + x[,"n_33"])*log10(0.5^2*(1 + 0.5)) -
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(0.5*(2 + 2*0.5 - 3*0.5^2)) - (x[,"n_12"] + x[,"n_21"])*log10(10 - 13*0.5 + 16*0.5^2 - 9*0.5^3) - 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(1 + 3*0.5 - 7*0.5^2 + 3*0.5^3) - (x[,"n_11"] + x[,"n_22"])*log10(4 + 8*0.5 - 11*0.5^2 + 9*0.5^3)
  
  ## Record the logL also:
  logL_cc <- (-x[,"n_00"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_21"] - x[,"n_22"] - x[,"n_30"] - x[,"n_33"])*log(3) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(4) + 
    (x[,"n_00"] + x[,"n_33"])*log(-(-1 + r_cc)^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log((-1 + r_cc)^2*r_cc) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(-((-1 + r_cc)*r_cc^2)) + 
    (x[,"n_03"] + x[,"n_30"])*log(r_cc^3) + (x[,"n_12"] + x[,"n_21"])*log(r_cc*(8 - 12*r_cc + 9*r_cc^2)) + (x[,"n_11"] + x[,"n_22"])*log(5 - 11*r_cc + 15*r_cc^2 - 9*r_cc^3)
  
  logL_cm <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(24) + 
    (x[,"n_00"] + x[,"n_33"])*log((-1 + r_cm)^2*r_cm) + (x[,"n_03"] + x[,"n_30"])*log(-((-1 + r_cm)*r_cm^2)) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(r_cm*(3 - 4*r_cm + 3*r_cm^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log(4 - 4*r_cm + 13*r_cm^2 - 9*r_cm^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(2 - 4*r_cm + 5*r_cm^2 - 3*r_cm^3) + 
    (x[,"n_11"] + x[,"n_22"])*log(4 + 5*r_cm - 14*r_cm^2 + 9*r_cm^3)
  
  logL_cr <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(12) + 
    (x[,"n_03"] + x[,"n_30"])*log((-1 + r_cr)^2*r_cr) + (x[,"n_00"] + x[,"n_33"])*log(-((-1 + r_cr)*r_cr^2)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(r_cr*(2 - 4*r_cr + 3*r_cr^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log(r_cr*(9 - 14*r_cr + 9*r_cr^2)) + (x[,"n_11"] + x[,"n_22"])*log(4 - 8*r_cr + 13*r_cr^2 - 9*r_cr^3) + 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(1 - 3*r_cr + 5*r_cr^2 - 3*r_cr^3)
  
  logL_rc <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(36) + 
    (x[,"n_03"] + x[,"n_30"])*log(-((-2 + r_rc)*r_rc^2)) + (x[,"n_00"] + x[,"n_33"])*log((-1 + r_rc)^2*(1 + r_rc)) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(r_rc*(4 - 5*r_rc + 3*r_rc^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log(8 - 8*r_rc + 14*r_rc^2 - 9*r_rc^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(2 - 3*r_rc + 4*r_rc^2 - 3*r_rc^3) + 
    (x[,"n_11"] + x[,"n_22"])*log(5 + 7*r_rc - 13*r_rc^2 + 9*r_rc^3)
  
  logL_rm <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - 2*x[,"n_23"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"])*log(3) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(8) + 
    (x[,"n_03"] + x[,"n_30"])*log((-2 + r_rm)*(-1 + r_rm)*r_rm) + (x[,"n_00"] + x[,"n_33"])*log(-((-1 + r_rm)*r_rm*(1 + r_rm))) + (x[,"n_11"] + x[,"n_22"])*log(4 - r_rm + 4*r_rm^2 - 3*r_rm^3) + 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(4 - 5*r_rm + 6*r_rm^2 - 3*r_rm^3) + (x[,"n_12"] + x[,"n_21"])*log(4 + 2*r_rm - 5*r_rm^2 + 3*r_rm^3) + 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(2 + 2*r_rm - 3*r_rm^2 + 3*r_rm^3)
  
  logL_rr <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(36) + 
    (x[,"n_03"] + x[,"n_30"])*log(-((-2 + r_rr)*(-1 + r_rr)^2)) + (x[,"n_00"] + x[,"n_33"])*log(r_rr^2*(1 + r_rr)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(r_rr*(2 + 2*r_rr - 3*r_rr^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log(10 - 13*r_rr + 16*r_rr^2 - 9*r_rr^3) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(1 + 3*r_rr - 7*r_rr^2 + 3*r_rr^3) + 
    (x[,"n_11"] + x[,"n_22"])*log(4 + 8*r_rr - 11*r_rr^2 + 9*r_rr^3)
  
  return(list(r_mat=cbind(r_cc, r_cm, r_cr, r_rc, r_rm, r_rr),
              LOD_mat=cbind(LOD_cc, LOD_cm, LOD_cr, LOD_rc, LOD_rm, LOD_rr),
              logL_mat=cbind(logL_cc, logL_cm, logL_cr, logL_rc, logL_rm, logL_rr),
              phasing_strategy="MLL", 
              possible_phases=c("coupling coupling", 
                                "coupling mixed",
                                "coupling repulsion",
                                "repulsion coupling",
                                "repulsion mixed",
                                "repulsion repulsion")))
  
}

#' @rdname r4_functions
#' @noRd
r4_2.1_2.1<-function(x, ncores=1){ # copy of r_1.2_1.2.. Only with different phases
  #########################
  ## COUPLING COUPLING
  #########################
  
  logLcc<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n03 - n11 - n12 - n21 - n22 - n30 - n33)*log(3) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(4) + 
      (n00 + n33)*log(-(-1 + r)^3) + (n01 + n10 + n23 + n32)*log((-1 + r)^2*r) + (n02 + n13 + n20 + n31)*log(-((-1 + r)*r^2)) + 
      (n03 + n30)*log(r^3) + (n12 + n21)*log(r*(8 - 12*r + 9*r^2)) + (n11 + n22)*log(5 - 11*r + 15*r^2 - 9*r^3)
    return(L)}
  
  inter_logLcc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcc<-function(r){ 
  #     L <- 3^(-x[,"n_00"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_21"] - x[,"n_22"] - x[,"n_30"] - x[,"n_33"])*4^(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*(-(-1 + r)^3)^(x[,"n_00"] + x[,"n_33"])*((-1 + r)^2*r)^(x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(-((-1 + r)*r^2))^(x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(r^3)^(x[,"n_03"] + x[,"n_30"])*(r*(8 - 12*r + 9*r^2))^(x[,"n_12"] + x[,"n_21"])*(5 - 11*r + 15*r^2 - 9*r^3)^(x[,"n_11"] + x[,"n_22"])
  #     return(L)}
  
  #########################
  ## COUPLING MIXED
  #########################
  logLcm<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(24) + 
      (n00 + n33)*log((-1 + r)^2*r) + (n03 + n30)*log(-((-1 + r)*r^2)) + (n02 + n13 + n20 + n31)*log(r*(3 - 4*r + 3*r^2)) + 
      (n12 + n21)*log(4 - 4*r + 13*r^2 - 9*r^3) + (n01 + n10 + n23 + n32)*log(2 - 4*r + 5*r^2 - 3*r^3) + 
      (n11 + n22)*log(4 + 5*r - 14*r^2 + 9*r^3)
    return(L)}
  
  inter_logLcm <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcm,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #########################
  ## COUPLING REPULSION
  #########################
  logLcr<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(12) + 
      (n03 + n30)*log((-1 + r)^2*r) + (n00 + n33)*log(-((-1 + r)*r^2)) + (n01 + n10 + n23 + n32)*log(r*(2 - 4*r + 3*r^2)) + 
      (n12 + n21)*log(r*(9 - 14*r + 9*r^2)) + (n11 + n22)*log(4 - 8*r + 13*r^2 - 9*r^3) + 
      (n02 + n13 + n20 + n31)*log(1 - 3*r + 5*r^2 - 3*r^3)
    return(L)}
  
  inter_logLcr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLcr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #########################
  ## REPULSION COUPLING
  #########################
  logLrc<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(36) + 
      (n03 + n30)*log(-((-2 + r)*r^2)) + (n00 + n33)*log((-1 + r)^2*(1 + r)) + (n02 + n13 + n20 + n31)*log(r*(4 - 5*r + 3*r^2)) + 
      (n12 + n21)*log(8 - 8*r + 14*r^2 - 9*r^3) + (n01 + n10 + n23 + n32)*log(2 - 3*r + 4*r^2 - 3*r^3) + 
      (n11 + n22)*log(5 + 7*r - 13*r^2 + 9*r^3)
    return(L)}
  
  inter_logLrc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrc,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #########################
  ## REPULSION MIXED
  #########################
  logLrm<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - 2*n10 - n11 - n12 - 2*n13 - 2*n20 - n21 - n22 - 2*n23 - 2*n30 - 2*n31 - 2*n32 - 2*n33)*log(3) + 
      (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(8) + 
      (n03 + n30)*log((-2 + r)*(-1 + r)*r) + (n00 + n33)*log(-((-1 + r)*r*(1 + r))) + (n11 + n22)*log(4 - r + 4*r^2 - 3*r^3) + 
      (n02 + n13 + n20 + n31)*log(4 - 5*r + 6*r^2 - 3*r^3) + (n12 + n21)*log(4 + 2*r - 5*r^2 + 3*r^3) + 
      (n01 + n10 + n23 + n32)*log(2 + 2*r - 3*r^2 + 3*r^3)
    return(L)}
  
  inter_logLrm <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrm,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #########################
  ## REPULSION REPULSION
  #########################
  logLrr<-function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){ 
    L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(36) + 
      (n03 + n30)*log(-((-2 + r)*(-1 + r)^2)) + (n00 + n33)*log(r^2*(1 + r)) + (n01 + n10 + n23 + n32)*log(r*(2 + 2*r - 3*r^2)) + 
      (n12 + n21)*log(10 - 13*r + 16*r^2 - 9*r^3) + (n02 + n13 + n20 + n31)*log(1 + 3*r - 7*r^2 + 3*r^3) + 
      (n11 + n22)*log(4 + 8*r - 11*r^2 + 9*r^3)
    return(L)}
  
  inter_logLrr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33){
    optimize(logLrr,c(0,0.5),n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  r_cc <- parallel::mcmapply(inter_logLcc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_cm <- parallel::mcmapply(inter_logLcm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_cr <-  parallel::mcmapply(inter_logLcr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  r_rm <-  parallel::mcmapply(inter_logLrm,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  r_rr <-  parallel::mcmapply(inter_logLrr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  LOD_cc <- (x[,"n_00"] + x[,"n_33"])*log10(-(-1 + r_cc)^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10((-1 + r_cc)^2*r_cc) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(-((-1 + r_cc)*r_cc^2)) + 
    (x[,"n_03"] + x[,"n_30"])*log10(r_cc^3) + (x[,"n_12"] + x[,"n_21"])*log10(r_cc*(8 - 12*r_cc + 9*r_cc^2)) + (x[,"n_11"] + x[,"n_22"])*log10(5 - 11*r_cc + 15*r_cc^2 - 9*r_cc^3) - 
    (x[,"n_00"] + x[,"n_33"])*log10(-(-1 + 0.5)^3) - (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10((-1 + 0.5)^2*0.5) - (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(-((-1 + 0.5)*0.5^2)) - 
    (x[,"n_03"] + x[,"n_30"])*log10(0.5^3) - (x[,"n_12"] + x[,"n_21"])*log10(0.5*(8 - 12*0.5 + 9*0.5^2)) - (x[,"n_11"] + x[,"n_22"])*log10(5 - 11*0.5 + 15*0.5^2 - 9*0.5^3)
  
  LOD_cm <- (x[,"n_00"] + x[,"n_33"])*log10((-1 + r_cm)^2*r_cm) + (x[,"n_03"] + x[,"n_30"])*log10(-((-1 + r_cm)*r_cm^2)) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(r_cm*(3 - 4*r_cm + 3*r_cm^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log10(4 - 4*r_cm + 13*r_cm^2 - 9*r_cm^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 - 4*r_cm + 5*r_cm^2 - 3*r_cm^3) + 
    (x[,"n_11"] + x[,"n_22"])*log10(4 + 5*r_cm - 14*r_cm^2 + 9*r_cm^3) - (x[,"n_00"] + x[,"n_33"])*log10((-1 + 0.5)^2*0.5) - (x[,"n_03"] + x[,"n_30"])*log10(-((-1 + 0.5)*0.5^2)) - 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(0.5*(3 - 4*0.5 + 3*0.5^2)) - (x[,"n_12"] + x[,"n_21"])*log10(4 - 4*0.5 + 13*0.5^2 - 9*0.5^3) - 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 - 4*0.5 + 5*0.5^2 - 3*0.5^3) - (x[,"n_11"] + x[,"n_22"])*log10(4 + 5*0.5 - 14*0.5^2 + 9*0.5^3)
  
  LOD_cr <- (x[,"n_03"] + x[,"n_30"])*log10((-1 + r_cr)^2*r_cr) + (x[,"n_00"] + x[,"n_33"])*log10(-((-1 + r_cr)*r_cr^2)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(r_cr*(2 - 4*r_cr + 3*r_cr^2)) +
    (x[,"n_12"] + x[,"n_21"])*log10(r_cr*(9 - 14*r_cr + 9*r_cr^2)) + (x[,"n_11"] + x[,"n_22"])*log10(4 - 8*r_cr + 13*r_cr^2 - 9*r_cr^3) + 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(1 - 3*r_cr + 5*r_cr^2 - 3*r_cr^3) - (x[,"n_03"] + x[,"n_30"])*log10((-1 + 0.5)^2*0.5) - (x[,"n_00"] + x[,"n_33"])*log10(-((-1 + 0.5)*0.5^2)) - 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(0.5*(2 - 4*0.5 + 3*0.5^2)) - (x[,"n_12"] + x[,"n_21"])*log10(0.5*(9 - 14*0.5 + 9*0.5^2)) - 
    (x[,"n_11"] + x[,"n_22"])*log10(4 - 8*0.5 + 13*0.5^2 - 9*0.5^3) - (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(1 - 3*0.5 + 5*0.5^2 - 3*0.5^3)
  
  LOD_rc <- (x[,"n_03"] + x[,"n_30"])*log10(-((-2 + r_rc)*r_rc^2)) + (x[,"n_00"] + x[,"n_33"])*log10((-1 + r_rc)^2*(1 + r_rc)) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(r_rc*(4 - 5*r_rc + 3*r_rc^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log10(8 - 8*r_rc + 14*r_rc^2 - 9*r_rc^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 - 3*r_rc + 4*r_rc^2 - 3*r_rc^3) + 
    (x[,"n_11"] + x[,"n_22"])*log10(5 + 7*r_rc - 13*r_rc^2 + 9*r_rc^3) - (x[,"n_03"] + x[,"n_30"])*log10(-((-2 + 0.5)*0.5^2)) - (x[,"n_00"] + x[,"n_33"])*log10((-1 + 0.5)^2*(1 + 0.5)) -
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(0.5*(4 - 5*0.5 + 3*0.5^2)) - (x[,"n_12"] + x[,"n_21"])*log10(8 - 8*0.5 + 14*0.5^2 - 9*0.5^3) - 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 - 3*0.5 + 4*0.5^2 - 3*0.5^3) - (x[,"n_11"] + x[,"n_22"])*log10(5 + 7*0.5 - 13*0.5^2 + 9*0.5^3)
  
  LOD_rm <- (x[,"n_03"] + x[,"n_30"])*log10((-2 + r_rm)*(-1 + r_rm)*r_rm) + (x[,"n_00"] + x[,"n_33"])*log10(-((-1 + r_rm)*r_rm*(1 + r_rm))) + (x[,"n_11"] + x[,"n_22"])*log10(4 - r_rm + 4*r_rm^2 - 3*r_rm^3) +
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(4 - 5*r_rm + 6*r_rm^2 - 3*r_rm^3) + (x[,"n_12"] + x[,"n_21"])*log10(4 + 2*r_rm - 5*r_rm^2 + 3*r_rm^3) + 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 + 2*r_rm - 3*r_rm^2 + 3*r_rm^3) - (x[,"n_03"] + x[,"n_30"])*log10((-2 + 0.5)*(-1 + 0.5)*0.5) - 
    (x[,"n_00"] + x[,"n_33"])*log10(-((-1 + 0.5)*0.5*(1 + 0.5))) - (x[,"n_11"] + x[,"n_22"])*log10(4 - 0.5 + 4*0.5^2 - 3*0.5^3) - 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(4 - 5*0.5 + 6*0.5^2 - 3*0.5^3) - (x[,"n_12"] + x[,"n_21"])*log10(4 + 2*0.5 - 5*0.5^2 + 3*0.5^3) - 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 + 2*0.5 - 3*0.5^2 + 3*0.5^3)
  
  LOD_rr <- (x[,"n_03"] + x[,"n_30"])*log10(-((-2 + r_rr)*(-1 + r_rr)^2)) + (x[,"n_00"] + x[,"n_33"])*log10(r_rr^2*(1 + r_rr)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(r_rr*(2 + 2*r_rr - 3*r_rr^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log10(10 - 13*r_rr + 16*r_rr^2 - 9*r_rr^3) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(1 + 3*r_rr - 7*r_rr^2 + 3*r_rr^3) + 
    (x[,"n_11"] + x[,"n_22"])*log10(4 + 8*r_rr - 11*r_rr^2 + 9*r_rr^3) - (x[,"n_03"] + x[,"n_30"])*log10(-((-2 + 0.5)*(-1 + 0.5)^2)) - (x[,"n_00"] + x[,"n_33"])*log10(0.5^2*(1 + 0.5)) -
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(0.5*(2 + 2*0.5 - 3*0.5^2)) - (x[,"n_12"] + x[,"n_21"])*log10(10 - 13*0.5 + 16*0.5^2 - 9*0.5^3) - 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(1 + 3*0.5 - 7*0.5^2 + 3*0.5^3) - (x[,"n_11"] + x[,"n_22"])*log10(4 + 8*0.5 - 11*0.5^2 + 9*0.5^3)
  
  ## Record the logL also:
  logL_cc <- (-x[,"n_00"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_21"] - x[,"n_22"] - x[,"n_30"] - x[,"n_33"])*log(3) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(4) + 
    (x[,"n_00"] + x[,"n_33"])*log(-(-1 + r_cc)^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log((-1 + r_cc)^2*r_cc) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(-((-1 + r_cc)*r_cc^2)) + 
    (x[,"n_03"] + x[,"n_30"])*log(r_cc^3) + (x[,"n_12"] + x[,"n_21"])*log(r_cc*(8 - 12*r_cc + 9*r_cc^2)) + (x[,"n_11"] + x[,"n_22"])*log(5 - 11*r_cc + 15*r_cc^2 - 9*r_cc^3)
  
  logL_cm <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(24) + 
    (x[,"n_00"] + x[,"n_33"])*log((-1 + r_cm)^2*r_cm) + (x[,"n_03"] + x[,"n_30"])*log(-((-1 + r_cm)*r_cm^2)) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(r_cm*(3 - 4*r_cm + 3*r_cm^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log(4 - 4*r_cm + 13*r_cm^2 - 9*r_cm^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(2 - 4*r_cm + 5*r_cm^2 - 3*r_cm^3) + 
    (x[,"n_11"] + x[,"n_22"])*log(4 + 5*r_cm - 14*r_cm^2 + 9*r_cm^3)
  
  logL_cr <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(12) + 
    (x[,"n_03"] + x[,"n_30"])*log((-1 + r_cr)^2*r_cr) + (x[,"n_00"] + x[,"n_33"])*log(-((-1 + r_cr)*r_cr^2)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(r_cr*(2 - 4*r_cr + 3*r_cr^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log(r_cr*(9 - 14*r_cr + 9*r_cr^2)) + (x[,"n_11"] + x[,"n_22"])*log(4 - 8*r_cr + 13*r_cr^2 - 9*r_cr^3) + 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(1 - 3*r_cr + 5*r_cr^2 - 3*r_cr^3)
  
  logL_rc <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(36) + 
    (x[,"n_03"] + x[,"n_30"])*log(-((-2 + r_rc)*r_rc^2)) + (x[,"n_00"] + x[,"n_33"])*log((-1 + r_rc)^2*(1 + r_rc)) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(r_rc*(4 - 5*r_rc + 3*r_rc^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log(8 - 8*r_rc + 14*r_rc^2 - 9*r_rc^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(2 - 3*r_rc + 4*r_rc^2 - 3*r_rc^3) + 
    (x[,"n_11"] + x[,"n_22"])*log(5 + 7*r_rc - 13*r_rc^2 + 9*r_rc^3)
  
  logL_rm <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - 2*x[,"n_23"] - 2*x[,"n_30"] - 2*x[,"n_31"] - 2*x[,"n_32"] - 2*x[,"n_33"])*log(3) + 
    (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(8) + 
    (x[,"n_03"] + x[,"n_30"])*log((-2 + r_rm)*(-1 + r_rm)*r_rm) + (x[,"n_00"] + x[,"n_33"])*log(-((-1 + r_rm)*r_rm*(1 + r_rm))) + (x[,"n_11"] + x[,"n_22"])*log(4 - r_rm + 4*r_rm^2 - 3*r_rm^3) + 
    (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(4 - 5*r_rm + 6*r_rm^2 - 3*r_rm^3) + (x[,"n_12"] + x[,"n_21"])*log(4 + 2*r_rm - 5*r_rm^2 + 3*r_rm^3) + 
    (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(2 + 2*r_rm - 3*r_rm^2 + 3*r_rm^3)
  
  logL_rr <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(36) + 
    (x[,"n_03"] + x[,"n_30"])*log(-((-2 + r_rr)*(-1 + r_rr)^2)) + (x[,"n_00"] + x[,"n_33"])*log(r_rr^2*(1 + r_rr)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(r_rr*(2 + 2*r_rr - 3*r_rr^2)) + 
    (x[,"n_12"] + x[,"n_21"])*log(10 - 13*r_rr + 16*r_rr^2 - 9*r_rr^3) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(1 + 3*r_rr - 7*r_rr^2 + 3*r_rr^3) + 
    (x[,"n_11"] + x[,"n_22"])*log(4 + 8*r_rr - 11*r_rr^2 + 9*r_rr^3)
  
  return(list(r_mat=cbind(r_cc, r_cm, r_cr, r_rc, r_rm, r_rr),
              LOD_mat=cbind(LOD_cc, LOD_cm, LOD_cr, LOD_rc, LOD_rm, LOD_rr),
              logL_mat=cbind(logL_cc, logL_cm, logL_cr, logL_rc, logL_rm, logL_rr),
              phasing_strategy="MLL", 
              possible_phases=c("coupling coupling", 
                                "mixed coupling",
                                "repulsion coupling",
                                "coupling repulsion",
                                "mixed repulsion",
                                "repulsion repulsion")))
  
}

#' @rdname r4_functions
#' @noRd
r4_1.3_1.3<-function(x, ncores=1){
  
  #########################
  ## COUPLING COUPLING
  #########################
  
  logLcc<-function(r,n11,n12,n13,n21,n22,n23,n31,n32,n33){ 
    L <- (-2*n11 - n12 - 2*n13 - n21 - n23 - 2*n31 - n32 - 2*n33)*log(2) + (n11 + n33)*log((-1 + r)^2) + 
      (n12 + n21 + n23 + n32)*log(-((-1 + r)*r)) + (n13 + n31)*log(r^2) + n22*log(1/2 - r + r^2)
    return(L)}
  
  inter_logLcc <- function(n11,n12,n13,n21,n22,n23,n31,n32,n33){
    optimize(logLcc,c(0,0.5),n11,n12,n13,n21,n22,n23,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcc<-function(r){ #Coupling function
  #     L <- 2^(-2*x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - x[,"n_21"] - x[,"n_23"] - 2*x[,"n_31"] - x[,"n_32"] - 2*x[,"n_33"])*((-1 + r)^2)^(x[,"n_11"] + x[,"n_33"])*(-((-1 + r)*r))^(x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*(r^2)^(x[,"n_13"] + x[,"n_31"])*(1/2 - r + r^2)^x[,"n_22"]
  #     return(L)}
  #   
  #########################
  ## COUPLING MIXED
  #########################
  logLcm<-function(r,n11,n12,n13,n21,n22,n23,n31,n32,n33){ 
    L <- (-2*n11 - n12 - 2*n13 - n21 - n22 - n23 - 2*n31 - n32 - 2*n33)*log(2) + 
      (-n11 - n12 - n13 - n21 - n22 - n23 - n31 - n32 - n33)*log(3) + 
      (n13 + n31)*log(-((-2 + r)*r)) + (n11 + n33)*log(-((-1 + r)*(1 + r))) + n22*log(1 + 2*r - 2*r^2) + 
      (n12 + n21 + n23 + n32)*log(1 - r + r^2)
    return(L)}
  
  inter_logLcm <- function(n11,n12,n13,n21,n22,n23,n31,n32,n33){
    optimize(logLcm,c(0,0.5),n11,n12,n13,n21,n22,n23,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lcm<-function(r){ 
  #     L <- 2^(-2*x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - 2*x[,"n_31"] - x[,"n_32"] - 2*x[,"n_33"])*3^(-x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*(-((-2 + r)*r))^(x[,"n_13"] + x[,"n_31"])*(-((-1 + r)*(1 + r)))^(x[,"n_11"] + x[,"n_33"])*(1 + 2*r - 2*r^2)^x[,"n_22"]*(1 - r + r^2)^(x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])
  #     return(L)}
  
  #########################
  ## REPULSION COUPLING
  #########################
  
  logLrc<-function(r,n11,n12,n13,n21,n22,n23,n31,n32,n33){ 
    L <- (-2*n11 - n12 - 2*n13 - n21 - n22 - n23 - 2*n31 - n32 - 2*n33)*log(2) + 
      (-n11 - n12 - n13 - n21 - n22 - n23 - n31 - n32 - n33)*log(3) + 
      (n13 + n31)*log(-((-2 + r)*r)) + (n11 + n33)*log(-((-1 + r)*(1 + r))) + n22*log(1 + 2*r - 2*r^2) + (n12 + n21 + n23 + n32)*log(1 - r + r^2)
    return(L)}
  
  inter_logLrc <- function(n11,n12,n13,n21,n22,n23,n31,n32,n33){
    optimize(logLrc,c(0,0.5),n11,n12,n13,n21,n22,n23,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  
  #   Lrc<-function(r){ 
  #     L <- 2^(-2*x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - 2*x[,"n_31"] - x[,"n_32"] - 2*x[,"n_33"])*3^(-x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*(-((-2 + r)*r))^(x[,"n_13"] + x[,"n_31"])*(-((-1 + r)*(1 + r)))^(x[,"n_11"] + x[,"n_33"])*(1 + 2*r - 2*r^2)^x[,"n_22"]*(1 - r + r^2)^(x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])
  #     return(L)}
  
  #########################
  ## REPULSION MIXED
  #########################
  
  
  logLrm<-function(r,n11,n12,n13,n21,n22,n23,n31,n32,n33){
    L <- (-2*n11 - n12 - 2*n13 - n21 - n22 - n23 - 2*n31 - n32 - 2*n33)*log(2) + 
      (-n11 - n12 - n13 - n21 - n22 - n23 - n31 - n32 - n33)*log(9) + 
      (n13 + n31)*log((-2 + r)^2) + (n12 + n21 + n23 + n32)*log(-((-2 + r)*(1 + r))) + (n11 + n33)*log((1 + r)^2) + 
      n22*log(5 - 2*r + 2*r^2)
    return(L)}
  
  inter_logLrm <- function(n11,n12,n13,n21,n22,n23,n31,n32,n33){
    optimize(logLrm,c(0,0.5),n11,n12,n13,n21,n22,n23,n31,n32,n33,
             maximum=T,lower = 0, upper = 0.5)$maximum}
  
  #   Lrm<-function(r){ 
  #     L <- 2^(-2*x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - 2*x[,"n_31"] - x[,"n_32"] - 2*x[,"n_33"])*9^(-x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*((-2 + r)^2)^(x[,"n_13"] + x[,"n_31"])*(-((-2 + r)*(1 + r)))^(x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*((1 + r)^2)^(x[,"n_11"] + x[,"n_33"])*(5 - 2*r + 2*r^2)^x[,"n_22"]
  #     return(L)}
  
  r_cc <- parallel::mcmapply(inter_logLcc, x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_cm <- parallel::mcmapply(inter_logLcm,x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                             mc.cores = ncores)
  
  r_rc <-  parallel::mcmapply(inter_logLrc, x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  r_rm <-  parallel::mcmapply(inter_logLrm, x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_31"],x[,"n_32"],x[,"n_33"],
                              mc.cores = ncores)
  
  ## Estimate the log L values:
  logL_cc <- (-2*x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - x[,"n_21"] - x[,"n_23"] - 2*x[,"n_31"] - x[,"n_32"] - 2*x[,"n_33"])*log(2) + (x[,"n_11"] + x[,"n_33"])*log((-1 + r_cc)^2) + 
    (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log(-((-1 + r_cc)*r_cc)) + (x[,"n_13"] + x[,"n_31"])*log(r_cc^2) + x[,"n_22"]*log(1/2 - r_cc + r_cc^2)
  
  logL_cm <- (-2*x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - 2*x[,"n_31"] - x[,"n_32"] - 2*x[,"n_33"])*log(2) + 
    (-x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(3) + 
    (x[,"n_13"] + x[,"n_31"])*log(-((-2 + r_cm)*r_cm)) + (x[,"n_11"] + x[,"n_33"])*log(-((-1 + r_cm)*(1 + r_cm))) + x[,"n_22"]*log(1 + 2*r_cm - 2*r_cm^2) + 
    (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log(1 - r_cm + r_cm^2)
  
  logL_rc <- (-2*x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - 2*x[,"n_31"] - x[,"n_32"] - 2*x[,"n_33"])*log(2) + 
    (-x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(3) + 
    (x[,"n_13"] + x[,"n_31"])*log(-((-2 + r_rc)*r_rc)) + (x[,"n_11"] + x[,"n_33"])*log(-((-1 + r_rc)*(1 + r_rc))) + x[,"n_22"]*log(1 + 2*r_rc - 2*r_rc^2) + 
    (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log(1 - r_rc + r_rc^2)
  
  logL_rm <- (-2*x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - 2*x[,"n_31"] - x[,"n_32"] - 2*x[,"n_33"])*log(2) + 
    (-x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(9) + 
    (x[,"n_13"] + x[,"n_31"])*log((-2 + r_rm)^2) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log(-((-2 + r_rm)*(1 + r_rm))) + (x[,"n_11"] + x[,"n_33"])*log((1 + r_rm)^2) + 
    x[,"n_22"]*log(5 - 2*r_rm + 2*r_rm^2)
  
  #############################################################################################
  LOD_cc <- (x[,"n_11"] + x[,"n_33"])*log10((-1 + r_cc)^2) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(-((-1 + r_cc)*r_cc)) + (x[,"n_13"] + x[,"n_31"])*log10(r_cc^2) + 
    x[,"n_22"]*log10(1/2 - r_cc + r_cc^2) - (x[,"n_11"] + x[,"n_33"])*log10((-1 + 0.5)^2) - (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(-((-1 + 0.5)*0.5)) - 
    (x[,"n_13"] + x[,"n_31"])*log10(0.5^2) - x[,"n_22"]*log10(1/2 - 0.5 + 0.5^2)
  
  LOD_cm <- (x[,"n_13"] + x[,"n_31"])*log10(-((-2 + r_cm)*r_cm)) + (x[,"n_11"] + x[,"n_33"])*log10(-((-1 + r_cm)*(1 + r_cm))) + x[,"n_22"]*log10(1 + 2*r_cm - 2*r_cm^2) + 
    (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(1 - r_cm + r_cm^2) - (x[,"n_13"] + x[,"n_31"])*log10(-((-2 + 0.5)*0.5)) - (x[,"n_11"] + x[,"n_33"])*log10(-((-1 + 0.5)*(1 + 0.5))) - 
    x[,"n_22"]*log10(1 + 2*0.5 - 2*0.5^2) - (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(1 - 0.5 + 0.5^2)
  
  LOD_rc <- (x[,"n_13"] + x[,"n_31"])*log10(-((-2 + r_rc)*r_rc)) + (x[,"n_11"] + x[,"n_33"])*log10(-((-1 + r_rc)*(1 + r_rc))) + x[,"n_22"]*log10(1 + 2*r_rc - 2*r_rc^2) + 
    (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(1 - r_rc + r_rc^2) - (x[,"n_13"] + x[,"n_31"])*log10(-((-2 + 0.5)*0.5)) - (x[,"n_11"] + x[,"n_33"])*log10(-((-1 + 0.5)*(1 + 0.5))) - 
    x[,"n_22"]*log10(1 + 2*0.5 - 2*0.5^2) - (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(1 - 0.5 + 0.5^2)
  
  LOD_rm <- (x[,"n_13"] + x[,"n_31"])*log10((-2 + r_rm)^2) + (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(-((-2 + r_rm)*(1 + r_rm))) + (x[,"n_11"] + x[,"n_33"])*log10((1 + r_rm)^2) + 
    x[,"n_22"]*log10(5 - 2*r_rm + 2*r_rm^2) - (x[,"n_13"] + x[,"n_31"])*log10((-2 + 0.5)^2) - (x[,"n_12"] + x[,"n_21"] + x[,"n_23"] + x[,"n_32"])*log10(-((-2 + 0.5)*(1 + 0.5))) - 
    (x[,"n_11"] + x[,"n_33"])*log10((1 + 0.5)^2) - x[,"n_22"]*log10(5 - 2*0.5 + 2*0.5^2)
  
  return(list(r_mat=cbind(r_cc, r_cm, r_rc, r_rm),
              LOD_mat=cbind(LOD_cc, LOD_cm, LOD_rc, LOD_rm),
              logL_mat=cbind(logL_cc, logL_cm, logL_rc, logL_rm),
              phasing_strategy="MLL", 
              possible_phases=c("coupling coupling", 
                                "coupling mixed",
                                "repulsion coupling",
                                "repulsion mixed")))
}
