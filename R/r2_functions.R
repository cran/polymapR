############################################################################################
## Likelihood, LOD and recombination frequency functions for a diploid x diploid cross 
## Peter Bourke, Wageningen UR Plant Breeding. January 2019
############################################################################################

#' Calculate recombination frequency, LOD and log-likelihood from frequency tables in a random pairing diploid cross.
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
#' @name r2_functions
NULL

#' @rdname r2_functions
r2_1.0_1.0 <- function(x,ncores=1){
r_c <- (x[,"n_01"] + x[,"n_10"])/(x[,"n_00"] + x[,"n_01"] + x[,"n_10"] + x[,"n_11"])
logL_c <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_10"] - x[,"n_11"])*log(2) + (x[,"n_00"] + x[,"n_11"])*log(pmax(1e-6,1 - r_c)) + (x[,"n_01"] + x[,"n_10"])*log(pmax(1e-6,r_c))
LOD_c <- (x[,"n_01"] + x[,"n_10"])*log10(2) + (x[,"n_00"] + x[,"n_11"])*log10(2) + (x[,"n_00"] + x[,"n_11"])*log10(pmax(1e-6,1 - r_c)) + (x[,"n_01"] + x[,"n_10"])*log10(pmax(1e-6,r_c))


r_r <- (x[,"n_00"] + x[,"n_11"])/(x[,"n_00"] + x[,"n_01"] + x[,"n_10"] + x[,"n_11"])
logL_r <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_10"] - x[,"n_11"])*log(2) + (x[,"n_01"] + x[,"n_10"])*log(pmax(1e-6,1 - r_r)) + (x[,"n_00"] + x[,"n_11"])*log(pmax(1e-6,r_r))
LOD_r <- (x[,"n_01"] + x[,"n_10"])*log10(2) + (x[,"n_00"] + x[,"n_11"])*log10(2) + (x[,"n_01"] + x[,"n_10"])*log10(pmax(1e-6,1 - r_r)) + (x[,"n_00"] + x[,"n_11"])*log10(pmax(1e-6,r_r))


return(list(
r_mat = cbind(r_c,r_r,0.499),
LOD_mat = cbind(LOD_c,LOD_r,0),
logL_mat = cbind(logL_c,logL_r,-1e6),
phasing_strategy = "MLL",
possible_phases = c("coupling","repulsion","unknown")
)
)
}

#' @rdname r2_functions
r2_1.0_1.1 <- function(x,ncores=1){
r_c <- (x[,"n_02"] + x[,"n_10"])/(x[,"n_00"] + x[,"n_02"] + x[,"n_10"] + x[,"n_12"])
logL_c <- 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"])*log(2) + (x[,"n_00"] + x[,"n_12"])*log(pmax(1e-6,1 - r_c)) + (x[,"n_02"] + x[,"n_10"])*log(pmax(1e-6,r_c))
LOD_c <- (x[,"n_02"] + x[,"n_10"])*log10(2) + (x[,"n_00"] + x[,"n_12"])*log10(2) + (x[,"n_00"] + x[,"n_12"])*log10(pmax(1e-6,1 - r_c)) + (x[,"n_02"] + x[,"n_10"])*log10(pmax(1e-6,r_c))


r_r <- (x[,"n_00"] + x[,"n_12"])/(x[,"n_00"] + x[,"n_02"] + x[,"n_10"] + x[,"n_12"])
logL_r <- 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"])*log(2) + (x[,"n_02"] + x[,"n_10"])*log(pmax(1e-6,1 - r_r)) + (x[,"n_00"] + x[,"n_12"])*log(pmax(1e-6,r_r))
LOD_r <- (x[,"n_02"] + x[,"n_10"])*log10(2) + (x[,"n_00"] + x[,"n_12"])*log10(2) + (x[,"n_02"] + x[,"n_10"])*log10(pmax(1e-6,1 - r_r)) + (x[,"n_00"] + x[,"n_12"])*log10(pmax(1e-6,r_r))


return(list(
r_mat = cbind(r_c,r_r,0.499),
LOD_mat = cbind(LOD_c,LOD_r,0),
logL_mat = cbind(logL_c,logL_r,-1e6),
phasing_strategy = "MLL",
possible_phases = c("coupling","repulsion","unknown")
)
)
}

#' @rdname r2_functions
r2_1.1_1.1 <- function(x,ncores=1){
logL_cc <- function(r,n00,n01,n02,n10,n12,n20,n21,n22,n11) {
L <- (-2*n00 - n01 - 2*n02 - n10 - n12 - 2*n20 - n21 - 2*n22)*log(2) + 2*(n00 + n22)*log(pmax(1e-6,1 - r)) + 2*(n02 + n20)*log(pmax(1e-6,r)) + (n01 + n10 + n12 + n21)*(log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + n11*log(1/2 - r + r^2)
return(L)}
interlogL_cc <- function(n00,n01,n02,n10,n12,n20,n21,n22,n11) {
optimize(logL_cc,c(0,0.5), n00,n01,n02,n10,n12,n20,n21,n22,n11, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_cc <- parallel::mcmapply(interlogL_cc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_11"],mc.cores = ncores)


LOD_cc <- 2*x[,"n_11"]*log10(2) + 2*(x[,"n_02"] + x[,"n_20"])*log10(2) + 2*(x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(2) + 2*(x[,"n_00"] + x[,"n_22"])*log10(2) + 2*(x[,"n_00"] + x[,"n_22"])*log10(pmax(1e-6,1 - r_cc)) + 2*(x[,"n_02"] + x[,"n_20"])*log10(pmax(1e-6,r_cc)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*(log10(pmax(1e-6,1 - r_cc)) + log10(pmax(1e-6,r_cc))) + x[,"n_11"]*log10(1/2 - r_cc + r_cc^2)


logL_cc <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*log(2) + 2*(x[,"n_00"] + x[,"n_22"])*log(pmax(1e-6,1 - r_cc)) + 2*(x[,"n_02"] + x[,"n_20"])*log(pmax(1e-6,r_cc)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*(log(pmax(1e-6,1 - r_cc)) + log(pmax(1e-6,r_cc))) + x[,"n_11"]*log(1/2 - r_cc + r_cc^2)


logL_cr <- function(r,n00,n01,n02,n10,n12,n20,n21,n22,n11) {
L <- 2*(-n00 - n01 - n02 - n10 - n12 - n20 - n21 - n22)*log(2) + (n00 + n02 + n11 + n20 + n22)*(log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n01 + n10 + n12 + n21)*log(1 - 2*r + 2*r^2)
return(L)}
interlogL_cr <- function(n00,n01,n02,n10,n12,n20,n21,n22,n11) {
optimize(logL_cr,c(0,0.5), n00,n01,n02,n10,n12,n20,n21,n22,n11, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_cr <- parallel::mcmapply(interlogL_cr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_11"],mc.cores = ncores)


LOD_cr <- (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(2) + 2*(x[,"n_00"] + x[,"n_02"] + x[,"n_11"] + x[,"n_20"] + x[,"n_22"])*log10(2) + (x[,"n_00"] + x[,"n_02"] + x[,"n_11"] + x[,"n_20"] + x[,"n_22"])*(log10(pmax(1e-6,1 - r_cr)) + log10(pmax(1e-6,r_cr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(1 - 2*r_cr + 2*r_cr^2)


logL_cr <- 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*log(2) + (x[,"n_00"] + x[,"n_02"] + x[,"n_11"] + x[,"n_20"] + x[,"n_22"])*(log(pmax(1e-6,1 - r_cr)) + log(pmax(1e-6,r_cr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log(1 - 2*r_cr + 2*r_cr^2)


logL_rr <- function(r,n00,n01,n02,n10,n12,n20,n21,n22,n11) {
L <- (-2*n00 - n01 - 2*n02 - n10 - n12 - 2*n20 - n21 - 2*n22)*log(2) + 2*(n02 + n20)*log(pmax(1e-6,1 - r)) + 2*(n00 + n22)*log(pmax(1e-6,r)) + (n01 + n10 + n12 + n21)*(log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + n11*log(1/2 - r + r^2)
return(L)}
interlogL_rr <- function(n00,n01,n02,n10,n12,n20,n21,n22,n11) {
optimize(logL_rr,c(0,0.5), n00,n01,n02,n10,n12,n20,n21,n22,n11, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_rr <- parallel::mcmapply(interlogL_rr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_11"],mc.cores = ncores)


LOD_rr <- 2*x[,"n_11"]*log10(2) + 2*(x[,"n_02"] + x[,"n_20"])*log10(2) + 2*(x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(2) + 2*(x[,"n_00"] + x[,"n_22"])*log10(2) + 2*(x[,"n_02"] + x[,"n_20"])*log10(pmax(1e-6,1 - r_rr)) + 2*(x[,"n_00"] + x[,"n_22"])*log10(pmax(1e-6,r_rr)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*(log10(pmax(1e-6,1 - r_rr)) + log10(pmax(1e-6,r_rr))) + x[,"n_11"]*log10(1/2 - r_rr + r_rr^2)


logL_rr <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*log(2) + 2*(x[,"n_02"] + x[,"n_20"])*log(pmax(1e-6,1 - r_rr)) + 2*(x[,"n_00"] + x[,"n_22"])*log(pmax(1e-6,r_rr)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*(log(pmax(1e-6,1 - r_rr)) + log(pmax(1e-6,r_rr))) + x[,"n_11"]*log(1/2 - r_rr + r_rr^2)


return(list(
r_mat = cbind(r_cc,r_cr,r_rr,0.499),
LOD_mat = cbind(LOD_cc,LOD_cr,LOD_rr,0),
logL_mat = cbind(logL_cc,logL_cr,logL_rr,-1e6),
phasing_strategy = "MLL",
possible_phases = c("coupling coupling","mixed","repulsion repulsion","unknown")
)
)
}
