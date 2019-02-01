############################################################################################
## Likelihood, LOD and recombination frequency functions for a tetraploid x diploid cross 
## Peter Bourke, Wageningen UR Plant Breeding. January 2019
############################################################################################

#' Calculate recombination frequency, LOD and log-likelihood from frequency tables in a random pairing triploid from a tetraploid x diploid cross.
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
#' @name r3_functions
NULL


#' @rdname r3_functions
r3_0.1_0.1 <- function(x,ncores=1){
r_c <- (x[,"n_01"] + x[,"n_10"])/(x[,"n_00"] + x[,"n_01"] + x[,"n_10"] + x[,"n_11"])
logL_c <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_10"] - x[,"n_11"])*log(2) + (x[,"n_00"] + x[,"n_11"])*log(pmax(1e-6,1 - r_c)) + (x[,"n_01"] + x[,"n_10"])*log(pmax(1e-6,r_c))
LOD_c <- (x[,"n_01"] + x[,"n_10"])*log10(2) + (x[,"n_00"] + x[,"n_11"])*log10(2) + (x[,"n_00"] + x[,"n_11"])*log10(pmax(1e-6,1 - r_c)) + (x[,"n_01"] + x[,"n_10"])*log10(pmax(1e-6,r_c))


logL_r <- function(r,n00,n01,n10,n11) {
L <- (-n00 - n01 - n10 - n11)*log(2) + (n01 + n10)*log(pmax(1e-6,1 - r)) + (n00 + n11)*log(pmax(1e-6,r))
return(L)}
interlogL_r <- function(n00,n01,n10,n11) {
optimize(logL_r,c(0,0.5), n00,n01,n10,n11, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_r <- parallel::mcmapply(interlogL_r,x[,"n_00"],x[,"n_01"],x[,"n_10"],x[,"n_11"],mc.cores = ncores)


LOD_r <- (x[,"n_01"] + x[,"n_10"])*log10(2) + (x[,"n_00"] + x[,"n_11"])*log10(2) + (x[,"n_01"] + x[,"n_10"])*log10(pmax(1e-6,1 - r_r)) + (x[,"n_00"] + x[,"n_11"])*log10(pmax(1e-6,r_r))


logL_r <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_10"] - x[,"n_11"])*log(2) + (x[,"n_01"] + x[,"n_10"])*log(pmax(1e-6,1 - r_r)) + (x[,"n_00"] + x[,"n_11"])*log(pmax(1e-6,r_r))


return(list(
r_mat = cbind(r_c,r_r,0.499),
LOD_mat = cbind(LOD_c,LOD_r,0),
logL_mat = cbind(logL_c,logL_r,-1e6),
phasing_strategy = "MLL",
possible_phases = c("coupling","repulsion","unknown")
)
)
}

#' @rdname r3_functions
r3_0.1_1.1 <- function(x,ncores=1){
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

#' @rdname r3_functions
r3_0.1_2.1 <- function(x,ncores=1){
logL_c <- function(r,n00,n02,n03,n10,n11,n13,n01,n12) {
L <- (-n00 - n02 - n03 - n10 - n11 - n13)*(2*log(2) + log(3)) + (n00 + n13)*log(pmax(1e-6,1 - r)) + (n01 + n12)*log(1/3 - r/4) + (n03 + n10)*log(pmax(1e-6,r)) + (n02 + n11)*log(1 + 3*r)
return(L)}
interlogL_c <- function(n00,n02,n03,n10,n11,n13,n01,n12) {
optimize(logL_c,c(0,0.5), n00,n02,n03,n10,n11,n13,n01,n12, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_c <- parallel::mcmapply(interlogL_c,x[,"n_00"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_13"],x[,"n_01"],x[,"n_12"],mc.cores = ncores)


LOD_c <- (x[,"n_03"] + x[,"n_10"])*log10(2) + (x[,"n_00"] + x[,"n_13"])*log10(2) + (x[,"n_01"] + x[,"n_12"])*(3*log10(2) + log10(3) - log10(5)) - (x[,"n_02"] + x[,"n_11"])*(-log10(2) + log10(5)) + (x[,"n_00"] + x[,"n_13"])*log10(pmax(1e-6,1 - r_c)) + (x[,"n_01"] + x[,"n_12"])*log10(1/3 - r_c/4) + (x[,"n_03"] + x[,"n_10"])*log10(pmax(1e-6,r_c)) + (x[,"n_02"] + x[,"n_11"])*log10(1 + 3*r_c)


logL_c <- (-x[,"n_00"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_13"])*(2*log(2) + log(3)) + (x[,"n_00"] + x[,"n_13"])*log(pmax(1e-6,1 - r_c)) + (x[,"n_01"] + x[,"n_12"])*log(1/3 - r_c/4) + (x[,"n_03"] + x[,"n_10"])*log(pmax(1e-6,r_c)) + (x[,"n_02"] + x[,"n_11"])*log(1 + 3*r_c)


logL_r <- function(r,n00,n01,n03,n10,n12,n13,n02,n11) {
L <- (-n00 - n01 - n03 - n10 - n12 - n13)*(2*log(2) + log(3)) + (n03 + n10)*log(pmax(1e-6,1 - r)) + (n02 + n11)*log(1/3 - r/4) + (n00 + n13)*log(pmax(1e-6,r)) + (n01 + n12)*log(1 + 3*r)
return(L)}
interlogL_r <- function(n00,n01,n03,n10,n12,n13,n02,n11) {
optimize(logL_r,c(0,0.5), n00,n01,n03,n10,n12,n13,n02,n11, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_r <- parallel::mcmapply(interlogL_r,x[,"n_00"],x[,"n_01"],x[,"n_03"],x[,"n_10"],x[,"n_12"],x[,"n_13"],x[,"n_02"],x[,"n_11"],mc.cores = ncores)


LOD_r <- (x[,"n_03"] + x[,"n_10"])*log10(2) + (x[,"n_00"] + x[,"n_13"])*log10(2) + (x[,"n_02"] + x[,"n_11"])*(3*log10(2) + log10(3) - log10(5)) - (x[,"n_01"] + x[,"n_12"])*(-log10(2) + log10(5)) + (x[,"n_03"] + x[,"n_10"])*log10(pmax(1e-6,1 - r_r)) + (x[,"n_02"] + x[,"n_11"])*log10(1/3 - r_r/4) + (x[,"n_00"] + x[,"n_13"])*log10(pmax(1e-6,r_r)) + (x[,"n_01"] + x[,"n_12"])*log10(1 + 3*r_r)


logL_r <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_03"] - x[,"n_10"] - x[,"n_12"] - x[,"n_13"])*(2*log(2) + log(3)) + (x[,"n_03"] + x[,"n_10"])*log(pmax(1e-6,1 - r_r)) + (x[,"n_02"] + x[,"n_11"])*log(1/3 - r_r/4) + (x[,"n_00"] + x[,"n_13"])*log(pmax(1e-6,r_r)) + (x[,"n_01"] + x[,"n_12"])*log(1 + 3*r_r)


return(list(
r_mat = cbind(r_c,r_r,0.499),
LOD_mat = cbind(LOD_c,LOD_r,0),
logL_mat = cbind(logL_c,logL_r,-1e6),
phasing_strategy = "MLL",
possible_phases = c("coupling","repulsion","unknown")
)
)
}

#' @rdname r3_functions
#' @noRd
r3_1.0_1.0 <- function(x,ncores=1) r4_1.0_1.0(x,ncores) #Directly pass to r4 function

#' @rdname r3_functions
#' @noRd
r3_1.0_1.1 <- function(x,ncores=1) r4_1.0_1.1(x,ncores) #Directly pass to r4 function

#' @rdname r3_functions
#' @noRd
r3_1.0_2.0 <- function(x,ncores=1) r4_1.0_2.0(x,ncores) #Directly pass to r4 function

#' @rdname r3_functions
#' @noRd
r3_1.0_2.1 <- function(x,ncores=1) r4_1.0_2.1(x,ncores) #Directly pass to r4 function

#' @rdname r3_functions
#' @noRd
r3_1.1_1.1 <- function(x,ncores=1){
logL_cc <- function(r,n00,n01,n02,n10,n12,n20,n21,n22,n11) {
L <- (-2*n00 - n01 - 2*n02 - n10 - n12 - 2*n20 - n21 - 2*n22)*log(2) + (n00 + n22)*log((-1 + r)^2) + 2*(n02 + n20)*log(pmax(1e-6,r)) + (n01 + n10 + n12 + n21)*(log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + n11*log(1/2 - r + r^2)
return(L)}
interlogL_cc <- function(n00,n01,n02,n10,n12,n20,n21,n22,n11) {
optimize(logL_cc,c(0,0.5), n00,n01,n02,n10,n12,n20,n21,n22,n11, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_cc <- parallel::mcmapply(interlogL_cc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_11"],mc.cores = ncores)


LOD_cc <- 2*x[,"n_11"]*log10(2) + 2*(x[,"n_02"] + x[,"n_20"])*log10(2) + 2*(x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(2) + 2*(x[,"n_00"] + x[,"n_22"])*log10(2) + 2*(x[,"n_00"] + x[,"n_22"])*log10(-1 + r_cc) + 2*(x[,"n_02"] + x[,"n_20"])*log10(pmax(1e-6,r_cc)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*(log10(pmax(1e-6,1 - r_cc)) + log10(pmax(1e-6,r_cc))) + x[,"n_11"]*log10(1/2 - r_cc + r_cc^2)


logL_cc <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*log(2) + (x[,"n_00"] + x[,"n_22"])*log((-1 + r_cc)^2) + 2*(x[,"n_02"] + x[,"n_20"])*log(pmax(1e-6,r_cc)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*(log(pmax(1e-6,1 - r_cc)) + log(pmax(1e-6,r_cc))) + x[,"n_11"]*log(1/2 - r_cc + r_cc^2)


logL_cr <- function(r,n00,n01,n02,n10,n12,n20,n21,n22,n11) {
L <- 2*(-n00 - n01 - n02 - n10 - n12 - n20 - n21 - n22)*log(2) + (n00 + n02 + n11 + n20 + n22)*(log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n01 + n10 + n12 + n21)*log(1 - 2*r + 2*r^2)
return(L)}
interlogL_cr <- function(n00,n01,n02,n10,n12,n20,n21,n22,n11) {
optimize(logL_cr,c(0,0.5), n00,n01,n02,n10,n12,n20,n21,n22,n11, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_cr <- parallel::mcmapply(interlogL_cr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_11"],mc.cores = ncores)


LOD_cr <- (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(2) + 2*(x[,"n_00"] + x[,"n_02"] + x[,"n_11"] + x[,"n_20"] + x[,"n_22"])*log10(2) + (x[,"n_00"] + x[,"n_02"] + x[,"n_11"] + x[,"n_20"] + x[,"n_22"])*(log10(pmax(1e-6,1 - r_cr)) + log10(pmax(1e-6,r_cr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(1 - 2*r_cr + 2*r_cr^2)


logL_cr <- 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*log(2) + (x[,"n_00"] + x[,"n_02"] + x[,"n_11"] + x[,"n_20"] + x[,"n_22"])*(log(pmax(1e-6,1 - r_cr)) + log(pmax(1e-6,r_cr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log(1 - 2*r_cr + 2*r_cr^2)


logL_rc <- function(r,n00,n01,n02,n10,n11,n12,n20,n21,n22) {
L <- (-2*n00 - n01 - 2*n02 - n10 - n11 - n12 - 2*n20 - n21 - 2*n22)*log(2) + (-n00 - n01 - n02 - n10 - n11 - n12 - n20 - n21 - n22)*log(3) + (n02 + n20)*(log(2 - r) + log(pmax(1e-6,r))) + (n00 + n22)*(log(pmax(1e-6,1 - r)) + log(1 + r)) + n11*log(1 + 2*r - 2*r^2) + (n01 + n10 + n12 + n21)*log(1 - r + r^2)
return(L)}
interlogL_rc <- function(n00,n01,n02,n10,n11,n12,n20,n21,n22) {
optimize(logL_rc,c(0,0.5), n00,n01,n02,n10,n11,n12,n20,n21,n22, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_rc <- parallel::mcmapply(interlogL_rc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],mc.cores = ncores)


LOD_rc <- (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*(2*log10(2) - log10(3)) - (x[,"n_02"] + x[,"n_20"])*(-2*log10(2) + log10(3)) - (x[,"n_00"] + x[,"n_22"])*(-2*log10(2) + log10(3)) - x[,"n_11"]*(-log10(2) + log10(3)) + (x[,"n_02"] + x[,"n_20"])*(log10(2 - r_rc) + log10(pmax(1e-6,r_rc))) + (x[,"n_00"] + x[,"n_22"])*(log10(pmax(1e-6,1 - r_rc)) + log10(1 + r_rc)) + x[,"n_11"]*log10(1 + 2*r_rc - 2*r_rc^2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(1 - r_rc + r_rc^2)


logL_rc <- (-2*x[,"n_00"] - x[,"n_01"] - 2*x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_20"] - x[,"n_21"] - 2*x[,"n_22"])*log(2) + (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*log(3) + (x[,"n_02"] + x[,"n_20"])*(log(2 - r_rc) + log(pmax(1e-6,r_rc))) + (x[,"n_00"] + x[,"n_22"])*(log(pmax(1e-6,1 - r_rc)) + log(1 + r_rc)) + x[,"n_11"]*log(1 + 2*r_rc - 2*r_rc^2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log(1 - r_rc + r_rc^2)


logL_rr <- function(r,n00,n01,n02,n10,n12,n20,n21,n22,n11) {
L <- 2*(-n00 - n01 - n02 - n10 - n12 - n20 - n21 - n22)*log(2) + (-n00 - n01 - n02 - n10 - n11 - n12 - n20 - n21 - n22)*log(3) + (n02 + n20)*(log(pmax(1e-6,1 - r)) + log(2 - r)) + (n00 + n22)*(log(pmax(1e-6,r)) + log(1 + r)) + (n01 + n10 + n12 + n21)*log(1 + 2*r - 2*r^2) + n11*log(1 - r + r^2)
return(L)}
interlogL_rr <- function(n00,n01,n02,n10,n12,n20,n21,n22,n11) {
optimize(logL_rr,c(0,0.5), n00,n01,n02,n10,n12,n20,n21,n22,n11, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_rr <- parallel::mcmapply(interlogL_rr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_10"],x[,"n_12"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_11"],mc.cores = ncores)


LOD_rr <- x[,"n_11"]*(2*log10(2) - log10(3)) - (x[,"n_02"] + x[,"n_20"])*(-2*log10(2) + log10(3)) - (x[,"n_00"] + x[,"n_22"])*(-2*log10(2) + log10(3)) - (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*(-log10(2) + log10(3)) + (x[,"n_02"] + x[,"n_20"])*(log10(pmax(1e-6,1 - r_rr)) + log10(2 - r_rr)) + (x[,"n_00"] + x[,"n_22"])*(log10(pmax(1e-6,r_rr)) + log10(1 + r_rr)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log10(1 + 2*r_rr - 2*r_rr^2) + x[,"n_11"]*log10(1 - r_rr + r_rr^2)


logL_rr <- 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*log(2) + (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"])*log(3) + (x[,"n_02"] + x[,"n_20"])*(log(pmax(1e-6,1 - r_rr)) + log(2 - r_rr)) + (x[,"n_00"] + x[,"n_22"])*(log(pmax(1e-6,r_rr)) + log(1 + r_rr)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_12"] + x[,"n_21"])*log(1 + 2*r_rr - 2*r_rr^2) + x[,"n_11"]*log(1 - r_rr + r_rr^2)


return(list(
r_mat = cbind(r_cc,r_cr,r_rc,r_rr,0.499),
LOD_mat = cbind(LOD_cc,LOD_cr,LOD_rc,LOD_rr,0),
logL_mat = cbind(logL_cc,logL_cr,logL_rc,logL_rr,-1e6),
phasing_strategy = "MLL",
possible_phases = c("coupling coupling","coupling repulsion","repulsion coupling","repulsion repulsion","unknown")
)
)
}

#' @rdname r3_functions
#' @noRd
r3_1.1_2.1 <- function(x,ncores=1){
logL_cc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23) {
L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - n10 - 2*n11 - 2*n12 - n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(3) + (n00 + n23)*log((-1 + r)^2) + 2*(n03 + n20)*log(pmax(1e-6,r)) + (n10 + n13)*(log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n02 + n21)*(log(3 - r) + log(pmax(1e-6,r))) + (n01 + n22)*log(2 - r - r^2) + (n11 + n12)*log(3 - 2*r + 2*r^2)
return(L)}
interlogL_cc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23) {
optimize(logL_cc,c(0,0.5), n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_cc <- parallel::mcmapply(interlogL_cc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],mc.cores = ncores)


LOD_cc <- 2*(x[,"n_10"] + x[,"n_13"])*log10(2) + 2*(x[,"n_03"] + x[,"n_20"])*log10(2) + 2*(x[,"n_00"] + x[,"n_23"])*log10(2) - (x[,"n_02"] + x[,"n_21"])*(-2*log10(2) + log10(5)) - (x[,"n_01"] + x[,"n_22"])*(-2*log10(2) + log10(5)) - (x[,"n_11"] + x[,"n_12"])*(-log10(2) + log10(5)) + 2*(x[,"n_00"] + x[,"n_23"])*log10(-1 + r_cc) + 2*(x[,"n_03"] + x[,"n_20"])*log10(pmax(1e-6,r_cc)) + (x[,"n_10"] + x[,"n_13"])*(log10(pmax(1e-6,1 - r_cc)) + log10(pmax(1e-6,r_cc))) + (x[,"n_02"] + x[,"n_21"])*(log10(3 - r_cc) + log10(pmax(1e-6,r_cc))) + (x[,"n_01"] + x[,"n_22"])*log10(2 - r_cc - r_cc^2) + (x[,"n_11"] + x[,"n_12"])*log10(3 - 2*r_cc + 2*r_cc^2)


logL_cc <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"])*log(2) + (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(3) + (x[,"n_00"] + x[,"n_23"])*log((-1 + r_cc)^2) + 2*(x[,"n_03"] + x[,"n_20"])*log(pmax(1e-6,r_cc)) + (x[,"n_10"] + x[,"n_13"])*(log(pmax(1e-6,1 - r_cc)) + log(pmax(1e-6,r_cc))) + (x[,"n_02"] + x[,"n_21"])*(log(3 - r_cc) + log(pmax(1e-6,r_cc))) + (x[,"n_01"] + x[,"n_22"])*log(2 - r_cc - r_cc^2) + (x[,"n_11"] + x[,"n_12"])*log(3 - 2*r_cc + 2*r_cc^2)


logL_cr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23) {
L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - 2*n10 - n11 - n12 - 2*n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(3) + (n00 + n03 + n20 + n23)*(log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n11 + n12)*log(1 + r - r^2) + (n01 + n22)*log(1 + r^2) + (n02 + n21)*log(2 - 2*r + r^2) + (n10 + n13)*log(1 - 2*r + 2*r^2)
return(L)}
interlogL_cr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23) {
optimize(logL_cr,c(0,0.5), n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_cr <- parallel::mcmapply(interlogL_cr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],mc.cores = ncores)


LOD_cr <- (x[,"n_10"] + x[,"n_13"])*log10(2) + 2*(x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*log10(2) - (x[,"n_11"] + x[,"n_12"])*(-2*log10(2) + log10(5)) - (x[,"n_02"] + x[,"n_21"])*(-2*log10(2) + log10(5)) - (x[,"n_01"] + x[,"n_22"])*(-2*log10(2) + log10(5)) + (x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*(log10(pmax(1e-6,1 - r_cr)) + log10(pmax(1e-6,r_cr))) + (x[,"n_11"] + x[,"n_12"])*log10(1 + r_cr - r_cr^2) + (x[,"n_01"] + x[,"n_22"])*log10(1 + r_cr^2) + (x[,"n_02"] + x[,"n_21"])*log10(2 - 2*r_cr + r_cr^2) + (x[,"n_10"] + x[,"n_13"])*log10(1 - 2*r_cr + 2*r_cr^2)


logL_cr <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"])*log(2) + (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(3) + (x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*(log(pmax(1e-6,1 - r_cr)) + log(pmax(1e-6,r_cr))) + (x[,"n_11"] + x[,"n_12"])*log(1 + r_cr - r_cr^2) + (x[,"n_01"] + x[,"n_22"])*log(1 + r_cr^2) + (x[,"n_02"] + x[,"n_21"])*log(2 - 2*r_cr + r_cr^2) + (x[,"n_10"] + x[,"n_13"])*log(1 - 2*r_cr + 2*r_cr^2)


logL_rc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23) {
L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - 2*n10 - n11 - n12 - 2*n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(3) + (n00 + n03 + n20 + n23)*(log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n11 + n12)*log(1 + r - r^2) + (n02 + n21)*log(1 + r^2) + (n01 + n22)*log(2 - 2*r + r^2) + (n10 + n13)*log(1 - 2*r + 2*r^2)
return(L)}
interlogL_rc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23) {
optimize(logL_rc,c(0,0.5), n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_rc <- parallel::mcmapply(interlogL_rc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],mc.cores = ncores)


LOD_rc <- (x[,"n_10"] + x[,"n_13"])*log10(2) + 2*(x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*log10(2) - (x[,"n_11"] + x[,"n_12"])*(-2*log10(2) + log10(5)) - (x[,"n_02"] + x[,"n_21"])*(-2*log10(2) + log10(5)) - (x[,"n_01"] + x[,"n_22"])*(-2*log10(2) + log10(5)) + (x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*(log10(pmax(1e-6,1 - r_rc)) + log10(pmax(1e-6,r_rc))) + (x[,"n_11"] + x[,"n_12"])*log10(1 + r_rc - r_rc^2) + (x[,"n_02"] + x[,"n_21"])*log10(1 + r_rc^2) + (x[,"n_01"] + x[,"n_22"])*log10(2 - 2*r_rc + r_rc^2) + (x[,"n_10"] + x[,"n_13"])*log10(1 - 2*r_rc + 2*r_rc^2)


logL_rc <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - 2*x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - 2*x[,"n_13"] - 2*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"])*log(2) + (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(3) + (x[,"n_00"] + x[,"n_03"] + x[,"n_20"] + x[,"n_23"])*(log(pmax(1e-6,1 - r_rc)) + log(pmax(1e-6,r_rc))) + (x[,"n_11"] + x[,"n_12"])*log(1 + r_rc - r_rc^2) + (x[,"n_02"] + x[,"n_21"])*log(1 + r_rc^2) + (x[,"n_01"] + x[,"n_22"])*log(2 - 2*r_rc + r_rc^2) + (x[,"n_10"] + x[,"n_13"])*log(1 - 2*r_rc + 2*r_rc^2)


logL_rr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23) {
L <- (-2*n00 - 2*n01 - 2*n02 - 2*n03 - n10 - 2*n11 - 2*n12 - n13 - 2*n20 - 2*n21 - 2*n22 - 2*n23)*log(2) + (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23)*log(3) + (n03 + n20)*log((-1 + r)^2) + 2*(n00 + n23)*log(pmax(1e-6,r)) + (n10 + n13)*(log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n01 + n22)*(log(3 - r) + log(pmax(1e-6,r))) + (n02 + n21)*log(2 - r - r^2) + (n11 + n12)*log(3 - 2*r + 2*r^2)
return(L)}
interlogL_rr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23) {
optimize(logL_rr,c(0,0.5), n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_rr <- parallel::mcmapply(interlogL_rr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],mc.cores = ncores)


LOD_rr <- 2*(x[,"n_10"] + x[,"n_13"])*log10(2) + 2*(x[,"n_03"] + x[,"n_20"])*log10(2) + 2*(x[,"n_00"] + x[,"n_23"])*log10(2) - (x[,"n_02"] + x[,"n_21"])*(-2*log10(2) + log10(5)) - (x[,"n_01"] + x[,"n_22"])*(-2*log10(2) + log10(5)) - (x[,"n_11"] + x[,"n_12"])*(-log10(2) + log10(5)) + 2*(x[,"n_03"] + x[,"n_20"])*log10(-1 + r_rr) + 2*(x[,"n_00"] + x[,"n_23"])*log10(pmax(1e-6,r_rr)) + (x[,"n_10"] + x[,"n_13"])*(log10(pmax(1e-6,1 - r_rr)) + log10(pmax(1e-6,r_rr))) + (x[,"n_01"] + x[,"n_22"])*(log10(3 - r_rr) + log10(pmax(1e-6,r_rr))) + (x[,"n_02"] + x[,"n_21"])*log10(2 - r_rr - r_rr^2) + (x[,"n_11"] + x[,"n_12"])*log10(3 - 2*r_rr + 2*r_rr^2)


logL_rr <- (-2*x[,"n_00"] - 2*x[,"n_01"] - 2*x[,"n_02"] - 2*x[,"n_03"] - x[,"n_10"] - 2*x[,"n_11"] - 2*x[,"n_12"] - x[,"n_13"] - 2*x[,"n_20"] - 2*x[,"n_21"] - 2*x[,"n_22"] - 2*x[,"n_23"])*log(2) + (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"])*log(3) + (x[,"n_03"] + x[,"n_20"])*log((-1 + r_rr)^2) + 2*(x[,"n_00"] + x[,"n_23"])*log(pmax(1e-6,r_rr)) + (x[,"n_10"] + x[,"n_13"])*(log(pmax(1e-6,1 - r_rr)) + log(pmax(1e-6,r_rr))) + (x[,"n_01"] + x[,"n_22"])*(log(3 - r_rr) + log(pmax(1e-6,r_rr))) + (x[,"n_02"] + x[,"n_21"])*log(2 - r_rr - r_rr^2) + (x[,"n_11"] + x[,"n_12"])*log(3 - 2*r_rr + 2*r_rr^2)


return(list(
r_mat = cbind(r_cc,r_cr,r_rc,r_rr,0.499),
LOD_mat = cbind(LOD_cc,LOD_cr,LOD_rc,LOD_rr,0),
logL_mat = cbind(logL_cc,logL_cr,logL_rc,logL_rr,-1e6),
phasing_strategy = "MLL",
possible_phases = c("coupling coupling","coupling repulsion","repulsion coupling","repulsion repulsion","unknown")
)
)
}

#' @rdname r3_functions
#' @noRd
r3_2.0_2.0 <- function(x,ncores=1) r4_2.0_2.0(x,ncores) #Directly pass to r4 function

#' @rdname r3_functions
#' @noRd
r3_2.0_2.1 <- function(x,ncores=1) r4_2.0_2.1(x,ncores) #Directly pass to r4 function

#' @rdname r3_functions
#' @noRd
r3_2.0_1.1 <- function(x,ncores=1) r4_2.0_1.1(x,ncores) #Directly pass to r4 function

#' @rdname r3_functions
r3_2.1_2.1 <- function(x,ncores=1){
logL_cc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
L <- 2*(-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(2) + (-n00 - n03 - n11 - n12 - n21 - n22 - n30 - n33)*log(3) + 3*(n00 + n33)*log(pmax(1e-6,1 - r)) + 3*(n03 + n30)*log(pmax(1e-6,r)) + (n01 + n10 + n23 + n32)*(2*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n02 + n13 + n20 + n31)*(log(pmax(1e-6,1 - r)) + 2*log(pmax(1e-6,r))) + (n12 + n21)*(log(pmax(1e-6,r)) + log(8 - 12*r + 9*r^2)) + (n11 + n22)*log(5 - 11*r + 15*r^2 - 9*r^3)
return(L)}
interlogL_cc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
optimize(logL_cc,c(0,0.5), n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_cc <- parallel::mcmapply(interlogL_cc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],mc.cores = ncores)


LOD_cc <- 3*(x[,"n_03"] + x[,"n_30"])*log10(2) + 3*(x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(2) + 3*(x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2) + 3*(x[,"n_00"] + x[,"n_33"])*log10(2) - (x[,"n_12"] + x[,"n_21"])*(-3*log10(2) + log10(17)) - (x[,"n_11"] + x[,"n_22"])*(-3*log10(2) + log10(17)) + 3*(x[,"n_00"] + x[,"n_33"])*log10(pmax(1e-6,1 - r_cc)) + 3*(x[,"n_03"] + x[,"n_30"])*log10(pmax(1e-6,r_cc)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(2*log10(pmax(1e-6,1 - r_cc)) + log10(pmax(1e-6,r_cc))) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(log10(pmax(1e-6,1 - r_cc)) + 2*log10(pmax(1e-6,r_cc))) + (x[,"n_12"] + x[,"n_21"])*(log10(pmax(1e-6,r_cc)) + log10(8 - 12*r_cc + 9*r_cc^2)) + (x[,"n_11"] + x[,"n_22"])*log10(5 - 11*r_cc + 15*r_cc^2 - 9*r_cc^3)


logL_cc <- 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(2) + (-x[,"n_00"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_21"] - x[,"n_22"] - x[,"n_30"] - x[,"n_33"])*log(3) + 3*(x[,"n_00"] + x[,"n_33"])*log(pmax(1e-6,1 - r_cc)) + 3*(x[,"n_03"] + x[,"n_30"])*log(pmax(1e-6,r_cc)) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(2*log(pmax(1e-6,1 - r_cc)) + log(pmax(1e-6,r_cc))) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(log(pmax(1e-6,1 - r_cc)) + 2*log(pmax(1e-6,r_cc))) + (x[,"n_12"] + x[,"n_21"])*(log(pmax(1e-6,r_cc)) + log(8 - 12*r_cc + 9*r_cc^2)) + (x[,"n_11"] + x[,"n_22"])*log(5 - 11*r_cc + 15*r_cc^2 - 9*r_cc^3)


logL_cr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*(2*log(2) + log(3)) + (n00 + n33)*(2*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n03 + n30)*(log(pmax(1e-6,1 - r)) + 2*log(pmax(1e-6,r))) + (n02 + n13 + n20 + n31)*(log(pmax(1e-6,r)) + log(2 - 4*r + 3*r^2)) + (n11 + n22)*(log(pmax(1e-6,r)) + log(9 - 14*r + 9*r^2)) + (n12 + n21)*log(4 - 8*r + 13*r^2 - 9*r^3) + (n01 + n10 + n23 + n32)*log(1 - 3*r + 5*r^2 - 3*r^3)
return(L)}
interlogL_cr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
optimize(logL_cr,c(0,0.5), n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_cr <- parallel::mcmapply(interlogL_cr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],mc.cores = ncores)


LOD_cr <- 3*(x[,"n_03"] + x[,"n_30"])*log10(2) + 3*(x[,"n_00"] + x[,"n_33"])*log10(2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(3*log10(2) - log10(3)) - (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(-3*log10(2) + log10(3)) - (x[,"n_12"] + x[,"n_21"])*(-3*log10(2) + log10(17)) - (x[,"n_11"] + x[,"n_22"])*(-3*log10(2) + log10(17)) + (x[,"n_00"] + x[,"n_33"])*(2*log10(pmax(1e-6,1 - r_cr)) + log10(pmax(1e-6,r_cr))) + (x[,"n_03"] + x[,"n_30"])*(log10(pmax(1e-6,1 - r_cr)) + 2*log10(pmax(1e-6,r_cr))) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(log10(pmax(1e-6,r_cr)) + log10(2 - 4*r_cr + 3*r_cr^2)) + (x[,"n_11"] + x[,"n_22"])*(log10(pmax(1e-6,r_cr)) + log10(9 - 14*r_cr + 9*r_cr^2)) + (x[,"n_12"] + x[,"n_21"])*log10(4 - 8*r_cr + 13*r_cr^2 - 9*r_cr^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(1 - 3*r_cr + 5*r_cr^2 - 3*r_cr^3)


logL_cr <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*(2*log(2) + log(3)) + (x[,"n_00"] + x[,"n_33"])*(2*log(pmax(1e-6,1 - r_cr)) + log(pmax(1e-6,r_cr))) + (x[,"n_03"] + x[,"n_30"])*(log(pmax(1e-6,1 - r_cr)) + 2*log(pmax(1e-6,r_cr))) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(log(pmax(1e-6,r_cr)) + log(2 - 4*r_cr + 3*r_cr^2)) + (x[,"n_11"] + x[,"n_22"])*(log(pmax(1e-6,r_cr)) + log(9 - 14*r_cr + 9*r_cr^2)) + (x[,"n_12"] + x[,"n_21"])*log(4 - 8*r_cr + 13*r_cr^2 - 9*r_cr^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(1 - 3*r_cr + 5*r_cr^2 - 3*r_cr^3)


logL_mc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*(3*log(2) + log(3)) + (n00 + n33)*(2*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n03 + n30)*(log(pmax(1e-6,1 - r)) + 2*log(pmax(1e-6,r))) + (n02 + n13 + n20 + n31)*(log(pmax(1e-6,r)) + log(3 - 4*r + 3*r^2)) + (n12 + n21)*log(4 - 4*r + 13*r^2 - 9*r^3) + (n01 + n10 + n23 + n32)*log(2 - 4*r + 5*r^2 - 3*r^3) + (n11 + n22)*log(4 + 5*r - 14*r^2 + 9*r^3)
return(L)}
interlogL_mc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
optimize(logL_mc,c(0,0.5), n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_mc <- parallel::mcmapply(interlogL_mc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],mc.cores = ncores)


LOD_mc <- 3*(x[,"n_03"] + x[,"n_30"])*log10(2) + 3*(x[,"n_00"] + x[,"n_33"])*log10(2) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(3*log10(2) - log10(7)) - (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(-3*log10(2) + log10(7)) - (x[,"n_12"] + x[,"n_21"])*(-3*log10(2) + log10(3) + log10(11)) - (x[,"n_11"] + x[,"n_22"])*(-3*log10(2) + log10(3) + log10(11)) + (x[,"n_00"] + x[,"n_33"])*(2*log10(pmax(1e-6,1 - r_mc)) + log10(pmax(1e-6,r_mc))) + (x[,"n_03"] + x[,"n_30"])*(log10(pmax(1e-6,1 - r_mc)) + 2*log10(pmax(1e-6,r_mc))) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(log10(pmax(1e-6,r_mc)) + log10(3 - 4*r_mc + 3*r_mc^2)) + (x[,"n_12"] + x[,"n_21"])*log10(4 - 4*r_mc + 13*r_mc^2 - 9*r_mc^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2 - 4*r_mc + 5*r_mc^2 - 3*r_mc^3) + (x[,"n_11"] + x[,"n_22"])*log10(4 + 5*r_mc - 14*r_mc^2 + 9*r_mc^3)


logL_mc <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*(3*log(2) + log(3)) + (x[,"n_00"] + x[,"n_33"])*(2*log(pmax(1e-6,1 - r_mc)) + log(pmax(1e-6,r_mc))) + (x[,"n_03"] + x[,"n_30"])*(log(pmax(1e-6,1 - r_mc)) + 2*log(pmax(1e-6,r_mc))) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(log(pmax(1e-6,r_mc)) + log(3 - 4*r_mc + 3*r_mc^2)) + (x[,"n_12"] + x[,"n_21"])*log(4 - 4*r_mc + 13*r_mc^2 - 9*r_mc^3) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log(2 - 4*r_mc + 5*r_mc^2 - 3*r_mc^3) + (x[,"n_11"] + x[,"n_22"])*log(4 + 5*r_mc - 14*r_mc^2 + 9*r_mc^3)


logL_mr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*(3*log(2) + log(3)) + (n03 + n30)*(2*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n00 + n33)*(log(pmax(1e-6,1 - r)) + 2*log(pmax(1e-6,r))) + (n01 + n10 + n23 + n32)*(log(pmax(1e-6,r)) + log(3 - 4*r + 3*r^2)) + (n11 + n22)*log(4 - 4*r + 13*r^2 - 9*r^3) + (n02 + n13 + n20 + n31)*log(2 - 4*r + 5*r^2 - 3*r^3) + (n12 + n21)*log(4 + 5*r - 14*r^2 + 9*r^3)
return(L)}
interlogL_mr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
optimize(logL_mr,c(0,0.5), n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_mr <- parallel::mcmapply(interlogL_mr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],mc.cores = ncores)


LOD_mr <- 3*(x[,"n_03"] + x[,"n_30"])*log10(2) + 3*(x[,"n_00"] + x[,"n_33"])*log10(2) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(3*log10(2) - log10(7)) - (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(-3*log10(2) + log10(7)) - (x[,"n_12"] + x[,"n_21"])*(-3*log10(2) + log10(3) + log10(11)) - (x[,"n_11"] + x[,"n_22"])*(-3*log10(2) + log10(3) + log10(11)) + (x[,"n_03"] + x[,"n_30"])*(2*log10(pmax(1e-6,1 - r_mr)) + log10(pmax(1e-6,r_mr))) + (x[,"n_00"] + x[,"n_33"])*(log10(pmax(1e-6,1 - r_mr)) + 2*log10(pmax(1e-6,r_mr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(log10(pmax(1e-6,r_mr)) + log10(3 - 4*r_mr + 3*r_mr^2)) + (x[,"n_11"] + x[,"n_22"])*log10(4 - 4*r_mr + 13*r_mr^2 - 9*r_mr^3) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(2 - 4*r_mr + 5*r_mr^2 - 3*r_mr^3) + (x[,"n_12"] + x[,"n_21"])*log10(4 + 5*r_mr - 14*r_mr^2 + 9*r_mr^3)


logL_mr <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*(3*log(2) + log(3)) + (x[,"n_03"] + x[,"n_30"])*(2*log(pmax(1e-6,1 - r_mr)) + log(pmax(1e-6,r_mr))) + (x[,"n_00"] + x[,"n_33"])*(log(pmax(1e-6,1 - r_mr)) + 2*log(pmax(1e-6,r_mr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(log(pmax(1e-6,r_mr)) + log(3 - 4*r_mr + 3*r_mr^2)) + (x[,"n_11"] + x[,"n_22"])*log(4 - 4*r_mr + 13*r_mr^2 - 9*r_mr^3) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(2 - 4*r_mr + 5*r_mr^2 - 3*r_mr^3) + (x[,"n_12"] + x[,"n_21"])*log(4 + 5*r_mr - 14*r_mr^2 + 9*r_mr^3)


logL_rc <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
L <- (-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*(2*log(2) + log(3)) + (n03 + n30)*(2*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n00 + n33)*(log(pmax(1e-6,1 - r)) + 2*log(pmax(1e-6,r))) + (n01 + n10 + n23 + n32)*(log(pmax(1e-6,r)) + log(2 - 4*r + 3*r^2)) + (n12 + n21)*(log(pmax(1e-6,r)) + log(9 - 14*r + 9*r^2)) + (n11 + n22)*log(4 - 8*r + 13*r^2 - 9*r^3) + (n02 + n13 + n20 + n31)*log(1 - 3*r + 5*r^2 - 3*r^3)
return(L)}
interlogL_rc <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
optimize(logL_rc,c(0,0.5), n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_rc <- parallel::mcmapply(interlogL_rc,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],mc.cores = ncores)


LOD_rc <- 3*(x[,"n_03"] + x[,"n_30"])*log10(2) + 3*(x[,"n_00"] + x[,"n_33"])*log10(2) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(3*log10(2) - log10(3)) - (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(-3*log10(2) + log10(3)) - (x[,"n_12"] + x[,"n_21"])*(-3*log10(2) + log10(17)) - (x[,"n_11"] + x[,"n_22"])*(-3*log10(2) + log10(17)) + (x[,"n_03"] + x[,"n_30"])*(2*log10(pmax(1e-6,1 - r_rc)) + log10(pmax(1e-6,r_rc))) + (x[,"n_00"] + x[,"n_33"])*(log10(pmax(1e-6,1 - r_rc)) + 2*log10(pmax(1e-6,r_rc))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(log10(pmax(1e-6,r_rc)) + log10(2 - 4*r_rc + 3*r_rc^2)) + (x[,"n_12"] + x[,"n_21"])*(log10(pmax(1e-6,r_rc)) + log10(9 - 14*r_rc + 9*r_rc^2)) + (x[,"n_11"] + x[,"n_22"])*log10(4 - 8*r_rc + 13*r_rc^2 - 9*r_rc^3) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(1 - 3*r_rc + 5*r_rc^2 - 3*r_rc^3)


logL_rc <- (-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*(2*log(2) + log(3)) + (x[,"n_03"] + x[,"n_30"])*(2*log(pmax(1e-6,1 - r_rc)) + log(pmax(1e-6,r_rc))) + (x[,"n_00"] + x[,"n_33"])*(log(pmax(1e-6,1 - r_rc)) + 2*log(pmax(1e-6,r_rc))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(log(pmax(1e-6,r_rc)) + log(2 - 4*r_rc + 3*r_rc^2)) + (x[,"n_12"] + x[,"n_21"])*(log(pmax(1e-6,r_rc)) + log(9 - 14*r_rc + 9*r_rc^2)) + (x[,"n_11"] + x[,"n_22"])*log(4 - 8*r_rc + 13*r_rc^2 - 9*r_rc^3) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log(1 - 3*r_rc + 5*r_rc^2 - 3*r_rc^3)


logL_rr <- function(r,n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
L <- 2*(-n00 - n01 - n02 - n03 - n10 - n11 - n12 - n13 - n20 - n21 - n22 - n23 - n30 - n31 - n32 - n33)*log(2) + (-n00 - n03 - n11 - n12 - n21 - n22 - n30 - n33)*log(3) + 3*(n03 + n30)*log(pmax(1e-6,1 - r)) + 3*(n00 + n33)*log(pmax(1e-6,r)) + (n02 + n13 + n20 + n31)*(2*log(pmax(1e-6,1 - r)) + log(pmax(1e-6,r))) + (n01 + n10 + n23 + n32)*(log(pmax(1e-6,1 - r)) + 2*log(pmax(1e-6,r))) + (n11 + n22)*(log(pmax(1e-6,r)) + log(8 - 12*r + 9*r^2)) + (n12 + n21)*log(5 - 11*r + 15*r^2 - 9*r^3)
return(L)}
interlogL_rr <- function(n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33) {
optimize(logL_rr,c(0,0.5), n00,n01,n02,n03,n10,n11,n12,n13,n20,n21,n22,n23,n30,n31,n32,n33, maximum=TRUE, lower=0, upper=0.5)$maximum}


r_rr <- parallel::mcmapply(interlogL_rr,x[,"n_00"],x[,"n_01"],x[,"n_02"],x[,"n_03"],x[,"n_10"],x[,"n_11"],x[,"n_12"],x[,"n_13"],x[,"n_20"],x[,"n_21"],x[,"n_22"],x[,"n_23"],x[,"n_30"],x[,"n_31"],x[,"n_32"],x[,"n_33"],mc.cores = ncores)


LOD_rr <- 3*(x[,"n_03"] + x[,"n_30"])*log10(2) + 3*(x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*log10(2) + 3*(x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*log10(2) + 3*(x[,"n_00"] + x[,"n_33"])*log10(2) - (x[,"n_12"] + x[,"n_21"])*(-3*log10(2) + log10(17)) - (x[,"n_11"] + x[,"n_22"])*(-3*log10(2) + log10(17)) + 3*(x[,"n_03"] + x[,"n_30"])*log10(pmax(1e-6,1 - r_rr)) + 3*(x[,"n_00"] + x[,"n_33"])*log10(pmax(1e-6,r_rr)) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(2*log10(pmax(1e-6,1 - r_rr)) + log10(pmax(1e-6,r_rr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(log10(pmax(1e-6,1 - r_rr)) + 2*log10(pmax(1e-6,r_rr))) + (x[,"n_11"] + x[,"n_22"])*(log10(pmax(1e-6,r_rr)) + log10(8 - 12*r_rr + 9*r_rr^2)) + (x[,"n_12"] + x[,"n_21"])*log10(5 - 11*r_rr + 15*r_rr^2 - 9*r_rr^3)


logL_rr <- 2*(-x[,"n_00"] - x[,"n_01"] - x[,"n_02"] - x[,"n_03"] - x[,"n_10"] - x[,"n_11"] - x[,"n_12"] - x[,"n_13"] - x[,"n_20"] - x[,"n_21"] - x[,"n_22"] - x[,"n_23"] - x[,"n_30"] - x[,"n_31"] - x[,"n_32"] - x[,"n_33"])*log(2) + (-x[,"n_00"] - x[,"n_03"] - x[,"n_11"] - x[,"n_12"] - x[,"n_21"] - x[,"n_22"] - x[,"n_30"] - x[,"n_33"])*log(3) + 3*(x[,"n_03"] + x[,"n_30"])*log(pmax(1e-6,1 - r_rr)) + 3*(x[,"n_00"] + x[,"n_33"])*log(pmax(1e-6,r_rr)) + (x[,"n_02"] + x[,"n_13"] + x[,"n_20"] + x[,"n_31"])*(2*log(pmax(1e-6,1 - r_rr)) + log(pmax(1e-6,r_rr))) + (x[,"n_01"] + x[,"n_10"] + x[,"n_23"] + x[,"n_32"])*(log(pmax(1e-6,1 - r_rr)) + 2*log(pmax(1e-6,r_rr))) + (x[,"n_11"] + x[,"n_22"])*(log(pmax(1e-6,r_rr)) + log(8 - 12*r_rr + 9*r_rr^2)) + (x[,"n_12"] + x[,"n_21"])*log(5 - 11*r_rr + 15*r_rr^2 - 9*r_rr^3)


return(list(
r_mat = cbind(r_cc,r_cr,r_mc,r_mr,r_rc,r_rr,0.499),
LOD_mat = cbind(LOD_cc,LOD_cr,LOD_mc,LOD_mr,LOD_rc,LOD_rr,0),
logL_mat = cbind(logL_cc,logL_cr,logL_mc,logL_mr,logL_rc,logL_rr,-1e6),
phasing_strategy = "MLL",
possible_phases = c("coupling coupling","coupling repulsion","mixed coupling","mixed repulsion","repulsion coupling","repulsion repulsion","unknown")
)
)
}
