#  Assign dosage scores to parents when using probabilistic genotype data
#' @param chk Output list as returned by function \code{\link{checkF1}}
#' @param probgeno_df A data.frame as read from the scores file produced by function
#' \code{saveMarkerModels} of R package \code{fitPoly}, or alternatively, a data.frame containing the following columns:
#' \itemize{
#' \item{SampleName}{
#' Name of the sample (individual)
#' }
#' \item{MarkerName}{
#' Name of the marker
#' }
#' \item{P0}{
#' Probabilities of dosage score '0'
#' }
#' \item{P1...}{
#' Probabilities of dosage score '1' etc. (up to max offspring dosage, e.g. P4 for tetraploid population)
#' }
#' \item{maxP}{
#' Maximum genotype probability identified for a particular individual and marker combination
#' }
#' \item{maxgeno}{
#' Most probable dosage for a particular individual and marker combination
#' }
#' \item{geno}{
#' Most probable dosage for a particular individual and marker combination, if \code{maxP} exceeds a user-defined threshold (e.g. 0.9), otherwise \code{NA}}
#' }
#' @param scorefile filename for tab-separated text file with the dosages. If NA no file is written
#' @return A data.frame with three columns: "MarkerName", "parent1" and "parent2".
#' @examples 
#' data("chk1","gp_df")
#' pardose<-assign_parental_dosage(chk=chk1,probgeno_df=gp_df)
#' @noRd
assign_parental_dosage <- function(chk,
                                   probgeno_df,
                                   scorefile = NA){
  
  probgeno_df <- test_probgeno_df(probgeno_df)
  
  result <- compareProbes.gp(chk = chk, 
                             scores = probgeno_df, 
                             probe.suffix = NA,
                             fracdiff.threshold = 0.04, 
                             qall_flavor = "qall_mult",
                             compfile = NA, 
                             combscorefile = scorefile)
  r <- result$combscores
  r <- r[,c("MarkerName","parent1","parent2")]
  r <- r[!is.na(r$parent1) & !is.na(r$parent2),]
  return(invisible(r))
} #assign_parental_dosage




#' calc_binning_thresholds
#' @description Calculate LOD and recombination frequency thresholds for marker binning
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @noRd
calc_binning_thresholds <- function(dosage_matrix){
  Ntot <- ncol(dosage_matrix) #exclude marker name and 2 parent cols
  the_na_rate <-
    sum(is.na(dosage_matrix)) / length(dosage_matrix) #0.02054193
  N_e <- Ntot * (1 - the_na_rate)
  r_thresh <- 1 / N_e
  lod_thresh <- 23.42072 + 0.115817 * N_e # See Bourke et al. (2016) for details
  
  out <- c(r_thresh, lod_thresh)
  names(out)<- c("r_thresh", "lod_thresh")
  
  return(out)
} #calc_binning_thresholds



#'@title Calculate potential shifts in marker dosages
#'@description Internal function for \code{correctDosages} function
#'@param segtypes Vector of segregation types
#'@param ploidy The ploidy of the species
#'@return Data.frame with columns:
#'\item{segtype}{ The segtype code before the underscore i.e. without the first dosage class
#'}
#'\item{mindos}{ The first dosage class
#'}
#'\item{shift}{ -1, 0 or 1: potental shifts to try out for each segtype whether a shift of +/- 1 
#'should be tried or not (0) (only when number of dosages in segtype in <= ploidy-2 and 
#'min.dosage==1 or max.dosage==ploidy-1)
#'}
#'@noRd
calcPotShifts <- function(segtypes, ploidy) {
  shift <- rep(0, length(segtypes))
  segtypes[is.na(segtypes) | segtypes==""] <- "x_0" #avoid problems with strsplit
  segtypes <- do.call(rbind, strsplit(as.character(segtypes), split="_"))
  numdos <- nchar(segtypes[,1]) #number of dosage classes in segtype
  compl <- substr(segtypes[,1],1,1) == "c"
  lencompl <- as.integer(substr(segtypes[compl, 1], 2, numdos[compl]-1))
  #nr of dosages in complex segtype
  numdos[compl] <- lencompl
  mindos <- as.integer(segtypes[,2])
  #shift down: numdos <= ploidy-2 and lowest dosage = 1:
  shift[numdos <= ploidy-2 & mindos == 1] <- -1
  #shift up: numdos <= ploidy-2 and highest dosage = ploidy-1:
  maxdos <- mindos + numdos - 1
  shift[numdos <= ploidy-2 & maxdos == ploidy-1] <- 1
  data.frame(seg=segtypes[,1], mindos=mindos, shift=shift,
             stringsAsFactors=FALSE)
} #calcPotShifts




#'calc_q1 ***********************************************************
#'@description Calculate q1 quality score. q1 is a measure of the fit of selfit. 
#'It is based on Pvalue and the fraction of valid F1 scores. If either is below their 
#'threshold the q1 value is 0. Also, if bestParentfit is different from bestfit this indicates a segregation distortion.
#'@param Pvalue P-value of bestParentfit segtype
#'@param fracInvalid Fraction invalid scores of bestParentfit segtype
#'@param bestfit segtype of best fit of population
#'@param bestParentfit segtype of best fit of parents
#'@param Pvalue_threshold Threshold P-value of bestParentfit segtype
#'@param fracInvalid_threshold a maximum threshold for the fracInvalid of the bestParentfit segtype
#' (with a larger fraction of invalid dosages in the F1 the q1 quality parameter will be set to 0)
#' @noRd
calc_q1 <- function(Pvalue, 
                    fracInvalid, 
                    bestfit, 
                    bestParentfit,
                    Pvalue_threshold, 
                    fracInvalid_threshold) {
  
  #x compares bestfit and bestParentfit:
  if (bestfit == bestParentfit) {
    x <- 1
  } else {
    x <- 0.5
  }
  #Version 20150220: nonzero but low y values also below P=0.001 depending
  #on threshold:
  if (Pvalue < Pvalue_threshold) y <- 0 else {
    if (Pvalue < 0.01) {
      if (Pvalue < 0.001) y <- linterpol(Pvalue, c(0,0), c(0.001, 0.05)) else
        y <- linterpol(Pvalue, c(0.001, 0.05), c(0.01, 0.6))
    } else if (Pvalue < 0.15) {
      if (Pvalue < 0.05) y <- linterpol(Pvalue, c(0.01, 0.6), c(0.05, 0.8)) else
        y <- linterpol(Pvalue, c(0.05, 0.8), c(0.15, 1))
    } else y <- 1
  }
  #z is determined by fraction F1 scores in valid categories for the selfit
  #segtype. z is equal to 1-fracInvalid, but below 1-frqinvalid.threshold z=0.
  z <- ifelse(fracInvalid > fracInvalid_threshold, 0, 1 - fracInvalid)
  x * y * z
} #calc_q1

#'calc_q2 ***********************************************************
#'@description q2 is about how well the parents are scored and fit with the selfit segtype
#'@param par.conflicts Amount of conflicts in parental scores
#'@param par.NAfrac Amount of missing values in the parents
#'@param matchParents Should parents match F1? One of the following: "No", "Unknown", "OneOK" or "Yes"
#'@param parentsScoredWithF1 Logical, if \code{TRUE} then q2 is determined not by matchParents but 
#'by the amount of conflicts and missing values in the parents. If \code{FALSE}, par.conflicts and 
#'par.NAfrac are not used and the information on the match between parents and F1 segtype has to come 
#'from matchParents.
#'@noRd
calc_q2 <- function(par.conflicts, 
                    par.NAfrac, 
                    matchParents,
                    parentsScoredWithF1) {
  
  if (parentsScoredWithF1) {
    q2 <- max(0, 1 - 0.5 * sum(par.conflicts) - 0.5 * sum(par.NAfrac))
    #In the version of 20150220 matchParents cannot be "No": selfit is always
    #bestParentfit (sometimes after promoting lowconf parental scores).
    #Therefore the remaining matchParents codes Unknown, OneOK and Yes are
    #correlated with the q2 score, no point to use the matchParents
  } else {
    # !parentsScoredWithF1
    # Note that matchParents=="No" cannot occur in the versions since 20150220.
    q2 <- match(matchParents,
                c("No", "Unknown", "OneOK", "Yes")) #1:4 or NA
    #q2 <-  c(0, 0.5, 0.9, 1)[q2] #translate into 0..1 values
    q2 <-  c(0, 0.65, 0.9, 1)[q2] #translate into 0..1 values
    #version of 20160324: the penalty for Unknown is less if the parents
    #are not scored with the F1: there is still the lack of information, but
    #the lack of scores does not mean that the (F1) genotyping is less reliable.
  }
  q2
} #calc_q2

#'calc_q3 ***********************************************************
#'@description   #q3 is about the fraction missing F1 scores. If that is larger than NA.threshold, q3=0;
#'else q3 depends on the fraction scored F1's
#'@param F1.NAfrac Fraction of missing values in the F1
#'@param fracNA_threshold Threshold for fraction missing values
#'@noRd
calc_q3 <- function(F1.NAfrac, 
                    fracNA_threshold) {
  
  if (F1.NAfrac > fracNA_threshold) {
    q3 <- 0
  } else if (fracNA_threshold > 0.2) {
    if (F1.NAfrac >= 0.2) {
      # 0.2<=F1.NAfrac<threshold
      q3 <- linterpol(F1.NAfrac, c(0.2, 0.5), c(fracNA_threshold, 0))
    } else {
      # 0<=F1.NAfrac<0.2
      q3 <- linterpol(F1.NAfrac, c(0, 1), c(0.2, 0.5))
    }
  } else {
    # fracNA_threshold <= 0.2, 0<=F1.NAfrac<threshold
    #: let q3 go from 1 at frac=0 to 0.5 at frac=threshold:
    q3 <- linterpol(F1.NAfrac, c(0, 1), c(fracNA_threshold, 0.5))
  }
  q3
} #calc_q3

#'calc_qall ***********************************************************
#'@param Pvalue_threshold Threshold can be disabled by setting to 0
#'@param fracInvalid_threshold Threshold can be disabled by setting to 1
#'@param fracNA_threshold Threshold can be disabled by setting to 1
#'@param Pvalue P-value of bestParentfit segtype
#'@param fracInvalid Fraction invalid scores of the bestParentfit segtype
#'@param F1.NAfrac fraction missing scores among F1 samples
#'@param matchParents Should parents match F1? One of the following: "No", "Unknown", "OneOK" or "Yes"
#'@param bestfit segtype of best fit of population
#'@param bestParentfit segtype of best fit of parents
#'@param par.conflicts Logical vector: is there a conflict for P1, P2
#'@param par.NAfrac number of NA scores for the samples of P1, P2
#'@param critweight numeric vector of length 3 (3 quality criteria) with their weights, 
#'or NA (in that case qall = q1*q2*q3, not a weighted average)
#'@param parentsScoredWithF1 Were parents scored with F1?
#'@return vector of length 4 or 5. The last (two) components of the vector are qall_mult 
#'(obtained by multiplication of q1..q3), and if critweights is supplied, qall_weights 
#'(a weighted average of q1..q3)
#'@noRd
calc_qall <- function(Pvalue_threshold, 
                      fracInvalid_threshold,
                      fracNA_threshold,
                      Pvalue, 
                      fracInvalid, 
                      F1.NAfrac,
                      matchParents,
                      bestfit, 
                      bestParentfit,
                      par.conflicts, 
                      par.NAfrac,
                      critweight,
                      parentsScoredWithF1) {
  
  q <- numeric(ifelse (length(critweight) == 3, 5, 4)) #last (two) are qall
  q[1] <- calc_q1(Pvalue, fracInvalid, bestfit, bestParentfit,
                  Pvalue_threshold, fracInvalid_threshold)
  q[2] <- calc_q2(par.conflicts, par.NAfrac, matchParents,
                  parentsScoredWithF1)
  q[3] <- calc_q3(F1.NAfrac, fracNA_threshold)
  q[4] <- prod(q[1:3]) #q[4] is qall_mult
  if (length(critweight) == 3) {
    q[5] <- sum(q[1:3] * critweight) / sum(critweight) #q[5] is qall_weights
  }
  q
} #calc_qall


#'checkFilename ***********************************************************
#'@description Check if a filename is valid for writing by trying to create the file.
#'If the filename contains a path, result will only be TRUE if the whole path already exists
#'@param filename The filename to check
#'@param overwrite if TRUE (default) an existing file of that name will be deleted (if it is not
#' protected by being locked, read-only, owned by another user etc)
#' @noRd
checkFilename <- function(filename, overwrite=TRUE) {
  filename <- filename[1] # we test only the first filename
  if (!overwrite && file.exists(filename)) return(FALSE)
  tryCatch(
    { suppressWarnings(
      if (file.create(filename)) {
        file.remove(filename)
        TRUE
      } else FALSE)
    }, error = function(e) { FALSE }
  )
} #checkFilename



#'chk2integer ***********************************************************
#'@description Convert data from checkF1 to integer
#'@param chk a data frame as returned by checkF1
#'@return the same data frame, but with all columns from Parent1 to last ancestor 
#'converted to integer (as these are sometimes factors with levels not equal to values 
#'after reading checkF1 from file)
#'@noRd
chk2integer <- function(chk) {
  #find out the ploidy:
  F1cap <- names(chk)[which(names(chk) == "F1_NA") - 1]
  ploidyF1 <- as.integer(sub("F1_", "", F1cap))
  #find the column range to convert to integers:
  firstcol <- which(names(chk) == "parent1")
  #determine which column has the last parent or ancestor sample:
  #if 3 columns for each segtype are present (showAll TRUE), then
  #the first column after the last sample has name frqInvalid_<segtype> (where
  #<segtype> depends on the ploidy).
  #else (showAll FALSE) the next column has name "bestfit", followed
  #by "frqInvalid_bestfit"
  lastcol <- which(leftstr(names(chk), 10) == "frqInvalid")[1]
  if (rightstr(names(chk)[lastcol], 7) == "bestfit") {
    lastcol <- lastcol - 2
  } else {
    lastcol <- lastcol - 1
  }
  if (length(firstcol) != 1 || is.na(lastcol) ||
      lastcol - firstcol < 3 + ploidyF1) {
    stop("chk2integer: column names of chk incorrect")
  }
  for (i in firstcol:lastcol) {
    suppressWarnings(chk[[i]] <- as.integer(as.character(chk[[i]])))
    #warnings caused by NAs
  }
  chk
} #chk2integer


#' Add colour bar scale to heatplot
#' @description Function to generate a scale for heatplots
#' @param col.data vector of colours
#' @param min minimum colour
#' @param max maximum colour
#' @param cex.ticks size of ticks on colour bar
#' @param nticks number of ticks on colour bar
#' @param ticks vector of positions of ticks on colour bar
#' @param title optional title for colour bar
#' @param ylab optional y-axis label for colour bar
#' @param cex.lab size of labels on colour bar
#' @noRd
colour.bar <- function(col.data, min, max=-min, cex.ticks = 1.2, nticks=11, 
                       ticks=seq(min, max, len=nticks), title='', ylab = '',
                       cex.lab = 1) {
  scale <- length(col.data)/(max-min)
  
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=ylab, main=title, cex.lab = cex.lab)
  axis(2, ticks, las=1, cex.axis = cex.ticks)
  for (i in 1:length(col.data)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=col.data[i], border=NA)
  }
} #colour.bar

#' Function to compare probes of the same marker and eventually assist merge
#' @description compareProbes.gp is used to compare probes of the same marker
#' @param chk Output list as returned by function \code{\link{checkF1}}
#' @param scores Input data, probabilistic genotype scores
#' @param probe.suffix Codes used to identify probes targetting same SNP
#' @param fracdiff.threshold Acceptable proportion of differences in probes to still be considered identical
#' @param parent1 Parent 1 identifier
#' @param parent2 Parent 2 identifier
#' @param qall_flavor Which quality criterion to use? By default qall_mult
#' @param shiftParents Logical, parental scores shifted?
#' @param compfile optional filename to write comparison info to
#' @param combscorefile Optional file to write comparison scores
#' @noRd
compareProbes.gp <- function (chk, 
                              scores, 
                              probe.suffix = c("P", "Q", "R"), 
                              fracdiff.threshold = 0.04,
                              qall_flavor = "qall_mult",
                              compfile, 
                              combscorefile){
  parent1 <- chk$meta$parent1 
  parent2 <- chk$meta$parent2
  F1 <- chk$meta$F1
  ancestors <- chk$meta$ancestors
  polysomic <- chk$meta$polysomi
  disomic <- chk$meta$disomic 
  mixed <- chk$meta$mixed
  ploidy <- chk$meta$ploidy
  ploidy2 <- ifelse(is.null(chk$meta$ploidy2),chk$meta$ploidy,chk$meta$ploidy2)
  shiftParents <- chk$meta$shiftParents
  chk <- chk$checked_F1
  
  noprobes <- length(probe.suffix) == 1 && is.na(probe.suffix)
  
  if (noprobes){
    funcname <- "writeScoresfile"
  } else{
    funcname <- "compareProbes"
  }
    
  if (!noprobes && (length(probe.suffix) != 3 || sum(is.na(probe.suffix)) >
                    0 || max(nchar(probe.suffix)) != min(nchar(probe.suffix)) ||
                    sum(abs(match(probe.suffix, probe.suffix) - c(1, 2, 3))) >
                    0))
    stop("compareProbes: probe.suffix invalid")
  
  if (ploidy%%2 != 0)
    stop(paste(funcname, ": odd parental ploidy not allowed",
               sep = ""))
  if (missing(ploidy2) || is.na(ploidy2))
    ploidy2 <- ploidy
  else if (ploidy2%%2 != 0)
    stop(paste(funcname, ": odd parental ploidy2 not allowed",
               sep = ""))
  ploidyF1 <- (ploidy + ploidy2)/2
  segtypecol <- which(tolower(names(chk)) == "segtype")
  if (length(segtypecol) != 1)
    segtypecol <- which(tolower(names(chk)) == "selfit")
  if (length(segtypecol) != 1)
    segtypecol <- which(tolower(names(chk)) == "bestparentfit")
  if (length(segtypecol) != 1)
    stop(paste(funcname, ": chk should have one column segtype or bestParentfit",
               sep = ""))
  if (class(chk[, segtypecol]) != "character")
    chk[, segtypecol] <- as.character(chk[, segtypecol])
  if (class(chk$MarkerName) != "character")
    chk$MarkerName <- as.character(chk$MarkerName)
  chk <- chk[!is.na(chk[, segtypecol]) & chk[, segtypecol] !=
               "", ]
  qallcol <- which(tolower(names(chk)) == qall_flavor)
  if (length(qallcol) > 1)
    qallcol <- qallcol[1]
  shiftspresent <- FALSE
  shiftcol <- which(tolower(names(chk)) == "shift")
  if (length(shiftcol) > 1) {
    stop(paste(funcname, ": chk should have at most one column 'shift'"),
         sep = "")
  }else if (length(shiftcol) == 0) {
    chk$shift <- rep(0, nrow(chk))
    shiftcol <- which(tolower(names(chk)) == "shift")
    shiftParents <- TRUE
  }else {
    if (missing(shiftParents) || is.na(shiftParents)) {
      if (ploidy == ploidy2)
        shiftParents <- TRUE
      else stop(paste(funcname, ": if ploidy2 != ploidy, shiftParents must be specified (FALSE/TRUE)",
                      sep = ""))
    }
    if (sum(is.na(chk[, shiftcol])) > 0)
      stop(paste(funcname, ": column 'shift' in chk may not contain missing values",
                 sep = ""))
    shiftspresent <- TRUE
  }
  progenitors <- c(parent1, parent2, ancestors)
  progenitors <- progenitors[!is.na(progenitors)]
  if (length(progenitors) == 0)
    shiftParents <- FALSE
  allsamp <- c(progenitors, F1)
  if (shiftParents)
    shfsamp <- rep(TRUE, length(allsamp))
  else shfsamp <- c(rep(FALSE, length(progenitors)), rep(TRUE, length(F1)))
  missamp <- allsamp[!(allsamp %in% unique(as.character(scores$SampleName)))]
  if (length(missamp) > 0)
    stop(paste(funcname, ": some samples not in scores:",
               paste(missamp, collapse = " ")))
  # missamp <- progenitors[!(progenitors %in% names(chk))]
  # if (length(missamp) > 0)
  #   stop(paste(funcname, ": some samples not in chk:", paste(missamp,
  #                                                            collapse = " ")))
  if (length(F1) != round(sum(chk[1, 5:(6 + ploidyF1)])))
    stop(paste(funcname, ": F1 has", length(F1), "samples but in chk",
               sum(chk[1, 5:(6 + ploidyF1)]), "values counted"))
  which_shf <- which(rightstr(chk$MarkerName, 4) == "_shf")
  chk$SNPname <- chk$MarkerName
  chk$SNPname[which_shf] <- leftstr(chk$MarkerName[which_shf],
                                    -4)
  if (sum(is.na(match(chk$SNPname, unique(as.character(scores$MarkerName))))) >
      0)
    stop(paste(funcname, ": not all markers from chk appear in scores",
               sep = ""))
  if (noprobes) {
    chk$probe <- rep("", nrow(chk))
  }else {
    mnl <- nchar(chk$SNPname)
    chk$probe <- substr(chk$SNPname, mnl - nchar(probe.suffix[1]) +
                          1, mnl)
    chk$SNPname <- substr(chk$SNPname, 1, mnl - nchar(probe.suffix[1]))
  }
  if (!noprobes && sum(!(chk$probe %in% probe.suffix[1:2])) >
      0)
    stop("compareProbes: not all marker names have a specified probe suffix")
  if (!noprobes && sum(nchar(chk$SNPname) == 0) > 0)
    stop("compareProbes: not all marker names have a specified probe suffix")
  if (!is.na(compfile) && !checkFilename(compfile))
    stop(paste(funcname, ": compfile not valid or not writable",
               sep = ""))
  if (!is.na(combscorefile) && !checkFilename(combscorefile))
    stop(paste(funcname, ": combscorefile not valid or not writable",
               sep = ""))
  if (is.factor(chk$parent1))
    chk$parent1 <- as.integer(as.character(chk$parent1))
  if (is.factor(chk$parent2))
    chk$parent2 <- as.integer(as.character(chk$parent2))
  if (is.factor(chk[, shiftcol]))
    chk[, shiftcol] <- as.integer(as.character(chk[, shiftcol]))
  if (is.factor(parent1))
    parent1 <- as.character(parent1)
  if (is.factor(parent2))
    parent2 <- as.character(parent2)
  if (is.factor(F1))
    F1 <- as.character(F1)
  if (is.factor(ancestors))
    ancestors <- as.character(ancestors)
  snpnames <- sort(unique(chk$SNPname))
  segtype <- array(NA, dim = c(length(snpnames), 2, 2), dimnames = list(snpname = snpnames,
                                                                        probe = c("probe1", "probe2"), shift = c("false", "true")))
  chkshift <- segtype
  qall <- segtype
  cons.parent1 <- segtype
  cons.parent2 <- segtype
  Rsegtype <- segtype
  for (sh in 1:2) {
    if (sh == 1)
      sel.sh <- chk$shift == 0
    else sel.sh <- chk$shift != 0
    if (noprobes)
      proberange <- 1
    else proberange <- 1:2
    for (prb in proberange) {
      if (noprobes) {
        prbdat <- chk[sel.sh, ]
      }
      else {
        prbdat <- chk[chk$probe == probe.suffix[prb] &
                        sel.sh, ]
      }
      rows <- match(prbdat$SNPname, snpnames)
      segtype[rows, prb, sh] <- prbdat[, segtypecol]
      chkshift[rows, prb, sh] <- prbdat[, shiftcol]
      if (length(qallcol) == 1)
        qall[rows, prb, sh] <- prbdat[, qallcol]
      cons.parent1[rows, prb, sh] <- prbdat$parent1
      cons.parent2[rows, prb, sh] <- prbdat$parent2
    }
  }
  seginfo <- selSegtypeInfo(calcSegtypeInfo(ploidy, ploidy2),
                            polysomic, disomic, mixed)
  if (sum(is.na(match(unique(chk[, segtypecol]), names(seginfo)))) >
      0)
    stop("compareProbes: chk contains segtypes not expected with the settings\nof polysomic/disomic/mixed/ploidy/ploidy2")
  scores <- scores[scores$SampleName %in% allsamp, ]
  combscores <- makeCombscoresDf(parent1, parent2, F1, ancestors,
                                 nrow = nrow(chk))
  combrow <- 0
  batchsize <- 100
  batchnr <- 1
  batchscores <- list()
  while (batchsize * (batchnr - 1) < length(snpnames)) {
    cat(paste("batch", batchnr, "\n"))
    minmrk <- batchsize * (batchnr - 1) + 1
    maxmrk <- min(length(snpnames), batchsize * batchnr)
    batchmarkers <- snpnames[minmrk:maxmrk]
    if (!noprobes) {
      batchmarkers <- c(paste(batchmarkers, probe.suffix[1],
                              sep = ""), paste(batchmarkers, probe.suffix[2],
                                               sep = ""))
    }
    batchscores <- scores[scores$MarkerName %in% batchmarkers,
                          ]
    for (mrk in minmrk:maxmrk) {
      sc <- list()
      progen.sc <- array(NA, dim = c(length(progenitors),
                                     2, 2), dimnames = list(sample = progenitors,
                                                            probe = c("probe1", "probe2"), shift = c("FALSE",
                                                                                                     "TRUE")))
      F1.sc <- array(NA, dim = c(length(F1), 2, 2), dimnames = list(sample = F1,
                                                                    probe = c("probe1", "probe2"), shift = c("FALSE",
                                                                                                             "TRUE")))
      parent.sc <- array(NA, dim = c(2, 2, 2), dimnames = list(sample = c("parent1",
                                                                          "parent2"), probe = c("probe1", "probe2"), shift = c("FALSE",
                                                                                                                               "TRUE")))
      P0col <- which(names(scores) == "P0")
      if (length(P0col) != 1)
        stop("compareProbes: scores should contain columns P0 .. P<ploidy>")
      for (prb in 1:2) if (sum(!is.na(segtype[mrk, prb,
                                              ])) > 0) {
        if (noprobes) {
          mrkname_scores <- snpnames[mrk]
        }
        else {
          mrkname_scores <- paste(snpnames[mrk], probe.suffix[prb],
                                  sep = "")
        }
        sc[[prb]] <- list()
        sc[[prb]][[1]] <- batchscores[batchscores$MarkerName ==
                                        mrkname_scores & batchscores$SampleName %in%
                                        allsamp, ]
        sc[[prb]][[1]] <- sc[[prb]][[1]][match(allsamp,
                                               sc[[prb]][[1]]$SampleName), ]
        if (!is.na(segtype[mrk, prb, 2])) {
          sc[[prb]][[2]] <- sc[[prb]][[1]]
          sv <- chkshift[mrk, prb, 2]
          sc[[prb]][[2]]$geno[shfsamp] <- sc[[prb]][[2]]$geno[shfsamp] +
            sv
          inval <- shfsamp & !(sc[[prb]][[2]]$geno %in%
                                 0:ploidyF1)
          sc[[prb]][[2]]$geno[inval] <- NA
          if (sv > 0) {
            sc[[prb]][[2]][shfsamp, P0col + ploidyF1] <- sum(sc[[prb]][[2]][shfsamp,
                                                                            P0col + ((ploidyF1 - sv):ploidyF1)])
            if (sv < ploidy)
              sc[[prb]][[2]][shfsamp, P0col + (sv:(ploidyF1 -
                                                     1))] <- sc[[prb]][[2]][shfsamp, P0col +
                                                                              (0:(ploidyF1 - 1 - sv))]
            sc[[prb]][[2]][shfsamp, P0col + (0:(sv -
                                                  1))] <- 0
          }
          else {
            sc[[prb]][[2]][shfsamp, P0col] <- sum(sc[[prb]][[2]][shfsamp,
                                                                 P0col + (0:(-sv))])
            if ((-sv) < ploidyF1)
              sc[[prb]][[2]][shfsamp, P0col + (1:(ploidyF1 +
                                                    sv))] <- sc[[prb]][[2]][shfsamp, P0col +
                                                                              ((1 - sv):ploidyF1)]
            sc[[prb]][[2]][shfsamp, P0col + (ploidyF1 +
                                               sv + 1):ploidyF1] <- 0
          }
          suppressWarnings(sc[[prb]][[2]]$maxP <- apply(sc[[prb]][[2]][P0col +
                                                                         (0:ploidyF1)], 1, max, na.rm = TRUE))
        }
        for (sh in 1:2) {
          if (!is.na(segtype[mrk, prb, sh])) {
            mrkname_out <- mrkname_scores
            if (shiftspresent)
              mrkname_out <- paste(mrkname_out, c("n",
                                                  "s")[sh], sep = "")
            progensc <- sc[[prb]][[sh]][sc[[prb]][[sh]]$SampleName %in%
                                          progenitors, ]
            progen.sc[match(progensc$SampleName, progenitors),
                      prb, sh] <- progensc$geno
            F1.sc[, prb, sh] <- sc[[prb]][[sh]]$geno[sc[[prb]][[sh]]$SampleName %in%
                                                       F1]
            parent.sc[, prb, sh] <- getPargeno(cons.parent1[mrk,
                                                            prb, sh], cons.parent2[mrk, prb, sh], segtype = segtype[mrk,
                                                                                                                    prb, sh], seginfo = seginfo)
            combrow <- combrow + 1
            if (combrow > nrow(combscores))
              combscores <- rbind(combscores, makeCombscoresDf(parent1,
                                                               parent2, F1, ancestors, nrow = 5000))
            combscores$MarkerName[combrow] <- mrkname_out
            combscores$segtype[combrow] <- segtype[mrk,
                                                   prb, sh]
            combscores[combrow, 3:length(combscores)] <- c(progen.sc[,
                                                                     prb, sh], parent.sc[, prb, sh], F1.sc[,
                                                                                                           prb, sh])
          }
        }
      }
      if (!noprobes) {
        numcomb <- 0
        for (sh1 in 1:2) for (sh2 in 1:2) if (!is.na(segtype[mrk,
                                                             1, sh1]) && !is.na(segtype[mrk, 2, sh2]) &&
                                              segtype[mrk, 1, sh1] == segtype[mrk, 2, sh2]) {
          bothscored <- sum(!is.na(sc[[1]][[sh1]]$geno) &
                              !is.na(sc[[2]][[sh2]]$geno))
          diff <- sc[[1]][[sh1]]$geno != sc[[2]][[sh2]]$geno
          Ndifferent <- sum(diff, na.rm = TRUE)
          fracdiff <- Ndifferent/bothscored
          if (!is.na(fracdiff) && fracdiff <= fracdiff.threshold) {
            Rsegtype[mrk, sh1, sh2] <- segtype[mrk, 1,
                                               sh1]
            numcomb <- numcomb + 1
            combgeno <- sc[[1]][[sh1]]$geno
            w <- which(is.na(sc[[1]][[sh1]]$geno) & !is.na(sc[[2]][[sh2]]$geno))
            combgeno[w] <- sc[[2]][[sh2]]$geno[w]
            w <- which(!is.na(sc[[1]][[sh1]]$geno) &
                         !is.na(sc[[2]][[sh2]]$geno) & sc[[2]][[sh2]]$geno !=
                         sc[[1]][[sh1]]$geno & sc[[2]][[sh2]]$maxP >
                         sc[[1]][[sh1]]$maxP)
            combgeno[w] <- sc[[2]][[sh2]]$geno[w]
            parent.na <- matrix(NA, nrow = 2, ncol = 2)
            parent.na[1, ] <- is.na(parent.sc[, 1, sh1])
            parent.na[2, ] <- is.na(parent.sc[, 2, sh2])
            combparent <- parent.sc[, 1, sh1]
            w <- which(parent.na[1, ] & !parent.na[2,
                                                   ])
            combparent[w] <- parent.sc[w, 2, sh2]
            w <- which(!parent.na[1, ] & !parent.na[2,
                                                    ] & parent.sc[, 1, sh1] != parent.sc[,
                                                                                         2, sh2])
            combparent[w] <- NA
            if (xor(parent.na[1, 1], parent.na[1, 2]) &&
                xor(parent.na[2, 1], parent.na[2, 2]) &&
                xor(parent.na[1, 1], parent.na[2, 1]) &&
                getMatchParents(combparent, seginfo[[segtype[mrk,
                                                             prb]]]) == "No") {
              combparent <- c(NA, NA)
            }
            if (shiftspresent) {
              mrkname_cmb <- paste(snpnames[mrk], "R",
                                   paste(c("n", "s")[c(sh1, sh2)], collapse = ""),
                                   sep = "")
            }
            else {
              mrkname_cmb <- paste(snpnames[mrk], "R",
                                   sep = "")
            }
            combrow <- combrow + 1
            if (combrow > nrow(combscores))
              combscores <- rbind(combscores, makeCombscoresDf(parent1,
                                                               parent2, F1, ancestors, nrow = 5000))
            combscores$MarkerName[combrow] <- mrkname_cmb
            combscores$segtype[combrow] <- segtype[mrk,
                                                   1, sh1]
            combscores[combrow, 3:length(combscores)] <- c(combgeno[1:length(progenitors)],
                                                           combparent, combgeno[(length(progenitors) +
                                                                                   1):length(combgeno)])
          }
        }
      }
    }
    batchnr <- batchnr + 1
  }
  if (combrow == 0)
    combscores <- combscores[0, ]
  else combscores <- combscores[1:combrow, ]
  if (!is.na(combscorefile))
    write.table(combscores, combscorefile, quote = FALSE,
                sep = "\t", na = "", row.names = FALSE, col.names = TRUE)
  if (noprobes) {
    compstat <- NA
  }else {
    compstat <- data.frame(SNPname = snpnames)
    sufx <- matrix(c(paste(probe.suffix[1], "n", sep = ""),
                     paste(probe.suffix[1], "s", sep = ""), paste(probe.suffix[2],
                                                                  "n", sep = ""), paste(probe.suffix[2], "s", sep = "")),
                   2, 2, byrow = TRUE)
    max_sh <- 2
    Rsufx <- matrix(c("nn", "ns", "sn", "ss"), 2, 2, byrow = TRUE)
    if (!shiftspresent) {
      max_sh <- 1
      sufx <- substr(sufx, 1, 1)
      Rsufx <- matrix("")
    }
    for (prb in 1:2) for (sh in 1:max_sh) {
      compstat[[paste("segtype", sufx[prb, sh], sep = "")]] <- segtype[,
                                                                       prb, sh]
      if (length(qallcol) == 1) {
        compstat[[paste("qall", sufx[prb, sh], sep = "")]] <- qall[,
                                                                   prb, sh]
      }
    }
    for (sh1 in 1:max_sh) for (sh2 in 1:max_sh) {
      compstat[[paste("segtype", probe.suffix[3], Rsufx[sh1,
                                                        sh2], sep = "")]] <- Rsegtype[, sh1, sh2]
    }
    for (prb in 1:2) {
      compstat[[paste("count", probe.suffix[prb], sep = "")]] <- rowSums(!is.na(segtype[,
                                                                                        prb, ]))
    }
    compstat[[paste("count", probe.suffix[3], sep = "")]] <- apply(!is.na(Rsegtype),
                                                                   1, sum)
    rownames(compstat) <- NULL
    if (!is.na(compfile))
      write.table(compstat, compfile, row.names = FALSE,
                  col.names = TRUE, sep = "\t", quote = FALSE,
                  na = "")
  }
  return(invisible(list(combscores = combscores, compstat = compstat)))
} #compareProbes.gp

#' User interface for specifying linkage group for a homologue
#' @param LG_hom_stack A data.frame defining clusters
#' @param LG_number Number of chromosomes (linkage groups)
#' @noRd
createChmHomList <- function(LG_hom_stack, LG_number=12){
  ## This is intended to allow the user to define which cluster elements should be
  ## put together into a single file. The user defines how many chromosomes are identified
  ## in the data:
  homolog_list <- sapply(1:LG_number, function(x) NULL)
  
  for(c in 1:LG_number){
    cat(paste("Linkage group ",c," :\n"))
    cat("_______________________________________\n")
    
    cluster_tracker <- NULL #To make sure unique clusters are entered
    
    temp_hom <- readline("Enter first cluster number:   ")
    temp_group <- NULL
    
    while(temp_hom != "x") {
      if(temp_hom %in% cluster_tracker){
        cat("Error, this cluster was already assigned\n")
        temp_hom <- readline("Please re-enter next cluster number (or x to exit):   ")
      } else if(temp_hom %in% names(LG_hom_stack) == FALSE){
        cat("Error, impossible cluster assignment\n")
        temp_hom <- readline("Please re-enter next cluster number (or x to exit):   ")
      }else{
        cluster_tracker <- c(cluster_tracker,temp_hom)
        temp_group <- c(temp_group,temp_hom)
        temp_hom <- readline("Enter next cluster number (or x to exit):   ")
      }
    }
    homolog_list[[c]] <- temp_group
    cat("_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n") 
  }
  
  cat("Clusters defined as: \n")
  print(homolog_list)
  
  user_abort <- "n"
  
  suppressWarnings(
    if(length(setdiff(1:length(LG_hom_stack),as.numeric(unlist(homolog_list)))) != 0){
      cat("The following cluster(s) have not been assigned to any linkage group:\n")
      print(setdiff(1:length(LG_hom_stack),as.numeric(unlist(homolog_list))))
      user_abort <- readline("Do you wish to abort and try again? (y/n)   ")
    })
  
  if(user_abort != "y"){
    
    ## Now we should go and collect the marker names from each cluster, add chromosome numbers
    names(homolog_list)<-as.character(1:length(homolog_list))
    homolog_m<-stack(homolog_list)
    colnames(homolog_m)<-c("homologue", "LG")
    marker_cluster<-stack(LG_hom_stack)
    colnames(marker_cluster)<-c("SxN_Marker", "homologue")
    output_df<-merge(marker_cluster, homolog_m, by="homologue")
    output_df<-output_df[,c("SxN_Marker","LG","homologue")]
    
    } else{
    cat("Attempt to define linkage groups has been aborted. Please re-run function.\n")
  }
  return(output_df)
} #createChmHomList



#'digamfrq ***********************************************************
#'@description Get the integer ratios of gametes with dosage 0:(ploidy/2)
#'from an allopolyploid parent (disomic inheritance) with given ploidy and dosage
#'@param ploidy The ploidy level of the parent
#'@param dosage The parental dosage
#'@noRd
digamfrq <- function(ploidy, dosage) {
  gp <- ploidy/2 #gamete ploidy
  #first step: calculate all possible compositions of diploid parental genomes
  #i.e. the numbers of subgenomes with nulliplex, simplex and duplex dosages
  maxdup <- dosage %/% 2 #max nr of subgenomes with duplex dosage
  mindup <- max(dosage-gp, 0) #min nr of subgenomes with duplex dosage
  genomes <- matrix(integer((maxdup-mindup+1)*3), ncol=3)
  colnames(genomes) <- c("nulli", "sim", "du")
  for (i in 1:nrow(genomes)) {
    genomes[i,3] <- mindup + i - 1
    genomes[i,2] <- dosage - 2*genomes[i,3]
    genomes[i,1] <- gp - sum(genomes[i, 2:3])
  }
  #next step: per genomic composition calculate the gamete distributions:
  gam <- matrix(numeric(nrow(genomes)*(gp+1)), nrow=nrow(genomes)) #all 0
  colnames(gam) <- 0:gp
  for (i in 1:nrow(genomes)) {
    mingamdos <- genomes[i,3] #each duplex genome contributes 1
    maxgamdos <- gp-genomes[i,1] #each nulliplex genome contributes 0
    bincoeff <- choose(genomes[i,2], 0:genomes[i,2])
    
    # to return fractions:
    #gam[i, (mingamdos+1):(maxgamdos+1)] <- bincoeff / sum(bincoeff)
    
    #to return integer ratios:
    gam[i, (mingamdos+1):(maxgamdos+1)] <- bincoeff
    #the smallest bincoeff is always 1 so division by the gcd is not needed
  }
  #list(genomes, gam)
  gam
} #digamfrq



#'gcd ***********************************************************
#'@description  Return the greatest common divisor of two integers (vectorized: x and y may be vectors)
#'@param x integer (or vector of integers)
#'@param y integer (or vector of integers)
#'@noRd
gcd <- function(x,y) {
  r <- x %% y;
  return(ifelse(r, gcd(y, r), y))
} #gcd



#'gcd_all ***********************************************************
#'@description Return the greatest common divisor of all elements of vector x
#'@param x Input vector
#'@noRd
gcd_all <- function(x) {
  if (length(x) == 1) x else {
    if (length(x) == 2) gcd(x[1], x[2]) else
      gcd_all(c(gcd(x[1], x[2]), x[3:length(x)]))
  }
} #gcd_all



#'getConsensusGeno ***********************************************************
#'@description Get consensus genotype
#'@param geno A vector with scores (0..ploidy or NA) for all samples of one parent.
#'@param maxNAfrac The maximum fraction of missing scores allowed to assign a high-confidence 
#'consensus geno score; if more are missing a missing geno score is returned (the default requires 
#'MORE than half of the samples to be scored, so one missing out of 2 samples already causes a missing geno)
#'@param lowconf.NAfrac if the fraction missing scores is more than \code{maxNAfrac} but less than \code{lowconf.NAfrac}, 
#'or if there is only one sample, a \code{lowconf.geno} score is assigned
#'@return Returns a list with 4 components:
#'\item{geno}{ 
#'The consensus of the dosages in vector geno. NA if geno is empty,
#'if all elements of geno are NA, if there are different non-NA
#'values in geno, or if the fraction of NA values in geno larger than \code{maxNAfrac}
#'}
#'\item{lowconf.geno}{
#'If there is no conflict and the fraction of missing scores
#'is between \code{maxNAfrac} and \code{lowconf.NAfrac}, the consensus is assigned
#'as a "low-confidence" geno - this needs confirmation from the
#'F1 segregation to be accepted}
#'\item{conflict}{ 
#'\code{TRUE} if there are different non-NA values in geno, else \code{FALSE}
#'}
#'\item{NAfrac}{ 
#'The fraction NA in geno (0.5 if length(geno)==0)
#'}
#'@noRd
getConsensusGeno <- function(geno, 
                             maxNAfrac=0.499,
                             lowconf.NAfrac=0.751) {
  nwgeno <- NA; nwlogeno <- NA; conflict <- FALSE; NAfrac <- 0
  if (length(geno) == 0) NAfrac <- 0.5 else {   #we need a NAfrac value for q2
    NAfrac <- sum(is.na(geno) / length(geno))
    if (NAfrac < 1.0) {
      if (length(geno) == 1) nwlogeno <- geno else {
        if (max(geno, na.rm=TRUE) != min(geno, na.rm=TRUE)) {
          conflict <- TRUE
        } else  if (NAfrac <= maxNAfrac) {
          nwgeno <- min(geno, na.rm=TRUE)
        } else if (NAfrac <= lowconf.NAfrac) {
          nwlogeno <- min(geno, na.rm=TRUE)
        }
      }
    }
  }
  list(geno=nwgeno, lowconf.geno=nwlogeno, conflict=conflict, NAfrac=NAfrac)
} #getConsensusGeno



#'getMatchParents ***********************************************************
#'@description Does the specified segtype match the parental genotypes?
#'@param parGeno integer vector of length 2 with the two parental (consensus) dosages (can each be NA)
#'@param seginfoItem one of the items (segtypes) of a list as returned by \code{\link{calcSegtypeInfo}}
#'@noRd
getMatchParents <- function(parGeno, 
                            seginfoItem) {
  p <- which(!is.na(parGeno))
  if (length(p) == 0) {
    #both parental genotypes unknown
    return("Unknown")
  } else if (length(p) == 1) {
    #only one of the parental genotypes known
    if (sum(seginfoItem$pardosage[, p] == parGeno[p]) > 0)
      return("OneOK") else return("No")
  } else {
    #both parental genotypes known
    i <- which(seginfoItem$pardosage[, 1] == parGeno[1])
    if (length(i) == 1 && parGeno[2] == seginfoItem$pardosage[i, 2])
      return("Yes") else return("No")
  }
} #getMatchParents

#' getPargeno
#' @description Get parental genotypes
#' @param P1consensus the consensus scores of parent1
#' @param P2consensus the consensus scores of parent2
#' @param segtype the name of the selected segtype (selfit or (in the new polyploid version) bestParentfit from \code{checkF1};
#' should occur in names(seginfo); OR the number of the segtype in seginfo
#' @param matchParents "Yes" if both consensus parental dosages match segtype
#' (bestParentfit.matchParents from \code{checkF1})
#' @param seginfo a list as returned by \code{calcSegtypeInfo}, with the par.dosages
#' limited to those fitting the auto, allo and mixed parameters used to calculate \code{checkF1}
#' @param allparentscores if FALSE (default) and more than one combination would be
#' possible, c(NA, NA) is returned. If TRUE then a matrix with
#' one row for each possibility and two columns for the 2 parental
#' dosages is returned.
#' @return an integer vector of two elements: the inferred dosages
#' of parent 1 and 2; or (if allparentscores is TRUE and there are
#' multiple combinations) a 2-column matrix with one row per parental combination
#' @noRd
getPargeno <- function(P1consensus,
                       P2consensus,
                       segtype, 
                       matchParents=NULL, 
                       seginfo,
                       allparentscores=FALSE) {

  # we try to fill in the expected parental dosage
  # based on the segtype and the consensus of the observed parental scores:
  parent.sc <- c(P1consensus, P2consensus)
  if (is.character(segtype)) {
    segtype <- which(names(seginfo) == segtype)
    if (length(segtype) != 1) stop("internal error in getPargeno")
  } #else segtype is already the number of the segtype
  exp.par <- seginfo[[segtype]]$pardosage #already selected to match
  #                                           auto/allo/mixed
  # now get the final parental dosage from consensus and segtype:
  if (nrow(exp.par) == 1) {
    # both parents must have the same genotype so we can fill these in,
    # irrespective of the observed parental genotype and the matchParent status:
    parent.sc <- exp.par[1,]
  } else {
    if (is.null(matchParents)) {
      matchParents <- getMatchParents(parGeno=parent.sc,
                                      seginfoItem=seginfo[[segtype]])
    }
    if (matchParents != "Yes") {
      # there is more than one possible parental combination.
      # if matchParents=="Yes" we don't need to do anything;
      # else (matchParents= Unknown, OneOK or No) we can only fill in the parents
      # if one parent matches one combination and the other matches none
      pmatch <- rep(NA, 2)
      for (p in 1:2) pmatch[p] <- match(parent.sc[p], exp.par[, p])
      if (sum(!is.na(pmatch)) == 1) {
        #one parent matches one expected parental combination and the other
        #doesn't match any
        r <- min(pmatch, na.rm=TRUE)
        parent.sc <- exp.par[r,]
      }
      else {
        # two possibilities:
        # - both parents each match a different row of exp.par
        # - none of the parents matches a row of exp.par
        # (if both parents would match the same row, matchParents would be Yes)
        # In both cases we cannot assign the parental genotypes:
        if (allparentscores) {
          parent.sc <- exp.par #a 2-column matrix with 2 or more rows
        } else parent.sc <- rep(NA, 2)
      }
    } #matchParents!="Yes"
  } #nrow(exp.par) == 1 else
  parent.sc
} #getPargeno



#'@title Get substrings from the lefthand side
#'@description Get substrings from the lefthand side
#'@usage leftstr(x, n)
#'@param x a character vector (or something having an as.character method)
#'@param n a single  number: if n>=0, the leftmost n characters of each element
#'of x are selected, if n<0 the (-n) rightmost characters are omitted
#'@return character vector with leftside substrings of x
#'@noRd
leftstr <- function(x, n) {
  #vectorized for x, not n
  #n>=0: take the leftmost n characters
  #n<0: take all but the rightmost (-n) characters
  if (n >= 0) {
    substr(x, 1, n) #automatically converts x to character if needed
  } else {
    #n < 0
    if (!is.character(x)) x <- as.character(x)
    len <- nchar(x)
    substr(x, 1, len + n)
  }
} #leftstr



#'linterpol ***********************************************************
#'@description linear interpolation: returns the y values matching the x values (vector)
#'on the line through points pnt1 and pnt2 (both are vectors of length 2)
#'possible error if both X-coordinates equal not checked or caught
#'@param x x-coordinate for which corresponding y-coordinate is sought. 
#'@param pnt1 Point 1
#'@param pnt2 Point 2
#'@noRd
linterpol <- function(x, pnt1, pnt2) {
  a <- (pnt2[2] - pnt1[2]) / (pnt2[1] - pnt1[1])
  b <- pnt1[2] - a * pnt1[1]
  a * x + b
} #linterpol


#' makeCombscoresDf
#' @description Create an empty data frame with the correct size and column names and types
#' @param parent1 character vector with sampleIDs for parent 1. (0, 1 or multiple samples per 
#' parent allowed, 0 samples are specified as character(0))
#' @param parent2 character vector with sampleIDs for parent 2. (0, 1 or multiple samples per 
#' parent allowed, 0 samples are specified as character(0)
#' @param F1 character vector with sampleIDs for F1 individuals (one sample per F1 individual)
#' @param ancestors character vector of other ancestor samples listed in the \code{checkF1} output; character(0) if none
#' @param other character vector with sampleIDs of any other samples; character(0) if none
#' @param nrow the number of rows to create
#' @noRd
makeCombscoresDf <- function(parent1, parent2, F1, ancestors=character(0),
                             other=character(0), nrow){
  
  ncol <- length(parent1) + length(parent2) + length(ancestors) + 2 +
    length(F1) + length(other)
  scores <- matrix(rep(NA_integer_, nrow * ncol), nrow=nrow)
  df <- data.frame(
    MarkerName=rep("", nrow),
    segtype=rep("", nrow),
    scores,
    stringsAsFactors=FALSE
  )
  names(df) <- c("MarkerName", "segtype", parent1, parent2, ancestors,
                 "parent1", "parent2", F1, other)
  return(df)
} #makeCombscoresDf


#'makeProgeny ***********************************************************
#'@description Returns a vector with the integer ratios of the F1 progeny
#'@param gametes1 vector of the integer ratios of the gametes produced by parent 1
#'@param gametes2 vector of the integer ratios of the gametes produced by parent 2
#'@noRd
makeProgeny <- function(gametes1, gametes2) {
  gp1 <- length(gametes1)-1 #gametes1 ploidy
  gp2 <- length(gametes2)-1 #gametes2 ploidy
  m <- matrix(0, nrow=gp1+gp2+1, ncol=gp2+1)
  for (i in 1:(gp2+1))
    m[i:(i+gp1), i] <- gametes1 * gametes2[i]
  D <- round(rowSums(m)) #the integer ratios of the F1 generation
  names(D) <- 0:(length(D)-1)
  d <- gcd_all(D[D>0])
  round(D/d) #the integer ratios of the F1 generation simplified
} #makeProgeny


merge_marker_assignments.gp <- function(MarkerType,
                                        target_parent = "P1",
                                        other_parent = "P2",
                                        LG_hom_stack,
                                        SN_linked_markers,
                                        ploidy,
                                        LG_number,
                                        log = NULL) {
  LG_hom_stack <- test_LG_hom_stack(LG_hom_stack)
  markers_LG_hom_stack <- as.character(LG_hom_stack[, "SxN_Marker"])
  LG_hom_stack <- LG_hom_stack[, c("LG", "homologue")]
  rownames(LG_hom_stack) <- markers_LG_hom_stack
  colnames(LG_hom_stack) <- c("Assigned_LG", "Assigned_Homolog")
  comb <- as.matrix(do.call(rbind, SN_linked_markers))
  
  #Add SxN markers
  SN_LG_mat <- t(matrix(sapply(LG_hom_stack[, "Assigned_LG"], function(x) {
    a <- rep(0, LG_number)
    a[x] <- 1
    return(a)
  }),nrow = LG_number, dimnames = list(paste0("LG", 1:LG_number),markers_LG_hom_stack)
  ))

  LG_mat <- rbind(SN_LG_mat, comb[, paste0("LG", 1:LG_number), drop = FALSE])
  LG_mat <-
    cbind(c(LG_hom_stack[, "Assigned_LG"], comb[, "Assigned_LG"]), LG_mat)
  rownames(LG_mat) <- c(markers_LG_hom_stack, rownames(comb))
  colnames(LG_mat)[1] <- "Assigned_LG"
  
  SN_hom_mat <- t(matrix(sapply(LG_hom_stack[, "Assigned_Homolog"], function(x) {
    a <- rep(0, ploidy)
    a[x] <- 1
    return(a)
  }),nrow = ploidy, dimnames = list(paste0("Hom", 1:ploidy),markers_LG_hom_stack)
  ))

  counts_hom_mat <- rbind(SN_hom_mat, comb[, paste0("Hom", 1:ploidy)])
  
  assigned_hom_mat <- matrix(c(LG_hom_stack[, "Assigned_Homolog"],
                               rep(NA, (ploidy - 1) * nrow(LG_hom_stack))),
                             ncol = ploidy)

  ##incorporate number of linkages per homologue
  assigned_hom_mat <-
    rbind(assigned_hom_mat, comb[, paste0("Assigned_hom", 1:ploidy)])
  Assigned_LG_hom <- cbind(LG_mat, counts_hom_mat, assigned_hom_mat)

  if(target_parent == "P1"){
    colNme <- c("parent1","parent2")
  }else{
    colNme <- c("parent2","parent1")
  }
  parental_dosages <- MarkerType[match(rownames(Assigned_LG_hom),MarkerType$MarkerName),colNme]

  colnames(parental_dosages) <- c(target_parent, other_parent)
  rownames(parental_dosages) <- rownames(Assigned_LG_hom)
  Assigned_LG_hom <-
    as.matrix(cbind(parental_dosages, Assigned_LG_hom))
  class(Assigned_LG_hom) <- "integer"
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  write(paste0("\n####Marker numbers  parent ", target_parent, "\n"),
        log.conn)
  count_table <-
    table(
      paste0(Assigned_LG_hom[, target_parent], "x", Assigned_LG_hom[, other_parent]),
      Assigned_LG_hom[, "Assigned_LG"],
      dnn = list("markertype", "linkage group")
    )
  
  write(knitr::kable(count_table),
        log.conn)
  
  if (!is.null(log))
    close(log.conn)
  
  return(Assigned_LG_hom)
} #merge_marker_assignments.gp


#' Calculate the mode of a vector
#' @param x Vector input, can be numeric, character or factor
#' @noRd
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
} #Mode



#' Plot links between 1.0 markers
#' @description Make a plot with recombination frequency versus LOD score.
#' @param linkage_df A linkage data.frame.
#' @param LG_hom_stack A \code{data.frame} as a result of \code{\link{bridgeHomologues}}
#' @param LG Linkage group (LG) number.
#' @param h1 Homologue to be compared.
#' @param h2 Homologue to be compared against h1
#' @param ymax Maximum limit of y-axis
#' @noRd
plot_SNlinks <- function(linkage_df,LG_hom_stack,
                         LG, h1, h2, ymax=NULL){
  
  h1Markers <- LG_hom_stack[LG_hom_stack[,"LG"] == LG & LG_hom_stack[,"homologue"] == h1,"SxN_Marker"]
  # find markers in cluster to combine
  h2Markers <- LG_hom_stack[LG_hom_stack[,"LG"] == LG & LG_hom_stack[,"homologue"] == h2,"SxN_Marker"]
  
  # find markers that are in both clusters 
  subdata1 <- linkage_df[linkage_df[,"marker_a"] %in% h1Markers & linkage_df[,"marker_b"] %in% h2Markers,]
  subdata2 <- linkage_df[linkage_df[,"marker_a"] %in% h2Markers & linkage_df[,"marker_b"] %in% h1Markers,]
  
  plotdata <- rbind(subdata1,subdata2)
  
  if(is.null(ymax)){ ymax<- max(plotdata$LOD)}
  
  if(nrow(plotdata) > 0){
    
    plot_linkage_df(linkage_df = plotdata, #there was droplevels() here..removed July 2020
                    main=paste0("LG",LG," h",h1," & h",h2), 
                    add_legend = F)
    
    # plot(NULL,xlab="r",ylab="LOD",
    #      xlim=c(0,0.5),ylim=c(0,ymax),
    #      main=paste0("LG",LG," h",h1," & h",h2))
    # with(plotdata[plotdata$phase=="coupling",],
    #      points(r,LOD,pch=19,col="limegreen"))
    # with(plotdata[plotdata$phase=="repulsion",],
    #      points(r,LOD,pch=19,col="dodgerblue"))
  } else{
    warning(paste("No SxN linkage found between h",h1," and h",h2," (LOD too low?)",sep=""))
  }
} #plot_SNlinks



#'polygamfrq ***********************************************************
#'@description Get the integer ratios of gametes with dosage 0:(ploidy/2)
#'from an autopolyploid parent (polysomic inheritance) with given ploidy and dosage
#'the result can be converted to fractions by y <- x/(sum(x))
#'@param ploidy The ploidy level of the parent
#'@param dosage The parental dosage
#'@noRd
polygamfrq <- function(ploidy, dosage) {
  gp <- ploidy/2 #gamete ploidy
  
  #using hypergeometric distribution:
  #dhyper(0:gp, dosage, ploidy-dosage, gp)
  
  #using n-over-k: function choose, in order to return integer ratios:
  A <- choose(dosage, 0:gp)
  B <- choose(ploidy-dosage, gp-(0:gp))
  #C <- choose(ploidy, gp)
  #A * B / C this would give the answer as fractions,
  # like dhyper(0:gp, dosage, ploidy-dosage, gp)
  D <- round(A*B)
  d <- gcd_all(D[D>0])
  round(D/d)
} #polygamfrq



#' Prepare pwd for mapping scripts e.g. MDSmap
#' @description Prepare a dataframe with pairwise recombination frequency and LOD score for mapping by re-ordering
#' @param pwd A pwd data.frame
#' @noRd
prepare_pwd <- function(pwd) {
  pwd[,c("marker_a", "marker_b")]  <-
    t(apply(pwd[,c("marker_a", "marker_b")], 1, function(x)
      x[order(x)]))
  
  switched.pwd<-rbind(pwd[,c("marker_a","marker_b")],
                      data.frame("marker_a"=pwd$marker_b,"marker_b"=pwd$marker_a))
  dupes<-which(duplicated(switched.pwd,fromLast = T))
  if(any(dupes <= nrow(pwd))) pwd <- pwd[-dupes[dupes <= nrow(pwd)],]
  pwd <- pwd[order(pwd$marker_a, pwd$marker_b),]
  pwd <- pwd[!duplicated(paste0(pwd$marker_a, pwd$marker_b)),]
  return(pwd[,c("marker_a", "marker_b", "r", "LOD")])
} #prepare_pwd



#'@title Get substrings from the righthand side
#'@description Get substrings from the righthand side
#'@usage rightstr(x, n)
#'@param x a character vector (or something having an as.character method)
#'@param n a single  number: if n>=0, the rightmost n characters of each element
#'of x are selected, if n<0 the (-n) leftmost characters are omitted
#'@return character vector with rightside substrings of x
#'@noRd
rightstr <- function(x, n) {
  #vectorized for x, not n
  #n>=0: take the rightmost n characters
  #n<0: take all but the leftmost (-n) characters
  if (!is.character(x)) x <- as.character(x)
  len <- nchar(x)
  if (n >= 0) {
    start <- len - n + 1
    substr(x, start, len)
  } else {
    #n < 0
    substr(x, -n + 1, len)
  }
} #rightstr



# segtypeInfoSummary *********************************************************
#'@title Summarize the segtypeInfo list
#'@description From a list of segregation types as produced by calcSegtypeInfo
#'or selSegtypeInfo, produce a data frame that only lists the parental
#'dosage combinations for each segtype and whether these produce the
#'segtype under polysomic, disomic and/or mixed inheritance.
#'Useful to quickly look up which segtypes match a given parental dosage
#'combination.
#'@usage segtypeInfoSummary(segtypeInfo)
#'
#'@param segtypeInfo a list as returned by calcSegtypeInfo or selSegtypeInfo
#'@return A data frame summarizing the segtypeInfo list, with columns:
#'\itemize{
#'\item{segtype: the name of the segregation type (see details of
#'calcSegtypeInfo)}
#'\item{segtypenr: the sequential number of the segtype in parameter segtypeInfo}
#'\item{parent1, parent2: dosages of the two parents}
#'\item{par.poly, par.di, par.mixed: whether these parental dosages produce this
#'segtype under polysomic, disomic and/or mixed inheritance}
#\item{segtype.poly, segtype.di, segtype.mixed: whether this segtype does occur
#under polysomic, disomic and/or mixed inheritance}
#'}
#'@noRd
segtypeInfoSummary <- function(segtypeInfo) {
  totrows <- 0
  for (i in 1:length(segtypeInfo))
    totrows <- totrows + nrow(segtypeInfo[[i]]$pardosage)
  result <- data.frame (segtype=character(totrows),
                        segtypenr=integer(totrows),
                        parent1=integer(totrows),
                        parent2=integer(totrows),
                        par.poly=integer(totrows),
                        par.di=integer(totrows),
                        par.mixed=integer(totrows),
                        #segtype.poly=integer(totrows),
                        #segtype.di=integer(totrows),
                        #segtype.mixed=integer(totrows),
                        stringsAsFactors=FALSE)
  r <- 1
  for (i in 1:length(segtypeInfo)) {
    for (p in 1:nrow(segtypeInfo[[i]]$pardosage)) {
      result$segtype[r] = names(segtypeInfo)[i]
      result$segtypenr[r] = i
      result$parent1[r] = segtypeInfo[[i]]$pardosage[p,1]
      result$parent2[r] = segtypeInfo[[i]]$pardosage[p,2]
      result$par.poly[r] = segtypeInfo[[i]]$parmode[p,1]
      result$par.di[r] = segtypeInfo[[i]]$parmode[p,2]
      result$par.mixed[r] = segtypeInfo[[i]]$parmode[p,3]
      #result$segtype.poly[r] = 0 + segtypeInfo[[i]]$polysomic
      #result$segtype.di[r] = 0 + segtypeInfo[[i]]$disomic
      #result$segtype.mixed[r] = 0 + segtypeInfo[[i]]$mixed
      r <- r + 1
    }
  }
  result
} #segtypeInfoSummary


# selSegtypeInfo ***********************************************************
#'@title Restrict a list of segregation types to specified inheritance modes
#'@description From a list of segregation types as produced by calcSegtypeInfo,
#'this function selects only those segtypes that occur with polysomic,
#'disomic and/or mixed inheritance if the corresponding parameters are set to
#'TRUE, and from those segtypes only the parental dosages with the same
#'restrictions are retained.
#'@usage selSegtypeInfo(segtypeInfo, polysomic, disomic, mixed)
#'@param segtypeInfo a list as returned by calcSegtypeInfo
#'@param polysomic If TRUE all segtypes with poly TRUE are retained, and from
#'those segtypes all parental dosage combinations with parmode[,1] TRUE
#'@param disomic If TRUE all segtypes with di TRUE are retained, and from
#'those segtypes all parental dosage combinations with parmode[,2] TRUE
#'@param mixed If TRUE all segtypes with mixed TRUE are retained, and from
#'those segtypes all parental dosage combinations with parmode[,3] TRUE
#'@return A list like segtypeInfo, modified as specified by parameters polysomic,
#'disomic and mixed
#'@noRd
selSegtypeInfo <- function(segtypeInfo, polysomic, disomic, mixed) {
  result <- list()
  selected <- integer(0)
  for (s in seq_along(segtypeInfo)) {
    if ((polysomic && segtypeInfo[[s]]$polysomic) ||
        (disomic && segtypeInfo[[s]]$disomic) ||
        (mixed && segtypeInfo[[s]]$mixed)) {
      selected <- c(selected, s)
      r <- length(result) + 1
      result[[r]] <- segtypeInfo[[s]]
      selpar <- (polysomic & result[[r]]$parmode[,1]) |
        (disomic & result[[r]]$parmode[,2]) |
        (mixed & result[[r]]$parmode[,3])
      result[[r]]$pardosage <-  result[[r]]$pardosage[selpar,, drop=FALSE]
      result[[r]]$parmode <-  result[[r]]$parmode[selpar,, drop=FALSE]
    }
  }
  names(result) <- names(segtypeInfo)[selected]
  result
} #selSegtypeInfo



#' Error and warning handling for cluster_stack 
#' @param cluster_stack A data.frame with a column "marker" specifying markernames, 
#' and a column "cluster" specifying marker cluster.
#' @noRd
test_cluster_stack <- function(cluster_stack){
  cn <- colnames(cluster_stack)
  cnw <- c("marker", "cluster")
  if(!all(cnw %in% cn)){
    warning("The names \"marker\" and \"cluster\" should be part of the columnames of cluster_stack")
    message("Trying to change columnames..")
    colnames(cluster_stack) <- cnw 
  }
  return(cluster_stack)
} #test_cluster_stack



#' Error and warning handling for dosage_matrix 
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @noRd
test_dosage_matrix <- function(dosage_matrix){
  if(inherits(dosage_matrix,"data.frame")){
    warning("dosage_matrix should be a matrix, now it's a data.frame.")
    message("Trying to convert it to matrix, assuming markernames are in the first column..")
    rownames(dosage_matrix) <- dosage_matrix[,1]
    dosage_matrix <- as.matrix(dosage_matrix[,-1])
    class(dosage_matrix) <- "integer"
  } else if(inherits(dosage_matrix,"matrix")){
    rn <- rownames(dosage_matrix)
    cn <- colnames(dosage_matrix)
    if(is.null(rn)) stop("The rownames of dosage_matrix should contain markernames. Now NULL")
    if(is.null(cn)) stop("The columnnames of dosage_matrix should contain genotype names. Now NULL")
    if(!(typeof(dosage_matrix) =="integer" | typeof(dosage_matrix) =="double")){
      warning("dosage_matrix should be integer or numeric. Trying to convert it.. ")
      class(dosage_matrix) <- "integer"
    }
  } else {
    stop("dosage_matrix should be a matrix of integers. 
         See the manual of this function for more information.")
  }
  if(any(duplicated(rownames(dosage_matrix)))){
    warning("Duplicated marker names detected. Removing duplicates as quick-fix.. (but advise to re-check your input also!)")
    dosage_matrix <- dosage_matrix[!duplicated(dosage_matrix),]
    }
  return(dosage_matrix)
} #test_dosage_matrix



#' Error and warning handling for LG_hom_stack 
#' @param LG_hom_stack A data.frame with a column "SxN_Marker" specifying markernames, 
#' a column "homologue" specifying homologue cluster and "LG" specifying linkage group.
#' @noRd
test_LG_hom_stack <- function(LG_hom_stack){
  cn <- colnames(LG_hom_stack)
  cnw <- c("SxN_Marker", "homologue", "LG")
  if(!all(cnw %in% cn)){
    warning("The names \"SxN_Marker\", \"homologue\" and \"LG\" should be part of the columnames of LG_hom_stack")
    message("Trying to change columnames..")
    colnames(LG_hom_stack) <- cnw 
  }
  if(!is.factor(LG_hom_stack$LG)) LG_hom_stack$LG <- as.factor(LG_hom_stack$LG)
  return(LG_hom_stack)
} #test_LG_hom_stack



#' Error and warning handling for linkage_df 
#' @param linkage_df A linkage data.frame as output of \code{\link{linkage}} calculating linkage between 1.0 markers.
#' @noRd
test_linkage_df <- function(linkage_df){
  if(!is.data.frame(linkage_df)){ 
    stop("linkage_df should be an object of class data.frame")
  } else {
    cn <- colnames(linkage_df)
    cnw <- c("marker_a","marker_b","r","LOD","phase")
    if(!identical(cn[1:5], cnw)){
      warning(
        "The first five columnnames of linkage_df are not of format c(\"marker_a\",\"marker_b\",\"r\",\"LOD\",\"phase\")")
      message("Trying to convert colnames..")
      colnames(linkage_df)[1:5] <- cnw
    }
    
    linkage_df[,c("marker_a", "marker_b")]<- 
      sapply(linkage_df[,c("marker_a", "marker_b")], as.character)
    linkage_df$phase <- as.factor(linkage_df$phase)
    classes <- sapply(linkage_df, class)
    names(classes) <- NULL
    classesw <- c("character", "character", "numeric", "numeric", "factor")
    if(!identical(classes[1:5], classesw)){
      stop(paste("The first five columns of linkage_df should be of classes:
                 ", paste(classesw, collapse = " ")))
    }
  }
  return(linkage_df)
} #test_linkage_df



#' Error and warning handling for probgeno_df data-frame of probabilistic genotypes (scores)
#' @param probgeno_df A data frame as read from the scores file produced by function
#' \code{saveMarkerModels} of R package \code{fitPoly}, or alternatively, a data frame containing the following columns:
#' \itemize{
#' \item{SampleName}{
#' Name of the sample (individual)
#' }
#' \item{MarkerName}{
#' Name of the marker
#' }
#' \item{P0}{
#' Probabilities of dosage score '0'
#' }
#' \item{P1...}{
#' Probabilities of dosage score '1' etc. (up to max offspring dosage, e.g. P4 for tetraploid population)
#' }
#' \item{maxP}{
#' Maximum genotype probability identified for a particular individual and marker combination
#' }
#' \item{maxgeno}{
#' Most probable dosage for a particular individual and marker combination
#' }
#' \item{geno}{
#' Most probable dosage for a particular individual and marker combination, if \code{maxP} exceeds a user-defined threshold (e.g. 0.9), otherwise \code{NA}}
#' }
#' @noRd
test_probgeno_df <- function(probgeno_df){
  if(!inherits(probgeno_df, "data.frame")) {
    warning("probgeno_df should be a data-frame. Attempting to coerce...")
    probgeno_df <- as.data.frame(probgeno_df)
  }
  if(!all(c("MarkerName", "SampleName", "geno","maxgeno","maxP") %in% colnames(probgeno_df))) stop("colnames MarkerName, SampleName, maxgeno, geno, and maxP are required!")
  if(!is.factor(probgeno_df$SampleName)) probgeno_df$SampleName <- as.factor(probgeno_df$SampleName)
  
  return(probgeno_df)
  } #test_probgeno_df


#' Large vector or in standardized matrix
#' @description Turns a vector in a matrix with fixed number of columns
#' @param x A vector
#' @param n.columns Integer, number of columns
#' @noRd
vector.to.matrix <- function(x, n.columns){
  if(length(x)>n.columns){
    x<-c(x, rep("", n.columns-length(x)%%n.columns))
  } else {
    n.columns <- length(x)
  }
  x.m <- matrix(x, ncol=n.columns, byrow=T)
  colnames(x.m)<-rep("_", n.columns)
  return(x.m)
} #vector.to.matrix


#' Write a header for the log file
#' @description Functionalized writing of function name and arguments as start for log paragraph.
#' @param matc A object of class \code{call}
#' @param log A character string specifying the log file
#' @noRd
write.logheader <- function(matc, log){
  args <- as.list(matc)
  mod <- "w"
  if(file.exists(log)) mod <- "a"
  log.conn <- file(log, mod)
  if(mod=="w") write(c("<style type=\"text/css\">.table {width: 40%;}</style>",
                       "## polymapR log file",
                       "This file is written in [markdown](https://en.wikipedia.org/wiki/Markdown).",
                       "Use knitr or [Markable](https://markable.in) to export it as a nicely formatted html or word file."),
                     log.conn)
  close(log.conn)
  mod <- "a"
  log.conn <- file(log, mod)
  write(c(paste0("\n### Log for function call of `", args[[1]], "`"), 
          as.character(Sys.time())),file = log.conn)
  close(log.conn)
  log.conn <- file(log, "a")
  write("\nWith the following call:\n", log.conn)
  write("```", log.conn)
  sink(log.conn)
  print(matc)
  sink()
  write("```\n", log.conn)
  close(log.conn)
} #write.logheader

