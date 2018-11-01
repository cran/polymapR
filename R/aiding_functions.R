#' Write a header for the log file
#' @description Functionalized writing of function name and arguments as start for log paragraph.
#' @param matc A object of class \code{call}
#' @param log A character string specifying the log file
#' @noRd
write.logheader <- function(matc, log){
  args <- as.list(matc)
  #args.df <- data.frame(names=names(args), values=as.character(args))
  #args.df <- args.df[-1,]
  mod <- "w"
  if(file.exists(log)) mod <- "a"
  log.conn <- file(log, mod)
  if(mod=="w") write(c("<style type=\"text/css\">.table {width: 40%;}</style>",
                       "##polymapR log file",
                       "This file is written in [markdown](https://en.wikipedia.org/wiki/Markdown).",
                       "Use knitr or [Markable](https://markable.in) to export it as a nicely formatted html or word file."),
                     log.conn)
  close(log.conn)
  mod <- "a"
  log.conn <- file(log, mod)
  write(c(#"\n-----------------------------------------------------------",
    paste0("\n###Log for function call of `", args[[1]], "`"), 
    as.character(Sys.time())
  ),
  #"-----------------------------------------------------------"),
  file = log.conn)
  close(log.conn)
  log.conn <- file(log, "a")
  write("\nWith the following call:\n", log.conn)
  write("```", log.conn)
  sink(log.conn)
  print(matc)
  sink()
  write("```\n", log.conn)
  #write("With the following arguments:", file=log.conn)
  #write.table(args.df, quote=F, sep=" = ", file=log.conn, row.names=F, col.names=F)
  close(log.conn)
}

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
}

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
} 


#' Plot links between 1.0 markers
#' @description Make a plot with recombination frequency versus LOD score.
#' @param linkage_df A linkage data.frame.
#' @param LG_hom_stack A \code{data.frame} as a result of \code{\link{bridgeHomologues}}
#' @param LG_number Chromosome (LG) number.
#' @param h1 Homologue to be compared.
#' @param h2 Homologue to be compared against h1
#' @param ymax Maximum limit of y-axis
#' @noRd
plot_SNlinks <- function(linkage_df,LG_hom_stack,
                         LG_number, h1, h2, ymax=NULL){
  
  h1Markers <- LG_hom_stack[LG_hom_stack[,"LG"] == LG_number & LG_hom_stack[,"homologue"] == h1,"SxN_Marker"]
  # find markers in cluster to combine
  h2Markers <- LG_hom_stack[LG_hom_stack[,"LG"] == LG_number & LG_hom_stack[,"homologue"] == h2,"SxN_Marker"]
  
  # find markers that are in both clusters 
  subdata1 <- linkage_df[linkage_df[,"marker_a"] %in% h1Markers & linkage_df[,"marker_b"] %in% h2Markers,]
  subdata2 <- linkage_df[linkage_df[,"marker_a"] %in% h2Markers & linkage_df[,"marker_b"] %in% h1Markers,]
  
  plotdata <- rbind(subdata1,subdata2)
  
  if(is.null(ymax)){ ymax<- max(plotdata$LOD)}
  
  if(nrow(plotdata) > 0){
    plot(NULL,xlab="r",ylab="LOD",
         xlim=c(0,0.5),ylim=c(0,ymax),main=paste0("LG",LG_number," h",h1," & h",h2))
    with(plotdata[plotdata$phase=="coupling",],
         points(r,LOD,pch=19,col="limegreen"))
    with(plotdata[plotdata$phase=="repulsion",],
         points(r,LOD,pch=19,col="dodgerblue"))
  } else{
    warning(paste("No SxN linkage found between h",h1," and h",h2," (LOD too low?)",sep=""))
  }
} 

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
}

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
}

#' Error and warning handling for dosage_matrix 
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @noRd
test_dosage_matrix <- function(dosage_matrix){
  if(class(dosage_matrix) == "data.frame"){
    warning("dosage_matrix should be a matrix, now it's a data.frame.")
    message("Trying to convert it to matrix, assuming markernames are in the first column..")
    rownames(dosage_matrix) <- dosage_matrix[,1]
    dosage_matrix <- as.matrix(dosage_matrix[,-1])
    class(dosage_matrix) <- "integer"
  } else if(class(dosage_matrix) == "matrix"){
    rn <- rownames(dosage_matrix)
    cn <- colnames(dosage_matrix)
    if(is.null(rn)) stop("The rownames of dosage_matrix should contain markernames. Now NULL")
    if(is.null(cn)) stop("The columnnames of dosage_matrix should contain genotype names. Now NULL")
    if(!(typeof(dosage_matrix)=="integer" | typeof(dosage_matrix)=="double")){
      warning("dosage_matrix should be integer or numeric. Trying to convert it.. ")
      class(dosage_matrix) <- "integer"
    }
  } else {
    stop("dosage_matrix should be a matrix of integers. 
         See the manual of this function for more information.")
  }
  return(dosage_matrix)
}

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
}

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
}

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
  return(LG_hom_stack)
}

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
}

#' Calculate the mode of a vector
#' @param x Vector input, can be numeric, character or factor
#' @noRd
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

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
} #check.filename

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
}

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
}

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