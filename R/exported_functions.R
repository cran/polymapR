#' Linkage analysis in polyploids
#'
#' This package uses dosage-scored (or probabilistic) bi-allelic markers from an F1 cross to perform linkage analysis in polyploids
#' @docType package
#' @name polymapR
#' @import utils stats foreach grDevices graphics
NULL

#' Add back duplicate markers after mapping
#' @description Often there will be duplicate markers that can be put aside to speed up mapping. These may be added back to the maps afterwards.
#' @param maplist A list of maps. Output of MDSMap_from_list.
#' @param bin_list A list of marker bins containing marker duplicates. One of the list outputs of \code{\link{screen_for_duplicate_markers}}
#' @param marker_assignments Optional argument to include the marker_assignments (output of \code{\link{check_marker_assignment}}). If included, marker assignment information will also be copied.
#' @return
#' A list with the following items:
#' \itemize{
#' \item{maplist}{List of maps, now with duplicate markers added}
#' \item{marker_assignments}{If required, marker assignment list with duplicate markers added}
#' }
#' @export
add_dup_markers <- function(maplist, 
                            bin_list, 
                            marker_assignments = NULL){
  #bin_list <- screened_data_list$bin_list
  maplist.out <- list()
  
  for(lg in names(maplist)){
    map <- maplist[[lg]]
    dupmark <- as.character(map$marker[map$marker %in% names(bin_list)])
    
    for(m in dupmark){
      newmark <- bin_list[[m]]
      ind <- which(map$marker == m)
      newinfo <- map[rep(ind, length(newmark)),]
      newinfo$marker <- newmark
      map <- rbind(map, newinfo)
    }
    
    maplist.out[[lg]] <- map[order(map$position),]
  }
  
  if(!is.null(marker_assignments)){
    
    if(!all(names(marker_assignments) == c("P1","P2"))) stop("marker_assignments should be the output of check_marker_assignment function, with 2 list elements named P1 and P2.")
    
    dupmark1 <- rownames(marker_assignments$P1)[rownames(marker_assignments$P1) %in% names(bin_list)]
    dupmark2 <- rownames(marker_assignments$P2)[rownames(marker_assignments$P2) %in% names(bin_list)]
    
    for(m in dupmark1){
      newmark <- bin_list[[m]]
      ind <- which(rownames(marker_assignments$P1) == m)
      newinfo <- matrix(marker_assignments$P1[rep(ind, length(newmark)),,drop=FALSE],
                        nrow = length(newmark),
                        dimnames = list(newmark,colnames(marker_assignments$P1)))
      marker_assignments$P1 <- rbind(marker_assignments$P1, newinfo)
    }
    
    for(m in dupmark2){
      newmark <- bin_list[[m]]
      ind <- which(rownames(marker_assignments$P2) == m)
      newinfo <- matrix(marker_assignments$P2[rep(ind, length(newmark)),,drop=FALSE],
                        nrow = length(newmark),
                        dimnames = list(newmark,colnames(marker_assignments$P2)))
      marker_assignments$P2 <- rbind(marker_assignments$P2, newinfo)
    }
    
    marker_assignments$P1 <- marker_assignments$P1[order(marker_assignments$P1[,"Assigned_LG"]),]
    marker_assignments$P2 <- marker_assignments$P2[order(marker_assignments$P2[,"Assigned_LG"]),]
    
    outlist <- list("maplist" = maplist.out,
                    "marker_assignments" = marker_assignments)
    
  } else{
    outlist <- list("maplist" = maplist.out,
                    "marker_assignments" = NULL)
  }
  
  return(outlist)
} #add_dup_markers

#' Assign non-SN markers to a linkage group and homologue(s).
#' @description \code{assign_linkage_group} quantifies per marker number of linkages to a linkage group and evaluates to which linkage group (and homologue(s)) the marker belongs.
#' @param linkage_df A linkage \code{data.frame} as output of \code{\link{linkage}}.
#' @param LG_hom_stack A \code{data.frame} with markernames (\code{"SxN_Marker"}), linkage group (\code{"LG"}) and homologue (\code{"homologue"})
#' @param SN_colname The name of the column in linkage_df harbouring the 1.0 markers
#' @param unassigned_marker_name The name of the column in linkage_df harbouring the marker that are to be assigned.
#' @param phase_considered The phase that is used to assign the markers (deprecated)
#' @param LG_number The number of chromosomes (linkage groups) in the species.
#' @param LOD_threshold The LOD score at which a linkage to a linkage group is significant.
#' @param ploidy The ploidy of the plant species.
#' @param assign_homologue Logical. Should markers be assigned to homologues? If \code{FALSE} markers will be assigned to all homologues
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return
#' Output is a data.frame with at least the following columns:
#' \item{Assigned_LG}{The assigned linkage group}
#' \item{Assigned_hom1}{The homologue with most linkages}
#' The columns LG1 - LGn and Hom1 - Homn give the number of hits per marker for that linkage group/homologue. Assigned_hom2 .. gives the nth homologue with most linkages.
#' @examples
#' data("SN_DN_P1", "LGHomDf_P1_1")
#' assigned_df<-assign_linkage_group(linkage_df = SN_DN_P1,
#'                      LG_hom_stack = LGHomDf_P1_1,
#'                      LG_number = 5, ploidy = 4)
#' @export
assign_linkage_group <- function(linkage_df,
                                 LG_hom_stack,
                                 SN_colname = "marker_a",
                                 unassigned_marker_name = "marker_b",
                                 phase_considered = "coupling",
                                 LG_number,
                                 LOD_threshold = 3,
                                 ploidy,
                                 assign_homologue = T,
                                 log = NULL) {
  #linkage_df <- test_linkage_df(linkage_df)
  LG_hom_stack <- test_LG_hom_stack(LG_hom_stack) #polymapR:::test_LG_hom_stack(LG_hom_stack)
  # filter for LOD threshold and phase:
  
  if(length(levels(factor(LG_hom_stack$homologue))) > ploidy){
    stop("The number of homologues per chromosome should not exceed ploidy.")
  }
  
  if(length(unique(LG_hom_stack$LG)) != LG_number){
    stop(paste("Only",length(unique(LG_hom_stack$LG)),"linkage groups were identified in LG_hom_stack. Please revise LG_number accordingly."))
  }
  
  linkage_df <-
    linkage_df[linkage_df$LOD > LOD_threshold &
                 linkage_df$phase == phase_considered,,drop = FALSE]
  
  if (is.null(linkage_df)) {
    message("There were no linkage groups the marker could be assigned to")
    return(NULL)
  }
  
  if (nrow(linkage_df) == 0) {
    message("There were no linkage groups the marker could be assigned to")
    return(NULL)
  }
  
  # get linked SNxSN markers
  SN_markers <- levels(as.factor(linkage_df[, SN_colname]))
  
  # get vector of unique markernames of unassigned markers
  unassigned_markers <-
    levels(as.factor(as.character((linkage_df[, unassigned_marker_name]))))
  
  # make sure the LG and homologue columns are factors
  LG_hom_stack$LG <- as.factor(LG_hom_stack$LG)
  LG_hom_stack$homologue <- as.factor(LG_hom_stack$homologue)
  
  # merge linkage_df and LG_hom_stack by SNxSN markers to assign unassigned markers to LG and homolog
  comb_df <-
    merge(
      linkage_df[, c(SN_colname, unassigned_marker_name)],
      LG_hom_stack,
      by.x = SN_colname,
      by.y = "SxN_Marker",
      all.x = T
    )
  
  # for every unassigned marker make a count table for markers assigned to which homolog in which LG
  count_tables <-
    tapply(1:nrow(comb_df), as.character(comb_df[, unassigned_marker_name]),
           function(x) {
             table(comb_df[x, "LG"], comb_df[x, "homologue"])
           })
  
  # get counts per chromosome (LG)
  if(LG_number > 1){
    chm.counts <- sapply(count_tables, rowSums)
  } else{
    chm.counts <- matrix(sapply(count_tables, rowSums),nrow = 1, dimnames = list(1,names(count_tables)))
  }
  
  counts_chm <- t(matrix(chm.counts,nrow = nrow(chm.counts),
                         dimnames = list(paste0("LG", rownames(chm.counts)),names(count_tables))))
  
  unlinked_markers <- rownames(counts_chm)[rowSums(counts_chm) == 0]
  
  counts_chm <-
    counts_chm[!rownames(counts_chm) %in% unlinked_markers, , drop = FALSE]
  count_tables <- "["(count_tables, rownames(counts_chm))
  
  # give warning message if markers are assigned to more than one LG:
  warn_lg <- apply(counts_chm, 1, function(x) {
    m <- max(x, na.rm = T)
    a <- m / x < 2
    return(sum(a, na.rm = T) > 1)
  })
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  if (sum(warn_lg) > 0) { #Possible bug fixed, strange it wasn't giving errors!
    write("\n####Marker(s) showing ambiguous linkage to more than one LG:\n",
          log.conn)
    amb.m <- vector.to.matrix(unassigned_markers[warn_lg], 4)
    write(knitr::kable(amb.m), log.conn)
  }
  
  # linkage group with most counts (shouldn't we get rid of markers that link to >1 LG?)
  Assigned_LG <- as.numeric(rownames(chm.counts))[apply(counts_chm, 1, which.max)]
  
  # get homolog count data for the selected linkage group
  counts_hom <- t(sapply(count_tables, function(x) {
    x[which.max(rowSums(x)),]
  }))
  
  if(ncol(counts_hom) != ploidy){
    counts_hom <- cbind(counts_hom,matrix(0,ncol = ploidy - ncol(counts_hom),
                                          nrow = nrow(counts_hom)))
  }
  
  colnames(counts_hom) <- paste0("Hom", 1:ploidy)
  
  if (assign_homologue) {
    # define assigned homologs in order of counts
    assigned_hom <- t(apply(counts_hom, 1, function(x) {
      s <- sum(x < 1)
      o <- order(x)
      o[0:s] <- NA
      return(rev(o))
    }))
    
    # get differential homolog set if phase is repulsion (deprecated?)
    if (phase_considered == "repulsion") {
      assigned_hom <- t(apply(assigned_hom, 1, function(x) {
        d <- setdiff(1:ploidy, x)
        return(c(d, rep(NA, sum(!is.na(
          x
        )))))
      }))
    }
    
  } else {
    nmark <- length(Assigned_LG)
    assigned_hom <-
      matrix(rep(1:ploidy, nmark), ncol = ploidy, byrow =
               T)
  }
  
  colnames(assigned_hom) <- paste0("Assigned_hom", 1:ploidy)
  
  output <- cbind(Assigned_LG, counts_chm, counts_hom, assigned_hom)
  
  write(paste("\n In total,",length(Assigned_LG),"out of",
              length(unassigned_markers),"markers were assigned."),
        log.conn)
  
  # write unlinked markers to standard out
  if (length(unlinked_markers) > 0) {
    write("\n####Marker(s) not assigned:\n", log.conn)
    unl.m <- vector.to.matrix(unlinked_markers, 4)
    write(knitr::kable(unl.m), log.conn)
  }
  if (!is.null(log))
    close(log.conn)
  
  return(output)
} #assign_linkage_group


#' Assign (leftover) 1.0 markers
#' @description Some 1.0 markers might have had ambiguous linkages, or linkages with low LOD scores leaving them unlinked to a linkage group.
#' \code{assign_SN_SN} finds 1.0 markers unlinked to a linkage group and tries to assign them.
#' @param linkage_df A \code{data.frame} as output of \code{\link{linkage}} with arguments markertype1=c(1,0) and markertype2=NULL.
#' @param LG_hom_stack A \code{data.frame} with markernames (\code{"SxN_Marker"}), linkage group (\code{"LG"}) and homologue (\code{"homologue"})
#' @param LOD_threshold A LOD score at which linkages between markers are significant.
#' @param ploidy Integer. The ploidy level of the plant species.
#' @param LG_number Integer. Number of chromosomes (linkage groups)
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return Returns a \code{data.frame} with the following columns:
#' \item{SxN_Marker}{The markername}
#' \item{Assigned_hom1}{The assigned homologue}
#' \item{Assigned_LG}{The assigned linkage group}
#' @examples
#' data("SN_SN_P1", "LGHomDf_P1_1")
#' SN_assigned<-assign_SN_SN(linkage_df = SN_SN_P1,
#'              LG_hom_stack = LGHomDf_P1_1,
#'              LOD_threshold= 4,
#'              ploidy=4,
#'              LG_number=5)
#' @export
assign_SN_SN <- function(linkage_df,
                         LG_hom_stack,
                         LOD_threshold,
                         ploidy,
                         LG_number,
                         log = NULL) {
  LG_hom_stack <- test_LG_hom_stack(LG_hom_stack)
  linkage_df <- test_linkage_df(linkage_df)
  all_SN <-
    unique(c(
      as.character(linkage_df$marker_a),
      as.character(linkage_df$marker_b)
    ))
  assigned_SN <- LG_hom_stack$SxN_Marker
  unassigned_SN <- all_SN[!all_SN %in% assigned_SN]
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  write(paste("Number of unassigned markers:", length(unassigned_SN)), log.conn)
  SN_SN_unass <-
    linkage_df[linkage_df$marker_a %in% unassigned_SN |
                 linkage_df$marker_b %in% unassigned_SN,]
  SN_SN_unass <-
    SN_SN_unass[!(SN_SN_unass$marker_a %in% unassigned_SN &
                    SN_SN_unass$marker_b %in% unassigned_SN),]
  change_order <- SN_SN_unass$marker_a %in% unassigned_SN
  change <- SN_SN_unass[change_order,]
  change_markerb <- change$marker_b
  change_markera <- change$marker_a
  change$marker_a <- change_markerb
  change$marker_b <- change_markera
  not_change <- SN_SN_unass[!change_order,]
  SN_SN_unass <- rbind(not_change, change)
  
  SN_SN_assigned <- assign_linkage_group(
    SN_SN_unass,
    LG_hom_stack = LG_hom_stack,
    SN_colname = "marker_a",
    unassigned_marker_name = "marker_b",
    LOD_threshold = LOD_threshold,
    ploidy = ploidy,
    LG_number = LG_number
  )
  
  if (is.null(SN_SN_assigned))
    return(NULL)
  
  ambiguous_assignments <-
    apply(SN_SN_assigned[, paste0("Hom", 1:ploidy), drop = FALSE],
          1,
          function(x) {
            m <- max(x, na.rm = T)
            a <- m / x < 2
            return(sum(a, na.rm = T) > 1)
            #sum(x>0)>1
          })
  
  unamb_ass <-
    SN_SN_assigned[!ambiguous_assignments, , drop = FALSE]
  write(paste("Number of extra assigned markers", nrow(unamb_ass)), log.conn)
  if (!is.null(log))
    close(log.conn)
  return(data.frame(
    SxN_Marker = rownames(unamb_ass),
    homologue = unamb_ass[, "Assigned_hom1"],
    LG = unamb_ass[, "Assigned_LG"]
  ))
} #assign_SN_SN

#' Use bridge markers to cluster homologues into linkage groups
#' @aliases bridgeHomologues assembleDuplexLinks
#' @description Clustering at high LOD scores results in marker clusters representing homologues.
#' \code{bridgeHomologues} clusters these (pseudo)homologues to linkage groups using linkage information between 1.0 and
#' bridge markers within a parent (e.g. 2.0 for a tetraploid).
#' If parent-specific bridge markers (e.g. 2.0) cannot be used, biparental markers can also be used (e.g. 1.1, 1.2, 2.1, 2.2 and 1.3 markers).
#' The linkage information between 1.0 and biparental markers can be combined.
#' @param cluster_stack A \code{data.frame} with a column \code{"marker"} specifying markernames,
#' and a column \code{"cluster"} specifying marker cluster
#' @param cluster_stack2 Optional. A \code{cluster_stack} for the other parent.
#' Use this argument if cross-parent markers are used (e.g. when using 1.1 markers).
#' @param linkage_df A linkage \code{data.frame} as output of \code{\link{linkage}} between bridge (e.g. 1.0 and 2.0) markers.
#' @param linkage_df2 Optional. A \code{linkage_df} specifying linkages between 1.0 and cross-parent markers in the other parent.
#' Use this argument if cross-parent markers are used (e.g. when using 1.1, 2.1, 1.2 and/or 2.2 markers).
#' The use of multiple types of cross-parent markers is allowed.
#' @param LOD_threshold Integer. The LOD threshold specifying at which LOD score a link between 1.0 and bridge (e.g. 2.0) markers is used for clustering homologues.
#' @param automatic_clustering Logical. Should clustering be executed without user input?
#' @param LG_number Integer. Expected number of chromosomes (linkage groups)
#' @param parentname Name of the parent. Used in the main title of the plot.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @param min_links The minimum number of cross-parent linkages for a marker to be considered.
#' Make this number higher if there are a lot of spurious links.
#' @param min_bridges The minimum number of linking markers to link two homologues together. 
#' @param only_coupling Logical, should only coupling linkages be used in the process? By default \code{FALSE}
#' @return A data.frame with markers classified by homologue and linkage group.
#' @examples
#' data("P1_homologues", "P2_homologues", "SN_DN_P1", "SN_SS_P1", "SN_SS_P2")
#' ChHomDf<-bridgeHomologues(cluster_stack = P1_homologues[["5"]],
#'                  linkage_df=SN_DN_P1,
#'                  LOD_threshold=4,
#'                  automatic_clustering=TRUE,
#'                  LG_number=5,
#'                  parentname="P1")
#'
#' ChHomDf<-bridgeHomologues(cluster_stack = P1_homologues[["5"]],
#'                            cluster_stack2 = P2_homologues[["5"]],
#'                  linkage_df=SN_SS_P1,
#'                  linkage_df2=SN_SS_P2,
#'                  LOD_threshold=4,
#'                  automatic_clustering=TRUE,
#'                  LG_number=5,
#'                  parentname="P1")
#' @export
bridgeHomologues <- function(
  cluster_stack,
  cluster_stack2 = NULL,
  linkage_df,
  linkage_df2 = NULL,
  LOD_threshold = 5,
  automatic_clustering = TRUE,
  LG_number,
  parentname = "",
  min_links = 1,
  min_bridges = 1,
  only_coupling = FALSE,
  log = NULL
) {
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  cluster_stack <- test_cluster_stack(cluster_stack)
  linkage_df <- test_linkage_df(linkage_df)
  if(!is.null(linkage_df2)) linkage_df2 <- test_linkage_df(linkage_df2)
  if(!is.null(cluster_stack2)) cluster_stack2 <- test_cluster_stack(cluster_stack2)
  
  linkage_df <- linkage_df[linkage_df$LOD >= LOD_threshold,]
  if(only_coupling) linkage_df <- linkage_df[linkage_df$phase == "coupling",]
  
  if (is.null(cluster_stack2) & is.null(linkage_df2)) {
    # define edges between SN and DN markers
    linkage_df <- 
      edges_linkage_df <- linkage_df[,c("marker_a", "marker_b")]
    
    homvec <- cluster_stack$cluster
    names(homvec) <- cluster_stack$marker
    
    combs <- tapply(linkage_df$marker_a, linkage_df$marker_b, function(x) x)
    
    dupcombs <- t(sapply(combs, function(x) table(homvec[x])))
    
    # Possible issue with this part from GvG May 2020, seems can be removed without issue (not fully checked)
    # for(i in seq(nrow(dupcombs))){
    #   tmp <- dupcombs[i,]
    #   dupcombs[i, -order(tmp, decreasing = T)[1:2]] <- 0
    # }
    
    pc <- t(combn(levels(homvec), 2))
    class(pc) <- "integer"
    
    weights <- rep(0, nrow(pc))
    for(i in seq(nrow(pc))){
      sub <- dupcombs[,as.character(pc[i,])]
      weights[i] <- sum(sub[,1] >= min_links & sub[,2] >= min_links)
    }
    
    edges <- as.data.frame(cbind(pc, weights))
    colnames(edges) <- c("h1", "h2", "weights")
    
  } else if (!is.null(cluster_stack2) & !is.null(linkage_df2)) {
    linkage_df2 <- linkage_df2[linkage_df2$LOD >= LOD_threshold,]
    if(only_coupling) linkage_df2 <- linkage_df2[linkage_df2$phase == "coupling",]
    
    assign_cluster <- function(linkage_df,
                               cluster_stack) {
      
      SN_markers <- levels(as.factor(linkage_df[, "marker_a"]))
      
      # get vector of unique markernames of unassigned markers
      unassigned_markers <-
        levels(as.factor(as.character((linkage_df[, "marker_b"]))))
      
      # make sure the LG and homologue colums are factors
      cluster_stack$cluster <-
        as.factor(cluster_stack$cluster)
      
      # merge linkage_df and cluster_stack by SNxSN markers to assign unassigned markers to LG and homolog
      comb_df <-
        merge(
          linkage_df[, c("marker_a", "marker_b")],
          cluster_stack,
          by.x = "marker_a",
          by.y = "marker",
          all.x = T
        )
      
      # for every unassigned marker make a count table for markers assigned to which homolog in which LG
      max_cluster <-
        tapply(1:nrow(comb_df), as.character(comb_df[, "marker_b"]),
               function(x) {
                 t <- table(comb_df[x, "cluster"])
                 m <- max(t)
                 if (m >= min_links) {
                   ass <- names(t)[which.max(t)]
                 } else {
                   ass <- NA
                 }
                 return(ass)
               })
      class(max_cluster) <- "integer"
      return(max_cluster)
    }
    
    assigned1 <- assign_cluster(linkage_df,
                                cluster_stack)
    assigned2 <- assign_cluster(linkage_df2,
                                cluster_stack2)
    mboth <- intersect(names(assigned1), names(assigned2))
    assigned_both <- cbind(assigned1[mboth], assigned2[mboth])
    assigned_both <- assigned_both[complete.cases(assigned_both),]
    
    edges <-
      tapply(assigned_both[, 1], assigned_both[, 2], function(x) {
        u <- unique(x) # add weights for number of links
        if(length(u) < 2) return(NULL)
        g <- t(combn(u,2))
        weights <- rep(0, nrow(g))
        for(i in seq(nrow(g))){
          weights[i] <- min(sum(x == g[i,1]), sum (x == g[i,2]))
        }
        out <- as.data.frame(cbind(g, weights))
        colnames(out) <- c("h1", "h2", "weights")
        return(out)
      })
    # rbind edgelist
    edges <- do.call(rbind, edges)
    
    for(i in seq(nrow(edges))){
      edges[i, 1:2] <- edges[i,1:2][order(edges[i,1:2])]
    }
    edges <- aggregate(edges$weights, by = list(edges[,1], edges[,2]), FUN = sum)
    colnames(edges) <- c("h1", "h2", "weights")
    
  } else {
    stop("Both cluster_stack2 and linkage_df2 should be specified or none of the two.")
  }
  
  edges <- edges[edges$weights >= min_bridges,]
  # make network between SN x SN clusters
  nw <- igraph::graph.data.frame(edges[,c("h1", "h2")], directed = F)
  nw <- igraph::set_edge_attr(nw, "weight", value = edges[,"weights"])
  gcl <- igraph::groups(igraph::clusters(nw))
  
  # plot network
  igraph::plot.igraph(
    nw,
    mark.groups = gcl,
    mark.col = "white",
    vertex.size = 10,
    layout = igraph::layout_nicely(nw),
    vertex.color = "white",
    edge.width = igraph::E(nw)$weight,
    main = paste(
      parentname,"SNxSN homologue links with LOD >",
      LOD_threshold
    )
  )
  
  repr_cluster_stack <-
    levels(cluster_stack$cluster)[!levels(cluster_stack$cluster) %in% unlist(gcl)]
  
  if (length(repr_cluster_stack) > 0) {
    sub <- cluster_stack[cluster_stack$cluster %in% repr_cluster_stack,]
    t <- table(droplevels(sub$cluster))
    message("Not all homologue clusters are represented in the graph.")
    write(paste0("\n####Unrepresented homologues\n"),
          file = log.conn)
    #sink(log.conn)
    write(knitr::kable(data.frame(cluster = names(t), size = as.vector(t))),
          log.conn)
    #suppressWarnings(sink())
  }
  
  if (automatic_clustering) {
    # automatic clustering proceeds only if number of clusters are expected chromosome number
    
    # list to df
    hom_cl <- stack(gcl)
    colnames(hom_cl) <- c("homologue", "LG")
    colnames(cluster_stack) <- c("SxN_Marker", "homologue")
    # merge markers to homolog df
    SN_SN_chm_cl <- merge(cluster_stack, hom_cl, by = "homologue")
    SN_SN_chm_cl <-
      SN_SN_chm_cl[, c("SxN_Marker","LG","homologue")]
    if (length(gcl) != LG_number) {
      warning(
        "Number of clusters do not represent number of expected chromosomes.
        If you're not okay with that, try again with different LOD_cutoff, min_links, min_bridges or automatic_clustering = FALSE"
      )
    }
  } else {
    # user input which homologs per chromosome
    SN_SN_cluster_list <- unstack(cluster_stack)
    SN_SN_chm_cl <-
      createChmHomList(SN_SN_cluster_list, LG_number = LG_number)
  }
  
  for (chm in levels(SN_SN_chm_cl$LG)) {
    clusters_chm <-
      as.factor(as.numeric(SN_SN_chm_cl$homologue[SN_SN_chm_cl$LG == chm]))
    levels(clusters_chm) <- 1:length(levels(clusters_chm))
    SN_SN_chm_cl$homologue[SN_SN_chm_cl$LG == chm] <-
      as.numeric(clusters_chm)
  }
  SN_SN_chm_cl$homologue <-
    as.factor(as.character(SN_SN_chm_cl$homologue))
  
  if (!is.null(log))
    close(log.conn)
  
  return(SN_SN_chm_cl)
} #bridgeHomologues

# calcSegtypeInfo ***********************************************************
#'@title Build a list of segregation types
#'@description For each possible segregation type in an F1 progeny with given
#'parental ploidy (and ploidy2, if parent2 has a different ploidy than parent1)
#'information is given on the segregation ratios, parental dosages and whether
#'the segregation is expected under polysomic, disomic and/or mixed inheritance.
#'
#'@usage calcSegtypeInfo(ploidy, ploidy2=NULL)
#'
#'@param ploidy The ploidy of parent 1 (must be even, 2 (diploid) or larger).
#'@param ploidy2 The ploidy of parent 2. If omitted (default=NULL) it is
#'assumed to be equal to ploidy.
#'@details The names of the segregation types consist of a short sequence of
#'digits (and sometimes letters), an underscore and a final number. This is
#'interpreted as follows, for example segtype 121_0: 121 means that there
#'are three consecutive dosages in the F1 population with frequency ratios 1:2:1,
#'and the 0 after the underscore means that the lowest of these dosages is
#'nulliplex. So 121_0 means a segregation of 1 nulliplex : 2 simplex : 1 duplex.
#'A monomorphic F1 (one single dosage) is indicated as e.g. 1_4 (only one
#'dosage, the 4 after the underscore means that this is monomorphic quadruplex).
#'If UPPERCASE letters occur in the first part of the name these are interpreted
#'as additional digits with values of A=10 to Z=35, e.g. 18I81_0 means a
#'segregation of 1:8:18:8:1 (using the I as 18), with the lowest dosage being
#'nulliplex.\cr
#'With higher ploidy levels higher numbers (above 35) may be required.
#'In that case each unique ratio number above 35 is assigned a lowercase letter.
#'E.g. one segregation type in octaploids is 9bcb9_2: a 9:48:82:48:9
#'segregation where the lowest dosage is duplex.\cr
#'Segregation types with more than 5 dosage classes are considered "complex"
#'and get codes like c7e_1 (again in octoploids): this means a complex type
#'(the first c) with 7 dosage classes; the e means that this is the fifth
#'type with 7 classes. Again the _1 means that the lowest dosage is simplex.
#'It is always possible (and for all segtype names with lowercase letters it is
#'necessary) to look up the actual segregation ratios in the intratio item
#'of the segtype. For octoploid segtype c7e_1 this shows 0:1:18:69:104:69:18:1:0
#'(the two 0's mean that nulli- and octoplexes do not occur).
#'
#'@return A list with for each different segregation type (segtype) one item.
#'The names of the items are the names of the segtypes.
#'Each item is itself a list with components:
# (alternatives for itemize: enumerate and describe,
# see https://cran.r-project.org/web/packages/roxygen2/vignettes/formatting.html)
#'\itemize{
#'\item{freq}{a vector of the ploidy+1 fractions of the dosages in the F1}
#'\item{intratios}{an integer vector with the ratios as the simplest integers}
#'\item{expgeno}{a vector with the dosages present in this segtype}
#'\item{allfrq}{the allele frequency of the dosage allele in the F1}
#'\item{polysomic}{boolean: does this segtype occur with polysomic inheritance?}
#'\item{disomic}{boolean: does this segtype occur with disomic inheritance?}
#'\item{mixed}{boolean: does this segtype occur with mixed inheritance (i.e. with
#'polysomic inheritance in one parent and disomic inheritance in the other)?}
#'\item{pardosage}{integer matrix with 2 columns and as many rows as there
#'are parental dosage combinations for this segtype;
#'each row has one possible combination of dosages for
#'parent 1 (1st column) and parent 2 (2nd column)}
#'\item{parmode}{logical matrix with 3 columns and the same number of rows as
#'pardosage. The 3 columns are named polysomic, disomic and mixed and
#'tell if this parental dosage combination will generate this
#'segtype under polysomic, disomic and mixed inheritance}
#'}
#'@examples
#'si4 <- calcSegtypeInfo(ploidy=4) # two 4x parents: a 4x F1 progeny
#'print(si4[["11_0"]])
#'
#'si3 <- calcSegtypeInfo(ploidy=4, ploidy2=2) # a 4x and a diplo parent: a 3x progeny
#'print(si3[["11_0"]])
#'@export
calcSegtypeInfo <- function(ploidy, ploidy2=NULL) {
  if (ploidy %% 2 != 0) stop("calcSegtypeInfo: odd ploidy not allowed")
  if (is.null(ploidy2)) ploidy2 <- ploidy else
    if (ploidy2 %% 2 != 0) stop("calcSegtypeInfo: odd ploidy2 not allowed")
  ploidyF1 <- (ploidy + ploidy2) / 2
  
  result <- list()
  #we check all parental combinations and enumerate the segtypes:
  for (p1 in 0:ploidy) {
    if (p1 %in% c(0, 1, ploidy-1, ploidy)) {
      p1ir <- matrix(polygamfrq(ploidy, p1), nrow=1)
      p1mode <- 0 #0=don't care
    } else {
      p1ir <- rbind(polygamfrq(ploidy, p1), digamfrq(ploidy, p1))
      p1mode <- c(1, rep(2, nrow(p1ir)-1)) #1=polysomic, 2=disomic
    }
    if (ploidy2 == ploidy) maxp2 <- p1 else maxp2 <- ploidy2
    for (p2 in 0:maxp2) {
      if (p2 %in% c(0, 1, ploidy2-1, ploidy2)) {
        p2ir <- matrix(polygamfrq(ploidy2, p2), nrow=1)
        p2mode <- 0 #0=don't care
      } else {
        p2ir <- rbind(polygamfrq(ploidy2, p2), digamfrq(ploidy2, p2))
        p2mode <- c(1, rep(2, nrow(p2ir)-1)) #1=polysomic, 2=disomic
      }
      for (i1 in 1:nrow(p1ir)) for (i2 in 1:nrow(p2ir)) {
        ir <- makeProgeny(p1ir[i1,], p2ir[i2,])
        #check if this segtype already exists in result:
        r <- length(result)
        while (r > 0 && sum(ir != result[[r]]$intratios) > 0) r <- r-1
        if (r == 0) {
          #new segtype:
          r <- length(result) + 1
          result[[r]] <- list()
          result[[r]]$freq <- ir/sum(ir)
          result[[r]]$intratios <- ir
          result[[r]]$expgeno <- which(ir > 0) - 1; names(result[[r]]$expgeno) <- NULL
          result[[r]]$allfrq <- sum(result[[r]]$freq * 0:ploidyF1) / ploidyF1
          result[[r]]$polysomic <- FALSE
          result[[r]]$disomic <- FALSE
          result[[r]]$mixed <- FALSE
          result[[r]]$pardosage <- matrix(c(p1, p2), nrow=1)
          colnames(result[[r]]$pardosage) <- c("P1", "P2")
          result[[r]]$parmode <- matrix(FALSE, nrow=1, ncol=3)
          colnames(result[[r]]$parmode) <- c("polysomic", "disomic", "mixed")
          p <- 1
        } else {
          #segtype already existed, see if parental combination
          #already existed:
          p <- nrow(result[[r]]$pardosage)
          while (p > 0 &&
                 sum(c(p1,p2) != result[[r]]$pardosage[p,]) > 0) p <- p-1
          if (p == 0) {
            #this parental dosage combination is new
            p <- nrow(result[[r]]$pardosage) + 1
            result[[r]]$pardosage <- rbind(result[[r]]$pardosage, c(p1, p2))
            result[[r]]$parmode <- rbind(result[[r]]$parmode,
                                         rep(FALSE, 3))
          }
        }
        result[[r]]$parmode[p, 1] <- result[[r]]$parmode[p, 1] ||
          (p1mode[i1] != 2 && p2mode[i2] != 2) #is polysomic?
        result[[r]]$parmode[p, 2] <- result[[r]]$parmode[p, 2] ||
          (p1mode[i1] != 1 && p2mode[i2] != 1) #is disomic?
        result[[r]]$parmode[p, 3] <- result[[r]]$parmode[p, 3] ||
          (p1mode[i1] == 0 || p2mode[i2] == 0 || p1mode[i1] != p2mode[i2]) #mixed?
      } #for i1 for i2
    } #for p2
  } #for p1
  #next we check if each segtype occurs with poly- or disomic or
  #mixed inheritance,
  #we add all reciprocal parental dosages,
  #we order all pardosage and parmode in order of increasing p1,
  #we check which intratios occur in segtypes with less than complex classes,
  #and we compile the numdosages and firstdosage statistics:
  numdosages <- firstdosage <- maxratio <- integer(length(result))
  complex <- 6 #the number of dosage classes from which we consider the segtype
  #             complex, i.e. we don't show the full segregation in the name
  simpleratios <- integer(0)
  for (r in seq_along(result)) {
    #calculate poly- and disomic and mixed:
    result[[r]]$polysomic <- sum(result[[r]]$parmode[, 1]) > 0
    result[[r]]$disomic <- sum(result[[r]]$parmode[, 2]) > 0
    result[[r]]$mixed <- sum(result[[r]]$parmode[, 3]) > 0
    #add the reverse of parents with unequal dosage:
    if (ploidy2 == ploidy) {
      #if ploidy2 not null, all parental combinations already done
      uneq <- result[[r]]$pardosage[,1] != result[[r]]$pardosage[,2]
      if (sum(uneq) > 0) {
        #there are rows with unequal parental dosage
        result[[r]]$pardosage <- rbind(result[[r]]$pardosage,
                                       result[[r]]$pardosage[uneq, 2:1])
        result[[r]]$parmode <- rbind(result[[r]]$parmode,
                                     result[[r]]$parmode[uneq,])
      }
    }
    #sort pardosage and parmode in order of increasing P1 dosage:
    or <- order(result[[r]]$pardosage[,1])
    result[[r]]$pardosage <- result[[r]]$pardosage[or,, drop=FALSE]
    result[[r]]$parmode <- result[[r]]$parmode[or,, drop=FALSE]
    #add the intratios in not complex:
    if (length(result[[r]]$expgeno) < complex)
      simpleratios <- union(simpleratios, result[[r]]$intratios)
    #calculate statistics for later ordering of result list:
    numdosages[r] <- length(result[[r]]$expgeno)
    firstdosage[r] <- result[[r]]$expgeno[1]
    maxratio[r] <- max(result[[r]]$intratios)
  }
  #next we reorder result from simple to complex segtypes
  or <- order(numdosages, maxratio, firstdosage)
  result <- result[or]
  #finally we assign names according to numdosages and firstdosage:
  numdosages <- numdosages[or]
  firstdosage <- firstdosage[or]
  maxratio <- maxratio[or]
  simple <- max(which(numdosages < complex))
  alldigits <- c("0","1","2","3","4","5","6","7","8","9", LETTERS)
  #we translate the intratios of non-complex segtypes to a sequence of digits
  #and uppercase letters (for ratios 10:35); if any ratios > 35 occur (as in
  #octo- and dodecaploids) we assign one of the lowercase letters for each of
  #these dosages (and assume there are less than 27 unique ratios above 35).
  #Just using the lowercase letter for ratios 36:61 does not work as already
  #with octoploids some intratios of non-complex segtypes are above 61.
  simpleratios <- sort(simpleratios[simpleratios > 35])
  for (r in 1:simple) {
    chr <- character(complex-1)
    name <- ""
    ir <- result[[r]]$intratios
    for (i in 1:(ploidyF1+1)) {
      if (ir[i] > 0) {
        if (ir[i] <= 35) {
          name <- paste(name, alldigits[ir[i] + 1], sep="")
        } else name <- paste(name, letters[which(simpleratios == ir[i])], sep="")
      }
    }
    names(result)[r] <- paste(name, firstdosage[r], sep="_")
  }
  if (complex <= max(numdosages)) for (n in complex:max(numdosages)) {
    #we want to assign codes like c6a_1, where:
    #c just means complex
    #6 is the number of different dosages in the segtype
    #c6a, c6b etc are segtypes with different ratio combinations, all with 6 dosages
    #_1 means that the first dosage is simplex
    selnum <- numdosages == n
    maxrat <- sort(unique(maxratio[selnum]))
    for (m in 1:length(maxrat)) {
      seltype <- selnum & maxratio == maxrat[m]
      if (length(maxrat) == 1) name <- paste("c", n, sep="") else
        name <- paste("c", n, letters[m], sep="")
      firstrange <- sort(unique(firstdosage[selnum]))
      for (f in firstrange) {
        sel <- seltype & firstdosage == f
        names(result)[sel] <- paste(name, f, sep="_")
      }
    }
  }
  result
} #calcSegtypeInfo

#  checkF1 *******************************************************************
#' @title Identify the best-fitting F1 segregation types 
#' @description For a given set of F1 and parental samples, this function
#' finds the best-fitting segregation type using either discrete or probabilistic input data. 
#' It can also perform a dosage shift prior to selecting the segregation type.
#' @param input_type Can be either one of 'discrete' or 'probabilistic'. For the former (default), a \code{dosage_matrix} must be supplied,
#' while for the latter a \code{probgeno_df} must be supplied. 
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
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
#' Most probable dosage for a particular individual and marker combination, if \code{maxP} exceeds a user-defined threshold (e.g. 0.9), otherwise \code{NA}
#' }
#' }
#' @param parent1 character vector with the sample names of parent 1
#' @param parent2 character vector with the sample names of parent 2
#' @param F1 character vector with the sample names of the F1 individuals
#' @param ancestors character vector with the sample names of any other
#' ancestors or other samples of interest. The dosages of these samples will
#' be shown in the output (shifted if shiftParents \code{TRUE}) but they are not used
#' in the selection of the segregation type.
#' @param polysomic if \code{TRUE} at least all polysomic segtypes are considered;
#' if \code{FALSE} these are not specifically selected (but if e.g. disomic is \code{TRUE},
#' any polysomic segtypes that are also disomic will still be considered)
#' @param disomic if \code{TRUE} at least all disomic segtypes are considered (see
#' \code{polysomic})
#' @param mixed if \code{TRUE} at least all mixed segtypes are considered (see
#' \code{polysomic}). A mixed segtype occurs when inheritance in one parent is
#' polysomic (random chromosome pairing) and in the other parent disomic (fully
#' preferential chromosome pairing)
#' @param ploidy The ploidy of parent 1 (must be even, 2 (diploid) or larger).
#' @param ploidy2 The ploidy of parent 2. If omitted it is
#' assumed to be equal to ploidy.
#' @param outfile the tab-separated text file to write the output to; if NA a temporary file
#' checkF1.tmp is created in the current working directory and deleted at end
#' @param critweight NA or a numeric vector containing the weights of three quality
#' criteria; do not need to sum to 1. If NA, the output will not contain a
#' column qall_weights. Else the weights specify how qall_weights will be
#' calculated from quality parameters q1, q2 and q3.
#' @param Pvalue_threshold a minimum threshold value for the Pvalue of the
#' bestParentfit segtype (with a smaller Pvalue the q1 quality parameter will
#' be set to 0)
#' @param fracInvalid_threshold a maximum threshold for the fracInvalid of the
#' bestParentfit segtype (with a larger fraction of invalid dosages in the F1
#' the q1 quality parameter will be set to 0)
#' @param fracNA_threshold a maximum threshold for the fraction of unscored F1
#' samples (with a larger fraction of unscored samples in the F1
#' the q3 quality parameter will be set to 0)
#' @param shiftmarkers if specified, shiftmarkers must be a data frame with
#' columns MarkerName and shift; for the markernames that match exactly
#' (upper/lowercase etc) those in the input (either \code{dosage_matrix} or \code{probgeno_df}), the dosages are increased by the
#' amount specified in column shift,
#' e.g. if shift is -1, dosages 2..ploidy are converted to 1..(ploidy-1)
#' and dosage 0 is a combination of old dosages 0 and 1, for all samples.
#' The segregation check is then performed with the shifted dosages.
#' A shift=NA is allowed, these markers will not be shifted.
#' The sets of markers in the input (either \code{dosage_matrix} or \code{probgeno_df}) and shiftmarkers
#' may be different, but markers may occur only once in shiftmarkers.
#' A column shift is added at the end of the returned data frame.\cr
#' If parameter shiftParents is \code{TRUE}, the parental and ancestor scores are
#' shifted as the F1 scores, if \code{FALSE} they are not shifted.
#' @param parentsScoredWithF1 \code{TRUE} if parents are scored in the same experiment
#' and the same \code{fitPoly} run as the F1, else \code{FALSE}.
#' If \code{TRUE}, their fraction missing scores
#' and conflicts tell something about the quality of the scoring. If \code{FALSE}
#' (e.g. when the F1 is triploid and the parents are diploid and tetraploid) the
#' quality of the F1 scores can be independent of that of the parents.\cr
#' If not specified, \code{TRUE} is assumed if ploidy2 == ploidy and \code{FALSE} if
#' ploidy2 != ploidy
#' @param shiftParents only used if parameter shiftmarkers is specified. If \code{TRUE},
#' apply the shifts also to the parental and ancestor scores.
#' By default \code{TRUE} if \code{parentsScoredWithF1} is \code{TRUE}
#' @param showAll (default \code{FALSE}) if \code{TRUE}, for each segtype 3 columns
#' are added to the returned data frame with the frqInvalid, Pvalue and
#' matchParents values for these segtype (see the description of the return value)
#' @param append_shf if \code{TRUE} and parameter shiftmarkers is specified, _shf is
#' appended to all marker names where shift is not 0. This is not required for
#' any of the functions in this package but may prevent duplicated marker names
#' when using other software. 
#' @details For each marker is tested how well the different segregation types
#' fit with the observed parental and F1 dosages. The results are summarized
#' by columns bestParentfit (which is the best fitting segregation type,
#' taking into account the F1 and parental dosages) and columns qall_mult
#' and/or qall_weights (how good is the fit of the bestParentfit segtype: 0=bad,
#' 1=good).\cr
#' Column bestfit in the results gives the segtype best fitting the F1
#' segregation without taking account of the parents. This bestfit segtype is
#' used by function correctDosages, which tests for possible "shifts" in
#' the marker models.\cr
#' In case the parents are not scored together with the F1 (e.g. if the F1 is
#' triploid and the parents are diploid and tetraploid) \code{dosage_matrix}
#' should be edited to contain the parental as well as the F1 scores.
#' In case the diploid and tetraploid parent are scored in the same run of
#' function \code{saveMarkerModels} (from package \code{fitPoly})
#' the diploid is initially scored as nulliplex-duplex-quadruplex (dosage 0, 2
#' or 4); that must be converted to the true diploid dosage scores (0, 1 or 2).
#' Similar corrections are needed with other combinations, such as a diploid
#' parent scored together with a hexaploid population etc.
#'
#' @return A list containing two elements, \code{checked_F1} and \code{meta}. \code{meta} is itself
#' a list that stores the parameter settings used in running \code{checkF1} which can 
#' be useful for later reference. The first element (\code{checked_F1}) contains the actual results: a data
#' frame with one row per marker, with the following columns:
#' \itemize{
#' \item{m: the sequential number of the marker (as assigned by \code{fitPoly})}
#' \item{MarkerName: the name of the marker, with _shf appended if the marker
#' is shifted and append_shf is \code{TRUE}}
#' \item{parent1: consensus dosage score of the samples of parent 1}
#' \item{parent2: consensus dosage score of the samples of parent 2}
#' \item{F1_0 ...	F1_<ploidy>: the number of F1 samples with dosage scores
#' 0 ... <ploidy>}
#' \item{F1_NA: the number of F1 samples with a missing dosage score}
#' \item{sample names of parents and ancestors: the dosage scores for those
#' samples}
#' \item{bestfit: the best fitting segtype, considering only the F1 samples}
#' \item{frqInvalid_bestfit: for the bestfit segtype, the frequency of F1 samples
#' with a dosage score that is invalid (that should not occur). The frequency is
#' calculated as the number of invalid samples divided by the number of non-NA
#' samples}
#' \item{Pvalue_bestfit: the chisquare test P-value for the observed
#' distribution of dosage scores vs the expected fractions. For segtypes
#' where only one dosage is expected (1_0, 1_1 etc) the binomial probability of
#' the number of invalid scores is given, assuming an error
#' rate of seg_invalidrate (hard-coded as 0.03)}
#' \item{matchParent_bestfit: indication how the bestfit segtype matches the
#' consensus dosages of parent 1 and 2: "Unknown"=both parental
#' dosages unknown; "No"=one or both parental dosages known
#' and conflicting with the segtype; "OneOK"= only one parental
#' dosage known, not conflicting with the segtype; "Yes"=both
#' parental dosages known and combination matching with
#' the segtype. This score is initially assigned based on
#' only high-confidence parental consensus scores; if
#' low-confidence dosages are confirmed by the F1, the
#' matchParent for (only) the selected segtype is
#' updated, as are the parental consensus scores.}
#' \item{bestParentfit: the best fitting segtype that does not conflict with
#' the parental consensus scores}
#' \item{frqInvalid_bestParentfit, Pvalue_bestParentfit,
#' matchParent_bestParentfit: same as the corresponding columns for bestfit.
#' Note that matchParent_bestParentfit cannot be "No".}
#' \item{q1_segtypefit: a value from 0 (bad) to 1 (good), a measure of the fit of
#' the bestParentfit segtype based on Pvalue, invalidP and whether bestfit is
#' equal to bestParentfit}
#' \item{q2_parents: a value from 0 (bad) to 1 (good), based either on the
#' quality of the parental scores (the number of missing scores and of
#' conflicting scores, if parentsScoredWithF1 is TRUE) or on matchParents
#' (No=0, Unknown=0.65, OneOK=0.9, Yes=1, if parentsScoredWithF1 is FALSE)}
#' \item{q3_fracscored: a value from 0 (bad) to 1 (good), based on the fraction
#' of F1 samples that have a non-missing dosage score}
#' \item{qall_mult: a value from 0 (bad) to 1 (good), a summary quality score
#' equal to the product q1*q2*q3. Equal to 0 if any of these is 0, hence
#' sensitive to thresholds; a natural selection criterion would be to accept
#' all markers with qall_mult > 0}
#' \item{qall_weights: a value from 0 (bad) to 1 (good), a weighted average of
#' q1, q2 and q3, with weights as specified in parameter critweight. This column is
#' present only if critweight is specified. In this case there is no "natural"
#' threshold; a threshold for selection of markers must be obtained by inspecting
#' XY-plots of markers over a range of qall_weights values}
#' \item{shift: if shiftmarkers is specified a column shift is added with
#' for all markers the applied shift (for the unshifted markers the shift value
#' is 0)}
#' }
#' qall_mult and/or qall_weights can be used to compare the quality
#' of the SNPs within one analysis and one F1 population but not between analyses
#' or between different F1 populations.\cr
#' If parameter showAll is \code{TRUE} there are 3 additional columns for each
#' segtype with names frqInvalid_<segtype>, Pvalue_<segtype> and
#' matchParent_<segtype>; see the corresponding columns for bestfit for an
#' explanation. These extra columns are inserted directly before the bestfit
#' column.
#' @examples 
#' \dontrun{
#' data("ALL_dosages")
#' chk1<-checkF1(input_type="discrete",dosage_matrix=ALL_dosages,parent1="P1",parent2="P2",
#' F1=setdiff(colnames(ALL_dosages),c("P1","P2")),polysomic=T,disomic=F,mixed=F,
#' ploidy=4)
#' data("gp_df")
#' chk1<-checkF1(input_type="probabilistic",probgeno_df=gp_df,parent1="P1",parent2="P2",
#' F1=setdiff(levels(gp_df$SampleName),c("P1","P2")),polysomic=T,disomic=F,mixed=F,
#' ploidy=4)
#' }
#' @export
checkF1 <- function(input_type = "discrete",
                    dosage_matrix,
                    probgeno_df,
                    parent1,
                    parent2,
                    F1,
                    ancestors=character(0),
                    polysomic, 
                    disomic, 
                    mixed,
                    ploidy, 
                    ploidy2,
                    outfile = "",
                    critweight=c(1.0, 0.4, 0.4),
                    Pvalue_threshold=0.0001,
                    fracInvalid_threshold=0.05,
                    fracNA_threshold=0.25,
                    shiftmarkers,
                    parentsScoredWithF1=TRUE,
                    shiftParents=parentsScoredWithF1,
                    showAll=FALSE,
                    append_shf=FALSE) {
  input_type <- match.arg(input_type, choices = c("discrete","probabilistic"))
  
  if(input_type == "discrete"){
    dosage_matrix <- test_dosage_matrix(dosage_matrix)
  } else{
    probgeno_df <- test_probgeno_df(probgeno_df)
  }
  
  if (ploidy %% 2 != 0) stop("checkF1: odd ploidy not allowed")
  if (missing(ploidy2) || is.na(ploidy2)) ploidy2 <- ploidy else
    if (ploidy2 %% 2 != 0) stop("checkF1: odd ploidy2 not allowed")
  
  ploidyF1 <- (ploidy + ploidy2) / 2
  
  if (!polysomic && !disomic && !mixed)
    stop("checkF1: at least one of polysomic, disomic and mixed must be TRUE")
  seginfo <- calcSegtypeInfo(ploidy, ploidy2)
  allsegtypenames <- names(seginfo)
  seginfo <- selSegtypeInfo(seginfo, polysomic, disomic, mixed)
  seginfoSummary <- segtypeInfoSummary(seginfo)
  #check critweight
  if (is.null(critweight) || is.na(critweight[1])) critweight <- NA else {
    if (!is.numeric(critweight) || length(critweight) != 3 ||
        sum(is.na(critweight)) > 0 || sum(critweight) == 0) {
      stop("invalid critweight")
    }
  }
  #check shiftmarkers:
  if (missing(shiftmarkers)) shiftmarkers <- NA
  if (is.data.frame(shiftmarkers)) {
    if (sum(is.na(match(c("MarkerName", "shift"), names(shiftmarkers)))) > 0)
      stop("checkF1: shiftmarkers must have columns MarkerName and shift")
    shiftmarkers$shift[is.na(shiftmarkers$shift)] <- 0
    if (nrow(shiftmarkers) == 0) shiftmarkers <- NA else {
      if (sum(is.na(shiftmarkers$shift)) > 0 ||
          sum(!(shiftmarkers$shift %in% ((-ploidyF1):ploidyF1))) > 0)
        stop("checkF1: shiftmarkers contains invalid shift values")
      shiftmarkers$MarkerName <- as.character(shiftmarkers$MarkerName)
      if (length(unique(shiftmarkers$MarkerName)) < nrow(shiftmarkers))
        stop("checkF1: some markernames occur more than once in shiftmarkers")
    }
  }
  
  seg_invalidrate <- 0.03 #assumed percentage of scoring errors leading to invalid dosages
  file.del <- is.na(outfile) || outfile == ""
  if (file.del) outfile <- "checkF1.tmp"
  if (!checkFilename(outfile))
    stop(paste("checkF1: cannot write file", outfile))
  #fill in missing vectors and/or clean up:
  if (missing(parent1) || is.logical(parent1)) parent1 <- character(0)
  if (missing(parent2) || is.logical(parent2)) parent2 <- character(0)
  if (missing(ancestors) || is.logical(ancestors)) ancestors <- character(0)
  parent1 <- as.character(parent1[!is.na(parent1)])
  parent2 <- as.character(parent2[!is.na(parent2)])
  ancestors <- as.character(ancestors[!is.na(ancestors)])
  # names of quality parameters:
  qnames <- c("q1_segtypefit"," q2_parents","q3_fracscored","qall_mult")
  if (length(critweight) == 3) qnames[5] <- "qall_weights"
  if(input_type == "discrete"){
    mrknames <- rownames(dosage_matrix)
  } else{
    # mrknames <- as.character(probgeno_df$MarkerName)
    mrknames <- sort(as.character(unique(probgeno_df$MarkerName)))
  }
  
  createResultsdf <- function(mrkcount) {
    #function within checkF1
    saf <- getOption("stringsAsFactors")
    options(stringsAsFactors = FALSE)
    #create the results data frame for this batch:
    mat <- matrix(
      integer((2+ploidyF1+2+length(parent1)+length(parent2)+length(ancestors)) *
                mrkcount), nrow=mrkcount)
    colnames(mat) <- c("parent1", "parent2", paste("F1", 0:ploidyF1, sep="_"),
                       "F1_NA", parent1, parent2, ancestors)
    #segdf is the template for all individual segtype data frames
    #(sets of 3 columns, and one row for each marker):
    segdf <- data.frame(frqInvalid=character(mrkcount), #printed as formatted strings
                        Pvalue=character(mrkcount),
                        matchParent=factor(rep(NA,mrkcount),
                                           levels=c("No","OneOK","Unknown","Yes")))
    #bres is the data frame that will hold the results (one row per marker)
    bres <- data.frame(
      m=integer(mrkcount),
      MarkerName=character(mrkcount),
      mat
    )
    if (showAll) for (sg in 1:length(seginfo)) {
      df <- segdf
      names(df) <- paste(names(df), names(seginfo)[sg], sep="_")
      bres <- data.frame(bres, df)
    }
    segdf <- data.frame(
      fit = factor(rep(NA, mrkcount), levels=names(seginfo)),
      segdf
    )
    for (fi in c("bestfit", "bestParentfit")) {
      df <- segdf
      names(df) <- c(fi,
                     paste(names(df)[2:4], fi, sep="_"))
      bres <- data.frame(bres, df)
    }
    bres <- data.frame (
      bres,
      q1_segtypefit=character(mrkcount),
      q2_parents=character(mrkcount),
      q3_fracscored=character(mrkcount),
      qall_mult=character(mrkcount)
    )
    if (length(critweight) == 3)
      bres <- data.frame(
        bres,
        qall_weights=character(mrkcount)
      )
    if (is.data.frame(shiftmarkers))
      bres <- data.frame(
        bres,
        shift=integer(mrkcount))
    options(stringsAsFactors = saf)
    bres
  } #createResultsdf within checkF1
  
  segtypeBestSelcrit <- function(candidates) {
    #function within checkF1
    #candidates is a vector indexing the segtypes in results from which
    #           the one with the best selcrit must be selected
    #return value: index to that segtype, or 0 if none
    if (length(candidates) == 0) return(0)
    candSelcrit <- selcrit[candidates]
    candidates[which.max(candSelcrit)]
  } # segtypeBestSelcrit within checkF1
  
  compareFit <- function(newsegtype, oldsegtype) {
    #function within checkF1
    #newsegtype, oldsegtype: two indices into results
    #return value: TRUE if newsegtype is selected,
    #              FALSE if oldsegtype is selected
    (results$fracInvalid[newsegtype] <=
       max(0.05, 1.5 * results$fracInvalid[oldsegtype])) &&
      (results$Pvalue[newsegtype] >=
         min(0.01, 0.1 * results$Pvalue[oldsegtype]))
  } #compareFit within checkF1
  
  shiftdosages <- function(dosages, shift, ploidyF1) {
    #dosages: vector of integer dosages, shift: one shift value
    dosages <- dosages + shift
    below <- !is.na(dosages) & dosages < 0
    dosages[below] <- 0
    above <- !is.na(dosages) & dosages > ploidyF1
    dosages[above] <- ploidyF1
    dosages
  }
  
  #subsetting from a big dosage_matrix takes a lot of time
  #(approx 2 sec for each marker, in a file of ~26000 markers and 480 samples)
  #Therefore we subdivide it in batches of 100 markers
  batchsize <- 100
  batchnr <- 1
  while (batchsize * (batchnr-1) < length(mrknames)) {
    minmrk <- batchsize*(batchnr-1) + 1
    maxmrk <- min(length(mrknames), batchsize*batchnr)
    
    if(input_type == "discrete"){
      batchscores <- dosage_matrix[mrknames[minmrk:maxmrk],] #only 260 * 2 sec
    } else{
      batchscores <- probgeno_df[probgeno_df$MarkerName %in% mrknames[minmrk:maxmrk],]
    }
    
    bres <- createResultsdf(maxmrk - minmrk + 1) #bres: batchresults
    
    ## If we use genotype probabilities, the most different part is the count of the probability
    #function to count in each dosages, what is the sum of the probabilities
    count_probabi <- function(ploidy,sc){
      geno.count <- c()
      for(i in 0:ploidy){
        geno.ea <- sum(as.numeric(as.character(sc[[paste0("P",i)]])),na.rm=TRUE)
        geno.count <- c(geno.count, geno.ea)
      }
      return(geno.count)
    }
    
    
    for (mrk in minmrk:maxmrk) {
      
      if(input_type == "discrete"){
        mrknr <- which(rownames(batchscores) == mrknames[mrk])
        sc <- data.frame(SampleName=colnames(batchscores),
                         geno=as.integer(batchscores[mrknr,]))
        # mrknr <- rownames(batchscores)[mrknr] #marker number in fitPoly input
      } else{
        # mrknr <- which(batchscores$MarkerName == mrknames[mrk])
        # #possible bug here if batchscores has a different column order than expected...
        # sc <- data.frame(SampleName = names(batchscores)[3:length(batchscores)],
        #                  geno = as.integer(batchscores[mrknr, 3:length(batchscores)]))
        # # mrknr <- batchscores$marker[mrknr]
        
        sc <- batchscores[batchscores$MarkerName == mrknames[mrk],]
        if ("marker" %in% names(sc)) mrknr <- sc$marker[1] else mrknr <- mrk
      }
      
      parent.geno <- list()
      
      if(input_type == "discrete"){
        #get all the parent and ancestor genotypes in the order of the given vectors:
        parent.geno[[1]] <- sc$geno[match(parent1, sc$SampleName)]
        parent.geno[[2]] <- sc$geno[match(parent2, sc$SampleName)]
        ancestors.geno <- sc$geno[match(ancestors, sc$SampleName)]
        F1.geno <- sc$geno[sc$SampleName %in% F1] #F1 does not have to be ordered
      } else{
        #probability use the maxgeno of parent
        parent.geno[[1]] <- sc$maxgeno[match(parent1, sc$SampleName)]
        parent.geno[[2]] <- sc$maxgeno[match(parent2, sc$SampleName)]
        ancestors.geno <- sc$maxgeno[match(ancestors, sc$SampleName)]
        #take the colsum of each category
        chosen_col <- paste0("P",seq(0,ploidy,1))
        F1.geno <- sc[sc$SampleName %in% F1,chosen_col]
      }
      
      
      shift <- 0
      if (is.data.frame(shiftmarkers)) {
        whichshift <- which(shiftmarkers$MarkerName == mrknames[mrk])
        if (length(whichshift) == 1) {
          shift <- shiftmarkers$shift[whichshift]
          if (shift != 0) {
            #changed: now all shifted dosages outside 0..ploidyF1 are set
            #to 0 or ploidyF1 (i.e. dosages that were separated are now merged),
            #before 20161130 they were set to NA
            #F1.geno <- ancestors.geno + shift
            #F1.geno[!(F1.geno %in% 0:ploidyF1)] <- NA
            F1.geno <- shiftdosages(F1.geno, shift, ploidyF1)
            if (shiftParents) {
              #also here new shiftdosages is used:
              parent.geno[[1]] <- shiftdosages(parent.geno[[1]], shift, ploidyF1)
              parent.geno[[2]] <- shiftdosages(parent.geno[[2]], shift, ploidyF1)
              ancestors.geno <- shiftdosages(ancestors.geno, shift, ploidyF1)
            }
          }
        }
      }
      
      #get the consensus genotype and conflict status of each of the two parents:
      par.geno <- c(0, 0) #integer
      par.lowconf.geno <- c(0, 0) #integer
      par.conflicts <- c(FALSE, FALSE) #logical
      par.NAfrac <- c(0.5, 0.5) #floatingpoint
      for (parent in 1:2) {
        #if (strict.parents) maxNAfrac <- 0.499 else maxNAfrac <- 1
        parresult <- getConsensusGeno(geno=parent.geno[[parent]],
                                      maxNAfrac=0.499,
                                      lowconf.NAfrac=0.751)
        par.geno[parent] <- parresult$geno
        # Update: following line was converted to numeric in PG version - necessary?
        par.lowconf.geno[parent] <- as.numeric(as.character(parresult$lowconf.geno))
        par.conflicts[parent] <- parresult$conflict
        par.NAfrac[parent] <- parresult$NAfrac
      }
      
      if(input_type == "discrete"){
        #get the freq. distribution of the F1 samples:
        F1.naCount <- sum(is.na(F1.geno))
        F1.nobs <- length(F1.geno) - F1.naCount #the number of non-NA F1 samples
        F1.counts <- tabulate(bin=F1.geno + 1, nbins=ploidyF1 + 1)
      } else{
        F1.naCount <- sum(rowSums(is.na(F1.geno)) == ploidyF1 + 1)
        F1.nobs <- nrow(F1.geno) - F1.naCount
        F1.counts <- count_probabi(ploidy = ploidy,
                                   sc = F1.geno)
        #correct the errors
        proba_correct <- 0.05 * nrow(F1.geno)/(ploidy + 1)
        F1.counts[F1.counts < proba_correct] <- 0
        F1.counts <- (F1.counts/sum(F1.counts))*F1.nobs  #re-normalize to the same individual number
      }
      
      #note that F1.counts[1] has the number of F1 samples with geno==0 etc !
      bestfit <- NA; bestParentfit <- NA;
      q <- rep(NA, length(qnames))
      results <- data.frame(segtype=names(seginfo),
                            fracInvalid=rep(1.0, length(seginfo)),
                            invalidP=rep(0.0, length(seginfo)),
                            Pvalue=rep(0.0, length(seginfo)),
                            matchParents=I(as.character(rep(NA,length(seginfo)))))
      if (F1.nobs > 10) {
        for (s in 1:length(seginfo)) {
          #does the current segtype match the parental genotypes?
          results$matchParents[s] <-
            getMatchParents(parGeno=par.geno, seginfoItem=seginfo[[s]])
          #calculate fraction invalid F1 scores given the current segtype:
          exp.geno <- seginfo[[s]]$expgeno #the genotypes that can occur given the segtype
          if(input_type == "discrete"){
            F1.invalid <- length(F1.geno[!(F1.geno %in% exp.geno)]) - F1.naCount #number of invalid scores 
          } else{
            F1.invalid <- sum(F1.counts) - sum(F1.counts[exp.geno+1])
          }
          
          results$fracInvalid[s] <- F1.invalid / F1.nobs
          #calculate chi-square P-value given the current segtype and invalidP value:
          if (F1.nobs - F1.invalid > 0) {
            #at least some F1 samples in expected peaks
            #(else Pvalue and invalidP remain 0)
            results$invalidP[s] <-
              pbinom(q=F1.nobs - F1.invalid,   #nr of valid scores
                     size=F1.nobs,             #total nr of scores
                     prob=1 - seg_invalidrate) #expected fraction of valid scores
            if (length(exp.geno) == 1) {
              results$Pvalue[s] <- 1.0
            } else {
              #more than one genotype expected, chisquare test
              if (sum(F1.counts[exp.geno + 1]) == 0) {
                #all observations are in invalid classes
                results$Pvalue[s] <- 0.0
                results$invalidP[s] <- 0.0
              } else {
                suppressWarnings(
                  results$Pvalue[s] <-
                    chisq.test(F1.counts[exp.geno+1],
                               p=seginfo[[s]]$freq[exp.geno+1])$p.value)
                #warnings probably generated by class counts < 5;
                #these do not worry us
                #as those P-values are likely to be very significant anyway
                #also the reported warnings are not informative
              }
            }
          }
        } # for (s in 1:length(seginfo))
        #select best fitting segtype based on fracInvalid and Pvalue
        #selcrit: selection criterion
        selcrit <- results$invalidP * results$Pvalue
        bestfit <- which.max(selcrit) #discards NA's
        if (bestfit == 0)
          stop(paste("Error in checkF1: bestfit is 0 at marker", mrknames[mrk]))
        
        #select which segtype best fits the parental genotype(s)
        ParentFit <- which(results$matchParents %in% c("Yes","OneOK","Unknown"))
        bestParentfit <- segtypeBestSelcrit(ParentFit)
        if (bestParentfit == 0)
          stop(paste("Error in checkF1: bestParentfit is 0 at marker", mrknames[mrk]))
        #bestParentfit is 0 if all segtypes have matchParents=="No" (should never occur)
        
        #Can we improve on bestParentfit by using low-confidence parental scores?
        lowParentFit <- NA
        lowc <- which(!is.na(par.lowconf.geno))
        if (length(lowc) > 0) {
          if (length(lowc) == 1) {
            if (is.na(par.geno[3 - lowc])) {
              #one parent missing, other lowconf
              low.segtypes <- seginfoSummary$segtypenr[
                seginfoSummary[, 2+lowc] == par.lowconf.geno[lowc]]
            } else {
              #one parent high conf, other lowconf score
              low.segtypes <- seginfoSummary$segtypenr[
                (seginfoSummary[, 2+lowc] == par.lowconf.geno[lowc]) &
                  (seginfoSummary[, 5-lowc] == par.geno[3-lowc])]
            }
            if (length(low.segtypes) > 0) {
              lowParentfit <- segtypeBestSelcrit(low.segtypes)
              #now we must check if lowParentfit is equal to or not
              #too much worse than bestParentfit: if so, we promote the
              #par.lowconf.geno to a true geno and lowParentfit to
              #bestParentfit
              if (compareFit(lowParentfit, bestParentfit)) {
                par.geno[lowc] <- par.lowconf.geno[lowc]
                bestParentfit <- lowParentfit
                results$matchParents[bestParentfit] <-
                  getMatchParents(parGeno=par.geno,
                                  seginfoItem=seginfo[[bestParentfit]])
              }
            }
          } else {
            #both parents have a low-conf score:
            #we should check for both scores if they are acceptable
            #(i.e. equal to, or not too much worse than bestfit)
            #If one of them is ok and not the other, we promote it and its lowParentfit;
            #if both pass and their combination is ok we promote both
            #and their combined lowParentfit;
            #if both pass but the combination does not pass we don't
            #promote either and keep bestfit as bestParentfit
            #First the combination:
            low.segtypes <- seginfoSummary$segtypenr[
              (seginfoSummary[, 3] == par.lowconf.geno[1]) &
                (seginfoSummary[, 4] == par.lowconf.geno[2])]
            lowParentfit <- segtypeBestSelcrit(low.segtypes)
            if (compareFit(lowParentfit, bestParentfit)) {
              #the combination of both par.lowconf.geno is acceptable
              par.geno <- par.lowconf.geno
              bestParentfit <- lowParentfit
              results$matchParents[bestParentfit] <-
                getMatchParents(parGeno=par.geno,
                                seginfoItem=seginfo[[bestParentfit]])
              #always "Yes"?
            } else {
              #the combination of both par.lowconf.geno is not acceptable,
              #check the two parents separately:
              lowParentfit <- c(0,0)
              for (p in 1:2) {
                low.segtypes <- seginfoSummary$segtypenr[
                  seginfoSummary[, 2+p] == par.lowconf.geno[p]]
                lowParentfit[p] <- segtypeBestSelcrit(low.segtypes)
                if (!compareFit(lowParentfit[p], bestParentfit))
                  lowParentfit[p] <- 0
              }
              p <- which(lowParentfit != 0)
              if (length(p) == 1) {
                #just one parent gives an acceptable lowParentfit, use that one:
                par.geno[p] <- par.lowconf.geno[p]
                bestParentfit <- lowParentfit[p]
                results$matchParents[bestParentfit] <-
                  getMatchParents(parGeno=par.geno,
                                  seginfoItem=seginfo[[bestParentfit]])
                #always "OneOK"?
              }
            } #two par.lowconf.geno, the combination was not acceptable
          } #two par.lowconf.geno
        } #at least one par.lowconf.geno
        
        #Note that we don't attempt here to resolve mismatches between the final
        #par.geno and bestParentfit, nor to fill in missing par.geno if the
        #bestParentfit segtype would allow only one solution.
        #That is done in compareProbes (and therefore also in writeDosagefile).
        #A further expansion to generate markers with all possible par.geno
        #allowed by the segtype and the known par.geno is done in
        #expandUnknownParents
        
        if(input_type == "discrete"){
          F1_NAfrac <- F1.naCount/length(F1.geno)
        } else {
          F1_NAfrac <- F1.naCount/nrow(F1.geno)
        }
        
        # generate a set of quality parameters of the marker in general and
        # the fit of bestParentfit:
        q <- calc_qall(Pvalue_threshold, fracInvalid_threshold, fracNA_threshold,
                       Pvalue=results$Pvalue[bestParentfit],
                       fracInvalid=results$fracInvalid[bestParentfit],
                       F1.NAfrac=F1_NAfrac,
                       matchParents=results$matchParents[bestParentfit],
                       bestfit=bestfit, bestParentfit=bestParentfit,
                       par.conflicts=par.conflicts,
                       par.NAfrac=par.NAfrac,
                       critweight=critweight,
                       parentsScoredWithF1=parentsScoredWithF1)
      } # F1.nobs>10
      
      #get the output line:
      bix <- mrk - minmrk + 1 #index to row in bres
      bres$m[bix] <- mrknr
      if (shift==0 || !append_shf) bres$MarkerName[bix] <- mrknames[mrk] else
        bres$MarkerName[bix] <- paste(mrknames[mrk], "shf", sep="_")
      bres[bix, 3:4] <- par.geno
      bres[bix, 5:(5+ploidyF1)] <- F1.counts
      bres$F1_NA[bix] <- F1.naCount
      startcol <- ploidyF1 + 7
      if (length(parent1) > 0) {
        bres[bix, startcol:(startcol-1+length(parent1))] <- parent.geno[[1]]
        startcol <- startcol + length(parent1)
      }
      if (length(parent2) > 0) {
        bres[bix, startcol:(startcol-1+length(parent2))] <- parent.geno[[2]]
        startcol <- startcol + length(parent2)
      }
      if (length(ancestors) > 0) {
        bres[bix, startcol:(startcol-1+length(ancestors))] <- ancestors.geno[[1]]
        startcol <- startcol + length(ancestors)
      }
      if (showAll) {
        bres[bix, startcol+seq(0, by=3, length.out=length(seginfo))] <-
          sprintf("%.4f", results$fracInvalid)
        bres[bix, startcol+seq(1, by=3, length.out=length(seginfo))] <-
          sprintf("%.4f", results$Pvalue)
        bres[bix, startcol+seq(2, by=3, length.out=length(seginfo))] <-
          results$matchParents
        startcol <- startcol + 3 * length(seginfo)
      }
      if (is.na(bestfit)) {
        bres[bix, startcol:(startcol+3)] <- rep(NA, 4)
      } else {
        bres[bix, startcol] <- results$segtype[bestfit]
        bres[bix, startcol+1] <- sprintf("%.4f", results$fracInvalid[bestfit])
        bres[bix, startcol+2] <- sprintf("%.4f", results$Pvalue[bestfit])
        bres[bix, startcol+3] <- results$matchParents[bestfit]
      }
      startcol <- startcol + 4
      if (is.na(bestParentfit)) {
        bres[bix, startcol:(startcol+3)] <- rep(NA, 4)
      } else {
        bres[bix, startcol] <- results$segtype[bestParentfit]
        bres[bix, startcol+1] <- sprintf("%.4f", results$fracInvalid[bestParentfit])
        bres[bix, startcol+2] <- sprintf("%.4f", results$Pvalue[bestParentfit])
        bres[bix, startcol+3] <- results$matchParents[bestParentfit]
      }
      startcol <- startcol + 4
      if (is.na(q[1])) {
        bres[bix, startcol:(startcol-1+length(q))] <- rep(NA, length(q))
      } else {
        bres[bix, startcol:(startcol-1+length(q))] <- sprintf("%.4f", q)
      }
      startcol <- startcol + length(q)
      if (is.data.frame(shiftmarkers)) bres[bix, startcol] <- shift
    } # for mrk
    if (batchnr ==1) {
      write.table(bres, file=outfile, quote=FALSE, sep="\t",
                  na="", row.names=FALSE, col.names=TRUE)
    } else {
      write.table(bres, file=outfile, append=TRUE, quote=FALSE, sep="\t",
                  na="", row.names=FALSE, col.names=FALSE)
    }
    batchnr <- batchnr + 1
  } #while batch
  output <- read.table(outfile, header=TRUE, sep="\t",
                       na.strings="", check.names=FALSE)
  if(input_type == "discrete"){
    output <- chk2integer(output)
  }
  
  
  output$bestfit <- factor(as.character(output$bestfit), levels=allsegtypenames)
  output$bestParentfit <- factor(as.character(output$bestParentfit), levels=allsegtypenames)
  if (file.del) file.remove(outfile)
  
  return(invisible(list("checked_F1" = output,
                        "meta" = list(parent1 = parent1,
                                      parent2 = parent2,
                                      F1 = F1,
                                      ancestors = ancestors,
                                      polysomic = polysomic, 
                                      disomic = disomic, 
                                      mixed = mixed,
                                      ploidy = ploidy, 
                                      ploidy2 = ploidy2,
                                      outfile = outfile,
                                      critweight = critweight,
                                      Pvalue_threshold = Pvalue_threshold,
                                      fracInvalid_threshold = fracInvalid_threshold,
                                      fracNA_threshold = fracNA_threshold,
                                      shiftmarkers = shiftmarkers,
                                      parentsScoredWithF1 = parentsScoredWithF1))))
} #checkF1


#' Check the quality of a linkage map using heatplots
#' @description Perform a series of checks on a linkage map and visualise the results using heatplots. Also shows the discrepency between
#' the pairwise and multi-point r estimates, plotted against the LOD of the pairwise estimate.
#' @param linkage_list A named \code{list} with r and LOD of markers within linkage groups.
#' @param maplist A list of maps. In the first column marker names and in the second their position.
#' @param mapfn The map function used in generating the maps, either one of "haldane" or "kosambi". By default "haldane" is assumed.
#' @param lod.thresh Numeric. Threshold for the LOD values to be displayed in heatmap, by default 5 (set at 0 to display all values)
#' @param tidyplot If \code{TRUE}, an attempt is made to reduce the plot density, using the \code{hexbin} package. 
#' This can have a considerable performance impact for high-density maps
#' @param detail Level of detail for heatmaps, by default 1 cM. Values less than 0.5 cM can have serious performance implications.
#' @param sortmarkers If \code{TRUE} (by default) the markers in the linkage_list are sorted: first the lower position, then the higher. 
#' Results in an averaging (in plot B) over all markers at a gives pair of positions 
#' @param plottype Option to specify graphical device for plotting, (either png or pdf), or by default "", in which case plots are directly plotted within R
#' @param prefix Optional prefix appended to plot names if outputting plots.
#' @examples
#' \dontrun{
#' data("maplist_P1","all_linkages_list_P1")
#' check_map(linkage_list = all_linkages_list_P1, maplist = maplist_P1)
#' }
#' @export
check_map <- function (linkage_list, 
                       maplist, 
                       mapfn = "haldane", 
                       lod.thresh = 5, 
                       tidyplot = TRUE, 
                       detail = 1, 
                       sortmarkers=TRUE, 
                       plottype=c("", "pdf", "png")[1],
                       prefix="") { 
  if (length(linkage_list) != length(maplist) ||
      !all(names(linkage_list) == names(maplist))) 
    stop("linkage_list and maplist do not correspond.")
  mapfn <- match.arg(mapfn, choices = c("haldane", "kosambi"))
  rev.haldane <- function(d) (1 - exp(-d/50))/2
  rev.kosambi <- function(d) ((exp(d/25) - 1)/(exp(d/25) + 
                                                 1))/2
  orig.mar <- c(5.1, 4.1, 4.1, 2.1)
  colbar.mar <- c(5.1, 2, 4.1, 0.5)
  
  lgfunc <- function(l) {
    lgname <- names(linkage_list)[l]
    posmat <- matrix(c(maplist[[l]][match(linkage_list[[l]]$marker_a, 
                                          maplist[[l]]$marker), ]$position, 
                       maplist[[l]][match(linkage_list[[l]]$marker_b, 
                                          maplist[[l]]$marker), ]$position, 
                       linkage_list[[l]]$r, 
                       linkage_list[[l]]$LOD), 
                     ncol = 4)
    # sort the rows to have the positions in ascending order
    # (so for plot B the LOD and r of all marker pairs at the same positions
    #  are averaged):
    if (sortmarkers)
      posmat[,1:2] <- t(apply(posmat[,1:2], 1, function(x) x[order(x)]))
    if (mapfn == "haldane") {
      expected.recom <- rev.haldane(abs(posmat[, 1] - posmat[, 2]))
    } else if (mapfn == "kosambi") {
      expected.recom <- rev.kosambi(abs(posmat[, 1] - posmat[, 2]))
    } else stop("incorrect mapfn")
    dev <- linkage_list[[l]]$r - expected.recom
    wRMSE <- sqrt(mean((abs(dev) * linkage_list[[l]]$LOD)^2))
    # Plot A: LOD vs delta_r
    if (plottype=="pdf") {
      pdf(paste0(prefix, lgname, "_check_map_plotA.pdf"), 
          height = 5, width = 5)
    } else if (plottype=="png") {
      png(paste0(prefix, lgname, "_check_map_plotA.png"), 
          height = 500, width = 500)
      
    }
    if (tidyplot) {
      hbin <- hexbin::hexbin(dev, linkage_list[[l]]$LOD, 
                             xbins = 150)
      suppressWarnings(hexbin::plot(hbin, ylab = "LOD", 
                                    xlab = expression(delta(r)), legend = FALSE, 
                                    main = bquote(.(lgname) ~ ": r"["pairwise"] ~ 
                                                    "- r"["map"])))
    }
    else {
      plot(dev, linkage_list[[l]]$LOD, ylab = "LOD", 
           xlab = expression(delta(r)), cex.lab = 1.25, 
           main = expression("|r"["pairwise"] ~ 
                               "- r"["map"] ~ "|"))
    }
    message(paste0(lgname, " weighted RMSE = ", round(wRMSE, 3)))
    if (plottype %in% c("pdf", "png")) {
      dev.off()
      Sys.sleep(0.5)
    }  
    if (plottype=="pdf") {
      pdf(paste0(prefix, lgname, "_check_map_plotB.pdf"), 
          height = 5, width = 8)
    } else if (plottype=="png") {
      png(paste0(prefix, lgname, "_check_map_plotB.png"), 
          height = 500, width = 1000)
      
    }
    layout(matrix(1:4, ncol = 4, byrow = TRUE), widths = c(1, 
                                                           0.2, 1, 0.2))
    par(oma = c(0, 0, 3, 0))
    #ds1 <- seq(min(posmat[, 1]), max(posmat[, 1]), detail) # position marker_a steps
    #ds2 <- seq(min(posmat[, 2]), max(posmat[, 2]), detail) # position marker_b steps
    ds <- seq(min(posmat[,1]), max(posmat[,2]), detail) # position steps; position1 always <= position2
    fi1 <- findInterval(posmat[, 1], ds) # in which interval is the position of the lowest marker
    fi2 <- findInterval(posmat[, 2], ds) # in which interval is the position of the highest marker
    fi1[is.na(fi1)] <- max(fi1)
    fi2[is.na(fi2)] <- max(fi2)
    fi1[fi1 == 0] <- 1
    fi2[fi2 == 0] <- 1
    L1 <- tapply(posmat[, 4], INDEX = list(fi1, fi2), mean) # mean LOD over all marker pairs in this pairwise position-bin
    # "tapply calls FUN for each cell that has any data in it"; therefore:
    L1[which(is.na(L1))] <- 0
    r1 <- tapply(posmat[, 3], INDEX = list(fi1, fi2), mean) # mean r
    r1[which(is.na(r1))] <- 0
    
    expandmatrix <- function(m, size, default, symmetric) {
      # if m has missing rows or columns, insert these and fill with default value
      if (nrow(m) < size | ncol(m) < size) {
        x <- matrix(0, nrow=size, ncol=ncol(m))
        x[as.integer(rownames(m)),] <- m 
        cn <- as.integer(colnames(m))
        m <- matrix(0, nrow=nrow(x), ncol=size)
        m[, cn] <- x
      } 
      if (symmetric) {
        # convert the triangular data to symmetric data
        m <- m + t(m)
        # works because the empty cells, including one triangle and the diagonal, 
        # are zeroes
      }
      m
    }
    
    L1 <- expandmatrix(L1, length(ds), 0, symmetric=sortmarkers)
    r1 <- expandmatrix(r1, length(ds), 0, symmetric=sortmarkers)
    colours <- colorRampPalette(c("green", "yellow", 
                                  "orange", "red"))(100)
    rcolmin <- min(r1)
    rcolmax <- max(r1)
    LODcolmin <- lod.thresh
    LODcolmax <- max(L1)
    LODcolbreaks <- seq(LODcolmin, LODcolmax, (LODcolmax - LODcolmin)/100)
    LODhits <- findInterval(x = L1, vec = LODcolbreaks)
    LODcols <- rep("white", length(LODhits))
    LODcols[LODhits != 0] <- colours[LODhits[LODhits != 0]]
    LODcols[is.na(LODcols)] <- colours[100]
    LODcols <- matrix(LODcols, ncol = length(ds))
    rcolbreaks <- seq(rcolmin, rcolmax, (rcolmax - rcolmin)/100)
    rhits <- findInterval(x = r1, vec = rcolbreaks)
    rcols <- rep("white", length(rhits))
    rcols[rhits != 0] <- colours[rhits[rhits != 0]]
    rcols[is.na(rcols)] <- colours[100]
    rcols <- matrix(rcols, ncol = length(ds))
    plot(NULL, xlim = range(ds), ylim = range(ds), main = "r", 
         cex = 3, xlab = "cM", ylab = "cM")
    for (i in 1:(length(ds) - 1)) {
      for (j in 1:(length(ds) - 1)) {
        rect(xleft = ds[i], ybottom = ds[j], 
             xright = ds[i + 1], ytop = ds[j + 1], 
             col = rcols[i, j], border = NA)
      }
    }
    par(mar = colbar.mar)
    colour.bar(col.data = colours, min = rcolmin, max = rcolmax, 
               nticks = 8, 
               ticks = round(seq(rcolmin, rcolmax, len = 8), 2), 
               cex.ticks = 0.8)
    par(mar = orig.mar)
    plot(NULL, xlim = range(ds), ylim = range(ds), main = "LOD", 
         cex = 3, xlab = "cM", ylab = "cM")
    for (i in 1:(length(ds) - 1)) {
      for (j in 1:(length(ds) - 1)) {
        rect(xleft = ds[i], ybottom = ds[j], 
             xright = ds[i + 1], ytop = ds[j + 1], col = LODcols[i, j], 
             border = NA)
      }
    }
    par(mar = colbar.mar)
    colour.bar(col.data = colours, min = LODcolmin, max = LODcolmax, 
               nticks = 8, 
               ticks = round(seq(LODcolmin, LODcolmax, len = 8), 2), 
               cex.ticks = 0.8)
    mtext(text = paste0(lgname, " map diagnostics"), 
          side = 3, outer = TRUE, cex = 2)
    if (plottype %in% c("pdf", "png")) {
      dev.off()
    }
  }
  
  sapply(seq(length(linkage_list)), lgfunc)
  par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = orig.mar)
} #check_map


#' Check for consistent marker assignment between both parents
#' @description Function to ensure there is consistent marker assignment to chromosomal linkage groups
#' for biparental markers
#' @param marker_assignment.P1 A marker assignment matrix for parent 1 with markernames as rownames and at least containing the column \code{"Assigned_LG"}; the output of \code{\link{homologue_lg_assignment}}.
#' @param marker_assignment.P2 A marker assignment matrix for parent 2 with markernames as rownames and at least containing the column \code{"Assigned_LG"}; the output of \code{\link{homologue_lg_assignment}}.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL (by default) log is send to stdout.
#' @param verbose Should messages be sent to stdout or log?
#' @examples
#' data("marker_assignments_P1"); data("marker_assignments_P2")
#' check_marker_assignment(marker_assignments_P1,marker_assignments_P2)
#' @return Returns a list of matrices with corrected marker assignments.
#' @export
check_marker_assignment <- function(marker_assignment.P1,
                                    marker_assignment.P2,
                                    log = NULL,
                                    verbose = TRUE){
  
  outlist <- list("P1" = marker_assignment.P1,
                  "P2" = marker_assignment.P2)
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  biparentals <- intersect(rownames(marker_assignment.P1),rownames(marker_assignment.P2))
  
  if(verbose) message(paste("Checking consistency of LG assignment of", length(biparentals), "bi-parental markers..."))
  
  mismatch <- which(marker_assignment.P1[match(biparentals,rownames(marker_assignment.P1)),"Assigned_LG"] !=
                      marker_assignment.P2[match(biparentals,rownames(marker_assignment.P2)),"Assigned_LG"])
  
  if(length(mismatch) > 0){
    if(verbose) message(paste("Inconsistent assignment detected for",length(mismatch),"markers. These will be removed.."))
    
    mismatch.df <- data.frame("marker" = biparentals[mismatch],
                              "P1" = marker_assignment.P1[match(biparentals[mismatch],rownames(marker_assignment.P1)),"Assigned_LG"],
                              "P2" = marker_assignment.P2[match(biparentals[mismatch],rownames(marker_assignment.P2)),"Assigned_LG"])
    rownames(mismatch.df) <- 1:nrow(mismatch.df)
    
    if(verbose) {
      write("The following bi-parental markers were assigned to different linkage groups across parents and were subsequently removed:\n", log.conn)
      write(knitr::kable(mismatch.df), log.conn)
    }
    
    outlist$P1 <- outlist$P1[-which(rownames(outlist$P1) %in% mismatch.df$marker),]
    outlist$P2 <- outlist$P2[-which(rownames(outlist$P2) %in% mismatch.df$marker),]
  } else{
    if(verbose) message("\nFull consistency in LG assignment found!")
  }
  
  if (!is.null(log))
    close(log.conn)
  
  return(outlist)
} #check_marker_assignment



#' check your dataset's maxP distribution
#' @description Function to assess the distribution of maximum genotype probabilities (\code{maxP}), if these are available. The function
#' plots a violin graph showing the distribution of the samples' \code{maxP}.
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
#' Most probable dosage for a particular individual and marker combination, if \code{maxP} exceeds a user-defined threshold (e.g. 0.9), otherwise \code{NA}
#' }
#' }
#' @examples
#' data("gp_df")
#' check_maxP(gp_df)
#' @return This function does not return any value, is simply a visualisation tool to help assess data quality.
#' @export
check_maxP <- function(probgeno_df){
  
  probgeno_df <- test_probgeno_df(probgeno_df)
  
  #1. Missing
  proba <- probgeno_df$maxP
  proba_NA <- sum(is.na(proba))  #nbr of missing
  proba_nonNA <- proba[complete.cases(proba)]
  writeLines(paste0(round(proba_NA/length(proba)*100,2),'% samples do not fit mixture model and result in missing values.'))
  
  #2. Random pick numbers for plot (otherwise this step takes too long)
  sample_size <- 5000
  if(length(proba)/10 < sample_size){
    sample_size <- round(length(proba)/10)
  }
  plot_sample <- sample(proba,sample_size)
  probability <- data.frame('maxP' = plot_sample,'Sample' = as.factor(""))
  
  p <- ggplot2::ggplot(probability, ggplot2::aes_string("Sample","maxP"))
  suppressWarnings(print(p + ggplot2::geom_violin(fill = "dodgerblue",colour = "navyblue")))
  
  #3. each threshold number count
  window_choice <- seq(0,1,0.05)
  count_summary <- data.frame()
  for(i in 1:(length(window_choice)-1)){
    min <- window_choice[i]
    max <- window_choice[i+1]
    count <- length(proba_nonNA[proba_nonNA > min & proba_nonNA <= max])
    count_summary <- rbind(count_summary,data.frame(paste0(min,' ~ ',max),count))
  }
  colnames(count_summary)[1] <- 'interval'
  count_summary <- rbind(count_summary,data.frame('interval' = 'NA', 'count' = proba_NA))
  count_summary$percent <- paste0(round(count_summary$count/length(proba)*100,2),'%')
  writeLines('\nSummary of the number of samples within each probability interval:')
  print.data.frame(count_summary)
  
} #check_maxP


#' Cluster 1.0 markers into correct homologues per linkage group
#' @description Clustering at one LOD score for all markers does usually not result in correct classification of homologues. Usually there are more clusters of (pseudo)homologues than expected. This function lets you inspect every linkage group separately and allows for clustering at a different LOD threshold per LG.
#' @param LG Integer. Linkage group to investigate.
#' @param linkage_df A data.frame as output of \code{\link{linkage}} with arguments \code{markertype1 = c(1,0)} and \code{markertype2=NULL}.
#' @param LG_hom_stack A \code{data.frame} with columns \code{"SxN_Marker"} providing 1.0 markernames and \code{"LG"}
#' and \code{"homologue"} providing linkage group and homologue respectively.
#' @param LOD_sequence A numeric or vector of numerics giving LOD threshold(s) at which clustering should be performed.
#' @param modify_LG_hom_stack Logical. Should \code{LG_hom_stack} be modified and returned?
#' @param nclust_out Number of clusters in the output. If there are more clusters than this number only the nclust_out largest clusters are returned.
#' @param network.layout Network layout: \code{"circular"} or \code{"stacked"}. If \code{"n"} no network is plotted.
#' @param device Function of the graphics device to plot to (e.g. \code{\link{pdf}}, \code{\link{png}}, \code{\link{jpeg}}). The active device is used when \code{NULL}
#' @param label.offset Offset of labels. Only used if \code{network.layout="circular"}.
#' @param cex.lab label character expansion. Only for \code{network.layout="circular"}.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @param \dots Arguments passed to \code{device}.
#' @return A modified LG_hom_stack \code{data.frame} if \code{modify_LG_hom_stack = TRUE}
#' @examples
#' data("SN_SN_P2", "LGHomDf_P2_1")
#' #take only markers in coupling:
#' SN_SN_P2_coupl <- SN_SN_P2[SN_SN_P2$phase=="coupling",]
#' cluster_per_LG(LG = 2,
#'                linkage_df=SN_SN_P2_coupl,
#'                LG_hom_stack=LGHomDf_P2_1,
#'                LOD_sequence=seq(4,10,2),
#'                modify_LG_hom_stack=FALSE,
#'                nclust_out=4,
#'                network.layout="circular",
#'                device=NULL,
#'                label.offset=1.2,
#'                cex.lab=0.75)
#' @export
cluster_per_LG <- function(LG,
                           linkage_df,
                           LG_hom_stack,
                           LOD_sequence,
                           modify_LG_hom_stack = FALSE,
                           nclust_out = NULL,
                           network.layout = c("circular", "stacked", "n"),
                           device = NULL,
                           label.offset = 1,
                           cex.lab = 0.7,
                           log = NULL,
                           ...) {
  linkage_df <- test_linkage_df(linkage_df)
  LG_hom_stack <- test_LG_hom_stack(LG_hom_stack)
  network.layout <- match.arg(network.layout)
  
  LG_hom_stack_sub <- LG_hom_stack[LG_hom_stack$LG == LG,]
  markers <- LG_hom_stack_sub$SxN_Marker
  SN_sub <-
    linkage_df[linkage_df$marker_a %in% markers &
                 linkage_df$marker_b %in% markers,]
  edges <- SN_sub[SN_sub$LOD >= LOD_sequence[1]
                  ,
                  c("marker_a", "marker_b")]
  
  write(paste("Total number of edges:", nrow(edges)), stdout())
  
  # make a network (graph)
  nw <- igraph::graph.data.frame(edges, directed = F)
  gcl <- igraph::groups(igraph::clusters(nw))
  len_init_gcl <- length(gcl)
  size_init_gcl <- lengths(gcl)
  
  if (modify_LG_hom_stack) {
    groups <- stack(gcl)
    colnames(groups) <- c("SxN_Marker", "homologue")
    if (!is.null(nclust_out)) {
      t <- table(groups$homologue)
      topt <- t[order(t, decreasing = T)][1:nclust_out]
      groups <- groups[groups$homologue %in% names(topt),]
      groups$homologue <- as.factor(as.character(groups$homologue))
      levels(groups$homologue) <- c(1:nclust_out)
    }
    groups$LG <- LG
    LG_hom_stack <- LG_hom_stack[LG_hom_stack$LG != LG,]
    LG_hom_stack <- rbind(LG_hom_stack, groups)
    LG_hom_stack$homologue <-
      as.factor(as.character(LG_hom_stack$homologue))
  }
  
  
  if (network.layout == "circular") {
    if (!is.null(device))
      device(...)
    
    group_mem <- igraph::clusters(nw)$membership
    l <- igraph::layout.circle(nw, order = order(group_mem))
    rownames(l) <- igraph::V(nw)$name
    names.layout <- igraph::V(nw)$name
    
    for (cut in LOD_sequence) {
      edgesi <- SN_sub[SN_sub$LOD > cut, c("marker_a", "marker_b")]
      check.layout <-
        edgesi$marker_a %in% names.layout &
        edgesi$marker_b %in% names.layout
      if (sum(!check.layout) > 0)
        warning(
          paste(
            "At LOD",
            cut,
            ": the number of clustering markers is higher than layout LOD. You might want to use a lower layout LOD"
          )
        )
      edgesi <- edgesi[check.layout,]
      nwi <- igraph::graph.data.frame(edgesi, directed = F)
      names <- igraph::V(nwi)$name
      igraph::plot.igraph(
        nwi,
        layout = l[names,],
        edge.color = rgb(0, 0, 0, alpha = 0.1),
        edge.width = 2,
        vertex.label = "",
        vertex.size = 0.1
      )
      par(new = T)
    }
    
    n <- length(igraph::V(nw)$name)
    s <- c(1:n) - 1
    alpha <- (s / n) * 2 * pi
    fincox <- label.offset * (cos(alpha))
    fincoy <- label.offset * (sin(alpha))
    g <- group_mem
    rot <- alpha * 180 / pi
    labnames <- names(g[order(g)])
    cols <- rainbow(length(unique(g)))
    textcols <- cols[g[labnames]]
    
    for (i in seq(n)) {
      text(
        fincox[i],
        fincoy[i],
        labnames[i],
        adj = 0,
        col = textcols[i],
        srt = rot[i],
        cex = cex.lab
      )
    }
    legend(
      -2,
      1.5,
      legend = paste("hom", 1:len_init_gcl, "n=", size_init_gcl),
      pch = 15,
      col = rainbow(len_init_gcl),
      bty = "n"
    )
    text(-2,-1.35, paste("Layout at LOD", LOD_sequence[1]), adj = 0)
    text(-2,-1.5, paste(
      "Ranging from LOD",
      min(LOD_sequence),
      "to",
      max(LOD_sequence)
    ), adj = 0)
    if (!is.null(device))
      dev.off()
  }
  
  
  if (network.layout == "stacked") {
    if (!is.null(device))
      device(...)
    
    if (length(LOD_sequence) > 1) {
      for (cut in LOD_sequence) {
        edgesi <- SN_sub[SN_sub$LOD > cut, c("marker_a", "marker_b")]
        nwi <- igraph::graph.data.frame(edgesi, directed = F)
        g <- igraph::groups(igraph::clusters(nwi))
        gcl <- c(gcl, g)
      }
      
      margin.def <- par("mar")
      par(mar = c(1, 1, 1, 1))
      init_cols <- rainbow(len_init_gcl, alpha = 0.5)
      cols <-
        c(init_cols, rep(rgb(0, 0, 0, alpha = 0.1), length(gcl) - len_init_gcl))
      igraph::plot.igraph(
        nw,
        mark.groups = gcl,
        vertex.size = 2.5,
        vertex.color = "red",
        vertex.label = "",
        mark.col = cols,
        mark.border = "black",
        mark.expand = c(rep(15, len_init_gcl), rep(8, length(gcl) -
                                                     len_init_gcl)),
        edge.color = "orange",
        edge.width = 0.5
      )
      legend(
        "topright",
        legend = paste("hom", 1:len_init_gcl, "n=", size_init_gcl),
        pch = 15,
        col = rainbow(len_init_gcl, alpha = 0.5),
        bty = "n"
      )
      par(mar = margin.def)
      
      if (!is.null(device))
        dev.off()
    }
    
  }
  
  if (!is.null(log)) {
    matc <- match.call()
    write.logheader(matc, log)
  }
  
  if (modify_LG_hom_stack)
    return(LG_hom_stack)
}


#' Cluster 1.0 markers
#' @description \code{cluster_SN_markers} clusters simplex nulliplex at different LOD scores.
#' @param linkage_df A linkage data.frame as output of \code{\link{linkage}} calculating linkage between 1.0 markers.
#' @param LOD_sequence A numeric vector. Specifying a sequence of LOD thresholds at which clustering is performed.
#' @param independence_LOD Logical. Should the LOD of independence be used for clustering? (by default, \code{FALSE}.)
#' @param LG_number Expected number of chromosomes (linkage groups)
#' @param ploidy Ploidy level of the plant species
#' @param parentname Name of parent
#' @param plot_network Logical. Should a network be plotted. Recommended FALSE with large number of marker combinations.
#' @param min_clust_size Integer. The minimum cluster size to be plotted. This does not delete clusters. All clusters are returned.
#' @param plot_clust_size Logical. Should exact cluster size be plotted as vertex labels?
#' @param max_vertex_size Integer. The maximum vertex size. Only used if \code{plot_clust_size=FALSE}.
#' @param min_vertex_size Integer. The minimum vertex size. Only used if \code{plot_clust_size=FALSE}.
#' @param phase_considered Character string. By default all phases are used, but "coupling" or "repulsion" are also allowed.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout (console).
#' @return A list with cluster data.frames.
#' @examples
#' data("SN_SN_P1")
#' cluster_list<-cluster_SN_markers(SN_SN_P1,LOD_sequence=c(4:10),parentname="P1",ploidy=4,LG_number=5)
#' @export
cluster_SN_markers <- function(linkage_df,
                               LOD_sequence = 7,
                               independence_LOD = FALSE,
                               LG_number,
                               ploidy,
                               parentname = "",
                               plot_network = F,
                               min_clust_size = 1,
                               plot_clust_size = TRUE,
                               max_vertex_size = 5,
                               min_vertex_size = 2,
                               phase_considered = "All",
                               log = NULL) {
  linkage_df <- test_linkage_df(linkage_df)
  
  total_marker <-
    unique(c(
      as.character(linkage_df$marker_a),
      as.character(linkage_df$marker_b)
    ))
  
  total_markernr <- length(total_marker)
  
  LODscore <- "LOD"
  if(independence_LOD){
    if(! "LOD_independence" %in% colnames(linkage_df))
      stop("The column LOD_independence should be part of linkage_df when clustering with LOD of independence.
           To obtain the LOD of independence, re-run linkage() with G2_test = TRUE")
    LODscore <- "LOD_independence"
  }
  
  if(phase_considered != "All") {
    if(!phase_considered%in%c("coupling","repulsion")) stop("Unknown phase type")
    linkage_df <- linkage_df[linkage_df[,"phase"]==phase_considered,]
  }
  
  if(max(linkage_df[,LODscore]) < max(LOD_sequence)){
    LOD_sequence <- LOD_sequence[LOD_sequence >= max(linkage_df[,LODscore])]
    if(length(LOD_sequence) == 0) {
      stop(paste("LOD_sequence only contains LOD values higher than present in the data. 
           Clustering only possible up to",max(linkage_df[,LODscore])))
    } else{
      warning(paste("LOD_sequence contains values higher than present in the data. Removing values higher than",
                    max(linkage_df[,LODscore])))
    }
  }
  
  # define edges between SN markers
  edges <- linkage_df[linkage_df[,LODscore] >= min(LOD_sequence),
                      c("marker_a", "marker_b")]
  
  write(paste("Total number of edges:", nrow(edges)), stdout())
  
  # make a network (graph)
  nw <- igraph::graph.data.frame(edges, directed = F)
  gcl <- igraph::groups(igraph::clusters(nw))
  clsizes <- sapply(gcl,length)
  gcl <- gcl[clsizes >= min_clust_size] 
  
  ngroups <- length(gcl)
  groups <- stack(gcl)
  colnames(groups) <- c("marker", "cluster")
  
  # open file for writing if !is.null(log)
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  if (length(LOD_sequence) > 1) {
    groups <- list()
    
    for (i in 1:length(LOD_sequence)) {
      edgesi <- linkage_df[linkage_df[,LODscore] >= LOD_sequence[i],
                           c("marker_a", "marker_b")]
      
      # make a network (graph)
      nwi <- igraph::graph.data.frame(edgesi, directed = F)
      
      # define groups in network
      g <- igraph::groups(igraph::clusters(nwi))
      
      gcl <- c(gcl, g)
      if (LOD_sequence[i] == max(LOD_sequence))
        ngroups <- length(g)
      groupsi <- stack(g)
      colnames(groupsi) <- c("marker", "cluster")
      gname <- as.character(i)
      groups[[gname]] <- groupsi
    }
    
    names(groups) <- LOD_sequence
    groups.out <- groups
    
    nonlinked_markers <- c()
    homologue_numbers <- c()
    for (i in names(groups)) {
      homNr <- length(unique(groups[[i]][, "cluster"]))
      homologue_numbers <- c(homologue_numbers, homNr)
      nonLinked <- total_markernr - nrow(groups[[i]])
      nonlinked_markers <- c(nonlinked_markers, nonLinked)
    }
    
    par(mar = c(5, 5, 2, 5))
    plot(
      LOD_sequence,
      homologue_numbers,
      type = "l",
      col = "black",
      xlab = "LOD score",
      ylab = "Number of clusters",
      cex.lab = 1.2,
      lwd = 2
    )
    segments(par("usr")[1],ploidy*LG_number,
             par("usr")[2],ploidy*LG_number,
             lty=2,lwd=2,col="red2")
    expected <- ifelse(ploidy==2,1,ploidy)*LG_number
    text(x=mean(LOD_sequence),
         y=expected+1,
         paste(expected,"clusters expected"),
         col="red2")
    par(new = TRUE)
    plot(
      LOD_sequence,
      nonlinked_markers,
      type = "l",
      axes = F,
      xlab = NA,
      ylab = NA,
      col = "cadetblue4",
      lwd = 2,
      lty=3
    )
    axis(side = 4,at = seq(from=range(nonlinked_markers)[1],
                           to=range(nonlinked_markers)[2],
                           by=1))
    mtext(
      side = 4,
      line = 3,
      "Number of unlinked markers",
      col = "cadetblue4",
      cex = 1.2
    )
    
    #unique cluster names for all thresholds
    for (i in seq(length(groups))) {
      groups[[i]]$LOD <- names(groups)[i]
      groups[[i]]$unique_cluster <-
        paste(groups[[i]]$cluster, groups[[i]]$LOD, sep =
                "_")
    }
    
    groups.df <- do.call(rbind, groups)
    
    vertex_size <-
      tapply(groups.df$marker, groups.df$unique_cluster, length)
    
    vertex_size <- vertex_size[vertex_size >= min_clust_size]
    
    groups.df <-
      groups.df[groups.df$unique_cluster %in% names(vertex_size),]
    
    edges <-
      tapply(groups.df$unique_cluster, groups.df$marker, function(x) {
        cbind(x[-length(x)], x[-1])
      })
    
    # rbind edgelist
    edges <- do.call(rbind, edges)
    
    # only unique rows
    edges <- unique(edges)
    
    nw.cl <- igraph::graph.data.frame(edges, directed = FALSE)
    
    minLOD <- min(LOD_sequence)
    
    roots <- unique(groups[[as.character(minLOD)]]$unique_cluster)
    roots <- roots[roots %in% igraph::V(nw.cl)$name]
    
    vertex_order <- igraph::V(nw.cl)
    
    if (plot_clust_size) {
      vertex_label <- vertex_size[vertex_order$name]
      vertex_color = "white"
      vertex_size = max_vertex_size
    } else {
      vertex_label <- ""
      vertex_color <- "red"
      sizes <- vertex_size[vertex_order$name]
      size_range <- max_vertex_size - min_vertex_size
      max_cluster_size <- max(sizes)
      vertex_size <-
        min_vertex_size + (sizes / max_cluster_size * size_range)
    }
    par(mar=c(2,4.1,1.5,2))
    plot(
      0,
      xlim = c(-1, 1),
      ylim = c(-1, 1),
      xaxt = "n",
      yaxt = "n",
      bty = "n",
      type = "n",
      xlab = "",
      ylab = "LOD"
    )
    nLOD <- length(LOD_sequence)
    y.coord <- c(0, cumsum(rep(2 / (nLOD - 1), nLOD - 1)))
    
    axis(2, 1 - y.coord, labels = LOD_sequence)
    for (i in y.coord) {
      lines(c(-1, 1), c(1 - i, 1 - i), col = "grey", lty = 2)
    }
    
    l <- igraph::layout_as_tree(nw.cl, root = roots)
    igraph::plot.igraph(
      nw.cl,
      layout = l,
      add = TRUE,
      vertex.size = vertex_size,
      vertex.label.cex = 0.7,
      vertex.label.family = "sans",
      vertex.color = vertex_color,
      edge.color = "black",
      edge.arrow.size = 0.5,
      vertex.label = vertex_label,
      asp = 1,
      vertex.frame.color = "white",
      margin = 0
    )
    
    # order groups.out, so it has the same order as the layout
    dimnames(l) <- list(vertex_order$name, c("x", "y"))
    
    l[,"y"] <- max(l[,"y"]) - l[,"y"]
    l <- l[order(l[,"x"], l[,"y"]),]
    
    for(i in seq(length(groups.out))){
      
      cluster_names <- rownames(l[l[,"y"] == i-1,,drop = FALSE])
      m <- matrix(unlist(strsplit(cluster_names, split = "_")), byrow = TRUE, ncol = 2)
      cluster_order <- m[,1]
      
      clmark <- groups.df[groups.df$LOD == LOD_sequence[i],"marker"]
      groups.out[[i]] <- groups.out[[i]][groups.out[[i]]$marker %in% clmark,]
      
      ps <- as.character(groups.out[[i]][,"cluster"])
      psm <- rep(NA, length(ps))
      for(j in seq(length(cluster_order))){
        psm[ps == cluster_order[j]] <- j
      }
      groups.out[[i]][,"cluster"] <- as.factor(psm)
      groups.out[[i]] <- groups.out[[i]][order(psm),]
      class(groups.out[[i]]) <- c(class(groups.out[[i]]), "cluster_stack")
    }
    
  } else {
    nonLinked <- total_markernr - nrow(groups)
    write(
      paste(
        "Number of clusters at LOD >",
        LOD_sequence,
        ":",
        ngroups,
        "\n",
        nonLinked,
        "unlinked marker(s)",
        sep = " "
      ),
      stdout()
    )
    groups.out <- list(groups)
    names(groups.out) <- LOD_sequence
  }
  
  if (plot_network) {
    # plot network and groups
    igraph::plot.igraph(
      nw,
      layout = igraph::layout_with_graphopt,
      mark.groups = gcl,
      vertex.size = 2.5,
      vertex.color = "red",
      vertex.label = "",
      mark.col = rgb(0, 0, 0, alpha = 0.01),
      mark.border = "black",
      mark.expand = 8,
      edge.color = "orange",
      edge.width = 0.5,
      main = paste0(
        ngroups,
        " SN x SN clusters at LOD score >",
        max(LOD_sequence),
        " of ",
        parentname
      )
    )
  }
  
  # Indicate whether all chromosomes / homologues are (potentially) separated
  if (ngroups < expected) {
    message("Please note: the number of clusters was smaller than expected.")
  }
  
  if (!is.null(log))
    close(log.conn)
  
  return(groups.out)
}


#' Compare linkage maps, showing links between connecting markers common to neighbouring maps
#' @description This function allows the visualisation of connections between different maps, showing them side by side. 
#' @param maplist A list of maps. This is probably most conveniently built on-the-fly in the function call itself.
#' If names are assigned to different maps (list items) these will appear above
#' the maps. In cases of multiple comparisons, for example comparing 1 map of interest to 3 others, the map of interest can
#' be supplied multiple times in the list, interspersed between the other maps. See the example below for details.
#' @param chm.wd The width in inches that linkage groups should be drawn. By default 0.2 inches is used.
#' @param bg.col The background colour of the maps, by default white. It can be useful to use a different background colour for the maps.
#' In this case, supply \code{bg.col} as a vector of colour identifiers, with the same length as \code{maplist} and corresponding to its elements in 
#' the same order. See the example below for details.
#' @param links.col The colour with which links between maps are drawn, by default grey.
#' @param thin.links Option to thin the plotting of links between maps, which might be useful if there are very many shared markers in a 
#' small genetic region. By default \code{NULL}, otherwise supply a value (in cM) for the minimum genetic distance between linking-lines 
#' (e.g. 0.5).
#' @param type Plot type, by default "karyotype". If "scatter" is requested a scatter plot is drawn, but only if the comparison is between 2 maps.
#' @param \dots option to supply arguments to the \code{plot} function (e.g. \code{main =} to add a title to the plot)
#' @return \code{NULL}
#' @examples 
#' data("map1","map2","map3")
#' compare_maps(maplist=list("1a"=map1,"c08"=map2,"1b"=map3),bg.col=c("thistle","white","skyblue"))
#' @export
compare_maps <- function(maplist,
                         chm.wd = 0.2,
                         bg.col = "white",
                         links.col = "grey42",
                         thin.links = NULL,
                         type = "karyotype",
                         ...){
  
  ## Check input:
  if(!is.list(maplist)) stop("Input maplist must be a list!")
  if(length(maplist) < 2 ) stop("Input maplist must contain at least 2 maps for comparison!")
  if(any(!sapply(maplist, is.data.frame))) stop("Input maps must be data frames!")
  if(any(!sapply(maplist, function(map) all(c("marker","position") %in% colnames(map))))) 
    stop("Colnames marker and position are required in input maps!")
  type <- match.arg(type, c("karyotype", "scatter"))
  
  if(type == "scatter" & length(maplist) != 2) stop("Currently scatter plot option is only available for comparison of two maps.")
  
  ##As links are drawn between neighbouring maps only, generate a list of common markers between neighbouring maps
  common_marks <- lapply(1:(length(maplist) - 1), function(n) 
    intersect(as.character(maplist[[n]]$marker),as.character(maplist[[n + 1]]$marker)))
  
  if(any(sapply(common_marks,length) == 0)) stop(paste0("Input had the following number of linking markers: ",
                                                        paste0(sapply(common_marks,length), collapse = ", "),
                                                        ".\nWithout linking markers, a comparison is impossible"))
  
  maxy <- max(sapply(maplist, function(x) max(x$position)))
  
  if(type == "karyotype"){
    
    ## Generalise bg.col into a vector to handle vector input on bg.col:
    if(length(bg.col) != 1){
      if(length(bg.col) != length(maplist)) stop("If supplying multiple background colours, please specify only the required number (same number of elements as maplist)!")
    } else{
      bg.col <- rep(bg.col, length(maplist))
    }
    
    ## Set up the plot area:
    plot(NULL,xlim = c(0,2*length(maplist)), 
         ylim = c(1.1*maxy,-0.1*maxy),
         xlab = "", ylab = "cM", axes = FALSE, ...)
    
    axis(2)
    
    ## Draw chromosome outlines:
    for(i in seq(length(maplist))){
      symbols(x = 2*i - 1, y = min(maplist[[i]]$position), circles= chm.wd, bg=bg.col[i], add=TRUE, inches = FALSE)
      symbols(x = 2*i - 1, y = max(maplist[[i]]$position), circles= chm.wd, bg=bg.col[i], add=TRUE, inches = FALSE)
      rect(2*i - 1 - chm.wd, min(maplist[[i]]$position), 2*i - 1 + chm.wd, max(maplist[[i]]$position), col = bg.col[i])
      
      ## Add the rest of the marker positions:
      for(j in 1:nrow(maplist[[i]])){
        segments(x0 = 2*i - 1 - chm.wd,y0 = maplist[[i]]$position[j],
                 x1 = 2*i - 1 + chm.wd,y1 = maplist[[i]]$position[j])
      }
      
      ## Add chromosome names:
      if(!is.null(names(maplist)[i])) {
        par("xpd" = NA)
        text(x = 2*i - 1, y = par("usr")[4],
             labels = names(maplist)[i], font = 2)
        par("xpd" = FALSE)
      }
      
    }
    
    ## Add links between markers
    for(i in seq(length(common_marks))){
      common.mrks <- common_marks[[i]]
      
      if(!is.null(thin.links)){
        ## thin out the links 
        mp1.pos <- maplist[[i]][match(common.mrks,maplist[[i]]$marker),"position"]
        mp2.pos <- maplist[[i+1]][match(common.mrks,maplist[[i+1]]$marker),"position"]
        
        mp1.range <- range(mp1.pos)
        mp2.range <- range(mp2.pos)
        
        mp1.ints <- seq(mp1.range[1],mp1.range[2],thin.links)
        mp2.ints <- seq(mp2.range[1],mp2.range[2],thin.links)
        
        mp1.dist <- findInterval(maplist[[i]][maplist[[i]]$marker %in% common.mrks,"position"],mp1.ints)
        mp2.dist <- findInterval(maplist[[i+1]][maplist[[i+1]]$marker %in% common.mrks,"position"],mp2.ints)
        
        common.mrks <- union(common.mrks[!duplicated(mp1.dist)],
                             common.mrks[!duplicated(mp2.dist)])
        
      }
      
      for(mark in common.mrks){
        segments(x0 = 2*i - 1 + chm.wd,y0 = maplist[[i]][maplist[[i]]$marker == mark,"position"],
                 x1 = 2*i + 1 - chm.wd,y1 = maplist[[i+1]][maplist[[i+1]]$marker == mark,"position"],
                 col = links.col, lty = 3)
      }
    }
    
  } else{ #scatter plot instead
    
    common.mrks <- common_marks[[1]]
    
    mp1.pos <- maplist[[1]][match(common.mrks,maplist[[1]]$marker),"position"]
    mp2.pos <- maplist[[2]][match(common.mrks,maplist[[2]]$marker),"position"]
    
    plot(mp1.pos,mp2.pos,xlab = names(maplist)[1],
         ylab = names(maplist)[2],
         ...)
    
  }
  return(NULL)
} #compare_maps

#' Consensus LG assignment
#' @description Assign markers to an LG based on consensus between two parents.
#' @param P1_assigned A marker assignment file of the first parent. Should contain the number of linkages per LG per marker.
#' @param P2_assigned A marker assignment file of the second parent. Should be the same markertype as first parent and contain the number of linkages per LG per marker.
#' @param LG_number Number of linkage groups (chromosomes).
#' @param ploidy Ploidy level of plant species.
#' @param consensus_file Filename of consensus output. No output is written if NULL.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return Returns a list containing the following components:
#' \item{P1_assigned}{
#'   A (modified) marker assignment matrix of the first parent.
#' }
#' \item{P2_assigned}{
#'   A (modified) marker assignment matrix of the second parent.
#' }
#' @examples
#' data("P1_SxS_Assigned", "P2_SxS_Assigned_2")
#' SxS_Assigned_list <- consensus_LG_assignment(P1_SxS_Assigned,P2_SxS_Assigned_2,5,4)
#' @export
consensus_LG_assignment <- function(P1_assigned,
                                    P2_assigned,
                                    LG_number,
                                    ploidy,
                                    consensus_file = NULL,
                                    log = NULL) {
  intersect.markers <-
    intersect(rownames(P1_assigned), rownames(P2_assigned))
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  write(paste("Files have", length(intersect.markers), "markers in common"), log.conn)
  LG_col <- c(paste0("LG", 1:LG_number))
  P1_LG <- P1_assigned[intersect.markers, LG_col]
  P2_LG <- P2_assigned[intersect.markers, LG_col]
  cons_LG <- P1_LG + P2_LG
  Assigned_LG_new <- apply(cons_LG, 1,
                           function(x) {
                             m <- max(x, na.rm = T)
                             a <- m / x < 2
                             if (sum(a, na.rm = T) > 1)
                               return(NA)
                             return(which.max(x))
                           })
  write(paste(
    "There are",
    sum(is.na(Assigned_LG_new)),
    "markers with ambiguous assignments."
  ),
  log.conn)
  P1_assigned[intersect.markers, "Assigned_LG"] <- Assigned_LG_new
  P2_assigned[intersect.markers, "Assigned_LG"] <- Assigned_LG_new
  if (!is.null(consensus_file)) {
    Assigned_LG <- Assigned_LG_new
    cons_LG <- cbind(Assigned_LG, cons_LG)
    write.table(cons_LG, consensus_file)
  }
  
  if (!is.null(log))
    close(log.conn)
  
  return(list(P1_assigned = P1_assigned, P2_assigned = P2_assigned))
}

#' Find consensus linkage group names
#' @description Chromosomes that should have same number, might have gotten different numbers between parents during clustering.
#' \code{consensus_LG_names} uses markers present in both parents (usually 1.1 markers) to modify the linkage group numbers in one parent with the other as template
#' @param modify_LG A \code{data.frame} with markernames, linkage group (\code{"LG"}) and homologue (\code{"homologue"}), in which the linkage group numbers will be modified
#' @param template_SxS A file with assigned markers of which (at least) part is present in both parents of the template parent.
#' @param modify_SxS A file with assigned markers of which (at least) part is present in both parents of the parent of which linkage group number are modified.
#' @param merge_LGs Logical, by default \code{TRUE}. If \code{FALSE}, any discrepency in the number of linkage groups will not be merged, but removed instead.
#' This can be needed if the number of chromosomes identified is not equal between parents, and the user wishes to proceed with a core set.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return A modified modified_LG according to the template_SxS linkage group numbering
#' @examples
#' data("LGHomDf_P2_2", "P1_SxS_Assigned", "P2_SxS_Assigned")
#' consensus_LGHomDf<-consensus_LG_names(LGHomDf_P2_2, P1_SxS_Assigned, P2_SxS_Assigned)
#' @export
consensus_LG_names <- function(modify_LG,
                               template_SxS,
                               modify_SxS,
                               merge_LGs = TRUE,
                               log = NULL) {
  # might want to incorporate some warning messages in this function
  # in order to check whether lg assignment is clear
  # find linkage groups of common markers and make a table
  modify_LG <- test_LG_hom_stack(modify_LG)
  
  common_markers <-
    intersect(rownames(template_SxS), rownames(modify_SxS))
  LG_P1 <- template_SxS[common_markers, "Assigned_LG"]
  LG_P2 <- modify_SxS[common_markers, "Assigned_LG"]
  
  cons_table <- table(LG_P1, LG_P2)
  
  eq.chms <- TRUE
  
  if(nrow(cons_table) < ncol(cons_table)) { #A merge problem arises when modify has more chms than template. Warning about this feature.
    warning(paste("modify_LG parent possesses more chromosomes than template parent. Extras will be",
                  ifelse(merge_LGs,"merged.","dropped.")))
    eq.chms <- FALSE
  }
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  write("\n####Original LG names\n", log.conn)
  write(knitr::kable(cons_table,
                     row.names = TRUE), log.conn)
  
  if(!eq.chms & !merge_LGs){
    ## Remove excess cols from cons_table:
    savecols <- sort(apply(cons_table, 1, which.max))
    remcols <- setdiff(1:ncol(cons_table), savecols)
    
    modify_LG <- modify_LG[-which(modify_LG$LG %in% levels(modify_LG$LG)[remcols]),]
    modify_LG$LG <- droplevels(modify_LG$LG) #drop these levels also
    
    cons_table <- cons_table[,savecols]
  }
  
  # get the clusters that have most assignments in P1
  changed_LGs <- apply(cons_table, 2, which.max)
  
  # Make sure there are no conflicts here - if there are, resolve these also:
  counter <- 1
  
  while(any(duplicated(changed_LGs)) & counter < 10){
    if(counter == 1)
      warning(paste0("Cannot simply rename LG ",
                     paste0(changed_LGs[duplicated(changed_LGs)],collapse = ", "),
                     ". Attempting a resolution...."))
    
    for(i in which(duplicated(changed_LGs))){
      
      #First determine which is the duplicate (if there are multiple, we are in trouble!)
      max.counts <- apply(cons_table[,which(changed_LGs == changed_LGs[i])],2,max)
      
      ## Duplicated is the one with the least counts. If there is a tie we are in trouble...
      if(length(max.counts) != length(unique(max.counts)))
        warning(paste("Cannot decide a duplicate concerning LG",changed_LGs[i]))
      
      i <- names(which.min(max.counts)) #this allows to change to a different duplicate within the loop. Possibly buggy..
      
      x <- sort(cons_table[,i],decreasing = T) > 0
      
      for(nam in as.numeric(names(x[x])[2:length(x[x])])){
        if(!nam %in% changed_LGs){
          changed_LGs[i] <- nam
          break()
        }
      }
      
    }
    counter <- counter + 1
  }
  
  if(counter == 10 | any(is.na(changed_LGs))) stop("Unable to resolve re-numbering puzzle. SxN cluster analysis may need to be revisited!")
  rm(counter)
  
  modify_LG$LG <- as.factor(modify_LG$LG)
  levels(modify_LG$LG) <- as.numeric(changed_LGs)
  
  # level ordering
  modify_LG$LG <-
    as.factor(as.numeric(as.character(modify_LG$LG)))
  
  
  colnames(cons_table) <- changed_LGs
  write("\n####Modified LG names\n", log.conn)
  reorder_cons_table <-
    cons_table[, order(as.numeric(colnames(cons_table))), drop = FALSE]
  write(knitr::kable(reorder_cons_table,
                     row.names = TRUE), log.conn)
  
  if (!is.null(log))
    close(log.conn)
  
  return(modify_LG)
}


#' Convert marker dosages to the basic types.
#' @description Convert marker dosages to the basic types which hold the same information and for which linkage calculations can be performed.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @param ploidy ploidy level of the plant species. If parents have different ploidy level, ploidy of parent1.
#' @param ploidy2 ploidy level of the second parent. NULL if both parents have the same ploidy level.
#' @param parent1 Character string specifying the first (usually maternal) parentname.
#' @param parent2 Character string specifying the second (usually paternal) parentname.
#' @param marker_conversion_info Logical, by default \code{FALSE}. Should marker conversion information be returned? This output can be useful for later map phasing step,
#' if original marker coding is desired (which is most likely the case).
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return A modified dosage matrix. If \code{marker_conversion_info = TRUE}, this function returns a list, with both the converted dosage_matrix, and 
#' information on the marker conversions performed per marker.
#' @examples
#' data("ALL_dosages")
#' conv<-convert_marker_dosages(dosage_matrix=ALL_dosages, ploidy = 4)
#' @export
convert_marker_dosages <- function(dosage_matrix,
                                   ploidy,
                                   ploidy2 = NULL,
                                   parent1 = "P1",
                                   parent2 = "P2",
                                   marker_conversion_info = FALSE,
                                   log = NULL) {
  dosage_matrix <- test_dosage_matrix(dosage_matrix)
  # data checks
  unexpected_dosages <- dosage_matrix > max(c(ploidy, ploidy2)) | dosage_matrix < 0
  sum_unexpected <- sum(unexpected_dosages, na.rm = TRUE)
  if (sum_unexpected != 0) {
    warning(
      paste("There are dosages greater than",ploidy,"or less than 0 in the dataset.
            If they include parental dosages, markers are removed from the output.
            Otherwise the dosage is made missing.")
    )
    dosage_matrix[unexpected_dosages] <- NA
    rm_markers <-
      apply(dosage_matrix[,c(parent1, parent2)], 1, function(x)
        any(is.na(x)))
    dosage_matrix <- dosage_matrix[!rm_markers,]
  }
  
  if(is.null(ploidy2)) ploidy2 <- ploidy
  
  prog_ploidy <- (ploidy + ploidy2)/2
  
  # data definition
  marker <- rownames(dosage_matrix)
  par_orig <- par <- dosage_matrix[, c(parent1, parent2)]
  progeny <-
    which(!colnames(dosage_matrix) %in% c(parent1, parent2))
  
  # sum and difference parents (for definition convertable markers)
  sumpar <- rowSums(par)
  difpar <- par[, parent1] - par[, parent2]
  absdifpar <- abs(difpar)
  
  # find homozygous
  homoP1 <- par[,parent1] == 0 | par[,parent1] == ploidy
  homoP2 <- par[,parent2] == 0 | par[,parent2] == ploidy2
  onehomo <- apply(cbind(homoP1, homoP2), 1, function(x)
    sum(x) == 1)
  
  mincovP1 <- par[,parent2] == 0 & par[,parent1] > ploidy/2 & par[,parent1] < ploidy
  mincovP2 <- par[,parent1] == 0 & par[,parent2] > ploidy2/2 & par[,parent2] < ploidy2
  
  P1_conversion <- par[,parent1] == ploidy & par[,parent2] < ploidy2/2
  P2_conversion <- par[,parent2] == ploidy2 & par[,parent1] < ploidy/2
  
  onepar_conversion <- P1_conversion | P2_conversion
  
  minus_conversion <- sumpar > prog_ploidy & !onepar_conversion | (mincovP1 | mincovP2)
  
  # unique(apply(dosage_matrix[minus_conversion, c(parent1,parent2)],1,function(x) paste0(x, collapse="")))
  
  dosage_matrix[minus_conversion, parent1] <-
    ploidy - dosage_matrix[minus_conversion, parent1]
  dosage_matrix[minus_conversion, parent2] <-
    ploidy2 - dosage_matrix[minus_conversion, parent2]
  dosage_matrix[minus_conversion,progeny] <-
    prog_ploidy - dosage_matrix[minus_conversion,progeny]
  
  par <- dosage_matrix[, c(parent1, parent2)]
  
  # find markers of which one parent can be converted
  P1_conversion <- par[,parent1] == ploidy & par[,parent2] < ploidy2/2
  P2_conversion <- par[,parent2] == ploidy2 & par[,parent1] < ploidy/2
  onepar_conversion <- P1_conversion | P2_conversion
  
  # find markers of which both parents can be converted
  
  # do onepar conversions parents
  par[P1_conversion, parent1] <- ploidy - par[P1_conversion, parent1]
  par[P2_conversion, parent2] <- ploidy2 - par[P2_conversion, parent2]
  
  # define parent corrections
  corr_P <- dosage_matrix[, c(parent1, parent2)] - par
  corr_P <- rowSums(corr_P)
  corr <- corr_P > 0
  
  # replace parental genotypes
  dosage_matrix[, c(parent1, parent2)] <- par
  
  # do ploidy minus conversion
  if(sum(corr) > 0){
    # correct progeny with 0.5*dosage correction
    dosage_matrix[corr, progeny] <-
      dosage_matrix[corr, progeny] - (0.5 * corr_P[corr])
    
    # no negative values (in case of double reduction)
    # dosage_matrix[corr, progeny] <-
    #   apply(dosage_matrix[corr, progeny, drop = FALSE], 1, function(x)
    #     max(x, na.rm = T)) - dosage_matrix[corr, progeny]
  }
  dosage_matrix[dosage_matrix < 0] <- NA
  
  # convert palindrome markers
  if(ploidy == ploidy2){
    palindrome_markers <-
      (dosage_matrix[, parent1] > dosage_matrix[, parent2]) &
      (abs(dosage_matrix[, parent1] - (0.5 * ploidy)) == abs(dosage_matrix[, parent2] -
                                                               (0.5 * ploidy)))
    dosage_matrix[palindrome_markers,] <-
      ploidy - dosage_matrix[palindrome_markers,]
  }
  
  # get stats of marker conversions
  par <- dosage_matrix[, c(parent1, parent2)]
  # minus_conversion
  par_swapped <- par != par_orig
  #Bug May 2020: these were incorrectly recorded as having only single parent conversion.
  par_swapped[minus_conversion & !onepar_conversion,] <- TRUE 
  par_swapped_num <- apply(par_swapped, 1, sum)
  
  # find non-segregating markers
  no_seg <- homoP1 & homoP2
  
  rownames(dosage_matrix) <- marker
  
  # remove non-segragating markers
  dosage_matrix <- dosage_matrix[!no_seg,]
  
  marker <- rownames(dosage_matrix) #Record rownames again..
  
  ## Identify any impossible dosage scores (excluding double reduction) and make missing:
  p1split <- t(apply(dosage_matrix[,parent1,drop=FALSE],1,function(x) c(rep(1,x),rep(0,ploidy-x))))
  p2split <- t(apply(dosage_matrix[,parent2,drop=FALSE],1,function(x) c(rep(1,x),rep(0,ploidy2-x))))
  
  possible.dosages <- lapply(1:nrow(p1split),
                             function(r) unique(rowSums(expand.grid(
                               unique(rowSums(expand.grid(p1split[r,], p1split[r,]))),
                               unique(rowSums(expand.grid(p2split[r,], p2split[r,])))
                             )
                             )
                             )
  )
  
  dosage_matrix <- do.call("rbind", lapply(1:nrow(dosage_matrix), function(r) {
    temp.prog <- dosage_matrix[r,progeny]
    temp.prog[!temp.prog %in% possible.dosages[[r]]] <- NA
    c(dosage_matrix[r,c(parent1,parent2)],temp.prog)
  }
  )
  )
  
  rownames(dosage_matrix) <- marker #Re-assign rownames...
  
  parental_info <-
    table(dosage_matrix[, parent1], dosage_matrix[, parent2])
  colnames(parental_info) <-
    paste0(parent2, "_", colnames(parental_info))
  rownames(parental_info) <-
    paste0(parent1, "_", rownames(parental_info))
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  write("\n####Marker dosage frequencies:\n",
        file = log.conn)

  write(knitr::kable(parental_info),
        log.conn)

  write(paste0("\nmarkers not converted: ", sum(par_swapped_num == 0)),
        file = log.conn)
  write(paste0("\nmarkers 1 parent converted: ", sum(par_swapped_num == 1)),
        file = log.conn)
  write(paste0("\nmarkers 2 parents converted: ", sum(par_swapped_num == 2)),
        file = log.conn)
  write(paste0("\nnon-segregating markers deleted: ", sum(no_seg)),
        file = log.conn)
  if (!is.null(log))
    close(log.conn)
  
  if (marker_conversion_info)
    return(list(dosage_matrix = dosage_matrix, marker_conversion_info = par_swapped))
  return(dosage_matrix)
  
} #convert_marker_dosages

#' Convert (probabilistic) genotype calling results from polyRAD to input compatible with polymapR
#' @param RADdata An RADdata (S3 class) object; output of the function \link[polyRAD]{PipelineMapping2Parents} having followed
#' the prior steps needed in the polyRAD pipeline. See the polyRAD vignette for details.
#' @return A data frame which include columns: 
#' MarkerName, SampleName,P0 ~ Pploidy (e.g. P0 ~ P4 for tetraploid, which represents
#' the probability assigning to this dosage), maxgeno (the most likely dosage),
#' and maxP (the maximum probability)
#' @examples
#' data("exampleRAD_mapping")
#' convert_polyRAD(RADdata = exampleRAD_mapping)
#' @export
convert_polyRAD <- function(RADdata){
  genotype_tmp <- RADdata$genotypeLikelihood[[1]]
  loci <- dimnames(genotype_tmp)[[3]] 
  output <- do.call(rbind,lapply(loci, function(l){
    tmp <- t(genotype_tmp[,,l])
    colnames(tmp) <- paste0('P',colnames(tmp))
    tmp_new <- data.frame("MarkerName" = replicate(nrow(tmp),l),
                          "SampleName" = rownames(tmp),
                          tmp,
                          'maxgeno' = as.numeric(apply(tmp,1,which.max)) - 1, #class - 1 become the dosage
                          "maxP" = as.numeric(apply(tmp,1,max)))
    rownames(tmp_new) <- NULL
    tmp_new
  }))
  return(output)
} #convert_polyRAD

#' Convert (probabilistic) genotype calling results from updog to input compatible with polymapR. 
#' @param mout An object of class multidog; output of the function \link[updog]{multidog}.
#' @return A data frame which include columns: 
#' MarkerName, SampleName,P0 ~ Pploidy (e.g. P0 ~ P4 for tetraploid, which represents
#' the probability assigning to this dosage), maxgeno (the most likely dosage),
#' and maxP (the maximum probability)
#' @examples
#' data("mout")
#' convert_updog(mout)
#' @export
convert_updog <- function(mout){
  probaset <- mout$inddf[,grep("Pr_",colnames(mout$inddf))]
  ploidy <- mout$snpdf[1,"ploidy"]
  colnames(probaset) <- c(paste0("P",seq(0,ploidy)))
  
  output <- data.frame("MarkerName" = mout$inddf[,"snp"],
                       "SampleName" = mout$inddf[,"ind"],
                       probaset,
                       "maxgeno" = mout$inddf[,"geno"],#fout$geno,
                       "maxP" = mout$inddf[,"maxpostprob"])#fout$maxpostprob)
  return(output)
} #convert_updog

#'@title Check if dosage scores may have to be shifted
#'@description fitPoly sometimes uses a "shifted" model to assign dosage
#'scores (e.g. all samples are assigned a dosage one higher than the true
#'dosage). This happens mostly when there are only few dosages present
#'among the samples. This function checks if a shift of +/-1 is possible.
#'@usage correctDosages(chk, dosage_matrix, parent1, parent2, ploidy,
#'polysomic=TRUE, disomic=FALSE, mixed=FALSE,
#'absent.threshold=0.04)
#'@param chk data frame returned by function checkF1 when called without
#'shiftmarkers
#'@param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#'@param parent1 character vector with names of the samples of parent 1
#'@param parent2 character vector with names of the samples of parent 2
#'@param ploidy ploidy of parents and F1 (correctDosages must not be used for
#'F1 populations where the parents have a different ploidy, or where the
#'parental genotypes are not scored together with the F1);
#'same as used in the call to checkF1 that generated data.frame chk
#'@param polysomic if TRUE at least all polysomic segtypes are considered;
#'if FALSE these are not specifically selected (but if e.g. disomic is TRUE,
#'any polysomic segtypes that are also disomic will still be considered);
#'same as used in the call to checkF1 that generated data.frame chk
#'@param disomic if TRUE at least all disomic segtypes are considered (see
#'param polysomic); same as used in the call to checkF1 that generated
#'data.frame chk
#'@param mixed if TRUE at least all mixed segtypes are considered (see
#'param polysomic). A mixed segtype occurs when inheritance in one parent is
#'polysomic (random chromosome pairing) and in the other parent disomic (fully
#'preferential chromosome pairing); same as used in the call to checkF1 that
#'generated data.frame chk
#'@param absent.threshold the threshold for the fraction of ALL samples
#'that has the dosage that is assumed to be absent due to mis-fitting of
#'fitPoly; should be at least the assumed error rate of the fitPoly scoring
#'assuming the fitted model is correct
#'
#'@return a data frame with columns
#'\itemize{
#'\item{markername}
#'\item{segtype: the bestfit (not bestParentfit!) segtype from chk}
#'\item{parent1, parent2: the consensus parental dosages; possibly
#'low-confidence, so may be different from those reported in chk}
#'\item{shift: -1, 0 or 1: the amount by which this marker should be shifted}
#'}
#'The next fields are only calculated if shift is not 0:
#'\itemize{
#'\item{fracNotOk: the fraction of ALL samples that are in the dosage
#'(0 or ploidy) that should be empty if the marker is indeed shifted.}
#'\item{parNA: the number of parental dosages that is missing (0, 1 or 2)}
#'}
#'@details A shift of -1 (or +1) is proposed when (1) the fraction of all
#'samples with dosage 0 (or ploidy) is below absent.threshold, (2) the
#'bestfit (not bestParentfit!) segtype in chk has one empty dosage on the
#'low (or high) side and more than one empty dosage at the high (or low) side,
#'and (3) the shifted consensus parental dosages do not conflict with the
#'shifted segregation type.\cr
#'The returned data.frame (or a subset, e.g. based on the values in the
#'fracNotOk and parNA columns) can serve as parameter shiftmarkers in a
#'new call to checkF1.\cr
#'Based on the quality scores assigned by checkF1 to
#'the original and shifted versions of each marker the user can decide if
#'either or both should be kept. A data.frame combining selected rows
#'of the original and shifted versions of the checkF1 output (which may
#'contain both a shifted and an unshifted version of some markers) can then be
#'used as input to compareProbes or writeDosagefile.
#'@export
correctDosages <- function(chk,
                           dosage_matrix,
                           parent1 = "P1",
                           parent2 = "P2",
                           ploidy = 4,
                           polysomic = TRUE,
                           disomic = FALSE,
                           mixed = FALSE,
                           absent.threshold = 0.04) {
  #This function identifies markers that are probably misscored by fitPoly
  #(bestfit in 1_1, 1_3, 11_1, 11_2 (if ploidy==4) and parents in same dosage
  #classes;
  # note that 1_1 and 1_3 can only occur if disomic=TRUE in check.F1)
  #and produce data to suggests a revised scoring for F1 and parent, based on
  #the segregation in an F1 population and the scores in all other samples
  if (ploidy %% 2 != 0)
    stop("correctDosages: only even ploidy allowed (parents must have been genotyped together with F1")
  if (sum(is.na(match(c(parent1, parent2), names(chk)))) > 0)
    stop("correctDosages: not all parental samples in names of chk")
  if ("shift" %in% names(chk))
    stop("correctDosages: chk must be calculated without shift")
  seginfo <- segtypeInfoSummary(calcSegtypeInfo(ploidy)) #seginfo <- polymapR:::segtypeInfoSummary(calcSegtypeInfo(ploidy))
  seginfo <- seginfo[(polysomic & seginfo$par.poly==1) |
                       (disomic & seginfo$par.di==1) |
                       (mixed & seginfo$par.mixed==1),]
  chk <- chk2integer(chk) #chk <- polymapR:::chk2integer(chk)
  if (is.factor(chk$bestfit)) chk$bestfit <- as.character(chk$bestfit)
  NAcol <- which(names(chk) == "F1_NA") # last before parental samples
  
  #calculate a (possibly very low conf) new parental consensus:
  #(the one in chk may be based partly on bestParentfit)
  for (p in 1:2) {
    if (p==1) parent <- parent1 else parent <- parent2
    pc <- which(names(chk) == "parent1") - 1 + p
    chk[,pc] <- rep(NA, nrow(chk))
    if (length(parent) > 0) {
      #for calculating the consensus score (pcons) we need to calculate the
      #parallel min and max over the columns of the parental samples;
      pcons <- apply(chk[,(NAcol+1):(NAcol+length(parent)), drop=FALSE],
                     1, min) #no na.rm=TRUE !!
      rw <- !is.na(pcons) &  pcons ==
        apply(chk[,(NAcol+1):(NAcol+length(parent)), drop=FALSE], 1, max)
      chk[rw, pc] <- pcons[rw] #so the "old" parent1 and parent2 columns now
      #                         have a new type of consensus score
      NAcol <- NAcol + length(parent) #startpoint for parent2
    }
  }
  dosage_matrix <- dosage_matrix[rownames(dosage_matrix) %in% as.character(chk$MarkerName),]
  # scores$MarkerName <- factor(scores$MarkerName) #drop unused levels
  nonNA <- function(x) sum(!is.na(x))
  is0 <- function(x) sum(x==0, na.rm=TRUE)
  isploidy <- function(x) sum(x==ploidy, na.rm=TRUE)
  
  scores.nonNA <- apply(dosage_matrix[,3:ncol(dosage_matrix)],1,nonNA)
  # scores.frc0 <- tapply(scores$geno, scores$MarkerName, FUN=is0)
  scores.frc0 <- apply(dosage_matrix[,3:ncol(dosage_matrix)],1,is0)
  scores.frc0 <- scores.frc0 / scores.nonNA
  # scores.frcploidy <- tapply(scores$geno, scores$MarkerName, FUN=isploidy)
  scores.frcploidy <- apply(dosage_matrix[,3:ncol(dosage_matrix)],1,isploidy)
  scores.frcploidy <- scores.frcploidy / scores.nonNA
  scores.frc0 <- scores.frc0[order(names(scores.frc0))] #needed!
  scores.frcploidy <- scores.frcploidy[order(names(scores.frcploidy))] #needed!
  #these are now in alphabetical order; do the same with the chk data frame:
  ch.ord <- order(chk$MarkerName) #we need this later
  chk <- chk[ch.ord,]
  
  #potential shift value for each marker:
  potShifts <- calcPotShifts(chk$bestfit, ploidy) #potShifts <- polymapR:::calcPotShifts(chk$bestfit, ploidy)
  #for which rows do we need to check a potential shift?
  mf <- !is.na(potShifts$shift) &
    ( (potShifts$shift == -1 & scores.frc0 <= absent.threshold) |
        (potShifts$shift == 1 & scores.frcploidy <= absent.threshold) )
  shfpar <- (cbind(chk$parent1 + potShifts$shift, chk$parent2 + potShifts$shift))[mf,,drop = FALSE]
  #    shifted parental consensus dosages
  shf_segtype <- paste(potShifts[mf ,1],
                       (potShifts[mf ,2] + potShifts[mf ,3]), sep= "_")
  #    shifted segtype: original first dosage + shift
  
  match_p_seg <- function(segtype, p1, p2) {
    # for one segtype: do the parental consensus dosages (p12) match it?
    if (is.na(p1) && is.na(p2)) return(TRUE)
    si <- seginfo[seginfo$segtype==segtype, 3:4] #parental dosage columns
    if (is.na(p1)) return(p2 %in% si[,2])
    #p1 now not NA
    if (is.na(p2)) return(p1 %in% si[,1])
    return(any(si[,1]==p1 & si[,2]==p2))
  }
  
  goodshifts <- mapply(FUN=match_p_seg, shf_segtype, shfpar[,1,drop = FALSE], shfpar[,2,drop = FALSE])
  shiftmf <- rep(0, sum(mf))
  shiftmf[goodshifts] <- potShifts$shift[mf][goodshifts]
  p1namf <- p2namf <- onepokmf <- bothpokmf <- bothpnamf <- rep(NA, sum(mf))
  p1namf[goodshifts] <- is.na(shfpar[goodshifts, 1])
  p2namf[goodshifts] <- is.na(shfpar[goodshifts, 2])
  #sum_isna <- p1namf[goodshifts] + p2namf[goodshifts]
  #bothpokmf[goodshifts] <- sum_isna == 0
  #onepokmf[goodshifts] <- sum_isna == 1
  #bothpnamf[goodshifts] <- sum_isna == 2
  
  nmrk <- nrow(chk)
  P1na = rep(NA, nmrk)
  P2na = rep(NA, nmrk)
  result <- data.frame(
    MarkerName = chk$MarkerName,
    segtype = chk$bestfit,
    parent1 = chk$parent1,
    parent2 = chk$parent2,
    shift = rep(0, nmrk),
    #the next columns will get non-NA values only for accepted shifts
    fracNotOk = rep(NA, nmrk) #fraction of all samples in the forbidden class
    #                           (0 or ploidy)
    #bothPok = rep(NA, nmrk), #TRUE if the combination matches the "misfit" parental scores
    #onePok = rep(NA, nmrk), #TRUE if one P missing and the other has one of the "misfit" parental scores
    #bothPna = rep(NA, nmrk)) #TRUE if both parents NA
  )
  result$shift[mf] <- shiftmf
  result$fracNotOk[result$shift == -1] <- scores.frc0[result$shift == -1]
  result$fracNotOk[result$shift == 1] <- scores.frcploidy[result$shift == 1]
  P1na[mf] <- p1namf
  P2na[mf] <- p2namf
  result$ParNA <- P1na + P2na #0, 1 or 2
  
  #sort result into same order as original chk:
  result <- result[order(ch.ord),]
  invisible(result)
} #correctDosages


#' Create input files for TetraOrigin using an integrated linkage map list and marker dosage matrix
#' @description \code{createTetraOriginInput} is a function for creating an input file for TetraOrigin, combining
#' map positions with marker dosages.
#' @param maplist A list of maps. In the first column marker names and in the second their position.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns. Either provide the unconverted dosages (i.e.
#'  before using the \code{\link{convert_marker_dosages}} function), or converted dosages (i.e. screened data), in matrix form.
#'  The analysis and results are unaffected by this choice, but it may be simpler to understand the results if converted dosages
#'  are used. Conversely, it may be advantageous to use the original unconverted dosages if particular marker alleles are being
#'  tracked for (e.g.) the development of selectable markers afterwards.
#' @param bin_size Numeric. Size (in cM) of the bins to include. If \code{NULL} (by default) then all markers are used (no binning).
#' @param bounds Numeric vector. If \code{NULL} (by default) then all positions are included, however if specified then output
#' is limited to a specific region, which is useful for later fine-mapping work.
#' @param remove_markers Optional vector of marker names to remove from the maps. Default is \code{NULL}.
#' @param outdir Output directory to which input files for TetraOrigin are written.
#' @param output_stem Character prefix to add to the .csv output filename.
#' @param plot_maps Logical. Plot the marker positions of the selected markers using \code{\link{plot_map}}.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @examples
#' \dontrun{
#' data("integrated.maplist","ALL_dosages")
#' createTetraOriginInput(maplist=integrated.maplist,dosage_matrix=ALL_dosages,bin_size=10)}
#' @export
createTetraOriginInput <- function(maplist,
                                   dosage_matrix,
                                   bin_size = NULL,
                                   bounds = NULL,
                                   remove_markers = NULL,
                                   outdir = "TetraOrigin",
                                   output_stem = "TetraOrigin_input",
                                   plot_maps = TRUE,
                                   log = NULL) {
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  dosage_matrix <- test_dosage_matrix(dosage_matrix)
  
  if(!file.exists(outdir)) dir.create(outdir)
  
  if(is.null(names(maplist))){
    if(length(maplist) == 1) stop("maplist provided was an unnamed list of length 1. Chromosome number not found.\n Assign names to maplist using names() function first!")
    warning("maplist provided was an unnamed list. Proceeding with automatic names")
    names(maplist) <- paste0("LG",seq(length(maplist)))
  }
  
  outlist <- lapply(seq(length(maplist)),
                    function(mapn){
                      
                      map <- maplist[[mapn]]
                      
                      chrn <- mapn
                      if(length(maplist) == 1) chrn <- gsub("(.*)([0-9]+)(.*)","\\2", names(maplist)[mapn])
                      
                      if(!is.null(remove_markers)) map <- map[!map$marker%in%remove_markers,]
                      if(!is.null(bounds)){
                        if(length(bounds) != 2 | bounds[2]<=bounds[1]) stop("Incorrect bounds specified. This should be a vector of two distinct cM positions.")
                        map <- map[map[,"position"] >= bounds[1] & map[,"position"] <= bounds[2],]
                      }
                      
                      if(!is.null(bin_size)){
                        bins<-seq(0,ceiling(max(map["position"])) + bin_size,by=bin_size)
                        binned.markers<-lapply(1:(length(bins)-1), function(n) as.character(map[,"marker"])[map[,"position"]>=bins[n] & map[,"position"]<bins[n+1]])
                        
                        remove.marks<-function(merkers){
                          numNA <- sapply(merkers, function(m) length(which(is.na(dosage_matrix[match(m,rownames(dosage_matrix)),3:ncol(dosage_matrix)]))))
                          ParSeg <- as.factor(sapply(merkers, function(m) paste0(dosage_matrix[match(m,rownames(dosage_matrix)),1],
                                                                                 dosage_matrix[match(m,rownames(dosage_matrix)),2]))
                          )
                          whittled<-do.call("c",lapply(levels(ParSeg), function(ps) setdiff(names(ParSeg[ParSeg==ps]),names(which.min(numNA[ParSeg==ps])))
                          ))
                          return(whittled)
                        }
                        
                        removed.markers<-do.call("c",
                                                 lapply(seq(length(binned.markers)),function(l)
                                                   if(length(binned.markers[[l]])>0) remove.marks(binned.markers[[l]]))
                        )
                        
                        map <- map[!map$marker%in%removed.markers,]
                        
                      }
                      
                      if(length(which(map[,"marker"] %in% rownames(dosage_matrix))) != nrow(map)) {
                        warning("Not all mapped markers have corresponding dosages. Attempting to continue with those that do...")
                        map <- map[map[,"marker"] %in% rownames(dosage_matrix),]
                      }
                      
                      if(nrow(map) < 2) stop("Insufficient map information to proceed. Check maplist and/or dosage_matrix.")
                      
                      tetraOrigin_out <- as.data.frame(matrix(0,ncol=nrow(map)+1,nrow=ncol(dosage_matrix)+3))
                      tetraOrigin_out[1,] <- c("marker",as.character(map[,"marker"]))
                      tetraOrigin_out[2,] <- c("chromosome",rep(chrn,nrow(map)))
                      tetraOrigin_out[3,] <- c("position",round(map[,"position"],2))
                      tetraOrigin_out[4:nrow(tetraOrigin_out),1] <- colnames(dosage_matrix)
                      
                      for(r in 1:nrow(map))
                        tetraOrigin_out[4:nrow(tetraOrigin_out),r+1] <- as.numeric(dosage_matrix[match(map[r,"marker"],
                                                                                                       rownames(dosage_matrix)),])
                      
                      write.table(tetraOrigin_out,
                                  file.path(outdir,paste0(output_stem,"_",names(maplist)[mapn],".csv")),
                                  sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)
                      
                      write(paste0("\n",nrow(map)," markers from a possible ",nrow(maplist[[mapn]])," on LG",mapn," written to '",
                                   outdir,"/",output_stem,"_",names(maplist)[mapn],".csv'"),log.conn)
                      
                      return(map)
                    }
  )
  
  if(!is.null(log)) close(log.conn)
  
  if(plot_maps){
    plot_map(maplist=outlist,
             color_by_type = TRUE,
             dosage_matrix = dosage_matrix,
             bg_col = "white")
  }
  
  return(outlist)
}

#' Create a phased homologue map list using the original dosages
#' @description \code{create_phased_maplist} is a function for creating a phased maplist, using
#' integrated map positions and original marker dosages.
#' @param input_type Can be either one of 'discrete' or 'probabilistic'. For the former (default), at least \code{dosage_matrix.conv} must be supplied,
#' while for the latter \code{chk} must be supplied. 
#' @param maplist A list of maps. In the first column marker names and in the second their position.
#' @param dosage_matrix.conv Matrix of marker dosage scores with markers in rows and individuals in columns. Note that dosages must be
#' in converted form, i.e. after having run the \code{\link{convert_marker_dosages}} function. Errors may result otherwise.
#' @param dosage_matrix.orig Optional, by default \code{NULL}.The unconverted dosages (i.e. raw dosage data before using
#' the \code{\link{convert_marker_dosages}} function). Required if \code{original_coding} is \code{TRUE}.
#' @param probgeno_df Probabilistic genotypes, for description see e.g. \code{\link{gp_overview}}. Required if probabilistic genotypes are used.
#' @param chk Output list as returned by function \code{\link{checkF1}}. Required if probabilistic genotypes are used.
#' @param remove_markers Optional vector of marker names to remove from the maps. Default is \code{NULL}.
#' @param original_coding Logical. Should the phased map use the original marker coding or not? By default \code{FALSE}.
#' @param N_linkages Number of significant linkages (as defined in \code{\link{homologue_lg_assignment}}) required for high-confidence linkage group assignment.
#' @param lower_bound Numeric. Lower bound for the rate at which homologue linkages (fraction of total for that marker) are recognised.
#' @param ploidy Integer. Ploidy of the organism.
#' @param ploidy2 Optional integer, by default \code{NULL}. Ploidy of parent 2, if different from parent 1.
#' @param marker_assignment.1 A marker assignment matrix for parent 1 with markernames as rownames and at least containing the column \code{"Assigned_LG"}.
#' @param marker_assignment.2 A marker assignment matrix for parent 2 with markernames as rownames and at least containing the column \code{"Assigned_LG"}.
#' @param parent1 character vector with names of the samples of parent 1
#' @param parent2 character vector with names of the samples of parent 2
#' @param marker_conversion_info One of the list elements generated by the function \code{\link{convert_marker_dosages}}. Required if \code{original_coding} is \code{TRUE}.
#' @param log Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.
#' @param verbose Logical, by default \code{TRUE}. Should details of the phasing process be given?
#' @examples
#' \dontrun{
#' data("integrated.maplist", "screened_data3", "marker_assignments_P1","marker_assignments_P2")
#' create_phased_maplist(maplist = integrated.maplist,
#'                      dosage_matrix.conv = screened_data3,
#'                      marker_assignment.1=marker_assignments_P1,
#'                      marker_assignment.2=marker_assignments_P2,
#'                      ploidy = 4)}
#' @export
create_phased_maplist <- function(input_type = "discrete",
                                  maplist,
                                  dosage_matrix.conv,
                                  dosage_matrix.orig = NULL,
                                  probgeno_df,
                                  chk,
                                  remove_markers = NULL,
                                  original_coding = FALSE,
                                  N_linkages = 2,
                                  lower_bound = 0.05,
                                  ploidy,
                                  ploidy2 = NULL,
                                  marker_assignment.1,
                                  marker_assignment.2,
                                  parent1 = "P1",
                                  parent2 = "P2",
                                  marker_conversion_info = NULL,
                                  log = NULL,
                                  verbose = TRUE) {
  input_type <- match.arg(input_type, choices = c("discrete", "probabilistic"))
  
  if(input_type == "discrete"){
    if(original_coding & is.null(dosage_matrix.orig)) stop("Unconverted dosage matrix should also be specified if original_coding = TRUE")
    if(original_coding & is.null(marker_conversion_info)) stop("marker_conversion_info should also be specified if original_coding = TRUE. \nThis can be obtained by re-running convert_marker_dosages with parameter marker_conversion_info = TRUE, and is then part of the list output of that function") 
  }
  
  mapped_markers <- unlist(lapply(maplist, function(x) as.character(x$marker)))
  
  if(input_type == "discrete"){
    if(!all(mapped_markers %in% rownames(dosage_matrix.conv))) stop("Not all markers on map have corresponding dosages! If duplicated markers were added back to maps, make sure to use an appropriate dosage matrix!")
  } else{
    if(!all(mapped_markers %in% chk$checked_F1$MarkerName)) stop("Not all markers on map were processed by checkF1! If duplicated markers were added back to maps, make sure to use appropriate input chk!")
    # ploidy <- chk1$meta$ploidy
    }
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  if(is.null(ploidy2)) ploidy2 <- ploidy
  
  if(input_type == "probabilistic"){ 
    # ALL_MT <- chk$checked_F1
    ALL_MT <- assign_parental_dosage(chk = chk,probgeno_df = probgeno_df)
    rownames(ALL_MT) <- ALL_MT$MarkerName
    ALL_MT <- ALL_MT[,c("parent1", "parent2")] 
    colnames(ALL_MT) <- c("P1","P2")
  }
  
  if(ploidy == ploidy2){
    if(input_type == "discrete"){
      palindromes <- rownames(dosage_matrix.conv)[which(dosage_matrix.conv[,parent1] != dosage_matrix.conv[,parent2] &
                                                          abs(dosage_matrix.conv[,parent1] - (0.5*ploidy)) == abs(dosage_matrix.conv[,parent2]-(0.5*ploidy2)))]
      
      ## If there are any unconverted palindromes, convert them:
      if(any(dosage_matrix.conv[palindromes,parent1] > dosage_matrix.conv[palindromes,parent2]))
        dosage_matrix.conv[palindromes[dosage_matrix.conv[palindromes,parent1] > dosage_matrix.conv[palindromes,parent2]],] <-
          ploidy - dosage_matrix.conv[palindromes[dosage_matrix.conv[palindromes,parent1] > dosage_matrix.conv[palindromes,parent2]],]
      
    } else{
      palindromes <- rownames(ALL_MT)[which(ALL_MT[,parent1] != ALL_MT[,parent2] &
                                              abs(ALL_MT[,parent1] - (0.5*ploidy)) == abs(ALL_MT[,parent2]-(0.5*ploidy2)))]
      
      ## If there are any unconverted palindromes, convert them:
      if(any(ALL_MT[palindromes,parent1] > ALL_MT[palindromes,parent2]))
        ALL_MT[palindromes[ALL_MT[palindromes,parent1] > ALL_MT[palindromes,parent2]],] <-
          ploidy - ALL_MT[palindromes[ALL_MT[palindromes,parent1] > ALL_MT[palindromes,parent2]],]
      
    }
      }
  
  # Begin by separating the SxN and NxS linkages:
  SxN_assigned <- marker_assignment.1[marker_assignment.1[,parent1]==1 &
                                        marker_assignment.1[,parent2]==0,]
  if(nrow(SxN_assigned) > 0){
    p1_assigned <- marker_assignment.1[-match(rownames(SxN_assigned),rownames(marker_assignment.1)),]
    noSN <- FALSE
  } else{
    p1_assigned <- marker_assignment.1
    noSN <- TRUE
  }
  
  NxS_assigned <- marker_assignment.2[marker_assignment.2[,parent1]==0 &
                                        marker_assignment.2[,parent2]==1,]
  if(nrow(NxS_assigned) > 0){
    p2_assigned <- marker_assignment.2[-match(rownames(NxS_assigned),rownames(marker_assignment.2)),]
    noNS <- FALSE
  } else{
    p2_assigned <- marker_assignment.2
    noNS <- TRUE
  }
  
  #Use only the markers with at least N_linkages significant linkages
  P1unlinked <- rownames(p1_assigned)[apply(p1_assigned[,3+grep("LG",colnames(p1_assigned)[4:ncol(p1_assigned)]),drop = FALSE],1,max)<N_linkages]
  P2unlinked <- rownames(p2_assigned)[apply(p2_assigned[,3+grep("LG",colnames(p2_assigned)[4:ncol(p2_assigned)]),drop = FALSE],1,max)<N_linkages]
  
  if(verbose) {
    removed.m1 <- vector.to.matrix(P1unlinked, n.columns = 4)
    removed.m2 <- vector.to.matrix(P2unlinked, n.columns = 4)
    
    if(nrow(removed.m1) > 0){
      write(paste("\nThe following P1 markers had less than", N_linkages,"significant linkages:\n_______________________________________\n"),log.conn)
      write(knitr::kable(removed.m1,format="markdown"), log.conn)
    }
    
    if(nrow(removed.m2) > 0){
      write(paste("\n\nThe following P2 markers had less than", N_linkages,"significant linkages:\n_______________________________________\n"),log.conn)
      write(knitr::kable(removed.m2,format="markdown"), log.conn)
      write("\n", log.conn)
    }
  }
  
  if(length(P1unlinked) > 0) p1_assigned <- p1_assigned[-match(P1unlinked,rownames(p1_assigned)),]
  if(length(P2unlinked) > 0) p2_assigned <- p2_assigned[-match(P2unlinked,rownames(p2_assigned)),]
  
  # Only select markers for which the number of homologue assignments match the seg type:
  p1cols <- 3+grep("Hom",colnames(p1_assigned)[4:ncol(p1_assigned)])
  p2cols <- 3+grep("Hom",colnames(p2_assigned)[4:ncol(p2_assigned)])
  
  P1rates <- p1_assigned[,p1cols]/rowSums(p1_assigned[,p1cols], na.rm = TRUE)
  P2rates <- p2_assigned[,p2cols]/rowSums(p2_assigned[,p2cols], na.rm = TRUE)
  
  P1rates[P1rates < lower_bound] <- 0
  P2rates[P2rates < lower_bound] <- 0
  
  P1linked <- apply(P1rates,1,function(x) length(which(x!=0)))
  P2linked <- apply(P2rates,1,function(x) length(which(x!=0)))
  
  p1.markers <- rownames(p1_assigned[p1_assigned[,parent1]!=0,])
  p2.markers <- rownames(p2_assigned[p2_assigned[,parent2]!=0,])
  
  ## Assuming markers are converted here; have to treat palindrome markers carefully:
  P1different <- rownames(p1_assigned[rownames(p1_assigned) %in% p1.markers & p1_assigned[,parent1] != P1linked,])
  P2different <- rownames(p2_assigned[setdiff(which(rownames(p2_assigned) %in% p2.markers & p2_assigned[,parent2] != P2linked),
                                              which(rownames(p2_assigned) %in% palindromes & ploidy2 - p2_assigned[,parent2] == P2linked)),])
  
  
  if(verbose) {
    
    removed.m1 <- if(!is.null(P1different)) {
      vector.to.matrix(P1different, n.columns = 4)
    } else matrix(,nrow=0,ncol=1) #catching error
    
    removed.m2 <- if(!is.null(P2different)){
      vector.to.matrix(P2different, n.columns = 4)
    } else matrix(,nrow=0,ncol=1) #catching error
    
    if(nrow(removed.m1) > 0){
      write(paste("\nThe following markers did not have the expected assignment in P1:\n_______________________________________\n"),log.conn)
      write(knitr::kable(removed.m1,format="markdown"), log.conn)
    }
    if(nrow(removed.m2) > 0){
      write(paste("\n\nThe following markers did not have the expected assignment in P2:\n_______________________________________\n"),log.conn)
      write(knitr::kable(removed.m2,format="markdown"), log.conn)
      write("\n", log.conn)
    }
  }
  
  P1rates <- P1rates[!rownames(p1_assigned) %in% P1different,]
  P2rates <- P2rates[!rownames(p2_assigned) %in% P2different,]
  
  #Update p1_assigned and p2_assigned
  p1_assigned <- p1_assigned[!rownames(p1_assigned) %in% P1different,]
  p2_assigned <- p2_assigned[!rownames(p2_assigned) %in% P2different,]
  
  ## Why is this step needed? Seems redundant:
  # rownames(P1rates) <- rownames(p1_assigned)
  # rownames(P2rates) <- rownames(p2_assigned)
  
  # return simplex x nulliplex markers
  if(!noSN) p1_assigned <- rbind(SxN_assigned,p1_assigned)
  if(!noNS) p2_assigned <- rbind(NxS_assigned,p2_assigned)
  
  P1rates <- rbind(SxN_assigned[,p1cols],P1rates)
  P2rates <- rbind(NxS_assigned[,p2cols],P2rates)
  
  # Remove the bi-parental markers that are not assigned in both parents (what about unconverted markers here? Logical test is only looks for a nulliplex parent.)
  bip1 <- rownames(p1_assigned[rowSums(p1_assigned[,c(parent1,parent2)]!=0)==2,])
  bip2 <- rownames(p2_assigned[rowSums(p2_assigned[,c(parent1,parent2)]!=0)==2,])
  
  BiP_different <- c(setdiff(bip1,intersect(bip1,bip2)),setdiff(bip2,intersect(bip1,bip2)))
  
  if (verbose & !is.null(BiP_different)) {
    removed.m <- vector.to.matrix(BiP_different, n.columns = 4)
    
    if(nrow(removed.m) > 0){
      write(paste("\nThe following markers did not have the expected assignment across both parents:\n_______________________________________\n"),log.conn)
      write(knitr::kable(removed.m,format="markdown"), log.conn)
      write("\n", log.conn)
    }
  }
  
  P1rates <- P1rates[!rownames(p1_assigned) %in% setdiff(bip1,intersect(bip1,bip2)),]
  P2rates <- P2rates[!rownames(p2_assigned) %in% setdiff(bip2,intersect(bip1,bip2)),]
  
  #Update p1_assigned and p2_assigned
  p1_assigned <- p1_assigned[!rownames(p1_assigned) %in% setdiff(bip1,intersect(bip1,bip2)),]
  p2_assigned <- p2_assigned[!rownames(p2_assigned) %in% setdiff(bip2,intersect(bip1,bip2)),]
  
  ALL_assigned <- unique(c(rownames(p1_assigned),rownames(p2_assigned)))
  
  # Make up the output
  maplist.out <- lapply(seq(length(maplist)),function(mapn) {
    
    map <- maplist[[mapn]]
    map <- map[map$marker%in%ALL_assigned,]
    
    outmap <- map[,c("marker","position")]
    
    hom_mat <- sapply(1:nrow(outmap), function(r){
      a <- rep(0, ploidy+ploidy2)
      
      temp <- P1rates[match(as.character(outmap$marker[r]),rownames(P1rates)),]
      if(length(which(temp!=0)) > 0) a[(1:ploidy)[which(temp!=0)]] <- 1
      
      temp <- P2rates[match(as.character(outmap$marker[r]),rownames(P2rates)),]
      if(length(which(temp!=0)) > 0) a[((ploidy+1):(ploidy+ploidy2))[which(temp!=0)]] <- 1
      
      return(a)
    })
    
    hom_mat <- t(hom_mat)
    colnames(hom_mat) <- paste0("h",seq(1,ploidy+ploidy2))
    
    # correct palindrome markers:
    if(any(outmap$marker %in% palindromes)){
      hom_mat[outmap$marker %in% palindromes,(ploidy+1):(ploidy+ploidy2)] <-
        (hom_mat[outmap$marker %in% palindromes,(ploidy+1):(ploidy+ploidy2)] + 1) %% 2
    }
    
    # recode using the original coding:
    if(original_coding){
      
      # orig_parents <- dosage_matrix.orig[match(as.character(outmap$marker),rownames(dosage_matrix.orig)),c(parent1,parent2)]
      if(input_type == "discrete"){
        orig_parents <- dosage_matrix.orig[as.character(outmap$marker),c(parent1,parent2)]
      } else{
        orig_parents <- ALL_MT[as.character(outmap$marker),c(parent1,parent2)]
      }
      
      orig_mat <- hom_mat
      ci <- marker_conversion_info[as.character(outmap$marker),]
      
      for(r in 1:nrow(orig_mat)){
        # if(sum(hom_mat[r,1:ploidy]) != orig_parents[r,1]) orig_mat[r,1:ploidy] <- (hom_mat[r,1:ploidy]+1)%%2
        #   
        #   if(sum(hom_mat[r,(ploidy+1):(ploidy+ploidy2)]) != orig_parents[r,2])
        #     orig_mat[r,(ploidy+1):(ploidy+ploidy2)] <- (hom_mat[r,(ploidy+1):(ploidy+ploidy2)]+1)%%2
        
        if(ci[r,1]) orig_mat[r,1:ploidy] <- (hom_mat[r,1:ploidy]+1)%%2
        
        if(ci[r,2]) orig_mat[r,(ploidy+1):(ploidy+ploidy2)] <- (hom_mat[r,(ploidy+1):(ploidy+ploidy2)]+1)%%2
        
      }
      outmap <- cbind(outmap,orig_mat)
    } else{
      outmap <- cbind(outmap,hom_mat)
    }
    
    return(outmap)
  }
  )
  
  names(maplist.out) <- names(maplist)
  
  phased_markers <- unlist(lapply(maplist.out, function(x) as.character(x$marker)))
  
  if(original_coding){
    if(input_type == "discrete"){
      mapped.dosages <- dosage_matrix.orig[mapped_markers,]
    } else{
      mapped.dosages <- ALL_MT[mapped_markers,] #I don't understand why this is used in both instances..
    }
    
  } else{ #converted coding used..
    if(input_type == "discrete"){
      mapped.dosages <- dosage_matrix.conv[mapped_markers,]
    } else{
      mapped.dosages <- as.matrix(ALL_MT[mapped_markers,]) #I don't understand why this is used in both instances..
    }
    
  }
  
  if(verbose){

    mds <- marker_data_summary(dosage_matrix = mapped.dosages,
                                  ploidy = (ploidy+ploidy2)/2,
                                  parent1 = parent1,
                                  parent2 = parent2,
                                  pairing = "random",
                                  shortform = TRUE,
                                  verbose = FALSE)
    
    # mds.aft <- marker_data_summary(dosage_matrix = dosage_matrix.conv[phased_markers,],
    #                                ploidy = (ploidy+ploidy2)/2,
    #                                parent1 = parent1,
    #                                parent2 = parent2,
    #                                pairing = "random", 
    #                                shortform = TRUE,
    #                                verbose = FALSE)

    
    # write(paste("\nMapped marker breakdown before phasing:\n_______________________________________\n"),log.conn)
    # write(knitr::kable(mds.b4$parental_info,format="markdown"), log.conn)
    # write("\n", log.conn)
    
    write(paste("\nPhased marker breakdown:\n_______________________________________\n"),log.conn)
    write(knitr::kable(mds$parental_info,format="markdown"), log.conn)
    write("\n", log.conn)
    
  }
  
  ## Run a final check to make sure that the phased marker dosages equal the original marker dosages:
  phased.dose <- do.call(rbind,lapply(maplist.out, function(x) {
    temp <- cbind(rowSums(x[,paste0("h",1:ploidy)]),
                  rowSums(x[,paste0("h",(ploidy + 1):(ploidy + ploidy2))]))
    rownames(temp) <- x[,"marker"]
    return(temp)
  }))
  
  orig.dose <- mapped.dosages[rownames(phased.dose),c(parent1,parent2)]
  
  conflicting <- which(rowSums(phased.dose == orig.dose) != 2)
  
  if(length(conflicting) > 0){
    warning("Not all phased markers matched original parental dosage. \nPerhaps unconverted marker dosages were supplied as converted dosages by mistake? \nThe following conflicts were detected and removed:")
    warn.df <- cbind(orig.dose[conflicting,,drop = FALSE],phased.dose[conflicting,,drop = FALSE])
    colnames(warn.df) <- c("P1_original","P2_original","P1_phased","P2_phased")
    write(knitr::kable(warn.df,format="markdown"), log.conn)
    
    ## Simply remove these markers from the output:
    rem.markers <- rownames(phased.dose)[conflicting]
    
    maplist.out <- lapply(maplist.out, function(x) x[!x$marker %in% rem.markers,])
  }
  
  if(!is.null(log)) close(log.conn)
  
  return(maplist.out)
}


#' Generate linkage group and homologue structure of SxN markers
#' @description Function which organises the output of \code{cluster_SN_markers} into a data frame of numbered linkage groups and homologues.
#' Only use this function if it is clear from the graphical output of \code{cluster_SN_markers} that there are LOD scores present which define both chromosomes (lower LOD)
#' and homologues (higher LOD).
#' @param cluster_list A list of cluster_stacks, the output of \code{cluster_SN_markers}.
#' @param LOD_chm Integer. The LOD threshold specifying at which LOD score the markers divide into chromosomal groups
#' @param LOD_hom Integer. The LOD threshold specifying at which LOD score the markers divide into homologue groups
#' @param LG_number Integer. Expected number of chromosomes (linkage groups). Note that if this number of clusters are not
#' present at LOD_chm, the function will abort.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return A data.frame with markers classified by homologue and linkage group.
#' @examples
#' data("P1_homologues")
#' ChHomDf<-define_LG_structure(cluster_list=P1_homologues,LOD_chm=3.5,LOD_hom=5,LG_number=5)
#' @export
define_LG_structure <- function(cluster_list,
                                LOD_chm,
                                LOD_hom,
                                LG_number,
                                log = NULL){
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  
  if(LG_number == 1){
    warning(paste("Since only 1 LG expected, all markers at LOD",LOD_hom,"will by default be included in chromosomal cluster!"))
    LGhom.df <- data.frame("SxN_Marker" = cluster_list[[as.character(LOD_hom)]]$marker,
                           "LG" = 1,
                           "homologue" = cluster_list[[as.character(LOD_hom)]]$cluster)
  } else if(length(levels(cluster_list[[as.character(LOD_chm)]]$cluster)) != LG_number){
    stop(paste("define_LG_structure: Incorrect number of linkage groups identified at LOD",LOD_chm))
  } else{
    
    LGhom.df <- merge(cluster_list[[as.character(LOD_chm)]],
                      cluster_list[[as.character(LOD_hom)]],
                      by = "marker")
    colnames(LGhom.df) <- c("SxN_Marker","LG","homologue")
  }
  
  LGhom.df$LG <-  as.numeric(LGhom.df$LG)
  LGhom.df$homologue <-  as.numeric(LGhom.df$homologue)
  
  LGhom.df <- LGhom.df[order(LGhom.df$LG,LGhom.df$homologue),]
  
  unlinked_markers <- setdiff(cluster_list[[as.character(LOD_chm)]]$marker, LGhom.df$SxN_Marker)
  
  LGhom.df$homologue <- do.call("c", lapply(unique(LGhom.df$LG), function(ch) {
    delta <- min(LGhom.df[LGhom.df$LG == ch,]$homologue) - 1 #This assumes related homologues are consecutive
    LGhom.df[LGhom.df$LG == ch,]$homologue <- LGhom.df[LGhom.df$LG == ch,]$homologue - delta
  }))
  
  LGhom.df <- LGhom.df[order(LGhom.df$LG, LGhom.df$homologue),]
  LGhom.df$homologue <-  as.factor(LGhom.df$homologue)
  LGhom.df$LG <-  as.factor(LGhom.df$LG)
  
  # write unlinked markers to standard out
  if (length(unlinked_markers) > 0) {
    write(paste0("\n####SxN Marker(s) lost in clustering step at LOD ",LOD_hom,":\n"), log.conn)
    unl.m <- vector.to.matrix(unlinked_markers, 4)
    write(knitr::kable(unl.m), log.conn)
  }
  if (!is.null(log))
    close(log.conn)
  
  return(LGhom.df)
}


#' Linkage analysis between all markertypes within LG.
#' @description \code{finish_linkage_analysis} is a wrapper for \code{\link{linkage}}, or in the case of probabilistic genotypes, \code{\link{linkage.gp}}.
#' The function performs linkage calculations between all markertypes within a linkage group.
#' @param input_type Can be either one of 'discrete' or 'probabilistic'. For the former (default), \code{dosage_matrix} must be supplied,
#' while for the latter \code{probgeno_df} and \code{chk} must be supplied. 
#' @param marker_assignment A marker assignment matrix with markernames as rownames and at least containing the column \code{"Assigned_LG"}.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
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
#' @param chk Output list as returned by function \code{\link{checkF1}}. This argument is only needed if probabilistic genotypes are used.

#' @param marker_combinations A matrix with four columns specifying marker combinations to calculate linkage.
#' If NULL all combinations are used for which there are rf functions.
#' Dosages of markers should be in the same order as specified in the names of rf functions.
#' E.g. if using 1.0_2.0 and 1.0_3.0 types use: \code{matrix(c(1,0,2,0,1,0,3,0), byrow = TRUE, ncol = 4)}
#' @param target_parent Character string specifying target parent.
#' @param other_parent Character string specifying other parent.
#' @param convert_palindrome_markers Logical. Should markers that behave the same for both parents be converted to a workable format for that parent? E.g.: should 3.1 markers be converted to 1.3?
#' @param ploidy Ploidy level of parent 1. If parent 2 has the same ploidy level, then also the ploidy level of parent 2.
#' @param ploidy2 Integer, by default \code{NULL}. If parental ploidies differ, use this to specify the ploidy of parent 2. Note that in cross-ploidy situations, ploidy2 must be smaller than ploidy.
#' @param pairing Type of pairing at meiosis, with options \code{"random"} or \code{"preferential"}.
#' @param prefPars The estimates for preferential pairing parameters for parent 1 and 2, in range 0 <= p < 2/3. By default this is c(0,0) (so, no preferential pairing).
#' See the function \code{\link{test_prefpairing}} and the vignette for more details.
#' @param LG_number Number of linkage groups (chromosomes).
#' @param verbose Should messages be sent to stdout or log?
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @param \dots (Other) arguments passed to \code{\link{linkage}}
#' @return Returns a matrix with marker assignments. Number of linkages of 1.0 markers are artificial.
#' @examples
#' \dontrun{
#' data("screened_data3", "marker_assignments_P1")
#' linkages_list_P1<-finish_linkage_analysis(marker_assignment=marker_assignments_P1,
#'                                           dosage_matrix=screened_data3,
#'                                           target_parent="P1",
#'                                           other_parent="P2",
#'                                           convert_palindrome_markers=FALSE,
#'                                           ploidy=4,
#'                                           pairing="random",
#'                                           LG_number=5)
#'                                           }
#' @export
finish_linkage_analysis <- function(input_type = "discrete",
                                    marker_assignment,
                                    dosage_matrix,
                                    probgeno_df,
                                    chk,
                                    marker_combinations = NULL,
                                    target_parent = "P1",
                                    other_parent = "P2",
                                    convert_palindrome_markers = TRUE,
                                    ploidy,
                                    ploidy2 = NULL,
                                    pairing = c("random", "preferential"),
                                    prefPars = c(0,0),
                                    LG_number,
                                    verbose = TRUE,
                                    log = NULL,
                                    ...) {
  input_type <- match.arg(input_type, choices = c("discrete","probabilistic"))
  
  if(input_type == "discrete"){
    dosage_matrix <- test_dosage_matrix(dosage_matrix)
    
    if(!target_parent %in% colnames(dosage_matrix) | !other_parent %in% colnames(dosage_matrix))
      stop("Incorrect column name identifiers supplied for parents (target_parents and/or other_parent). Please check!")
    
  } else{
    probgeno_df <- test_probgeno_df(probgeno_df)
    pardose <- assign_parental_dosage(chk = chk,probgeno_df = probgeno_df)
  }
  
  pairing <- match.arg(pairing)
  

  #initialise:
  ploidy.F1 <- ploidy
  sn.ignore <- "_0.1_"
  dip.ignore <- "h3k0"
  
  ## Currently only relevant for triploids:
  if(!is.null(ploidy2)){
    
    if(ploidy2 == ploidy) stop("ploidy2 only needs to be specified if it differs from ploidy, otherwise leave as default (NULL).")
    if(ploidy2 > ploidy) stop("Currently in cross-ploidy mapping, parent1 has to have the higher ploidy level.")
    if(!target_parent %in% c("P1","P2")) stop("To use cross-ploidy functionality, please rename parent1 as P1 and parent2 as P2 (e.g. in colnames() of your dosage_matrix).")
    
    ploidy.F1 <- (ploidy + ploidy2)/2
    
    if(ploidy2 == 2 & target_parent != "P1") {
      sn.ignore <- "_1.0_"
      dip.ignore <- "_2.0_"
    }
  }
  
  if (is.null(marker_combinations)) {
    avail_funs <- ls(getNamespace("polymapR"))
    rfuns <- avail_funs[grep(paste0(substr(pairing, 1, 1),
                                    ploidy.F1, "_"), avail_funs)]
    seg.hits <- grep("seg",rfuns)
    if(length(seg.hits) > 0) rfuns <- rfuns[-seg.hits]
    
    ignore.funs <- grep(sn.ignore,rfuns)
    if(length(ignore.funs) > 0) rfuns <- rfuns[-ignore.funs] #takes out SN funs for other parent
    
    ## Added this to prevent unnecessary calculations for triploids:
    ignore.funs2 <- grep(dip.ignore,rfuns)
    if(length(ignore.funs2) > 0) rfuns <- rfuns[-ignore.funs2] #takes out DN funs (not needed for diploid parent)
    
    marker_combinations <- do.call(rbind, strsplit(rfuns,
                                                   "[_.]"))
    marker_combinations <- marker_combinations[, -1]
    class(marker_combinations) <- "integer"
  }
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  if(verbose){
    write(paste0("There are ", nrow(marker_combinations), " marker combinations\n"),
          log.conn)
  }
  
  r_LOD_lists <- c()
  lgs <- unique(marker_assignment[, "Assigned_LG"])
  lgs <- lgs[order(lgs)]
  
  if (is.null(log) & length(marker_combinations) > 0)
    pb <- txtProgressBar(0, nrow(marker_combinations), style = 3)
  for (i in seq(nrow(marker_combinations))) {
    mtype1 <- marker_combinations[i, 1:2]
    if (marker_combinations[i, 1] == marker_combinations[i,3] & 
        marker_combinations[i, 2] == marker_combinations[i,4]) {
      mtype2 <- NULL
      mtype2_name <- marker_combinations[i, 3:4]
    } else {
      mtype2 <- marker_combinations[i, 3:4]
      mtype2_name <- mtype2
    }
    
    if(verbose){
      write(paste0(
        "\nCalculating r and LOD between ",
        paste(marker_combinations[i,1:2], collapse = "."),
        " and ",
        paste(marker_combinations[i,3:4], collapse = "."),
        " markers.\n"), log.conn)
    }
    
    r_LOD_list <- list()
    for (lg in lgs) {
      if(is.null(log) & verbose){
        write(paste0("Running LG", lg, "..."), log.conn)
      }
      markers <- rownames(marker_assignment)[marker_assignment[,"Assigned_LG"] == lg]
      
      if(input_type == "discrete"){
        r_LOD_list[[paste0("LG", lg)]] <- linkage(dosage_matrix[markers,],
                                                  markertype1 = mtype1,
                                                  markertype2 = mtype2,
                                                  target_parent = target_parent,
                                                  other_parent = other_parent,
                                                  convert_palindrome_markers = convert_palindrome_markers,
                                                  ploidy = ploidy,
                                                  ploidy2 = ploidy2,
                                                  pairing = pairing,
                                                  prefPars = prefPars,
                                                  verbose = F,
                                                  ...)
      } else{
        probgeno.lg <- probgeno_df[probgeno_df$MarkerName %in% markers,]
        chk.lg <- list(checked_F1 = chk$checked_F1[as.character(chk$checked_F1$MarkerName) %in% markers,],
                        meta = chk$meta)
        pardose.lg <- pardose[pardose$MarkerName %in% markers,]
        
        if(nrow(chk.lg$checked_F1) > 1){
        
        r_LOD_list[[paste0("LG", lg)]] <- linkage.gp(probgeno_df = probgeno.lg,
                                                     chk = chk.lg,
                                                     pardose = pardose.lg,
                                                     markertype1 = mtype1,
                                                     markertype2 = mtype2,
                                                     target_parent = target_parent,
                                                     # ploidy = ploidy,
                                                     # ploidy2 = ploidy2,
                                                     prefPars = prefPars,
                                                     verbose = F,
                                                     ...)
        } else{
          r_LOD_list[[paste0("LG", lg)]] <- NULL
        }
      }
    } #for(lg in lgs)...
    
    r_LOD_list_name <- paste0("r_LOD_list_", paste(mtype1,
                                                   collapse = "."), "_", paste(mtype2_name, collapse = "."))
    assign(r_LOD_list_name, get("r_LOD_list"))
    r_LOD_lists <- c(r_LOD_lists, r_LOD_list_name)
    if (!is.null(log))
      setTxtProgressBar(pb, i)
  }
  
  for (i in lgs) { 
    combined_mat <- data.frame(marker_a = character(), marker_b = character(),
                               r = numeric(), LOD = numeric(), phase = character())
    for (list in r_LOD_lists) {
      assign("l", get(list))
      iname <- paste0("LG", i)
      lgnames <- names(l)
      if (iname %in% lgnames) #no longer needed (deprecated)
        combined_mat <- rbind(combined_mat, l[[iname]])
    }
    assign(paste0("LG", i), get("combined_mat"))
  }
  if (!is.null(log))
    close(log.conn)
  return(mget(paste0("LG", lgs))) 
} #finish_linkage_analysis

#' Visualize and get all markertype combinations for which there are functions in polymapR
#'
#' @param ploidy Ploidy level
#' @param pairing Type of pairing. Either "random" or "preferential".
#' @param nonavailable_combinations Logical. Should nonavailable combinations be plotted with grey lines?
#' @return
#' A matrix with two columns. Each row represents a function with the first and second markertype.
#' @export
#' @examples
#' get_markertype_combinations(ploidy = 4, pairing = "random")
get_markertype_combinations <- function(ploidy,
                                        pairing,
                                        nonavailable_combinations = TRUE){
  fnames <- ls(getNamespace("polymapR"))
  prefix <- substr(pairing, 1,1)
  prefix <- paste0(prefix, ploidy, "_")
  fnames <- fnames[grepl(prefix, fnames)]
  write(paste0("There are ", length(fnames), " markercombinations"), stdout())
  fnames.df <- do.call(rbind, strsplit(fnames, "_"))
  fnames.df <- fnames.df[,c(2,3)]
  fnames.df <- gsub(".", "x", fnames.df, fixed = TRUE)
  
  mtypes <- unique(as.vector(fnames.df))
  allcombs <- t(combn(mtypes,2))
  
  fnames.df <- t(apply(fnames.df, 1, function(x) x[order(x)]))
  allcombs <- t(apply(allcombs, 1, function(x) x[order(x)]))
  
  if(nonavailable_combinations){
    g_all <- igraph::graph_from_data_frame(allcombs, directed = F)
    g_fn <- igraph::graph_from_data_frame(fnames.df, directed = F)
    intsct <- attributes(igraph::E(g_all))$vnames %in% attributes(igraph::E(g_fn))$vnames
    igraph::E(g_all)$color <- "grey"
    igraph::E(g_all)$color[intsct] <- "black"
    l <- igraph::layout_in_circle(g_all)
    igraph::plot.igraph(g_all, layout = l, vertex.size = 20,
                        vertex.label.color = "black",
                        vertex.label.family = "sans")
  } else {
    g <- igraph::graph_from_data_frame(fnames.df, directed = F)
    l <- igraph::layout_in_circle(g)
    igraph::plot.igraph(g, layout = l, vertex.size = 20, vertex.color = "black",
                        vertex.label.color = "black",
                        vertex.label.family = "sans")
  }
  
  return(fnames.df)
}


#' gp_overview
#' @description Function to generate an overview of genotype probabilities across a population
#' @param probgeno_df A data frame as read from the scores file produced by function
#' \code{saveMarkerModels} of R package \code{fitPoly}, or equivalently, a data frame containing the following columns:
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
#' Most probable dosage for a particular individual and marker combination, if \code{maxP} exceeds a user-defined threshold (e.g. 0.9), otherwise \code{NA}
#' }
#' }
#' @param cutoff a filtering threshold, by default 0.7, to identify individuals with more than \code{alpha} 
#' non-missing (maximum) genotype probabilities falling below this cut-off. In other words, by using this
#' default settings (\code{cutoff} = 0.7 and \code{alpha} = 0.1), you require that 90% of markers for an individual were assigned a probability of more than 0.7
#' in one of the possible genotype dosage classes. This can help identify problematic individuals with many examples of 
#' diffuse genotype calls. Lowering the threshold allows more diffuse calls to be accepted.
#' @param alpha Option to specify the quantile of an individuals' scores that will be used to test against \code{cutoff}, by default 0.1.
#' @return a list with the following elements:
#' \item{probgeno_df}{
#' Input data, filtered based on chosen \code{cutoff}
#' }
#' \item{population_overview}{
#' data.frame containing summary statistics of each individual's genotyping scores
#' }
#' @examples
#' \dontrun{
#' data("gp_df")
#' gp_overview(gp_df)
#' }
#' @export
gp_overview <- function(probgeno_df,
                        cutoff = 0.7,
                        alpha = 0.1){
  probgeno_df <- test_probgeno_df(probgeno_df)
  
  individual_temp <- as.data.frame(
    do.call(rbind,lapply(levels(probgeno_df$SampleName) ,function(i){
      temp <- probgeno_df[probgeno_df$SampleName %in% i,]
      na.count <- sum(is.na(temp$maxP))
      mean <- mean(temp$maxP,na.rm = TRUE)
      sd <- sd(temp$maxP, na.rm = TRUE)
      min <- min(temp$maxP, na.rm = TRUE)
      less_than_cutoff <- sum(temp$maxP <= cutoff,na.rm = TRUE)
      Qtle <- as.numeric(quantile(temp$maxP, na.rm = TRUE, probs = alpha))
      c(mean,sd,min,round(na.count/nrow(temp),2),Qtle,less_than_cutoff/nrow(temp))
    }))
  )
  colnames(individual_temp) <- c('mean','SD','min','NA',paste(alpha,'quantile'),paste('less than',cutoff))
  rownames(individual_temp) <- as.character(unique(probgeno_df$SampleName))
  
  ##boxplot to show each standard's number
  par(mfrow = c(2,3))
  sapply(colnames(individual_temp),function(c){
    boxplot(individual_temp[[c]],main = c, xlab = c,
            cex.axis = 1.2, cex.main = 2, cex.lab = 1.5)
  })
  
  ##set cutoff
  if(any(individual_temp[,5] <= cutoff)){
    remove_ind <- rownames(individual_temp)[individual_temp[,5] <= cutoff]
    writeLines(paste0(paste0(remove_ind,collapse = ','),' have been removed.'))
    probgeno_df <- droplevels(probgeno_df[!(probgeno_df$SampleName %in% remove_ind),])
  }
  
  par(mfrow = c(1,1)) #return to original par
  
  return(list('probgeno_df' = probgeno_df,
              'population_overview' = individual_temp))
} #gp_overview

#' Assign markers to linkage groups and homologues.
#' @description This is a wrapper combining \code{\link{linkage}} (or \code{\link{linkage.gp}}) and \code{\link{assign_linkage_group}}. 
#' It is used to assign all marker types to linkage groups by using linkage information with 1.0 markers. It allows for input of marker assignments for which this analysis has already been performed.
#' @param input_type Can be either one of 'discrete' or 'probabilistic'. For the former (default), \code{dosage_matrix} must be supplied,
#' while for the latter \code{probgeno_df} and \code{chk} must be supplied. 
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
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
#' @param chk Output list as returned by function \code{\link{checkF1}}. This argument is only needed if probabilistic genotypes are used.
#' @param assigned_list List of \code{data.frames} with marker assignments for which the assignment analysis is already performed.
#' @param assigned_markertypes List of integer vectors of length 2. Specifying the markertypes in the same order as assigned_list.
#' @param SN_functions A vector of function names to be used. If NULL all remaining linkage functions with SN markers are used.
#' @param LG_hom_stack A \code{data.frame} with markernames (\code{"SxN_Marker"}), linkage group (\code{"LG"}) and homologue (\code{"homologue"})
#' @param target_parent A character string specifying the target parent.
#' @param other_parent A character string specifying the other parent.
#' @param convert_palindrome_markers Logical. Should markers that behave the same for both parents be converted to a workable format for that parent? E.g.: should 3.1 markers be converted to 1.3?
#' @param ploidy Ploidy level of parent 1. If parent 2 has the same ploidy level, then also the ploidy level of parent 2.
#' @param ploidy2 Integer, by default \code{NULL}. If parental ploidies differ, use this to specify the ploidy of parent 2. Note that in cross-ploidy situations, ploidy2 must be smaller than ploidy.
#' @param pairing Type of pairing. Either \code{"random"} or \code{"preferential"}.
#' @param LG_number Expected number of chromosomes (linkage groups).
#' @param LOD_threshold LOD threshold at which a linkage is considered significant.
#' @param write_intermediate_files Logical. Write intermediate linkage files to working directory?
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @param \dots Arguments passed to \code{\link{linkage}}
#' @return A \code{data.frame} specifying marker assignments to linkage group and homologue.
#' @examples
#' \dontrun{
#' data("screened_data3", "P1_SxS_Assigned", "P1_DxN_Assigned", "LGHomDf_P1_1")
#' Assigned_markers<-homologue_lg_assignment(dosage_matrix = screened_data3,
#'                                           assigned_list = list(P1_SxS_Assigned, P1_DxN_Assigned),
#'                                           assigned_markertypes = list(c(1,1), c(2,0)),
#'                                           LG_hom_stack = LGHomDf_P1_1,ploidy=4,LG_number = 5,
#'                                           write_intermediate_files=FALSE)
#'                          }
#' @export
homologue_lg_assignment <- function(input_type = "discrete",
                                    dosage_matrix,
                                    probgeno_df,
                                    chk,
                                    assigned_list,
                                    assigned_markertypes,
                                    SN_functions = NULL,
                                    LG_hom_stack,
                                    target_parent = "P1",
                                    other_parent = "P2",
                                    convert_palindrome_markers = TRUE,
                                    ploidy,
                                    ploidy2 = NULL,
                                    pairing = "random",
                                    LG_number,
                                    LOD_threshold = 3,
                                    write_intermediate_files = TRUE,
                                    log = NULL,
                                    ...) {
  input_type <- match.arg(input_type, choices = c("discrete","probabilistic"))
  LG_hom_stack <- test_LG_hom_stack(LG_hom_stack)
  
  if(input_type == "discrete"){
    dosage_matrix <- test_dosage_matrix(dosage_matrix)
    if(!target_parent %in% colnames(dosage_matrix) | !other_parent %in% colnames(dosage_matrix))
      stop("Incorrect column name identifiers supplied for parents (target_parents and/or other_parent). Please check!")
  } else{
    probgeno_df <- test_probgeno_df(probgeno_df)
    pardose <- assign_parental_dosage(chk = chk, 
                                      probgeno_df = probgeno_df)
  }
  
  assigned_markers <- lapply(assigned_list, rownames)
  assigned_markers <- do.call("c", assigned_markers)
  
  if(input_type == "discrete"){
    filt_dosdat <-
      dosage_matrix[!rownames(dosage_matrix) %in% assigned_markers,]
  } else{
    filt_score <- probgeno_df[!probgeno_df$MarkerName %in% assigned_markers,] #this was logically incorrect in Yanlin's version (lacked !)
  }
  
  if(pairing == "random") {
    pairing_abbr <- "r"
  } else if(pairing == "preferential"){
    pairing_abbr <- "p"
  }
  
  sn.grep1 <- "_1.0_"
  sn.grep2 <- "_1.0_1.0"
  ploidy.F1 <- ploidy
  
  if(!is.null(ploidy2)){ # Currently only for triploids
    if(ploidy2 == ploidy) stop("ploidy2 only needs to be specified if it differs from ploidy, otherwise leave as default (NULL).")
    if(ploidy2 > ploidy) stop("Currently in cross-ploidy mapping, parent1 has to have the higher ploidy level.")
    if(!target_parent %in% c("P1","P2")) stop("To use cross-ploidy functionality, please rename parent1 as P1 and parent2 as P2 in colnames() of your dosage_matrix.")
    
    ploidy.F1 <- (ploidy + ploidy2)/2
    
    if(ploidy2 == 2 & target_parent != "P1") {
      sn.grep1 <- "_0.1_"
      sn.grep2 <- "_0.1_0.1"
    }
    
  }
  
  avail_funs <- ls(getNamespace("polymapR"))
  linkage_functions <- avail_funs[grep(paste0(pairing_abbr, ploidy.F1, "_"), avail_funs)]
  
  already_assigned_functions <-
    sapply(assigned_markertypes, function(x) {
      paste0(pairing_abbr, ploidy.F1, sn.grep1, paste(x, collapse = "."))
    })
  
  linkage_functions <-
    linkage_functions[!linkage_functions %in% c(
      paste0(pairing_abbr, ploidy.F1, sn.grep2),
      already_assigned_functions)]
  
  if(is.null(SN_functions)){
    SN_functions <-
      linkage_functions[grep(paste0(pairing_abbr, ploidy.F1, sn.grep1), linkage_functions)]
  }
  
  marker_combinations <-
    do.call(rbind, strsplit(SN_functions, "[_.]"))
  marker_combinations <- marker_combinations[,-1, drop = FALSE]
  class(marker_combinations) <- "integer"
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  #add a progress bar
  if(!is.null(log) & length(marker_combinations) > 0){
    pb <-
      txtProgressBar(
        min = 0,
        max = nrow(marker_combinations),
        style = 3
      )
  }
  
  target.ploidy <- ploidy
  
  if(!is.null(ploidy2) & target_parent == "P2") target.ploidy <- ploidy2
  
  for (i in seq(nrow(marker_combinations))) {
    if(!is.null(log)) sink(log.conn, append = TRUE)
    write(
      paste0(
        "Calculating r and LOD between ",
        paste(marker_combinations[i, 1:2], collapse = "."),
        " and ",
        paste(marker_combinations[i, 3:4], collapse = "."),
        " markers..."
      ),
      stdout()
    )
    
    mtype1 <- marker_combinations[i, 1:2]
    mtype2 <- marker_combinations[i, 3:4]
    
    if(input_type == "discrete"){
      linkage_df <- linkage(
        dosage_matrix = filt_dosdat,
        markertype1 = mtype1,
        markertype2 = mtype2,
        target_parent = target_parent,
        other_parent = other_parent,
        convert_palindrome_markers = convert_palindrome_markers,
        LOD_threshold = 0,
        ploidy = ploidy,
        ploidy2 = ploidy2,
        pairing = pairing,
        verbose = FALSE,
        ...
      )
    } else{
      linkage_df <- linkage.gp(
        probgeno_df = filt_score,
        chk = chk,
        pardose = pardose,
        markertype1 = mtype1,
        markertype2 = mtype2,
        target_parent = target_parent,
        LOD_threshold = 0,
        # ploidy = ploidy,
        # ploidy2 = ploidy2,
        verbose = FALSE,
        ...
      )
    }

    if (write_intermediate_files) {
      mname1 <- paste(marker_combinations[i, 1:2], collapse = "x")
      mname2 <- paste(marker_combinations[i, 3:4], collapse = "x")
      saveRDS(linkage_df,paste0(target_parent, "_", mname1, "_", mname2, ".RDS"))
    }
  
    assignedData <- assign_linkage_group(
      linkage_df = linkage_df,
      LG_hom_stack = LG_hom_stack,
      phase_considered = "coupling",
      LG_number = LG_number,
      LOD_threshold = LOD_threshold,
      ploidy = target.ploidy,
      assign_homologue = TRUE
    )
    
    if (write_intermediate_files) {
      write.table(assignedData,
                  paste0(target_parent, "_", mname2, "_Assigned.txt"),
                  sep = "\t")
    }
    
    assigned_name <-
      paste0("assigned_",
             paste(mtype1, collapse = "."),
             "_",
             paste(mtype2, collapse = "."))
    assign(assigned_name, get("assignedData"))
    assigned_list[[assigned_name]] <- get(assigned_name)
    write("\n________________________________________\n" ,
          stdout())
    if(!is.null(log)) sink()
    if(!is.null(log)) setTxtProgressBar(pb, i)
  }
  
  if(!is.null(log)) sink(log.conn, append = TRUE)
  
  if(input_type == "discrete"){
    marker_assignments <-
      merge_marker_assignments(
        dosage_matrix = dosage_matrix,
        target_parent = target_parent,
        other_parent = other_parent,
        LG_hom_stack = LG_hom_stack,
        SN_linked_markers = assigned_list,
        ploidy = target.ploidy,
        LG_number = LG_number
      )
  } else{
    #Note: this gp version of the original function is not exported 
    marker_assignments <-
      merge_marker_assignments.gp(MarkerType = pardose,
                                  target_parent = target_parent,#.new,
                                  other_parent = other_parent,#.new,
                                  LG_hom_stack = LG_hom_stack,
                                  SN_linked_markers = assigned_list,
                                  ploidy = target.ploidy,
                                  LG_number = LG_number)
  }
  
  if(!is.null(log)) sink()
  
  if (!is.null(log))
    close(log.conn)
  return(marker_assignments)
} #homologue_lg_assignment


#' Calculate recombination frequency, LOD and phase
#' @description \code{linkage} is used to calculate recombination frequency, LOD and phase within one type of marker or between two types of markers.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @param markertype1 A vector of length 2 specifying the first markertype to compare. The first element specifies the dosage in \code{target_parent}, the second in \code{other_parent}.
#' @param markertype2 A vector of length 2 specifying the first markertype to compare. This argument is optional. If not specified, the function will calculate
#' linkage within the markertype as specified by \code{markertype1}.
#' The first element specifies the dosage in \code{target_parent}, the second in \code{other_parent}.
#' @param target_parent Character string specifying the target parent as provided in the columnnames of dosage_matrix
#' @param other_parent Character string specifying the other parent as provided in the columnnames of dosage_matrix
#' @param G2_test Apply a G2 test (LOD of independence) in addition to the LOD of linkage.
#' @param convert_palindrome_markers Logical. Should markers that behave the same for both parents be converted to a workable format for that parent? E.g.: should 3.1 markers be converted to 1.3? If unsure, set to TRUE.
#' @param LOD_threshold Minimum LOD score of linkages to report. Recommended to use for large number (> millions) of marker comparisons in order to reduce memory usage.
#' @param ploidy Integer. The ploidy of parent 1. If parent 2 has the same ploidy level, then also the ploidy level of parent 2.
#' @param ploidy2 Integer, by default \code{NULL}. If parental ploidies differ, use this to specify the ploidy of parent 2. Note that in cross-ploidy situations, ploidy2 must be smaller than ploidy.
#' @param pairing Type of pairing. \code{"random"} or \code{"preferential"}.
#' @param prefPars The estimates for preferential pairing parameters for parent 1 and 2, in range 0 <= p < 2/3. By default this is c(0,0) (so, no preferential pairing).
#' See the function \code{\link{test_prefpairing}} and the vignette for more details.
#' @param combinations_per_iter Optional integer. Number of marker combinations per iteration.
#' @param iter_RAM A (very) conservative estimate of working memory in megabytes used per core. It only takes the size frequency matrices into account. Actual usage is more, especially with large number of linkages that are reported. Reduce memory usage by using a higher LOD_threshold.
#' @param ncores Number of cores to use. Works both for Windows and UNIX (using \code{doParallel}). Use \code{parallel::detectCores()} to find out how many cores you have available.
#' @param verbose Should messages be sent to stdout?
#' @param full_output Logical, by default \code{FALSE}. If \code{TRUE}, the complete output over all phases and showing marker combination counts is returned.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#'
#' @return
#' Returns a data.frame with columns:
#' \itemize{
#' \item{marker_a}{
#'   first marker of comparison. If markertype2 is specified, it has the type of markertype1.
#' }
#' \item{marker_b}{
#'   second marker of comparison. It has the type of markertype2 if specified.
#' }
#' \item{r}{
#'   (estimated) recombinations frequency
#' }
#' \item{LOD}{
#'   (estimated) LOD score
#' }
#' \item{phase}{
#'   phase between markers
#' }
#' }
#' @examples
#' data("screened_data3")
#' SN_SN_P1 <- linkage(dosage_matrix = screened_data3,
#'                    markertype1 = c(1,0),
#'                    target_parent = "P1",
#'                    other_parent = "P2",
#'                    ploidy = 4,
#'                    pairing = "random",
#'                    ncores = 1
#'                    )
#' @export
linkage <- function(dosage_matrix,
                    markertype1 = c(1,0),
                    markertype2 = NULL,
                    target_parent = "P1",
                    other_parent = "P2",
                    G2_test = FALSE,
                    convert_palindrome_markers = TRUE,
                    LOD_threshold = 0,
                    ploidy,
                    ploidy2 = NULL,
                    pairing = c("random", "preferential"),
                    prefPars = c(0, 0),
                    combinations_per_iter = NULL,
                    iter_RAM = 500,
                    ncores = 1,
                    verbose = TRUE,
                    full_output = FALSE,
                    log = NULL) {
  time_start <- Sys.time()
  dosage_matrix <- test_dosage_matrix(dosage_matrix)
  if(!target_parent %in% colnames(dosage_matrix) | !other_parent %in% colnames(dosage_matrix))
    stop("Incorrect column name identifiers supplied for parents (target_parents and/or other_parent). Please check!")
  
  if(identical(markertype1, markertype2)) markertype2 <- NULL
  
  pairing <- match.arg(pairing)
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  win <- Sys.info()["sysname"] == "Windows"
  if (win) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  } else {
    doParallel::registerDoParallel(cores = ncores)
  }
  
  if(is.null(ploidy2)) ploidy2 <- ploidy
  prog_ploidy <- (ploidy + ploidy2)/2
  
  if(!prog_ploidy %in% c(2,3,4,6)) stop(paste("F1 populations of ploidy",prog_ploidy,"are not currently catered for."))
  
  # Below does conversion from 3.1 to 1.3 in tetraploids
  # for hexaploids there are two of these cases: 5.1 and 4.2
  # this function could always be run, but is only needed in special cases.
  # if you want to be safe, make it TRUE
  if (convert_palindrome_markers) {
    if(prog_ploidy %% 2 == 0){ # No "palindrome markers" in populations of odd ploidy
      palindrome_markers <-
        (dosage_matrix[, target_parent] > dosage_matrix[, other_parent]) &
        (abs(dosage_matrix[, target_parent] - (0.5 * ploidy)) == abs(dosage_matrix[, other_parent] -
                                                                       (0.5 * ploidy)))
      if(length(which(palindrome_markers)) > 0)
        dosage_matrix[palindrome_markers,] <-
          ploidy - dosage_matrix[palindrome_markers,]
    }
  }
  
  # Current workaround for odd ploidy populations:
  if(ploidy != ploidy2 & target_parent == "P2") {
    target_parent <- "P1"
    other_parent <- "P2"
  }
  
  # Based on whether we look at within marker combinations or between
  # Different subsets of the dosage_data need to be taken and different types of combinations need to be made.
  if (is.null(markertype2)) {
    # find markertypes and remove parents
    dosage_data <-
      dosage_matrix[dosage_matrix[, target_parent] == markertype1[1] &
                      dosage_matrix[, other_parent] == markertype1[2],-which(
                        colnames(dosage_matrix) == target_parent |
                          colnames(dosage_matrix) == other_parent
                      )]
    n <- nrow(dosage_data)
    
    nomarkers <- FALSE
    if(is.null(n)){
      nomarkers <- TRUE
    } else if (n == 0) {
      nomarkers <- TRUE
    }
    
    if(nomarkers){
      if (verbose){
        message(paste0("No markers in dosage_matrix of type ",
                       markertype1[1], "x", markertype1[2]))
      }
      
      if(win) parallel::stopCluster(cl)
      return(NULL)
    }
    
    combinations <- as.data.frame(t(combn(n, 2)))
    
    markertype2 <- markertype1
    
  } else {
    dosage_data1 <-
      dosage_matrix[dosage_matrix[, target_parent] == markertype1[1] &
                      dosage_matrix[, other_parent] == markertype1[2],-which(
                        colnames(dosage_matrix) == target_parent |
                          colnames(dosage_matrix) == other_parent
                      ), drop = FALSE]
    dosage_data2 <-
      dosage_matrix[dosage_matrix[, target_parent] == markertype2[1] &
                      dosage_matrix[, other_parent] == markertype2[2],-which(
                        colnames(dosage_matrix) == target_parent |
                          colnames(dosage_matrix) == other_parent
                      ), drop = FALSE]
    n1 <- nrow(dosage_data1)
    n2 <- nrow(dosage_data2)
    
    testNULLorZero <- function(x, y = 1){
      if(is.null(x) | is.null(y)){
        return(TRUE)
      } else if (x == 0 | y == 0) {
        return(TRUE)
      }
      return(FALSE)
    }
    
    if(testNULLorZero(n1,n2)){
      if(verbose){
        if(testNULLorZero(n1)) message(paste0("No markers in dosage_matrix of type ",
                                              markertype1[1], "x", markertype1[2]))
        if(testNULLorZero(n2)) message(paste0("No markers in dosage_matrix of type ",
                                              markertype2[1], "x", markertype2[2]))
      }
      
      if(win) parallel::stopCluster(cl)
      return(NULL)
    }
    
    combinations <- expand.grid(1:n1, (n1 + 1):(n1 + n2))
    dosage_data <- rbind(dosage_data1, dosage_data2)
    rm(dosage_data1, dosage_data2)
  }
  
  dosage_data <- t(dosage_data)
  
  # read the table with segregation data
  seg.fname <- paste0("seg_p", prog_ploidy, "_", pairing)
  seg <- get(seg.fname)#,envir = getNamespace("polymapR"))
  segpar <- seg[, c("dosage1", "dosage2")]
  segoff <- seg[, 3:ncol(seg)]
  segpos <- c(0:prog_ploidy)
  
  # get the expected offspring dosages
  offspring_dosage1 <-
    segpos[segoff[segpar$dosage1 == markertype1[1] &
                    segpar$dosage2 == markertype1[2],] > 0]
  offspring_dosage2 <-
    segpos[segoff[segpar$dosage1 == markertype2[1] &
                    segpar$dosage2 == markertype2[2],] > 0]
  
  # get all possible dosage combinations
  dosage_combinations <-
    expand.grid(offspring_dosage1, offspring_dosage2)
  dosage_levels <-
    paste0("n_", dosage_combinations[, 1], dosage_combinations[, 2])
  rownames(dosage_combinations) <- dosage_levels
  
  udc_template <- unique(dosage_combinations[, 1])
  udc_compare <- unique(dosage_combinations[, 2])
  
  ###### combinations per iter calculations here: ###########
  if (is.null(combinations_per_iter)) {
    expected_nr_comparisons <-
      (length(udc_template) + length(udc_compare)) * nrow(dosage_data) * as.numeric(nrow(combinations))
    # 14 bites per logical value
    bites_needed <- 14 * expected_nr_comparisons
    reserve_RAM <- iter_RAM * 1e6
    number_of_iterations <- ceiling(bites_needed / reserve_RAM)
    combinations_per_iter <-
      ceiling(nrow(combinations) / number_of_iterations)
    if (number_of_iterations == 1)
      reserve_RAM <- bites_needed
    if (verbose){
      message(
        paste0(
          "Number of combinations per iteration: ",
          combinations_per_iter,
          "\nReserving approximately ",
          round((
            reserve_RAM + (110 * combinations_per_iter)
          ) / 1e6),
          "MB RAM per iteration"
        )
      )
    }
  }
  split_factor <-
    ceiling((1:nrow(combinations)) / combinations_per_iter)
  combinations_list <- split(combinations, split_factor)
  
  if (verbose){
    message(
      paste(
        "In total",
        nrow(combinations),
        "combinations, which will run in",
        length(unique(split_factor)),
        "iteration(s)...",
        sep = " "
      )
    )
  }
  
  rm(combinations)
  rm(split_factor)
  
  if(pairing == "random") {
    pairing_abbr <- "r"
  } else if(pairing == "preferential"){
    pairing_abbr <- "p"
  }
  
  fname <- paste0(pairing_abbr,
                  prog_ploidy,
                  "_",
                  paste0(markertype1, collapse = "."),
                  "_",
                  paste0(markertype2, collapse="."))
  
  rfun <- get(fname)#,envir = getNamespace("polymapR"))
  
  r_LOD_tot <-
    foreach::foreach(
      i = 1:length(combinations_list),
      .combine = rbind,
      .inorder = F
    ) %dopar% {
      #.inorder=F
      combs <- as.matrix(combinations_list[[i]])
      colnames(combs) <- c("marker_a", "marker_b")
      # get all markers from first vector in the two lists
      template_matrix <- dosage_data[, combs[, 1], drop = F]
      compare_matrix <- dosage_data[, combs[, 2], drop = F]
      
      # initiate count matrix
      count_matrix <-
        matrix(integer(),
               nrow = ncol(template_matrix),
               ncol = length(dosage_levels))
      colnames(count_matrix) <- dosage_levels
      
      # make logical matrices of all possible dosages
      matrix_list_template <- lapply(udc_template, function(x) {
        template_matrix == x
      })
      names(matrix_list_template) <- udc_template
      
      matrix_list_compare <- lapply(udc_compare, function(x) {
        compare_matrix == x
      })
      names(matrix_list_compare) <- udc_compare
      
      # fill in count matrix for each dosage level
      for (level in dosage_levels) {
        # find all progeny that have a certain dosage for the first marker
        markera <-
          matrix_list_template[[as.character(dosage_combinations[level, 1])]]
        
        # find all progeny that have a certain dosage for the second marker
        markerb <-
          matrix_list_compare[[as.character(dosage_combinations[level, 2])]]
        
        # calculate the number of progeny that meet both criteria and fill in matrix
        compare_vec <-
          as.integer(colSums(markera & markerb, na.rm = T))
        count_matrix[, level] <- compare_vec
      }
      
      #####G2 test######
      if(G2_test){
        off1coords <- sapply(offspring_dosage1, function(n) which(n == dosage_combinations[,1]))
        off2coords <- sapply(offspring_dosage2, function(n) which(n == dosage_combinations[,2]))
        
        offspring1sums <- matrix(sapply(1:ncol(off1coords),function(c) rowSums(count_matrix[,off1coords[,c],drop=F])),ncol=ncol(off1coords))
        offspring2sums <- matrix(sapply(1:ncol(off2coords),function(c) rowSums(count_matrix[,off2coords[,c],drop=F])),ncol=ncol(off2coords))
        Totals <- rowSums(count_matrix)
        Sumcounts <- cbind(offspring1sums,offspring2sums)
        
        combos <- expand.grid(1:length(offspring_dosage1),(length(offspring_dosage1) + 1):(length(offspring_dosage1) + length(offspring_dosage2)))
        expected_counts <- matrix(sapply(1:nrow(combos),function(r) Sumcounts[,combos[r,1]]*Sumcounts[,combos[r,2]]/Totals),ncol=nrow(combos))
        G<-2*rowSums(matrix(sapply(1:ncol(count_matrix), function(c) 
          count_matrix[,c,drop=F]*(log(pmax(count_matrix[,c,drop=F],1)) - log(pmax(expected_counts[,c,drop=F],1)))),
          ncol = ncol(count_matrix)))
        df<-(length(offspring_dosage1)-1)*(length(offspring_dosage2)-1) #degrees of freedom
        
        if(df > 1){
          e<-exp(-G/(2*(df-1)))
          G1<-((4-e)*e -3)*(df-1)+G
          LOD_independence<-G1/(2*log(10))
        } else{
          LOD_independence<-G/(2*log(10))
        }
      }
      ################
      
      count_mat <- cbind(combs, count_matrix)
      
      if (pairing == "random") {
        r_list <- rfun(count_mat)
      } else {
        r_list <- rfun(count_mat,p1=prefPars[1],p2=prefPars[2])
      }
      markernames <- count_mat[, c("marker_a", "marker_b"),drop = FALSE]
      rm(count_mat)
      
      # make sure r>0.5 are not chosen
      r_over_0.5 <- r_list$r_mat >= 0.5
      r_list$r_mat[r_over_0.5] <- 1
      
      # based on phasing strategy chose phasing
      if (r_list$phasing_strategy == "MLL") {
        if (!any(is.na(r_list$logL_mat))) {
          # normal case
          r_list$logL_mat[r_over_0.5] <- -1e4
          phase_num <- apply(r_list$logL_mat, 1, which.max)
          
        } else {
          # exceptional case (happens only for few linkage functions)
          # works with missing values
          r_list$logL_mat[r_over_0.5] <- NA
          phase_num <- vector(length = nrow(r_list$logL_mat))
          for (i in seq(nrow(r_list$logL_mat))) {
            # print(i)
            max_logL <- which.max(r_list$logL_mat[i,])
            if (length(max_logL) == 0) {
              if(length(which.min(r_list$r_mat[i,])) == 0){ #Error, no estimation possible
                write(c("\nWARNING: likelihood singularity encountered (possible division by zero). \nAffected marker pair:\n",
                        paste(colnames(dosage_data)[combs[i,]], collapse=" and "),
                        "\nMarker dosage combinations encountered:\n",
                        knitr::kable(count_matrix[i,,drop=FALSE]),
                        "\nSuggest to run checkF1 and remove any markers with q_mult = 0\n"),
                      log.conn)
                phase_num[i] <- NA
              } else{
                phase_num[i] <- which.min(r_list$r_mat[i,])
              }
            } else if (is.na(r_list$r_mat[i, max_logL])) {
              phase_num[i] <- which.min(r_list$r_mat[i,])
            } else {
              phase_num[i] <- max_logL
            }
          }
        }
        
      } else if (r_list$phasing_strategy == "MINR") {
        
        ## Remove any possible undefined estimates:
        r_list$LOD_mat[is.na(r_list$r_mat)] <- 0
        r_list$r_mat[is.na(r_list$r_mat)] <- 0.499
        
        phase_num <- apply(r_list$r_mat, 1, function(x) {
          which.min(x)[1]
        }
        )
      } else {
        stop("Unknown phasing strategy. Check rf functions.")
      }
      
      # fix for NA values of r
      if(any(is.na(phase_num))){
        if(!"unknown" %in% r_list$possible_phases)
          r_list$possible_phases <- c(r_list$possible_phases,"unknown")
        phase_num[is.na(phase_num)] <- length(r_list$possible_phases)
      }
      
      phase <- r_list$possible_phases[phase_num]
      
      # r based on choice of phase
      r_mat_phase <- cbind(r_list$r_mat, phase_num)
      r <- apply(r_mat_phase, 1, function(x) {
        x[x["phase_num"]]
      })
      
      # LOD based on choice of phase
      LOD_mat_phase <- cbind(r_list$LOD_mat, phase_num)
      LOD <- apply(LOD_mat_phase, 1, function(x) {
        x[x["phase_num"]]
      })
      
      if(!full_output) rm(r_list)
      
      r_LOD <- cbind(as.data.frame(markernames), r, LOD, phase)
      negative_r <- which(r_LOD$r < 0 | r_LOD$r == 1)
      r_LOD$r[negative_r] <- 0.499
      r_LOD$LOD[negative_r] <- 0
      levels(r_LOD$phase) <- c(levels(r_LOD$phase), "unknown")
      r_LOD$phase[negative_r] <- "unknown"
      
      r_LOD$LOD[r_LOD$LOD < 0] <- 0
      
      if(G2_test){
        r_LOD <-
          cbind(r_LOD[, c("marker_a", "marker_b", "r", "LOD", "phase")], LOD_independence)
      } else{
        r_LOD <-
          r_LOD[, c("marker_a", "marker_b", "r", "LOD", "phase")]
      }
      
      r_LOD <- r_LOD[r_LOD$LOD >= LOD_threshold,]
      
      if(full_output)
        r_LOD <-
        cbind(r_LOD,
              count_matrix,
              r_list$r_mat[,1:(ncol(r_list$r_mat)-1)],
              r_list$LOD_mat[,1:(ncol(r_list$LOD_mat)-1)],
              r_list$logL_mat[,1:(ncol(r_list$logL_mat)-1)])
      
      return(r_LOD)
    }
  
  # if (win) snow::stopCluster(cl)
  if (win) parallel::stopCluster(cl)
  
  if (verbose){
    message(
      paste0(
        "\nTotal markercombinations: ",
        nrow(r_LOD_tot),
        " at LOD cutoff of ",
        LOD_threshold
      )
    )
  }
  
  # change marker numbers into markernames
  r_LOD_tot$marker_a <- colnames(dosage_data)[r_LOD_tot$marker_a]
  r_LOD_tot$marker_b <- colnames(dosage_data)[r_LOD_tot$marker_b]
  
  time_end <- Sys.time()
  
  if (verbose) {
    timediff <- as.numeric(time_end - time_start, units = "mins")
    
    write(
      paste(
        "\nFor",
        nrow(r_LOD_tot),
        "marker combinations r, LOD and phase were written to output."
      ),
      log.conn
    )
    write(
      paste(
        "\nRun on a",
        Sys.info()["sysname"],
        "machine using",
        ncores,
        ifelse(ncores>1,"cores,","core,"),
        "taking",
        round(timediff, 2),
        "minutes."
      ),
      log.conn
    )
    
  }
  
  if (!is.null(log)) close(log.conn)
  
  # class(r_LOD_tot) <- c("linkage_df", class(r_LOD_tot))
  return(r_LOD_tot)
  
} #linkage

#' Calculate recombination frequency, LOD and phase using genotype probabilities
#' @description \code{linkage.gp} is used to calculate recombination frequency, LOD and phase within one type of marker or between two types of markers.
#' @param probgeno_df A data frame as read from the scores file produced by function \code{saveMarkerModels} of R package \code{fitPoly}, or alternatively, a data frame containing the following columns:
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
#' @param chk Output list as returned by function \code{\link{checkF1}}
#' @param pardose Option to include the most likely (discrete) parental dosage scores, used mainly for internal calls of this function. By default \code{NULL}
#' @param markertype1 A vector of length 2 specifying the first markertype to compare. The first element specifies the dosage in \code{target_parent} (and the second in the other parent).
#' @param markertype2 A vector of length 2 specifying the first markertype to compare. This argument is optional. If not specified, the function will calculate
#' linkage within the markertype as specified by \code{markertype1}.
#' The first element specifies the dosage in \code{target_parent} (and the second in the other parent).
#' @param target_parent Which parent is being targeted (only acceptable options are "P1" or "P2"), ie. which parent is of specific interest? 
#' If this is the maternal parent, please specify as "P1". If the paternal parent, please use "P2". The actual identifiers of the two parents are
#' entered using the arguments \code{parent1_replicates} and \code{parent2_replicates}.
#' @param G2_test Apply a G2 test (LOD of independence) in addition to the LOD of linkage.
#' @param LOD_threshold Minimum LOD score of linkages to report. Recommended to use for large number (> millions) of marker comparisons in order to reduce memory usage.
#' @param prefPars The estimates for preferential pairing parameters for parent 1 and 2, in range 0 <= p < 2/3. By default this is c(0,0) (so, no preferential pairing).
#' See the function \code{\link{test_prefpairing}} and the vignette for more details.
#' @param combinations_per_iter Optional integer. Number of marker combinations per iteration.
#' @param iter_RAM A (very) conservative estimate of working memory in megabytes used per core. It only takes the size frequency matrices into account. Actual usage is more, especially with large number of linkages that are reported. Reduce memory usage by using a higher LOD_threshold.
#' @param ncores Number of cores to use. Works both for Windows and UNIX (using \code{doParallel}). Use \code{parallel::detectCores()} to find out how many cores you have available.
#' @param verbose Should messages be sent to stdout?
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return
#' Returns a data.frame with columns:
#' \itemize{
#' \item{marker_a}{
#'   first marker of comparison. If markertype2 is specified, it has the type of markertype1.
#' }
#' \item{marker_b}{
#'   second marker of comparison. It has the type of markertype2 if specified.
#' }
#' \item{r}{
#'   (estimated) recombinations frequency
#' }
#' \item{LOD}{
#'   (estimated) LOD score
#' }
#' \item{phase}{
#'   phase between markers
#' }
#' }
#' @examples
#' data("gp_df","chk1")
#' SN_SN_P1.gp <- linkage.gp(probgeno_df = gp_df,
#'                           chk = chk1,
#'                           markertype1 = c(1,0),
#'                           target_parent = "P1")
#' @export
linkage.gp <- function(probgeno_df,
                       chk,
                       pardose = NULL,
                       markertype1 = c(1,0),
                       markertype2 = NULL,
                       target_parent = match.arg(c("P1","P2")),
                       G2_test = FALSE,
                       LOD_threshold = 0,
                       prefPars = c(0, 0),
                       combinations_per_iter = NULL,
                       iter_RAM = 500,
                       ncores = 2,
                       verbose = TRUE,
                       log = NULL){
  time_start <- Sys.time()
  
  probgeno_df <- test_probgeno_df(probgeno_df)
  
  if(is.null(pardose)) {
    pardose <- assign_parental_dosage(chk = chk,probgeno_df = probgeno_df)
  } else{
    ## Subset pardose:
    pardose <- pardose[pardose$MarkerName %in% unique(probgeno_df$MarkerName),]
  }
  
  
  #Extract meta data
  F1 <- chk$meta$F1
  parent1_replicates <- chk$meta$parent1
  parent2_replicates <- chk$meta$parent2
  ploidy <- chk$meta$ploidy
  ploidy2 <- ifelse(is.null(chk$meta$ploidy2),chk$meta$ploidy,chk$meta$ploidy2)
  
  if (is.null(log)) {
    log.conn <- stdout()
  }else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  win <- Sys.info()["sysname"] == "Windows"
  
  if (win) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  } else {
    doParallel::registerDoParallel(cores = ncores)
  }

  if(length(unique(probgeno_df$MarkerName)) < 3){
    stop("Please check probgeno_df: not enough markers for linkage analysis.")
  }
  if(!any(unique(probgeno_df$SampleName) %in% parent1_replicates) || !any(unique(probgeno_df$SampleName) %in% parent2_replicates)){
    stop("Please check probgeno_df: not enough information for parents.")
  }
  if(length(which(unique(probgeno_df$SampleName) %in% F1)) < 11){
    stop("Please check the probgeno_df. There is no enough individuals for linkage anlaysis.")
  }
  if (identical(markertype1, markertype2)) {
    markertype2 <- NULL
  }
  # if(is.null(ploidy2)) ploidy2 <- ploidy
  
  prog_ploidy <- (ploidy + ploidy2)/2
  
  if (!prog_ploidy %in% c(3, 4, 6)) {
    stop(paste("F1 populations of ploidy", prog_ploidy, "are not currently catered for."))
  }
  
  ## Here we still use dosage_matrix instead of dosage_data because we do not filter for the marker type
  ## Find all marker's P1 and P2
  score_P <- probgeno_df[probgeno_df$SampleName %in% parent1_replicates | probgeno_df$SampleName %in% parent2_replicates,]
  markers <- as.character(unique(score_P$MarkerName))
  
  #make sure markers appear in the PossibleMarkerType
  markers <- markers[markers %in% pardose$MarkerName]
  
  convert_probability_score <- function(MT_mat = RealMT,
                                        expected_MT = markertype1,
                                        target_parent = "P1",
                                        scores_mat = probgeno_df,
                                        ploidy = ploidy){
    if(target_parent == "P1"){
      expected_MT <- paste0(expected_MT,collapse = ".")
      Tpart <- 1
      NTpart <- 2
    }else{
      expected_MT <- paste0(c(expected_MT[2],expected_MT[1]),collapse = ".")
      Tpart <- 2
      NTpart <- 1
    }
    
    MT_mat$MT <- paste0(MT_mat$parent1, ".",MT_mat$parent2)
    
    if(any(MT_mat$MT != expected_MT)){
      change_markers <- MT_mat[which(MT_mat$MT != expected_MT),]
      sub_scores_mat <- scores_mat[scores_mat$MarkerName %in% as.character(change_markers$MarkerName),]
      #Step 1: reverse all: 0 to 4, 1 to 3, 2 keep as 2
      dosc <- list()
      for(c in 0:ploidy){
        dosc[[as.character(c)]] <- sub_scores_mat[[paste0("P",c)]]
      }
      for(c in 0:ploidy){
        sub_scores_mat[paste0("P",c)] <- dosc[[as.character(ploidy-c)]]
      }
      change_markers_new <- ploidy - change_markers[,c(2,3)]
      change_markers_new$MT <-  paste0(change_markers_new$parent1, ".",change_markers_new$parent2)
      change_markers_new <- cbind("MarkerName" = change_markers$MarkerName,change_markers_new)
      scores_mat[scores_mat$MarkerName %in% as.character(change_markers$MarkerName),] <- sub_scores_mat
      
      #Step 2: make the 2 - 0
      if(any(change_markers[paste0("parent",Tpart)] < change_markers_new[paste0("parent",Tpart)])){
        change_markers_new1 <- change_markers_new[which(change_markers[paste0("parent",Tpart)] < change_markers_new[paste0("parent",Tpart)]),]
        sub_scores_mat <- scores_mat[scores_mat$MarkerName %in% as.character(change_markers_new1$MarkerName),]
        for(t in 0:(ploidy/2-1)){
          sub_scores_mat[paste0("P",t)] <-   sub_scores_mat[paste0("P",(ploidy/2 + t))]+  sub_scores_mat[paste0("P",t)]
        }
        change_markers_new1[paste0("parent",Tpart)]<- change_markers[which(change_markers[paste0("parent",Tpart)] < change_markers_new[paste0("parent",Tpart)]),
                                                                     paste0("parent",Tpart)]
        change_markers_new1[paste0("parent",NTpart)] <- replicate(nrow(change_markers_new1),0)
        
        change_markers_new[which(change_markers[paste0("parent",Tpart)] < change_markers_new[paste0("parent",Tpart)]),]<- change_markers_new1
        change_markers_new$MT <- paste0(change_markers_new$parent1,".",change_markers_new$parent2)
        scores_mat[scores_mat$MarkerName %in% as.character(change_markers_new1$MarkerName),] <- sub_scores_mat
        
      }
      
      #Step 3: substract by ploidy/2
      if(any(change_markers_new$MT != expected_MT)){
        change_markers_new2 <- change_markers_new[which(change_markers_new$MT != expected_MT),]
        sub_scores_mat <- scores_mat[scores_mat$MarkerName %in% as.character(change_markers_new2$MarkerName),]
        
        dosc <- list()
        for(c in 0:ploidy){
          dosc[[as.character(c)]] <- sub_scores_mat[[paste0("P",c)]]
        }
        for(c in 0:ploidy){
          sub_scores_mat[paste0("P",c)] <- dosc[[as.character(ploidy-c)]]
        }
        
        change_markers_new3 <- ploidy - change_markers_new2[,c(2,3)]
        change_markers_new3$MT <-  paste0(change_markers_new3$parent1, ".",change_markers_new3$parent2)
        change_markers_new3 <- cbind("MarkerName" = change_markers_new2$MarkerName,change_markers_new3)
        change_markers_new[which(change_markers_new$MT != expected_MT),] <- change_markers_new3
        scores_mat[scores_mat$MarkerName %in% as.character(change_markers_new3$MarkerName),] <- sub_scores_mat
        
        if(any(change_markers_new2[paste0("parent",Tpart)] < change_markers_new3[paste0("parent",Tpart)])){
          change_markers_new4 <- change_markers_new3[which(change_markers_new2[paste0("parent",Tpart)] < change_markers_new3[paste0("parent",Tpart)]),]
          sub_scores_mat <- scores_mat[scores_mat$MarkerName %in% as.character(change_markers_new4$MarkerName),]
          
          for(t in 0:(ploidy/2-1)){
            sub_scores_mat[paste0("P",t)] <-   sub_scores_mat[paste0("P",(ploidy/2 + t))]+  sub_scores_mat[paste0("P",t)]
          }
          change_markers_new4[paste0("parent",Tpart)]<- change_markers_new2[which(change_markers_new2[paste0("parent",Tpart)] < change_markers_new3[paste0("parent",Tpart)]),
                                                                            paste0("parent",Tpart)]
          change_markers_new4[paste0("parent",NTpart)] <- replicate(nrow(change_markers_new4),0)
          change_markers_new[change_markers_new$MarkerName %in% as.character(change_markers_new4$MarkerName),] <- change_markers_new4
          
          change_markers_new$MT <- paste0(change_markers_new$parent1,".",change_markers_new$parent2)
          scores_mat[scores_mat$MarkerName %in% as.character(change_markers_new4$MarkerName),] <- sub_scores_mat
        }
      }
      
      if(any(change_markers_new$MT != expected_MT)){
        writeLines("There is some errors")
      }
      
    }
    return(scores_mat)
  } #convert_probability_score
  
  pardose$markertype <- paste0(pardose$parent1, ".",pardose$parent2)
  
  #filter the possiblemarkertype
  if(nrow(pardose) != 0){
    #if these is no markertype2, only pick out the markers of markertype1
    if(is.null(markertype2)){
      #check whether there is a need to convert dosages
      if(target_parent == "P1"){
        MT <- markertype1
      }else{
        MT <- c(markertype1[2],markertype1[1])
      }
      #find all possible marker type
      
      if(any(MT %in% c(0,ploidy))){
        zero_loc <- which(MT %in% c(0,ploidy))
        non_zero_loc <- which(!MT %in% c(0,ploidy))
        
        zero <- c(0,4)
        non_zero <- unique(c(MT[non_zero_loc],ploidy-MT[non_zero_loc]))
        
        if(length(non_zero) == 0){
          MT_matrix <- data.frame("Var1" = zero,
                                  "Var2" = zero)
          MT_sum <- paste0(MT_matrix$Var1,".", MT_matrix$Var2)  #all possible MT - character
        }else{
          if(zero_loc < non_zero_loc){
            MT_matrix <- expand.grid(zero,non_zero)
            MT_sum <- paste0(MT_matrix$Var1,".", MT_matrix$Var2) #all possible MT - character
          }else{
            MT_matrix <- expand.grid(non_zero,zero)
            MT_sum <- paste0(MT_matrix$Var1,".", MT_matrix$Var2)
          }
        }
        
        
      }else{
        MT_c <- ploidy - MT
        MT_sum <- unique(c(paste0(MT, collapse = "."), paste0(MT_c, collapse = ".")))
      }
      
      PossibleMarker <- as.character(pardose[pardose$markertype %in% MT_sum,]$MarkerName)
      
      if(length(PossibleMarker) != 0){
        if(length(PossibleMarker) == 1){
          combinations <- data.frame(PossibleMarker,PossibleMarker)
        }else{
          combinations <- t(combn(PossibleMarker,2))
        }
        colnames(combinations) <- c("markera","markerb")
        rownames(combinations) <- seq(1, nrow(combinations),1)
        
        #Do the marker type conversion
        RealMT <- pardose[pardose$MarkerName %in% PossibleMarker,c("MarkerName","parent1","parent2")]
        
        probgeno_df1 <- convert_probability_score(MT_mat = RealMT,
                                                  expected_MT = markertype1,
                                                  target_parent = target_parent,
                                                  scores_mat = probgeno_df,
                                                  ploidy = ploidy)
        probgeno_df <- probgeno_df1
      }else{
        message("There is not enough markers")
        return(NULL)
      }
    }else{ #there is markertype1 and markertype2
      if(target_parent == "P1"){
        MT1 <- markertype1
        MT2 <- markertype2
      }else{
        MT1 <- c(markertype1[2], markertype1[1])
        MT2 <-  c(markertype2[2], markertype2[1])
      }
      #find all possibilities of MT
      if(any(MT1  %in% c(0,ploidy))){
        zero_loc <- which(MT1  %in% c(0,ploidy))
        non_zero_loc <- which(!MT1  %in% c(0,ploidy))
        
        zero <- c(0,ploidy) #BUG? replaced 4 with ploidy..
        non_zero <- unique(c(MT1[non_zero_loc],ploidy-MT1[non_zero_loc]))
        
        if(length(non_zero) == 0){
          MT1_matrix <- data.frame("Var1" = zero,
                                   "Var2" = zero)
          MT1_sum <- paste0(MT1_matrix$Var1,".", MT1_matrix$Var2)  #all possible MT - character
        }else{
          if(zero_loc < non_zero_loc){
            MT1_matrix <- expand.grid(zero,non_zero)
            MT1_sum <- paste0(MT1_matrix$Var1,".", MT1_matrix$Var2) #all possible MT - character
          }else{
            MT1_matrix <- expand.grid(non_zero,zero)
            MT1_sum <- paste0(MT1_matrix$Var1,".", MT1_matrix$Var2) #all possible MT - character
          }
        }
        
        
      }else{
        MT1_c <- ploidy - MT1
        MT1_sum <- unique(c(paste0(MT1, collapse = "."), paste0(MT1_c, collapse = ".")))
      }
      
      if(any(MT2  %in% c(0,ploidy))){
        zero_loc <- which(MT2  %in% c(0,ploidy))
        non_zero_loc <- which(!MT2  %in% c(0,ploidy))
        
        zero <- c(0,ploidy)
        non_zero <- unique(c(MT2[non_zero_loc],ploidy-MT2[non_zero_loc]))
        
        if(length(non_zero) == 0){
          MT2_matrix <- data.frame("Var1" = zero,
                                   "Var2" = zero)
          MT2_sum <- paste0(MT2_matrix$Var1,".", MT2_matrix$Var2)  #all possible MT - character
        }else{
          if(zero_loc < non_zero_loc){
            MT2_matrix <- expand.grid(zero,non_zero)
            MT2_sum <- paste0(MT2_matrix$Var1,".", MT2_matrix$Var2) #all possible MT - character
          }else{
            MT2_matrix <- expand.grid(non_zero,zero)
            MT2_sum <- paste0(MT2_matrix$Var1,".", MT2_matrix$Var2) #all possible MT - character
          }
        }
        
        
      }else{
        MT2_c <- ploidy - MT2
        MT2_sum <- unique(c(paste0(MT2, collapse = "."), paste0(MT2_c, collapse = ".")))
      }
      
      Possiblemarker1 <-  as.character(pardose[pardose$markertype %in% MT1_sum,]$MarkerName)
      Possiblemarker2 <-  as.character(pardose[pardose$markertype %in% MT2_sum,]$MarkerName)
      
      if(length(Possiblemarker1) != 0 & length(Possiblemarker2 != 0)){
        combinations <- expand.grid(Possiblemarker1,Possiblemarker2)
        colnames(combinations) <- c("markera","markerb")
        rownames(combinations) <- seq(1, nrow(combinations),1)
        
        RealMT1 <- pardose[pardose$MarkerName %in% Possiblemarker1,c("MarkerName","parent1","parent2")]
        RealMT2 <- pardose[pardose$MarkerName %in% Possiblemarker2,c("MarkerName","parent1","parent2")]
        
        probgeno_df1 <- convert_probability_score(MT_mat = RealMT1,
                                                  expected_MT = markertype1,
                                                  target_parent = target_parent,
                                                  scores_mat = probgeno_df,
                                                  ploidy = ploidy)
        
        probgeno_df2 <- convert_probability_score(MT_mat = RealMT2,
                                                  expected_MT = markertype2,
                                                  target_parent = target_parent,
                                                  scores_mat = probgeno_df1,
                                                  ploidy = ploidy)
        probgeno_df <- probgeno_df2
      }else{
        combinations <- data.frame()
        message("There is not enough markers")
        
        return(NULL)
      }
      
    }
  } #the return is after filtering, the already made marker combinations
  
  seginfo <- calcSegtypeInfo(ploidy, ploidy2) #obtain the seginfo for further use
  polysomic <- chk$meta$polysomic
  disomic <- chk$meta$disomic 
  mixed <- chk$meta$mixed
  
  seginfo <- selSegtypeInfo(seginfo, polysomic, disomic, mixed)
  seginfoSummary <- segtypeInfoSummary(seginfo)
  
  if(polysomic){
    pairing <- "random"
    pairing_abbr <- "r"
  }
  if(disomic){
    if(polysomic) stop("Cannot run linkage.gp function if both polysomic and disomic are TRUE. Please re-run checkF1 with one of these options FALSE")
    pairing <- "preferential"
    pairing_abbr <- "p"
  }
  
  ##Here we need to call the segregation table according to the ploidy
  if(is.null(markertype2)){
    fname <- paste0(pairing_abbr, prog_ploidy, "_", paste0(markertype1, collapse = "."), "_", paste0(markertype1, collapse = "."))
    
  }else{
    fname <- paste0(pairing_abbr, prog_ploidy, "_", paste0(markertype1, collapse = "."), "_", paste0(markertype2, collapse = "."))
    
  }
  seg.fname <- paste0("seg_p", ploidy, "_", pairing)
  seg <- get(seg.fname)
  segpos <- c(0:prog_ploidy)
  segoff <- seg[, 3:ncol(seg)]
  segpar <- seg[, c("dosage1", "dosage2")]
  
  if (is.null(combinations_per_iter)) {
    expected_nr_comparisons <- nrow(combinations)
    bites_needed <- 14 * expected_nr_comparisons
    reserve_RAM <- iter_RAM * 1e+06
    number_of_iterations <- ceiling(bites_needed/reserve_RAM)
    combinations_per_iter <- ceiling(nrow(combinations)/number_of_iterations)
    if (number_of_iterations == 1)
      reserve_RAM <- bites_needed
    if (verbose) {
      message(paste0("Number of combinations per iteration: ",
                     combinations_per_iter, "\nReserving approximately ",
                     round((reserve_RAM + (110 * combinations_per_iter))/1e+06),
                     "MB RAM per iteration"))
    }
  }
  
  split_factor <- ceiling((1:nrow(combinations))/combinations_per_iter)
  combinations_list <- split(as.data.frame(combinations), split_factor,drop = TRUE)
  if (verbose) {
    message(paste("In total", nrow(combinations), "combinations, which will run in",
                  length(unique(split_factor)), "iteration(s)...",
                  sep = " "))
  }
  
  
  if(is.null(markertype2)){
    offspring_dosage1 <- offspring_dosage2 <- segpos[segoff[segpar$dosage1 == markertype1[1] &
                                                              segpar$dosage2 == markertype1[2], ] > 0]
  }else{
    offspring_dosage1 <- segpos[segoff[segpar$dosage1 == markertype1[1] &
                                         segpar$dosage2 == markertype1[2], ] > 0]
    
    offspring_dosage2 <- segpos[segoff[segpar$dosage1 == markertype2[1] &
                                         segpar$dosage2 == markertype2[2], ] > 0]
  }
  
  dosage_combinations <- expand.grid(offspring_dosage1, offspring_dosage2)
  dosage_levels <- paste0("n_", dosage_combinations[, 1], dosage_combinations[, 2])  #make the name looks like: n_00, n_10, n_01, n_11
  rownames(dosage_combinations) <- dosage_levels
  
  convert_scores <- function(scores,
                             markername,
                             samplename){
    TempArray <- array(0, dim = c(length(samplename),ploidy+1,length(markername)),
                       dimnames = list(c(samplename),
                                       paste0("P",seq(0,ploidy)),
                                       c(markername)))
    for(m in markername){
      O_temp <- scores[scores$MarkerName == m,]
      rownames(O_temp) <- O_temp$SampleName
      TempArray[,,m] <- as.matrix(O_temp[samplename,paste0("P",seq(0,ploidy))])
    }
    return(TempArray)
  }
  
  array <- convert_scores(scores = probgeno_df,
                          markername = unique(c(as.character(combinations[,1]),as.character(combinations[,2]))),
                          samplename = as.character(unique(probgeno_df$SampleName))[as.character(unique(probgeno_df$SampleName)) %in% F1])
  
  rfun <- get(fname)
  
  r_LOD_tot <- foreach::foreach(i = 1:length(combinations_list),
                                .combine = rbind, .inorder = F) %dopar% {
                                  
                                  combs <- as.matrix(combinations_list[[i]])
                                  
                                  if(nrow(combs) > 1){
                                    a <- array[,,combs[,1]]
                                    b <- array[,,combs[,2]]
                                    compare_vec <- matrix(nrow = nrow(combs),ncol = length(dosage_levels))
                                    for (l in 1:length(dosage_levels)){
                                      proba <- a[,paste0("P",as.character(dosage_combinations[l, 1])),]
                                      probb <- b[,paste0("P",as.character(dosage_combinations[l, 2])),]
                                      colnames(proba) <- colnames(probb)  <- NULL
                                      sum <- colSums(proba * probb, na.rm = TRUE)
                                      compare_vec[,l] <- sum
                                    }
                                  }else if(nrow(combs) == 1){
                                    a <- array[,,combs[,1]]
                                    b <- array[,,combs[,2]]
                                    compare_vec <- matrix(nrow = nrow(combs),ncol = length(dosage_levels))
                                    for (l in 1:length(dosage_levels)){
                                      proba <- a[,paste0("P",as.character(dosage_combinations[l, 1]))]
                                      probb <- b[,paste0("P",as.character(dosage_combinations[l, 2]))]
                                      colnames(proba) <- colnames(probb)  <- NULL
                                      sum <- sum(proba * probb, na.rm = TRUE)
                                      compare_vec[,l] <- sum
                                    }
                                  }else{
                                    compare_vec <- NULL
                                  }
                                  
                                  count_matrix <- compare_vec
                                  
                                  if (G2_test) {
                                    off1coords <- sapply(offspring_dosage1, function(n) which(n == dosage_combinations[, 1]))
                                    off2coords <- sapply(offspring_dosage2, function(n) which(n == dosage_combinations[, 2]))
                                    offspring1sums <- matrix(sapply(1:ncol(off1coords),
                                                                    function(c) rowSums(count_matrix[, off1coords[,c], drop = F])), ncol = ncol(off1coords))
                                    offspring2sums <- matrix(sapply(1:ncol(off2coords),
                                                                    function(c) rowSums(count_matrix[, off2coords[,c], drop = F])), ncol = ncol(off2coords))
                                    Totals <- rowSums(count_matrix)
                                    Sumcounts <- cbind(offspring1sums, offspring2sums)
                                    combos <- expand.grid(1:length(offspring_dosage1),
                                                          (length(offspring_dosage1) + 1):(length(offspring_dosage1) + length(offspring_dosage2)))
                                    expected_counts <- matrix(sapply(1:nrow(combos),
                                                                     function(r) Sumcounts[, combos[r, 1]] * Sumcounts[,combos[r, 2]]/Totals), ncol = nrow(combos))
                                    G <- 2 * rowSums(matrix(sapply(1:ncol(count_matrix),
                                                                   function(c) count_matrix[, c, drop = F] * (log(pmax(count_matrix[,c, drop = F], 1)) - log(pmax(expected_counts[,
                                                                                                                                                                                  c, drop = F], 1)))), ncol = ncol(count_matrix)))
                                    df <- (length(offspring_dosage1) - 1) * (length(offspring_dosage2) - 1)
                                    if (df > 1) {
                                      e <- exp(-G/(2 * (df - 1)))
                                      G1 <- ((4 - e) * e - 3) * (df - 1) + G
                                      LOD_independence <- G1/(2 * log(10))
                                    }else {
                                      LOD_independence <- G/(2 * log(10))
                                    }
                                  }
                                  
                                  summary_linkage <- data.frame()
                                  
                                  if(!is.null(compare_vec)){
                                    colnames(count_matrix) <- dosage_levels
                                    count <- count_matrix
                                    if(pairing == "random"){
                                      r_list <- rfun(count)
                                      r_list1 <- rfun(round(count))
                                      if(any(r_list$r_mat < 0) & all(r_list1$r_mat > 0)) {
                                        r_list$r_mat[r_list$r_mat < 0] <- r_list1$r_mat[r_list$r_mat < 0]
                                      }
                                    }
                                    if(pairing == "preferential"){
                                      r_list <- rfun(count, p1 = prefPars[1], p2 = prefPars[2])
                                      r_list1 <- rfun(round(count), p1 = prefPars[1], p2 = prefPars[2])
                                      if(any(r_list$r_mat < 0) & all(r_list1$r_mat > 0)) {
                                        r_list$r_mat[which(r_list$r_mat < 0),] <- r_list1$r_mat[which(r_list$r_mat < 0),]
                                      }
                                    }
                                    
                                    
                                    r_over_0.5 <- r_list$r_mat >= 0.5
                                    r_list$r_mat[r_over_0.5] <- 1
                                    
                                    # based on phasing strategy chose phasing
                                    if (r_list$phasing_strategy == "MLL") {
                                      if (!any(is.na(r_list$logL_mat))) {
                                        # normal case
                                        r_list$logL_mat[r_over_0.5] <- -1e4
                                        phase_num <- apply(r_list$logL_mat, 1, which.max)
                                        
                                      } else {
                                        # exceptional case (happens only for few linkage functions)
                                        # works with missing values
                                        r_list$logL_mat[r_over_0.5] <- NA
                                        phase_num <- vector(length = nrow(r_list$logL_mat))
                                        for (i in seq(nrow(r_list$logL_mat))) {
                                          max_logL <- which.max(r_list$logL_mat[i,])
                                          if (length(max_logL) == 0) {
                                            phase_num[i] <- which.min(r_list$r_mat[i,])
                                          } else if (is.na(r_list$r_mat[i, max_logL])) {
                                            phase_num[i] <- which.min(r_list$r_mat[i,])
                                          } else {
                                            phase_num[i] <- max_logL
                                          }
                                        }
                                      }
                                      
                                    } else if (r_list$phasing_strategy == "MINR") {
                                      phase_num <- apply(r_list$r_mat, 1, function(x) {
                                        which.min(x)[1]
                                      }
                                      )
                                    } else {
                                      stop("Unknown phasing strategy. Check rf functions.")
                                    }
                                    
                                    # fix for NA values of r
                                    if(any(is.na(phase_num))){
                                      if(!"unknown" %in% r_list$possible_phases)
                                        r_list$possible_phases <- c(r_list$possible_phases,"unknown")
                                      phase_num[is.na(phase_num)] <- length(r_list$possible_phases)
                                    }
                                    
                                    phase <- r_list$possible_phases[phase_num]
                                    
                                    # r based on choice of phase
                                    r_mat_phase <- cbind(r_list$r_mat, phase_num)
                                    r <- apply(r_mat_phase, 1, function(x) {
                                      x[x["phase_num"]]
                                    })
                                    
                                    # LOD based on choice of phase
                                    LOD_mat_phase <- cbind(r_list$LOD_mat, phase_num)
                                    LOD <- apply(LOD_mat_phase, 1, function(x) {
                                      x[x["phase_num"]]
                                    })
                                    
                                    rm(r_list)
                                    
                                    r_LOD <- data.frame(combs, r, LOD, phase)
                                    
                                    negative_r <- which(r_LOD$r < 0 | r_LOD$r == 1)
                                    r_LOD$r[negative_r] <- 0.499
                                    r_LOD$LOD[negative_r] <- 0
                                    levels(r_LOD$phase) <- c(levels(r_LOD$phase), "unknown")
                                    r_LOD$phase[negative_r] <- "unknown"
                                    
                                    r_LOD$LOD[r_LOD$LOD < 0] <- 0
                                    
                                    r_LOD <- r_LOD[r_LOD$LOD >= LOD_threshold,]
                                    
                                    colnames(r_LOD)[1:2] <- c("marker_a", "marker_b")
                                    if(G2_test){
                                      r_LOD <-
                                        cbind(r_LOD[, c("marker_a", "marker_b", "r", "LOD", "phase")], LOD_independence)
                                    } else {
                                      r_LOD <-
                                        r_LOD[, c("marker_a", "marker_b", "r", "LOD", "phase")]
                                    }
                                    
                                    summary_linkage <- rbind(summary_linkage, r_LOD)
                                    
                                  }
                                  return(summary_linkage)
                                }
  
  #print the information out
  if (win)
    parallel::stopCluster(cl)
  
  if (verbose) {
    message(paste0("\nTotal marker combinations: ", nrow(r_LOD_tot),
                   " at LOD cutoff of ", LOD_threshold))
    
    time_end <- Sys.time()
  
    # matc <- match.call()
    # write.logheader(matc, log)
    # log.conn <- file(log, "a")
    
    timediff <- as.numeric(time_end - time_start, units = "mins")
    
    write(
      paste(
        "\nFor",
        nrow(r_LOD_tot),
        "marker combinations LOD, r and phase written to output. Calculated on a",
        Sys.info()["sysname"],
        "machine using",
        ncores,
        "cores, taking",
        round(timediff, 2),
        "minutes"
      ),
      log.conn
    )
    
  }
  
  if (!is.null(log)) close(log.conn)
  
  return(r_LOD_tot)
  
} #linkage.gp


#' Perform binning of markers.
#' @description \code{marker_binning} allows for binning of very closely linked markers and choses one representative.
#' @param dosage_matrix A dosage \code{matrix}.
#' @param linkage_df A linkage \code{data.frame}.
#' @param r_thresh Numeric. Threshold at which markers are binned. Is calculated if NA.
#' @param lod_thresh Numeric. Threshold at which markers are binned. Is calculated if NA.
#' @param target_parent A character string specifying the name of the target parent.
#' @param other_parent A character string specifying the name of the other parent.
#' @param max_marker_nr The maximum number of markers per homologue. If specified, LOD threshold is optimized based on this number.
#' @param max_iter Maximum number of iterations to find optimum LOD threshold. Only used if \code{max_marker_nr} is specified.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return A list with the following components:
#' \item{binned_df}{
#'   A linkage data.frame with binned markers removed.
#' }
#' \item{removed}{
#'   A data.frame containing binned markers and their representatives.
#' }
#' \item{left}{
#'   Integer. Number markers left.
#' }
#' @examples
#' data("screened_data3", "all_linkages_list_P1_split")
#' binned_markers<-marker_binning(screened_data3, all_linkages_list_P1_split[["LG2"]][["homologue3"]])
#' @export
marker_binning <-
  function(dosage_matrix,
           linkage_df,
           r_thresh = NA,
           lod_thresh = NA,
           target_parent = "P1",
           other_parent = "P2",
           max_marker_nr = NULL,
           max_iter = 10,
           log = NULL) {
    linkage_df <- test_linkage_df(linkage_df)
    dosage_matrix <- test_dosage_matrix(dosage_matrix)
    if(!target_parent %in% colnames(dosage_matrix) | !other_parent %in% colnames(dosage_matrix))
      stop("Incorrect column name identifiers supplied for parents (target_parents and/or other_parent). Please check!")
    
    if (is.null(log)) {
      log.conn <- stdout()
    } else {
      matc <- match.call()
      write.logheader(matc, log)
      log.conn <- file(log, "a")
    }
    
    rt <- r_thresh
    lt <- lod_thresh
    
    perform_binning <- function(dosage_matrix,
                                linkage_df,
                                r_thresh,
                                lod_thresh,
                                target_parent,
                                other_parent) {
      edges <- linkage_df[linkage_df$LOD > lod_thresh &
                            linkage_df$r < r_thresh,
                          c("marker_a", "marker_b")]
      
      nw <- igraph::graph.data.frame(edges, directed = F)
      gcl <- igraph::groups(igraph::clusters(nw))
      
      count_NA_values <- function(marker_genotypes) {
        sum(is.na(marker_genotypes))
      }
      
      filtered_list <- lapply(gcl, function(markers) {
        sub_dosages <- dosage_matrix[markers,]
        NA_counts <- apply(sub_dosages, 1, count_NA_values)
        SN_marker_subset <- sub_dosages[, target_parent] == 1 &
          sub_dosages[, other_parent] == 0
        SN_markers <- markers[SN_marker_subset]
        SN_markers_NA_counts <- NA_counts[SN_marker_subset]
        
        representing_marker <- ifelse(sum(SN_marker_subset) > 0,
                                      SN_markers[which.min(SN_markers_NA_counts)],
                                      markers[which.min(NA_counts)])
        remove_markers <-
          cbind(markers[markers != representing_marker], representing_marker)
        return(remove_markers)
      })
      
      if (length(filtered_list) == 0) {
        remove_markers <-
          data.frame(removed_marker = character(), representing_marker = character())
      } else {
        remove_markers <- do.call(rbind, filtered_list)
        colnames(remove_markers)[1] <- "removed_marker"
      }
      
      filter_markerA <-
        linkage_df$marker_a %in% remove_markers[,"removed_marker"]
      filter_markerB <-
        linkage_df$marker_b %in% remove_markers[,"removed_marker"]
      
      linkage_filter <- !(filter_markerA | filter_markerB)
      
      out_df <- linkage_df[linkage_filter,]
      nmark <-
        length(unique(as.vector(as.matrix(linkage_df[, c("marker_a", "marker_b")]))))
      nremoved <- nrow(remove_markers)
      if (nrow(remove_markers) == 0) {
        remove_markers <- NULL
      }
      
      return(list(
        binned_df = out_df,
        removed = remove_markers,
        left = nmark - nremoved
      ))
    }
    
    markers <- unique(c(linkage_df$marker_a, linkage_df$marker_b))
    
    thresholds <- calc_binning_thresholds(dosage_matrix)
    
    do_iter <- !is.null(max_marker_nr)
    if (do_iter)
      do_iter <- max_marker_nr < length(markers)
    
    if (do_iter) {
      #write("LOD threshold is optimized", log.conn)
      r_thresh <- thresholds["r_thresh"]
      minLOD <- min(linkage_df$LOD)
      maxLOD <- max(linkage_df$LOD)
      iter <- 0
      mleft <- length(markers)
      
      while (minLOD < maxLOD &
             (iter < max_iter) & (mleft != max_marker_nr)) {
        iter <- iter + 1
        lod_thresh <- (maxLOD + minLOD) / 2
        b <- perform_binning(dosage_matrix,
                             linkage_df,
                             r_thresh,
                             lod_thresh,
                             target_parent,
                             other_parent)
        
        if (b$left < max_marker_nr)
          minLOD <- lod_thresh
        if (b$left > max_marker_nr)
          maxLOD <- lod_thresh
        mleft <- b$left
        # print(paste("markers left", b$left))
        
      }
      write(paste("\nOptimized threshold:", lod_thresh), log.conn)
    }
    
    if (is.na(rt)) {
      r_thresh <-
        thresholds["r_thresh"]#Lets user over-rule this calculation if desired
      write(paste(
        "\nRecombination frequency threshold:", signif(r_thresh, 3)
      ), log.conn)
    }
    if (is.na(lt)) {
      if (!do_iter)
        lod_thresh <- thresholds["lod_thresh"]
      write(paste("\nLOD threshold:", signif(lod_thresh,3)), log.conn)
    }
    
    if (!is.null(log))
      close(log.conn)
    
    return(
      perform_binning(
        dosage_matrix,
        linkage_df,
        r_thresh,
        lod_thresh,
        target_parent,
        other_parent
      )
    )
    
  }




#' Summarize marker data
#' @description Gives a frequency table of different markertypes, relative frequency per markertype of incompatible offspring and the names of incompatible progeny.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @param ploidy Integer. Ploidy of plant species.
#' @param pairing Type of pairing. "random" or "preferential".
#' @param parent1 Name of first parent. Usually maternal parent.
#' @param parent2 Name of second parent. Usually paternal parent.
#' @param progeny_incompat_cutoff The relative number of incompatible dosages per genotype that results in reporting
#'  this genotype as incompatible. Incompatible dosages are greater than maximum number of alleles than can be inherited or
#'  smaller than the minimum number of alleles that can be inherited.
#' @param verbose Logical, by default \code{TRUE} - should intermediate messages be written to stout?
#' @param shortform Logical, by default \code{FALSE}. Returns only a shortened output with parental dosage summary, used internally by some functions.
#' @param log Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.
#' @return Returns a list containing the following components:
#' \item{parental_info}{
#'   frequency table of different markertypes. Names start with parentnames, and behind that the dosage score.
#' }
#' \item{offspring_incompatible}{
#'   Rate of incompatible ("impossible") marker scores (given as percentages of the total number of observed marker scores per marker class)
#' }
#' \item{progeny_incompatible}{
#'   progeny names having incompatible dosage scores higher than threshold at progeny_incompat_cutoff.
#' }
#' @examples
#' data("ALL_dosages")
#' summary_list<-marker_data_summary(dosage_matrix = ALL_dosages, ploidy = 4)
#' @export
marker_data_summary <- function(dosage_matrix,
                                ploidy,
                                pairing = c("random", "preferential"),
                                parent1 = "P1",
                                parent2 = "P2",
                                progeny_incompat_cutoff = 0.1,
                                verbose = TRUE,
                                shortform = FALSE,
                                log = NULL) {
  
  dosage_matrix <- test_dosage_matrix(dosage_matrix)
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  pardos <- dosage_matrix[, c(parent1, parent2)]
  
  if(any(is.na(pardos))){
    NAmark <- rownames(pardos)[is.na(pardos[,parent1]) | is.na(pardos[,parent2])]
    warning("There are parental scores with missing values. These are not considered in the analysis.
            It is recommended to remove those before proceeding to further steps.")
    dosage_matrix <- dosage_matrix[!rownames(dosage_matrix) %in% NAmark, ]
    if(verbose) write(paste(c("\nThe following marker have missing values in their parental scores:",
                              NAmark, "\n"), collapse = "\n\n"), file = log.conn)
  }
  
  
  pairing <- match.arg(pairing)
  
  add_P_to_table <- function(table) {
    #add parent info to table
    colnames(table) <- paste0("P2_", colnames(table))
    rownames(table) <- paste0("P1_", rownames(table))
    return(table)
  }
  
  test_row <- function(x, lu, parpos = c(1, 2)) {
    #analyse offspring incompatibility for a marker
    #with lu as lookup table for maximum and minimum offspring dosages
    progeny <- x[-parpos]
    partype <- lu$pmin == min(x[parpos]) & lu$pmax == max(x[parpos])
    min <- lu[partype, "min"]
    max <- lu[partype, "max"]
    return(!is.na(progeny) & progeny >= min & progeny <= max)
  }
  #######################################
  nm <- nrow(dosage_matrix)
  end_col <- ncol(dosage_matrix)
  
  if(verbose) write("Calculating parental info...", stdout())
  
  # contingency table number of markers
  
  parental_info <-
    table(as.factor(dosage_matrix[,	parent1]), as.factor(dosage_matrix[,	parent2]))
  
  parental_info <- add_P_to_table(parental_info)
  
  if(!shortform){
  
  #Checking offspring compatability
  
  if(verbose) write("Checking compatability between parental and offspring scores...",
                    stdout())
  
  parpos <- which(colnames(dosage_matrix) %in% c(parent1, parent2))
  
  progeny <- dosage_matrix[,-parpos]
  
  nr_offspring <- ncol(progeny)
  
  seg.fname <- paste0("seg_p", ploidy, "_", pairing)
  seg <- get(seg.fname)#,envir=getNamespace("polymapR"))
  segpar <- seg[, c("dosage1", "dosage2")]
  colnames(segpar) <- c("pmax", "pmin")
  segoff <- seg[, 3:ncol(seg)]
  segoff <- segoff > 0
  segpos <- c(0:ploidy)
  
  lu_min_max <- apply(segoff, 1, function(x) {
    a <- segpos[x]
    min <- min(a)
    max <- max(a)
    return(c(min, max))
  })
  
  rownames(lu_min_max) <- c("min", "max")
  lu <- cbind(segpar, t(lu_min_max))
  
  expected_dosage <-
    apply(dosage_matrix, 1, test_row, lu = lu, parpos = parpos)
  
  #NA should be "TRUE", now "FALSE"
  expected_dosage <- t(expected_dosage)
  if(length(which(is.na(progeny))) > 0) expected_dosage[is.na(progeny)] <- TRUE
  
  #two factorial table of parental dosages with percentage of "FALSE" per factor combination
  progeny_incompat <- colSums(!expected_dosage)
  na_progeny <- colSums(is.na(progeny))
  perc_incompat <-
    progeny_incompat / (nrow(expected_dosage) - na_progeny)
  progeny_incompat <-
    colnames(progeny)[perc_incompat > progeny_incompat_cutoff]
  
  nr_incompat <- rowSums(!expected_dosage)
  offspring_incompat <- tapply(
    nr_incompat,
    list(dosage_matrix[, parent1], dosage_matrix[, parent2]),
    FUN = function(x)
      sum(x) / (length(x) * nr_offspring) * 100
  )
  offspring_incompat <- round(offspring_incompat, 2)
  offspring_incompat <- add_P_to_table(offspring_incompat)
  
  summary <-
    list(parental_info, offspring_incompat, progeny_incompat)
  names(summary) <-
    c("parental_info",
      "offspring_incompatible",
      "progeny_incompatible")
  
  for (i in c(1, 2)) {
    if(verbose) {
      write(paste0("\n####", names(summary)[i], "\n"),
            file = log.conn)
      #sink(log.conn)
      write(knitr::kable(summary[[i]]),
            log.conn)
    }
    #suppressWarnings(sink())
  }
  
  if(verbose) write("\n####Incompatible individuals:\n", log.conn)
  if (length(progeny_incompat) == 0 & verbose)
    write("None\n", log.conn)
  
  if(verbose) write(summary$progeny_incompatible, log.conn)
  
  } else{
    summary <- list()
    summary$parental_info <- parental_info
  }
  
  if (!is.null(log))
    close(log.conn)
  
  return(summary)
} #marker_data_summary()


#' Wrapper function for MDSMap to generate linkage maps from list of pairwise linkage estimates
#' @description Create multidimensional scaling maps from a list of linkages
#' @param linkage_list A named \code{list} with r and LOD of markers within linkage groups.
#' @param write_to_file Should output be written to a file? By default \code{FALSE}, if \code{TRUE} then output,
#' including plots from \code{MDSMap} are saved in the same directory as the one used for input files. These
#' plots are currently saved as pdf images. If a different plot format is required (e.g. for publications),
#' then run the \code{MDSMap} function \code{\link[MDSMap]{estimate.map}} (or similar) directly and save the output
#' with a different plotting function as wrapper around the map function call.
#' @param mapdir Directory to which map input files are initially written. Also used for output if \code{write_to_file=TRUE}
#' @param plot_prefix prefix for the filenames of output plots.
#' @param log Character string specifying the log filename to which standard output should be written.
#' If NULL log is send to stdout.
#' @param \dots Arguments passed to \code{\link[MDSMap]{estimate.map}}.
#' @examples
#' \dontrun{
#' data("all_linkages_list_P1")
#' maplist_P1 <- MDSMap_from_list(all_linkages_list_P1[1])
#' }
#' @export
MDSMap_from_list <- function(linkage_list,
                             write_to_file = FALSE,
                             mapdir = "mapping_files_MDSMap",
                             plot_prefix = "",
                             log = NULL,
                             ...) {
  if(class(linkage_list) != "list"){
    stop(paste("linkage_list should be a list, now it's a", class(linkage_list)))
  }
  
  if(is.null(names(linkage_list))){
    stop("linkage_list should be named.
         Apply names to the list (e.g. homologue or LG names) using the function names")
  }
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  if(!file.exists(mapdir)) dir.create(mapdir)
  
  maplist <- list()
  
  for(lg in names(linkage_list)){
    pwd <- prepare_pwd(linkage_list[[lg]])
    
    write(length(unique(c(pwd$marker_a,pwd$marker_b))),
          file = file.path(mapdir,paste0("pwd_",lg,".txt")))
    
    write.table(pwd,file=file.path(mapdir,paste0("pwd_",lg,".txt")),
                row.names = FALSE,col.names = FALSE,
                quote=FALSE,append = TRUE)
    
    if(write_to_file) pdf(file.path(mapdir,paste0("MDSplots_",lg,".pdf")),
                          width=9,height=5)
    
    map <- MDSMap::estimate.map(
      file.path(mapdir,paste0("pwd_",lg,".txt")), ...)
    
    if(write_to_file) dev.off()
    
    maplist[[lg]] <- map
    
  }
  
  if(write_to_file){
    write_nested_list(maplist, directory = mapdir)
  }
  
  logmat <- matrix(nrow = length(linkage_list), ncol = 3)
  rownames(logmat) <- names(linkage_list)
  colnames(logmat) <-
    c("map_length", "stress", "nnfit")
  
  maplist_simpl <- list()
  for(lg in names(maplist)){
    maplist_simpl[[lg]] <-
      maplist[[lg]][["locimap"]][,c("locus", "position", "confplotno", "nnfit")]
    colnames(maplist_simpl[[lg]])[1] <- "marker"
    logmat[lg, "map_length"] <- maplist[[lg]]$length
    logmat[lg, "stress"] <- maplist[[lg]]$smacofsym$stress
    logmat[lg, "nnfit"] <-  maplist[[lg]]$meannnfit
  }
  
  write(knitr::kable(logmat),
        log.conn)
  
  if (!is.null(log))
    close(log.conn)
  return(maplist_simpl)
}


#' Merge homologues
#' @description Based on additional information, homologue fragments, separated during clustered should be merged again.
#' \code{merge_homologues} allows to merge homologues per linkage group based on user input.
#' @param LG_hom_stack A \code{data.frame} with markernames, linkage group (\code{"LG"}) and homologue (\code{"homologue"})
#' @param ploidy The ploidy level of the plant species.
#' @param LG The linkage group where the to be merged homologue fragments are in.
#' @param mergeList A list of vectors of length 2, specifying the numbers of the homologue fragments to be merged. User input is asked if \code{NULL}.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return A modified LG_hom_stack
#' @examples
#' data("LGHomDf_P2_1")
#' merged<-merge_homologues(LGHomDf_P2_1,ploidy=4,LG=2,mergeList=list(c(1,5)))
#' @export
merge_homologues <-
  function(LG_hom_stack,
           ploidy,
           LG,
           mergeList = NULL,
           log = NULL) {
    LG_hom_stack <- test_LG_hom_stack(LG_hom_stack)
    if (is.null(mergeList)) {
      user_input <-
        readline(
          paste0(
            "Provide homologue numbers to be merged for linkage group ",
            LG,
            " \n(dash for combinations, space between combinations e.g: 1-6 3-5):"
          )
        )
      user_input <- as.list(strsplit(user_input, " ")[[1]])
      mergeList <-
        lapply(user_input, function(x)
          as.numeric(strsplit(x, "-")[[1]]))
    }
    edges <- do.call(rbind, mergeList)
    nw <- igraph::graph.data.frame(edges, directed = F)
    gcl <- igraph::groups(igraph::clusters(nw))
    
    # built in error message if two combinations harbor the same homolog?
    for (combination in gcl) {
      # change the homolog number of the second homolog of the combination
      # into the number of the first homolog of the combination
      LG_hom_stack$homologue[LG_hom_stack$LG == LG &
                               (LG_hom_stack$homologue %in% combination)] <-
        combination[1]
    }
    
    # make tapply?:
    for (chm in levels(LG_hom_stack$LG)) {
      # rename factors in ascending order:
      clusters_chm <-
        as.factor(as.numeric(LG_hom_stack$homologue[LG_hom_stack$LG == chm]))
      levels(clusters_chm) <- 1:length(levels(clusters_chm))
      LG_hom_stack$homologue[LG_hom_stack$LG == chm] <-
        as.numeric(clusters_chm)
    }
    LG_hom_stack$homologue <-
      as.factor(as.character(LG_hom_stack$homologue))
    #colnames(LG_hom_stack)[which(colnames(LG_hom_stack)=="homologue")]<-"homologue"
    
    if (is.null(log)) {
      log.conn <- stdout()
    } else {
      matc <- match.call()
      write.logheader(matc, log)
      log.conn <- file(log, "a")
    }
    
    write("####Number of markers per homologue:\n", log.conn)
    write(knitr::kable(
      table(
        LG_hom_stack$homologue,
        LG_hom_stack$LG,
        dnn = list("homologue", "LG")
      ),
      row.names = TRUE
    ),
    log.conn)
    if (!is.null(log))
      close(log.conn)
    
    return(LG_hom_stack)
  }

#' Merge marker assignments
#' @description \code{merge_marker_assignments} Merges 1.0 backbone object with marker assignment objects
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @param target_parent Character string specifying target parent.
#' @param other_parent Character string specifying other parent.
#' @param LG_hom_stack data.frame specifying 1.0 marker assignments to linkage groups and homologues.
#' @param SN_linked_markers a list of marker assignment objects
#' @param ploidy Ploidy level of plant species.
#' @param LG_number Number of linkage groups (chromosomes).
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return Returns a matrix with marker assignments. Number of linkages of 1.0 markers are artificial.
#' @examples
#' data("screened_data3", "LGHomDf_P1_1", "P1_SxS_Assigned", "P1_DxN_Assigned")
#' merged_assignment<-merge_marker_assignments(screened_data3, target_parent="P1",
#'                          other_parent="P2",
#'                          LG_hom_stack=LGHomDf_P1_1,
#'                          SN_linked_markers=list(P1_SxS_Assigned, P1_DxN_Assigned),
#'                          ploidy=4,
#'                          LG_number=5)
#' @export
merge_marker_assignments <- function(dosage_matrix,
                                     target_parent = "P1",
                                     other_parent = "P2",
                                     LG_hom_stack,
                                     SN_linked_markers,
                                     ploidy,
                                     LG_number,
                                     log = NULL) {
  LG_hom_stack <- test_LG_hom_stack(LG_hom_stack)
  dosage_matrix <- test_dosage_matrix(dosage_matrix)
  if(!target_parent %in% colnames(dosage_matrix) | !other_parent %in% colnames(dosage_matrix))
    stop("Incorrect column name identifiers supplied for parents (target_parents and/or other_parent). Please check!")
  
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
  }),nrow = LG_number, dimnames = list(paste0("LG", levels(LG_hom_stack$Assigned_LG)),markers_LG_hom_stack)
  ))
  
  LG_mat <- rbind(SN_LG_mat, comb[, paste0("LG", levels(LG_hom_stack$Assigned_LG)), drop = FALSE])
  LG_mat <-
    cbind(c(as.numeric(as.character(LG_hom_stack[, "Assigned_LG"])), comb[, "Assigned_LG"]), LG_mat)
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
  
  matched_rows <- match(rownames(Assigned_LG_hom),rownames(dosage_matrix))
  
  if(any(is.na(matched_rows))){
    stop("Could not find all assigned markers in dosage_matrix. Please check supplied dosage_matrix is correct.")
  }
  
  parental_dosages <-
    dosage_matrix[matched_rows, c(target_parent, other_parent)]
  
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
}



#' Plotting 1.0 links between homologues
#' @description \code{overviewSNlinks} is written to enable merging of homologue fractions.
#' Fractions of homologues will have more markers in coupling than in repulsion, whereas separate homologues will only have markers in repulsion.
#' @param linkage_df A data.frame as output of \code{\link{linkage}} with arguments markertype1=c(1,0) and markertype2=NULL.
#' @param LG_hom_stack A data.frame with a column "SxN_Marker" specifying markernames,
#' a column "homologue" specifying homologue cluster and "LG" specifying linkage group.
#' @param LG Integer. Linkage group number of interest.
#' @param LOD_threshold Numeric. LOD threshold of linkages which are plotted.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @param ymax Maximum y-limit of the plots.
#' @examples
#' data("SN_SN_P1", "LGHomDf_P1_1")
#' overviewSNlinks(linkage_df=SN_SN_P1,
#'                LG_hom_stack=LGHomDf_P1_1,
#'                LG=5,
#'                LOD_threshold=3)
#' @export
overviewSNlinks <- function(linkage_df,
                            LG_hom_stack,
                            LG,
                            LOD_threshold,
                            ymax = NULL,
                            log = NULL) {
  LG_hom_stack <- test_LG_hom_stack(LG_hom_stack)
  linkage_df <- test_linkage_df(linkage_df)
  ## First, get all unique pairwise comparisons (there are n*(n-1)/2 combinations)
  linkage_df <- linkage_df[linkage_df$LOD > LOD_threshold,]
  mA_filtered <- LG_hom_stack[LG_hom_stack[, "LG"] == LG,]
  homs <- sort(unique(as.character(mA_filtered[, "homologue"])))
  homCombs <- t(combn(homs, 2))
  
  # all plots have the same height width ratio as "Plots" tab:
  numCols <- ceiling(sqrt(nrow(homCombs) + 1))
  
  # save default par
  default.mar <- par("mar")
  
  # new par and layout
  layout(matrix(1:numCols ^ 2, ncol = numCols, byrow = TRUE))
  
  par(mar = c(2.5, 2.5, 2.5, 0.5))
  
  plot.new()
  legend(
    "topleft",
    pch = 19,
    col = c("limegreen", "red3"),
    legend = c("coupling", "repulsion"),
    cex = 1.25,
    bty = "n"
  )
  
  for (i in 1:nrow(homCombs)) {
    plot_SNlinks(
      linkage_df = linkage_df,
      LG_hom_stack = LG_hom_stack,
      LG = LG,
      h1 = homCombs[i, 1],
      h2 = homCombs[i, 2],
      ymax = ymax
    ) #All repulsion. Don't join
  }
  
  # set back layout and graphical parameters
  layout(matrix(1))
  par(mar = default.mar, mfrow=c(1,1))
  
  if (!is.null(log)) {
    matc <- match.call()
    write.logheader(matc, log)
  }
  
} #overviewSNlinks()


#' Calculate frequency of each markertype.
#' @description Plots and returns frequency information for each markertype.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @param parent1 Character string specifying the first (usually maternal) parentname.
#' @param parent2 Character string specifying the second (usually paternal) parentname.
#' @param \dots Arguments passed to \code{\link{barplot}}
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return A named vector containing the frequency of each markertype in the dataset.
#' @examples
#' data("ALL_dosages","screened_data")
#' parental_quantities(dosage_matrix=ALL_dosages)
#' parental_quantities(dosage_matrix=screened_data)
#' @export
parental_quantities <-
  function(dosage_matrix,
           parent1 = "P1",
           parent2 = "P2",
           log = NULL,
           ...) {
    dosage_matrix <- test_dosage_matrix(dosage_matrix)
    parental_info <-
      table(as.factor(dosage_matrix[, parent1]), as.factor(dosage_matrix[, parent2]))
    markertypes <-
      expand.grid(list(rownames(parental_info), colnames(parental_info)))
    parental_quantities <- as.vector(parental_info)
    names(parental_quantities) <-
      markertypes <- paste0(markertypes[, 1], "x", markertypes[, 2])
    pq_plot <- parental_quantities[parental_quantities > 0]
    barplot(
      pq_plot,
      col = c("darkorchid"),
      main = "Marker segregation summary",
      ylab = "Nr. markers",
      cex.names = 0.8,
      ...
    )
    mtext(paste("Total number markers:", sum(parental_quantities)), 3, font = 3)
    
    if (is.null(log)) {
      log.conn <- stdout()
    } else {
      matc <- match.call()
      write.logheader(matc, log)
      log.conn <- file(log, "a")
    }
    freqs <- t(as.data.frame(pq_plot))
    rownames(freqs) <- "frequency"
    write(knitr::kable(freqs), file = log.conn)
    
    if (!is.null(log))
      close(log.conn)
    
    #print(pq_plot)
    return(parental_quantities)
  }

#' Perform a PCA on progeny
#' @description Principal component analysis in order to identify individuals that deviate from the population.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @param highlight A list of character vectors specifying individual names that should be highlighted
#' @param colors Highlight colors. Vector of the same length as \code{highlight}.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @details Missing values are imputed by taking the mean of marker dosages per marker.
#' @examples
#' data("ALL_dosages")
#' PCA_progeny(dosage_matrix=ALL_dosages, highlight=list(c("P1", "P2")), colors="red")
#' @export
PCA_progeny <-
  function(dosage_matrix,
           highlight = NULL,
           colors = NULL,
           log = NULL) {
    
    dosage_matrix <- test_dosage_matrix(dosage_matrix)
    
    if (length(highlight) != length(colors))
      stop("Arguments highlight and colors should have the same length")
    
    class(dosage_matrix) <- "numeric"
    
    na.imputed <- apply(dosage_matrix, 1, function(x) {
      m <- mean(x, na.rm = T)
      x[is.na(x)] <- m
      return(x)
    })
    
    fit <- prcomp(na.imputed)
    cox <- 1
    coy <- 2
    
    # get scores
    scores <- fit$x[, c(cox, coy)]
    
    # calculate variance explained
    ve1 <- fit$sdev[cox] ^ 2 / sum(fit$sdev ^ 2)
    ve2 <- fit$sdev[coy] ^ 2 / sum(fit$sdev ^ 2)
    
    # score plot
    plot(
      0,
      xlim = c(min(scores[, 1]), max(scores[, 1])),
      ylim = c(min(scores[, 2]), max(scores[, 2])),
      xlab = paste('PC', cox, ' (', round(ve1 * 100, 1), '%)', sep = ''),
      ylab = paste('PC', coy, ' (', round(ve2 * 100, 1), '%)', sep = ''),
      bty = 'n',
      type = "n",
      cex.lab = 1.2,
      cex.axis = 1.2
    )
    text(scores[, 1], scores[, 2], rownames(scores), cex = 0.5)
    
    if (!is.null(highlight)) {
      for (i in seq(highlight)) {
        if(length(highlight[[i]]) > 0)
          text(scores[highlight[[i]], 1],
               scores[highlight[[i]], 2],
               highlight[[i]],
               cex = 0.5,
               col = colors[i])
      }
    }
    abline(0, 0, lty = 2)
    abline(v = 0, lty = 2)
    
    if (!is.null(log)) {
      matc <- match.call()
      write.logheader(matc, log)
    }
  }


#' Phase 1.0 markers at the diploid level
#' @description \code{phase_SN_diploid} phases simplex x nulliplex markers for a diploid parent.
#' @param linkage_df A linkage data.frame as output of \code{\link{linkage}} calculating linkage between 1.0 markers.
#' @param cluster_list A list of cluster_stacks, the output of \code{cluster_SN_markers}.
#' @param LOD_chm Integer. The LOD threshold specifying at which LOD score the markers divide into chromosomal groups
#' @param LG_number Expected number of chromosomes (linkage groups)
#' @param independence_LOD Logical. Should the LOD of independence be used for clustering? (by default, \code{FALSE}.)
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout (console).
#' @return  A data.frame with markers classified by homologue and linkage group.
#' @examples
#' data("SN_SN_P2_triploid","P2_homologues_triploid")
#' cluster_list2<-phase_SN_diploid(SN_SN_P2_triploid,P2_homologues_triploid,LOD_chm=5,LG_number = 3)
#' @export
phase_SN_diploid <- function(linkage_df,
                             cluster_list,
                             LOD_chm = 3.5,
                             LG_number,
                             independence_LOD = FALSE,
                             log = NULL) {
  if(length(unique(cluster_list[[as.character(LOD_chm)]]$cluster)) != LG_number) {
    warning(paste0("Unexpected number (",length(unique(cluster_list[[as.character(LOD_chm)]]$cluster)),
                   ") of linkage groups found at LOD ",LOD_chm))
    LG_number <- length(unique(cluster_list[[as.character(LOD_chm)]]$cluster))
  }
  
  chm.clusters <- cluster_list[[as.character(LOD_chm)]]
  
  if(any(is.na(chm.clusters$cluster))){
    warning(paste("The following markers had no cluster assigned and were removed:",
                  paste0(chm.clusters$marker[is.na(chm.clusters$cluster)], collapse = ", ")))
    
    #Is 2 rows enough here? Seems a bit low..
    if(nrow(chm.clusters) < 2) stop("Insufficient marker information to proceed. Suggest to re-check cluster_list!")
  }
  
  linkage_df <- test_linkage_df(linkage_df)
  
  # total_marker <-
  # unique(c(
  # as.character(linkage_df$marker_a),
  # as.character(linkage_df$marker_b)
  # ))
  
  # total_markernr <- length(total_marker)
  
  LODscore <- "LOD"
  if(independence_LOD){
    if(! "LOD_independence" %in% colnames(linkage_df))
      stop("The column LOD_independence should be part of linkage_df when clustering with LOD of independence.
           To obtain the LOD of independence, re-run linkage() with G2_test = TRUE")
    LODscore <- "LOD_independence"
  }
  
  linkage_df.c <- linkage_df[linkage_df[,"phase"]=="coupling",]
  
  # define edges between SN markers
  edges <- linkage_df.c[linkage_df.c[,LODscore] >= LOD_chm, c("marker_a", "marker_b")]
  
  write(paste("Total number of edges:", nrow(edges)), stdout())
  
  # make a network (graph)
  nw <- igraph::graph.data.frame(edges, directed = F)
  gcl <- igraph::groups(igraph::clusters(nw))
  ngroups <- length(gcl)
  groups <- stack(gcl)
  colnames(groups) <- c("marker", "cluster")
  
  # open file for writing if !is.null(log)
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  ## Add any informative repulsion linkages over LOD_chm by symmetry:
  linkage_df.r <- linkage_df[linkage_df[,"phase"]=="repulsion",]
  edges.r <- linkage_df.r[linkage_df.r[,LODscore] >= LOD_chm, c("marker_a", "marker_b")]
  
  repul.linked <- unique(c(edges.r$marker_a, edges.r$marker_b))
  
  if(length(repul.linked) > length(groups$marker)){ #Some repulsion info to use
    unlinked.c <- setdiff(repul.linked,groups$marker)
    
    assigned.clusters <- sapply(unlinked.c, function(merker) {
      replinks.temp <- unique(c(edges.r[edges.r$marker_a == merker,"marker_b"],edges.r[edges.r$marker_b == merker,"marker_a"]))
      chm.temp <- Mode(chm.clusters[chm.clusters$marker %in% replinks.temp,"cluster"]) #polymapR:::Mode
      chm.clust <- groups[groups$marker %in% chm.clusters[chm.clusters$cluster == chm.temp,"marker"],]
      rep.clust <- Mode(chm.clust[chm.clust$marker %in% replinks.temp,"cluster"]) #polymapR:::Mode
      setdiff(unique(chm.clust$cluster),rep.clust)
    }
    )
    
    groups <- rbind(groups,data.frame("marker" = unlinked.c,
                                      "cluster" = assigned.clusters))
  } else{
    write(paste("Complete phase assignment possible using only coupling information at LOD", LOD_chm),log.conn)
  }
  
  LGhom.df <- merge(chm.clusters, groups, by = "marker")
  
  colnames(LGhom.df) <- c("SxN_Marker","LG","homologue")
  LGhom.df$homologue <-  as.numeric(LGhom.df$homologue)
  LGhom.df$LG <-  as.numeric(LGhom.df$LG)
  LGhom.df <- LGhom.df[order(LGhom.df$LG,LGhom.df$homologue),]
  
  unlinked_markers <- setdiff(cluster_list[[as.character(LOD_chm)]]$marker, LGhom.df$SxN_Marker)
  
  LGhom.df[,"homologue"] <- do.call("c", lapply(unique(LGhom.df$LG), function(ch) {
    oldhom.vect <- LGhom.df[LGhom.df$LG == ch,]$homologue
    old.homs <- unique(oldhom.vect)
    new.homs <- seq(length(old.homs))
    
    ## Relabel using the new homologue names:
    new.homs[match(oldhom.vect,old.homs)]
  }))
  
  LGhom.df <- LGhom.df[order(LGhom.df$LG, LGhom.df$homologue),]
  LGhom.df$homologue <-  as.factor(LGhom.df$homologue)
  LGhom.df$LG <-  as.factor(LGhom.df$LG)
  
  if(ncol(table(LGhom.df$LG,LGhom.df$homologue)) > 2){
    warning("More than 2 couplings-groupings detected in some linkage groups. Deleting these...\nPerhaps try cluster_SN_markers at a higher LOD score?")
    
    message("Incorrect solution:")
    print(table(LGhom.df$LG,LGhom.df$homologue))
    
    affected_lg <- unique(as.character((LGhom.df[!LGhom.df$homologue %in% c(1,2),"LG"])))
    
    for(lg in affected_lg){
      temp <- LGhom.df[LGhom.df$LG == lg,]
      save_homs <- names(sort(table(temp$homologue),decreasing = T)[1:2])
      unlinked_markers <- c(unlinked_markers,temp[!temp$homologue %in% save_homs,"SxN_Marker"])
      LGhom.df <- LGhom.df[-which(LGhom.df$LG == lg & !LGhom.df$homologue %in% save_homs),]
    }
    
    LGhom.df <- droplevels(LGhom.df)
  }
  
  # write unlinked markers to standard out
  if (length(unlinked_markers) > 0) {
    write(paste0("\n####SxN Marker(s) lost in clustering step at LOD ",LOD_chm,":\n"), log.conn)
    unl.m <- vector.to.matrix(unlinked_markers, 4)
    write(knitr::kable(unl.m), log.conn)
  }
  if (!is.null(log))
    close(log.conn)
  
  return(LGhom.df)
} #phase_SN_diploid


#' Plot homologue position versus integrated positions
#' @param map_df A dataframe of a map that defines a linkage group.
#' @param maplist_homologue A list of maps were each item represents a homoloogue.
#' @examples
#' data("integrated.maplist", "maplist_P1_subset")
#' colnames(integrated.maplist[["LG2"]]) <- c("marker", "position", "QTL_LOD")
#' plot_hom_vs_LG(map_df = integrated.maplist[["LG2"]],
#'                maplist_homologue = maplist_P1_subset[["LG2"]])
#' @export
plot_hom_vs_LG <- function(map_df,
                           maplist_homologue){
  colors <- RColorBrewer::brewer.pal(6, "Dark2")
  homlengths <- sapply(maplist_homologue, function(x) max(x$position))
  plot(0, type = "n", xlim = c(0, max(map_df$position)), ylim = c(0, max(homlengths)),
       xlab = "Position on linkage group", ylab = "Position on homologue",
       cex.lab = 1.3, cex.axis = 1.2)
  legend("topleft", legend = names(maplist_homologue), col = colors, pch = 19, bty = "n")
  
  for(i in seq(length(maplist_homologue))){
    merged <- merge(map_df, maplist_homologue[[i]], by = "marker")
    cor_pos <- cor(merged$position.x, merged$position.y)
    if(cor_pos < 0) merged$position.y <- max(merged$position.y)-merged$position.y
    points(merged$position.x, merged$position.y, col = colors[i])
  }
}


#' Plot r versus LOD grouped by phase
#' @description \code{plot_linkage_df} plots r versus LOD, colour separated for different phases.
#' @param linkage_df A linkage data.frame as output of \code{\link{linkage}}.
#' @param r_max Maximum r value to plot
#' @param stepsize Size of window in recombination frequency to generate mean LOD score. Larger values will reduce plot density.
#' @param add_legend Logical, should a legend be added to the plot?
#' @param \dots Arguments passed to base plot function
#' @examples
#' data("SN_SN_P1")
#' plot_linkage_df(SN_SN_P1)
#' @export
plot_linkage_df <- function(linkage_df,
                            r_max = 0.5,
                            stepsize = 0.001,
                            add_legend = TRUE,
                            ...){
  
  all_phases <- levels(as.factor(linkage_df$phase))
  
  ## Assume there are at most 9 phasings..
  colours <-
    c(
      "limegreen",
      "red3",
      "darkgoldenrod2",
      "darkorchid4",
      "dodgerblue2",
      "gray25",
      "darkmagenta",
      "darkorange1",
      "cyan"
    )
  
  ## There can be up to 9 phases with DxD + DxD pairs
  linkage_df <- linkage_df[linkage_df[,"r"] <= r_max,]
  
  with(linkage_df,
       plot(
         NULL,
         xlim = c(0, 0.5),
         ylim = c(0, max(LOD[complete.cases(LOD)])),
         xlab = "r",
         ylab = "LOD",
         ...
       ))
  
  
  for (p in 1:length(all_phases)) {
    phase_level <- as.character(all_phases[p])
    
    temp.df <- linkage_df[linkage_df$phase == phase_level & complete.cases(linkage_df$LOD),]
    
    rseq <- seq(-0.01,0.51,stepsize)
    bins <- findInterval(temp.df$r,rseq)
    
    LOD.means <- tapply(temp.df$LOD,bins,mean)
    
    points(rseq[as.numeric(names(LOD.means))], LOD.means, 
           col = colours[p],...)
    
  }
  
  if(add_legend) legend("topright",
                        col = colours[1:length(all_phases)],
                        pch = 1,
                        as.character(all_phases),
                        cex=0.8)
  
} #plot_linkage_df 


#' Plot linkage maps
#' @description Makes a simple plot of a list of generated linkage maps
#' @param maplist A list of maps. In the first column marker names and in the second their position.
#' @param highlight A list of the same length of maplist with vectors of length 2 that specifies the
#' limits in cM from and to which the plotted chromosomes should be highlighted.
#' @param bg_col The background colour of the map.
#' @param highlight_col The color of the highlight. Only used if \code{highlight} is specified.
#' @param colname_in_mark Optional. The column name of the value to be plotted as marker color.
#' @param colname_beside_mark Optional. The column name of the value to be plotted beside the markers.
#' @param palette_in_mark,palette_beside_mark Color palette used to plot values. Only used if colnames of the values are specified.
#' @param color_by_type Logical. Should the markers be coloured by type? If TRUE, dosage_matrix should be specified.
#' @param dosage_matrix Optional (by default \code{NULL}). Dosage matrix of marker genotypes, input of \code{\link{linkage}}
#' @param parent1 Character string specifying the first (usually maternal) parentname.
#' @param parent2 Character string specifying the second (usually paternal) parentname.
#' @param legend Logical. Should a legend be drawn?
#' @param legend.x Optional. The x value of the coordinates of the legend.
#' @param legend.y Optional. The y value of the coordinates of the legend.
#' @param \dots Arguments passed to \code{\link{plot}}
#' @examples
#' data("maplist_P1")
#' plot_map(maplist = maplist_P1, colname_in_mark = "nnfit", bg_col = "white",
#'          palette_in_mark = colorRampPalette(c("blue", "purple", "red")),
#'          highlight = list(c(20, 60),
#'          c(60,80),
#'          c(20,30),
#'          c(40,70),
#'          c(60,80)))
#' @export
plot_map <- function(maplist,
                     highlight = NULL,
                     bg_col = "grey",
                     highlight_col="yellow",
                     colname_in_mark = NULL,
                     colname_beside_mark = NULL,
                     palette_in_mark = colorRampPalette(c("white", "purple")),
                     palette_beside_mark = colorRampPalette(c("white", "green")),
                     color_by_type = FALSE,
                     dosage_matrix = NULL,
                     parent1 = "P1",
                     parent2 = "P2",
                     legend = FALSE,
                     legend.x = 1,
                     legend.y = 120,
                     ...) {
  if(!is.null(dosage_matrix)) dosage_matrix <- test_dosage_matrix(dosage_matrix)
  
  if(color_by_type & is.null(dosage_matrix)) stop("If color_by_type = TRUE, dosage_matrix should be specified")
  map_lengths <- sapply(maplist, function(x)
    max(x[,2]))
  
  if(color_by_type){
    all_markers <- unlist(sapply(maplist, function(x) as.character(x[,1])))
    markpres <- all_markers %in% rownames(dosage_matrix)
    types <- dosage_matrix[all_markers[markpres], c(parent1, parent2)]
    types <- paste(types[, parent1], types[, parent2], sep = "x")
    if(sum(!markpres) > 0){
      types <- c(types, rep("nxn", sum(!markpres)))
      warning(paste0(all_markers[!markpres], " not in dosage_matrix so type could not be determined"))
    }
    types <- as.factor(types)
    discrete_cols <- rainbow(length(levels(types)))
    markercols <- discrete_cols[types]
    names(markercols) <- all_markers
    
  }
  
  plot(
    0, type = "n", xlim = c(0, length(maplist)+2),
    ylim = c(1.1*max(map_lengths), -0.05*max(map_lengths)),
    xaxt = "n", bty = "n", ylab = "Location (cM)", xlab = ""#, ...
  )
  axis(side = 1, at = seq(length(maplist)), labels = names(maplist))
  Hmisc::minor.tick(nx=1,ny=4)
  
  if(!is.null(highlight) & (length(highlight) != length(maplist))){
    warning("Highlight is not the same length as maplist. Chromosomes are not highlighted")
    highlight <- NULL
  }
  
  get_colors_and_legend <- function(maplist, colname, palette, map_lengths, col_colname, y_place){
    
    if(length(maplist)>1){
      scores <- sapply(X = maplist, simplify = FALSE, FUN = function(x)
        x[,c(colname)])
      stacked_scores <- stack(scores)
      maxscore <- max(stacked_scores$values, na.rm = TRUE)
      
    } else {
      stacked_scores <- maplist[[1]][,colname]
      stacked_scores <- data.frame(ind = names(maplist), values = stacked_scores)
      maxscore <- max(stacked_scores$values)
    }
    
    maxl <- max(map_lengths)
    place <- length(maplist)+1
    colors_scale <- palette(6)
    stacked_scores$col_scale <- colors_scale[as.numeric(cut(stacked_scores$values,breaks = 6,
                                                            include.lowest = TRUE))]
    for(mapn in names(maplist)){
      maplist[[mapn]][,col_colname] <- stacked_scores$col_scale[stacked_scores$ind==mapn]
    }
    for(i in 1:5){
      rect(place - 0.15, (maxl/25*(i+y_place)), place + 0.15, maxl/25*(i-1+y_place),
           col = colors_scale[i+1])
      text(place+0.1,  maxl/25*(i+y_place), round(i/5*maxscore, 1), pos = 4, cex = 0.8)
    }
    return(maplist)
  }
  
  
  
  if(!is.null(colname_beside_mark)){
    maplist <- get_colors_and_legend(maplist,
                                     colname_beside_mark,
                                     palette_beside_mark,
                                     map_lengths,
                                     col_colname = "col_beside",
                                     y_place = 0)
  }
  
  if(!is.null(colname_in_mark)){
    maplist <- get_colors_and_legend(maplist,
                                     colname_in_mark,
                                     palette_in_mark,
                                     map_lengths,
                                     col_colname = "col_in",
                                     y_place = 6)
  }
  
  for (mapn in seq(length(maplist))) {
    
    map <- maplist[[mapn]]
    
    if(color_by_type){
      markercols_mapn <- markercols[as.character(map[,1])]
    } else if(!is.null(colname_in_mark)){
      markercols_mapn <- map$col_in
    } else {
      markercols_mapn <- rep("black", nrow(map))
    }
    
    y <- map[,2]
    symbols(x = mapn, y = 0, circles= 0.25, bg=bg_col, add=TRUE, inches = FALSE)
    symbols(x = mapn, y = max(y), circles= 0.25, bg=bg_col, add=TRUE, inches = FALSE)
    rect(-0.25 + mapn, max(y), 0.25 + mapn, 0, col = bg_col)
    if(!is.null(highlight)){
      rect(-0.25 + mapn, highlight[[mapn]][1], 0.25 + mapn, highlight[[mapn]][2], col=highlight_col)
    }
    for (i in seq(length(y))){
      lines(c(-0.25 + mapn, 0.25 + mapn), c(y[i], y[i]), col=markercols_mapn[i])
    }
    lines(c(rep(-0.25+mapn,2)), c(max(y),0))
    lines(c(rep(0.25+mapn,2)), c(max(y),0))
    
    if(!is.null(colname_beside_mark)){
      
      for (i in seq(length(y))){
        lines(c(0.3 + mapn, 0.5+mapn), c(y[i], y[i]), col = map$col_beside[i], lwd = 2)
      }
    }
  }
  
  if(legend) {
    if(!color_by_type) stop("Doesn't make sense to generate legend while color_by_type = FALSE")
    # par(xpd=NA) #Allow legend to be outside the plot
    legend(legend.x, legend.y, legend = levels(types), col = discrete_cols,#ncol = length(levels(types)),
           lty = 1, lwd =2)
    # par(xpd=FALSE)
  }
}


#' Visualise the phased homologue maplist
#' @description \code{plot_phased_maplist} is a function for visualising a phased maplist, the output of
#' \code{\link{create_phased_maplist}}
#' @param phased.maplist A list of phased linkage maps, the output of \code{\link{create_phased_maplist}}
#' @param ploidy Integer. Ploidy of the organism.
#' @param ploidy2 Optional integer, by default \code{NULL}. Ploidy of parent 2, if different from parent 1.
#' @param cols Vector of colours for the integrated, parent1 and parent2 maps, respectively.
#' @param width Width of the linkage maps, by default 0.2
#' @param mapTitles Optional vector of titles for maps, by default names of maplist, or titles LG1, LG2 etc. are used.
#' @examples
#' data("phased.maplist")
#' plot_phased_maplist(phased.maplist, ploidy = 4)
#' @export
plot_phased_maplist <- function(phased.maplist,
                                ploidy,
                                ploidy2 = NULL,
                                cols = c("black","darkred","navyblue"),
                                width = 0.2,
                                mapTitles = NULL) {
  
  if(is.null(ploidy2)) ploidy2 <- ploidy
  
  if(is.null(mapTitles)) {
    if(is.null(names(phased.maplist))) {
      mapTitles <- paste0("LG",seq(length(phased.maplist)))
    } else{
      mapTitles <- names(phased.maplist)
    }
  }
  
  for(lg in seq(length(phased.maplist))){
    phased.maplist[[lg]][1,]
    
    plot(NULL,xlim = c(0.5,ploidy+ploidy2+1),
         ylim=c(max(phased.maplist[[lg]]$position)+width, -width),
         xlab="",ylab="cM",axes=FALSE, main = mapTitles[lg])
    
    axis(1,at = c(1:ploidy,(ploidy+2):(ploidy+ploidy2+1)),
         labels = c(paste0("h",1:ploidy),paste0("h",(ploidy+1):(ploidy+ploidy2))))
    axis(2,las=1)
    Hmisc::minor.tick(nx=1,ny=4)
    
    ## Add the central chromosome:
    maxcM <- max(phased.maplist[[lg]]$position)
    mincM <- min(phased.maplist[[lg]]$position)
    symbols(x = ploidy + 1, y = mincM, circles= width, bg="white", add=TRUE, inches = FALSE)
    symbols(x = ploidy + 1, y = maxcM, circles= width, bg="white", add=TRUE, inches = FALSE)
    rect(-width + ploidy + 1, maxcM, width + ploidy + 1, mincM, col = "white")
    # segments(5-width+0.02,0,5+width-0.02,0,col="white")
    
    for (i in 1:nrow(phased.maplist[[lg]])){
      lines(c(-width + ploidy + 1, width + ploidy + 1),
            c(phased.maplist[[lg]]$position[i], phased.maplist[[lg]]$position[i]),
            col=cols[1])
    }
    
    ## Fill in the homologue marker distribution:
    for(cl in 1:ploidy){
      hom.pts <- which(phased.maplist[[lg]][,paste0("h",cl)] == 1)
      
      for(i in hom.pts) lines(c(-width + cl ,width + cl),
                              c(phased.maplist[[lg]][i,"position"],
                                phased.maplist[[lg]][i,"position"]),
                              col = cols[2])
      
    }
    
    for(cl in (ploidy + 1):(ploidy + ploidy2)){
      hom.pts <- which(phased.maplist[[lg]][,paste0("h",cl)] == 1)
      
      for(i in hom.pts) lines(c(-width + cl + 1, width + cl + 1),
                              c(phased.maplist[[lg]][i,"position"],
                                phased.maplist[[lg]][i,"position"]), col = cols[3])
    }
    
  }
  
} #plot_phased_maplist


#' Plot r versus LOD
#' @description \code{r_LOD_plot} plots r versus LOD, colour separated for different phases.
#' @param linkage_df A linkage data.frame as output of \code{\link{linkage}}.
#' @param plot_main A character string specifying the main title
#' @param chm Integer specifying chromosome
#' @param r_max Maximum r value to plot
#' @examples
#' data("SN_SN_P1")
#' r_LOD_plot(SN_SN_P1)
#' @export
r_LOD_plot <- function(linkage_df,
                       plot_main = "",
                       chm = NA,
                       r_max = 0.5) {
  all_phases <- levels(as.factor(linkage_df$phase))
  
  ## Assume there are at most 6 phasings..
  colours <-
    c(
      "limegreen",
      "red3",
      "darkgoldenrod2",
      "darkorchid4",
      "dodgerblue2",
      "gray25",
      "darkmagenta",
      "darkorange1",
      "cyan"
    )
  ## There can be up to 9 phases with DxD + DxD pairs
  
  
  if (!is.na(chm)) {
    plot_title <- paste(plot_main, "   -    chm", chm)
  } else{
    plot_title <- plot_main
  }
  linkage_df <- linkage_df[linkage_df[,"r"]<=r_max,]
  
  with(linkage_df,
       plot(
         NULL,
         xlim = c(0, 0.5),
         ylim = c(0, max(LOD[complete.cases(LOD)])),
         xlab = "r",
         ylab = "LOD",
         main = plot_title
       ))
  for (p in 1:length(all_phases)) {
    phase_level <- as.character(all_phases[p])
    phase_col <- colours[p]
    
    temp_data <- linkage_df[linkage_df$phase == phase_level,]
    
    with(temp_data[complete.cases(temp_data$LOD),],
         points(r, LOD, col = phase_col))
  }
  legend("topright",
         col = colours[1:length(all_phases)],
         pch = 1,
         as.character(all_phases))
  
}

#' Screen for duplicate individuals
#' @description \code{screen_for_duplicate_individuals} identifies and merges duplicate individuals.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @param cutoff Correlation coefficient cut off. At this correlation coefficient, individuals are merged. If NULL user input will be asked after plotting.
#' @param plot_cor Logical. Should correlation coefficients be plotted? Can be memory/CPU intensive with high number of individuals.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @param saveCOV = "PearsonCC"
#' @return A matrix similar to dosage_matrix, with merged duplicate individuals.
#' @examples
#' \dontrun{
#' #user input:
#' data("segregating_data")
#' screen_for_duplicate_individuals(dosage_matrix=segregating_data,cutoff=0.9,plot_cor=TRUE)
#' }
#' @export
screen_for_duplicate_individuals <-
  function(dosage_matrix,
           cutoff = NULL,
           plot_cor = T,
           log = NULL,
           saveCOV = "PearsonCC") {
    dosage_matrix <- test_dosage_matrix(dosage_matrix)
    cor.dosage_matrix <-
      cor(dosage_matrix, use = "pairwise.complete.obs")
    diag(cor.dosage_matrix) <- NA
    cor.dosage_matrix[lower.tri(cor.dosage_matrix)] <- NA
    
    if (is.null(log)) {
      log.conn <- stdout()
    } else {
      matc <- match.call()
      write.logheader(matc, log)
      log.conn <- file(log, "a")
    }
    
    if (plot_cor) {
      corvec <- c(cor.dosage_matrix[!is.na(cor.dosage_matrix)])
      
      plot(
        corvec,
        ylab = "Pearson's correlation coefficient",
        xlab = "",
        ylim = c(0, 1),
        pch = 20,
        bty = "n",
        col = rgb(0.5, 0.5, 0.5, 0.4),
        xaxt = "n",
        cex.lab = 1.2,
        cex.axis = 1.2
      )
      if (!is.null(cutoff)) {
        abline(
          h = cutoff,
          col = "grey",
          lwd = 2,
          lty = 2
        )
        points(which(corvec > cutoff),
               corvec[corvec > cutoff],
               col = rgb(1, 0, 0, 0.4),
               pch = 20)
        p1 <- colnames(cor.dosage_matrix)[1]
        p2 <- colnames(cor.dosage_matrix)[2]
        if(p1 == "P1" & p2 == "P2")
          points(1, corvec[1], col = rgb(0.5, 0, 0.5, 0.4), pch = 20)
      }
    }
    
    if (is.null(cutoff)) {
      cutoff <- as.numeric(readline("Input similarity threshold:  "))
      while (cutoff > 1 | is.na(cutoff) | cutoff < 0) {
        cutoff <-
          as.numeric(readline("Invalid threshold. Re-enter similarity threshold:  "))
      }
      write(paste0("\nInput similarity threshold: ", cutoff, "\n"),
            log.conn)
      if (plot_cor) {
        abline(
          h = cutoff,
          col = "grey",
          lwd = 1.2,
          lty = 2
        )
        points(which(corvec > cutoff),
               corvec[corvec > cutoff],
               col = rgb(1, 0, 0, 0.4),
               pch = 20)
      }
    }
    
    combsAboveCutoff <-
      which(cor.dosage_matrix > cutoff, arr.ind = T)
    nw <- igraph::graph.data.frame(combsAboveCutoff, directed = F)
    gcl_nw <- igraph::groups(igraph::clusters(nw))
    
    if(length(unique(unlist(gcl_nw))) >= ncol(dosage_matrix) - 2) stop("At this threshold, whole population would be merged!")
    
    if (length(gcl_nw) > 0) {
      remove <- c()
      for (group in gcl_nw) {
        group <- colnames(dosage_matrix)[as.numeric(group)]
        
        if(length(group) > 2) {
          warning("Multiple duplicates of single genotype identified at this threshold. Attempting to merge...")
          
          for(r in 2:length(group)){
            write(paste(c(
              "\nCombining",group[1],"&",group[r],"into",group[1]
            ), collapse = " "), log.conn)
            
            comp_mat <- dosage_matrix[, group[c(1,r)]]
            
            merged <- apply(comp_mat, 1, function(x) {
              ifelse(all(is.na(x)), NA,
                     ifelse(length(unique(x[!is.na(x)])) == 1, x[!is.na(x)][1], NA))
            })
            
            remove <- c(remove, group[r])
            
          }
          
        } else{ #only 2 members in this group, so a true duplicate example
          write(paste(c(
            "\nCombining",group[1],"&",group[2],"into",group[1]
          ), collapse = " "), log.conn)
          
          comp_mat <- dosage_matrix[, group]
          merged <- apply(comp_mat, 1, function(x) {
            ifelse(all(is.na(x)), NA,
                   ifelse(length(unique(x[!is.na(x)])) == 1, x[!is.na(x)][1], NA))
          })
          remove <- c(remove, group[2:length(group)])
        }
        
      }
      write(paste("\n####",length(remove),"individuals removed:\n"), log.conn)
      remove.m <- vector.to.matrix(remove, 4)
      if(!is.null(log)) message(paste("\n####",length(remove),"individuals removed:\n"));print(remove.m)
      write(knitr::kable(remove.m), log.conn)
      dosage_matrix <-
        dosage_matrix[,!colnames(dosage_matrix) %in% remove]
    } else {
      write("\nNo duplicates found\n", log.conn)
    }
    
    if (!is.null(log))
      close(log.conn)
    
    if (!is.null(saveCOV)){
      number <- length(corvec)
      variation <- var(corvec)
      rangec <- paste0(round(range(corvec),digits = 2), collapse = "_")
      CovFile <- data.frame("Number" = number,
                            "Variation" = variation,
                            "Range" = rangec)
      write.csv(CovFile, paste0(saveCOV,".csv"))
    }
    
    return(dosage_matrix)
  } #screen_for_duplicate_individuals


#' Screen for duplicate individuals using weighted genotype probabilities
#' @description \code{screen_for_duplicate_individuals.gp} identifies and merges duplicate individuals based on probabilistic genotypes.
#' See \code{\link{screen_for_duplicate_individuals}} for the original function.
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
#' Most probable dosage for a particular individual and marker combination, if \code{maxP} exceeds a user-defined threshold (e.g. 0.9), otherwise \code{NA}
#' }
#' }
#' @param ploidy The ploidy of parent 1
#' @param parent1 character vector with the sample names of parent 1
#' @param parent2 character vector with the sample names of parent 2
#' @param F1 character vector with the sample names of the F1 individuals
#' @param cutoff Correlation coefficient cut off to declare duplicates. At this correlation coefficient, individuals are merged. If \code{NULL} user input will be asked after plotting.
#' @param plot_cor Logical. Should correlation coefficients be plotted? Can be memory/CPU intensive with high number of individuals.
#' @param log Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.
#' @param saveCOV A file name where the Pearson's correlation coefficient's variation, number, and mean can be saved
#' @return A data frame similar to input \code{probgeno_df}, but with duplicate individuals merged. 
#' @export
screen_for_duplicate_individuals.gp <- function(probgeno_df,
                                                ploidy,
                                                parent1 = "P1",
                                                parent2 = "P2",
                                                F1,
                                                cutoff = 0.95,
                                                plot_cor = TRUE,
                                                saveCOV = "PeasonCC",
                                                log = NULL){
  probgeno_df <- test_probgeno_df(probgeno_df)
  
  #make all threshold add up to 1
  sum <- rep(0, nrow(probgeno_df))
  for(p in 0:ploidy){
    poly <- paste0("P",p)
    s <- probgeno_df[[poly]]
    sum <- sum + s
  }
  
  #get the proba probgeno_df
  t <- rep(0, nrow(probgeno_df))
  for(p in 0:ploidy){
    poly <- paste0("P",p)
    s <- probgeno_df[[poly]]*p/sum
    t <- t + s
  }
  probgeno_df$probascores <- t
  
  dosages <- do.call(rbind, lapply(as.character(unique(probgeno_df$MarkerName)),
                                   function(m) c(
                                       mean(probgeno_df[probgeno_df$MarkerName == m & probgeno_df$SampleName %in% parent1,"probascores"]),
                                       mean(probgeno_df[probgeno_df$MarkerName == m & probgeno_df$SampleName %in% parent2,"probascores"]),
                                       probgeno_df[probgeno_df$MarkerName == m & probgeno_df$SampleName %in% F1,"probascores"]
                                     )))
 
  colnames(dosages) <- c("P1", "P2", F1)
  rownames(dosages) <- as.character(unique(probgeno_df$MarkerName))
  
  check2 <- screen_for_duplicate_individuals(dosage_matrix = as.matrix(dosages),
                                                  cutoff = cutoff, 
                                                  plot_cor = plot_cor, 
                                                  saveCOV = saveCOV,
                                                  log = log)
  ind <- colnames(check2)
  
  if(all(c("P1", "P2") %in% ind)){
    pg1 <- probgeno_df[probgeno_df$SampleName %in% c(ind,parent1,parent2),]
  }else{
    writeLines("There is something wrong with the parents. Please check it!")
    pg1 <- NULL
  }
  return(pg1)
} #screen_for_duplicate_individuals.gp



#' Screen for and remove duplicated markers
#' @description \code{screen_for_duplicate_markers} identifies and merges duplicate markers.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @param merge_NA Logical. Should missing values be imputed if non-NA in duplicated marker? By default, \code{TRUE}.
#' If \code{FALSE} the dosage scores of representing marker are represented in the filtered_dosage_matrix.
#' @param plot_cluster_size Logical. Should an informative plot about duplicate cluster size be given? By default, \code{TRUE}.
#' @param ploidy Ploidy level of parent 1. Only needed if \code{estimate_bin_size} is \code{TRUE}
#' @param ploidy2 Integer, by default \code{NULL}. If parental ploidies differ, use this to specify the ploidy of parent 2. 
#' Only needed if \code{estimate_bin_size} is \code{TRUE}
#' @param LG_number Expected number of chromosomes (linkage groups). Only needed if \code{estimate_bin_size} is \code{TRUE}
#' @param estimate_bin_size Logical, by default \code{FALSE}. If \code{TRUE}, a very rudimentary calculation is made to estimate
#' the average size of a marker bin, assuming a uniform distribution of cross-over events and on average one cross-over per bivalent.  
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return
#' A list containing:
#' \itemize{
#' \item{bin_list} {list of binned markers. The list names are the representing markers.
#' This information can later be used to enrich the map with binned markers.}
#' \item{filtered_dosage_matrix} {dosage_matrix with merged duplicated markers.
#' The markers will be given the name of the marker with least missing values.}
#' }
#' @examples
#' data("screened_data3")
#' dupmscreened <- screen_for_duplicate_markers(screened_data3)
#' @export
screen_for_duplicate_markers <- function(dosage_matrix,
                                         merge_NA = TRUE,
                                         plot_cluster_size = TRUE,
                                         ploidy,
                                         ploidy2 = NULL,
                                         LG_number,
                                         estimate_bin_size = FALSE,
                                         log = NULL){
  dosage_matrix <- test_dosage_matrix(dosage_matrix)
  
  rn_orig <- rownames(dosage_matrix)
  ddo <- dosage_matrix[do.call(order, as.data.frame(dosage_matrix)),]
  clmat <- matrix(integer(), ncol =2)
  for(i in seq(nrow(ddo)-1)){
    diff <- ddo[i,] - ddo[i+1,]
    if(all(diff[!is.na(diff)]==0))
      clmat <- rbind(clmat, c(i, i+1))
  }
  
  nw <- igraph::graph_from_data_frame(d= as.data.frame(clmat), directed = FALSE)
  gcl_nw <- igraph::groups(igraph::clusters(nw))
  
  if(length(gcl_nw) == 0) { #No groups detected (no duplicates), effective abort:
    return(list(bin_list = NULL, filtered_dosage_matrix = dosage_matrix))
  } else{
    
    if(plot_cluster_size) plot(sapply(gcl_nw, length),
                               ylab = "Cluster size",
                               cex.lab = 1.2, cex.axis = 1.2,
                               col = rgb(1,0,0, alpha = 0.5),
                               pch = 20)
    
    
    dupnums <- as.integer(do.call(c, gcl_nw))
    rn <- rownames(ddo)
    ddo_wo_dups <- ddo[-dupnums,]
    bin_list <- list()
    
    mergeFun1 <- function(sub, repr){
      nwdos <- colMeans(sub, na.rm = TRUE)
      return(nwdos)
    }
    
    mergeFun <- function(sub, repr){
      nwdos <- sub[repr,]
      return(nwdos)
    }
    
    if(merge_NA){
      mergeFun <- mergeFun1
    }
    
    nwdos_list <- list()
    
    if(length(gcl_nw) > 1) pb <- txtProgressBar(min = 1, max = length(gcl_nw), style = 3)
    
    for(i in seq(length(gcl_nw))){
      sub <- ddo[as.integer(gcl_nw[[i]]),]
      narate <- apply(sub, 1, function(x) sum(is.na(x)))
      repr <- names(narate)[which.min(narate)]
      nwdos <- mergeFun(sub, repr)
      nwdos_list[[repr]] <- nwdos
      bin_list[[repr]] <- names(narate)[-(which.min(narate))]
      if(length(gcl_nw) > 1) setTxtProgressBar(pb, i)
    }
    
    nwdos_mat <- do.call(rbind, nwdos_list)
    outmat <- rbind(nwdos_mat, ddo_wo_dups)
    leftm <- rn_orig[rn_orig %in% rownames(outmat)]
    outmat <- outmat[leftm,]
    class(outmat) <- "integer"
    
    if (is.null(log)) {
      log.conn <- stdout()
    } else {
      matc <- match.call()
      write.logheader(matc, log)
      log.conn <- file(log, "a")
    }
    
    write(paste0("\nMarked and merged ", nrow(dosage_matrix)-nrow(outmat), " duplicated markers"),
          file = log.conn)
    

    #Estimate the expected distribution of bin sizes
    if(estimate_bin_size){
      
      if(is.null(ploidy2)) ploidy2 <- ploidy
          
      pop.size <- ncol(dosage_matrix) - 2
      nmark <- nrow(dosage_matrix)
      
      #estimate number of cross-over breaks across the genome and the whole population
      nbreaks <- LG_number*pop.size*(ploidy/2 + ploidy2/2)
      
      mds <- marker_data_summary(dosage_matrix = dosage_matrix,
                                 parent1 = colnames(dosage_matrix)[1],
                                 parent2 = colnames(dosage_matrix)[2],
                                 pairing = "random",#doesn't matter
                                 ploidy = ploidy,
                                 shortform = T, 
                                 verbose = F)
      
      ntype <- length(mds$parental_info[mds$parental_info!= 0]) #total nr of marker segregation types in the data
      
      bin.size <- round(nmark/(ntype*nbreaks),2)
      
      if(bin.size > 2) {
        write(paste("At a very rough estimate, there should be about",round(bin.size) - 1,"markers in each bin, although this assumes the distribution of markers as well as the distribution of recombination break-points are uniform across the genome."),file = log.conn)
      } else if (bin.size > 1 & bin.size <= 2){
        write("Your population size and number of markers appear to be relatively well-balanced.", file = log.conn)
      } else{
        write("Your marker set appears to be somewhat small to fully exploit the genetic resolution afforded by your mapping population!", file = log.conn)
      }
      
      bs <- data.frame("pop.size" = pop.size,
                       "nmark" = nmark,
                       "N_LG" = LG_number,
                       "P1_xo" = LG_number*ploidy*pop.size/2,
                       "P2_xo" = LG_number*ploidy2*pop.size/2,
                       "Tot_xo" = nbreaks,
                       "n_marker_types" = ntype,
                       "bin.size" = bin.size)
      
      write(knitr::kable(bs), file = log.conn)
    }
    
    
    if (!is.null(log))
      close(log.conn)
    
    
    return(list(bin_list = bin_list, filtered_dosage_matrix = outmat))
  }
  
}


#' Screen marker data for NA values
#' @description \code{screen_for_NA_values} identifies and can remove rows or columns of a marker dataset based on the relative frequency of missing values.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @param margin An integer at which margin the missing value frequency will be calculated. A value of 1 means rows (markers), 2 means columns (individuals)
#' @param cutoff Missing value frequency cut off. At this frequency, rows or columns are removed from the dataset. If NULL user input will be asked after plotting the missing value frequency histogram.
#' @param parentnames A character vector of length 2, specifying the parent names.
#' @param plot_breakdown Logical. Should the percentage of markers removed as breakdown per markertype be plotted? Can only be used if margin = 1.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @param print.removed Logical. Should removed instances be printed?
#' @return A matrix similar to dosage_matrix, with rows or columns removed that had a higher missing value frequency than specified.
#' @examples
#' data("segregating_data","screened_data")
#' screened_markers<-screen_for_NA_values(dosage_matrix=segregating_data, margin=1, cutoff=0.1)
#' screened_indiv<-screen_for_NA_values(dosage_matrix=screened_data, margin=2, cutoff=0.1)
#' @export
screen_for_NA_values <- function(dosage_matrix,
                                 margin = 1,
                                 cutoff = NULL,
                                 parentnames = c("P1", "P2"),
                                 plot_breakdown = FALSE,
                                 log = NULL,
                                 print.removed = TRUE) {
  dosage_matrix <- test_dosage_matrix(dosage_matrix)
  parent_cols <- which(colnames(dosage_matrix) %in% parentnames)
  parent_data <- dosage_matrix[, parent_cols]
  input_data <- dosage_matrix[,-parent_cols]
  names <- dimnames(input_data)[[margin]]
  nx <- dim(input_data)[margin]
  ny <- dim(input_data)[3 - margin]
  
  if(margin==1){
    missing <- rowSums(is.na(input_data))
    object <- "markers"
  } else if (margin ==2){
    missing <- colSums(is.na(input_data))
    object <- "individuals"
  } else {
    stop("Margin should be either 1 (markers) or 2 (individuals).")
  }
  
  fraction_missing <- missing / ny
  
  ## Finding an empirical cut-off for the allowable rate of missing values:
  cut_offs <- seq(0.001, ifelse(is.null(cutoff),0.2,max(0.2,cutoff)), by = 0.001)
  
  percent_cut <- sapply(cut_offs,
                        function(x) {
                          sum(fraction_missing >= x)
                        })
  
  plot(
    cut_offs,
    percent_cut,
    type = "l",
    col = "darkgrey",
    lwd = 2,
    ylab = paste("number of", object,  "removed"),
    xlab = "Threshold for acceptable 'NA' rate"
  )
  five_proc_line <- sum(fraction_missing >= 0.05)
  ten_proc_line <- sum(fraction_missing >= 0.1)
  lines(x = c(0.05, 0.05, 0), y = c(0, rep(five_proc_line, 2)), col = 'darkolivegreen',
        lty = 2)
  lines(x = c(0.1, 0.1, 0), y = c(0, rep(ten_proc_line, 2)), col = 'firebrick',
        lty = 2)
  
  text(x = 0.08,
       y = five_proc_line,
       paste(five_proc_line, object),
       col = "darkolivegreen")
  text(x = 0.12,
       y = ten_proc_line,
       paste(ten_proc_line, object),
       col = "firebrick")
  
  if (is.null(cutoff)) {
    ## Take user input on the cutoff threshold (max rate of NA allowed)...
    print("Please specify the allowable cutoff rate...")
    acceptable_cutoff_rate <-
      as.numeric(readline("(e.g. by typing 0.1 it means 10% NA is the max allowed rate):   "))
  } else {
    acceptable_cutoff_rate <- cutoff
  }
  
  if(acceptable_cutoff_rate >= 0.5) stop("Suggest to choose a smaller cutoff rate.")
  
  user_defined_missing <- fraction_missing >= acceptable_cutoff_rate
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  if (sum(user_defined_missing) > 0) {
    # some markers are removed..
    write(paste(
      "\nNumber of removed",
      object,
      ": ",
      sum(user_defined_missing),
      "\n"
    ),
    log.conn)
    
    if (!is.null(log)) message(paste("\nNumber of removed",object,": ",sum(user_defined_missing),"\n"))
    
    if (margin == 1) {
      output_data <- input_data[!user_defined_missing,]
      parent_data <- parent_data[!user_defined_missing,]
      if (print.removed) {
        removed <- rownames(input_data)[user_defined_missing]
        removed.m <- vector.to.matrix(removed, n.columns = 4)
        write(knitr::kable(removed.m), log.conn)
        write("\n", log.conn)
      }
    } else {
      output_data <- input_data[,!user_defined_missing]
      if (print.removed) {
        removed <- colnames(input_data)[user_defined_missing]
        removed.m <- vector.to.matrix(removed, n.columns = 4)
        write(knitr::kable(removed.m), log.conn)
        write("\n", log.conn)
      }
    }
    
    n <- dim(output_data)[margin]
    
    write(paste("\nThere are now", n, object,"leftover in the dataset."),
          log.conn)
    
    if(!is.null(log)){
      message(paste("\nThere are now", n, object,"leftover in the dataset."))
      close(log.conn)
    }
    
    output_data <- cbind(parent_data, output_data)
    
    if(margin == 1 & plot_breakdown){
      bp1 <- factor(dosage_matrix[,parentnames[1]])
      bp2 <- factor(dosage_matrix[,parentnames[2]])
      ap1 <- factor(output_data[,parentnames[1]])
      ap2 <- factor(output_data[,parentnames[2]])
      levels(ap1) <- levels(bp1)
      levels(ap2) <- levels(bp2)
      table_before_NArm <- table(bp1, bp2)
      table_after_NArm <- table(ap1, ap2)
      percentage_removed <- (table_before_NArm - table_after_NArm)/table_before_NArm*100
      markertypes <-
        expand.grid(list(rownames(percentage_removed), colnames(percentage_removed)))
      parental_quantities <- as.vector(percentage_removed)
      names(parental_quantities) <-
        markertypes <- paste0(markertypes[, 1], "x", markertypes[, 2])
      parental_quantities[is.na(parental_quantities)] <- 0
      pq_plot <- parental_quantities[parental_quantities > 0]
      barplot(
        pq_plot,
        col = c("darkorchid"),
        ylab = "Percentage > 10% NA",
        cex.names = 0.8,
        las=2
      )} else {
        if (margin != 1 & plot_breakdown) warning("Plotting the breakdown per markertype can only be done when margin = 1, so is skipped")
      }
  } else{
    print("Nothing removed at this rate.")
    
    #write.csv(input_data,paste(outname,".csv",sep=""),row.names=F)
    output_data <- dosage_matrix
  }
  return(output_data)
}

#' Identify deviations in LOD scores between pairs of simplex x nulliplex markers
#' @description \code{SNSN_LOD_deviations} checks whether the LOD scores obtained in the case of pairs of simplex x nulliple
#' markers are compatible with expectation. This can help identify problematic linkage estimates which can adversely affect
#' marker clustering.
#' @param linkage_df A linkage data.frame as output of \code{\link{linkage}}.
#' @param ploidy Integer. The ploidy level of the species.
#' @param N Numeric. The number of F1 individuals in the mapping population.
#' @param plot_expected Logical. Plot the observed and expected relationship between r and LOD.
#' @param alpha Numeric. Vector of upper and lower tolerances around expected line.
#' @param phase Character string. Specify which phase to examine for deviations (usually this is "coupling" phase).
#' @return A vector of deviations in LOD scores outside the range defined by tolerances input \code{alpha}
#' @examples
#' data("SN_SN_P1")
#' SNSN_LOD_deviations(SN_SN_P1,ploidy = 4, N = 198)
#' @export
SNSN_LOD_deviations <- function(linkage_df,
                                ploidy,
                                N,
                                plot_expected=TRUE,
                                alpha = c(0.05,0.2),
                                phase=c("coupling","repulsion")){
  
  linkage_df<-linkage_df[linkage_df[,"phase"]==phase,]
  phase <- match.arg(phase)
  
  LODc <- function(r,N_o) N_o*((1-r)*log10(1-r) + r*log10(pmax(r,1e-6)) + log10(2))
  LODr <- function(r,N_o) N_o*(((ploidy/2-1+r)*log10(ploidy/2-1+r) +
                                  (ploidy/2-r)*log10(ploidy/2-r))/(ploidy-1) +
                                 log10(2/(ploidy-1)))
  
  upr <- if(phase=="coupling"){
    LODc(linkage_df$r,N*(1+alpha[1]))
  }else{
    LODr(linkage_df$r,N*(1+alpha[1]))
  }
  
  lwr <- if(phase=="coupling"){
    LODc(linkage_df$r,N*(1-alpha[2]))
  } else{
    LODr(linkage_df$r,N*(1-alpha[2]))
  }
  
  devs<-rep(0,nrow(linkage_df))
  upr.outliers <- linkage_df$LOD>upr
  lwr.outliers <- linkage_df$LOD<lwr
  outliers <- upr.outliers+lwr.outliers>0
  
  devs[upr.outliers] <- (linkage_df$LOD-upr)[upr.outliers]
  devs[lwr.outliers] <- (lwr-linkage_df$LOD)[lwr.outliers]
  
  devs[is.na(devs)] <- 0
  
  if(plot_expected){
    r_LOD_plot(linkage_df)
    s<-seq(0,0.5,0.01)
    
    if(phase=="coupling"){
      lines(s,LODc(s,N*(1+alpha[1])),lwd=2,lty=3)
      lines(s,LODc(s,N*(1-alpha[2])),lwd=2,lty=3)
    } else{
      lines(s,LODr(s,N*(1+alpha[1])),lwd=2,lty=3)
      lines(s,LODr(s,N*(1-alpha[2])),lwd=2,lty=3)
    }
    points(linkage_df$r[outliers],linkage_df$LOD[outliers],pch=8)
  }
  
  return(devs)
}



#' Check for and estimate preferential pairing
#' @description Identify closely-mapped repulsion-phase simplex x nulliplex markers and test these
#' for preferential pairing, including estimating a preferential pairing parameter.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @param maplist A list of integrated chromosomal maps, as generated by e.g. \code{\link{MDSMap_from_list}}. In the first column marker names and in the second their position.
#' @param LG_hom_stack A \code{data.frame} with markernames (\code{"SxN_Marker"}), linkage group (\code{"LG"}) and homologue (\code{"homologue"}),
#' the output of \code{\link{define_LG_structure}} or \code{\link{bridgeHomologues}} usually.
#' @param target_parent Character string specifying the parent to be tested for preferential pairing as provided in the columnnames of dosage_matrix, by default "P1".
#' @param other_parent The other parent, by default "P2"
#' @param ploidy The ploidy level of the species, by default 4 (tetraploid) is assumed.
#' @param min_cM The smallest distance to be considered a true distance on the linkage map, by default distances less than 0.5 cM are considered essentially zero.
#' @param adj.method Method to correct p values of Binomial test for multiple testing, by default the FDR correction is used, other options are available, inherited from \code{\link{p.adjust}}
#' @param verbose Should messages be sent to stdout? If \code{NULL} log is send to stdout.
#' @examples
#' data("ALL_dosages","integrated.maplist","LGHomDf_P1_1")
#' P1pp <- test_prefpairing(ALL_dosages,integrated.maplist,LGHomDf_P1_1,ploidy=4)
#' @export
test_prefpairing <- function(dosage_matrix,
                             maplist,
                             LG_hom_stack,
                             target_parent = "P1",
                             other_parent = "P2",
                             ploidy,
                             min_cM = 0.5,
                             adj.method = "fdr",
                             verbose = TRUE) {
  
  dosage_matrix <- test_dosage_matrix(dosage_matrix)
  LG_hom_stack <- LG_hom_stack[order(LG_hom_stack[,2],LG_hom_stack[,3]),]
  
  if(!target_parent %in% colnames(dosage_matrix) | !other_parent %in% colnames(dosage_matrix))
    stop("Incorrect column name identifiers supplied for parents (target_parents and/or other_parent). Please check!")
  
  maplist.df <- do.call("rbind", lapply(seq(length(maplist)), function(LG) cbind(maplist[[LG]][,c("marker","position")],LG)))
  
  dosage_data.1 <-
    dosage_matrix[dosage_matrix[, target_parent] == 1 & dosage_matrix[,other_parent] == 0,
                  -which(colnames(dosage_matrix) %in% c(target_parent,other_parent))]
  
  ## Use only marker data that has been assigned and mapped:
  temp.dosage <- as.data.frame(cbind(rownames(dosage_data.1),dosage_data.1))
  colnames(temp.dosage)[1] <- colnames(LG_hom_stack)[1] <- "marker"
  dosage_data.x <- merge(LG_hom_stack, temp.dosage, by = "marker", sort = FALSE)
  dosage_data.y <- merge(dosage_data.x, maplist.df, by = "marker", sort = FALSE)
  
  assignment_data <- dosage_data.y[,c(1,2,3,ncol(dosage_data.y)-1)]
  dosage_data.2 <- as.matrix(dosage_data.x[,4:ncol(dosage_data.x)])
  
  dosage_data <- matrix(as.numeric(dosage_data.2), ncol = ncol(dosage_data.2),
                        dimnames = list(assignment_data[,1],colnames(dosage_data.2)))
  
  n <- nrow(dosage_data)
  
  if(n == 0){
    if (verbose) message("Insufficient data to complete analysis - check whether LG_hom_stack is consistent with dosage_matrix and/or target_parent.")
    return(NULL)
  }
  
  ## Only use repulsion-phase linkages within min_cM distance:
  combinations <- do.call("rbind", lapply(unique(assignment_data$LG), function(ch) t(combn(which(assignment_data$LG == ch),2))))
  combinations <- as.data.frame(combinations[-union(
    which(assignment_data[combinations[,1],"homologue"] == assignment_data[combinations[,2],"homologue"]),
    which(abs(assignment_data[combinations[,1],"position"] - assignment_data[combinations[,2],"position"]) >= min_cM)
  ),])
  
  if(nrow(combinations) == 0) stop(paste("No suitable repulsion-phase pairs found. Suggest to increase min_cM (current value =", min_cM,"cM) slightly."))
  
  dosage_data <- t(dosage_data)
  
  dosage_combinations <- expand.grid(c(0,1), c(0,1))
  dosage_levels <- paste0("n_", dosage_combinations[, 1], dosage_combinations[, 2])
  rownames(dosage_combinations) <- dosage_levels
  
  udc_template <- unique(dosage_combinations[, 1])
  udc_compare <- unique(dosage_combinations[, 2])
  
  if (verbose){message(paste("In total",
                             nrow(combinations),
                             "repulsion-phase combinations to be tested for preferential pairing...\n"))}
  
  combs <- as.matrix(combinations)
  colnames(combs) <- c("marker_a", "marker_b")
  
  template_matrix <- dosage_data[, combs[, 1], drop = F]
  compare_matrix <- dosage_data[, combs[, 2], drop = F]
  
  count_mat <-
    matrix(integer(),
           nrow = ncol(template_matrix),
           ncol = length(dosage_levels))
  colnames(count_mat) <- dosage_levels
  
  matrix_list_template <- lapply(udc_template, function(x) {
    template_matrix == x
  })
  names(matrix_list_template) <- udc_template
  
  matrix_list_compare <- lapply(udc_compare, function(x) {
    compare_matrix == x
  })
  names(matrix_list_compare) <- udc_compare
  
  for (level in dosage_levels) {
    markera <-
      matrix_list_template[[as.character(dosage_combinations[level, 1])]]
    
    markerb <-
      matrix_list_compare[[as.character(dosage_combinations[level, 2])]]
    
    compare_vec <-
      as.integer(colSums(markera & markerb, na.rm = T))
    count_mat[, level] <- compare_vec
  }
  
  n_tot <- count_mat[,"n_00"] + count_mat[,"n_01"] + count_mat[,"n_10"] + count_mat[,"n_11"]
  r_disom <- (count_mat[,"n_00"] + count_mat[,"n_11"])/n_tot
  
  LOD <- n_tot*log10(2/n_tot) + (count_mat[,"n_01"]+count_mat[,"n_10"])*log10(count_mat[,"n_01"] + count_mat[,"n_10"]) +
    (count_mat[,"n_00"] + count_mat[,"n_11"])*log10(count_mat[,"n_00"] + count_mat[,"n_11"])
  
  p.val <- sapply(1:nrow(count_mat), function(l) binom.test((count_mat[,"n_00"][l] + count_mat[,"n_11"][l]), n_tot[l], 1/3, "less")$p.value)
  p.adj <- p.adjust(p.val, method = adj.method)
  
  ## Assume the repulsion alleles are on homologues (within subgenome)
  pref.p1 <- (2*(count_mat[,"n_01"] + count_mat[,"n_10"]) - 4*(count_mat[,"n_00"] + count_mat[,"n_11"])) / (3*n_tot)
  ## The estimator across homoeologues is slightly different (between subgenomes)
  pref.p2 <- (8*(count_mat[,"n_00"] + count_mat[,"n_11"]) - 4*(count_mat[,"n_01"] + count_mat[,"n_10"])) / (3*n_tot)
  
  ## Only take the positive values in the final column:
  pref.p <- pref.p1
  pref.p[pref.p<0] <- pref.p2[pref.p<0]
  if(!all(pref.p >= 0)) warning("Unexpected preferential pairing parameters estimated. Suggest to contact the package developer.")
  
  output <- cbind(assignment_data[combs[,1],],
                  assignment_data[combs[,2],],
                  count_mat,
                  data.frame("n" = n_tot, "r_disomic" = r_disom,"LOD" = LOD,
                             "P_value"=p.val,
                             "P_value.adj" = p.adj,"pref.p_homol" = pref.p1,"pref.p_homoeol" = pref.p2, "pref.p" = pref.p)
  )
  colnames(output)[1:8] <- c("marker_a","LG_a","Hom_a","pos_a","marker_b","LG_b","Hom_b","pos_b")
  
  return(output)
} #test_prefpairing()


#' Write MapChart file
#' @description Write a .mct file of a maplist for external plotting with MapChart software (Voorrips ).
#' @param maplist A list of maps. In the first column marker names and in the second their position. All map data are
#' compiled into a single MapChart file.
#' @param mapdir Directory to which .mct files are written, by default the same directory
#' as for \code{\link{MDSMap_from_list}}
#' @param file_info A character string added to the first lines of the .mct file, by default a datestamp is recorded.
#' @param filename Character string of filename to write the .mct file to, by default "MapFile"
#' @param precision To how many decimal places should marker positions be specified (default = 2)?
#' @param showMarkerNames Logical, by default \code{FALSE}, if \code{TRUE}, the marker names will be diplayed in the
#' MapChart output as well.
#' @examples
#' \dontrun{
#' data("integrated.maplist")
#' write.mct(integrated.maplist)}
#' @export
write.mct <- function(maplist,
                      mapdir = "mapping_files_MDSMap",
                      file_info = paste("; MapChart file created on",Sys.Date()),
                      filename = "MapFile",
                      precision = 2,
                      showMarkerNames = FALSE){
  
  if(!file.exists(mapdir)) dir.create(mapdir)
  
  if(is.null(names(maplist))){
    if(length(maplist) == 1) stop("maplist provided was an unnamed list of length 1. Chromosome number not found.\n Assign names to maplist using names() function first!")
    warning("maplist provided was an unnamed list. Proceeding with automatic names")
    names(maplist) <- paste0("LG",seq(length(maplist)))
  }
  
  write(c(file_info,"\n"), file = file.path(mapdir,paste0(filename,".mct")))
  
  affix <- ifelse(showMarkerNames,"","O")
  
  for(lg in names(maplist)){
    mapdat <- maplist[[lg]]
    
    write(c("\n",
            paste(paste("GROUP",lg),
                  paste("S=",min(mapdat$position),sep=""),
                  paste("E=",round(max(mapdat$position),precision),sep=""),
                  sep="\t")
    ), file = file.path(mapdir,paste0(filename,".mct")),
    append=TRUE)
    
    for(r in 1:nrow(mapdat)){
      write(paste(mapdat[r,"marker"],round(mapdat[r,"position"],precision),affix,sep="  "),
            file = file.path(mapdir,paste0(filename,".mct")),
            append=TRUE)
    }
    
  }
}

#' Write a JoinMap compatible .pwd file from linkage data.frame.
#' @description Output of this function allows to use JoinMap to perform the marker ordering step.
#' @param linkage_df A linkage \code{data.frame}.
#' @param pwd_file A character string specifying a file open for writing.
#' @param file_info A character string added to the first lines of the .pwd file.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @examples
#' \dontrun{
#' data("all_linkages_list_P1_split")
#' write.pwd(all_linkages_list_P1_split[["LG3"]][["homologue1"]],
#'            "LG3_homologue1_P1.pwd",
#'            "Please feed me to JoinMap")}
#' @export
write.pwd <-
  function(linkage_df,
           pwd_file,
           file_info,
           log = NULL) {
    linkage_df <- test_linkage_df(linkage_df)
    fileConn <- file(pwd_file)
    writeLines(file_info, fileConn)
    close(fileConn)
    linkage_df$phase <- NULL
    write.table(
      linkage_df,
      file = pwd_file,
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE,
      append = TRUE,
      quote = FALSE
    )
    if (!is.null(log)) {
      matc <- match.call()
      write.logheader(matc, log)
    }
  }


#' Write TetraploidSNPMap input file
#' @description Output the phased linkage map files into format readable by TetraploidSNPMap (Hackett et al. 2017) to perform QTL analysis.
#' @param phased.maplist Phased maps in list format, the output of \code{\link{create_phased_maplist}}
#' @param outputdir Directory to which TetraploidSNPMap files are written, by default written to "TetraploidSNPMap_QTLfiles" folder
#' @param filename Character string of filename stem to write the output files to, by default "TSNPM" with linkage groups names appended
#' @param ploidy The ploidy of the species, currently only 4 is supported by TetraploidSNPMap
#' @param verbose Should messages be sent to stdout?
#' @return \code{NULL}
#' @examples
#' \dontrun{
#' data("phased.maplist")
#' write.TSNPM(phased.maplist,ploidy=4)}
#' @export
write.TSNPM <- function(phased.maplist,
                        outputdir = "TetraploidSNPMap_QTLfiles",
                        filename = "TSNPM",
                        ploidy,
                        verbose=FALSE) {
  
  if(ploidy != 4) stop("Currently only tetraploid data supported by TetraploidSNPMap!")
  
  if(!dir.exists(outputdir)) dir.create(outputdir)
  
  for(lg in seq(length(phased.maplist))){
    
    nr.mapped <- nrow(phased.maplist[[lg]])
    
    fileConn <- file(file.path(outputdir,paste0(filename,"_LG",lg,".txt")))
    
    temp.phase <- phased.maplist[[lg]][,paste0("h",1:(2*ploidy))] + 1
    
    ## Very specific spacing requirements for the input, has to be decimal point of position at col 31. Use writelines function
    outlines <- rep("n",nr.mapped)
    
    for(i in 1:nr.mapped){
      spacers <- 30 - nchar(as.character(phased.maplist[[lg]][i,"marker"])) - nchar(floor(phased.maplist[[lg]][i,"position"]))
      
      if(spacers <= 0) stop(paste("Currently marker names as long as:",
                                  as.character(phased.maplist[[lg]][i,"marker"]),
                                  "are unacceptable in TetraploidSNPMap. Please shorten!"))
      
      outlines[i] <- paste0(phased.maplist[[lg]][i,"marker"],
                            paste0(rep(" ",spacers),collapse=""),
                            format(round(phased.maplist[[lg]][i,"position"],2),nsmall=2),
                            "  ",
                            paste0(temp.phase[i,1:ploidy],collapse=""),
                            " ",
                            paste0(temp.phase[i,(ploidy+1):(2*ploidy)],collapse=""),
                            collapse="")
    }
    
    if(verbose) print(outlines)
    
    writeLines(c(nr.mapped,outlines),con = fileConn)
    
    close(fileConn)
  }
  
  return(NULL)
}

#' Write out a nested list
#' @description Write a nested list into a directory structure
#' @param nested_list A nested list.
#' @param directory Character string. Directory name to which to write the structure.
#' @param save_as_object Logical. Save as R object?
#' @param object_prefix Character. Prefix of R object. Only used if \code{save_as_object = TRUE}.
#' @param extension Character. File extension. Default is ".txt".
#' @param \dots Arguments passed to \code{\link{write.table}}
#' @examples
#' \dontrun{
#' data("all_linkages_list_P1_subset")
#' write_nested_list(nested_list = all_linkages_list_P1_subset,
#'                   directory = "all_linkages_P1",
#'                   sep="\t")}
#' @export
write_nested_list <-
  function(nested_list, directory, save_as_object = FALSE, object_prefix = directory,
           extension = if(save_as_object) ".Rdata" else ".txt", ...) {
    ff = function(x) {
      if (class(x)[1] == "list" | class(x)[1] == "array")
        lapply(x, ff)
      else
        TRUE
    }
    lnames = names(unlist(lapply(nested_list, ff)))
    fnames = file.path(file.path(directory),
                       paste0(gsub(".", "/", lnames, fixed = TRUE), extension))
    dirnames = unlist(lapply(fnames, dirname))
    varnames = strsplit(lnames, split = ".", fixed = TRUE)
    suppressWarnings(lapply(dirnames, dir.create, recursive = TRUE))
    if(save_as_object){
      for(i in seq_along(varnames)){
        obj_name <- paste(c(object_prefix, varnames[[i]]), collapse = "_")
        assign(obj_name, nested_list[[varnames[[i]]]])
        save(list = obj_name, file = fnames[[i]])
      }
    } else {
      for(i in seq_along(varnames)){
        obj <- nested_list[[varnames[[i]]]]
        trywrite <- try(as.data.frame(obj), silent = TRUE)
        if(class(trywrite) != "try-error")
          write.table(obj, file = fnames[[i]], ...)
      }
    }
  }


#' Write pwd files from a nested list
#' @description A wrapper for \code{\link{write.pwd}}, which allows to write multiple pwd files with a directory structure according to the nested linkage list.
#' @param linkages_list A nested \code{list} with linkage group on the first level and homologue on the second.
#' @param target_parent A  character string specifying the name of the target parent.
#' @param binned Logical. Are the markers binned? This information is used in the pwd header.
#' @param dir A character string specifying the directory in which the files are written. Defaults to working directory.
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @examples
#' \dontrun{
#' data("all_linkages_list_P1_split")
#' write_pwd_list(all_linkages_list_P1_split, target_parent="P1", binned=FALSE)}
#' @export
write_pwd_list <-
  function(linkages_list,
           target_parent,
           binned = FALSE,
           dir = getwd(),
           log = NULL) {
    bin <- ifelse(binned, "binned", "not binned")
    
    for (LG_name in names(linkages_list)) {
      LG_dir <- file.path(dir, LG_name)
      dir.create(LG_dir, showWarnings = FALSE)
      
      for (homologue_name in names(linkages_list[[LG_name]])) {
        linkage_df <- linkages_list[[LG_name]][[homologue_name]]
        pwd_file <-
          file.path(
            LG_dir,
            paste0(
              target_parent,
              "_",
              LG_name,
              "_",
              homologue_name,
              "_",
              bin,
              ".pwd"
            )
          )
        file_info <-
          c(
            paste(
              ";",
              target_parent,
              LG_name,
              homologue_name,
              "using JOINMAP",
              paste0("(", bin, ")") ,
              sep = " "
            ),
            "; Peter Bourke and Geert van Geest, WUR Plant Breeding",
            paste(
              "name =",
              paste(target_parent, LG_name, homologue_name, bin, sep =
                      "_"),
              sep = " "
            )
          )
        
        write.pwd(linkage_df, pwd_file, file_info)
        
      }
      
    }
    if (!is.null(log)) {
      matc <- match.call()
      write.logheader(matc, log)
    }
  }
