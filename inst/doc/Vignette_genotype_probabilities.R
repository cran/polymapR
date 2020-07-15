## ---- eval = FALSE------------------------------------------------------------
#  install.packages("polymapR")
#  library(polymapR)

## ---- eval = FALSE------------------------------------------------------------
#  geno <- read.csv("fitPoly_4x_output_2343_SNPs.csv")

## ---- echo = FALSE------------------------------------------------------------
library(polymapR)
load("genoprobsdata.RData")

## ---- eval = FALSE------------------------------------------------------------
#  knitr::kable(head(geno))

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(head(geno.sub))

## ---- eval = FALSE------------------------------------------------------------
#  parent1.reps <- c("P1","P1a")
#  parent2.reps <- "P2"
#  individuals <- setdiff(unique(geno$SampleName),c(parent1.reps,parent2.reps))

## ---- eval = FALSE------------------------------------------------------------
#  chk1 <- checkF1(input_type = "probabilistic",
#                  probgeno_df = geno,
#                  parent1 = parent1.reps,
#                  parent2 = parent2.reps,
#                  F1 = individuals,
#                  polysomic = TRUE,
#                  disomic = FALSE,
#                  mixed = FALSE,
#                  ploidy = 4)

## ---- eval = FALSE------------------------------------------------------------
#  chk1$checked_F1[1:5,]

## ---- echo = FALSE------------------------------------------------------------
chk1.sub[1:5,]

## ---- eval = FALSE------------------------------------------------------------
#  pardose <- polymapR:::assign_parental_dosage(chk = chk1,
#                                               probgeno_df = geno)
#  knitr::kable(head(pardose))

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(head(pardose.sub))

## ---- eval = FALSE------------------------------------------------------------
#  length(which(chk1$checked_F1$qall_mult == 0)) #859 markers

## ---- eval = FALSE------------------------------------------------------------
#  remove.mark <- chk1$checked_F1[chk1$checked_F1$qall_mult==0,"MarkerName"]
#  geno1 <- geno[!geno$MarkerName %in% remove.mark,]

## ---- eval = FALSE------------------------------------------------------------
#  gpo <- gp_overview(probgeno_df = geno1)

## ---- out.width = "500px", echo = FALSE, fig.align="center"-------------------
knitr::include_graphics("figures/gp_overview.png")

## ---- eval = FALSE------------------------------------------------------------
#  geno1 <- gpo$probgeno_df

## ---- eval = FALSE------------------------------------------------------------
#  maxP.chk <- check_maxP(probgeno_df = geno1)

## ---- echo = FALSE------------------------------------------------------------
print(maxP.chk)

## ---- out.width = "550px", echo = FALSE, fig.align="center"-------------------
knitr::include_graphics("figures/maxPdist.png")

## ---- eval = FALSE------------------------------------------------------------
#  geno2 <- screen_for_duplicate_individuals.gp(probgeno_df = geno2,
#                                               ploidy = 4,
#                                               parent1 = parent1.reps,
#                                               parent2 = parent2.reps,
#                                               F1 = individuals)

## ---- out.width = "550px", echo = FALSE, fig.align="center"-------------------
knitr::include_graphics("figures/dup_indivs.png")

## ---- echo = FALSE------------------------------------------------------------
 write("\nNo duplicates found\n",stdout())

## -----------------------------------------------------------------------------
nc <- parallel::detectCores() - 2

## ---- eval = FALSE------------------------------------------------------------
#  chk1 <- checkF1(input_type = "probabilistic",
#                  probgeno_df = geno2,
#                  parent1 = parent1.reps,
#                  parent2 = parent2.reps,
#                  F1 = individuals,
#                  polysomic = TRUE,
#                  disomic = FALSE,
#                  mixed = FALSE,
#                  ploidy = 4)

## ---- eval = FALSE------------------------------------------------------------
#  SN_SN_P1 <- linkage.gp(probgeno_df = geno2,
#                         chk = chk1,
#                         markertype1 = c(1,0),
#                         target_parent = "P1",
#                         LOD_threshold = 3,
#                         ncores = nc)

## ---- eval = FALSE------------------------------------------------------------
#  head(SN_SN_P1)

## ---- echo = FALSE------------------------------------------------------------
head(SN_SN_P1.sub)

## ---- eval = FALSE------------------------------------------------------------
#  par(mfrow = c(2,1))
#  P1_homologues <- cluster_SN_markers(linkage_df = SN_SN_P1,
#                                      LOD_sequence = seq(3,12,0.5),
#                                      LG_number = 5,
#                                      ploidy = 4,
#                                      parentname = "P1",
#                                      plot_clust_size = FALSE,
#                                      min_clust_size = 3)

## ---- out.width = "600px", echo = FALSE, fig.align="center"-------------------
knitr::include_graphics("figures/p1_clustering_gp.png")

## -----------------------------------------------------------------------------
sort(table(P1_homologues[['7']]$cluster),decreasing = T)
length(table(P1_homologues[['7']]$cluster)) 

## ---- eval = FALSE------------------------------------------------------------
#  SN_SN_P2 <- linkage.gp(probgeno_df = geno2,
#                         chk = chk1,
#                         markertype1 = c(1,0),
#                         target_parent = "P2",
#                         LOD_threshold = 3,
#                         ncores = nc)

## ---- eval = FALSE------------------------------------------------------------
#  P2_homologues <- cluster_SN_markers(linkage_df = SN_SN_P2,
#                                      LOD_sequence = seq(3,12,0.5),
#                                      LG_number = 5,
#                                      ploidy = 4,
#                                      parentname = "P2",
#                                      plot_clust_size = F,
#                                      min_clust_size = 3)

## ---- out.width = "600px", echo = FALSE, fig.align="center"-------------------
knitr::include_graphics("figures/p2_clustering_gp.png")

## -----------------------------------------------------------------------------
sort(table(P2_homologues[['6']]$cluster),decreasing = T)
length(table(P2_homologues[['6']]$cluster)) 

## ---- eval = FALSE------------------------------------------------------------
#  SN_SS_P1 <- linkage.gp(probgeno_df = geno2,
#                         chk = chk1,
#                         markertype1 = c(1,0),
#                         markertype2 = c(1,1),
#                         target_parent = "P1",
#                         ncores = nc)
#  
#  SN_SS_P2 <- linkage.gp(probgeno_df = geno2,
#                         chk = chk1,
#                         markertype1 = c(1,0),
#                         markertype2 = c(1,1),
#                         target_parent = "P2",
#                         ncores = nc)

## ---- eval = FALSE------------------------------------------------------------
#  LGHomDf_P1 <- bridgeHomologues(cluster_stack = P1_homologues[["7"]],
#                                 cluster_stack2 = P2_homologues[["7"]],
#                                 linkage_df = SN_SS_P1,
#                                 linkage_df2 = SN_SS_P2,
#                                 LOD_threshold = 5,
#                                 LG_number = 5)

## ---- out.width = "450px", echo = FALSE, fig.align="center"-------------------
knitr::include_graphics("figures/P1_bridges.png")

## -----------------------------------------------------------------------------
table(LGHomDf_P1$LG,LGHomDf_P1$homologue)

## ---- eval = FALSE------------------------------------------------------------
#  LGHomDf_P2 <- bridgeHomologues(cluster_stack = P2_homologues[["6"]],
#                                 cluster_stack2 = P1_homologues[["6"]],
#                                 linkage_df = SN_SS_P2,
#                                 linkage_df2 = SN_SS_P1,
#                                 LOD_threshold = 5,
#                                 LG_number = 5)

## ---- out.width = "450px", echo = FALSE, fig.align="center"-------------------
knitr::include_graphics("figures/P2_bridges.png")

## -----------------------------------------------------------------------------
table(LGHomDf_P2$LG,LGHomDf_P2$homologue)

## ---- eval = FALSE------------------------------------------------------------
#  LGHomDf_P1a <-cluster_per_LG(LG = 3,
#                               linkage_df = SN_SN_P1[SN_SN_P1$phase == "coupling",],
#                               LG_hom_stack = LGHomDf_P1,
#                               LOD_sequence = 3:10,
#                               modify_LG_hom_stack = TRUE,
#                               network.layout = "stacked",
#                               nclust_out = 4,
#                               label.offset=1.2)

## ---- out.width = "450px", echo = FALSE, fig.align="center"-------------------
knitr::include_graphics("figures/cluster_per_LG_P1.3.png")

## -----------------------------------------------------------------------------
table(LGHomDf_P1a$LG,LGHomDf_P1a$homologue)

## -----------------------------------------------------------------------------
head(LGHomDf_P1a)

## ---- eval = FALSE------------------------------------------------------------
#  P1_SxS_Assigned <- assign_linkage_group(linkage_df = SN_SS_P1,
#                                          LG_hom_stack = LGHomDf_P1a,
#                                          SN_colname = "marker_a",
#                                          unassigned_marker_name = "marker_b",
#                                          phase_considered = "coupling",
#                                          LG_number = 5,
#                                          LOD_threshold = 3,
#                                          ploidy = 4)
#  
#  P2_SxS_Assigned <- assign_linkage_group(linkage_df = SN_SS_P2,
#                                          LG_hom_stack = LGHomDf_P2,
#                                          SN_colname = "marker_a",
#                                          unassigned_marker_name = "marker_b",
#                                          phase_considered = "coupling",
#                                          LG_number = 5,
#                                          LOD_threshold = 3,
#                                          ploidy = 4)

## -----------------------------------------------------------------------------
LGHomDf_P2c <- consensus_LG_names(modify_LG = LGHomDf_P2, 
                                 template_SxS = P1_SxS_Assigned, 
                                 modify_SxS = P2_SxS_Assigned)

## ---- eval = FALSE------------------------------------------------------------
#  save(LGHomDf_P1a,LGHomDf_P2c, file = "LGHomDf_stacks.Rdata") #for example..

## ---- eval = FALSE------------------------------------------------------------
#  P2_SxS_Assigned <- assign_linkage_group(linkage_df = SN_SS_P2,
#                                          LG_hom_stack = LGHomDf_P2c, #this is changed
#                                          SN_colname = "marker_a",
#                                          unassigned_marker_name = "marker_b",
#                                          phase_considered = "coupling",
#                                          LG_number = 5,
#                                          LOD_threshold = 3,
#                                          ploidy = 4)

## ---- eval = FALSE------------------------------------------------------------
#  marker_assignments_P1 <- homologue_lg_assignment(input_type = "probabilistic",
#                                                   probgeno_df = geno2,
#                                                   chk = chk1,
#                                                   assigned_list = list(P1_SxS_Assigned),
#                                                   assigned_markertypes = list(c(1,1)),
#                                                   LG_hom_stack = LGHomDf_P1a,
#                                                   target_parent = "P1",
#                                                   other_parent = "P2",
#                                                   ploidy = 4,
#                                                   pairing = "random",
#                                                   convert_palindrome_markers = FALSE,
#                                                   LG_number = 5,
#                                                   LOD_threshold = 3,
#                                                   write_intermediate_files = FALSE)
#  
#  marker_assignments_P2 <- homologue_lg_assignment(input_type = "probabilistic",
#                                                   probgeno_df = geno2,
#                                                   chk = chk1,
#                                                   assigned_list = list(P2_SxS_Assigned),
#                                                   assigned_markertypes = list(c(1,1)),
#                                                   LG_hom_stack = LGHomDf_P2c,
#                                                   target_parent = "P2",
#                                                   other_parent = "P1",
#                                                   ploidy = 4,
#                                                   pairing = "random",
#                                                   convert_palindrome_markers = FALSE,
#                                                   LG_number = 5,
#                                                   LOD_threshold = 3,
#                                                   write_intermediate_files = FALSE)

## ---- eval = FALSE------------------------------------------------------------
#  marker_assignments <- check_marker_assignment(marker_assignments_P1,marker_assignments_P2)

## ---- eval = FALSE------------------------------------------------------------
#  saveRDS(marker_assignments, file = "marker_assignments.RDS")

## ---- eval = FALSE------------------------------------------------------------
#  all_linkages_list_P1 <- finish_linkage_analysis(input_type = "probabilistic",
#                                                  marker_assignment = marker_assignments$P1,
#                                                  probgeno_df = geno2,
#                                                  chk = chk1,
#                                                  target_parent = "P1",
#                                                  other_parent = "P2",
#                                                  convert_palindrome_markers = FALSE,
#                                                  ploidy = 4,
#                                                  pairing = "random",
#                                                  LG_number = 5,
#                                                  ncores = nc)
#  
#  all_linkages_list_P2 <- finish_linkage_analysis(input_type = "probabilistic",
#                                                  marker_assignment = marker_assignments$P2,
#                                                  probgeno_df = geno2,
#                                                  chk = chk1,
#                                                  target_parent = "P2",
#                                                  other_parent = "P1",
#                                                  convert_palindrome_markers = TRUE,
#                                                  ploidy = 4,
#                                                  pairing = "random",
#                                                  LG_number = 5,
#                                                  ncores = nc)

## ---- eval = FALSE------------------------------------------------------------
#  linkages <- list()
#  for(lg in names(all_linkages_list_P1)){
#    linkages[[lg]] <- rbind(all_linkages_list_P1[[lg]], all_linkages_list_P2[[lg]])
#  }

## ---- eval = FALSE------------------------------------------------------------
#  saveRDS(linkages, file = "linkages.RDS")

## ---- eval = FALSE------------------------------------------------------------
#  linkages <- readRDS("linkages.RDS")

## ---- eval = FALSE------------------------------------------------------------
#  integrated.maplist <- MDSMap_from_list(linkages)

## ---- eval = FALSE------------------------------------------------------------
#  phased.maplist <- create_phased_maplist(input_type = "probabilistic",
#                                          maplist = integrated.maplist,
#                                          chk = chk1,
#                                          ploidy = 4,
#                                          marker_assignment.1 = marker_assignments$P1,
#                                          marker_assignment.2 = marker_assignments$P2)

