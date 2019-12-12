#' A dosage matrix for a random pairing tetraploid with five linkage groups.
#' @format A matrix
"ALL_dosages"

#' @rdname ALL_dosages
"segregating_data"

#' @rdname ALL_dosages
"screened_data"

#' @rdname ALL_dosages
"screened_data2"

#' @rdname ALL_dosages
"screened_data3"

#' @rdname ALL_dosages
"TRI_dosages"

#' A linkage \code{data.frame}.
#' @format 
#' \itemize{
#'   \item marker_a. First marker in comparison
#'   \item marker_b. Second marker in comparison
#'   \item r. recombination frequency
#'   \item LOD. LOD score
#'   \item phase. The phase between markers
#' }
"SN_SN_P1"

#' @rdname SN_SN_P1
"SN_SN_P2"

#' @rdname SN_SN_P1
"SN_SS_P1"

#' @rdname SN_SN_P1
"SN_SS_P2"

#' @rdname SN_SN_P1
"SN_DN_P1"

#' @rdname SN_SN_P1
"SN_DN_P2"

#' @rdname SN_SN_P1
"SN_SN_P2_triploid"

#' A list of cluster stacks at different LOD scores
#' @format A list with with LOD thresholds as names. The list contains dataframes with the following format:
#' \itemize{
#'  \item marker. markername
#'  \item pseudohomologue. name of (pseudo)homologue
#'  }
"P1_homologues"

#' @rdname P1_homologues
"P2_homologues"

#' @rdname P1_homologues
"P2_homologues_triploid"

#' A \code{data.frame} with marker assignments
#' @format A data.frame with at least the following columns:
#' \itemize{
#'  \item Assigned_LG. The assigned linkage group
#'  \item Assigend_hom1. The homologue with most linkages
#' }
#' The columns LG1 - LGn and Hom1 - Homn give the number of hits per marker for that linkage group/homologue. Assigned_hom2 .. gives the nth homologue with most linkages. 
"P1_SxS_Assigned"

#' @rdname P1_SxS_Assigned
"P2_SxS_Assigned"

#' @rdname P1_SxS_Assigned
"P2_SxS_Assigned_2"

#' @rdname P1_SxS_Assigned
"P1_DxN_Assigned"

#' @rdname P1_SxS_Assigned
"P2_DxN_Assigned"

#' @rdname P1_SxS_Assigned
"marker_assignments_P1"

#' @rdname P1_SxS_Assigned
"marker_assignments_P2"

#' A \code{data.frame} specifying the assigned homologue and linkage group number per SxN marker
#' @format 
#' \itemize{
#'  \item SxN_Marker. Markername of simplex nulliplex marker
#'  \item homologue. Assigned homologue number
#'  \item LG Assigned. linkage group number
#' }
"LGHomDf_P1_1"

#' @rdname LGHomDf_P1_1
"LGHomDf_P2_1"

#' @rdname LGHomDf_P1_1
"LGHomDf_P2_2"

#' A (nested) list of linkage data frames classified per linkage group and homologue
"all_linkages_list_P1"

#' @rdname all_linkages_list_P1
"all_linkages_list_P1_split"

#' @rdname all_linkages_list_P1
"all_linkages_list_P1_subset"

#' A nested list with integrated maps
"integrated.maplist"

#' A list of phased maps
"phased.maplist"

#' A list of maps of one parent
"maplist_P1"

#' @rdname maplist_P1
"maplist_P1_subset"

#' @rdname maplist_P1
"maplist_P2_subset"

#' A sample map
"map1"

#' A sample map
"map2"

#' A sample map
"map3"
