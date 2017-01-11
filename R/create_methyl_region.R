#' Create methylation regions for each gene promoter.
#'
#' \code{create_methyl_region} creates methylation regions using BS-Seq and
#' annotated gene promoter regions. BS-Seq data give information for the
#' methylation of CpGs individually, and annotated data are used to locate the
#' TSS of each gene and its promoter region.
#'
#' @param bs_data \code{\link[GenomicRanges]{GRanges}} object containing the
#'   BS-Seq data. The GRanges object should also have two additional metadata
#'   columns named \code{total_reads} and \code{meth_reads}. A GRanges object
#'   used in this function can be the output of
#'   \code{\link{read_bs_encode_haib}} or its wrapper function
#'   \code{\link{preprocess_bs_seq}}.
#' @param prom_region \code{\link[GenomicRanges]{GRanges}} object containing
#'   promoter regions, i.e. N bp upstream and M bp downstream of TSS location.
#'   The GRanges object should also have one additional metadata column named
#'   \code{tss}. A GRanges object used in this function can be the output of
#'   \code{\link{create_prom_region}}.
#' @param cpg_density Optional integer defining the minimum number of CpGs that
#'   have to be in a methylated region. Regions with less than \code{n} CpGs are
#'   discarded.
#' @param sd_thresh Optional numeric defining the minimum standard deviation of
#'   the methylation change in a region. This is used to filter regions with no
#'   methylation change.
#' @param ignore_strand Logical, whether or not to ignore strand information.
#' @param is_single_cell Logical, if we work with single cell data, then we keep
#'   only the information of methylated or not for each CpG location.
#' @param fmin Optional minimum range value for region location scaling. Under
#'   this version, this parameter should be left to its default value.
#' @param fmax Optional maximum range value for region location scaling. Under
#'   this version, this parameter should be left to its default value.
#'
#' @return A \code{methyl_region} object containing the following information:
#'   \itemize{ \item{ \code{meth_data}: A list containing methylation data,
#'   where each entry in the list is an \eqn{L_{i} X 3} dimensional matrix,
#'   where \eqn{L_{i}} denotes the number of CpGs found in region \code{i}. The
#'   columns contain the following information: \enumerate{ \item{ 1st column:
#'   Contains the locations of CpGs relative to TSS. Note that the actual
#'   locations are scaled to the (fmin, fmax) region. } \item{ 2nd column:
#'   Contains the total reads of each CpG in the corresponding location. NOTE
#'   that for single-cell data this column is not stored in the object.} \item{
#'   3rd column: Contains the methylated reads each CpG in the corresponding
#'   location.} } } \item{ \code{prom_ind}: A vector storing the corresponding
#'   promoter indices, so as to map each methylation region with its
#'   corresponding gene promoter.} } The lengths of \code{meth_data} and
#'   \code{prom_ind} should be the same.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{preprocess_bs_seq}}, \code{\link{create_prom_region}}
#'
#' @examples
#' # Obtain the path to the BS file and then read it
#' bs_file <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
#' bs_data <- read_bs_encode_haib(bs_file)
#'
#' # Create promoter regions
#' rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "BPRMeth")
#' annot_data <- read_rna_encode_caltech(rnaseq_file)
#' prom_region <- create_prom_region(annot_data)
#'
#' # Finally, create methylation regions
#' meth_region <- create_methyl_region(bs_data, prom_region)
#'
#' @importFrom methods is
#' @export
create_methyl_region <- function(bs_data, prom_region, cpg_density = 10,
                                 sd_thresh = 10e-02, ignore_strand = TRUE,
                                 is_single_cell = FALSE, fmin = -1, fmax = 1){

    message("Creating methylation regions ...")
    assertthat::assert_that(methods::is(bs_data, "GRanges"))
    assertthat::assert_that(methods::is(prom_region, "GRanges"))

    # Find overlaps between promoter regions and BS-Seq data -------
    overlaps <- GenomicRanges::findOverlaps(query   = prom_region,
                                            subject = bs_data,
                                            ignore.strand = ignore_strand)
    if (length(overlaps) < 2){
        stop("Not enough matches between the BS-Seq and RNA-Seq data.")
    }

    # Convert data in vector format for faster lookup --------------
    query_hits <- S4Vectors::queryHits(overlaps)
    subj_hits  <- S4Vectors::subjectHits(overlaps)

    prom_loc   <- unique(query_hits)     # Indices of promoter locations
    tss_loc    <- prom_region$tss        # TSS locations
    prom_id    <- prom_region$id         # (Ensembl) IDs
    tss_strand <- as.character(GenomicRanges::strand(prom_region))
    cpg_loc    <- GenomicRanges::ranges(bs_data)@start  # CpG locations
    meth_reads <- bs_data$meth_reads     # Methylated read

    # If we work with single cell data
    if (is_single_cell){
      D <- 2  # Number of columns for each matrix in the list
    } else{
      tot_reads <- bs_data$total_reads   # Total reads
      D <- 3
    }

    # Extract upstream and downstream lengths in bps
    width        <- GenomicRanges::ranges(prom_region)@width[1]
    if (identical(tss_strand[1], "+")){
        upstream   <- GenomicRanges::ranges(prom_region)@start[1] - tss_loc[1]
        downstream <- width + upstream - 1
    }else if (identical(tss_strand[1], "-")){
        downstream <- abs(GenomicRanges::ranges(prom_region)@start[1] -
                              tss_loc[1])
        upstream   <- downstream - width + 1
    }

    # Initialize variables -----------------------------------------
    LABEL        <- FALSE                     # Flag variable
    meth_data    <- list()                    # List where data will be stored
    prom_counter <- 0                         # Promoter counter
    cpg_ind      <- vector(mode = "integer")  # Vector of CpG indices
    cpg_ind      <- c(cpg_ind, subj_hits[1])  # Add the first subject hit

    for (i in 2:NROW(query_hits)){
        # If query hits is the same as the previous one
        if (query_hits[i] == query_hits[i - 1]){
            cpg_ind <- c(cpg_ind, subj_hits[i])  # Add subject hit
            # In case we have the last region
            if (i == NROW(query_hits)){
              prom_counter <- prom_counter + 1  # Increase promoter counter
              LABEL <- TRUE
            }
        }else{
            prom_counter <- prom_counter + 1  # Increase promoter counter
            LABEL <- TRUE
        }

        if (LABEL){
            # TSS location for promoter 'promCount'
            id <- prom_id[prom_loc[prom_counter]]
            meth_data[[id]] <- NA
            # Only keep regions that have more than 'n' CpGs
            if (length(cpg_ind) > cpg_density){
                # If sd of the methylation level is above threshold
                if(is_single_cell){
                    obs_var <- stats::sd(meth_reads[cpg_ind])
                }else{
                    obs_var <- stats::sd(meth_reads[cpg_ind] /
                                         tot_reads[cpg_ind])
                }
                if (obs_var > sd_thresh){
                    # Locations of CpGs in the genome
                    region <- cpg_loc[cpg_ind]
                    # TSS location for promoter 'promCount'
                    tss <- tss_loc[prom_loc[prom_counter]]
                    # Extract strand information, i.e. direction
                    strand_direction <- tss_strand[prom_loc[prom_counter]]
                    # Shift CpG locations relative to TSS
                    center_data  <- .center_loc(region = region,
                                        tss = tss,
                                        strand_direction = strand_direction)
                    # In the "-" strand the order of the locations should change
                    Order <- base::order(center_data)

                    meth_data[[id]] <- matrix(0, nrow = length(cpg_ind),
                                              ncol = D)

                    # Store normalized locations of methylated CpGs
                    meth_data[[id]][, 1] <- BPRMeth:::.minmax_scaling(
                                                data = center_data[Order],
                                                xmin = upstream,
                                                xmax = downstream,
                                                fmin = fmin,
                                                fmax = fmax)

                    if (is_single_cell){
                      # Store methylated reads in the corresponding locations
                      meth_data[[id]][, 2] <- meth_reads[cpg_ind][Order]
                    }else{
                      # Store total reads in the corresponding locations
                      meth_data[[id]][, 2] <- tot_reads[cpg_ind][Order]
                      # Store methylated reads in the corresponding locations
                      meth_data[[id]][, 3] <- meth_reads[cpg_ind][Order]
                    }
                }
            }
            LABEL   <- FALSE
            cpg_ind <- vector(mode = "integer")
            cpg_ind <- c(cpg_ind, subj_hits[i])
        }
    }
    # Obtain promoters that did not have CpG coverage
    ind <- which(is.na(match(seq_along(tss_loc), prom_loc)))
    id <- prom_id[ind]
    for (i in seq_along(ind)){
        meth_data[[id[i]]] <- NA
    }
    # Get the indices of matches between annot data and meth_data
    ind <- match(prom_id, names(meth_data))
    # Reorder methylation data to match annot_data
    meth_data <- meth_data[ind]

    message("Done!\n")
    return(meth_data)
}


# Center CpG locations relative to TSS
#
# \code{center_loc} centera CpG locations relative to TSS
#
# @param region CpG locations
# @param tss TSS location
# @param strand_direction Strand direction
#
# @return Centered location data relative to TSS
#
.center_loc <- function(region, tss, strand_direction){
    assertthat::assert_that(is.character(strand_direction))
    center <- region - tss
    if (identical(strand_direction, "-")){
        center  <- (-center)  # If '-' strand, swap CpG locations
    }
    return(center)
}
