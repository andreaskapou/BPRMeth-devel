#' Pre-process BS-Seq data in any given format
#'
#' \code{preprocess_bs_seq} is a general function for reading and preprocessing
#' BS-Seq data. If a vector of files is given, these are considered as
#' replicates and are pooled together. Finally, noisy reads are discarded.
#'
#' @param files A vector of filenames containing replicate experiments. This can
#'   also be just a single replicate.
#' @param file_format A string denoting the file format that the BS-Seq data are
#'   stored. Current version allows "\code{encode_rrbs}" or "\code{bismark_cov}"
#'   formats.
#' @param chr_discarded A vector with chromosome names to be discarded.
#' @param min_bs_cov The minimum number of reads mapping to each CpG site. CpGs
#'   with less reads will be considered as noise and will be discarded.
#' @param max_bs_cov The maximum number of reads mapping to each CpG site. CpGs
#'   with more reads will be considered as noise and will be discarded.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object. The GRanges object
#'   contains two additional metadata columns: \itemize{ \item
#'   \code{total_reads}: total reads mapped to each genomic location. \item
#'   \code{meth_reads}: methylated reads mapped to each genomic location. }
#'   These columns can be accessed as follows: \code{granges_object$total_reads}
#'
#' @section Additional Info: Information about the file formats can be found in
#'   the following links:
#'
#'   Encode RRBS format:
#'   \url{http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?db=
#'   hg19&hgta_group=regulation&hgta_track=wgEncodeHaibMethylRrbs&hgta_table=
#'   wgEncodeHaibMethylRrbsBcbreast0203015BiochainSitesRep2&hgta_doSchema=
#'   describe+table+schema}
#'
#'   Bismark Cov format: \url{http://rnbeads.mpi-inf.mpg.de/data/RnBeads.pdf}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_bs_bismark_cov}}, \code{\link{read_bs_encode_haib}}
#'   \code{\link{pool_bs_seq_rep}}
#'
#' @examples
#' # Obtain the path to the files
#' bs_file <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
#' bs_data <- preprocess_bs_seq(bs_file, file_format = "encode_rrbs")
#'
#' # Extract the total reads and methylated reads
#' total_reads <- bs_data$total_reads
#' meth_reads <- bs_data$meth_reads
#'
#' @export
preprocess_bs_seq <- function(files, file_format = "encode_rrbs",
                              chr_discarded = NULL, min_bs_cov = 4,
                              max_bs_cov = 1000){
    # If we have more than one replicates
    if (length(files) > 1){
        bs_data <- pool_bs_seq_rep(files         = files,
                                   file_format   = file_format,
                                   chr_discarded = chr_discarded)
    }else{
        if (file_format == "encode_rrbs"){
            bs_data <- read_bs_encode_haib(file          = files,
                                           chr_discarded = chr_discarded,
                                           is_GRanges    = TRUE)
        }else if (file_format == "bismark_cov"){
            bs_data <- read_bs_bismark_cov(file          = files,
                                           chr_discarded = chr_discarded,
                                           is_GRanges    = TRUE)
        }
        else{
            stop("Wrong file format. Please check the available file formats!")
        }
    }

    bs_data <- .discard_bs_noise_reads(bs_data     = bs_data,
                                       min_bs_cov  = min_bs_cov,
                                       max_bs_cov  = max_bs_cov)
    message("Done!\n")
    return(bs_data)
}
