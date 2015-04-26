#' Gets mRNA accession identifiers from desired Granges object.
#' 
#' @param gff3.filtered An already filtered Granges object
#' @return List of mRNA accession identifiers contained in gff3.filtered
#' @examples
#' mRNA.acc <- get.mRNA.acc(gff3.filtered)
get.mRNA.acc <- function(gff3.filtered) {
  # grabs only unique mRNA accession identifers from gff3.filtered
  mRNA.acc <- unique(sub(".+Genbank:(XM_\\d+\\.\\d).+",
                         "\\1",
                         gff3.filtered$tags,
                         perl=TRUE))
  return(mRNA.acc)
}