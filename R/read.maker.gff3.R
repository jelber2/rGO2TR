#' Reads a MAKER create gff3 file into R.
#' 
#' @param gff3.file MAKER created GFF3 file (can be gzipped) that you wish to read into R
#' @return A GRanges object of the gene annotations
#' @examples
#' gff3 <- read.maker.gff3(gff3.file)
#' @export
read.maker.gff3 <- function(gff3.file) {
  gff3 <- read.gff(gff3.file, locus.tags = FALSE)
  return(gff3)
}