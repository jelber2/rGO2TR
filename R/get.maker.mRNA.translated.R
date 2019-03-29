#' Get translated sequences in FASTA format from MAKER mRNA accession identifier input.
#' 
#' @param input.R.object Name of R object containing mRNA accessions
#' @param output.R.object Name of R object to contain mRNA translated sequences
#' @param input.file.name Name of file with translated proteins sequences
#' @param output.file.name Desired name of the output fasta file
#' @return An R object containing mRNA translated sequences
#' @examples
#' mRNA1.translated <- get.maker.mRNA.translated(mRNA1.acc,
#'                                               mRNA1.translated,
#'                                               "maker.translated.fasta",
#'                                               "sparrow.mRNA1.translated.fasta")
#' mRNA2.translated <- get.maker.mRNA.translated(mRNA2.acc,
#'                                               mRNA2.translated,
#'                                               "maker.translated.fasta",
#'                                               "sparrow.mRNA2.translated.fasta")
#' @export
get.maker.mRNA.translated <- function(input.R.object,
                                      output.R.object,
                                      input.file.name,
                                      output.file.name) {
  # read in input.file.name
  intermediate.file <- seqinr::read.fasta(file=input.file.name,seqtype="AA",as.string=T, set.attributes=F)
  output.R.object <- intermediate.file[c(which(names(intermediate.file) %in% input.R.object))]
  seqinr::write.fasta(sequences=output.R.object, file.out=output.file.name, names=names(output.R.object))
  return(output.R.object)
}