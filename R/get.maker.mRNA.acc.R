#' Gets mRNA accession identifiers from MAKER derived Granges object.
#' 
#' @param gff3.filtered An already filtered Granges object
#' @return List of mRNA accession identifiers contained in gff3.filtered
#' @examples
#' mRNA.acc <- get.maker.mRNA.acc(gff3.filtered)
#' @export
get.maker.mRNA.acc <- function(gff3.filtered) {
  # grabs only unique mRNA accession identifers from gff3.filtered
  mRNA.acc <- unique(unlist(strsplit(sub(";",
                                  "",
                                  sub(".+Parent=(\\w+_\\d+\\-R\\w)",
                                      "\\1",
                                      gff3.filtered$tags,
                                      perl=TRUE)),
                              ",")))
  return(mRNA.acc)
}