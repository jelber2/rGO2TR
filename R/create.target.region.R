#' Filters Granges object by retained mRNA list to create target region
#' 
#' @param retained.mRNA.list mRNA retained by \code{filter.mRNA.GO.list} function
#' @param gff3.filtered Granges object to filter by \code{retained.mRNA.list}
#' @return target.region in the form of a Granges object with overlapping intervals
#' @examples 
#' target.region <- create.target.region(retained.mRNA.list, gff3.filtered)
create.target.region <- function (retained.mRNA.list, gff3.filtered) {
  j=0
  target.region <- ""

  gff3.filtered$tags <- sub(".+Genbank:(XM_\\d+\\.\\d).+",
                            "\\1",
                            gff3.filtered$tags,
                            perl=TRUE)

  target.region <- gff3.filtered[0]

  for (i in retained.mRNA.list)
  {
    j = j+1
    target.region <- c(target.region,
                       gff3.filtered[grepl(i,
                                           gff3.filtered$tags),],
                       ignore.mcols=TRUE)
    cat("Finished filtering retained mRNA",
        j,
        "of",
        length(retained.mRNA.list),
        "finished",
        "\n")
  }
  return(target.region)
}