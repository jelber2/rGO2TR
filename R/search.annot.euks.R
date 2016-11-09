#' Searches Genomes Annotated by NCBI Eukaryotic Genome Annotation Pipeline.
#' 
#' @param query A species' latin name to search for
#' @return Returns a query result based on \code{query}
#' @examples
#' query.results <- search.annot.euks("Zonotrichia albicollis")
#' query.results <- search.annot.euks("Canis lupus")
search.annot.euks <- function(query) {
  query.results <- euks.filtered[agrepl(query, euks.filtered$Species),]
  query.results$Links <- as.character(query.results$Links)
  return(query.results)
}