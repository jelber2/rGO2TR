#' Filters a GRange object by desired annotation source.
#' 
#' @param file.name Name of GRange object to filter
#' @param type.to.filter Valid type such as CDS, gene, or exon to filter by
#' @return Returns a filtered gff3 file
#' @examples
#' gff3.filtered <- filter.gff3(gff3, 'exon')
#' @export
filter.gff3 <- function(file.name, type.to.filter) {
  # filters gff3 annotation file by desired feature
  file.name.filtered <- file.name[file.name$feature == type.to.filter, ]
  
  # filters gff3 annotation file by predicted mRNA
  file.name.filtered <- file.name.filtered[grepl("XM_", file.name.filtered$tags), ]
  return(file.name.filtered)
}