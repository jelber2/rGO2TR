#' Parse genomes annotated by NCBI Eukaryotic Genome Annotation Pipeline.
#'
#' @return Returns all eukaryotic genomes annotated by the NCBI Eukaryotic
#' Genome Annotation Pipeline.
#' @examples
#' euks.filtered <- get.ncbi.annot.euk.genomes()
get.ncbi.annot.euk.genomes <- function() {
  # parse the webpage
  theurl <- "https://www.ncbi.nlm.nih.gov/genome/annotation_euk/all/"
  webpage <- xml2::read_html(theurl)
  
  # converts tables in html into list of data.frames
  tables <- rvest::html_table(x=webpage)
  
  # combines all tables into one table
  tables2 <- do.call(rbind, tables)

  # gets rid of columns 7-10 (they aren't needed)
  tables2[6:10] <- list(NULL)
  
  # grabs the desired ftp links
  xmlwebpage <- XML::htmlParse(webpage)
  web.links <- XML::getHTMLLinks(xmlwebpage)
  
  # gets rid of ftp link at end of page
  web.links <- web.links[web.links != "ftp://ftp.ncbi.nlm.nih.gov/"]
  ftp.links <- data.frame("Links" = (web.links[grepl("ftp://ftp.ncbi.",
                                                     web.links)]))
  
  # adds ftp links to tables, renames tables3 to euks.filtered, converts ftp
  # links from factor to character
  tables3 <- cbind (tables2, ftp.links)
  euks.filtered <- tables3
  euks.filtered <- unique(euks.filtered)
  euks.filtered$Links <- as.character(euks.filtered$Links)
  return(euks.filtered)
}