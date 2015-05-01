#' Parse genomes annotated by NCBI Eukaryotic Genome Annotation Pipeline.
#'
#' @return Returns all eukaryotic genomes annotated by the NCBI Eukaryotic
#' Genome Annotation Pipeline.
#' @examples
#' euks.filtered <- get.ncbi.annot.euk.genomes()
get.ncbi.annot.euk.genomes <- function() {
  # parse the webpage
  theurl <- "http://www.ncbi.nlm.nih.gov/genome/annotation_euk/all/"
  webpage <- rvest::html(theurl)

  # grabs the name for the first table called "Featured"
  nn <- sub("(.+) \\(\\d+\\)",
            "\\1",
            XML::xpathSApply(webpage,
                        '//*[@class="jig-ncbitoggler-open ui-ncbitoggler-
                        group-mygroup1"]',
                        xmlValue),
            perl=TRUE)
  # grabs the rest of the table names
  nn <- c(nn, sub("(.+) \\(\\d+\\)",
                  "\\1",
                  XML::xpathSApply(webpage,
                              '//*[@class="jig-ncbitoggler ui-ncbitoggler-
                              group-mygroup1"]',
                              xmlValue),
                  perl=TRUE))
  
  # converts tables in html into list of data.frames
  tables <- rvest::html_table(x=webpage)
  
  # applies table names to each table
  names(tables) = nn
  
  # adds the group name to each table
  tables <- lapply(seq(length(tables)), function(i) {
    cbind(tables[[i]], Group = names(tables)[i])
  })
  
  # combines all tables into one table
  tables2 <- do.call(rbind, tables)
  
  # gets rid of columns 7-9 (they aren't needed)
  tables2[6:9] <- list(NULL)
  
  # grabs the desired ftp links
  web.links <- XML::getHTMLLinks(webpage)
  
  # gets rid of ftp link at end of page
  web.links <- web.links[web.links != "ftp://ftp.ncbi.nlm.nih.gov/"]
  ftp.links <- data.frame("Links" = (web.links[grepl("ftp://ftp.ncbi.",
                                                     web.links)]))
  
  # adds ftp links to tables, renames tables3 to euks.filtered, converts ftp
  # links from factor to character
  tables3 <- cbind (tables2, ftp.links)
  euks.filtered <- tables3
  euks.filtered$Links <- as.character(euks.filtered$Links)
  euks.filtered <- unique(euks.filtered)
  return(euks.filtered)
}