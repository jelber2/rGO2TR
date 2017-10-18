#' Downloads the desired gff3 gene annotation file.
#' 
#' @param query.results Name given to ouput of \code{search.annot.euks} function
#' @param output.file.name Desired name to give to output file
#' @return A GRanges object of the gene annotations
#' @examples
#' gff3 <- get.gff3(query.results, "sparrow.genome.gff3.gz")
#' @export
get.gff3 <- function(query.results, output.file.name) {
  # appends "GFF/" to the ftp link for query result
  ftp.search.path <- paste(query.results$Links, "GFF/", sep = "")
  if (length(ftp.search.path) > 1) stop('\n\n You have selected > 1 species!
                                        Select only 1 species in query.results!
                                        Rerun search.AnnotatedEuks function with
                                        more specific query!\n')
  
  # looks at files in GFF for query species
  ftp.search.list <- ftpList(ftp = ftp.search.path, fileonly = TRUE)
  ftp.search.num <- grep("_top_level.gff3.gz", ftp.search.list$name)
  ftp.desired.path <- ftp.search.list$name[ftp.search.num]

  # if there is more than 1 assembly, let the user choose the desired assembly
  if (length(ftp.desired.path) > 1)
     cat('\nThere is more than 1 assembly. Which one do you want?\n')
     for (i in 1:length(ftp.desired.path)) cat('\n',ftp.desired.path[i])
     cat('\n')
     cat('\n')
     file.number <- as.integer(readline(prompt='Enter the reference number (ex: 1 for the first, 2 for the second, etc): '))
     ftp.desired.path <- ftp.desired.path[file.number]

  # creates the path to download the GFF3 file from queried species
  gff3.url <- paste(ftp.search.path, ftp.desired.path, sep = "")
  
  # downloads the desired file to the current working directory
  if (file.exists(output.file.name)) 
    stop('\n\n Choose a different file name. You already have a downloaded file
         with the same name!\n')
  else 
    gff3.file.gz <- download.file(gff3.url, output.file.name)
  # read the gff3 file
  gff3 <- read.gff(output.file.name, locus.tags = FALSE)
  return(gff3)
}
