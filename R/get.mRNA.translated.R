#' Get translated sequences in FASTA format from mRNA accession identifier input.
#' 
#' @param input.R.object Name of R object containing mRNA accessions
#' @param output.R.object Name of R object to contain mRNA translated sequences
#' @param email.address Your email address
#' @param output.file.name Desired name of the output fasta file
#' @return An R object containing mRNA translated sequences
#' @examples
#' mRNA1.translated <- getmRNAtranslated(mRNA1.acc,
#'                                       mRNA1.translated,
#'                                      "jelber2@@lsu.edu",
#'                                      "sparrow.mRNA1.fasta")
#' mRNA2.translated <- getmRNAtranslated(mRNA2.acc,
#'                                       mRNA2.translated,
#'                                      "jelber2@@lsu.edu",
#'                                      "sparrow.mRNA2.fasta")

get.mRNA.translated <- function(input.R.object,
                              output.R.object,
                              email.address,
                              output.file.name) {
  # makes an empty list
  output.R.object <- ""
  # specifies email address
  Sys.setenv(email=email.address)
  
  # if >200 gi's, split the prot gi's into a data.frame with each column
  # having <200 accessions
  input.R.object.split <- split(input.R.object, ceiling(seq_along(input.R.object)/200))
  for (i in names(input.R.object.split))
  {
    # for each column of 200 accessions, fetch the fasta sequence using the
    # genomes::efetch function
    output.R.object <- c(output.R.object, genomes::efetch(id = input.R.object.split[[i]],
                                                          "nucleotide", "gb", "xml"))
    output.R.object <- output.R.object[grepl('[A-Z]{13,}|accession-version|>XM_\\w+.\\d translated',
                                             output.R.object,
                                             perl=TRUE)]
    output.R.object <- sub("\\s+<GBSeq_accession-version>(.+)</GBSeq_accession-version>",
                           ">\\1 translated",
                           output.R.object,
                           perl=TRUE)
    output.R.object <- sub("\\s+<GBQualifier_value>(\\w+)</GBQualifier_value>",
                           "\\1",
                           output.R.object,
                           perl=TRUE)
    
    # Pause the system for 0.25 seconds, so there aren't too many efetch calls
    # at any one time
    Sys.sleep(0.25)
    
    # This takes a while, so let the user know the progress
    cat("List", i, "of", length(names(input.R.object.split)), "finished", "\n")
  }
  
  # create an external output file
  writeLines(output.R.object, con = output.file.name, sep = "\n")
  return(output.R.object)
}