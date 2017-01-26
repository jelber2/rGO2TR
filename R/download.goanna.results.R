#' Download results from GOanna website.
#'
#' @param result Results file produced by upload.fasta.to.goanna()
#' @param goanna.zip.file.name Name to call Goanna output as a zip file 
#' @return goannaoutput
#' @examples
#' goannaoutput<- download.goanna.results(results = results,
#'                                        goanna.zip.file.name = "mRNA.test.zip")

download.goanna.results <- function (result,
                                     goanna.zip.file.name){
  result2 <- read_html(result)
  result3 <- as.character(result2)
  result4 <- sub(".+Your GOanna job has been submitted with job_id: (\\w+).+", "\\1", result3)
  download.file(paste0("http://www.agbase.msstate.edu/tmp/GOAL/", result4, ".zip"), destfile = goanna.zip.file.name)
  goannaoutput <- read.table(unz("mRNA.test.zip", paste0(result4,".sliminput.txt")))
  return(goannaoutput)
}