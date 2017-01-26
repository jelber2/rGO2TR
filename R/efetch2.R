#' Retrieve Entrez database records at NCBI in a variety of formats
#' code from Chris Stubben's genomes2 package modified by Jean Elbers
#' https://github.com/cstubben/genomes2/blob/master/R/efetch.R 
#' @param id An EntrezHistory object or vector of Ids
#' @param db An Entrez database, default pubmed
#' @param rettype Retrieval type, see note for details
#' @param retmode Retrieval mode, see note for details
#' @param showURL display URL string
#' @param destfile location to save downloaded file using download.file. If missing, the url is loaded into R using readLines
#' @param ... Other key-value pairs passed to the efetch url string, e.g seq_stop
#' @return Returns efetch query
#' @examples
#' efetch.out <- efetch2(input.R.object.split[[i]], "nucleotide", "gb", "xml")
efetch2 <- function (id, db = "pubmed", rettype = "", retmode = "text", 
          showURL = FALSE, destfile, ...) 
{
  email <- Sys.getenv("email")
  if (email == "") {
    print("WARNING: please set your email using Sys.setenv(email='name@email.com')")
  }
  if (class(id)[1] == "EntrezHistory") {
    opts <- c(db = id$db, query_key = id$query_key, WebEnv = id$WebEnv)
  }
  else {
    id <- gsub(" ", "", id)
    if (is.vector(id)) 
      id <- paste(id, collapse = ",")
    opts <- c(id = id, db = db)
  }
  opts <- c(email = email, tool = "efetch.R", opts, rettype = rettype, 
            retmode = retmode, ...)
  opts <- paste(paste(names(opts), opts, sep = "="), collapse = "&")
  if (any(duplicated(names(opts)))) {
    stop("Duplicated keys are not allowed in url strings")
  }
  fetch <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  fetch <- paste(fetch, opts, sep = "?")
  if (showURL) 
    print(fetch)
  if (missing(destfile)) {
    content <- getURL(fetch)
    return(content)
  }
  else {
    download.file(fetch, destfile)
  }
}