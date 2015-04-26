#' Uploads translated mRNA sequences to the iPlant's Discovery Environment
#' 
#' @param user.name Your iPlant user name
#' @param password You iPlant password
#' @param file.name.to.upload Name of file to upload
#' @param local.file.path Path of \code{file.name.to.upload}
#' @examples
#' uploadDataToiPlant ("jelber2",
#'                    "Gopherus2011",
#'                    "sparrow.mRNA1.fasta",
#'                    "C:/Users/jelber2/Dropbox/LSU/Dissertation/Manuscripts/GO2OME/rGO2TR/white_throated_sparrow/")
#' uploadDataToiPlant ("jelber2",
#'                    "Gopherus2011",
#'                    "sparrow.mRNA2.fasta",
#'                    "C:/Users/jelber2/Dropbox/LSU/Dissertation/Manuscripts/GO2OME/rGO2TR/white_throated_sparrow/")
upload.data.iplant <- function(user.name,
                               password,
                               file.name.to.upload,
                               local.file.path) {
  # logs into iPlant Agave API
  Validate(user.name, password, api="agave", print.curl=FALSE)
  
  # uploads desired file to iPlant Discovery Environment
  # uploaded file is available at the following directory
  # once logged into the Discovery Environment
  # /iplant/home/user.name/
  UploadFile(file.name.to.upload, 
             local.file.path,
             filetype="FASTA-0")
  
  # Let's user know upload is complete.
  cat("\n", file.name.to.upload, "was successfully uploaded", "\n")
}