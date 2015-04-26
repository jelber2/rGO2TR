#' Makes mRNA-GO.list
#' 
#' @param path.to.GOanna.output.file Full path to GOanna output file
#' @param output.file.name Name to call output file
#' @return mRNA GO list
#' @examples
#' mRNA.GO.list <- make.mRNA.GO.list("sparrow_annot.sliminput.txt",
#'                                   "sparrow.mRNA.GO.list.txt")
make.mRNA.GO.list <- function (path.to.GOanna.output.file, output.file.name) {
  # read in the data from GOanna as mRNA GO list
  mRNA.GO.list <- read.table(path.to.GOanna.output.file)
  
  # remove 3rd row from protein GO list
  mRNA.GO.list <- mRNA.GO.list[-3]
  
  # create the mRNA.GO.list
  write.table(mRNA.GO.list,
              file = output.file.name,
              append = FALSE,
              quote = FALSE,
              sep = "\t",
              eol = "\n",
              row.names = FALSE,
              col.names = FALSE)
  return(mRNA.GO.list)
}