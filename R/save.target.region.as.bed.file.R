#' Saves final.target.region as a bed file
#' 
#' @param final.target.region R object to save as BED file
#' @param output.file.name Name of BED file
#' @examples
#' save.target.region.as.bed.file(final.target.region, "sparrow.pigmentome.bed")
#' @export
save.target.region.as.bed.file <- function (final.target.region,
                                            output.file.name) {
  # creates the data.frame df to store seqnames,start,end,name,score,and strand
  # code modified from https://www.biostars.org/p/89341/
  df <- data.frame(seqnames=seqnames(final.target.region),
                   starts=start(final.target.region)-1,
                   ends=end(final.target.region),
                   names=c(rep(".", length(final.target.region))),
                   scores=c(rep(".", length(final.target.region))),
                   strands=strand(final.target.region))
  # write df to a tab-delimited output.file.name
  write.table(df,
              file=output.file.name,
              quote=FALSE,
              sep="\t",
              row.names=FALSE,
              col.names=FALSE)
  cat("\n",
      "Finished writing",
      output.file.name,
      "\n")
}