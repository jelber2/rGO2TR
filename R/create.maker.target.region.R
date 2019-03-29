#' Filters Granges object by retained mRNA list to create MAKER target region
#' 
#' @param retained.mRNA.list mRNA retained by \code{filter.mRNA.GO.list} function
#' @param gff3.filtered Granges object to filter by \code{retained.mRNA.list}
#' @return target.region in the form of a Granges object with overlapping intervals
#' @examples
#' target.region <- create.maker.target.region(retained.mRNA.list, gff3.filtered)
#' @export
create.maker.target.region <- function (retained.mRNA.list, gff3.filtered) {
  df <- data.frame(seqnames = GenomeInfoDb::seqnames(gff3.filtered), 
                   starts = GenomicRanges::start(gff3.filtered),
                   ends = GenomicRanges::end(gff3.filtered), 
                   strands = GenomicRanges::strand(gff3.filtered),
                   tags = gff3.filtered$tags)
  
  df$seqnames <- as.character(df$seqnames)
  df$strands <- paste(df$strands)
  df$tags <- as.character(df$tags)
  
  j = 0
  target.region <- data.frame()
  df$tags <- sub(";",
                 "",
                 sub(".+Parent=(\\w+_\\d+\\-R\\w)",
                     "\\1",
                     df$tags,
                     perl=TRUE))
  
  for (i in retained.mRNA.list) {
    j = j + 1
    target.region <- rbind(target.region, df[grepl(i,
                                                   df$tags), ])
    cat("Finished filtering retained mRNA", j, "of", length(retained.mRNA.list), 
        "finished", "\n")
  }
  target.region <- GenomicRanges::makeGRangesFromDataFrame(target.region,
                                            keep.extra.columns=TRUE,
                                            ignore.strand=FALSE,
                                            seqinfo=NULL,
                                            seqnames.field="seqnames",
                                            start.field="starts",
                                            end.field="ends",
                                            strand.field="strands",
                                            starts.in.df.are.0based=FALSE)
  return(target.region)
}