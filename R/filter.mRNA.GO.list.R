#' Filters mRNA-GO list by GO id list to create Retained mRNA list
#' 
#' @param mRNA.GO.list R Object name of mRNA-GO list to filter
#' @param GO.id.list R Object name of GO id list to filter by
#' @return retained.mRNA.list filtered by \code{GO.id.list}
#' @examples
#' retained.mRNA.list <- filter.mRNA.GO.list(mRNA.GO.list, GO.id.list)
filter.mRNA.GO.list <- function (mRNA.GO.list, GO.id.list) {
  retained.mRNA.list <- data.frame(cbind("V1" = "NA",
                                         "V2" = "NA"),
                                   stringsAsFactors=FALSE)
  j = 0
  for (i in GO.id.list)
  {
    j = j+1
    # for each GO id in GO.id.list, search through the mRNA.GO.list and keep
    # only matches, save them to retained.mRNA.list, and sequentially add rows
    retained.mRNA.list <- data.frame(rbind(retained.mRNA.list,
                                           mRNA.GO.list[grepl(i,
                                                              mRNA.GO.list$V2),
                                                        ]),
                                     stringsAsFactors=FALSE)
    cat("GO id", j, "of", length(GO.id.list), "finished", "\n")
  }
  retained.mRNA.list <- retained.mRNA.list$V1
  retained.mRNA.list <- retained.mRNA.list[-1]
  retained.mRNA.list <- unique(retained.mRNA.list)
  return(retained.mRNA.list)
}