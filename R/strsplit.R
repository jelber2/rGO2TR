#' @export
strsplit2 <-function(x, split=" ", n=1, ...)
 # from Chris Stubben's genomes2 package modified by Jean Elbers
 # https://github.com/cstubben/genomes2/blob/master/R/strsplit2.R
{
  y <- strsplit(as.character( x), split=split, ...)
  sapply(y, "[", n=n)
}
