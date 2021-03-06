% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.maker.mRNA.acc.R
\name{get.maker.mRNA.acc}
\alias{get.maker.mRNA.acc}
\title{Gets mRNA accession identifiers from MAKER derived Granges object.}
\usage{
get.maker.mRNA.acc(gff3.filtered)
}
\arguments{
\item{gff3.filtered}{An already filtered Granges object}
}
\value{
List of mRNA accession identifiers contained in gff3.filtered
}
\description{
Gets mRNA accession identifiers from MAKER derived Granges object.
}
\examples{
mRNA.acc <- get.maker.mRNA.acc(gff3.filtered)
}
