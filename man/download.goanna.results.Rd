% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.goanna.results.R
\name{download.goanna.results}
\alias{download.goanna.results}
\title{Download results from GOanna website.}
\usage{
download.goanna.results(result, goanna.zip.file.name)
}
\arguments{
\item{result}{Results file produced by upload.fasta.to.goanna()}

\item{goanna.zip.file.name}{Name to call Goanna output as a zip file}
}
\value{
goannaoutput
}
\description{
Download results from GOanna website.
}
\examples{
goannaoutput<- download.goanna.results(results = results,
                                       goanna.zip.file.name = "mRNA.test.zip")
}
