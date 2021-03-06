% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/upload.fasta.to.goanna.R
\name{upload.fasta.to.goanna}
\alias{upload.fasta.to.goanna}
\title{Upload protein FASTA files to GOanna website.}
\usage{
upload.fasta.to.goanna(email.address, file.to.upload, expected.value,
  word.size, max.target.sequences, percent.identity, query.coverage)
}
\arguments{
\item{email.address}{Email address to send the output of GOanna as a string}

\item{file.to.upload}{Path to protein FASTA file to upload to GOanna as a string}

\item{expected.value}{Expected value or evalue for BLAST search as a string}

\item{word.size}{Word size for the BLAST search as a string}

\item{max.target.sequences}{Maximum number of target sequences for BLAST search}

\item{percent.identity}{Percent Identity of the BLAST search as a string}

\item{query.coverage}{Query coverage of BLAST search as a string}
}
\value{
results = a status of the GOanna webpage request
}
\description{
Upload protein FASTA files to GOanna website.
}
\examples{
result <- upload.fasta.to.goanna(email.address ="email",
                       file.to.upload = "sparrow.mRNA1.fasta",
                       expected.value = "10e-20",
                       word.size = "3",
                       max.target.sequences = "3",
                       percent.identity = "20",
                       query.coverage = "20")
}
