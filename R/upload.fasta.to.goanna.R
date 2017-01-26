#' Upload protein FASTA files to GOanna website.
#'
#' @param email.address Email address to send the output of GOanna as a string
#' @param file.to.upload Path to protein FASTA file to upload to GOanna as a string
#' @param expected.value Expected value or evalue for BLAST search as a string
#' @param word.size Word size for the BLAST search as a string
#' @param max.target.sequences Maximum number of target sequences for BLAST search
#' @param percent.identity Percent Identity of the BLAST search as a string
#' @param query.coverage Query coverage of BLAST search as a string
#'
#' @examples
#' upload.fasta.to.goanna(email.address ="email",
#'                        file.to.upload = "sparrow.mRNA1.fasta",
#'                        expected.value = "10e-20",
#'                        word.size = "3",
#'                        max.target.sequences = "3",
#'                        percent.identity = "20",
#'                        query.coverage = "20")
upload.fasta.to.goanna <- function(email.address,
                                   file.to.upload,
                                   expected.value,
                                   word.size,
                                   max.target.sequences,
                                   percent.identity,
                                   query.coverage) {

  url       <-"http://agbase.msstate.edu/cgi-bin/tools/GOanna.cgi"   ## page to spider
  pgsession <-html_session(url)               ## create session
  pgform    <-html_form(pgsession)[[1]]       ## pull form from session
  filled_form <- set_values(pgform,
                            pgform$url <- "",
                            pgform$fields$EMAIL$value <- email.address,
                            pgform$fields$IDLIST$value <- upload_file(file.to.upload),
                            pgform$fields$EXPECT$value <- expected.value,
                            pgform$fields$WORD_SIZE$value <- word.size,
                            pgform$fields$MAX_TARGET_SEQS$value <- max.target.sequences,
                            pgform$fields$PCTID$value <- percent.identity,
                            pgform$fields$QRYCOV$value <- query.coverage,
                            pgform$fields$blastp_query_cov$value <- query.coverage)
  submit_form(pgsession,filled_form)
}