#' rGO2TR: filters NCBI Annotated Eukaryotic Genomes by gene onotology.
#' 
#' The rGO2TR package is based off of the GO2TR:
#' \strong{G}ene \strong{O}ntology \strong{to} \strong{T}arget \strong{R}egion
#' workflow and includes functions for manipulating gene annotations, 
#' acquiring mRNA accession identifiers, uploading data to an annoation service,
#' creating a gene ontology list, and filtering gene annotations by
#' gene ontology.
#' 
#' @section rGO2TR functions:
#' Note: All functions are listed in chronological order.
#' 
#' The gene annotation functions include: 
#' \code{\link{get.ncbi.annot.euk.genomes}} - gets all available genomes 
#' annotated by NCBI Eukaryotic Genome Annotation Pipeline, 
#' \code{\link{search.annot.euks}} - searches for a specified species, and 
#' \code{\link{get.gff3}} - downloads the gene annotations for the 
#' specified species.
#' 
#' The mRNA accession functions include: \code{\link{filter.gff3}} - filters 
#' gene annotations by desired source (i.e., CDS, mRNA, exon, etc.),
#' \code{\link{get.mRNA.acc}} - searches through filtered gene annotations for 
#' predicted mRNA accessions identifiers (i.e., "XM_"), and
#' \code{\link{get.mRNA.translated}} - returns the protein translated 
#' sequences for the mRNA in FASTA format.
#' 
#' The upload data function, \code{\link{upload.data.iplant}} uploads FASTA 
#' files to the iPlant Discovery Enivronment by using their Agave API.
#' 
#' The gene ontology list functions are \strong{not} coded as package
#' functions, rather they are simple commands such as the following example,
#' which creates a GO identitifer list (inlcudes descendant/offspring terms)
#' for the GO term pigment:
#' 
#' \code{GO.term <- "GO:0043473" # sets pigment as base GO term}
#' 
#' \code{GO.term.descendants <- GOBPOFFSPRING$"GO:0043473" # the  
#' offspring/descendant GO terms for pigment}
#' 
#' \code{GO.id.list <- c(GO.term, GO.term.descendants) # combines base 
#' and descendant terms}
#' 
#' 
#' The filtering gene annotation by gene onotlogy functions include:
#' \code{\link{make.mRNA.GO.list}} - processes output from the GOanna 
#' program, which is used to assign gene ontology terms to mRNA,
#' \code{\link{filter.mRNA.GO.list}} - filters themRNA.GO.list by GO.id.list
#' to create a retained mRNA list to eventually filter the gene annotations,
#' \code{\link{create.target.region}} - filters the gene annotations by the 
#' retained mRNA, effectively creating a pigmentome for example, and 
#' \code{\link{save.target.region.as.bed.file}} exports the generated target 
#' region into the bed (\strong{B}rowser \strong{E}xtensible \strong{D}ata) 
#' format.
#' 
#' @docType package
#' @name rGO2TR
NULL