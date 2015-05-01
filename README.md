# rGO2TR
### About
    The rGO2TR package is based off of the GO2TR: Gene Ontology to Target Region workflow
    and includes functions for manipulating gene annotations, acquiring mRNA accession
    identifiers, uploading data to an annoation service, creating a gene ontology list,
    and filtering gene annotations by gene ontology.

### Installation
#### To install rGO2TR package from within R, make sure you have the following packages:
    devtools (from CRAN)
    genomes (from Bioconductor)
    GOstats (from Bioconductor)
    rvest (from CRAN)
    rPlant (>= 2.10.5, not yet on CRAN but on rPlant's R-Forge repository)
    
##### To install devtools and rvest:
    install.packages("devtools")
    install.packages("rvest")
    
##### To install Bioconductor packages:
    source("http://bioconductor.org/biocLite.R")
    biocLite("genomes")
    biocLite("GOstats")

##### To install rPlant version 2.10.5 from R-Forge repository:
    install.packages("rPlant", repos="http://R-Forge.R-project.org", type = "source")

##### Then use the following command to actually install rGO2TR package
    devtools::install_github("jelber2/rGO2TR", auth_token = "8d5c8580c1dfe4aff9e99d271eebece7e53c6fc4")

### Usage
##### From within R, load the required libraries
    library("rvest")
    library("GOstats")
    library("genomes")
    library("rPlant")
    library("rGO2TR")
    
##### Now you are ready to use the package to filter genomes by gene ontology!

###### 1. Get available genomes with the get.ncbi.annot.euk.genomes function
    euks.filtered <- get.ncbi.annot.euk.genomes
###### 2. Search for your desired genome with search.annot.euks function
    query.results <- search.annot.euks("Zonotrichia albicollis")
###### 3. Download desired gff3 gene annotation file with get.gff3 function
    gff3 <- get.gff3(query.results, "sparrow.genome.gff3.gz")
###### 4. Get mRNA accession ids from gene annotations with get.mRNA.acc function
    mRNA.acc <- get.mRNA.acc(gff3.filtered)
###### 5. Translate the mRNA using get.mRNA.translated function
    mRNA1.translated <- getmRNAtranslated(mRNA1.acc,
                                          mRNA1.translated,
                                          "jelber2@lsu.edu",
                                          "sparrow.mRNA1.fasta")
###### 6. Upload mRNA translated sequences to iPlant using upload.data.iplant function
    upload.data.iplant("user.name",
                       "password",
                       "sparrow.mRNA1.fasta",
                       "C:/Users/jelber2/white_throated_sparrow/")

###### 7. Make GO id list for desired gene ontology term
    # example using GO term for pigment = GO:0043473
    GO.term <- "GO:0043473" # sets GO term as pigment
    GO.term.descendants <- GOBPOFFSPRING$"GO:0043473" # gets descendants for pigment
    GO.id.list <- c(GO.term, GO.term.descendants) # makes GO id list

###### 8.

###### 9.
    
###### 10.
        
![Image of Discovery Environment](https://github.com/jelber2/rGO2TR/blob/master/images/1.PNG)
