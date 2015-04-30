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
