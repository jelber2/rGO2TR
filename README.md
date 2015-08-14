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



##### To install devtools and rvest:

    install.packages("devtools")
    install.packages("rvest")



##### To install Bioconductor packages:

    source("http://bioconductor.org/biocLite.R")
    biocLite("genomes")
    biocLite("GOstats")



##### Then use the following command to actually install rGO2TR package

    devtools::install_github("jelber2/rGO2TR")



### Usage

##### From within R, load the required libraries

    library("rvest")
    library("GOstats")
    library("genomes")
    library("rGO2TR")


##### Now you are ready to use the package to filter genomes by gene ontology!



###### 1. Get available genomes with the get.ncbi.annot.euk.genomes function

    euks.filtered <- get.ncbi.annot.euk.genomes


###### 2. Search for your desired genome with search.annot.euks function

    query.results <- search.annot.euks("Zonotrichia albicollis")


###### 3. Download desired gff3 gene annotation file with get.gff3 function

    gff3 <- get.gff3(query.results, "sparrow.genome.gff3.gz")


###### 4a. Get mRNA accession ids from gene annotations with get.mRNA.acc function

    mRNA.acc <- get.mRNA.acc(gff3.filtered)


###### 4b. Determine the length of mRNA.acc and split if necessary
    # GOanna is limited in how many accessions it can process at one time
    # code to determine length of mRNA.acc
    length(mRNA.acc)

    # code to split mRNA.acc into smaller lists if >15,000 gi's
    mRNA.acc.split <- split(mRNA.acc,
                            ceiling(seq_along(mRNA.acc)/15000))
    mRNA1.acc <- mRNA.acc.split[[1]]
    mRNA2.acc <- mRNA.acc.split[[2]]


###### 5. Translate the mRNA using get.mRNA.translated function

    mRNA1.translated <- get.mRNA.translated(mRNA1.acc,
                                            mRNA1.translated,
                                            "jelber2@lsu.edu",
                                            "sparrow.mRNA1.fasta")

    mRNA2.translated <- get.mRNA.translated(mRNA2.acc,
                                            mRNA2.translated,
                                            "jelber2@lsu.edu",
                                            "sparrow.mRNA2.fasta")


###### 6a. Assign gene ontology using [GOanna](http://www.agbase.msstate.edu/cgi-bin/tools/GOanna.cgi)
    >Parameters for BLAST search

        Program	blastp

        Email Address	your_email_address

        Input File format	FASTA

        File to Upload	sparrow.mRNA1.fasta

        Database	AgBase-UniProt

        Filter sequences and annotations from the selected databases to remove sequences 
             with no GO annotations or with IEA or ND annotations only	Selected

        Expect	10e-20

        Matrix	Blosum62

        Gap Costs	Existence11 Extension1

        Word Size	3

        Filter	Low complexity not selected

        Nbr. Target Seqs	3

        Pct Id. Filter	10

        Query Coverage Fil	10

        Blast results format selection	TSV format

        Type of Evidence to Return	Experimental Evidence Codes(EXP,IDA,IPI,IMP,IGI,IEP)

    >Click on BLAST

    >Repeat search using above parameters for sparrow.mRNA2.fasta file

    >Can only submit 3 searches at a time with a single email account

    >Save raw results from email and place in desired working directory

    >Unzip each zip file. You want the "*.sliminput.txt" files


###### 6b. Read in the sliminput.txt files into R using the following commands:

    goannaoutput1 <- read.table("sparrow_set1.sliminput.txt")
    goannaoutput2 <- read.table("sparrow_set2.sliminput.txt")
    goannaoutput <- rbind(goannaoutput1, goannaoutput2)
    write.table(goannaoutput,
                file = "sparrow_annot.sliminput.txt",
                append = FALSE,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                row.names = FALSE,
                col.names = FALSE)


###### 7. Make GO id list for desired gene ontology term

    # example using GO term for pigment = GO:0043473
    GO.term <- "GO:0043473" # sets GO term as pigment
    GO.term.descendants <- GOBPOFFSPRING$"GO:0043473" # gets descendants for pigment
    GO.id.list <- c(GO.term, GO.term.descendants) # makes GO id list


###### 8. Make mRNA-GO id list with make.mRNA.GO.list function

    mRNA.GO.list <- make.mRNA.GO.list("sparrow_annot.sliminput.txt",
                                      "sparrow.mRNA.GO.list.txt")


###### 9. Filter the mRNA-GO id list with filter.mRNA.GO.list function

    retained.mRNA.list <- filter.mRNA.GO.list(mRNA.GO.list, GO.id.list)


###### 10. Create an initial target region with create.target.region function

    target.region <- create.target.region(retained.mRNA.list, gff3.filtered)


###### 11. Remove overlapping intervals

    final.target.region <- reduce(target.region)


###### 12. Estimate the size of the final non-overlapping target region

    cat("The target region after removing overlaps is",
    sum(width(final.target.region)),
    "bp.")
