# rGO2TR

### About

    The rGO2TR package is based off of the GO2TR: Gene Ontology to Target Region workflow
    and includes functions for manipulating gene annotations, acquiring mRNA accession
    identifiers, uploading data to an annoation service, creating a gene ontology list,
    and filtering gene annotations by gene ontology.



### Installation

#### To install rGO2TR package from within R, make sure you have the following packages:

    devtools (from CRAN)
    GenomicRanges (from Bioconductor)
    GOstats (from Bioconductor)
    GO.db (from Bioconductor)
    rvest (from CRAN)
    seqinr (from CRAN)
    tidyr (from CRAN)



##### To install devtools and rvest:

    install.packages("devtools")
    install.packages("rvest")
    install.packages("XML")
    install.packages("httr")
  	install.packages("seqinr")
  	install.packages("RCurl")
  	install.packages("tidyr")


##### To install Bioconductor packages:

    source("http://bioconductor.org/biocLite.R")
    biocLite("GenomicRanges")
    biocLite("GOstats")
    biocLite("GO.db")



##### Then use the following command to actually install rGO2TR package

    library("devtools")
    devtools::install_github("jelber2/rGO2TR")



### Usage

##### From within R, load the required libraries

    library("httr")
	library("RCurl")
	library("seqinr")
    library("rvest")
    library("GOstats")
    library("GenomicRanges")
    library("XML")
    library("GO.db")
    library("tidyr")
    library("rGO2TR")

### Citation
#### To cite the rGO2TR package in publications use:

    Elbers, J. P. and Taylor, S. S. (2015) GO2TR: a gene ontology-based workflow to generate target regions for target enrichment
    experiments Conservation Genetics Resources 7:851--857

#### A BibTeX entry for LaTeX users is

    @Article{,
    author = {Jean P. Elbers and Sabrina S. Taylor},
    title = {GO2TR: a gene ontology-based workflow to generate target regions for target enrichment experiments},
    journal = {Conservation Genetics Resources},
    year = {2015},
    volume = {7},
    pages = {851--857},
    doi = {10.1007/s12686-015-0487-6},
    }

#### As rGO2TR is continually evolving, you may want to cite its version number. Find it with
    help(package=rGO2TR) or
    sessionInfo()

##### Now you are ready to use the package to filter genomes by gene ontology!

### [If you don't have a license for GOanna because AgBase is now pay-to-use, follow these instructions](https://github.com/jelber2/rGO2TR/blob/master/README-UNIPROT.md)

### [FOR MAKER ANNOTATED GENOMES FOLLOW THESE INSTRUCTIONS](https://github.com/jelber2/rGO2TR/blob/master/README-MAKER.md)

### FOR NCBI ANNOTATED GENOMES FOLLOW THE INSTRUCTIONS BELOW
###### 1. Get available genomes with the get.ncbi.annot.euk.genomes function

    euks.filtered <- get.ncbi.annot.euk.genomes()


###### 2. Search for your desired genome with search.annot.euks function

    query.results <- search.annot.euks("Zonotrichia albicollis")


###### 3a. Download desired gff3 gene annotation file with get.gff3 function

    gff3 <- get.gff3(query.results, "sparrow.genome.gff3.gz")

###### 3b. Filter gene annotations by desired source (i.e., CDS, mRNA, exon, etc.)

    gff3.filtered <- filter.gff3(gff3, 'exon')

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
                                            "jelber2@lsu.edu", # put your email address here
                                            "sparrow.mRNA1.fasta")

    mRNA2.translated <- get.mRNA.translated(mRNA2.acc,
                                            mRNA2.translated,
                                            "jelber2@lsu.edu", # put your email address here
                                            "sparrow.mRNA2.fasta")


###### 6a. Assign gene ontology using [GOanna](http://www.agbase.arizona.edu/cgi-bin/tools/GOanna.cgi) In January 2019 AgBase will be switching to a non-profit, subscription model.
    Note1: the upload.fasta.to.goanna() function seems to be broken, so manual file uploading is necessary 
    Note2: as of rGO2TR version 1.0.9, there is now a function called
          upload.fasta.to.goanna() that allows you to upload protein
          fasta files generated by get.mRNA.translated() to GOanna
          It is called by using
          results1 <- upload.fasta.to.goanna(email.address ="email",
                                           file.to.upload = "sparrow.mRNA1.fasta",
                                           expected.value = "10e-20",
                                           word.size = "3",
                                           max.target.sequences = "3",
                                           percent.identity = "20",
                                           query.coverage = "20")

          results2 <- upload.fasta.to.goanna(email.address ="email",
                                           file.to.upload = "sparrow.mRNA2.fasta",
                                           expected.value = "10e-20",
                                           word.size = "3",
                                           max.target.sequences = "3",
                                           percent.identity = "20",
                                           query.coverage = "20")
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

        Pct Id. Filter	70 (70 is recommended by AgBase, but usually use 20)

        Query Coverage Fil	70 (70 is recommended by AgBase, but usually use 20)

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
    # Note as of rGO2TR version 1.0.10 there is now a download.goanna.results()
      function which is called using:
      goannaoutput1 <- download.goanna.results(result = results1,
                                               goanna.zip.file.name = "mRNA1.test.zip")
      goannaoutput2 <- download.goanna.results(result = results2,
                                               goanna.zip.file.name = "mRNA2.test.zip")

    
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


###### 10a. Create an initial target region with create.target.region function

    target.region <- create.target.region(retained.mRNA.list, gff3.filtered)


###### 10b. Estimate number of genes and exons

    # how many unique genes are there?
    cat("There are at least",
    length(unique(sub("XM_\\d+.\\d;(\\w+)","\\1",target.region$tags,perl=TRUE))),
    "unique genes in the target region, but possibly more because genes in this 
    target region might overlap genes in the gff3 file.")


    # how many exons are there?
    cat("There are at least",
    length(target.region$tags),
    "exons in the target region, but possibly more because exons in this 
    target region might overlap exons in the gff3 file.")


###### 11. Merge overlapping intervals

    final.target.region <- reduce(target.region)


###### 12. Estimate the size of the final non-overlapping target region

    cat("The target region after removing overlaps is",
    sum(width(final.target.region)),
    "bp.")


###### 13. Convert final target region to a BED file

    save.target.region.as.bed.file(final.target.region, "sparrow.pigmentome.bed")
