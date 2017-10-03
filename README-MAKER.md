# rGO2TR using MAKER Annotated Genomes

## FOR MAKER ANNOTATED GENOMES FOLLOW THESE INSTRUCTIONS

### 1. Read in MAKER produce gff3 file (can be gzipped)

    gff3 <-read.maker.gff3("camel.gff3")


### 2. Filter gene annotations by desired source (i.e., CDS, mRNA, exon, etc.)

    gff3.filtered <- filter.maker.gff3(gff3, 'exon')

### 3. Get mRNA accession ids from gene annotations with get.maker.mRNA.acc function

    mRNA.acc <- get.maker.mRNA.acc(gff3.filtered)


### 4. Determine the length of mRNA.acc and split if necessary
    # GOanna is limited in how many accessions it can process at one time
    # code to determine length of mRNA.acc
    length(mRNA.acc)

    # code to split mRNA.acc into smaller lists if >15,000 gi's
    mRNA.acc.split <- split(mRNA.acc,
                            ceiling(seq_along(mRNA.acc)/15000))
    mRNA1.acc <- mRNA.acc.split[[1]]
    mRNA2.acc <- mRNA.acc.split[[2]]


### 5. Read in protein translations

    mRNA1.translated <- get.maker.mRNA.translated(mRNA1.acc,
                                              mRNA1.translated,
                                              "maker3_proteins_0-0.74_AED.fasta",
                                              "mRNA1.translated.fasta")

    mRNA2.translated <- get.maker.mRNA.translated(mRNA2.acc,
                                                  mRNA2.translated,
                                                  "maker3_proteins_0-0.74_AED.fasta",
                                                  "mRNA2.translated.fasta")

### 6a. Assign gene ontology using [GOanna](http://www.agbase.msstate.edu/cgi-bin/tools/GOanna.cgi)

    results1 <- upload.fasta.to.goanna(email.address ="email",
                                       file.to.upload = "mRNA1.translated.fasta",
                                       expected.value = "10e-20",
                                       word.size = "3",
                                       max.target.sequences = "3",
                                       percent.identity = "20",
                                       query.coverage = "20")

    results2 <- upload.fasta.to.goanna(email.address ="email",
                                       file.to.upload = "mRNA2.translated.fasta",
                                       expected.value = "10e-20",
                                       word.size = "3",
                                       max.target.sequences = "3",
                                       percent.identity = "20",
                                       query.coverage = "20")


### 6b. Read in the sliminput.txt files into R using the following commands:

      goannaoutput1 <- download.goanna.results(result = results1,
                                               goanna.zip.file.name = "mRNA1.goanna.zip")
      goannaoutput2 <- download.goanna.results(result = results2,
                                               goanna.zip.file.name = "mRNA2.goanna.zip")

    goannaoutput <- rbind(goannaoutput1, goannaoutput2)
    write.table(goannaoutput,
                file = "camel_annot.sliminput.txt",
                append = FALSE,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                row.names = FALSE,
                col.names = FALSE)


### 7. Make GO id list for desired gene ontology term

    # example using GO term for immune response = GO:0006955
    GO.term <- "GO:0006955" # sets GO term as immune response
    GO.term.descendants <- GOBPOFFSPRING$"GO:0006955" # gets descendants for immune response
    GO.id.list <- c(GO.term, GO.term.descendants) # makes GO id list


### 8. Make mRNA-GO id list with make.mRNA.GO.list function

    mRNA.GO.list <- make.mRNA.GO.list("camel_annot.sliminput.txt",
                                      "camel.mRNA.GO.list.txt")


### 9. Filter the mRNA-GO id list with filter.mRNA.GO.list function

    retained.mRNA.list <- filter.mRNA.GO.list(mRNA.GO.list, GO.id.list)


### 10a. Create an initial target region with create.target.region function

    target.region <- create.maker.target.region(retained.mRNA.list, gff3.filtered)


### 10b. Estimate number of genes and exons

    # how many unique genes are there?
    cat("There are at least",
    length(unique(sub("(\\w+_\\d+)\\-R\\w",
                             "\\1",
                             unlist(strsplit(target.region$tags, ",")),
                             perl=T))),
    "unique genes in the target region, but possibly more because genes in this 
    target region might overlap genes in the gff3 file.")


    # how many exons are there?
    cat("There are at least",
    length(target.region$tags),
    "exons in the target region, but possibly more because exons in this 
    target region might overlap exons in the gff3 file.")


### 11. Merge overlapping intervals

    final.target.region <- reduce(target.region)


### 12. Estimate the size of the final non-overlapping target region

    cat("The target region after removing overlaps is",
    sum(width(final.target.region)),
    "bp.")


### 13. Convert final target region to a BED file

    save.target.region.as.bed.file(final.target.region, "camel.immunome.bed")