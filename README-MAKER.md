# rGO2TR using MAKER Annotated Genomes

## FOR MAKER ANNOTATED GENOMES FOLLOW THESE INSTRUCTIONS

### 1. Read in MAKER produce gff3 file

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

### 6. [Follow step 6 and onwards on the main README.MD page] (https://github.com/jelber2/rGO2TR/blob/master/README.md)