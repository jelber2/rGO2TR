# rGO2TR using UniProt/Swiss-Prot

## FOR NCBI ANNOTATED GENOMES FOLLOW THESE INSTRUCTIONS

### Work in progress

## FOR MAKER ANNOTATED GENOMES FOLLOW THESE INSTRUCTIONS

### NOTE ASSUMES THAT YOUR MAKER IDs are in the form Cadr_111111-RA (regular expression= \w+_\d+-R\w)

### 1. Read in MAKER produce gff3 file (can be gzipped)

    gff3 <-read.maker.gff3("camel.gff3")


### 2. Filter gene annotations by desired source (i.e., CDS, mRNA, exon, etc.)

    gff3.filtered <- filter.maker.gff3(gff3, 'exon')

### 3. On a Linux computer, BLASTP your MAKER predicted proteins with UniProt/Swiss-Prot proteins

    #A. Get UniProt/Swiss-Prot proteins
         wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    #B. Decompress the file
         gunzip uniprot_sprot.fasta.gz
    #C. Make BLAST DB
         makeblastdb -dbtype pep -in uniprot_sprot.fasta
    #D. Perform BLAST search, requiring at least an E value of 1e-20 or less
         blastp -num_threads enter-the-number-of-cores-you-have-on-your-computer-here -db uniprot_sprot.fasta \
         -evalue 1e-20 -query maker.proteins.fasta -outfmt 6 > blastp.result
    #E. Get the first (i.e., "best") BLAST hit for each protein
         cat blastp.result | awk '!seen[$1]++' > blastp.result.filtered
    #F. Get unique UniProtKB entry-ids from the BLAST search (you will use entry-ids to get gene ontology term ids (i.e., go-ids) for these later)
         cut -f 2 blastp.result.filtered |cut -f 2 -d "|" |sort -u > UniProtKB-entry-ids-unique.txt
    #G. Make a list of your proteins compared to UniProtKB ids from the BLAST search
         cut -f 1 blastp.result.filtered > UniProtKB-protein-ids.txt
         cut -f 2 blastp.result.filtered |cut -f 2 -d "|" > UniProtKB-entry-ids.txt
         paste UniProtKB-protein-ids.txt UniProtKB-entry-ids.txt  > UniProtKB-protein-entry-list.txt
    #H. Use the following code to fetch the go-ids for each entry-id
         while read i;do
         wget -q -O - "http://www.uniprot.org/uniprot/?query=${i}&format=tab&columns=id%2Cgo-id" |grep -v "Entry" >> UniProtKB-entry-ids-unique-go-ids.txt
         done < UniProtKB-entry-ids-unique.txt

### 4. Use rGO2TR to manipulate the UniProt output files

    #A. Read in the data     
         uniprot.output <- read.table("UniProtKB-entry-ids-unique-go-ids.txt",sep = "\t",header=F,na.strings = "")
    #B. Convert the column V1 as character, not as factors
         uniprot.output$V1 <- as.character(uniprot.output$V1)
    #C. Convert the Go-id column into as character, not as factors
         uniprot.output$V2 <- as.character(uniprot.output$V2)
    #D. Get rid of spaces in column2
         uniprot.output$V2 <- gsub(" ", "", uniprot.output$V2)
    #E. Use tidyr separate rows to  convert A1  GO:1,GO:2 to
      #                                     A1  GO:1
      #                                     A1  GO:2
         uniprot.output <- tidyr::separate_rows(data = uniprot.output,V2,sep = ";")

    #F. Get rid of uniprot proteins without any GO ids
         uniprot.output <- na.omit(uniprot.output)
    #G. Read in the BLAST results
    # column 1 is MAKER protein-id and column 2 is the UniProt entry-id
         protein.entry.df <- read.table("UniProtKB-protein-entry-list.txt")

    #H. Make these two columns as characters not factors
         protein.entry.df$V1 <- as.character(protein.entry.df$V1)
         protein.entry.df$V2 <- as.character(protein.entry.df$V2)

    #I. copy the uniprot data as a new dataframe called protein.go.id.df
         protein.go.id.df <- uniprot.output

    #J. Go through each pair of MAKER protein-ids and UniProt entry-ids and replace the UniProt entry-id with the MAKER protein-id
         for (i in 1:nrow(protein.entry.df)){
         protein.go.id.df$V1 <- sub(pattern=protein.entry.df$V2[i],
                             replacement=protein.entry.df$V1[i],
                             protein.go.id.df$V1)
         cat("Protein", i, "of", nrow(protein.entry.df), "\n")
         }


### 5. Make GO-id list for desired gene ontology term

    # example using GO term for immune response = GO:0006955
    GO.term <- "GO:0006955" # sets GO term as immune response
    GO.term.descendants <- GOBPOFFSPRING$"GO:0006955" # gets descendants for immune response
    GO.id.list <- c(GO.term, GO.term.descendants) # makes GO id list


### 6. Filter the protein.go.id.df with filter.mRNA.GO.list function

    retained.protein.list <- filter.mRNA.GO.list(protein.go.id.df, GO.id.list)


### 7a. Create an initial target region with create.target.region function

    target.region <- create.maker.target.region(retained.protein.list, gff3.filtered)


### 7b. Estimate number of genes and exons

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


### 7c. Merge overlapping intervals

    final.target.region <- reduce(target.region)


### 7d. Estimate the size of the final non-overlapping target region

    cat("The target region after removing overlaps is",
    sum(width(final.target.region)),
    "bp.")


### 8. Convert final target region to a BED file

    save.target.region.as.bed.file(final.target.region, "camel.immunome.bed")