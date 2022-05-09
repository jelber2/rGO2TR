RGO2TR update for NCBI ANNOTATED genomes using UniProt/Swiss-Prot

```bash
mkdir -p ~/cheetah

cd ~/cheetah

# see https://github.com/mamba-org/mamba about mamba installation

# create mamba/conda environment
mamba create -n cheetah -c conda-forge -c bioconda \
r-tidyr bioconductor-GO.db \
r-XML bioconductor-GenomicRanges \
bioconductor-GenomeInfoDb=1 \
bioconductor-GOstats bioconductor-graph bioconductor-Category \
r-base=3.6.3 \
r-Matrix \
bioconductor-AnnotationDbi \
bioconductor-IRanges \
bioconductor-Biobase \
bioconductor-BiocGenerics \
r-rvest \
r-seqinr \
r-RCurl \
r-httr \
conda-ecosystem-user-package-isolation \
r-devtools \
r-data.table \
blast \
parallel

# activate conda environment
conda activate cheetah

# enter R shell
R
```

```R
# install rGO2TR
library("devtools")
devtools::install_github("jelber2/rGO2TR",force=TRUE)

# type in "3" to not update any libraries

# load libraries
library("data.table")
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


# get cheetah gff file from NCBI
system("wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/32536/101/GCF_003709585.1_Aci_jub_2/GCF_003709585.1_Aci_jub_2_genomic.gff.gz")

# unzip
system("gunzip GCF_003709585.1_Aci_jub_2_genomic.gff.gz")

# read in gff3
gff3 <-read.maker.gff3("./GCF_003709585.1_Aci_jub_2_genomic.gff")

# filter gff3
gff3.filtered <- filter.gff3(gff3, 'exon')

# Get mRNA accession ids from gene annotations with get.mRNA.acc function
mRNA.acc <- get.mRNA.acc(gff3.filtered)

# this step takes a long time
# Get mRNA translated
mRNA.translated <- get.mRNA.translated(mRNA.acc,
                                        mRNA.translated,
                                        "jean.elbers@gmail.com", # put your email address here
                                        "mRNA.fasta")


#A. Get UniProt/Swiss-Prot proteins

system("wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz")

#B. Decompress the file

system("gunzip uniprot_sprot.fasta.gz")

#C. Make BLAST DB

system("makeblastdb -dbtype prot -in uniprot_sprot.fasta")

#D. Perform BLAST search in parallel, requiring at least an E value of 1e-20 or less
# for jobs, enter your number of maximum threads/cores

# this step takes a very long time
system("cat mRNA.fasta | parallel --no-notice --jobs 96 --block 10k --recstart '>' --pipe blastp -evalue 1e-20 -outfmt 6 -db uniprot_sprot.fasta -query - >> blastp.result")



#E. Get the first (i.e., "best") BLAST hit for each protein

system("cat blastp.result | awk '!seen[$1]++' > blastp.result.filtered")

#F. Get unique UniProtKB entry-ids from the BLAST search (you will use entry-ids to get gene ontology term ids (i.e., go-ids) for these later)

system("cut -f 2 blastp.result.filtered |cut -f 2 -d '|' |sort -u > UniProtKB-entry-ids-unique.txt")

#G. Make a list of your proteins compared to UniProtKB ids from the BLAST search

system("cut -f 1 blastp.result.filtered > UniProtKB-protein-ids.txt")

system("cut -f 2 blastp.result.filtered |cut -f 2 -d '|' > UniProtKB-entry-ids.txt")

system("paste UniProtKB-protein-ids.txt UniProtKB-entry-ids.txt  > UniProtKB-protein-entry-list.txt")

#H. Use the following code to fetch the go-ids for each entry-id

system('while read i;do wget -q -O - "http://www.uniprot.org/uniprot/?query=${i}&format=tab&columns=id%2Cgo-id" |grep -v "Entry" >> UniProtKB-entry-ids-unique-go-ids.txt ;done < UniProtKB-entry-ids-unique.txt')




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

# example using GO term for immune response = GO:0006955
GO.term <- "GO:0006955" # sets GO term as immune response
GO.term.descendants <- GOBPOFFSPRING$"GO:0006955" # gets descendants for immune response
GO.id.list <- c(GO.term, GO.term.descendants) # makes GO id list

retained.protein.list <- filter.mRNA.GO.list(protein.go.id.df, GO.id.list)

target.region <- create.maker.target.region(retained.protein.list, gff3.filtered)


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

final.target.region <- reduce(target.region)

cat("The target region after removing overlaps is",
    sum(width(final.target.region)),
    "bp.")

save.target.region.as.bed.file(final.target.region, "cheetah.immunome.bed")

