install.packages("R.utils")
library(R.utils)
BiocManager::install("Rsubread")
library(Rsubread)

install.packages("data.table")
library(data.table)

BiocManager::install("RUVSeq")
library(RUVSeq)

BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

install.packages("pheatmap")
library(pheatmap)

install.packages("RColorBrewer")
library(RColorBrewer)

install.packages("ggplot2")
library(ggplot2)

BiocManager::install("Rqc")
library(Rqc)

#Downloading the data

linkURL<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/006/SRR5924196/SRR5924196_1.fastq.gz"
file<-"SRR5924196_1.fastq.gz"
download.file(linkURL,file)

linkURL<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/006/SRR5924196/SRR5924196_2.fastq.gz"
file<-"SRR5924196_2.fastq.gz"
download.file(linkURL,file)

linkURL<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/008/SRR5924198/SRR5924198_1.fastq.gz"
file<-"SRR5924198_1.fastq.gz"
download.file(linkURL,file)

linkURL<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/008/SRR5924198/SRR5924198_2.fastq.gz"
file<-"SRR5924198_2.fastq.gz"
download.file(linkURL,file)


#Downloading the yeast genome (fasta and gtf file)
link<-"ftp://ftp.ensembl.org/pub/release-96/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
file<-"Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
download.file(link,file)
gunzip(file)

link<-"ftp://ftp.ensembl.org/pub/release-96/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.96.gtf.gz"
file<-"Saccharomyces_cerevisiae.R64-1-1.96.gtf.gz"
download.file(link,file)
gunzip(file)

# Quality control
install.packages("fastqcr")
library(fastqcr)
fastqc_install()
fastqc()
qc <- qc_aggregate("FASTQC/")
qc

# BUILDING genome Index

buildindex("Sc_full_index_rsubread",
           "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa",
           indexSplit=F)
reads1 <- list.files(pattern = "_1.fastq.gz$" )
reads2 <- list.files(pattern = "_2.fastq.gz$" )
all.equal(length(reads1),length(reads2))


#performing Alignment
align(index="Sc_full_index_rsubread",
      readfile1=reads1,
      readfile2=reads2,
      input_format="gzFASTQ",
      output_format="BAM",
      nthreads=10)

# checking the mapping quality
checkbam.files <- list.files(pattern = ".BAM$", full.names = TRUE)

check<- propmapped(files=checkbam.files)
check

#feature counts regeneration
featurec <- featureCounts(file= checkbam.files, annot.ext = "Saccharomyces_cerevisiae.R64-1-1.96.gtf", GTF.featureType = "exon",GTF.attrType = "gene_id",isGTFAnnotationFile = TRUE, isPairedEnd = TRUE)
featurdata<- data.frame(featurec[["counts"]])
colnames(featurdata) <- c("Normal","Tumor")

