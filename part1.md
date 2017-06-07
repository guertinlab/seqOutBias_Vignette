---
layout: docs
title: Processing of DNase-seq data from ENCODE
prev: index
next: part2
number: 1
---

# Processing of DNase-seq data from ENCODE

Enzymatic digestion sequence preference was characterized on the genome-wide scale in a joint publication
from Shirley Liu’s and Myles Brown’s groups using DNase-seq data (He *et al.*, 2014). We use
publicly available DNase-seq data from ENCODE (Stamatoyannopoulos Lab) and data deposited into
GEO (Lazarovici *et al.*, 2013) as examples of how to use `seqOutBias` to correct enzymatic accessibility
data.

## Retrieving raw data from ENCODE
Download the *fastq* files directly from ENCODE ( [http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/) ).
Use wget to retrieve the raw fastq files from ENCODE, `tar` the downloaded files, combine all the
replicates for each condition and compress the resultant file. Note that we are combining all the data
sets for the purpose of having more sequencing depth of coverage and for ease of analysis downstream.
These files result from DNase-nicking of crude nuclei isolations, so peaks represent regions of open
chromatin *in vivo*. A recent comprehensive review of outlines the molecular biology details of DNaseseq
(Vierstra and Stamatoyannopoulos, 2016). It is noteworthy that DNase does not cleave double
stranded DNA, instead DNase nicks the phosphodiester backbone of DNA and four nicking events are
needed to detect a DNA fragment by DNase-seq, although only two nicking events are detected per
fragment (Vierstra and Stamatoyannopoulos, 2016; Thomas, 1956).

```bash
mkdir ~/DNase_ENCODE
cd ~/DNase_ENCODE
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseMcf7Est100nm1hRawDataRep1.fastq.tgz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseMcf7Est100nm1hRawDataRep2.fastq.tgz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseMcf7Estctrl0hRawDataRep1.fastq.tgz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseMcf7Estctrl0hRawDataRep2.fastq.tgz
tar -xvf wgEncodeUwDnaseMcf7Estctrl0hRawDataRep2.fastq.tgz
tar -xvf wgEncodeUwDnaseMcf7Estctrl0hRawDataRep1.fastq.tgz
tar -xvf wgEncodeUwDnaseMcf7Est100nm1hRawDataRep2.fastq.tgz
tar -xvf wgEncodeUwDnaseMcf7Est100nm1hRawDataRep1.fastq.tgz
cat UwStam_MCF7-*fastq > UW_MCF7_both.fastq
gzip UW_MCF7_both.fastq
rm *.fastq.tgz
rm *fastq
```

Download short read archive data set (SRA accession SRX247626) of DNase-seq data from DNA purfied
from IMR90 cells (Lazarovici *et al.*, 2013), convert the *sra* to a *fastq* file using `fastq-dump` (herein we use
version: 2.7.0), change the name to be descriptive, and compress the file. Note that this is DNase-seq
data from naked DNA digestion (i.e. no bound proteins), which provided the most compeling evidence
that DNase signatures at the site of transcription factor (TF) binding are not a result of protein binding
(He *et al.*, 2014).

```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR769/SRR769954/SRR769954.sra
fastq-dump SRR769954.sra
mv SRR769954.fastq IMR90_Naked_DNase.fastq
gzip IMR90_Naked_DNase.fastq
rm *sra
```

## Index the Appropriate Genome File
Retrieve the relevant genome from UCSC (Karolchik *et al.*, 2014), we will use the latest assembly, hg38.
This is a zipped *fasta* file of the entire human genome. Bowtie 2 is an efficient tool for aligning sequencing
reads to long reference sequences. For this execution we used Bowtie2 version 2.2.6. The first task
is to build the genome index with [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) (Langmead *et al.*, 2009). This only has to be performed once per
genome.
```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bowtie2-build hg38.fa hg38
```

## Align fastq.gz files with bowtie2
This code will loop through files in a directory. The file name is split on the ’.fastq.gz’ string and the
variable ’name’ is assigned to the first string.
```bash
for fq in *.fastq.gz
do
    name=$(echo $fq | awk -F".fastq.gz" '{print $1}')
    echo $name
    bowtie2 -x hg38 -U $fq -S $name.sam
done
```

## Convert to bam file format
Sam files retain all the information from the *fastq* files, but include additional information, including
alignment coordinates ([http://samtools.github.io/hts-specs/SAMv1.pdf](http://samtools.github.io/hts-specs/SAMv1.pdf)). The header has all the
chromosome size information from the hg38.fa file.
Next convert the sam file to the compressed and **sorted** BAM format using `samtools` (version 1.2 used
herein) (Li *et al.*, 2009).
```bash
for sam in *.sam
do
    name=$(echo $sam | awk -F".sam" '{print $1}')
    echo $name
    samtools view -b $sam | samtools sort - $name
done
rm *sam
```