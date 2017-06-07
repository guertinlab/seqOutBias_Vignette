---
layout: docs
title: Correction of MNase sequence bias from MNase-seq data
prev: part3
next: part5
number: 4
---

# Correction of MNase sequence bias from MNase-seq data

We use the same process to correct MNase-seq data.

## Processing MNase-seq data with seqOutBias
The MNase-seq data is paired-end. Note the use of the pdist=100:400 to specifically process reads
that have insert sizes between 100 and 400 base pairs.

```bash
#completely new
mkdir MNase_Zhang
cd MNase_Zhang
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX564/SRX564203/SRR1323041/SRR1323041.sra 
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX564/SRX564204/SRR1323042/SRR1323042.sra 
fastq-dump --split-3 SRR1323041.sra
fastq-dump --split-3 SRR1323042.sra
mv SRR1323041_1.fastq MNase_MCF7_PE1_rep1.fastq
mv SRR1323042_1.fastq MNase_MCF7_PE1_rep2.fastq
mv SRR1323041_2.fastq MNase_MCF7_PE2_rep1.fastq
mv SRR1323042_2.fastq MNase_MCF7_PE2_rep2.fastq
gzip *fastq

for fq in *_PE1_rep1.fastq.gz
do
    name=$(echo $fq | awk -F"_PE1_rep1.fastq.gz" '{print $1}')
    echo $name
    bowtie2 -x ~/DNase_ENCODE/hg38 -1 $fq,${name}_PE1_rep2.fastq.gz -2 ${name}_PE2_rep1.fastq.gz,${name}_PE2_rep2.fastq.gz -S $name.sam
    samtools view -b $name.sam | samtools sort - $name
    seqOutBias ~/DNase_ENCODE/hg38.fa $name.bam --no-scale --bw=${name}_0-mer.bigWig --shift-counts --read-size=101 --pdist=100:400
    seqOutBias ~/DNase_ENCODE/hg38.fa $name.bam --kmer-mask=NNNNCNNNN --bw=${name}_8-mer.bigWig --shift-counts --read-size=101 --pdist=100:400
done

mkdir MNase_final
mv MNase_MCF7_0-mer.bigWig MNase_0-mer.bigWig
mv MNase_MCF7_8-mer.bigWig MNase_8-mer.bigWig
mv MNase_0-mer.bigWig MNase_final
mv MNase_8-mer.bigWig MNase_final
```

## Plotting MNase-seq composites using `R`
Plot the composite MNase profile using the MCF7 ChIP-seq peaks from Section 2.4. These MNase-seq data are relatively low coverage and MNase-seq reads are not enriched at TF binding sites, as in DNase-seq. Therefore, the sequence bias correction is more apparent when you average over all the motif instances in the genome that we identified by FIMO in Section 2.4.

```r
source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')

all.composites.mnase.mcf7 = cycle.fimo.new.not.hotspots(path.dir.mast = '~/DNase_ENCODE/',
    path.dir.bigWig = '/Users/guertinlab/MNase_Zhang/MNase_final', window = 30, exp = 'MCF7_MNase')
all.composites.mnase = cycle.fimo.new.not.hotspots(path.dir.fimo = '~/DNase_ENCODE/',
    path.dir.bigWig = '/Users/guertinlab/MNase_Zhang/MNase_final', window = 30, exp = 'MNase')

save(all.composites.mnase.mcf7, file = '~/MNase_Zhang/all.composites.mnase.mcf7.Rdata')
save(all.composites.mnase, file = '~/MNase_Zhang/all.composites.mnase.Rdata')
```
