---
layout: docs
title: Scaling Tissue Accessible Chromatin (TACh) Benzonase and Cyanase Digested DNA
prev: part4
next: part6
number: 5
---

# Scaling Tissue Accessible Chromatin (TACh) Benzonase and Cyanase Digested DNA

Many enzymes nick DNA with distinct sequence biases. Here we characterize the biases of Benzonase and Cyanase and show that seqOutBias can scale DNA digestion data resulting from Benzonase and Cyanase treatment (Grøntved *et al.*, 2012).

## Retrieving TACh-seq data from mouse liver

```bash
mkdir ~/TACh_Grontved
cd ~/TACh_Grontved 
url=ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX174/ 
wget ${url}SRX174756/SRR535737/SRR535737.sra
wget ${url}SRX174757/SRR535738/SRR535738.sra
wget ${url}SRX174757/SRR535739/SRR535739.sra
wget ${url}SRX174758/SRR535740/SRR535740.sra
wget ${url}SRX174761/SRR535744/SRR535744.sra
wget ${url}SRX174760/SRR535742/SRR535742.sra
wget ${url}SRX174760/SRR535743/SRR535743.sra
wget ${url}SRX174759/SRR535741/SRR535741.sra
wget ${url}SRX174755/SRR535735/SRR535735.sra
wget ${url}SRX174755/SRR535736/SRR535736.sra

for i in *sra
do
    fastq-dump $i
done

mv SRR535737.fastq mm10_liver_Benzonase0.25U.fastq
mv SRR535738.fastq mm10_liver_Benzonase1U_1.fastq
mv SRR535739.fastq mm10_liver_Benzonase1U_2.fastq
mv SRR535740.fastq mm10_liver_Benzonase4U.fastq
mv SRR535741.fastq mm10_liver_Cyanase0.25U.fastq
mv SRR535742.fastq mm10_liver_Cyanase1U_1.fastq
mv SRR535743.fastq mm10_liver_Cyanase1U_2.fastq
mv SRR535744.fastq mm10_liver_Cyanase4U.fastq
mv SRR535735.fastq DNaseI_a.fastq 
mv SRR535736.fastq DNaseI_b.fastq
cat *Benz* > mm10_liver_Benzonase.fastq 
cat *Cyan* > mm10_liver_Cyanase.fastq
cat DNaseI_*.fastq > mm10_liver_DNase.fastq 
gzip *ase.fastq
rm *fastq rm *sra
```

## Index the mm10 genome file and align to mm10

Retrieve the compressed mm10 genome from UCSC (Karolchik *et al.*, 2014). This is a *2bit* compressed file and needs to be converted to a *fasta* using `twoBitToFa` from [http://hgdownload.soe.ucsc.edu/admin/exe/](http://hgdownload.soe.ucsc.edu/admin/exe/).

```bash
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit
twoBitToFa mm10.2bit mm10.fa
bowtie2-build mm10.fa mm10

for fq in *.fastq.gz
do
  name=$(echo $fq | awk -F".fastq.gz" '{print $1}')
  echo $name
  bowtie2 -x mm10 -U $fq -S $name.sam
  samtools view -b $name.sam | samtools sort - $name
  rm $name.sam
done
```

## Using seqOutBias to determine the sequence preference for Cyanase and Benzonase

```bash
for bam in mm10_liver*.bam
do
    name=$(echo $bam | awk -F"mm10_liver_" '{print $NF}' | awk -F".bam" '{print $1}')
    echo $name
    seqOutBias mm10.fa $bam --no-scale --bw=${name}_0-mer.bigWig --shift-counts --skip-bed --read-size=35 
    seqOutBias mm10.fa $bam --kmer-mask=NNNCNNN --bw=${name}_6-mer.bigWig --shift-counts --read-size=35 
    seqOutBias mm10.fa $bam --kmer-mask=NNNNCNNNN --bw=${name}_8-mer.bigWig --shift-counts --read-size=35
done

mv DNase_0-mer.bigWig DNase_mouse_0-mer.bigWig 
mv DNase_6-mer.bigWig DNase_mouse_6-mer.bigWig 
mv DNase_8-mer.bigWig DNase_mouse_8-mer.bigWig

mkdir dnase
mkdir benzonase
mkdir cyanase
mv Benzonase*bigWig benzonase 
mv Cyanase*bigWig cyanase
mv DNase*bigWig dnase
```

## Retrieving ChIP-seq binding and sequence motif data for mouse liver

To look at composite footprints that are the result of transcription factors binding to DNA in the context of chromatin, we need to first find all the regions bound by the factor. We will retrieve TF binding data from several sources (Seo *et al.*, 2009; Stamatoyannopoulos *et al.*, 2012; Grøntved *et al.*, 2013). Then we convert the peak files to the latest genome assembly.

```bash
#CTCF
wget https://www.encodeproject.org/files/ENCFF001YAM/@@download/ENCFF001YAM.bed.gz
mv ENCFF001YAM.bed.gz CTCF.mm9.bed.gz
#FOXA2
url=ftp://ftp.ncbi.nlm.nih.gov/geo/series/
wget ${url}GSE25nnn/GSE25836/suppl/GSE25836_Mouse_Liver_FOXA2_GLITR_1p5_FDR.bed.gz
mv GSE25836_Mouse_Liver_FOXA2_GLITR_1p5_FDR.bed.gz FOXA2.mm8.bed.gz
#CEBP-beta
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46047/suppl/GSE46047%5FCEBPb%5Fpeaks%5Fveh%5Fmouse%5Fliver%5Fmm9%2Etxt%2Egz 
mv GSE46047_CEBPb_peaks_veh_mouse_liver_mm9.txt.gz CEBP-beta_temp.mm9.bed.gz

wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz 
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm8/liftOver/mm8ToMm10.over.chain.gz

gunzip *.over.chain.gz
gunzip *mm*bed.gz
tail +2 CEBP-beta_temp.mm9.bed > CEBP-beta.mm9.bed 
rm CEBP-beta_temp.mm9.bed

for peak in *mm8.bed
do
    name=$(echo $peak | awk -F".mm8" '{print $1}')
    echo $name
    liftOver $peak mm8ToMm10.over.chain $name.mm10.bed $name.mm10.unmapped.txt -bedPlus=6 
    fastaFromBed -fi mm10.fa -bed $name.mm10.bed -fo $name.mm10.fasta
done

for peak in *mm9.bed
do
    name=$(echo $peak | awk -F".mm9" '{print $1}')
    echo $name
    liftOver $peak mm9ToMm10.over.chain $name.mm10.bed $name.mm10.unmapped.txt -bedPlus=6 
    fastaFromBed -fi mm10.fa -bed $name.mm10.bed -fo $name.mm10.fasta
done

head -9 ~/DNase_ENCODE/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme > header_meme_temp.txt

grep -i -A 16 'MOTIF MA0047.2 FOXA2' ~/DNase_ENCODE/motif_databases/JASPAR/JASPAR_CORE_2016.meme > foxa2_temp.txt 
grep -i -A 23 'MOTIF MA0139.1 CTCF' ~/DNase_ENCODE/motif_databases/JASPAR/JASPAR_CORE_2016.meme > ctcf_temp.txt 
grep -i -A 15 'MOTIF MA0466.2 CEBPB' ~/DNase_ENCODE/motif_databases/JASPAR/JASPAR_CORE_2016.meme > cebpb_temp.txt 
cat header_meme_temp.txt foxa2_temp.txt > FOXA2_minimal_meme.txt
cat header_meme_temp.txt ctcf_temp.txt > CTCF_minimal_meme.txt
cat header_meme_temp.txt cebpb_temp.txt > CEBP-beta_minimal_meme.txt rm *temp.txt

for meme in *.mm10.fasta
do
    name=$(echo $meme | awk -F".mm10.fasta" '{print $1}')
    echo $name
    mast ${name}_minimal_meme.txt $meme -hit_list -mt 0.0001 > ${name}_mast.txt 
    fimo --thresh 0.0001 --text ${name}_minimal_meme.txt mm10.fa > ${name}_fimo.txt ceqlogo -i1 ${name}_minimal_meme.txt -o ${name}_logo.eps -N -Y
done
```

## Plotting Benzonase and Cyanase composites using `R`

```r
source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')

#note that the full path is needed to the directory containing the bigWigs
all.composites.cyanase = cycle.fimo.new.not.hotspots(path.dir.mast = '~/TACh_Grontved/', path.dir.bigWig = '/Users/guertinlab/TACh_Grontved/cyanase/', window = 30, exp = 'Cyanase')
all.composites.benzonase = cycle.fimo.new.not.hotspots(path.dir.mast = '~/TACh_Grontved/', path.dir.bigWig = '/Users/guertinlab/TACh_Grontved/benzonase/', window = 30, exp = 'Benzonase')
all.composites.dnase = cycle.fimo.new.not.hotspots(path.dir.mast = '~/TACh_Grontved/', path.dir.bigWig = '/Users/guertinlab/TACh_Grontved/dnase/', window = 30, exp = 'DNase_mm10')

save(all.composites.cyanase, all.composites.benzonase, all.composites.dnase, file = 'MOUSE_composites.Rdata')

composites.func.panels.naked.chromatin(all.composites.benzonase[all.composites.benzonase$cond == 'Benzonase_0-mer' | all.composites.benzonase$cond == 'Benzonase_8-mer',], fact= "Benzonase8", summit= "Motif",num = 24,
    col.lines = c(rgb(0,0,1,1/2), rgb(0,0,0,1/2)), fill.poly = c(rgb(0,0,1,1/4), rgb(0,0,0,1/4)))

composites.func.panels.naked.chromatin(all.composites.cyanase[all.composites.cyanase$cond == 'Cyanase_0-mer' | all.composites.cyanase$cond == 'Cyanase_8-mer',], fact= "Cyanase8", summit= "Motif",num = 24,
    col.lines = c(rgb(0,0,1,1/2), rgb(0,0,0,1/2)), fill.poly = c(rgb(0,0,1,1/4), rgb(0,0,0,1/4)))

composites.func.panels.naked.chromatin(all.composites.dnase[all.composites.dnase$cond == 'DNase_mouse_0-mer' | all.composites.dnase$cond == 'DNase_mouse_6-mer',], fact= "Dnase6", summit= "Motif",num = 24,
    col.lines = c(rgb(0,0,1,1/2), rgb(0,0,0,1/2)), fill.poly = c(rgb(0,0,1,1/4), rgb(0,0,0,1/4)))
```
