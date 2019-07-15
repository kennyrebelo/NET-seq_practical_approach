# Supplemental Material and Methods

In this section we provide the specific options used for each program to aid the reader setup their own analysis.

**Dataset used**: Long reads S5P mNET-seq dataset (Nojima et al., 2018). From [GSE106881](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106881) we used [GSM2856674](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2856674) and [GSM2856677](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2856677).

## S2.1. Quality control and Adapter trimming
Initial state of the data was verified using FastQC (version 0.11.7)
FastQC command:
```
fastqc mNET_Long_S5P_rep1_1.fastq
fastqc mNET_Long_S5P_rep1_2.fastq
```

Cutadapt (version 1.18) was used with an error rate of 0.05, and allow it to match ‘N’s in the reads to the adapter sequence; discard the reads that are shorter than 10 bases, and remove the adapter only once from each read, as described in (Nojima et al., 2016).
Cutadapt command:
```
cutadapt -a TGGAATTCTCGGGTGCCAAGG -A GATCGTCGGACTGTAGAACTCTGAAC -m 10 -e 0.05 --match-read-wildcards -n 1 -o mNET_Long_S5P_rep1_1_tr.fastq.gz -p mNET_Long_S5P_rep1_2_tr.fastq.gz mNET_Long_S5P_rep1_1.fastq mNET_Long_S5P_rep1_2.fastq
```


## S2.2. Mapping of reads to the reference genome 
STAR (version 2.6.0b) set to detect chimeric alignments with the minimum mapped length of at least 20nt on each end.
STAR index generation:
```
STAR --runMode genomeGenerate --genomeDir ./starIndex/ --genomeFastaFiles /genomes/human/hg38/GRCh38.primary.genome.fa –-sjdbGTFfile /genomes/human/hg38/gencode.v28.annotation.gtf
```

STAR command (paired reads):
```
STAR --runMode alignReads --genomeDir /genomes/human/hg38/star/ --readFilesIn ./mNET_Long_S5P_rep1_1_tr.fastq.gz ./mNET_Long_S5P_rep1_2_tr.fastq.gz --chimSegmentMin 20 --outSAMtype BAM Unsorted --readFilesCommand gunzip -c --outFileNamePrefix /alignments/mNET_Long_S5P_rep1_
```
To obtain the uniquely mapped reads using SAMtools (version 1.7):
```
samtools view -H mNET_Long_S5P_rep1_Aligned.out.bam > mNET_Long_S5P_rep1_header.sam
samtools view -q 255 mNET_Long_S5P_rep1_Aligned.out.bam > mNET_Long_S5P_rep1_unique.sam
cat mNET_Long_S5P_rep1_header.sam mNET_Long_S5P_rep1_unique.sam > mNET_Long_S5P_rep1_unique_H.sam
samtools view -Sb -h mNET_Long_S5P_rep1_unique_H.sam > mNET_Long_S5P_rep1_unique.bam
rm -f mNET_Long_S5P_rep1_unique_H.sam
rm -f mNET_Long_S5P_rep1_unique.sam
```

## S2.3. Identification of RNA 3' ends
Using python (version 2.7.12) script [get_SNR_bam_ignoreSoftClip.py](https://github.com/kennyrebelo/mNET_snr)
command used:
```
python get_SNR_bam_ignoreSoftClip.py -f mNET_Long_S5P_rep1_unique.bam -s mNET_Long_S5P_rep1 -d ./
```


## S2.4. Identification and removal of PCR internal priming 
Python script [Filter_InternalPriming.py](https://github.com/kennyrebelo/Filtering_InternalPriming)

This script removes reads that result from internal priming events. Reads that result from internal priming will have their downstream
nucleotides (after alignment) matching partially/fully to the adapter sequence that was connected to the 3'OH end (TGGAATTCTCGGGTGCCAAGG) during library preparation.
This script checks for the presence or absence of the 3'OH connected adapter sequence downstream of the aligned reads.

command used for paired reads:
```
python Filter_InternalPriming.py -f /alignments/mNET_Long_S5P_rep1_unique_sorted.bam -s paired -a TGG.. -g /genomes/human/hg38/GRCh38.primary.genome.fa
```


## S2.5 Distinguishing nascent RNA from splicing intermediates
NET-seq signal from co-transcriptional splicing corresponds to reads that map to the last nucleotide of exons/introns. If there is interest in co-transcriptional splicing then only use these genomic regions. If the interest is in nascent RNA then you should disregard these regions from the analysis.
Some example commands on how to exclude reads that map to the very last nucleotide of exons/introns:

```
intersectBed -a mNET_Long_S5P_rep1_SNR.bam -b exons_lastNT.bed -wa -v | samtools view - > mNET_Long_S5P_rep1_SNR_noLastNT_temp.sam
cat mNET_Long_S5P_rep1_header.sam mNET_Long_S5P_rep1_SNR_noLastNT_temp.sam > mNET_Long_S5P_rep1_SNR_noLastNT.sam
samtools view -bS mNET_Long_S5P_rep1_SNR_noLastNT.sam > mNET_Long_S5P_rep1_SNR_noLastNT.bam
rm -f mNET_Long_S5P_rep1_SNR_noLastNT_temp.sam
rm -f mNET_Long_S5P_rep1_SNR_noLastNT.sam
```

## S2.6. Selection of transcriptionally active genes
Read counts values in each 500kb window were obtained using [ ](https://bedtools.readthedocs.io/en/latest/index.html) (version  v2.27.1-1-gb87c465)

command used:
```
coverageBed -a intergenic_regions_500kb_Windows.bed -b Filter_IP/mNET_Long_S5P_rep1_noInternalPriming.bam -counts > intergenic_regions_500kb_Windows_cov.bed
```
