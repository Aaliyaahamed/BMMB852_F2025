## Week 13 Assignment: RNA-Seq Analysis
## RNA-Seq with Hisat2
## Aaliya • BMMB852 • 2025-12-09

Overview\
This week the assignment is to run an RNA-Seq analysis including the following steps:

Download the genome and sequence data.

1. Align/Classify reads relative to the genome

2. Quantify the expression levels

3. Compare the expression levels to find up- and down-regulated genes

4. Visualize and interpret the results

For this analysis, I utilized the Human Brain Reference (HBR) and Universal Human Reference (UHR) dataset. The analysis was restricted to Chromosome 22 to demonstrate the alignment and quantification pipeline. The goal is to replicate the biological expectation that brain-specific genes should be upregulated in the HBR samples compared to the UHR samples.

## Download the genome and sequence data
Source Data: Biostar Handbook (UHR-HBR subset).

Unlike the full SRA download, I utilized a curated subset of the data specific to Chromosome 22 to streamline the processing time. The data includes 3 replicates for the Control (HBR) and 3 replicates for the Treatment (UHR).

To download and unpack the data, I used the following commands:

```
wget -nc http://data.biostarhandbook.com/data/uhr-hbr.tar.gz
tar xzvf uhr-hbr.tar.gz
Reference and Annotations:
```
Genome: Homo sapiens Chromosome 22 (chr22.genome.fa)

Annotation: GRCh38 GTF (chr22.gtf)

## Align/Classify reads relative to the genome
To automate the alignment process, I wrote a Makefile that handles indexing, alignment, and sorting. I used Hisat2 for splice-aware alignment and Samtools for BAM processing.

Index the genome:

Makefile
```
hisat2-build refs/chr22.genome.fa refs/hisat2_index
```
Align samples: The alignment was performed for all 6 samples (HBR_1-3, UHR_1-3). The specific command used within the Makefile was:

Makefile
```
hisat2 -x refs/hisat2_index -U reads/HBR_1_R1.fq | samtools sort > bam/HBR_1.bam
```
IGV Verification of RNA-Seq Data: I generated BigWig files to visualize the alignments in IGV. As seen in the screenshots below, the data is confirmed to be RNA-Seq data because the read coverage is discontinuous; reads align to exons and "skip" introns, leaving gaps that match the gene models in the annotation track.

## Quantify the expression levels
I used featureCounts to quantify the number of reads mapping to each gene feature.

Command:
```
featureCounts -a refs/chr22.gtf -o counts.txt bam/HBR_1.bam bam/HBR_2.bam bam/HBR_3.bam bam/UHR_1.bam bam/UHR_2.bam bam/UHR_3.bam
```
Results: Unlike the classmate's example where few genes were found, my quantification was successful. The counts.txt matrix contains robust counts for thousands of genes.

Top expressed genes in HBR Replicate 1:

```
ENSG00000100201.18    9995
ENSG00000241973.10    11407
ENSG00000100321.14    762
```
## Compare the expression levels
To verify differential expression, I manually inspected the count matrix for biologically relevant genes. I selected SYN3 (Synapsin III), a gene known to be involved in neurotransmitter release.

Observed Counts (HBR vs UHR):

HBR (Brain) Samples: Consistently high expression (~700–1000 reads).

HBR_1: 762 | HBR_2: 996 | HBR_3: 865

UHR (Universal) Samples: Negligible expression (~50 reads).

UHR_1: 54 | UHR_2: 51 | UHR_3: 42

This manual inspection confirms that the SYN3 gene is significantly upregulated in the Brain Reference samples compared to the Universal Reference.

## Visualize and interpret the results
![alt text](image-1.png)
I visualized the SYN3 locus (ENSG00000100321) in IGV to confirm the values seen in the count matrix.

Top 3 Tracks (HBR): These tracks show dense, tall peaks of read coverage over the exons.

Bottom 3 Tracks (UHR): These tracks are nearly flat, corresponding to the low counts observed in the matrix.

The visual evidence in IGV perfectly aligns with the count matrix, validating the pipeline.

## Discussion
I was able to successfully replicate the expected biological signal for this dataset. By aligning the UHR and HBR samples to Chromosome 22, I generated a count matrix that accurately reflects tissue-specific gene expression.

In contrast to potential issues with annotation mismatch or over-downsampling that can occur in SRA projects, the use of the matched chr22.genome.fa and chr22.gtf allowed featureCounts to successfully assign reads to genes. The analysis of SYN3 served as a positive control, confirming that the pipeline correctly identified upregulation in the brain tissue samples. The splicing patterns observed in IGV further confirm the high quality of the RNA-Seq library preparation and alignment.