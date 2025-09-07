# Week 2 Assignment: Demonstrate data analysis at Unix command line
#### Aaliya Ahamed • BMMB852 • 2025-09-06

---
## 1) Working directory
Since, I haver to replicate this analysis, I am defining the path to my wokring directory. The name of my working directory is "ASSIGNMENT_2"

```bash
DIR=YOURPATH/02_Assignment2
mkdir -p $DIR
cd $DIR

```

## 2) GFF3 file from Ensemble
I am downloading the GGF3 file of the organism.\
(*Mola mola* is the world’s heaviest bony fish — [Wikipedia](https://en.wikipedia.org/wiki/Ocean_sunfish))

The assembly used here is ASM169857v1.\
Use `gunzip` to unzip the file.

```bash
GFF=Mola_mola.ASM169857v1.115.gff3
GFF_LINK=http://ftp.ensembl.org/pub/current_gff3/mola_mola/Mola_mola.ASM169857v1.115.gff3.gz
wget $GFF_LINK
gunzip $GFF.gz
```

## 3) Exploring the organism
To look at the content of the GFF file:
```bash
head $GFF
tail $GFF
```
The head confirms this is a valid GFF3 file with Ensembl metadata and sequence region declarations.

The tail shows how scaffold-level annotations are recorded, ending with small scaffolds.

## 4) Extracting sequence regions
```bash
grep "sequence-region" $GFF > sequence_regions.txt
``` 
The file lists many scaffolds (e.g., KV751321.1, KV751319.1, etc.), which indicates the Mola mola genome is scaffold-level rather than fully chromosome-level.

## 5) Count and Rank features
Count all features (excluding comment lines)
``` bash
grep -v "^#" $GFF | wc -l > all_features_count.txt
```
Result: 619,591 features

Rank feature types:
```bash
cut -f3 $GFF | grep -v "^#" | sort | uniq -c | sort -r > all_features_ranked.txt
```
Top 10 features:
```bash
head all_features_ranked.txt > chr_features_top10.txt
```
Top 10 annotated features:

- exon-270,282

- CDS-269,689

- mRNA-29,008

- biological_region-22,671

- gene-21,404

- region-5,552

- ncRNA_gene-444

- snoRNA-202

- snRNA-86

- miRNA-78

#### These results show that exons and CDS dominate the annotation, as expected for a protein-coding genome.
## 6) Total number of genes
```bash
grep -v "^#" $GFF | awk -F"\t" '$3=="gene"' | wc -l > gene_count.txt
```
Result: 21,404 genes

## 7) Feature Type?
only 2 lnc_RNA annotations were identified, suggesting these RNAs are annotated but rare in this assembly.

## 8) Confirming results

According to my Ensembl GFF3 analysis, the *Mola mola* genome (ASM169857v1) contains **21,404 genes**.

When I checked the NCBI Assembly page for *Mola mola* ([NCBI ASM169857v1](https://www.ncbi.nlm.nih.gov/assembly/169857)), I found:

- Genome size: **639.5 Mb**
- Number of scaffolds: **5,552**
- Number of contigs: **51,826**
- GC content: **41%**
- Assembly level: **Scaffold**
- Sample: adult female *Mola mola* (blood tissue)

NCBI explicitly notes:  
> *“This scaffold-level genome assembly includes 5,552 scaffolds and no assembled chromosomes.”*

### This confirms my observation from the GFF file that the genome is scaffold-based, not chromosome-based.  
---

### 9) Extra analysis
- Average exons per gene (from Ensembl GFF3): ~12.6 exons per gene, consistent with expectations for fish genomes.  
- Very few lncRNAs (2 annotations) are present, suggesting non-coding RNAs are under-annotated in this genome.  

Together, Ensembl and NCBI confirm that while *Mola mola* has a robust gene-level annotation, its genome assembly remains incomplete at the chromosome level.


## 10) My thoughts
The annotation appears fairly complete at the gene level with ~21,400 genes and nearly equal numbers of exons and CDS (~270k each). However, the genome is assembled at the scaffold-level (many KV-numbered scaffolds instead of chromosomes). Overall, the annotation quality is strong, but the incomplete assembly reduces its utility for higher-resolution genomic studies.