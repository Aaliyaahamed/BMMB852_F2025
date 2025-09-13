# Week 3 Assignment: Visualizing genomic data
#### Aaliya Ahamed • BMMB852 • 2025-09-12
------------------------------------

### 1. Navigate to working directory
``` bash
cd /mnt/d/BMMB852_F2025
mkdir 03_Assignment3
cd 03_Assignment3
```
### 2. Download genome (FASTA) and annotation (GFF3) for Aquifex aeolicus from Ensembl Bacteria FTP
>For this assignment, I worked with the bacterium **Aquifex aeolicus**, using its reference genome and annotation files downloaded from the Ensembl Bacteria FTP database. 


``` bash
# Genome FASTA
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/current/fasta/bacteria_0_collection/aquifex_aeolicus_vf5_gca_000008625/dna/Aquifex_aeolicus_vf5_gca_000008625.ASM862v1.dna.toplevel.fa.gz

# Annotation GFF3
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/current/gff3/bacteria_0_collection/aquifex_aeolicus_vf5_gca_000008625/Aquifex_aeolicus_vf5_gca_000008625.ASM862v1.62.gff3.gz

```
### 3. Unzip
``` bash
gunzip Aquifex_aeolicus_vf5_gca_000008625.ASM862v1.dna.toplevel.fa.gz
gunzip Aquifex_aeolicus_vf5_gca_000008625.ASM862v1.62.gff3.gz
```
### 4. Renaming files
``` bash
mv Aquifex_aeolicus_vf5_gca_000008625.ASM862v1.dna.toplevel.fa aquifex.fa
mv Aquifex_aeolicus_vf5_gca_000008625.ASM862v1.62.gff3 aquifex.gff
```

### 5. Check genome statistics
``` bash
seqkit stats aquifex.fa
```
Output:
type → nucleotide sequences (DNA).

num_seqs → there are 2 sequences (most likely the main chromosome and an additional contig or plasmid).

sum_len → the total genome length = 1,590,791 bp (~1.59 Mb).

min_len → the smallest sequence length = 39,456 bp.

avg_len → the average sequence length = 795,395 bp.

max_len → the largest sequence length = 1,551,335 bp.
> The Aquifex aeolicus genome is about 1.59 Mb in size, represented by 2 sequences (a main chromosome and a smaller contig)

### 6. Count features in GFF file
``` bash
grep -v "#" aquifex.gff | cut -f3 | sort | uniq -c
```
Output:
Output:
1552 CDS
1553 genes
1553 mRNA
1604 exons
44 tRNA
6 rRNA
51 ncRNA genes
1 tmRNA
1 chromosome, 1 region

### 7. Extract only genes and transcripts to a new GFF
``` bash
grep -v "#" aquifex.gff | awk '$3=="gene" || $3=="transcript"' > aquifex_genes_transcripts.gff
```
 
 ### 8. samtools faidx aquifex.fa
This step was recommended by perplexity bot as IGV was not permitting direct opening of .fa file 

### 9. Extract CDS only
This step will help in locating start and stop codon in IGV
```bash 
grep -v "#" aquifex.gff | awk '$3=="CDS"' > aquifex_cds.gff

# *Check first and last CDS*
head -n 1 aquifex_cds.gff
tail -n 1 aquifex_cds.gff
```
#### In IGV, the coding sequence for gene aq_aa39 ends with a stop codon, visualized as a red * in the translation track. This confirms that the gene terminates properly.we zoomed in on the gene fusA (aq_001). The first CDS begins with the canonical start codon ATG (translated to methionine, shown as green "M"). The coding sequence terminates with a stop codon (TAA/TAG/TGA, shown as a red "*" in IGV). This confirms that the gene annotations are consistent with expected biological rules of translation.
> Screenshots are included in the images/ folder