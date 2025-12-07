# Makefile for Assignment 13: Transcriptomic Profiling & Matrix Generation

SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# --- Global Variables ---
REF_BUILD=chr22.genome
REF_FA=refs/${REF_BUILD}.fa
REF_GTF=refs/chr22.gtf

# Input/Output Directories
READS_DIR=reads
BAM_DIR=bam
RES_DIR=results

# Sample Identifier 
SAMPLE=HBR_1_R1
FASTQ=${READS_DIR}/${SAMPLE}.fq

# --- Targets ---

# 1. Align reads to the reference 
# Usage: make align SAMPLE=HBR_1_R1
align:
	mkdir -p ${BAM_DIR}
	# Aligning with BWA MEM and sorting
	bwa mem -t 4 ${REF_FA} ${FASTQ} | samtools sort > ${BAM_DIR}/${SAMPLE}.bam
	# Indexing BAM for IGV visualization
	samtools index ${BAM_DIR}/${SAMPLE}.bam

# 2. Generate BigWig for Visualization
# Usage: make bigwig SAMPLE=HBR_1_R1
bigwig:
	# Create intermediate BedGraph
	bedtools genomecov -ibam ${BAM_DIR}/${SAMPLE}.bam -split -bg | sort -k1,1 -k2,2n > ${BAM_DIR}/${SAMPLE}.bedgraph
	# Convert to BigWig using reference index
	bedGraphToBigWig ${BAM_DIR}/${SAMPLE}.bedgraph ${REF_FA}.fai ${BAM_DIR}/${SAMPLE}.bw
	# Cleanup intermediate file
	rm ${BAM_DIR}/${SAMPLE}.bedgraph

# 3. Feature Counting
# Generates raw counts for downstream Differential Expression or Multi-omics integration
count_features:
	mkdir -p ${RES_DIR}
	featureCounts -a ${REF_GTF} -o ${RES_DIR}/counts.txt \
		${BAM_DIR}/HBR_1_R1.bam ${BAM_DIR}/HBR_2_R1.bam ${BAM_DIR}/HBR_3_R1.bam \
		${BAM_DIR}/UHR_1_R1.bam ${BAM_DIR}/UHR_2_R1.bam ${BAM_DIR}/UHR_3_R1.bam

# 4. Matrix formatting
matrix:
	# Convert counts to CSV
	Rscript ~/src/r/format_featurecounts.r -c ${RES_DIR}/counts.txt -o ${RES_DIR}/counts.csv
	# Fetch gene metadata (using Human Ensembl mapping)
	Rscript ~/src/r/create_tx2gene.r -d hsapiens_gene_ensembl
	# Merge metadata with counts
	Rscript ~/src/r/format_featurecounts.r -c ${RES_DIR}/counts.txt -t tx2gene.csv -o ${RES_DIR}/final_matrix.csv

# Setup target to fetch data 
setup:
	wget -nc http://data.biostarhandbook.com/data/uhr-hbr.tar.gz
	tar xzvf uhr-hbr.tar.gz
	# Index the provided reference immediately
	bwa index ${REF_FA}
	samtools faidx ${REF_FA}

.PHONY: align bigwig count_features matrix setup