# Makefile for a reproducible alignment pipeline.
# This Makefile is designed to be called by GNU Parallel,
# which provides the SRR, SAMPLE, GENOME_ACC, and GENOME_NAME variables.

SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# ===== Default Parameters (can be overridden from command line) =====
GENOME_ACC    ?= U00096.3
GENOME_NAME   ?= ecoli_k12
SRR    ?= SRR1972739
SAMPLE ?= sample_default

# ===== Directories =====
GENOME_DIR  := genome
READS_DIR   := reads
RESULTS_DIR := results
TMP_DIR     := .tmp

# ===== File Paths =====
# Genome files are named with the user-friendly GENOME_NAME
GENOME_FA  := $(GENOME_DIR)/$(GENOME_NAME).fna
GENOME_FAI := $(GENOME_FA).fai
GENOME_BWA_INDEX := $(GENOME_FA).bwt

# Read files are based on SRR
R1 := $(READS_DIR)/$(SRR)_1.fastq.gz
R2 := $(READS_DIR)/$(SRR)_2.fastq.gz
RS := $(READS_DIR)/$(SRR).fastq.gz

# Result files are based on the SAMPLE name from design.csv
BAM      := $(RESULTS_DIR)/$(SAMPLE).sorted.bam
BAI      := $(BAM).bai
BW       := $(RESULTS_DIR)/$(SAMPLE).bw
STATS    := $(RESULTS_DIR)/$(SAMPLE).stats.txt
FASTQC_DIR := $(RESULTS_DIR)/$(SAMPLE)_fastqc
FASTQC_OK  := $(RESULTS_DIR)/.$(SAMPLE).fastqc_ok


# ===== Main Pipeline Targets =====
.PHONY: all genome reads fastqc index align stats bigwig clean

all: align bigwig stats

genome: $(GENOME_FAI)
reads: $(READS_DIR)/.$(SRR).ok
fastqc: $(FASTQC_OK)
index: $(GENOME_BWA_INDEX)
align: $(BAM)
stats: $(STATS)
bigwig: $(BW)

# ===== Recipes =====

# Download Genome
$(GENOME_FA):
	@mkdir -p $(GENOME_DIR)
	@echo "[-->] Fetching reference genome $(GENOME_ACC) and naming it $(GENOME_NAME)..."
	@curl -sSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$(GENOME_ACC)&rettype=fasta&retmode=text" > $@
	@test -s $@ || { echo "ERROR: Genome download failed for $(GENOME_ACC)."; exit 1; }

$(GENOME_FAI): $(GENOME_FA)
	@echo "[-->] Indexing genome with samtools faidx..."
	samtools faidx $(GENOME_FA)

# Download Reads
$(READS_DIR)/.$(SRR).ok:
	@mkdir -p $(READS_DIR)
	@echo "[-->] Downloading fastq for $(SRR)..."
	@if command -v fasterq-dump >/dev/null 2>&1; then \
		echo "[i] Using SRA Toolkit (fasterq-dump) for $(SRR)..."; \
		fasterq-dump --split-files -O $(READS_DIR) --threads 4 $(SRR); \
		if ls $(READS_DIR)/$(SRR)_1.fastq 1>/dev/null 2>&1; then \
			gzip -f $(READS_DIR)/$(SRR)_1.fastq $(READS_DIR)/$(SRR)_2.fastq; \
		elif ls $(READS_DIR)/$(SRR).fastq 1>/dev/null 2>&1; then \
			gzip -f $(READS_DIR)/$(SRR).fastq; \
		else \
			echo "ERROR: SRA Toolkit did not produce FASTQ files for $(SRR)."; exit 1; \
		fi; \
	else \
		echo "ERROR: fasterq-dump not found. Please install SRA Toolkit (conda install -c bioconda sra-tools)"; exit 1; \
	fi
	@touch $@

# Run FastQC
$(FASTQC_OK): reads
	@mkdir -p $(RESULTS_DIR)
	@echo "[-->] Running FastQC for $(SAMPLE) (SRR: $(SRR))..."
	fastqc -o $(RESULTS_DIR) $(READS_DIR)/$(SRR)*.fastq.gz
	@touch $@

# BWA Index
$(GENOME_BWA_INDEX): $(GENOME_FA)
	@echo "[-->] Indexing reference genome with BWA..."
	bwa index $(GENOME_FA)

# Align Reads
$(BAM): index genome reads
	@mkdir -p $(RESULTS_DIR)
	@echo "[-->] Aligning $(SRR) to create $(SAMPLE).sorted.bam..."
	@READ_ARGS=""; \
	if [ -f "$(R1)" ]; then \
		echo "  - Paired-end alignment"; \
		READ_ARGS="$(R1) $(R2)"; \
	else \
		echo "  - Single-end alignment"; \
		READ_ARGS="$(RS)"; \
	fi; \
	bwa mem -t 4 $(GENOME_FA) $$READ_ARGS | samtools view -b - | samtools sort -o "$(BAM)";
	@echo "[-->] Indexing BAM for $(SAMPLE)..."
	samtools index "$(BAM)"

# Generate Alignment Stats
$(STATS): $(BAM)
	@echo "[-->] Computing alignment stats for $(SAMPLE)..."
	@( \
	  echo "--- Alignment Report for Sample: $(SAMPLE) (SRR: $(SRR)) ---"; \
	  samtools flagstat "$(BAM)"; \
	  echo; echo "--- Per-reference read counts ---"; \
	  samtools idxstats "$(BAM)"; \
	) > "$(STATS)"
	@cat "$(STATS)"

# Generate BigWig Coverage File
$(BW): $(BAM)
	@echo "[-->] Generating BigWig coverage file for $(SAMPLE)..."
	@if ! command -v bamCoverage >/dev/null 2>&1; then \
		echo "[!] WARNING: bamCoverage not found (conda install -c bioconda deeptools). Skipping BigWig."; \
		touch "$(BW)"; exit 0; \
	fi
	bamCoverage --binSize 10 -b "$(BAM)" -o "$(BW)"

# Clean up generated files
clean:
	@echo "[-->] Cleaning all generated files..."
	rm -rf $(GENOME_DIR) $(READS_DIR) $(RESULTS_DIR) $(TMP_DIR)

