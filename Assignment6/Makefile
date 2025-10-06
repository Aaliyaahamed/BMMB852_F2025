# Makefile for Ebola Makona 2014 Alignment (Assignment 6)
# Author: Aaliya Ahamed

SHELL := /bin/bash
# Execute all commands in a single shell instance for each rule, allowing variables to be shared.
.ONESHELL:
# Run the shell with strict error checking: exit on error, fail on pipe errors.
.SHELLFLAGS := -eu -o pipefail -c
# Automatically delete target files if a command in the rule fails.
.DELETE_ON_ERROR:
# Warn if a variable is used without being defined.
MAKEFLAGS += --warn-undefined-variables
# Disable built-in rules (like converting .c to .o files) to avoid unexpected behavior.
MAKEFLAGS += --no-builtin-rules

# --- Variables--
# Dataset from Gire et al., 2014, Science (unique to this project)
ACC          = KM034562.1
SRR          = SRR1553606
SPOTS        = 5000
THREADS      = 4

# Organized directories
GENOME_DIR   = genome
READS_DIR    = reads
RESULTS_DIR  = results

# File paths (the final targets of our workflow)
GENOME_FA    = $(GENOME_DIR)/ebola_genome.fna
READ1_GZ     = $(READS_DIR)/$(SRR)_1.fastq.gz
READ2_GZ     = $(READS_DIR)/$(SRR)_2.fastq.gz
BAM          = $(RESULTS_DIR)/$(SRR).sorted.bam
STATS_FILE   = $(RESULTS_DIR)/alignment_stats.txt

# --- Phony Targets (Commands that are not files) ---
.PHONY: all qc index align stats clean

# --- Main Pipeline ---
# The default target runs the entire workflow and produces the final stats file.
all: $(STATS_FILE)

# --- Workflow Steps ---

# 1. Download Reference Genome
# This rule creates the genome file if it doesn't exist.
$(GENOME_FA):
	@echo "[-->] Creating genome directory..."
	mkdir -p $(GENOME_DIR)
	@echo "[-->] Fetching reference genome ${ACC}..."
	efetch -db nuccore -id ${ACC} -format fasta > $@

# 2. Download, Subset, and Compress Sequencing Reads
# This is a grouped target. If either .gz file is missing, this recipe runs.
# It creates both target files in one step.
$(READ1_GZ) $(READ2_GZ):
	@echo "[-->] Creating reads directory..."
	mkdir -p $(READS_DIR)
	@echo "[-->] Downloading and subsetting ${SPOTS} spots for ${SRR}..."
	fastq-dump --split-files -X ${SPOTS} -O $(READS_DIR) ${SRR}
	@echo "[-->] Compressing FASTQ files with gzip..."
	# Gzip the uncompressed files to create the final target .gz files
	gzip -f $(READS_DIR)/$(SRR)_1.fastq
	gzip -f $(READS_DIR)/$(SRR)_2.fastq

# 3. Quality Control (Optional target, not run by default)
qc: $(READ1_GZ) $(READ2_GZ)
	@echo "[-->] Running FastQC on reads..."
	fastqc -o $(READS_DIR) $(READ1_GZ) $(READ2_GZ)

# 4. Index Reference Genome
# This is a prerequisite for the alignment step.
index: $(GENOME_FA).bwt
$(GENOME_FA).bwt: $(GENOME_FA)
	@echo "[-->] Indexing reference genome with BWA..."
	bwa index $(GENOME_FA)

# 5. Align Reads to Genome and Create a Sorted BAM
# The final alignment file depends on the indexed genome and the compressed reads.
align: $(BAM)
$(BAM): $(GENOME_FA).bwt $(READ1_GZ) $(READ2_GZ)
	@echo "[-->] Creating results directory..."
	mkdir -p $(RESULTS_DIR)
	@echo "[-->] Aligning reads with BWA-MEM and sorting..."
	bwa mem -t $(THREADS) $(GENOME_FA) $(READ1_GZ) $(READ2_GZ) | samtools sort -o $(BAM)
	@echo "[-->] Indexing sorted BAM file..."
	samtools index $(BAM)

# 6. Generate Final Statistics Report
# This creates the final text report and prints it to the screen.
stats: $(STATS_FILE)
$(STATS_FILE): $(BAM)
	@echo "[-->] Generating final statistics report..."
	# Use '>' to create the file and write the header
	@echo "--- Alignment and Coverage Report for SRR: $(SRR) ---" > $@
	# Use '>>' to append the rest of the stats to the file
	@echo "--- Mapping Statistics (samtools flagstat) ---" >> $@
	samtools flagstat $(BAM) >> $@
	@echo -e "\n--- Coverage Statistics (samtools depth) ---" >> $@
	samtools depth $(BAM) \
	  | awk '{s+=$$3; n++} END{if(n>0) printf "Average coverage = %.2fx\n", s/n; else print "Average coverage = 0x"}' \
	  >> $@
	@echo "[✓] Analysis complete. Report saved to: $(STATS_FILE)"
	@echo "----------------------------------------------------"
	@cat $(STATS_FILE)

# --- Utility Targets ---
# A rule to clean up all generated files and start fresh.
clean:
	@echo "[-->] Removing all generated directories and files..."
	rm -rf $(GENOME_DIR) $(READS_DIR) $(RESULTS_DIR)

