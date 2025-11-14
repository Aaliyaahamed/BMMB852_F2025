# This is a self-contained Makefile for Assignment 10.
# It downloads all necessary files, runs the full alignment
# and variant calling pipeline, and merges the VCFs.

# --- Shell Configuration ---
SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:

# --- Tool Configuration ---
BWA    := bwa
SAM    := samtools
SRA    := fastq-dump
BCF    := bcftools

# --- Directory and File Definitions ---
GENOME_DIR  := genome
DATA_DIR    := data
ALIGN_DIR   := align
RESULTS_DIR := results

# Reference Genome
REF         := $(GENOME_DIR)/reference.fasta
REF_FAI     := $(REF).fai
REF_BWA     := $(REF).bwt # BWA index file

# VCF for IGV
MERGED_VCF     := $(RESULTS_DIR)/all_samples.merged.vcf.gz
MERGED_VCF_IDX := $(MERGED_VCF).tbi

# --- PHONY TARGETS ---
.PHONY: help call_all merge_vcf clean run_samples

help:
	echo "Targets:"
	echo "  call      SRR=<SRR> SAMPLE=<SampleName> -> Creates a VCF for one sample"
	echo "  call_all                               -> Runs the entire pipeline for all samples."
	echo "  merge_vcf                            -> Merges all sample VCFs into one."
	echo "  clean                                -> Removes all generated files."

# --- Main Pipeline Target ---
# This runs everything: runs all samples, then merges the VCF.
call_all: $(MERGED_VCF)
	echo "\n[COMPLETE] I have finished the full pipeline!"
	echo "  Your final VCF is at: $(MERGED_VCF)"


# --- Main Variant Calling Pipeline ---
run_samples:
	echo "Step 1: Reading the design.csv file so I can process all your samples."
	# This awk command is fixed and will not corrupt
	awk -F, 'NR==1{for(i=1;i<=NF;i++){gsub(/\r/,"",$$i); if($$i=="SRR") srr_col=i; if($$i=="SampleName") sample_col=i}; next} {print "echo \"\nProcessing " $$sample_col " (SRR: " $$srr_col ") ...\"; make call SRR=" $$srr_col " SAMPLE=" $$sample_col}' design.csv | bash

# This rule merges the VCFs into our "master" file for IGV
$(MERGED_VCF): run_samples
	echo "\nStep 2: All individual samples are processed. I will now merge their VCF files."
	echo "[merge] I am collecting all the sample VCFs..."
	VCFS=$$(ls -1 $(RESULTS_DIR)/*.vcf.gz 2>/dev/null | grep -v "all_samples" || echo ""); \
	if [ -z "$$VCFS" ]; then echo "No VCFs found to merge."; exit 1; fi; \
	echo "I will now merge these files: $$VCFS"; \
	$(BCF) merge -Oz -o $(MERGED_VCF) $$VCFS
	echo "Merging complete. Now indexing the final VCF..."
	$(BCF) index -t $(MERGED_VCF)
	echo "Done! Your master VCF is: $(MERGED_VCF)"

# This is a hidden target to make merge_vcf work
merge_vcf: $(MERGED_VCF)

# --- Reference Genome Rules (FIXED) ---
# These rules build the indexes. They both depend on the reference file.
$(REF_BWA): $(REF)
	echo "Building BWA index..."
	$(BWA) index $(REF)

$(REF_FAI): $(REF)
	echo "Building FASTA index..."
	$(SAM) faidx $<

# This rule downloads the reference file.
# It explicitly depends on the GENOME_DIR target.
$(REF): $(GENOME_DIR)
	echo "I'm downloading the Ebola reference genome..."
	wget -O $(REF) "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.cgi?db=nuccore&id=KJ660346.2&rettype=fasta&retmode=text"

# This rule creates the directory. It is a prerequisite for $(REF).
$(GENOME_DIR):
	mkdir -p $(GENOME_DIR)

# --- Per-sample derived variables ---
FASTQ      := $(DATA_DIR)/$(SRR).fastq
BAM_FILE   := $(ALIGN_DIR)/$(SAMPLE).bam
SORTED_BAM := $(ALIGN_DIR)/$(SAMPLE).sorted.bam
BAM_INDEX  := $(ALIGN_DIR)/$(SAMPLE).sorted.bam.bai
VCF_FILE   := $(RESULTS_DIR)/$(SAMPLE).vcf.gz
VCF_INDEX  := $(VCF_FILE).tbi

# --- File-Generating Rules (Per-Sample) ---
# This is the target you run from the command line
call: $(VCF_FILE) $(VCF_INDEX)
	echo "Done: Finished call for $(SAMPLE) -> $(VCF_FILE)"

$(VCF_INDEX): $(VCF_FILE)
	$(BCF) index -t $<

$(VCF_FILE): $(BAM_INDEX) $(REF_FAI)
	mkdir -p $(RESULTS_DIR)
	$(BCF) mpileup -Ou -f $(REF) $(SORTED_BAM) \
	| $(BCF) call -mv -Ou --ploidy 1 \
	| $(BCF) norm -f $(REF) -Ou \
	| $(BCF) filter -s LOWCOV -e 'DP<10 || QUAL<30' \
	| $(BCF) view -Oz -o $@

$(BAM_INDEX): $(SORTED_BAM)
	$(SAM) index $<

$(SORTED_BAM): $(BAM_FILE)
	$(SAM) sort $< -o $@

# This rule now correctly depends on the reference indexes
$(BAM_FILE): $(FASTQ) $(REF_BWA)
	mkdir -p $(ALIGN_DIR)
	$(BWA) mem $(REF) $< | $(SAM) view -bS - > $@

# This rule depends on the DATA_DIR directory
$(FASTQ): $(DATA_DIR)
	$(SRA) -X 10000 -O $(DATA_DIR) $(SRR)

# This rule creates the data directory
$(DATA_DIR):
	mkdir -p $(DATA_DIR)

# --- Clean Rule ---
clean:
	echo "Cleaning up all generated directories and files..."
	rm -rf $(ALIGN_DIR) $(DATA_DIR) $(RESULTS_DIR)
	echo "Cleanup complete."