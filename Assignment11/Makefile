# This is the complete Makefile for Assignment 11.
# It runs the full variant pipeline and generates a
# variant effect report using snpEff.

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

# --- snpEff Configuration ---
SNPEFF_JAR     := snpEff/snpEff.jar
SNPEFF_CONFIG  := snpEff/snpEff.config
SNPEFF_DB_NAME := ebola_zaire
# This is the name snpEff's database uses for the chromosome
SNPEFF_CHR_NAME := KJ660346.1

# --- Directory and File Definitions ---
GENOME_DIR  := genome
DATA_DIR    := data
ALIGN_DIR   := align
RESULTS_DIR := results

# Reference Genome
REF         := $(GENOME_DIR)/reference.fasta
REF_FAI     := $(REF).fai
REF_BWA     := $(REF).bwt
# This is the name our FASTA file uses for the chromosome
FASTA_CHR_NAME := KJ660346.2

# Final VCF
MERGED_VCF     := $(RESULTS_DIR)/all_samples.merged.vcf.gz

# Files for snpEff
CHROMS_MAP   := chrom_map.txt
RENAMED_VCF  := $(RESULTS_DIR)/all_samples.for_snpeff.vcf.gz
SNPEFF_HTML  := $(RESULTS_DIR)/snpEff_report.html
# This is the annotated VCF file that `cat` can read
SNPEFF_VCF   := $(RESULTS_DIR)/all_samples.ann.vcf

# --- PHONY TARGETS ---
.PHONY: help call_all merge_vcf clean annotate run_samples

help:
	echo "Targets:"
	echo "  call      SRR=<SRR> SAMPLE=<SampleName> -> Creates a VCF for one sample"
	echo "  call_all                               -> Runs the entire pipeline for all samples."
	echo "  merge_vcf                            -> Merges all sample VCFs into one."
	echo "  annotate                             -> Creates the snpEff report."
	echo "  clean                                -> Removes all generated files."

# --- Main Pipeline Target ---
# This runs everything: runs all samples, merges the VCF,
# and finally runs the snpEff annotation.
call_all: $(SNPEFF_HTML)
	echo "\n[COMPLETE] I have finished the full pipeline!"
	echo "  Your snpEff report is at: $(SNPEFF_HTML)"
	echo "  Your annotated VCF is at: $(SNPEFF_VCF)"

# --- snpEff Annotation Pipeline ---
# This target creates the HTML report.
annotate: $(SNPEFF_HTML)

$(SNPEFF_HTML): $(RENAMED_VCF) $(SNPEFF_JAR)
	echo "I'm annotating the variants with snpEff..."
	# This command saves the HTML stats AND saves the text VCF file
	java -Xmx4g -jar $(SNPEFF_JAR) -v $(SNPEFF_DB_NAME) -c $(SNPEFF_CONFIG) \
		-stats $(SNPEFF_HTML) \
		$(RENAMED_VCF) > $(SNPEFF_VCF)
	echo "snpEff report is at: $(SNPEFF_HTML)"

# This rule creates the VCF for snpEff by renaming the chromosome
$(RENAMED_VCF): $(MERGED_VCF) $(CHROMS_MAP)
	echo "I'm creating a special VCF for snpEff (renaming chromosome)..."
	$(BCF) annotate --rename-chrs $(CHROMS_MAP) $(MERGED_VCF) -Oz -o $(RENAMED_VCF)

# This rule creates the map file snpEff needs
$(CHROMS_MAP):
	echo "Creating chromosome map for snpEff..."
	echo "$(FASTA_CHR_NAME) $(SNPEFF_CHR_NAME)" > $(CHROMS_MAP)

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

# --- Per-sample derived variables ---
FASTQ      := $(DATA_DIR)/$(SRR).fastq
BAM_FILE   := $(ALIGN_DIR)/$(SAMPLE).bam
SORTED_BAM := $(ALIGN_DIR)/$(SAMPLE).sorted.bam
BAM_INDEX  := $(ALIGN_DIR)/$(SAMPLE).sorted.bam.bai
VCF_FILE   := $(RESULTS_DIR)/$(SAMPLE).vcf.gz
VCF_INDEX  := $(VCF_FILE).tbi

# --- File-Generating Rules (Per-Sample) ---
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

$(BAM_FILE): $(FASTQ) $(REF_BWA)
	mkdir -p $(ALIGN_DIR)
	$(BWA) mem $(REF) $< | $(SAM) view -bS - > $@

$(FASTQ):
	mkdir -p $(DATA_DIR)
	$(SRA) -X 10000 -O $(DATA_DIR) $(SRR)

# --- Reference Genome and Annotation Rules ---
# This rule downloads the reference.fasta if it is missing.
$(REF):
	echo "I'm downloading the Ebola reference genome..."
	mkdir -p $(GENOME_DIR)
	wget -O $(REF) "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.cgi?db=nuccore&id=KJ660346.2&rettype=fasta&retmode=text"

$(REF_BWA): $(REF)
	echo "Building BWA index..."
	$(BWA) index $(REF)

$(REF_FAI): $(REF)
	echo "Building FASTA index..."
	$(SAM) faidx $<

# --- Clean Rule ---
clean:
	echo "Cleaning up all generated directories and files..."
	rm -rf $(ALIGN_DIR) $(DATA_DIR) $(RESULTS_DIR) $(CHROMS_MAP)
	echo "Cleanup complete."