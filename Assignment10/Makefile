# ==========================
# Assignment 10 Makefile
# ==========================

SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:

# ---- Configuration
GENOME_DIR := genome
DATA_DIR   := data
ALIGN_DIR  := align
RESULTS    := results
REF 	   := $(GENOME_DIR)/reference.fasta
REFFAI 	   := $(REF).fai

# ---- Tools
BWA := bwa
SAM := samtools
SRA := fastq-dump
BCF := bcftools

# ================ 
# PHONY TARGETS
# ================ 
.PHONY: help call_all call merge_vcf clean

help:
	@echo "Targets:"
	@echo "  call        SRR=<SRR> SAMPLE=<SampleName>    -> Creates a VCF for one sample"
	@echo "  call_all                                   -> Runs 'call' for all samples in design.csv, then merges."
	@echo "  merge_vcf                                  -> Merges all VCFs from design.csv"
	@echo "  clean"

# =========================
# FULL PIPELINE
# =========================
# This target reads the design.csv, then runs 'make call' for each sample,
# and finishes by running 'make merge_vcf'.
call_all:
	@echo "[STEP 1] Reading design.csv to call all samples..."
	@awk -F, 'NR==1{for(i=1;i<=NF;i++){gsub(/\r/,"",$$i); if($$i=="SRR") srr_col=i; if($$i=="SampleName") sample_col=i}; next} {print "echo \"\n--- Processing " $$sample_col " (SRR: " $$srr_col ") ---\"; make call SRR=" $$srr_col " SAMPLE=" $$sample_col}' design.csv | bash
	@echo "\n[STEP 2] All samples processed. Now merging VCFs..."
	@$(MAKE) merge_vcf
	@echo "\n[COMPLETE] Full pipeline finished."

# ----------------
# FILE DEFINITIONS 
# ----------------
FASTQ 	   := $(DATA_DIR)/$(SRR).fastq
BAM_FILE   := $(ALIGN_DIR)/$(SAMPLE).bam
SORTED_BAM := $(ALIGN_DIR)/$(SAMPLE).sorted.bam
BAM_INDEX  := $(ALIGN_DIR)/$(SAMPLE).sorted.bam.bai
VCF_FILE   := $(RESULTS)/$(SAMPLE).vcf.gz
VCF_INDEX  := $(RESULTS)/$(SAMPLE).vcf.gz.tbi



# Main target to create a single VCF
call: $(VCF_FILE) $(VCF_INDEX)
	@echo "Done: call for $(SAMPLE) -> $(VCF_FILE)"

# Rule to build VCF
$(VCF_FILE): $(BAM_INDEX) $(REFFAI)
	mkdir -p $(RESULTS)
	$(BCF) mpileup -Ou -f $(REF) $(SORTED_BAM) \
	| $(BCF) call -mv -Ou --ploidy 1 \
	| $(BCF) norm -f $(REF) -Ou \
	| $(BCF) filter -s LOWCOV -e 'DP<10 || QUAL<30' \
	| $(BCF) view -Oz -o $@
	$(BCF) index -t $@

# Rule to build BAM Index
$(BAM_INDEX): $(SORTED_BAM)
	$(SAM) index $<

# Rule to build Sorted BAM
$(SORTED_BAM): $(BAM_FILE)
	$(SAM) sort $< -o $@

# Rule to build BAM 
$(BAM_FILE): $(FASTQ) $(REFFAI)
	mkdir -p $(ALIGN_DIR)
	$(BWA) mem $(REF) $< | $(SAM) view -bS - > $@

# Rule to build FASTQ
$(FASTQ):
	mkdir -p $(DATA_DIR)
	$(SRA) -X 10000 -O $(DATA_DIR) $(SRR)

# Rule to build FASTA Index
$(REFFAI): $(REF)
	$(SAM) faidx $<

# ----------------
# MERGE & CLEAN
# ----------------
merge_vcf:
	@echo "[merge] collecting sample VCFs from design.csv ..."
	VCFS=$$(awk -F, 'NR==1{for(i=1;i<=NF;i++){gsub(/\r/,"",$$i); if($$i=="SampleName") c=i}; next} {print "results/"$$c".vcf.gz"}' design.csv | sort -u); \
	echo "Merging: $$VCFS"; \
	$(BCF) merge -Oz -o $(RESULTS)/all_samples.merged.vcf.gz $$VCFS
	$(BCF) index -t $(RESULTS)/all_samples.merged.vcf.gz
	@echo "Merged VCF: results/all_samples.merged.vcf.gz"

clean:
	rm -rf $(ALIGN_DIR)/*.sam $(ALIGN_DIR)/*.bam $(ALIGN_DIR)/*.bai $(RESULTS)/*.vcf.gz*

