# Makefile
# This is Assignment 11: Ebola Variant Calling Pipeline

SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables

# -Variables 

SRR     ?= SRR0000000
SAMPLE  ?= SampleX

GENOME  := genome
DATA    := data
ALIGN   := align
RESULTS := results

REF     := $(GENOME)/reference.fasta
REF_BWA := $(REF).bwt
REF_FAI := $(REF).fai

FASTQ       := $(DATA)/$(SRR).fastq
BAM         := $(ALIGN)/$(SAMPLE).bam
SORTED_BAM  := $(ALIGN)/$(SAMPLE).sorted.bam
BAM_IDX     := $(SORTED_BAM).bai

VCF         := $(RESULTS)/$(SAMPLE).vcf.gz
VCF_IDX     := $(VCF).tbi
MERGED_VCF  := $(RESULTS)/all_samples.merged.vcf.gz

CHROMS_MAP  := chrom_map.txt
RENAMED_VCF := $(RESULTS)/all_samples.for_snpeff.vcf.gz
SNPEFF_VCF  := $(RESULTS)/all_samples.ann.vcf
SNPEFF_HTML := $(RESULTS)/snpEff_report.html

# Tool paths
BWA    := bwa
SAM    := samtools
SRA    := fastq-dump
BCF    := bcftools

SNPEFF_JAR    := snpEff/snpEff.jar
SNPEFF_CONFIG := snpEff/snpEff.config
SNPEFF_DB     := ebola_zaire

FASTA_CHR_NAME := KJ660346.2
SNPEFF_CHR_NAME := KJ660346.1

#Pipeline 

all: call_all

call_all: run_samples merge_vcf annotate
	@echo "\n Pipeline complete. Outputs:"
	@echo " - $(MERGED_VCF)"
	@echo " - $(SNPEFF_VCF)"
	@echo " - $(SNPEFF_HTML)"

run_samples:
	@echo " Reading design.csv..."
	awk -F, 'NR==1{for(i=1;i<=NF;i++){gsub(/\r/,"",$$i); if($$i=="SRR") srr_col=i; if($$i=="SampleName") sample_col=i}; next} {print "make call SRR=" $$srr_col " SAMPLE=" $$sample_col}' design.csv | bash

call: $(VCF) $(VCF_IDX)
	@echo " Variant call complete for $(SAMPLE)"

$(FASTQ):
	mkdir -p $(DATA)
	$(SRA) -X 10000 -O $(DATA) $(SRR)

$(BAM): $(FASTQ) $(REF_BWA)
	mkdir -p $(ALIGN)
	$(BWA) mem $(REF) $< | $(SAM) view -bS - > $@

$(SORTED_BAM): $(BAM)
	$(SAM) sort $< -o $@

$(BAM_IDX): $(SORTED_BAM)
	$(SAM) index $<

$(VCF): $(BAM_IDX) $(REF_FAI)
	mkdir -p $(RESULTS)
	$(BCF) mpileup -Ou -f $(REF) $(SORTED_BAM) \
	| $(BCF) call -mv -Ou --ploidy 1 \
	| $(BCF) norm -f $(REF) -Ou \
	| $(BCF) filter -s LOWCOV -e 'DP<10 || QUAL<30' \
	| $(BCF) view -Oz -o $@

$(VCF_IDX): $(VCF)
	$(BCF) index -t $<

# Merge VCF 

merge_vcf: $(MERGED_VCF)

$(MERGED_VCF): run_samples
	VCFS=$$(ls -1 $(RESULTS)/*.vcf.gz 2>/dev/null | grep -v all_samples || echo ""); \
	if [ -z "$$VCFS" ]; then echo " No VCFs found to merge."; exit 1; fi; \
	echo " Merging VCFs..."; \
	$(BCF) merge -Oz -o $@ $$VCFS; \
	$(BCF) index -t $@

# Annotation

annotate: $(SNPEFF_HTML)

$(CHROMS_MAP):
	echo "$(FASTA_CHR_NAME) $(SNPEFF_CHR_NAME)" > $(CHROMS_MAP)

$(RENAMED_VCF): $(MERGED_VCF) $(CHROMS_MAP)
	$(BCF) annotate --rename-chrs $(CHROMS_MAP) $(MERGED_VCF) -Oz -o $@

$(SNPEFF_HTML): $(RENAMED_VCF) $(SNPEFF_JAR)
	java -Xmx4g -jar $(SNPEFF_JAR) -v $(SNPEFF_DB) -c $(SNPEFF_CONFIG) \
		-stats $@ \
		$(RENAMED_VCF) > $(SNPEFF_VCF)

# Reference

$(REF):
	mkdir -p $(GENOME)
	wget -O $@ "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.cgi?db=nuccore&id=KJ660346.2&rettype=fasta&retmode=text"

$(REF_BWA): $(REF)
	$(BWA) index $<

$(REF_FAI): $(REF)
	$(SAM) faidx $<

#Utilities 

clean:
	rm -rf $(ALIGN) $(DATA) $(RESULTS) $(CHROMS_MAP)

help:
	@echo "🛠  Makefile Commands:"
	@echo "  make call SRR=SRRXXXX SAMPLE=Name   # Call variants for a sample"
	@echo "  make call_all                       # Full pipeline"
	@echo "  make clean                          # Remove all outputs"
