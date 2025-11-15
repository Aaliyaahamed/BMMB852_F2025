# =====================================
# Somatic Variant Analysis (chr17:7M–9M)
# Custom Makefile for TP53 region
#
# VERSION 3: What did i fix?
#   - FIX: Downloads *only the region* for BAMs, not the full 90GB files.
#   - FIX: Corrects bcftools filter syntax ('&&' to '&').
#   - FIX: Creates tmp2 dir before touching files!!
# =====================================

SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --no-builtin-rules --warn-undefined-variables

# === Parameters ===
CHR := chr17
START := 7000000
END := 9000000
INTERVAL := $(CHR):$(START)-$(END)

# === Input URLs ===
REF_URL := https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
CTRL_URL := https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Element-AVITI-20241216/HG008-N-D_Element-StdInsert_77x_GRCh38-GIABv3.bam
TUMOR_URL := https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Element-AVITI-20241216/HG008-T_Element-StdInsert_111x_GRCh38-GIABv3.bam
GOLD_URL := https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/Ultima_DeepVariant-somatic-SNV-INDEL_20240822/HG008_edv_AF_recalibration.result.PASS.vcf.gz

# === Paths ===
DATA := data
REF := $(DATA)/ref.fasta.gz
CTRL := $(DATA)/ctrl.bam
TUMOR := $(DATA)/tumor.bam
VCF := results
OUT := report

# === VCF Outputs ===
CTRL_VCF := $(VCF)/ctrl.vcf.gz
TUMOR_VCF := $(VCF)/tumor.vcf.gz
UNIQUE := $(VCF)/tumor_unique.vcf.gz
FILTERED := $(VCF)/tumor_filtered.vcf.gz
GOLD := $(VCF)/gold.vcf.gz
GOLD_SUB := $(VCF)/gold_sub.vcf.gz

TP := $(VCF)/tp.vcf.gz
FP := $(VCF)/fp.vcf.gz
FN := $(VCF)/fn.vcf.gz

# === Targets ===
all: $(OUT)/summary.txt

# === Reference ===
# FIX: Added -C - to resume download if it gets interrupted.
$(REF):
	mkdir -p $(DATA)
	@echo "-> Downloading 3GB reference genome (will resume if interrupted)..."
	curl -L -C - $(REF_URL) -o $@
	@echo "-> Indexing reference genome..."
	samtools faidx $@

# === BAM Files (Region Extraction) ===
# FIX: This now downloads ONLY the region, not the full 70GB file.
$(CTRL):
	mkdir -p $(DATA)
	@echo "-> Downloading CONTROL BAM (chr17 region only)..."
	samtools view -b $(CTRL_URL) $(INTERVAL) > $@
	@echo "-> Indexing CONTROL BAM..."
	samtools index $@

$(TUMOR):
	mkdir -p $(DATA)
	@echo "-> Downloading TUMOR BAM (chr17 region only)..."
	samtools view -b $(TUMOR_URL) $(INTERVAL) > $@
	@echo "-> Indexing TUMOR BAM..."
	samtools index $@

# === Variant Calling ===
$(CTRL_VCF): $(CTRL) $(REF)
	@echo "-> Calling variants for CONTROL..."
	mkdir -p $(VCF)
	bcftools mpileup -f $(REF) -r $(INTERVAL) -d 100 $< | \
	bcftools call -mv -Oz -o $@
	tabix -p vcf $@

$(TUMOR_VCF): $(TUMOR) $(REF)
	@echo "-> Calling variants for TUMOR..."
	mkdir -p $(VCF)
	bcftools mpileup -f $(REF) -r $(INTERVAL) -d 100 $< | \
	bcftools call -mv -Oz -o $@
	tabix -p vcf $@

# === Tumor-Only Variants ===
$(UNIQUE): $(CTRL_VCF) $(TUMOR_VCF)
	@echo "-> Finding tumor-only variants..."
	mkdir -p $(VCF)/tmp1
	bcftools isec -p $(VCF)/tmp1 -n=2 $^
	bcftools view -Oz -o $@ $(VCF)/tmp1/0001.vcf
	tabix -p vcf $@ 2>/dev/null || true
	rm -rf $(VCF)/tmp1

# === Low-Threshold Filtering ===
# FIX: Changed '&&' to '&' for bcftools syntax.
$(FILTERED): $(UNIQUE)
	@echo "-> Filtering tumor-only variants..."
	bcftools filter -i 'QUAL > 10 & INFO/DP > 5' -Oz -o $@ $<
	tabix -p vcf $@ 2>/dev/null || true

# === Gold Standard Subset ===
$(GOLD):
	@echo "-> Downloading gold standard..."
	curl -L -C - $(GOLD_URL) -o $@
	curl -L -C - $(GOLD_URL).tbi -o $@.tbi

$(GOLD_SUB): $(GOLD)
	@echo "-> Subsetting gold standard..."
	bcftools view -r $(INTERVAL) -Oz -o $@ $<
	tabix -p vcf $@

# === Comparison ===
compare: $(FILTERED) $(GOLD_SUB)
	@echo "-> Comparing pipeline to gold standard..."
	@# FIX: Create the tmp2 directory *before* trying to use it.
	mkdir -p $(VCF)/tmp2
	@touch $(VCF)/tmp2/0000.vcf $(VCF)/tmp2/0001.vcf $(VCF)/tmp2/0002.vcf
	bcftools isec -p $(VCF)/tmp2 $^
	bcftools view -Oz -o $(TP) $(VCF)/tmp2/0002.vcf
	bcftools view -Oz -o $(FP) $(VCF)/tmp2/0000.vcf
	bcftools view -Oz -o $(FN) $(VCF)/tmp2/0001.vcf
	tabix -p vcf $(TP) 2>/dev/null || true
	tabix -p vcf $(FP) 2>/dev/null || true
	tabix -p vcf $(FN) 2>/dev/null || true
	rm -rf $(VCF)/tmp2

# === Report ===
$(OUT)/summary.txt: compare
	@echo "-> Generating final report..."
	mkdir -p $(OUT)
	TP_CT=$$(zgrep -vc "^#" $(TP) || true); \
	FP_CT=$$(zgrep -vc "^#" $(FP) || true); \
	FN_CT=$$(zgrep -vc "^#" $(FN) || true); \
	echo "==== TP53 Somatic Variant Report ====" > $@; \
	echo "Region: $(INTERVAL)" >> $@; \
	echo "Filter: 'QUAL > 10 & INFO/DP > 5'" >> $@; \
	echo "-----------------------------------" >> $@; \
	echo "True Positives:  $$TP_CT" >> $@; \
	echo "False Positives: $$FP_CT" >> $@; \
	echo "False Negatives: $$FN_CT" >> $@; \
	echo "-----------------------------------" >> $@; \
	if [ $$((TP_CT + FP_CT)) -gt 0 ]; then \
		P=$$(awk "BEGIN {printf \"%.2f\", 100 * $$TP_CT / ($$TP_CT + $$FP_CT)}"); \
		echo "Precision: $$P%" >> $@; \
	else echo "Precision: N/A" >> $@; fi; \
	if [ $$((TP_CT + FN_CT)) -gt 0 ]; then \
		R=$$(awk "BEGIN {printf \"%.2f\", 100 * $$TP_CT / ($$TP_CT + $$FN_CT)}"); \
		echo "Recall:    $$R%" >> $@; \
	else echo "Recall:    N/A" >> $@; fi
	@echo "==== REPORT COMPLETE ===="
	@cat $@

# === Cleanup ===
clean:
	rm -rf $(VCF) $(OUT)

clean-all:
	rm -rf $(VCF) $(OUT) $(DATA)