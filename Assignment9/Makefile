# ==============================================================================
# Simplified Makefile for BMMB852 Assignment 9:post feedback from week 8
#
# - No internal loops.
# - Processes one sample at a time using variables (SRR, SAMPLE).
# - Uses fastq-dump -X 10000 as requested.
# ==============================================================================

# --- Variables ---
GENOME_DIR = genome
DATA_DIR = data
ALIGN_DIR = align

# Assumes genome is manually placed here.
GENOME = $(GENOME_DIR)/reference.fasta

# --- Tool Definitions ---
BWA = bwa
SAMTOOLS = samtools
FASTQ_DUMP = fastq-dump

# --- Phony Targets ---
.PHONY: all align clean

# --- Main Rule ---
# This target now depends on the final, indexed BAM for the *single sample*
# passed via the command line.
align: $(ALIGN_DIR)/$(SAMPLE).sorted.bam.bai
	@echo "✅ Pipeline complete for sample: $(SAMPLE)"

# --- Step 1: Index Genome (One-time step) ---
# This rule is run manually *once*, not per-sample.
$(GENOME).bwt: $(GENOME)
	@echo "🧬 Indexing reference genome..."
	$(BWA) index $<

# --- Step 2: Download Reads (Per-sample) ---
# This rule uses the SRR variable and the course-specific command.
$(DATA_DIR)/$(SRR).fastq:
	@echo "⏬ Downloading 10,000 reads for $(SRR)..."
	mkdir -p $(DATA_DIR)
	$(FASTQ_DUMP) -X 10000 -O $(DATA_DIR)/ $(SRR)

# --- Step 3: Align (Per-sample) ---
# This rule creates the SAM file. It depends on the genome index and the reads.
$(ALIGN_DIR)/$(SAMPLE).sam: $(GENOME).bwt $(DATA_DIR)/$(SRR).fastq
	@echo "🔬 Aligning $(SAMPLE)..."
	mkdir -p $(ALIGN_DIR)
	$(BWA) mem $(GENOME) $< > $@

# --- Step 4: Convert SAM to BAM (Per-sample) ---
$(ALIGN_DIR)/$(SAMPLE).bam: $(ALIGN_DIR)/$(SAMPLE).sam
	@echo "📦 Converting $(SAMPLE) to BAM..."
	$(SAMTOOLS) view -bS $< > $@

# --- Step 5: Sort BAM (Per-sample) ---
$(ALIGN_DIR)/$(SAMPLE).sorted.bam: $(ALIGN_DIR)/$(SAMPLE).bam
	@echo "Sorting $(SAMPLE) BAM..."
	$(SAMTOOLS) sort $< -o $@

# --- Step 6: Index BAM (Per-sample) ---
$(ALIGN_DIR)/$(SAMPLE).sorted.bam.bai: $(ALIGN_DIR)/$(SAMPLE).sorted.bam
	@echo "Indexing $(SAMPLE) sorted BAM..."
	$(SAMTOOLS) index $<

# --- Cleanup ---
clean:
	@echo "🧹 Cleaning up intermediate files..."
	rm -rf $(ALIGN_DIR)/*.sam
	rm -rf $(ALIGN_DIR)/*.bam
	rm -rf $(DATA_DIR)/*.fastq