# ==============================================================================
# Simplified Makefile for BMMB852 Assignment 9
#
# - Fixes 'all' target by matching rules to the 'results/' directory.
# - Adds missing rules for .bw, .stats.txt, and _fastqc.html.
# ==============================================================================

# --- Variables ---
GENOME_DIR = genome
DATA_DIR = data
RESULTS_DIR = results # CHANGED: Was ALIGN_DIR, now matches 'all' target

# Assumes genome is manually placed here.
GENOME = $(GENOME_DIR)/reference.fasta

# --- Tool Definitions ---
BWA = bwa
SAMTOOLS = samtools
FASTQ_DUMP = fastq-dump
BAMCOVERAGE = bamCoverage # ADDED: For deeptools
FASTQC = fastqc         # ADDED: For FastQC

# --- Phony Targets ---
.PHONY: all clean
all: $(RESULTS_DIR)/$(SAMPLE).sorted.bam \
     $(RESULTS_DIR)/$(SAMPLE).sorted.bam.bai \
     $(RESULTS_DIR)/$(SAMPLE).bw \
     $(RESULTS_DIR)/$(SAMPLE).stats.txt \
     $(RESULTS_DIR)/$(SAMPLE)_fastqc.html

# --- Step 1: Index Genome (One-time step) ---
# This rule is run manually *once*, not per-sample.
$(GENOME).bwt: $(GENOME)
	@echo " Indexing reference genome..."
	$(BWA) index $<

# --- Step 2: Download Reads (Per-sample) ---
# This rule uses the SRR variable and the course-specific command.
$(DATA_DIR)/$(SRR).fastq:
	@echo " Downloading 10,000 reads for $(SRR)..."
	mkdir -p $(DATA_DIR)
	$(FASTQ_DUMP) -X 10000 -O $(DATA_DIR)/ $(SRR)

# --- Step 3: Align (Per-sample) ---
# This rule creates the SAM file. It depends on the genome index and the reads.
$(RESULTS_DIR)/$(SAMPLE).sam: $(GENOME).bwt $(DATA_DIR)/$(SRR).fastq
	@echo " Aligning $(SAMPLE)..."
	mkdir -p $(RESULTS_DIR)
	$(BWA) mem $(GENOME) $(word 2, $^) > $@ # Use 'word 2' for the .fastq dependency

# --- Step 4: Convert SAM to BAM (Per-sample) ---
$(RESULTS_DIR)/$(SAMPLE).bam: $(RESULTS_DIR)/$(SAMPLE).sam
	@echo " Converting $(SAMPLE) to BAM..."
	$(SAMTOOLS) view -bS $< > $@

# --- Step 5: Sort BAM (Per-sample) ---
$(RESULTS_DIR)/$(SAMPLE).sorted.bam: $(RESULTS_DIR)/$(SAMPLE).bam
	@echo "Sorting $(SAMPLE) BAM..."
	$(SAMTOOLS) sort $< -o $@

# --- Step 6: Index BAM (Per-sample) ---
$(RESULTS_DIR)/$(SAMPLE).sorted.bam.bai: $(RESULTS_DIR)/$(SAMPLE).sorted.bam
	@echo "Indexing $(SAMPLE) sorted BAM..."
	$(SAMTOOLS) index $<

# --- Step 7: Create Coverage Map (bamCoverage) --- [NEW RULE]
# Depends on the indexed sorted BAM.
$(RESULTS_DIR)/$(SAMPLE).bw: $(RESULTS_DIR)/$(SAMPLE).sorted.bam.bai
	@echo "Creating coverage map for $(SAMPLE)..."
	$(BAMCOVERAGE) -b $(RESULTS_DIR)/$(SAMPLE).sorted.bam -o $@

# --- Step 8: Get Alignment Stats (samtools stats) --- [NEW RULE]
# Depends on the sorted BAM.
$(RESULTS_DIR)/$(SAMPLE).stats.txt: $(RESULTS_DIR)/$(SAMPLE).sorted.bam
	@echo "Gathering stats for $(SAMPLE)..."
	$(SAMTOOLS) stats $< > $@

# --- Step 9: Run FastQC --- [NEW RULE]
# Depends on the raw reads.
$(RESULTS_DIR)/$(SAMPLE)_fastqc.html: $(DATA_DIR)/$(SRR).fastq
	@echo "Running FastQC on raw reads for $(SAMPLE)..."
	mkdir -p $(RESULTS_DIR)
	# FastQC creates outputs in the directory specified by -o
	$(FASTQC) $< -o $(RESULTS_DIR)

# --- Cleanup ---
clean:
	@echo "Cleaning up all generated files..."
	rm -rf $(RESULTS_DIR)
	rm -rf $(DATA_DIR)
	# Note: This does not remove the genome index.