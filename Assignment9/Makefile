# Variables
GENOME_DIR = genome
DATA_DIR = data
RESULTS_DIR = results

# Path to the reference genome fasta
GENOME = $(GENOME_DIR)/reference.fasta

BWA = bwa
SAMTOOLS = samtools
FASTQ_DUMP = fastq-dump
BAMCOVERAGE = bamCoverage 
FASTQC = fastqc
EFETCH = efetch

# Phony Targets
.PHONY: all clean

# after feedback 2
all: $(RESULTS_DIR)/$(SAMPLE).sorted.bam \
     $(RESULTS_DIR)/$(SAMPLE).sorted.bam.bai \
     $(RESULTS_DIR)/$(SAMPLE).bw \
     $(RESULTS_DIR)/$(SAMPLE).stats.txt \
     $(RESULTS_DIR)/$(SAMPLE)_fastqc.html

#I have fixed this from before
$(GENOME):
    @echo " Fetching genome $(GENOME_ACC)..."
    mkdir -p $(GENOME_DIR)
    $(EFETCH) -db nuccore -id $(GENOME_ACC) -format fasta > $@

# Index the genome
$(GENOME).bwt: $(GENOME)
    @echo " Indexing reference genome..."
    $(BWA) index $<

# Download the reads
$(DATA_DIR)/$(SRR).fastq:
    @echo " Downloading 10,000 reads for $(SRR)..."
    mkdir -p $(DATA_DIR)
    $(FASTQ_DUMP) -X 10000 -O $(DATA_DIR)/ $(SRR)

# Align reads to genome
$(RESULTS_DIR)/$(SAMPLE).sam: $(GENOME).bwt $(DATA_DIR)/$(SRR).fastq
    @echo " Aligning $(SAMPLE)..."
    mkdir -p $(RESULTS_DIR)
    # $(GENOME) is the ref, $(word 2, $^) is the .fastq file
    $(BWA) mem $(GENOME) $(word 2, $^) > $@ 

# Convert SAM to BAM
$(RESULTS_DIR)/$(SAMPLE).bam: $(RESULTS_DIR)/$(SAMPLE).sam
    @echo " Converting $(SAMPLE) to BAM..."
    $(SAMTOOLS) view -bS $< > $@

# Sort the BAM
$(RESULTS_DIR)/$(SAMPLE).sorted.bam: $(RESULTS_DIR)/$(SAMPLE).bam
    @echo "Sorting $(SAMPLE) BAM..."
    $(SAMTOOLS) sort $< -o $@

# Index the sorted BAM
$(RESULTS_DIR)/$(SAMPLE).sorted.bam.bai: $(RESULTS_DIR)/$(SAMPLE).sorted.bam
    @echo "Indexing $(SAMPLE) sorted BAM..."
    $(SAMTOOLS) index $<

# Create coverage map.
$(RESULTS_DIR)/$(SAMPLE).bw: $(RESULTS_DIR)/$(SAMPLE).sorted.bam.bai
    @echo "Creating coverage map for $(SAMPLE)..."
    $(BAMCOVERAGE) -b $(RESULTS_DIR)/$(SAMPLE).sorted.bam -o $@

# Get stats
#Deepends on the BAM
$(RESULTS_DIR)/$(SAMPLE).stats.txt: $(RESULTS_DIR)/$(SAMPLE).sorted.bam
    @echo "Gathering stats for $(SAMPLE)..."
    $(SAMTOOLS) stats $< > $@

# Run FastQC
$(RESULTS_DIR)/$(SAMPLE)_fastqc.html: $(DATA_DIR)/$(SRR).fastq
    @echo "Running FastQC on raw reads for $(SAMPLE)..."
    mkdir -p $(RESULTS_DIR)
    # FastQC creates outputs in the directory specified by -o
    $(FASTQC) $< -o $(RESULTS_DIR)

#Cleanup
clean:
    @echo "Cleaning up all generated files..."
    rm -rf $(RESULTS_DIR)
    rm -rf $(DATA_DIR)
    rm -rf $(GENOME_DIR)