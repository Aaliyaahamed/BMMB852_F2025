# Define the reference genome and annotation
REF = refs/chr22.genome.fa
GTF = refs/chr22.gtf
INDEX_BASE = refs/hisat2_index

# List of samples (matches the roots in your reads/ folder)
SAMPLES = HBR_1 HBR_2 HBR_3 UHR_1 UHR_2 UHR_3

# Define the final output files we want
BAMS = $(patsubst %, bam/%.bam, $(SAMPLES))
BIGWIGS = $(patsubst %, bam/%.bw, $(SAMPLES))
COUNTS = counts.txt

# Target: all
# This runs the whole pipeline when you type 'make'
all: index $(BAMS) $(BIGWIGS) $(COUNTS)

# Target: index
# Builds the Hisat2 index if it doesn't exist
index: $(INDEX_BASE).1.ht2
$(INDEX_BASE).1.ht2:
	hisat2-build $(REF) $(INDEX_BASE)

# Target: Alignment
# Aligns reads, sorts them, and indexes the BAM
bam/%.bam: reads/%_R1.fq
	mkdir -p bam
	# Align with Hisat2 and pipe to samtools sort
	hisat2 -x $(INDEX_BASE) -U $< | samtools sort > $@
	# Index the bam file
	samtools index $@

# Target: BigWig
# Creates coverage tracks for IGV
bam/%.bw: bam/%.bam
	# Create a genome index if missing (needed for bigwig conversion)
	samtools faidx $(REF)
	# Convert BAM to BedGraph
	bedtools genomecov -ibam $< -split -bg > bam/$*.bg
	# Convert BedGraph to BigWig
	bedGraphToBigWig bam/$*.bg $(REF).fai $@
	# Clean up the intermediate bedgraph
	rm bam/$*.bg

# Target: Counting
# Runs featureCounts to create the matrix
$(COUNTS): $(BAMS)
	featureCounts -a $(GTF) -o $@ $(BAMS)

# Target: clean
# Removes generated files
clean:
	rm -rf bam/*.bam bam/*.bai bam/*.bw counts.txt counts.txt.summary