# Week 9 Assignment: Revising and Improving Your Automation Code
### Aaliya Ahamed • BMMB852 • 2025-10-26

This week's assignment automated the process of downloading, processing, and aligning high-throughput sequencing data. It used a Makefile to define the workflow for a single sample and leveraged GNU parallel to efficiently scale the analysis across multiple samples defined in a design.csv file.

The pipeline was based on the data from the paper on the 2014 Ebola virus, BioProject PRJNA257197.

# Analysis Workflow

The core pipeline was a series of automated tasks performed on each sample:

Reference Genome Acquisition: The process began by fetching the specified reference genome sequence from NCBI.

Read Data Retrieval: Raw sequencing reads for each sample run were downloaded from the Sequence Read Archive (SRA) using the fasterq-dump utility.

Quality Assurance: The integrity and quality of the raw reads were assessed using FastQC to generate a QC report.

Sequence Alignment: The downloaded reads were mapped to the reference genome using the BWA-MEM algorithm.

BAM File Generation: The alignment output was converted into a sorted, indexed BAM file, the standard format for downstream analysis.

Coverage Map Creation: A BigWig (BW) file was generated from the alignment to provide a coverage track, which was useful for visualization in genome browsers like IGV.

Summary Statistics: A final report was generated containing key alignment metrics, such as mapping rates, from samtools.


# Design.csv

The pipeline's batch functionality was driven by a sample manifest file named design.csv. Proper formatting of this file was critical for successful execution. The file was saved in a standard comma-separated values (CSV) format.

Example design.csv Structure:
```bash
SRR,SampleName
SRR1972976,G1749.1
SRR1972977,G1750.1
SRR1972978,G1751.1
```

# Execution

The Makefile was executed from the terminal, where all necessary parameters for a given run were provided as arguments.

1. Initial Sanity Check: A Single Run

```bash
make all SRR=SRR1972976 SAMPLE=G1749.1 GENOME_ACC=AF086833 GENOME_NAME=ebola_zaire
```

_A successful execution resulted in the creation of genome/, reads/, and results/ directories containing the output for the G1749.1 sample_.

2. Full-Scale Analysis: Batch Processing

To process multiple samples, GNU parallel was used to iterate through the design.csv manifest and launch a make job for each entry.

Command to process the first 10 samples:
```bash

cat design.csv | parallel --colsep , --header : \
    make align SRR={SRR} SAMPLE={SampleName}
```

# Directory Structure and Key Artifacts

After a successful run, the project folder was organized as follows, containing the key generated files:
```bash

.
├── design.csv
├── Makefile
├── README.md
├── genome/
│   └── ebola_zaire.fna         # Reference genome
│   └── ... (BWA index files)
├── reads/
│   └── SRR1972976_1.fastq.gz   # Raw sequencing reads
│   └── SRR1972976_2.fastq.gz
│   └── ... (other read files)
└── results/
    ├── G1749.1.sorted.bam      # Main alignment output
    ├── G1749.1.sorted.bam.bai  # Index for the BAM file
    ├── G1749.1.bw              # Genome coverage map
    ├── G1749.1.stats.txt       # Text-based alignment metrics
    └── G1749.1_fastqc.html     # HTML quality control report
    └── ... (outputs for all other samples)
```