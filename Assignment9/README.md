# Week 9 Assignment: Automate a large scale analysis
### Aaliya Ahamed • BMMB852 • 2025-10-26
# redo it - merge this with the same branch

This week's assignment automated the process of downloading, processing, and aligning high-throughput sequencing data. It used a Makefile to define the workflow for a single sample and leveraged GNU parallel to efficiently scale the analysis across multiple samples defined in a design.csv file.

The pipeline was based on the data from the paper on the 2014 Ebola virus, BioProject PRJNA257197.

1. Initial Sanity Check: A Single Run

It was strongly recommended to perform a test run on a single sample. This verified that all software dependencies were met and the configuration was correct before launching a large batch job. The command was executed from the project's root directory.
```bash
make all SRR=SRR1972976 SAMPLE=G1749.1 GENOME_ACC=AF086833 GENOME_NAME=ebola_zaire
```

_A successful execution resulted in the creation of genome/, reads/, and results/ directories containing the output for the G1749.1 sample_.

2. Full-Scale Analysis: Batch Processing

To process multiple samples, GNU parallel was used to iterate through the design.csv manifest and launch a make job for each entry.

Command to process the first 10 samples:
```bash

head -n 11 design.csv | parallel --colsep ',' --header : --lb -j  "make all SRR={SRR} SAMPLE={SampleName} GENOME_ACC=AF086833 GENOME_NAME=ebola_zaire"
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