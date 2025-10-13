
# Reusable Alignment Makefile

SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# ===== Require params only for pipeline targets (not for clean/help) =====
NEED_PARAMS_FOR := all genome reads links index align bigwig stats primary maxdepth gene_forward
ifneq ($(filter $(NEED_PARAMS_FOR),$(MAKECMDGOALS)),)
  ifndef SRR
    $(error Provide SRR=..., e.g., make all SRR=SRR1972739 ACC=U00096.3)
  endif
  ifndef ACC
    $(error Provide ACC=..., e.g., make all SRR=SRR1972739 ACC=U00096.3)
  endif
endif

# ===== User-tunable knobs =====
# BW_METHOD can be: auto | deeptools | bedtools
BW_METHOD ?= auto
BIN_SIZE  ?= 10  # for bamCoverage

# Optional region for the "gene_forward" query:
#   make gene_forward SRR=... ACC=... CHR=NC_000913.3 START=190000 END=195000
CHR    ?=
START  ?=
END    ?=

# ===== Directories =====
GENOME_DIR  := genome
READS_DIR   := reads
RESULTS_DIR := results
TMP_DIR     := .tmp

# ===== Files =====
GENOME_FA  := $(GENOME_DIR)/$(ACC).fna
GENOME_FAI := $(GENOME_FA).fai
GENOME_BWA_INDEX := $(GENOME_FA).bwt

READS_FIELDS := fastq_https,fastq_ftp
READS_MANIFEST := $(TMP_DIR)/$(SRR).ena.tsv

# Normalized links (created regardless of original names)
R1 := $(READS_DIR)/$(SRR)_1.fastq.gz
R2 := $(READS_DIR)/$(SRR)_2.fastq.gz
RS := $(READS_DIR)/$(SRR).fastq.gz

BAM := $(RESULTS_DIR)/$(SRR).sorted.bam
BAI := $(BAM).bai

# Coverage/summary outputs
BEDGRAPH := $(RESULTS_DIR)/$(SRR).bedgraph
BW       := $(RESULTS_DIR)/$(SRR).bw
STATS    := $(RESULTS_DIR)/$(SRR).stats.txt
PRIMARY  := $(RESULTS_DIR)/$(SRR).primary.txt
MAXDEPTH := $(RESULTS_DIR)/$(SRR).maxdepth.txt
GENEFWD  := $(RESULTS_DIR)/$(SRR).gene_forward.txt

# ===== Phony =====
.PHONY: all tools genome reads links index align bigwig stats primary maxdepth gene_forward clean help

# ===== Top-level =====
all: tools $(BAM) $(BW) $(STATS)
	@echo
	@echo "[✓] All steps completed successfully for $(SRR). Outputs in $(RESULTS_DIR)/"

help:
	@echo "Examples:"
	@echo "  make all SRR=SRR519924 ACC=U00096.3"
	@echo "  make bigwig SRR=... ACC=... [BW_METHOD=deeptools|bedtools] [BIN_SIZE=10]"
	@echo "  make primary  SRR=... ACC=..."
	@echo "  make maxdepth SRR=... ACC=..."
	@echo "  make gene_forward SRR=... ACC=... CHR=NC_000913.3 START=190000 END=195000"

# ===== Tool checks =====
tools:
	@command -v bwa       >/dev/null 2>&1 || { echo "ERROR: bwa not found"; exit 1; }
	@command -v samtools  >/dev/null 2>&1 || { echo "ERROR: samtools not found"; exit 1; }
	@command -v wget      >/dev/null 2>&1 || command -v curl >/dev/null 2>&1 || { echo "ERROR: need wget or curl"; exit 1; }
	@true

# ===== Genome download + faidx =====
$(GENOME_FA):
	@mkdir -p $(GENOME_DIR)
	@echo "[-->] Fetching reference genome $(ACC) ..."
	@if command -v efetch >/dev/null 2>&1; then \
	  efetch -db nuccore -id $(ACC) -format fasta > $@; \
	else \
	  echo "[i] efetch not found; using NCBI E-utilities via curl."; \
	  curl -sSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$(ACC)&rettype=fasta&retmode=text" > $@; \
	fi
	@test -s $@ || { echo "ERROR: Genome download empty. Check ACC ($(ACC))."; exit 1; }

$(GENOME_FAI): $(GENOME_FA)
	@echo "[-->] samtools faidx $(ACC) ..."
	samtools faidx $(GENOME_FA)

genome: $(GENOME_FAI)

# ===== ENA manifest (best-effort; SRA fallback if empty) =====
$(READS_MANIFEST):
	@mkdir -p $(TMP_DIR)
	@echo "[-->] Querying ENA filereport for $(SRR) ..."
	( if command -v wget >/dev/null 2>&1; then \
	    wget -qO- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$(SRR)&result=read_run&fields=$(READS_FIELDS)&format=tsv" || true; \
	  else \
	    curl -sS --fail "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$(SRR)&result=read_run&fields=$(READS_FIELDS)&format=tsv" || true; \
	  fi ) > $@.tmp || true
	if grep -q "^fastq_https" $@.tmp; then \
	  mv $@.tmp $@; \
	  echo "[i] ENA manifest saved to $@"; \
	else \
	  echo "[i] ENA manifest unavailable or malformed; will use SRA Toolkit fallback."; \
	  echo -e "fastq_https\tfastq_ftp" > $@; \
	  rm -f $@.tmp; \
	fi

# ===== Read download: HTTPS/FTP -> fasterq-dump fallback =====
reads: $(READS_DIR)/.ok
$(READS_DIR)/.ok: $(READS_MANIFEST)
	@mkdir -p $(READS_DIR)
	@echo "[-->] Downloading fastq for $(SRR) ..."
	URLS=$$(tail -n +2 $(READS_MANIFEST) \
	  | awk -F'\t' '{print $$1" "$$2}' \
	  | sed 's/;/ /g' | tr ' ' '\n' | sed '/^$$/d'); \
	DL=0; \
	if [ -n "$$URLS" ]; then \
	  for U in $$URLS; do \
	    case "$$U" in http*|ftp*) FULL="$$U" ;; *) FULL="https://$$U" ;; esac; \
	    FNAME=$$(basename "$$FULL"); OUT="$(READS_DIR)/$$FNAME"; \
	    if [ ! -s "$$OUT" ]; then \
	      echo "   - $$FNAME"; \
	      if command -v wget >/dev/null 2>&1; then \
	        wget --quiet -O "$$OUT" "$$FULL" && DL=1 || echo "[i] wget failed for $$FULL"; \
	      else \
	        curl -sSL -o "$$OUT" "$$FULL" && DL=1 || echo "[i] curl failed for $$FULL"; \
	      fi; \
	    else DL=1; fi; \
	    if [ -s "$$OUT" ]; then gzip -t "$$OUT" || { echo "ERROR: Corrupt $$OUT"; exit 1; }; fi; \
	  done; \
	fi; \
	if [ "$$DL" -eq 0 ]; then \
	  echo "[i] Falling back to SRA Toolkit (fasterq-dump) for $(SRR) ..."; \
	  command -v fasterq-dump >/dev/null 2>&1 || { echo "ERROR: fasterq-dump not found (conda install -c bioconda sra-tools)"; exit 1; }; \
	  fasterq-dump --split-files -O $(READS_DIR) --threads 4 $(SRR); \
	  if ls $(READS_DIR)/$(SRR)_1.fastq 1>/dev/null 2>&1; then \
	    gzip -f $(READS_DIR)/$(SRR)_1.fastq $(READS_DIR)/$(SRR)_2.fastq; \
	  elif ls $(READS_DIR)/$(SRR).fastq 1>/dev/null 2>&1; then \
	    gzip -f $(READS_DIR)/$(SRR).fastq; \
	  else \
	    echo "ERROR: SRA Toolkit did not produce FASTQ files for $(SRR)."; exit 1; \
	  fi; \
	fi
	@touch $@

# ===== Create normalized symlinks (R1/R2/RS) =====
links: $(TMP_DIR)/.links_ok
$(TMP_DIR)/.links_ok: $(READS_DIR)/.ok
	@mkdir -p $(TMP_DIR)
	@echo "[-->] Normalizing read filenames for $(SRR) ..."
	PE1=$$(ls $(READS_DIR)/*_1.fastq.gz 2>/dev/null | head -n1 || true); \
	PE2=$$(ls $(READS_DIR)/*_2.fastq.gz 2>/dev/null | head -n1 || true); \
	SE=$$(ls $(READS_DIR)/*.fastq.gz 2>/dev/null | grep -v "_[12]\.fastq\.gz$$" | head -n1 || true); \
	if [ -n "$$PE1" ] && [ -n "$$PE2" ]; then \
	  ln -sf "$$PE1" "$(R1)"; ln -sf "$$PE2" "$(R2)"; echo "   - Paired-end detected."; \
	elif [ -n "$$SE" ]; then \
	  ln -sf "$$SE" "$(RS)"; echo "   - Single-end detected."; \
	else \
	  echo "ERROR: Could not classify reads as paired or single."; exit 1; \
	fi
	@touch $@

# ===== BWA index =====
$(GENOME_BWA_INDEX): $(GENOME_FA)
	@echo "[-->] Indexing reference genome with BWA ..."
	bwa index $(GENOME_FA)

index: $(GENOME_BWA_INDEX)

# ===== Alignment (produces BAM+BAI) =====
$(BAM): $(GENOME_BWA_INDEX) $(GENOME_FAI) links
	@mkdir -p $(RESULTS_DIR)
	@echo "[-->] Aligning $(SRR) ..."
	if [ -e "$(R1)" ] && [ -e "$(R2)" ]; then \
	  echo "   - Paired-end alignment"; \
	  bwa mem -t 4 $(GENOME_FA) "$(R1)" "$(R2)" | samtools view -b - | samtools sort -o "$(BAM)"; \
	elif [ -e "$(RS)" ]; then \
	  echo "   - Single-end alignment"; \
	  bwa mem -t 4 $(GENOME_FA) "$(RS)" | samtools view -b - | samtools sort -o "$(BAM)"; \
	else \
	  echo "ERROR: No normalized read links found (links target)."; exit 1; \
	fi
	@echo "[-->] Indexing BAM ..."
	samtools index "$(BAM)"

align: $(BAM)

# ===== BigWig creation =====
# Route A: deepTools (bamCoverage)
$(BW): $(BAM) $(GENOME_FAI)
ifeq ($(BW_METHOD),deeptools)
	@echo "[-->] Generating BigWig with deepTools (binSize=$(BIN_SIZE)) ..."
	@if command -v bamCoverage >/dev/null 2>&1; then \
	  bamCoverage -b "$(BAM)" -o "$(BW)" --binSize $(BIN_SIZE); \
	else \
	  echo "ERROR: bamCoverage not found. Set BW_METHOD=bedtools or install deepTools."; exit 1; \
	fi
else ifeq ($(BW_METHOD),bedtools)
	@$(MAKE) bedgraph BW=$(BW) BEDGRAPH=$(BEDGRAPH)
else
	# auto: prefer deepTools; else bedtools; else warn
	@if command -v bamCoverage >/dev/null 2>&1; then \
	  echo "[-->] Generating BigWig with deepTools (auto) ..."; \
	  bamCoverage -b "$(BAM)" -o "$(BW)" --binSize $(BIN_SIZE); \
	elif command -v bedtools >/dev/null 2>&1 && command -v bedGraphToBigWig >/dev/null 2>&1; then \
	  echo "[-->] Generating BigWig with bedtools/bedGraphToBigWig (auto) ..."; \
	  $(MAKE) bedgraph BW="$(BW)" BEDGRAPH="$(BEDGRAPH)"; \
	else \
	  echo "WARNING: Neither deepTools nor (bedtools+bedGraphToBigWig) found. Skipping BigWig."; \
	  : > "$(BW)"; \
	fi
endif

# Route B helper: bedtools -> bedGraphToBigWig (professor’s method)
.PHONY: bedgraph
bedgraph: $(BAM) $(GENOME_FAI)
	@command -v bedtools        >/dev/null 2>&1 || { echo "ERROR: bedtools not found"; exit 1; }
	@command -v bedGraphToBigWig>/dev/null 2>&1 || { echo "ERROR: bedGraphToBigWig not found"; exit 1; }
	@echo "[-->] Generating bedGraph and BigWig via bedtools ..."
	LC_ALL=C; bedtools genomecov -ibam "$(BAM)" -split -bg \
	  | sort -k1,1 -k2,2n > "$(BEDGRAPH)"
	bedGraphToBigWig "$(BEDGRAPH)" "$(GENOME_FAI)" "$(BW)"

bigwig: $(BW)

# ===== Core stats (for report) =====
$(STATS): $(BAM)
	@echo "[-->] Computing alignment stats (flagstat + idxstats) ..."
	{ \
	  echo "--- Alignment Report for SRR: $(SRR) ---"; \
	  samtools flagstat "$(BAM)"; \
	  echo; echo "--- Per-reference read counts (idxstats) ---"; \
	  samtools idxstats "$(BAM)"; \
	} > "$(STATS)"
	@cat "$(STATS)"

stats: $(STATS)

# ===== Helpers for the specific questions =====
# Primary alignments (exclude secondary 0x100 and supplementary 0x800)
$(PRIMARY): $(BAM)
	@echo "[-->] Counting PRIMARY alignments ..."
	@echo -n "primary_alignments: " > "$(PRIMARY)"
	@samtools view -c -F 0x900 "$(BAM)" >> "$(PRIMARY)"
	@cat "$(PRIMARY)"

primary: $(PRIMARY)

# Max observed coverage (position with largest depth)
$(MAXDEPTH): $(BAM)
	@echo "[-->] Finding position with MAX depth ..."
	@samtools depth -a "$(BAM)" | sort -k3,3nr | head -1 > "$(MAXDEPTH)" || true
	@if [ ! -s "$(MAXDEPTH)" ]; then echo "NA\t0\t0" > "$(MAXDEPTH)"; fi
	@echo -n "max_depth (chrom pos depth): "; cat "$(MAXDEPTH)"

maxdepth: $(MAXDEPTH)

# Forward-strand alignments covering a region (your chosen gene window)
$(GENEFWD): $(BAM)
	@if [ -z "$(CHR)" ] || [ -z "$(START)" ] || [ -z "$(END)" ]; then \
	  echo "ERROR: Provide CHR, START, END (e.g., CHR=NC_000913.3 START=190000 END=195000)"; exit 1; \
	fi
	@echo "[-->] Counting forward-strand alignments in $(CHR):$(START)-$(END) ..."
	# Forward strand => flag 0x10 NOT set, hence -F 0x10
	@echo -n "$(CHR):$(START)-$(END) forward_alignments: " > "$(GENEFWD)"
	@samtools view -c -F 0x10 "$(BAM)" "$(CHR):$(START)-$(END)" >> "$(GENEFWD)"
	@cat "$(GENEFWD)"

gene_forward: $(GENEFWD)

# ===== Utility =====
clean:
	@echo "[-->] Cleaning generated outputs ..."
	rm -rf "$(GENOME_DIR)" "$(READS_DIR)" "$(RESULTS_DIR)" "$(TMP_DIR)"

# --- override: disable normalize step ---
.tmp/.links_ok:
	@mkdir -p .tmp
	@echo "[skip] Normalizing read filenames disabled."
	@touch $$@
