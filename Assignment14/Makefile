# ============================================================
# Makefile: RNA-Seq Differential Expression & Enrichment
# ============================================================

# --- Standard Shell Settings (Safety First) ---
SHELL := /usr/bin/env bash
.SHELLFLAGS := -euo pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# --- Variables ---
# Input Files
COUNTS  := counts.txt
DESIGN  := design.csv
SCRIPT  := scripts/analyze.R

# Output Directory
OUTDIR  := results

# Expected Output Files
DE_RES      := $(OUTDIR)/edger_results.csv
PCA_PDF     := $(OUTDIR)/pca.pdf
HEATMAP_PDF := $(OUTDIR)/heatmap.pdf
ENRICH_CSV  := $(OUTDIR)/enrichment.csv

# --- Targets ---
.PHONY: all usage de pca heatmap enrichment clean

# Target: all (Runs everything)
all: de pca heatmap enrichment

# Target: usage (Displays help menu like Example 2)
usage:
	@echo "# =========================================="
	@echo "# RNA-Seq Analysis Pipeline"
	@echo "# =========================================="
	@echo "# usage: make [target]"
	@echo "#"
	@echo "# Targets:"
	@echo "#   all        Run the full pipeline"
	@echo "#   de         Perform differential expression"
	@echo "#   pca        Generate PCA plot"
	@echo "#   heatmap    Generate Heatmap"
	@echo "#   enrichment Run functional enrichment"
	@echo "#   clean      Remove all results"
	@echo "# =========================================="

# --- Analysis Rules ---

# Note: Since your R script (scripts/analyze.R) performs all steps
# in one go, these targets act as checkpoints.

# 1. Differential Expression
de: $(DE_RES)

$(DE_RES): $(COUNTS) $(DESIGN) $(SCRIPT)
	@echo "Running Differential Expression..."
	@mkdir -p $(OUTDIR)
	Rscript $(SCRIPT) $(COUNTS) $(DESIGN) $(OUTDIR)

# 2. PCA Plot
pca: $(PCA_PDF)

$(PCA_PDF): $(DE_RES)
	@echo "Checking for PCA plot..."
	@if [ ! -f $@ ]; then echo "Error: PCA plot not generated."; exit 1; fi

# 3. Heatmap
heatmap: $(HEATMAP_PDF)

$(HEATMAP_PDF): $(DE_RES)
	@echo "Checking for Heatmap..."
	@if [ ! -f $@ ]; then echo "Error: Heatmap not generated."; exit 1; fi

# 4. Functional Enrichment
enrichment: $(ENRICH_CSV)

$(ENRICH_CSV): $(DE_RES)
	@echo "Checking for Enrichment results..."
	@if [ ! -f $@ ]; then echo "Enrichment skipped (no significant genes found) or failed."; fi

# --- Clean Up ---
clean:
	@echo "Cleaning up results..."
	rm -rf $(OUTDIR)