# Assignment 12: Evaluate data from the Cancer Genome in a Bottle project
**Name:** Aaliya\
**Date:** 2025-11-14

---
##  Usage

```bash
# Run full analysis from scratch 
make all

# Or run step-by-step:
make ref        # Download reference genome
make bam        # Download tumor + normal BAMs
make vcf        # Variant calling
make somatic    # Find tumor-only variants
make filter     # Apply variant filters
make gold       # Download & subset benchmark VCF
make compare    # Evaluate pipeline vs gold standard
make summary    # Generate final performance report

Copy code
make clean       # Remove VCFs, comparisons, and report
make clean-all   # Delete everything (BAMs, reference, etc.)
```

Project Overview\
This pipeline compares tumor and normal samples to identify somatic mutations, focusing on a 2 Mb region surrounding TP53 on chromosome 17.

TP53 is a well-studied tumor suppressor gene, frequently mutated in cancers including breast, colorectal, and lung cancers.

Region of Analysis
Chromosome: chr17

Coordinates: 7,000,000–9,000,000

Gene of Interest: TP53

Reference Genome: GRCh38 (GIAB v3, no-alt)

Data Sources\
Sample Type	Sample Name	Coverage	Source\
Normal	HG008-N-D	~77×	GIAB (Element AVITI)\
Tumor	HG008-T	~111×	GIAB (Element AVITI)\
Benchmark	DeepVariant	Validated VCF	GIAB (Ultima Genomics)

Pipeline Summary
```bash
Download reference genome and alignment BAMs

Call variants for normal and tumor separately

Find tumor-only variants (somatic candidates)

Filter variants using relaxed criteria

Compare filtered calls against the benchmark
🧾 Report performance metrics (TP, FP, FN, precision, recall)
```

Filtering Criteria
Final filter used:
``` bash
QUAL > 10 & INFO/DP > 5
```
Rationale

I tested several filters:

Strict: QUAL > 25, DP > 10 →  no variants passed

Moderate: QUAL > 20, DP > 8 → very few true positives

Final choice: QUAL > 10 & INFO/DP > 5 offered the best trade-off for:

Detecting more somatic variants (increased recall)

Allowing benchmark comparison

Highlighting false positive rate of bcftools


Example Report
Generated with make summary:
```bash
Region: chr17:7000000-9000000
Filter: QUAL > 10 & INFO/DP > 5

True Positives:    0
False Positives: 2630
False Negatives:   6

Precision: 0.00%
Recall:    0.00%
```
Indicates that:

bcftools produced many low-confidence calls

None matched the validated gold standard

A machine learning caller (e.g., Mutect2, DeepVariant) would improve results

Ouput Metrics
```bash
### Performance Metrics Summary


| Metric               | Value | Interpretation                                            
| True Positives (TP)  |   0   | No variants called by `bcftools` matched the benchmark set.                      |
| False Positives (FP) | 2630  | Many variants were called, but none were validated — suggesting high noise.      |
| False Negatives (FN) |   6   | All 6 validated somatic variants from DeepVariant were missed.                   |
| Precision        |  0%   | None of the variants we called were correct.                                     |
| Recall           |  0%   | We missed 100% of the validated variants.                                        |
| F1 Score         |  0%   | No balance between precision and recall — performance was completely off.        |

```
False Positives (FP): 2630
Our simple pipeline identified 2,630 potential variants that are all noise. These are likely alignment errors or sequencing artifacts that look like variants but are not real. This results in a Precision of 0.00% (none of our calls were correct).

False Negatives (FN): 6
There were 6 real somatic variants in this region according to the gold standard, and our pipeline missed every single one. They were likely too subtle (e.g., low allele frequency) for our simple bcftools call to detect. This results in a Recall of 0.00% (we found none of the real variants).

Conclusion: This pipeline successfully highlights the challenges of somatic variant calling. It proves that a simple bcftools pipeline with lenient filters is not sufficient to distinguish real, low-frequency somatic mutations from background noise, leading to thousands of false positives and a failure to detect the true variants.

## My Thoughts 

```bash
At first, I thought the 0% precision was an error in my Makefile, but it was actually the correct result.

It's clear that a simple bcftools pipeline is not reliable for this kind of analysis, as it generated over 2,600 false positives while missing all 6 true variants.

This result shows why advanced tools like DeepVariant are necessary. They are much better at finding the "needle in the haystack" (the 6 true variants) and ignoring the "hay" (the 2,630 noisy artifacts).

If I were to improve this, the next step would be to add much stricter filters...
```

My File Structure
```bash
Assignment12/
├── Makefile               # Pipeline logic
├── README.md              # Project documentation
├── data/                  # Reference genome + BAMs
├── vcfs/                  # VCF files (raw, filtered)
├── evals/                 # Evaluation results (TP, FP, FN)
└── reports/               # Summary output
Tools Used
bcftools (variant calling, filtering, comparison)

samtools (indexing, subsetting BAMs)

tabix (VCF indexing)

make (workflow automation)

awk, bash (report generation)
```
# Extended Experiment:Relaxed Pipeline

Changes in relaxed analysis:

No filters applied and Smaller TP53‑focused window. 
```bash
chr17:7650000–7680000
```
Relaxed Pipeline Usage
```bash
make relaxed-all
```
Cleanup:
```bash
make relaxed-clean
make relaxed-clean-all
```
## Results
```bash
True Positives:  0
False Positives: 46
False Negatives: 1

Precision: 0.00%
Recall:    0.00%
```
## File Structue 
```bash
Assignment12/
├── Makefile               # Full pipeline (strict + relaxed)
├── README.md              # This file
├── data/                  # BAMs and reference
├── results/               # Strict filtered VCFs and comparison
├── report/                # Summary metrics (strict)
├── relaxed_data/          # Relaxed BAMs and reference
├── relaxed_vcf/           # Raw VCFs (relaxed)
├── relaxed_output/        # Relaxed metrics and comparison
```

### After getting 0 true positives and 6 false negatives in my original pipeline, I realized my filters were likely too strict for detecting low-frequency somatic mutations. So, I created a relaxed version that: 
used a smaller region focused around known TP53 variants

Removed all filtering (no QUAL or DP thresholds)

I hoped this would at least recover some true positives, even at the cost of more false positives.

What I got:

1. Fewer false positives (only 46, due to smaller region)

2. Still 0 true positives — even without filters

3. 1 false negative (gold standard missed)