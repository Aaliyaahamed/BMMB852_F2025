# Assignment 11: Variant Effect Prediction

**Name:** Aaliya\
**Date:** 2025-11-07

---

## 1. Process Summary

For this assignment, I extended my previous `Makefile` pipeline to include variant effect prediction. The goal was to generate a VCF file from our Ebola samples and then evaluate the functional effects of at least three different variants using a software tool.

Here is the process I followed:

1.  **Pipeline Setup:** I used a `Makefile` to automate the entire process. This pipeline first ensures all reference files (FASTA) and tools (`bwa`, `samtools`, `bcftools`, `htslib-tools`) are in place.
2.  **Sample Processing:** The `call_all` target looped through my `design.csv` file, processing 10 samples. For each sample, it downloaded 10,000 reads from the SRA using `fastq-dump`.
3.  **Alignment & Variant Calling:** Each sample's reads were aligned to the Ebola reference genome (`KJ660346.2`) using `bwa mem`. The resulting alignments were sorted (`samtools sort`) and indexed. Variants were called using `bcftools mpileup` and `bcftools call`.
4.  **Merging:** All 10 single-sample VCF files were merged into one master file: `all_samples.merged.vcf.gz`.
5.  **Annotation with snpEff:** This is where the new steps occurred.
    * To make `snpEff` work, I had to create a special copy of the VCF. The `snpEff` database (`ebola_zaire`) expected the chromosome to be named `KJ660346.1`, but my VCF used `KJ660346.2`.
    * The `Makefile` automatically created a "map file" (`chrom_map.txt`) to handle this.
    * It then used `bcftools annotate --rename-chrs` to create a new, properly-named VCF: `all_samples.for_snpeff.vcf.gz`.
    * Finally, it ran `snpEff` on this new file, which successfully generated the `snpEff_report.html` (for summary stats) and the final annotated VCF, `all_samples.ann.vcf`.
6.  **Analysis:** I analyzed the output from `snpEff` to find variants of different impacts.

## 2. Variant Analysis Results

The `snpEff` analysis of the 9 variants found in the merged VCF file revealed a total of 43 functional annotations. These were categorized by impact (based on my `snpEff_report.html`):

* **MODERATE:** 5 variants (11.6%)
* **LOW:** 5 variants (11.6%)
* **MODIFIER:** 33 variants (76.7%)

The 10 coding annotations were evenly split between 5 **missense** (Moderate) and 5 **synonymous** (Low) variants. I selected three of these from my `all_samples.ann.vcf` file for my final analysis.

---

### Variant 1: MODERATE Impact (Missense)

* **Gene:** `NP` (nucleoprotein)
* **Position:** 800
* **Change:** `C > T`
* **Effect:** `p.Arg111Cys` (Arginine to Cysteine)
* **Analysis:** This is a **missense variant**, which means the nucleotide change resulted in a *different* amino acid being coded. An Arginine (a basic, positively charged amino acid) was replaced by a Cysteine (a polar, uncharged amino acid). This change in chemical properties could potentially affect the protein's folding or function.

### Variant 2: LOW Impact (Synonymous)

* **Gene:** `NP` (nucleoprotein)
* **Position:** 1849
* **Change:** `T > C`
* **Effect:** `p.Asp460Asp` (Aspartic acid to Aspartic acid)
* **Analysis:** This is a **synonymous variant** (also called a "silent" mutation). Although the nucleotide changed from T to C, the resulting codon still codes for the same amino acid (Aspartic acid). This mutation has a "LOW" impact because it does not change the final protein sequence.

### Variant 3: MODERATE Impact (Different Gene)

* **Gene:** `GP` (glycoprotein)
* **Position:** 6283
* **Change:** `C > T`
* **Effect:** `p.Ala82Val` (Alanine to Valine)
* **Analysis:** This is another **missense variant**. It changes an Alanine to a Valine. Both are nonpolar, hydrophobic amino acids, but Valine is slightly larger. This change is less drastic than the first example but still modifies the protein sequence and is therefore considered a "MODERATE" impact.


Understanding My Final Outputs:
![alt text](image-1.png)
![alt text](image.png)
![alt text](image-2.png)

### Commands I ran on terminal
**Main Command:**

Just run `make call_all` to do everything from start to finish.

```bash
# This runs the entire pipeline
make call_all
```

**How I Wrote My Report (My Analysis)**\
To finish the assignment, I used the files my pipeline created.

I opened my results/snpEff_report.html file to see the summary charts (like 5 MODERATE, 5 LOW impact).

I used a command to pull the variant data from the text-based VCF file:

```Bash
cat results/all_samples.ann.vcf | grep -E "missense_variant|synonymous_variant"
```