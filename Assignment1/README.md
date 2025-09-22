# Assignment 1 – Setup and Unix

## 1. Samtools Version in bioinfo Environment
```bash
conda activate bioinfo
samtools --version
```
Result: samtools 1.17
Using htslib 1.17

## 2. Create Nested Directory Structure
```bash
mkdir -p project/data/raw project/data/processed project/results/figures
```
 
 ## 3. Create Files in Different Directories
 ```bash
 echo "raw data file" > project/data/raw/sample1.txt
echo "processed data file" > project/data/processed/sample1_processed.txt
echo "figure placeholder" > project/results/figures/figure1.txt
```

## 4.Access Files Using Relative Paths
``` bash
cat project/data/raw/sample1.txt
cat project/results/figures/figure1.txt
```

## 5. Access Files Using Absolute Paths
``` bash
cat /mnt/d/BMMB852_F2025/01_Assignment1/project/data/processed/sample1_processed.txt
``` 

