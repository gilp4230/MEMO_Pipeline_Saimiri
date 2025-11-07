# 🧬 Hard Filtering Process — Introduction and Details

## Information from the Sequencing Lab

1. The dataset may include **related** or **admixed** individuals.  
2. **Individuals with QC issues** (coverage, missing data, allele balance) were **not removed**.  
   - You can find this information in the provided *list of samples* and remove them if desired.  
   - **Note:** There are **no Saimiri** individuals in that list.  
3. The lab performed its **internal QC**, but errors may still exist.  
   - You should run **your own QC** before analyses and report any issues found.  
4. The **VCF files are unfiltered**.  
   - Apply your own filtering criteria depending on your analysis goals.

---

## Reference

**Nuria et al., 2024**  
📄 DOI: [10.1038/s42003-024-06901-3](https://doi.org/10.1038/s42003-024-06901-3)

---

## Hard Filtering Criteria (from the paper)

Hard filtering was applied to the resulting SNPs in each chunk according to the following steps:

- **Removed non-biallelic variants**  
  *Tool:* `vcflib vcfbiallelic` (v1.0.5)

- **Excluded variants based on depth (DP):**
  - Minimum: **5×** (one third of the median minimum coverage across samples)
  - Maximum: **70×** (double the median maximum coverage across samples)

- **Removed variants with >60% missingness**  
  *Tool:* `vcftools` (v0.1.12)

- **Applied GATK Best Practices filtering** using `GATK VariantFiltration` with the following expression: QD < 2 | FS > 60 | MQ < 40 | SOR > 3 | ReadPosRankSum < -8.0 | MQRankSum < -12.5

- **Merged** the remaining SNPs from each chunk into a single dataset.

- **Filtered for allele imbalance**, keeping variants with allele frequencies between **0.25–0.75**  
*(see Supplementary Methods, Section 1.2; Fig. S1)*

- **Removed SNPs** in:
- Scaffolds belonging to **sex chromosomes** (identified in *Kuderna et al.*)
- Scaffolds **shorter than 0.5 Mb**

---

✅ **Final dataset:** 

---

# 🧬 B. Basic (useful) Command Lines

## 1. Module Loading
```bash
module load BCFtools/1.14-GCC-11.2.0
module load VCFtools/0.1.16-GCC-12.3.0
```
## 2. Number of Variants
```bash
bcftools view -H {file-name} | wc -l
```
## 3. Creating Sample Files
```bash
bcftools query -l {input-file-name} | sort | uniq > Analysis_BCF/{output-file-name}.txt
```
## 4. Create Statistics to Compare Files
```bash
bcftools stats {input-file-name} > {output-file-name}.txt
```
---

# 🧬 C. BCFtools Analysis

## 1. Number of Variants

```bash
bcftools view -H Saimiri_boliviensis__Saimiri.pop.geno.snps.raw.reheader.vcf.gz | wc -l
```
→ 80,590,885 variants
→ 59 individuals

## 2. General Statistics from VCF

```bash
bcftools stats Saimiri_boliviensis__Saimiri.pop.geno.snps.raw.reheader.vcf.gz > Analysis_BCF/Saimiri_BCF_stats.txt
```
---

# 🧬 D. Removing Non-Biallelic SNPs

## 1. Remove Non-Biallelic SNPs

```bash
nohup bash -c 'bcftools view -m2 -M2 -v snps \
-o Saimiri_boliviensis__Saimiri_filtered_biallelic.vcf.gz \
Saimiri_boliviensis__Saimiri.pop.geno.snps.raw.reheader.vcf.gz' > biallelic_log.out 2>&1 &
```
## 2. Check If Only Biallelic SNPs Remain

```bash
bcftools view -m2 -M2 -v snps Saimiri_boliviensis__Saimiri.pop.geno.snps.raw.reheader.vcf.gz | wc -l
```
→ 74,698,589 variants

## 3. Indexing the Output File

```bash
bcftools index Saimiri_boliviensis__Saimiri_filtered_biallelic.vcf.gz
```
---

# 🧬 E. GATK Hard Filtering and Quality Control

This section applies GATK-based hard filters to remove low-quality variants, indexes files, and assesses the results.

## 1. Preparing Files for GATK

Before applying filters, ensure that:
- File names are correct.
- An index file exists in the format recognized by GATK.

### 1.1 Load Tabix

```bash
module load tabix/0.2.6-GCCcore-10.2.0
```
### 1.2 Index the VCF File for GATK

```bash
tabix -p vcf Saimiri_boliviensis__Saimiri_filtered_biallelic.vcf.gz
```

## 2. Load GATK Module

```bash
module load GATK/4.2.4.1-GCCcore-10.2.0-Java-1.8
```
## 3. Apply Hard Filtering with GATK

```bash
nohup gatk VariantFiltration \
  -V Saimiri_boliviensis__Saimiri_filtered_biallelic.vcf.gz \
  -O Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || ReadPosRankSum < -8.0 || MQRankSum < -12.5" \
  --filter-name "HardFilter" > HF_log_out 2>&1 &
```

## 4. Index the GATK-Filtered File

```bash
bcftools index Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz
```
## 5. Count Variants Before and After Filtering

```bash
bcftools view -H Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz | wc -l
```
→ Before filtering: xxx xxx xxx variants
→ After filtering: 74,698,589 variants

## 6. Inspect the FILTER Column

```bash
bcftools view -H Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz | cut -f 7 | sort | uniq -c
```
## 7. Remove Failed Variants and Re-Index

```bash 
bcftools view -f PASS -Oz -o Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz \
Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz
```
### 7.1 Index the Final File

```bash
tabix -p vcf Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz

bcftools index Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz
```
## 8. Check Missingness per Individual

```bash
vcftools --gzvcf Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz \
  --missing-indv \
  --out missingness
```
---

# 🧬 F. Checking and Filtering Missing Data per Sample

This section assesses missing genotype data per individual and, if necessary, removes variants or samples exceeding the allowed missingness threshold.

## 1. Assess Missingness with VCFtools

```bash
vcftools --gzvcf Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz --missing-indv --out missingness
```
### 1.1 Optional — Filter Individuals or Variants by Missingness

```bash
nohup bash -c "bcftools filter -i 'F_MISSING <= 0.6' -Oz -o Saimiri_boliviensis__Saimiri_filtered_final_missingindv.vcf.gz \
Saimiri_boliviensis__Saimiri.pop.geno.snps.raw.reheader.vcf.gz" > bcf_missing.out 2>&1 &
```
### 1.2 Filter Variants with >60% Missing Data (Alternative Example)

```bash
# Using vcftools (alternative approach)
nohup bash -c 'vcftools --vcf Pithecia_depth_filtered.vcf.gz --max-missing 0.4 --recode --out Pithecia_missing_filtered' > missing_log.out 2>&1 &

# Using bcftools (preferred for compressed VCFs)
nohup bash -c "bcftools filter -i 'F_MISSING <= 0.6' -Oz -o Pithecia_missing_filtered.vcf.gz \
Pithecia_depth_filtered.vcf.gz" > bcf_missing.out 2>&1 &
```

### 1.3 Count Variants After Filtering

```bash
bcftools view -H Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz | wc -l
```
→ Before filtering : 74,698,589 variants
→ After filetring : 61,993,478 variants

---

# 🧬 G. VCFtools Analysis and Depth-Based Filtering

Reference: [Filtering VCFs — Speciation Genomics](https://speciationgenomics.github.io/filtering_vcfs/)

This section calculates key quality metrics (allele frequencies, depth, missingness, heterozygosity, etc.), determines coverage thresholds, and applies hard filters using **VCFtools** and **BCFtools**.

## 1. Create VCF Variables

> ⚠️ Run these commands only on **biallelic VCF files**.

```bash
# Define input and output paths
INPUT=~/Primates/memo_saimiri/Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz
OUTPUT=~/Primates/memo_saimiri/Analysis_VCF/Saimiri_stats_VCF
```

### 1.1 Calculate Key Statistics

```bash
# Allele frequency
vcftools --gzvcf $INPUT --freq2 --out $OUTPUT --max-alleles 2

# Mean depth per individual
vcftools --gzvcf $INPUT --depth --out $OUTPUT

# Mean depth per site
vcftools --gzvcf $INPUT --site-mean-depth --out $OUTPUT

# Site quality
vcftools --gzvcf $INPUT --site-quality --out $OUTPUT

# Proportion of missing data per individual
vcftools --gzvcf $INPUT --missing-indv --out $OUTPUT

# Proportion of missing data per site
vcftools --gzvcf $INPUT --missing-site --out $OUTPUT

# Heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf $INPUT --het --out $OUTPUT
```

## 2. Calculating Depth Distribution
### 2.1 Extract Per-Sample Depths
```bash
bcftools query -f '[%DP\n]' Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz | grep -v '^0$' > Analysis_BCF/all_depths.txt
```

### 2.2 Clean the Depth File

```bash
awk '{for (i=1; i<=NF; i++) if ($i ~ /^[0-9]+$/) printf "%s ", $i; print ""}' all_depths.txt > Analysis_BCF/all_depths_cleaned.txt
```

### 2.3 Calculate the Per-Sample Median Depth

```bash
sort all_depths_cleaned.txt | uniq -c | sort -nr | head -1
```
### 2.4 Compute Minimum and Maximum Thresholds

```bash
mode=35  # replace with your actual mode value
min_thresh=$(echo "$mode / 3" | bc)  # one-third of the mode
max_thresh=$(echo "$mode * 2" | bc)  # double the mode
echo "Min DP: $min_thresh"
echo "Max DP: $max_thresh"
```

## 3. Setting Hard-Filter Parameters

```bash
MISS=0.75      # maximum allowed missingness (75%)
QUAL=30        # minimum site quality
MIN_DEPTH=11   # minimum depth threshold
MAX_DEPTH=70   # maximum depth threshold
```
## 4. Perform Filtering with VCFtools

```bash
vcftools --gzvcf Saimiri_boliviensis__Saimiri_biallelic.vcf.gz \
  --remove-indels \
  --max-missing $MISS \
  --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
  --minDP $MIN_DEPTH --maxDP $MAX_DEPTH \
  --recode --stdout | gzip -c > Saimiri_filtered.vcf.gz
```
### 4.1 

```bash
bcftools view -H Saimiri_filtered.vcf.gz | wc -l
```
→ Before filtering : xxx xxx xxx variants
→ After filetring : 67 420 702 variants

## 5. Visualize Depth Distribution (Python)

```bash
### PYTHON SCRIPT ###
import numpy as np
import matplotlib.pyplot as plt

# Load depths from file
depths = []
with open("all_depths_cleaned.txt") as f:
    for line in f:
        try:
            dp = int(line.strip())
            depths.append(dp)
        except:
            continue

depths = np.array(depths)

# Calculate statistics
mode = int(np.bincount(depths).argmax())
p95 = int(np.percentile(depths, 95))
print(f"Mode depth: {mode}")
print(f"95th percentile depth: {p95}")

# Plot histogram
plt.figure(figsize=(10, 6))
plt.hist(depths, bins=range(0, min(200, depths.max()+1)), color="skyblue", edgecolor="black")
plt.axvline(mode, color='green', linestyle='--', label=f"Mode = {mode}")
plt.axvline(p95, color='red', linestyle='--', label=f"95th percentile = {p95}")
plt.xlabel("Depth (DP)")
plt.ylabel("Frequency")
plt.title("Read Depth Distribution")
plt.legend()
plt.tight_layout()
plt.savefig("depth_distribution.png", dpi=300)
######
```
Mode depth: 35
95th percentile depth: 44
In this study, the 95% depth valus as a criterion for the maximum DP cutoff works well

---

# 🧬 H. Scaffold and Minor Allele Frequency (MAF) Filtering

This section filters SNPs by **scaffold length** and **minor allele frequency (MAF)** to retain only high-confidence variants suitable for downstream population genetic analyses.

## 1. Filter SNPs in Scaffolds Shorter than 0.5 Mb

### 1.1 Index the VCF (if not already indexed)

```bash
bcftools index Saimiri.vcf.gz

bcftools view -H Saimiri.vcf.gz | wc -l
```
### 1.2 Extract Scaffold Names and Lengths from the VCF Header

```bash
bcftools view -h Saimiri.vcf.gz | grep "^##contig" | sed -E 's/.*ID=([^,]+),length=([0-9]+).*/\1\t\2/' > scaffold_lengths.txt
```

### 1.3 Keep Only Scaffolds ≥ 500,000 bp

```bash
awk '$2 >= 500000 {print $1}' scaffold_lengths.txt > long_scaffolds.txt
```

### 1.4 Convert the List to a Comma-Separated String

```bash
paste -sd, long_scaffolds.txt > scaffold_list.txt
LONG_SCAFFOLDS=$(cat scaffold_list.txt)
```
### 1.5 Filter the VCF to Retain Only SNPs on Long Scaffolds

```bash
bcftools view -r "$LONG_SCAFFOLDS" Saimiri.vcf.gz -Oz -o Saimiri_scaffold_filtered.vcf.gz
```
### 1.6 Index the Filtered VCF

```bash
bcftools index Saimiri_scaffold_filtered.vcf.gz

bcftools view -H Saimiri_scaffold_filtered.vcf.gz | wc -l
```

## 2. Filter SNPs by Minor Allele Frequency (MAF ≥ 0.015)

Filtering variants based on Minor Allele Frequency (MAF) removes extremely rare alleles that may represent sequencing errors or artifacts.

### 2.1 Tag the VCF with MAF Values

```bash
bcftools +fill-tags Saimiri_scaffold_filtered.vcf.gz -Oz -o Saimiri_maf_tagged.vcf.gz -- -t MAF
```
### 2.2 Apply the MAF Threshold (≥ 0.015)

```bash 
bcftools view -q 0.015:minor Saimiri_maf_tagged.vcf.gz -Oz -o Saimiri_MAF0015.vcf.gz
```
### 2.3 Index the Result

```bash
bcftools index Saimiri_MAF0015.vcf.gz

bcftools view -H Saimiri_MAF0015.vcf.gz | wc -l
```

## 3. Filter SNPs by Minor Allele Frequency (MAF ≥ 0.05)

This step applies a more stringent filtering threshold (MAF ≥ 0.05), retaining only variants with sufficient representation in the population.

### 3.1 Tag VCF with Updated MAF Values

```bash
bcftools +fill-tags Saimiri_MAF0015.vcf.gz -Oz -o Saimiri_maf_tagged.vcf.gz -- -t MAF
```
### 3.2 Apply the MAF Threshold (≥ 0.05)

```bash
bcftools view -q 0.05:minor Saimiri_maf_tagged.vcf.gz -Oz -o Saimiri_MAF005.vcf.gz
```
### 3.3 Index the Final Filtered File

```bash
bcftools index Saimiri_MAF005.vcf.gz

bcftools view -H Saimiri_MAF005.vcf.gz | wc -l
``` 

---

# 🧩 End of Filtering Pipeline

All filtering steps are now complete.
You can proceed to downstream analyses (e.g., PCA, admixture, phylogenetics) using the final filtered dataset:

```bash
Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz
```