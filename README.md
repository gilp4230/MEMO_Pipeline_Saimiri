# MEMO_Pipeline_Saimiri
## 🧬 Saimiri VCF Filtering Workflow

This repository provides a complete workflow for **VCF quality control and hard filtering** of *Saimiri boliviensis* genomic data.  
The protocol follows best practices inspired by **GATK**, **vcftools**, and **bcftools** recommendations, as well as criteria from  
[Núria et al., 2024](https://doi.org/10.1038/s42003-024-06901-3).

The main guide, [`Saimiri_vcf_filtering_workflow.md`](./Saimiri_vcf_filtering_workflow.md), details each filtering step — from raw unfiltered VCFs to a finalized high-quality dataset ready for population genomic analyses (e.g., PCA, admixture, structure).

---

### 📋 Contents

- **Module loading and setup**
- **Hard filtering with bcftools and GATK**
- **Depth and missingness evaluation**
- **VCFtools-based filtering**
- **Allele frequency and heterozygosity statistics**
- **Final QC and visualization**

---

## 🧬 Saimiri Population Structure Analysis
