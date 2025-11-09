# 🧬 A. Population Structure Analysis

## I. Principal Component Analysis (all samples Identified)
Script adapted from https://speciationgenomics.github.io/pca/

### 1. Module Loading

```bash
module load PLINK/2.00a3.7-gfbf-2023a
```

### 2. Prepare Input File

```bash
VCF=/srv/home/gpierre/Primates/memo_saimiri/Saimiri_filtered_final.vcf.gz
```
### 3. Perform linkage disequilibrium (LD) pruning

Removes correlated SNPs (r^2 > 0.1 in 50-SNP windows sliding by 10)

```bash
plink --vcf /srv/home/gpierre/Primates/memo_saimiri/Saimiri_filtered_final.vcf.gz \
      --keep Samples_identified.txt \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 \
      --out Identified/Saimiri_identified
```

### 4. Perform PCA using only the pruned SNPs

```bash
plink --vcf /srv/home/gpierre/Primates/memo_saimiri/Saimiri_filtered_final.vcf.gz \
      --keep Samples_identified.txt \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --extract Identified/Saimiri_identified.prune.in \
      --make-bed --pca --out Identified/Saimiri_identified
```

### 5. Plotting the PCA (all samples identified) output with R

```bash
### R SCRIPT ###
#Load necessary packages
if (!require(tidyverse)) install.packages("tidyverse", repos = "https://cloud.r-project.org")
if (!require(ggrepel)) install.packages("ggrepel", repos="https://cloud.r-project.org")
library(ggrepel)
library(tidyverse)
library(dplyr)
#Read in data
pca <- read_table("Identified/Saimiri_identified.eigenvec", col_names = FALSE)
eigenval <- scan("Identified/Saimiri_identified.eigenval")
#Sort out the PCA data
pca <- pca[,-1]  #Remove nuisance column
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
print(pca$ind)  #View all individual names
#Define species mapping
species_mapping <- c(
  #Saimiri boliviensis
  "PD_2704" = "Saimiri_boliviensis", "PD_1233" = "Saimiri_boliviensis", "PD_1232" = "Saimiri_boliviensis", "PD_2703" = "Saimiri_boliviensis", "PD_2706" = "Saimiri_boliviensis", "PD_2705" = "Saimiri_boliviensis",
  "PD_0370" = "Saimiri_boliviensis",
  #Saimiri cassiquiarensis
  "PD_1251" = "Saimiri_cassiquiarensis", "PD_0359" = "Saimiri_cassiquiarensis", "PD_0358" = "Saimiri_cassiquiarensis",
  "PD_0354" = "Saimiri_cassiquiarensis", "PD_0356" = "Saimiri_cassiquiarensis", "PD_0355" = "Saimiri_cassiquiarensis",
  "PD_0360" = "Saimiri_cassiquiarensis", "PD_0353" = "Saimiri_cassiquiarensis", "PD_0373" = "Saimiri_cassiquiarensis",
  "PD_0357" = "Saimiri_cassiquiarensis", "PD_1226" = "Saimiri_cassiquiarensis",
  #Saimiri macrodon
  "PD_1238" = "Saimiri_macrodon", "PD_0362" = "Saimiri_macrodon", "PD_1243" = "Saimiri_macrodon", "PD_1230" = "Saimiri_macrodon",
  "PD_1240" = "Saimiri_macrodon", "PD_1242" = "Saimiri_macrodon", "PD_1237" = "Saimiri_macrodon", "PD_1234" = "Saimiri_macrodon",
  "PD_0361" = "Saimiri_macrodon", "PD_0372" = "Saimiri_macrodon", "PD_1239" = "Saimiri_macrodon", "PD_1246" = "Saimiri_macrodon",
  "PD_1231" = "Saimiri_macrodon", "PD_1236" = "Saimiri_macrodon", "PD_1235" = "Saimiri_macrodon", "PD_1244" = "Saimiri_macrodon",
  "PD_1245" = "Saimiri_macrodon", "PD_1247" = "Saimiri_macrodon",
  #Saimiri oerstedii
  "PD_0176" = "Saimiri_oerstedii",
  #Saimiri sciureus
  "WPB399" = "Saimiri_sciureus", "WPB374" = "Saimiri_sciureus", "PD_1250" = "Saimiri_sciureus", "PD_1249" = "Saimiri_sciureus",
  "PD_1248" = "Saimiri_sciureus", "PD_0364" = "Saimiri_sciureus", "PD_0365" = "Saimiri_sciureus",
  #Saimiri sp
  "PD_2174" = "Saimiri_sp", "PD_1310" = "Saimiri_sp", "PD_1275" = "Saimiri_sp",
  #Saimiri ustus
  "PD_1227" = "Saimiri_ustus", "PD_0371" = "Saimiri_ustus", "PD_0375" = "Saimiri_ustus", "PD_1241" = "Saimiri_ustus",
  "PD_1229" = "Saimiri_ustus", "PD_0363" = "Saimiri_ustus", "PD_0374" = "Saimiri_ustus",
  "PD_0368" = "Saimiri_ustus", "PD_0366" = "Saimiri_ustus", "PD_0367" = "Saimiri_ustus", "PD_1225" = "Saimiri_ustus",
  "PD_1228" = "Saimiri_ustus"
)

pop_colors <- c(
  Saimiri_boliviensis     = "#E69F00",  # orange
  Saimiri_cassiquiarensis = "#56B4E9",  # sky blue
  Saimiri_macrodon        = "#009E73",  # bluish green
  Saimiri_oerstedii       = "#F0E442",  # yellow
  Saimiri_sciureus        = "#0072B2",  # blue
  Saimiri_ustus           = "#D55E00"   # vermillion
)

#Assign species using your predefined mapping
spp <- rep(NA, length(pca$ind))  # Initialize empty vector
#Apply the species mapping
pca <- pca %>%
  mutate(
    original_id = ind,
    spp = species_mapping[ind],
    ind = paste0(species_mapping[ind], "_", original_id)
  )
#Convert to tibble
pca <- as_tibble(data.frame(pca, spp))
#Percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval / sum(eigenval) * 100)
#Scree plot
PCA_Saimiri <- ggplot(pve, aes(PC, pve)) +
  geom_bar(stat = "identity") +
  ylab("Percentage variance explained") +
  theme_light()
ggsave("Identified/PCA_Saimiri_Identified.png", plot = PCA_Saimiri, height = 1280, width = 2048, units = "px")
#Print cumulative variance
print(cumsum(pve$pve))

PCA_Saimiri_map_labelled <- ggplot(pca, aes(PC1, PC2, col = spp)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = original_id), size = 1.5, max.overlaps = 100, box.padding = 0.5, point.padding = 0.25) +
  scale_colour_manual(values = pop_colors) +
  coord_equal() + theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave("Identified/PCA_Saimiri_Identified__map_labelled.png", plot = PCA_Saimiri_map_labelled, height = 1280, width = 2048, units = "px")

PCA_Saimiri_map <- ggplot(pca, aes(PC1, PC2, col = spp)) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = pop_colors) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave("Identified/PCA_Saimiri_Identified_map.png", plot = PCA_Saimiri_map, height = 1280, width = 2048, units = "px")
######
```

## II. Principal Component Analysis (all samples identified - without S. oerstedii)
### 1. Perform linkage disequilibrium (LD) pruning

```bash
plink --vcf /srv/home/gpierre/Primates/memo_saimiri/Saimiri_filtered_final.vcf.gz \
      --keep Samples_identified_oerstedii_out.txt \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 \
      --out Identified_no_oerstedii/Saimiri_identified_no_oerstedii
```

### 2. Perform PCA using only the pruned SNPs

```bash
plink --vcf /srv/home/gpierre/Primates/memo_saimiri/Saimiri_filtered_final.vcf.gz \
      --keep Samples_identified_oerstedii_out.txt \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --extract Identified_no_oerstedii/Saimiri_identified_no_oerstedii.prune.in \
      --make-bed --pca --out Identified_no_oerstedii/Saimiri_identified_no_oerstedii
```

### 3. Plotting the PCA (all samples identified - without S. oersetdii) output with R

```bash
### R script ###
#Load necessary packages
if (!require(tidyverse)) install.packages("tidyverse", repos = "https://cloud.r-project.org")
if (!require(ggrepel)) install.packages("ggrepel", repos = "https://cloud.r-project.org")
library(tidyverse)
library(dplyr)
library(ggrepel)
#Read in PCA data
pca <- read_table("Identified_no_oerstedii/Saimiri_identified_no_oerstedii.eigenvec", col_names = FALSE)
eigenval <- scan("Identified_no_oerstedii/Saimiri_identified_no_oerstedii.eigenval")
#Preprocess PCA data
pca <- pca[, -1]  # Remove extra column
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
#Species mapping without Saimiri_oerstedii
species_mapping <- c(
  #Saimiri boliviensis
  "PD_2704" = "Saimiri_boliviensis", "PD_1233" = "Saimiri_boliviensis", "PD_1232" = "Saimiri_boliviensis", "PD_2703" = "Saimiri_boliviensis", "PD_2706" = "Saimiri_boliviensis", "PD_2705" = "Saimiri_boliviensis", "PD_0370" = "Saimiri_boliviensis",
  #Saimiri cassiquiarensis
  "PD_1251" = "Saimiri_cassiquiarensis", "PD_0359" = "Saimiri_cassiquiarensis", "PD_0358" = "Saimiri_cassiquiarensis",
  "PD_0354" = "Saimiri_cassiquiarensis", "PD_0356" = "Saimiri_cassiquiarensis", "PD_0355" = "Saimiri_cassiquiarensis",
  "PD_0360" = "Saimiri_cassiquiarensis", "PD_0353" = "Saimiri_cassiquiarensis", "PD_0373" = "Saimiri_cassiquiarensis",
  "PD_0357" = "Saimiri_cassiquiarensis", "PD_1226" = "Saimiri_cassiquiarensis",
  #Saimiri macrodon
  "PD_1238" = "Saimiri_macrodon", "PD_0362" = "Saimiri_macrodon", "PD_1243" = "Saimiri_macrodon", "PD_1230" = "Saimiri_macrodon",
  "PD_1240" = "Saimiri_macrodon", "PD_1242" = "Saimiri_macrodon", "PD_1237" = "Saimiri_macrodon", "PD_1234" = "Saimiri_macrodon",
  "PD_0361" = "Saimiri_macrodon", "PD_0372" = "Saimiri_macrodon", "PD_1239" = "Saimiri_macrodon", "PD_1246" = "Saimiri_macrodon",
  "PD_1231" = "Saimiri_macrodon", "PD_1236" = "Saimiri_macrodon", "PD_1235" = "Saimiri_macrodon", "PD_1244" = "Saimiri_macrodon",
  "PD_1245" = "Saimiri_macrodon", "PD_1247" = "Saimiri_macrodon",
  #Saimiri oerstedii
  "PD_0176" = "Saimiri_oerstedii",
  #Saimiri sciureus
  "WPB399" = "Saimiri_sciureus", "WPB374" = "Saimiri_sciureus", "PD_1250" = "Saimiri_sciureus", "PD_1249" = "Saimiri_sciureus",
  "PD_1248" = "Saimiri_sciureus", "PD_0364" = "Saimiri_sciureus", "PD_0365" = "Saimiri_sciureus",
  #Saimiri sp
  "PD_2174" = "Saimiri_sp", "PD_1310" = "Saimiri_sp", "PD_1275" = "Saimiri_sp",
  #Saimiri ustus
  "PD_1227" = "Saimiri_ustus", "PD_0371" = "Saimiri_ustus", "PD_0375" = "Saimiri_ustus", "PD_1241" = "Saimiri_ustus",
  "PD_1229" = "Saimiri_ustus", "PD_0363" = "Saimiri_ustus", "PD_0374" = "Saimiri_ustus",
  "PD_0368" = "Saimiri_ustus", "PD_0366" = "Saimiri_ustus", "PD_0367" = "Saimiri_ustus", "PD_1225" = "Saimiri_ustus",
  "PD_1228" = "Saimiri_ustus"
)
pop_colors <- c(
  Saimiri_boliviensis     = "#E69F00",  # orange
  Saimiri_cassiquiarensis = "#56B4E9",  # sky blue
  Saimiri_macrodon        = "#009E73",  # bluish green
  Saimiri_oerstedii       = "#F0E442",  # yellow
  Saimiri_sciureus        = "#0072B2",  # blue
  Saimiri_ustus           = "#D55E00"   # vermillion
)
#Assign species using your predefined mapping
spp <- rep(NA, length(pca$ind))  # Initialize empty vector
#Apply the species mapping
pca <- pca %>%
  mutate(
    original_id = ind,
    spp = species_mapping[ind],
    ind = paste0(species_mapping[ind], "_", original_id)
  )
#Convert to tibble
pca <- as_tibble(data.frame(pca, spp))
#Variance explained per PC
pve <- data.frame(PC = 1:20, pve = eigenval / sum(eigenval) * 100)
#Scree plot
PCA_Saimiri <- ggplot(pve, aes(PC, pve)) +
  geom_bar(stat = "identity") +
  ylab("Percentage variance explained") +
  theme_light()
ggsave("Identified_no_oerstedii/PCA_Saimiri_no_oerstedii_scree.png", plot = PCA_Saimiri, height = 1280, width = 2048, units = "px")
#Print cumulative variance
print(cumsum(pve$pve))
#PCA scatter plot
PCA_Saimiri_map <- ggplot(pca, aes(PC1, PC2, col = spp)) +
  geom_point(size = 3) +
  scale_colour_manual(values = pop_colors) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave("Identified_no_oerstedii/PCA_Saimiri_no_oerstedii_map.png", plot = PCA_Saimiri_map, height = 1280, width = 2048, units = "px")
PCA_Saimiri_map_label <- ggplot(pca, aes(PC1, PC2, col = spp, label = original_id)) +
  geom_point(size = 3) +
  geom_text_repel(size = 1.5, max.overlaps = 100) +  # add labels
  scale_colour_manual(values = pop_colors) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave("Identified_no_oerstedii/PCA_Saimiri_identified_no_oerstedii_map_with_labels.png", plot = PCA_Saimiri_map_label, width = 2048, height = 1280, units = "px")
###
```

## III. Principal Component Analysis 3D (all samples identified)
### 1.

```bash
### R Script ###
install.packages("scatterplot3d")  # if needed
library(scatterplot3d)

# keep only samples with species assigned
p3 <- dplyr::filter(pca, !is.na(spp))

# colors + shapes per species (tweak as you like)
pop_colors <- c(
  boliviensis     = "#E69F00",  # orange
  cassiquiarensis = "#56B4E9",  # sky blue
  macrodon        = "#009E73",  # bluish green
  oerstedii       = "#F0E442",  # yellow
  sciureus        = "#0072B2",  # blue
  ustus           = "#D55E00"   # vermillion
)
shps <- c(
  Saimiri_boliviensis = 16,   # solid circle
  Saimiri_cassiquiarensis = 17,# triangle
  Saimiri_macrodon = 15,      # square
  Saimiri_oerstedii = 18,     # diamond
  Saimiri_sciureus = 19,      # big circle
  Saimiri_ustus = 20          # hollow circle
)

pt_col <- cols[p3$spp]
pt_pch <- shps[p3$spp]

# axis labels with variance explained
xl <- paste0("PC1 (", round(pve$pve[1], 2), "%)")
yl <- paste0("PC2 (", round(pve$pve[2], 2), "%)")
zl <- paste0("PC3 (", round(pve$pve[3], 2), "%)")

# save a static PNG
png("Identified/PCA_Saimiri_Identified_3D_static.png",
    width = 2048, height = 1600, res = 300)
par(mar = c(4, 4, 2, 8), xpd = NA)  # extra right margin for legend

s3d <- scatterplot3d(
  x = p3$PC1, y = p3$PC2, z = p3$PC3,
  xlab = xl, ylab = yl, zlab = zl,
  color = pt_col, pch = pt_pch, cex.symbols = 1.0,
  grid = TRUE, box = FALSE,
  angle = 55,    # camera azimuth
  theta = 25     # elevation; remove if your version lacks 'theta'
)

# legend (on the right)
legend("right", inset = c(-0.12, 0),
       legend = names(cols),
       col = unname(cols),
       pch = unname(shps[names(cols)]),
       pt.cex = 1, bty = "n", y.intersp = 1.1)

dev.off()
###
```

# 🧬 B. Admixture Analysis
## I. Admixture (all samples identified)
*Script adapted from* [Speciation Genomics ADMIXTURE tutorial](https://speciationgenomics.github.io/ADMIXTURE/)  
*Manual reference:* [ADMIXTURE User Guide (UCLA)](http://software.genetics.ucla.edu/admixture/admixture-manual.pdf)

### 0. Initialize Environment
Admixture was installed via **conda**:
```bash
conda install bioconda::admixture
```

### 0.1 Initialize and load module
```bash
module load PLINK/2.00a3.7-gfbf-2023a
module load BCFtools/1.14-GCC-11.2.0
```
### 1. Extract Individuals of Interest
```bash
bcftools view \
  --samples-file Analysis_admixture/Identified/identified_individuals_species.txt \
  --output-type z \
  --output Analysis_admixture/Identified/Saimiri_admixture.vcf.gz \
  --threads 16 \
  Saimiri_filtered_final.vcf.gz
```
### 2. Creation of PLINK Input File
This step converts the filtered VCF into PLINK binary files required by ADMIXTURE.  
It also performs linkage pruning to remove correlated SNPs and assigns default variant IDs.

```bash
plink --vcf /srv/home/gpierre/Primates/memo_saimiri/Saimiri_filtered_final.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --make-bed --out /srv/home/gpierre/Primates/memo_saimiri/Analysis_admixture/Identified_2nd/Saimiri_identified_admx
``` 

### 2.1 Chromosome Name Correction
ADMIXTURE accepts only numeric chromosome names. This command replaces the first column of the BIM file with zero to make all chromosomes valid.

```bash
awk '{$1="0";print $0}' Saimiri_identified_admx.bim > Saimiri_identified_admx.bim.tmp

mv Saimiri_identified_admx.bim.tmp Saimiri_identified_admx.bim
```
### 3. Running ADMIXTURE
ADMIXTURE estimates ancestry proportions for different numbers of clusters.  
Each run tests a value of K, which represents the number of ancestral populations.  
The cross-validation error helps identify the best K value.

```bash
# test 1 : 
admixture -j16 --cv Saimiri_identified_admx.bed 5 | tee log5.out
# test 2 : 
nohup admixture --cv Saimiri_identified_admx.bed 5 > log5.out 2>&1 &
```
### 3.1 Running ADMIXTURE Loop
To test several K values, a loop is used to run ADMIXTURE from K equals two to K equals ten. Each run records time and saves cross-validation results in a summary file.

```bash
nohup bash -c '
set -euo pipefail
start=$(date +%s)   # record global start time
total=9             # total number of runs (K=2..10)
for i in {2..10}; do
  run_start=$(date +%s)
  echo "[$(date)] Running ADMIXTURE K=$i with 16 threads..."
  admixture -j16 --cv Saimiri_identified_admx.bed $i |& tee "log${i}.out"
  run_end=$(date +%s)
  elapsed=$((run_end - start))
  last_runtime=$((run_end - run_start))
  done_so_far=$((i - 1))
  remaining=$((total - done_so_far))
  #Estimate remaining time based on last run’s duration
  eta=$((remaining * last_runtime / 3600))  # in hours (rounded down)
  echo "[$(date)] Finished K=$i (last run ${last_runtime}s)."
  echo "    Total elapsed: $((elapsed/60)) min."
  echo "    Estimated time remaining: ~${eta}h"
done
#Extract CV error lines and make a clean table: K <tab> CV_error
grep -h "CV error" log*.out \
  | sed -E "s/.*\(K=([0-9]+)\):[[:space:]]*([0-9.]+)/\1\t\2/" \
  | sort -n > cv_summary.tsv
echo "[$(date)] Finished all runs. CV summary written to cv_summary.tsv"
' > admixture_nohup.out 2>&1 &
```

### 3.2 Extracting Cross-Validation Errors
These commands collect and clean CV error values from all log files into one summary.

```bash
awk '/CV/ {print $5,$9}' *out | cut -c 4,7-20 > saimiri_identified_admx.cv.error

grep "CV" *out | awk '{print $5,$9}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > saimiri_identified_admx.cv.error

grep "CV" *out | awk '{print $5,$9}' | cut -c 4,7-20 > saimiri_identified_admx.cv.error
```
### 4. Visualization and Plotting
This step uses an R script to visualize ADMIXTURE results and identify the best K value. The script creates bar plots showing ancestry proportions for each individual.

### 4.1 Download R Script
The following commands download the R script used to generate the ADMIXTURE plots. The script requires four input arguments: the prefix of the ADMIXTURE output files,  the file containing species information,  the maximum K value to plot,  and the list of populations in the order they should appear.

```bash
wget https://github.com/speciationgenomics/scripts/raw/master/plotADMIXTURE.r

chmod +x plotADMIXTURE.r
```
### 4.2 Run the R Script
These commands run the plotting script with different maximum K values. The populations are listed in the desired display order.

```bash
  Rscript plotADMIXTURE.r \
  -p Saimiri_identified_admx \
  -i identified_individuals_species.txt \
  -k 5 \
  -l Saimiri_boliviensis,Saimiri_cassiquiarensis,Saimiri_macrodon,Saimiri_sciureus,Saimiri_ustus,Saimiri_oerstedii
```
### 4.3 Plotting Cross-Validation Error
This R code plots the cross-validation error for each tested K value and highlights the lowest error.

```bash
### R SCRIPT ###
library(ggplot2)
library(readr)

cv <- read_table("saimiri.cv.error", col_names = c("K", "CVerror"))

bestK <- cv$K[which.min(cv$CVerror)]

p <- ggplot(cv, aes(x = K, y = CVerror)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  geom_point(data = cv[cv$K == bestK, ], aes(x = K, y = CVerror),
             color = "red", size = 4) +
  geom_text(data = cv[cv$K == bestK, ], aes(label = paste0("Best K = ", K)),
            vjust = -1, color = "red") +
  theme_bw(base_size = 14) +
  labs(title = "ADMIXTURE cross-validation error",
       x = "Number of clusters (K)",
       y = "CV error")

ggsave("CV_error_plot.png", plot = p, width = 6, height = 4, dpi = 300)
###
```

### 5. Excluding *Saimiri oerstedii* for a Finer Analysis
This step removes *Saimiri oerstedii* individuals from the dataset to focus on variation within the remaining species.  The filtered file is then reindexed and prepared for a new ADMIXTURE run.

### 5.1 Remove the Sample and Reindex
The following commands remove the individual PD_0176, which belongs to *S. oerstedii*, from the VCF file.  A new file is created without this sample, followed by indexing and verification.

```bash
INVCF=/srv/home/gpierre/Primates/memo_saimiri/Saimiri_filtering_process/Saimiri_MAF005.vcf.gz

bcftools query -l "$INVCF" | grep -n '^PD_0176$'

bcftools view -s ^PD_0176 -Oz -o /srv/home/gpierre/Primates/memo_saimiri/Saimiri_filtering_process/Saimiri_MAF005_no_oerstedii.vcf.gz "$INVCF"

bcftools index -t /srv/home/gpierre/Primates/memo_saimiri/Saimiri_filtering_process/Saimiri_MAF005_no_oerstedii.vcf.gz

bcftools query -l /srv/home/gpierre/Primates/memo_saimiri/Saimiri_filtering_process/Saimiri_MAF005_no_oerstedii.vcf.gz | grep PD_0176 || echo "✅ PD_0176 successfully removed"
```
### 5.1 Rerun ADMIXTURE Without S. oerstedii
After creating the new VCF file, repeat the ADMIXTURE workflow using this updated dataset. Follow the same steps as before (from 3. to 4.3), starting from PLINK file conversion and running ADMIXTURE for the desired K values.

