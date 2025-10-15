# A Hard filtering process - Introduction and precision   
#First, let's keep in mind what the lab where samples were sequenced said:
#1. There may be related individuals and admixed individuals in the dataset. 
#2. Individuals with problems in our QC (coverage, missing data and allele balance) were not removed, find the information in the list of samples and remove them if desired. 
#There are no Pithecia individuals in that list
#3. We have run our internal QC, but there may be mistakes in the dataset. Please before starting all the analyses run your QC and let us know if you find anything wrong. 
#4. The VCF files are not filtered. Apply your own filtering criteria depending on the analyses you want to do
#From Nuria's paper (https://doi.org/10.1038/s42003-024-06901-3): 
#Hard filtering was applied to the resulting SNPs in each chunk by removing those variants that were not biallelic (vcflib vcfbilallelic (v10.5)),
#those with depth values outside the range between 70 (double the median maximum coverage across samples) and 5 (one third of the median minimum
#coverage across samples), and missingness higher than 60% using vcftools (v0.1.12). 
#Furthermore, following GATK Best Practices Protocol, we excluded those SNPs not complying with the following expression:
#“QD< 2 | FS >60|MQ<40 | SOR> 3 |ReadPosRankSum<–8.0 | MQRankSum < –12.5” employing GATK VariantFiltration. 
#We merged the remaining SNPs in each chunk-VCF into a single dataset, which we then filtered for allele imbalance by keeping variants with frequencies within the
#range of 0.25–0.75 (Supplementary Methods (Section 1.2), Fig. S1). Lastly, we removed SNPs in scaffolds belonging to sex chromosomes (identified in Kuderna et al) and those in scaffolds shorter than 0.5Mb from the dataset, thus keeping 367 scaffolds.



# B Basic command lines
## 1 module loading
module load BCFtools/1.14-GCC-11.2.0
module load VCFtools/0.1.16-GCC-12.3.0

## 2 Number of variants
bcftools view -H {files-names} | wc -l 

## 3 Creating sample files
#Creating a sample_file with the codes. It may be useful later:
bcftools query -l {input - file - name}   | sort | uniq > Analysis_BCF/{output - file - name}.txt

## 4 Create statistics to compare files
#Compare Files to See Differences. Use bcftools stats to compare before and after filtering
bcftools stats {input - file - name} > {output - file - name}.txt

# C BCFtools Analysis
## 1 Number of variants
bcftools view -H Saimiri_boliviensis__Saimiri.pop.geno.snps.raw.reheader.vcf.gz| wc -l 
-> 80 590 885 variants
-> 59 individuals (see tables "Saimiri_samples" to have acces to id and species)

## 2 general statistics from vcf
#using bcf tools to extract some general statistics from vcf
bcftools stats Saimiri_boliviensis__Saimiri.pop.geno.snps.raw.reheader.vcf.gz  > Analysis_BCF/Saimiri_BCF_stats.txt

# D Removing non-biallelic SNPs 
## 1 remove non-biallelic SNPs
#Using bcftools to remove non-biallelic SNPs and mantain only biallelic SNPs. This is a necessary step for downstream analyses (PCA, Structure, etc) that will assume you are using biallelic SNPs
nohup bash -c 'bcftools view -m2 -M2 -v snps -o Saimiri_boliviensis__Saimiri_filtered_biallelic.vcf.gz Saimiri_boliviensis__Saimiri.pop.geno.snps.raw.reheader.vcf.gz ' > biallelic_log.out 2>&1 &

## 2 check if only biallelic remained
#Checking Non-Biallelic SNPs Removed
bcftools view -m2 -M2 -v snps Saimiri_boliviensis__Saimiri.pop.geno.snps.raw.reheader.vcf.gz | wc -l 
-> 74 700 345 biallelic SNPs
Saimiri_boliviensis__Saimiri_biallelic.vcf.gz -> 74 698 589 variants

## 3 Indexing the output file 
bcftools index Saimiri_boliviensis__Saimiri_filtered_biallelic.vcf.gz

# E GATK Hard Filtering
## 1.0 Applying quality-based filtering
#Double-check the name of files
#First, you have to create an index file in the format that GATK will recognize

## 1.1 Loading tabix
module load tabix/0.2.6-GCCcore-10.2.0

## 1.2 Running tabix
tabix -p vcf Saimiri_boliviensis__Saimiri_filtered_biallelic.vcf.gz

## 2 Loading GATK 
module load GATK/4.2.4.1-GCCcore-10.2.0-Java-1.8

## 3 Running GATK
nohup gatk VariantFiltration \
    -V Saimiri_boliviensis__Saimiri_filtered_biallelic.vcf.gz \
    -O Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || ReadPosRankSum < -8.0 || MQRankSum < -12.5" \
    --filter-name "HardFilter" > HF_log_out 2>&1 &
#This step tags variants that fail quality thresholds.

## 4 Indexing file with BCF
bcftools index Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz

## 5 Count variants 
bcftools view -H Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz| wc -l
#Before Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz ->  
#After Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz ->  74 698 589 variants

## 6 Check the filter colunm 
#This will show counts for PASS and filtered variants (e.g., HardFilter)
bcftools view -H Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz | cut -f 7 | sort | uniq -c
-> 5056154 HardFilter
-> 69642435 PASS

## 7.0  Remove failed variants:
bcftools view -f PASS -Oz -o Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz Saimiri_boliviensis__Saimiri_filtered_GATK.vcf.gz

## 7.1 Indexing
module load tabix/0.2.6-GCCcore-10.2.0 
tabix -p vcf Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz
module load BCFtools/1.14-GCC-11.2.0
bcftools index Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz

## 8 Count variants
#After filtering Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz -> 69 642 435 -> 10 948 450 variants removed from the original file

## 9 Check missingness
#This will generate missingness.imiss, where: N_DATA is the total number of genotypes per sample. F_MISS is the fraction of missing genotypes per sample
vcftools --gzvcf Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz --missing-indv --out missingness

# F Checking missing data per sample - not done
## 1.0 removing missingness 
#Using vcftools to check missing data
vcftools --gzvcf Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz --missing-indv --out missingness
#This will generate missingness.imiss, where: N_DATA is the total number of genotypes per sample. F_MISS is the fraction of missing genotypes per sample
#if some problems, try this one 
nohup bash -c "bcftools filter -i 'F_MISSING <= 0.6' -Oz -o Saimiri_boliviensis__Saimiri_filtered_final_missingindv.vcf.gz Saimiri_boliviensis__Saimiri.pop.geno.snps.raw.reheader.vcf.gz" > bcf_missing.out 2>&1 &
#Remove SNPs with >60% Missing Data
#nohup bash -c 'vcftools --vcf Pithecia_depth_filtered.vcf.gz --max-missing 0.4 --recode --out Pithecia_missing_filtered' > missing_log.out 2>&1 &
nohup bash -c "bcftools filter -i 'F_MISSING <= 0.6' -Oz -o Pithecia_missing_filtered.vcf.gz Pithecia_depth_filtered.vcf.gz" > bcf_missing.out 2>&1 &
## 1.1 Count variants
#After filtering Pithecia_final_filtered.vcf.gz -> 61 993 478 -> XXXX variants removed


# G VCFtools Analysis
-> https://speciationgenomics.github.io/filtering_vcfs/ 
## 1 Create VCF variables  
create VCF variables - ONLY BIALLELIC VCF FILES : 
INPUT=~/Primates/memo_saimiri/Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz
OUTPUT=~/Primates/memo_saimiri/Analysis_VCF/Saimiri_stats_VCF
#Calculate allele Frequency
vcftools --gzvcf $INPUT --freq2 --out $OUTPUT --max-alleles 2
#Calculate mean depth per individual
vcftools --gzvcf $INPUT --depth --out $OUTPUT
#Calculate mean depth per site
vcftools --gzvcf $INPUT --site-mean-depth --out $OUTPUT
#Calculate site quality
vcftools --gzvcf $INPUT --site-quality --out $OUTPUT
#Calculate proportion of missing data per individual
vcftools --gzvcf $INPUT --missing-indv --out $OUTPUT
#Calculate proportion of missing data per site
vcftools --gzvcf $INPUT --missing-site --out $OUTPUT
Calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf $INPUT --het --out $OUTPUT

## 2.0 Calculating depth 
#Extract Per-Sample Depths
#Calculating the depth values per samples. These are the values to calculate the one-third of the median minimum coverage and the double the median maximum coverage
bcftools query -f '[%DP\n]' Saimiri_boliviensis__Saimiri_filtered_final.vcf.gz| grep -v '^0$' > Analysis_BCF/all_depths.txt
#In the output file you will find a list of millions of numbers in one single colunm. These values are the depths (DP) for each sample at a specific SNP and will be used to calculate the mode         
#Check if there are only numeric values in the output file
awk '{for (i=1; i<=NF; i++) if ($i ~ /^[0-9]+$/) printf "%s ", $i; print ""}' all_depths.txt > Analysis_BCF/all_depths_cleaned.txt

## 2.1 Calculate the Per-Sample Median Depth
#Now, computing the median depth per sample
#Findging the mode:
#The mode reflects typical, consistent coverage
#The mode is the most frequent read depth across the genome.
#it represents the typical coverage per base per sample in your dataset, excluding outliers (unlike the mean, which can be skewed).
#So it's a stable reference point to set thresholds around.
sort all_depths_cleaned.txt | uniq -c | sort -nr | head -1
-> I got the following output 326856017 35, which means that the mode is 35

mode=35  # replace with your actual value
min_thresh=$(echo "$mode / 3" | bc) #bc is a command line calculator "Basic Calculator"
max_thresh=$(echo "$mode * 2" | bc)
echo "Min DP: $min_thresh"
echo "Max DP: $max_thresh"
#I got:
#Min DP = 11
#Max DP = 70

## 3.0 Set filters for VCF hard filtration
MISS=0.75
QUAL=30
MIN_DEPTH=11
MAX_DEPTH=70

## 3.1 Perform the filtering with vcftools
vcftools --gzvcf Saimiri_boliviensis__Saimiri_biallelic.vcf.gz \
--remove-indels --max-missing $MISS \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > Saimiri_filtered.vcf.gz

## 3.2 Checking variants remained
Saimiri_filtered.vcf.gz  -> 67 420 702 variants

## 3.3 Visualise depth 
#To visualise it with this python script: 
nano plot_depths_distribution.py
#######################################PYTHON SCRIPT#####################################
import numpy as np
import matplotlib.pyplot as plt
#Load depths from file
depths = []
with open("all_depths_cleaned.txt") as f:
    for line in f:
        try:
            dp = int(line.strip())
            depths.append(dp)
        except:
            continue

depths = np.array(depths)
#Calculate statistics
mode = int(np.bincount(depths).argmax())
p95 = int(np.percentile(depths, 95))
print(f"Mode depth: {mode}")
print(f"95th percentile depth: {p95}")
#plot histogram
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
#########################################################################################
#Run this file in a conda environment
#in my case, I have created the environment called wgs. Therefore, I used:
conda activate wgs
#Then, I installed numpy and matplotlib.
pip install numpy matplotlib
#...and I run the python script below to get a .png figure with the depth distributions
python plot_depth_distribution.py
Mode depth: 35
95th percentile depth: 44
#In my case, the 95% depth valus as a criterion for the maximum DP cutoff works well

#######################################################################################

# Z useful COMMAND LINES for terminal 

## 1 
#if you want the full command + time of what is running in the terminal:
top -u gpierre

## 2
#if you want to have acces to prompt again
ctrl + c

## 3
#pour visualiser un ficher dans un éditeur de txt dans le terminal
nano
ctrl + x pour sortir 

## 4
#Pour lancer R : r
#Pour quitter R : q()

## 5
#Pour lancer python : python
#Pour quitter python :  ctrl + z

## 6
#Pour visualiser un fichier ou le télécharger, utilise GIT bash sur l'ordinateur. ça ouvre un terminal capable d'aller rechercher un fichier sur le server ebe1.
gilp4@HP-GilP MINGW64 ~
$ scp gpierre@10.75.1.146:~/my_plot.png /c/Users/gilp4/Documents/Gil_ULB
The authenticity of host '10.75.1.146 (10.75.1.146)' can't be established.
ED25519 key fingerprint is SHA256:y6rxYfNKKiYt/aJK+JkMBjAHBG3Ki1N5TfPbAaOiZfw.
This host key is known by the following other names/addresses:
    ~/.ssh/known_hosts:1: ebe-gpu01.hpda.ulb.ac.be
Are you sure you want to continue connecting (yes/no/[fingerprint])? y
Please type 'yes', 'no' or the fingerprint: yes
Warning: Permanently added '10.75.1.146' (ED25519) to the list of known hosts.
(gpierre@10.75.1.146) Password:

## 7 
#meme chose que 6 mais pour faire tout un dossier cette fois il suffit de rajouter -r apres scp

## 8
#to unzip a file .gz with the function "-k" which keep the input file in a compressed way
gunzip -k file.gz

## 9
#When you want to transfer a file from yout local pc to the server -> use git bash 
gilp4@HP-GilP MINGW64 ~
$ scp "/c/Users/gilp4/OneDrive - Université Libre de Bruxelles/Memo_Saimiri/Saimiri_data/Saimiri_R/Saimiri_biallelic_stats_plots_code_server.R" gpierre@10.75.1.146:/srv/home/gpierre/Primates/Saimiri/Analysis_VCF/

## 10
#When you want to move several files. if they have the same extension you can use this formula : 
mv *.ext  *.xml *.txt /path/to/dest/folder/

## 11
#view header
bcftools view -h Saimiri_filtered.vcf.gz 

## 12
#view 10 first lines without header
bcftools view -H Saimiri_filtered.vcf.gz | head -n 5

## 13
#get samples names
bcftools query -l Saimiri_filtered.vcf.gz

## 14
#remove samples and variants 
bcftools view -S ^PD_2174 , PD_1310 , PD_1275 -Oz -o filtered_samples.vcf.gz 
Saimiri_filtered_samples.vcf.gz


echo -e "PD_2174\nPD_1310\nPD_1275" > exclude_samples.txt
bcftools view -S ^exclude_samples.txt -Oz -o filtered_samples.vcf.gz filtered_samples.vcf.gz