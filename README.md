Local ancestry PheWAS pipeline scripts and summary statistics
=================
NOTE: This file is a work in progress expected completion: 10/23/24

Citation: "Systematic comparison of phenome wide admixture mapping and genome-wide association in a diverse health system biobank"

Author: Sinead Cullina

##Pipeline Summary##
Calculate global ancestry, phase data and subset, run local ancestry with GNOMIX, convert GNOMIX output to VCF files, compare local ancestry output to global ancestry proportions, filter samples and variants, calculate GRM, run SAIGE step 1 and step 2 for admixture mapping, admixture mapping results post processing, pre GWAS QC of genotype data, run SAIGE step 1 and step 2 for GWAS.

This pipeline is composed of a series of scripts. This document describes the order, input and output of each of these scripts to conduct an admixture pheWAS.

v
## Pipeline Map: ##
#### 0.) Global ancestry inference  #####

* Overview
  * Merge files using plink
  * Filter for ADMIXTURE input
  * Run ADMIXTURE
  * Remove samples with complex admixture patterns.

##### 1.) Infer local ancestry #####
  * Downsample merged file to keep reference panels and samples passing QC
  * Phase using EAGLE
  * Convert to VCF format using SHAPEIT 
  * Normalize
  * Infer local ancestry using GNOMIX
  * Plot and filter GNOMIX output

##### 2.1) Run admixture mapping: Two-way and three-way #####
* Convert GNOMIX output to VCFs for SAIGE
* Create genotype files for SAIGE Step 1
* Run SAIGE Step 1
* Run SAIGE Step 2
* Process Results

##### 2.2) Run GWAS #####
* Create genotype files for SAIGE Step 1
* Run SAIGE Step 1
* Run SAIGE Step 2

##### 3.1) Plotting results #####


## 0.) Global ancestry inference ###
##### Overview #####
You can use PLINK v1.9 to merge your query genotype data and reference panels. First you should find the maximum number of overlapping variants between all files. Then you may need to flip and rename some of the variants so that the order and name of the variants matches between all genotype files. I usually use PLINKv2 for these steps. As of writing PLINKv2 does not have a suitable merging function.  

Inputs:
Query dataset (i.e. samples you are using for association testing)
Reference panels 

To flip variants:  

```
plink2  --alt1-allele gda_to_flip2.txt \
  --bfile GSA_GDA_V2_TOPMED_HG38_id_update \
  --make-bed \
  --out GSA_GDA_V2_TOPMED_HG38_id_update_flipped \
```
To rename variants: 

```
plink2  --bfile GSA_GDA_V2_TOPMED_HG38
  --extract range GDA_genotyped_sites.plink
  --keep-allele-order
  --make-bed
  --max-alleles 2
  --new-id-max-allele-len 66
  --out GSA_GDA_V2_TOPMED_HG38_id_update
  --set-all-var-ids @:#:$r:$a
 ```
 
To merge genotype files:

```
plink --bfile GSA_GDA_V2_TOPMED_HG38_id_update_flipped_ID
  --bmerge TGP_HGDP_WHI_SGDP_BIOME
  --extract gda_plus_refs_perfect_match.plink
  --make-bed
  --out GDA_GSA_imputed_TGP_HGDP_WHI_SGDP_BIOME

```

I then use this merged dataset for global ancestry inference. Firstly I apply some filtering steps:

#remove rare variants and missingnesss filters. Also remove regions of the genome known to be under strong selection (see exclusion_regions_under_selection.bed file) 
```
 plink2  --bfile ../GDA_GSA_imputed_TGP_HGDP_WHI_SGDP_BIOME
  --exclude range exclusion_regions_under_selection.bed
  --geno 0.05
  --mac 10
  --maf 0.01
  --make-bed
  --mind 0.05
  --out 
 
 ```
Remove palindromic sites, LD prune and remove 2nd degree relateds: 


Then launch ADMIXTURE: 
Run according to the documentation (https://dalexander.github.io/admixture/)
citation: Alexander, David H., John Novembre, and Kenneth Lange. "Fast model-based estimation of ancestry in unrelated individuals." Genome research 19.9 (2009): 1655-1664.

```
geno_path=/sc/arion/projects/kennylab/Sinead/LAI_GDA_GSA/HIS_AA_imputed_merge/ADMIXTURE/

for i in {2,3,4,5,6}

do
file=GDA_GSA_imputed_TGP_HGDP_WHI_SGDP_BIOME_mind_maf_geno_no_palindrome_exclusion_regions_no_duplicates_LDprune_cell_paper.bed
echo 'cd '${geno_path}'' > Admixture_cell_${i}.pbs
echo 'module load admixture' >> Admixture_cell_${i}.pbs
echo 'admixture '${geno_path}$file' '${i}' -j48 --cv' >> Admixture_cell_${i}.pbs

## if premium doesn't work for you, you may want to try changing to -q alloc

bsub -q premium -P acc_ipm2  -n 48 -W 20:00 -R rusage[mem=6000] -o Admixture_cell_${i}.log < Admixture_cell_${i}.pbs

done

```

The output will consist of files that end in .Q and .P for each admixture component tested.


Script to order admixture components, visualize, select reference panel samples and filter query samples:


#### 1)  Infer local ancestry ####

Use Eagle to phase genotype data. Run according to the documentation. Examples as follows:

Note: 1000 Genomes reference panels and genetic maps can be downloaded here: https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html.

```
geno_path=/sc/arion/projects/kennylab/Sinead/LAI_GDA_GSA/HIS_AA_imputed_merge/PHASING/
for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}

do
echo 'cd '${geno_path}'' > Launch_Eagle_${i}.pbs
echo 'module load eagle/2.4' >> Launch_Eagle_${i}.pbs
echo 'eagle --bfile HIS_AFR_plus_Ref_no_2nd_deg_rel --geneticMapFile /hpc/packages/minerva-common/eagle/2.4/Eagle_v2.4/tables/genetic_map_hg38_withX.txt.gz --chrom '${i}' --outPrefix Phased_HIS_AFR_plus_Ref_no_2nd_deg_rel'${i}'' >> Launch_Eagle_${i}.pbs
## if premium doesn't work for you, you may want to try changing to -q alloc

bsub -q premium -P acc_kennylab -n 10 -R span[ptile=10] -R rusage[mem=1200] -W 60:00 -o Launch_Eagle_${i}.log < Launch_Eagle_${i}.pbs

done

```

Use Eagle to phase genotype data. Run according to the documentation. Examples as follows:

```
geno_path=/sc/arion/projects/kennylab/Sinead/LAI_GDA_GSA/HIS_AA_imputed_merge/PHASING/

for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}

do
echo 'cd '${geno_path}'' > Launch_shapeit_${i}.pbs
echo 'module load shapeit/v2r900' >> Launch_shapeit_${i}.pbs
echo 'gunzip Phased_HIS_AFR_plus_Ref_no_2nd_deg_rel'${i}'.haps.gz' >> Launch_shapeit_${i}.pbs
echo 'shapeit -convert --input-haps Phased_HIS_AFR_plus_Ref_no_2nd_deg_rel'${i}' --output-vcf Phased_HIS_AFR_plus_Ref_no_2nd_deg_rel'${i}'.vcf' >> Launch_shapeit_${i}.pbs


bsub -q premium -P acc_kennylab -W 12:00 -M 70000 -o Launch_shapeit_${i}.log < Launch_shapeit_${i}.pbs

done

```

Bcftools index and normalize variants calls:

```
geno_path=/sc/arion/projects/kennylab/Sinead/LAI_GDA_GSA/HIS_AA_imputed_merge/PHASING/

for i in {1,2,3,4,5,6,7,8,9,11,10,12,13,14,15,16,17,18,19,20,21,22}

do
echo 'cd '${geno_path}'' > Launch_normalization_chr${i}.pbs
echo 'module load bcftools' >>  Launch_normalization_chr${i}.pbs
echo 'bcftools norm --check-ref -s -f /sc/arion/projects/ipm/data/GenomeReferences/GRCh38Ref/Converted/B38.primary_assembly.genome.fa Phased_HIS_AFR_plus_Ref_no_2nd_deg_rel'${i}'.vcf -Oz >  Phased_norm_HIS_AFR_plus_Ref_no_2nd_deg_rel'${i}'.vcf.gz' >>  Launch_normalization_chr${i}.pbs

bsub -q premium -P acc_kennylab -W 07:00 -M 50000 -o  Launch_normalization_chr${i}.log <  Launch_normalization_chr${i}.pbs

done

```

Separate out samples into query and reference panels for GNOMIX, in my case I had four files because I was training my own local ancestry model:
Hispanic Latino Reference VCF
Hispanic Latino Query VCF
African American Reference VCF
African American Query VCF


Run GNOMIX according to the documentation (https://github.com/AI-sandbox/gnomix). Examples as follows:

```
IN_PRE=$1 #/sc/arion/projects/kennylab/Sinead/LAI_GDA_GSA/HIS_AA_imputed_merge/GNOMIX_LAI/
OUT_PRE=$2 #launch_GNOMIX_

NOW=$(date +"%d%m%y")

CONFIG_PRE=${OUT_PRE}_
# mkdir $CONFIG_PRE
for i in {1,2,3,4,5,6,8,9,10,11,12,14,15,16,17,18,20,21,22}
do
CONFIG_ADDR=${CONFIG_PRE}${i}_subfile_
cat <<EOF > ${CONFIG_ADDR}
#!/bin/bash

#BSUB -J GNOMIX_AA_SHORT_$i
#BSUB -P acc_kennylab
#BSUB -q premium
#BSUB -n 3
#BSUB -R "span[hosts=1]"
#BSUB -R himem
#BSUB -R rusage[mem=120000]
#BSUB -W 15:00
#BSUB -o ${CONFIG_PRE}${NOW}_${i}_AA_SHORT.out
#BSUB -e ${CONFIG_PRE}${NOW}_${i}_AA_SHORT.err
#BSUB -L /bin/bash
ml anaconda3
ml bcftools
source activate /sc/arion/projects/kennylab/roohy/conda/envs/igh_gnomix/

python3 /sc/arion/projects/kennylab/roohy/utils/gnomix/gnomix.py ${IN_PRE}query_AA_Shapeit_normalized_Phased_GDA_Chromosome_${i}.vcf ${IN_PRE}gnomix_chr${i}_local_ancestry_short_ref_AA ${i} False ${IN_PRE}genetic_map_files/chr${i}.gmap ${IN_PRE}REF_SHORT_AA_Shapeit_normalized_Phased_GDA_Chromosome_${i}.vcf ${IN_PRE}AA_small_ref_panel.smap /sc/arion/projects/kennylab/roohy/utils/gnomix/configs/config_array.yaml

EOF
bsub < ${CONFIG_ADDR}
done

```

This will generate a directory for each chromosome with local ancestry calls summarized in .msp, .fb and .lai files. I then use R script  make_VCF_file_from_GNOMIX_AA.R to convert the local ancestry calls to VCF style format for two way local ancestry,   make_VCF_file_from_GNOMIX_HL.R does the same for three-way local ancestry calls. The script takes a .msp file as input, in the case of three way local ancestry it will output three different VCF files, one for each local ancestry background (i.e. AFR, EUR, NAT). For two-way local ancestry only one VCF file is output.

Merge vcf files of local ancestry calls and remove outlier samples and regions identified in QC for local ancestry inference.





