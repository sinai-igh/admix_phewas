Local ancestry PheWAS pipeline scripts and summary statistics
=================
NOTE: This file is a work in progress expected completion: 11/05/24

Citation: "Systematic comparison of phenome wide admixture mapping and genome-wide association in a diverse health system biobank"

Author: Sinead Cullina

## Pipeline Summary ##

This pipeline calculates global ancestry, phases genotype data and subsets samples, runs local ancestry with GNOMIX, converts GNOMIX output to VCF files, compares local ancestry output to global ancestry proportions, filters samples and variants, and runs SAIGE for admixture mapping and fine-mapping. It also includes some post-processing steps and result plotting.

The pipeline consists of a series of scripts and commands. Below is the order of operations, commands and scripts used to conduct an admixture PheWAS.

Toy dataset that you can use with these commands: https://drive.google.com/drive/folders/1J0xpkTBwyMI3MMCDRi24TPrpbuPQuQtn?usp=sharing

This toy dataset uses thousand genomes (1KGP) and human genome diversity project (HGDP) whole genome sequencing data published here: https://pubmed.ncbi.nlm.nih.gov/38749656/. The dataset has been downsampled to include samples with the following population labels at a selection of GDA panel sites:

##### Potential Reference Samples #####
* European: CEU, TSI
* African: YRI, LWK
* Native American: Pima, Surui, Colombian,  Karitiana, Maya (all from HGDP), MXL, PUR


##### Admixed populations for LAI #####
* CLM - Colombians from Medellín, Colombia (three-way admixture)
* ASW - African Ancestry in Southwest US (two-way admixture)

A mapper file with population labels for these samples called: sample_map_file_1kgp_hgdp.txt is also in the google drive.

## Pipeline Map: ##
#### 0.) Global ancestry inference  #####

* Overview
  * Merge files using plink
  * Filter for ADMIXTURE input
  * Run ADMIXTURE
  * Remove samples with complex admixture patterns

##### 1.) Infer local ancestry #####
  * Downsample merged files for QC-passed reference panels and samples
  * Phase using EAGLE
  * Convert to VCF format using SHAPEIT 
  * Normalize
  * Infer local ancestry using GNOMIX
  * Plot and filter GNOMIX output

##### 2.1) Run admixture mapping: Two-way and three-way #####
* Convert GNOMIX output to VCFs for SAIGE
* Create genotype files for SAIGE Step 1
* Run SAIGE Step 1 and Step 2
* Process Results

##### 2.2) Run GWAS #####
* Create genotype files for SAIGE Step 1
* Run SAIGE Step 1 and Step 2

##### 3.1) Plotting results #####
* Use plotting scripts to visualize results from both local ancestry and GWAS.

## 0.) Global ancestry inference ###
##### Overview #####

The first step in the analysis is to merge your query genotype data with the reference panels you are using.  You should find the maximum number of overlapping variants between all files and downsample your files to just these variants. You may need to flip and rename some of the variants so that the order and name of the variants is matching between all genotype files. I usually use PLINKv2 for these steps. As of writing PLINKv2 does not have a suitable merging function. You can use PLINK v1.9 to merge the filtered files at the end. If using the toy dataset this step is unnessesary as the query and reference samples are already merged.

Inputs:
*Query dataset (i.e. samples you are using for association testing)
*Reference panels 
*exclusion_regions_under_selection.bed  #provided in repository

To flip variants:  

```
plink2  --alt1-allele list_of_variants_to_flip.txt \
  --bfile query_or_reference_panel_file \
  --extract range position_intervals_to_keep.bed
  --make-bed \
  --out query_or_reference_panel_file_flipped_alleles \
```
To rename variants: 
```
plink2  --bfile query_or_reference_panel_file_flipped_alleles
  --keep-allele-order
  --make-bed
  --max-alleles 2
  --new-id-max-allele-len 66
  --set-all-var-ids @:#:\$r:\$a
  --out query_or_reference_panel_file_flipped_alleles_renamed
 ```
To merge genotype files:

```
plink --bfile query_or_reference_panel_file_flipped_alleles_renamed
  --bmerge list_of_other_genotype_files_to_merge
  --extract overlapping_snps.plink
  --make-bed
  --out merged_query_and_refs_dataset

```
I then use this merged dataset for global ancestry inference. Firstly I apply some filtering steps:

Remove rare variants and apply missingnesss filters. Also remove regions of the genome known to be under strong selection (see exclusion_regions_under_selection.bed file) 

```
 plink2  --bfile merged_query_and_refs_dataset
  --exclude range exclusion_regions.bed
  --geno 0.05
  --mac 10
  --maf 0.01
  --make-bed
  --mind 0.05
  --out merged_query_and_refs_dataset_admixture_qc
 
 ```
Remove palindromic sites, LD prune and remove 2nd degree relateds. Note that you can apply more stringent/relaxed ld filters but you should have at least ~120K variants in the end for input into ADMIXTURE for accurate admixture estimation.

```
awk '{ if (($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "C" && $6 == "G") || ($5 == "G" && $6 == "C")) print $2 }' merged_query_and_refs_dataset_admixture_qc.bim > palindromic_snps.txt #get list of palindromix sites
plink2 --bfile merged_query_and_refs_dataset_admixture_qc --indep-pairwise 50 5 0.2 --out ld_prunelist #ld prune
plink2 --bfile merged_query_and_refs_dataset_admixture_qc --king-cutoff 0.354 --out king_duplicates #remove duplicate samples (if neccesary)

cat palindromic_snps.txt ld_prunelist.out > palindromic_and_ld_snps.out
plink2 --bfile merged_query_and_refs_dataset_admixture_qc --exclude palindromic_and_ld_snps.out --keep king_duplicates.king.cutoff.in.id --make-bed --out admixture_input_file
```

Launch ADMIXTURE: 
Run according to the documentation (https://dalexander.github.io/admixture/)
citation: Alexander, David H., John Novembre, and Kenneth Lange. "Fast model-based estimation of ancestry in unrelated individuals." Genome research 19.9 (2009): 1655-1664.

```
geno_path=/path/to/your/files/

for i in {2,3,4,5,6} #set number of Ks you want to calculate

do
file=admixture_input_file.bed
echo 'cd '${geno_path}'' > admixture_${i}.pbs
echo 'module load admixture' >> admixture_${i}.pbs
echo 'admixture '${geno_path}$file' '${i}' -j48 --cv' >> admixture_${i}.pbs

bsub -q submission_parameter -P submission_parameter  -n 48 -W 20:00 -R rusage[mem=6000] -o admixture_${i}.log < admixture_${i}.pbs #if submitting to bsub compute cluster

done
```
The output will consist of files that end in .Q (admixture proportions) and .P (admixture probabilities) for each admixture component tested. 

Next we run a script to order admixture components, visualize, select reference panel samples and filter query samples
input: 
* admixture_input_file.fam
* admixture_input_file.Q 
* mapper_file.txt - file with three columns. One containing sample ID, one containing name of dataset sample came from, one containing population label

output:
* ordered_admixture_file.Q

There are two separate R scripts. One for ordering K3 ADMIXTURE output and one for ordering K4 ADMIXTURE output.
```
Rscript annotate_admixture_output.K4.R admixture_input_file.Q  admixture_input_file.fam mapper_file.txt  #for K3

Rscript annotate_admixture_output.K4.R admixture_input_file.Q  admixture_input_file.fam mapper_file.txt #for K4
```
Then plot admixture results using R script: 

```
Rscript plot_admixture.R 3 #for K3
Rscript plot_admixture.R 4 #for K4
```
Once you have matched ancestry components based on genetic similarity with reference panels you can apply sample filters. Cohort and reference panel selection based on global ancestry componenets is subjective and project-specific. Use the plotted admixture data to remove query samples that have considerable components from ancestries that are not going to be inferred in local ancestry inference. Subset reference samples to those with >90\% of the ancestry component being inferred. 

#### 1)  Infer local ancestry ####
* use merged file prior to ADMIXTURE filtering: merged_query_and_refs_dataset
* subset cohort to keep query and reference samples
* remove 2nd degree relateds
* apply missingness filter (Eagle automatically removes variants with >10\% missingness.)

Use Eagle to phase genotype data. Run according to the documentation https://alkesgroup.broadinstitute.org/Eagle/. Examples as follows:

Note: Eagle genetic maps downloaded from: https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/

```
#!/bin/bash

# Path to the genotype files
geno_path=/path/to/your/files/

# Loop through chromosomes 1 to 22
for i in {1..22}
do
  # Create a PBS script for each chromosome
  pbs_file="Launch_Eagle_${i}.pbs"

  # Write the PBS commands to the file
  echo "cd ${geno_path}" > ${pbs_file}
  echo "module load eagle/2.4" >> ${pbs_file}
  echo "eagle --bfile merged_query_and_refs_dataset_downsampled  \\
          --geneticMapFile /hpc/packages/minerva-common/eagle/2.4/Eagle_v2.4/tables/genetic_map_hg38_withX.txt.gz \\
          --chrom ${i} \\
          --outPrefix phased_query_and_refs_chr_${i}" >> ${pbs_file}

  # Submit the job to the cluster using bsub with specified resources
  bsub -q  submission_parameter  -P  submission_parameter   -n 10 -R "span[ptile=10]" -R "rusage[mem=1200]" -W 60:00 -o "Launch_Eagle_${i}.log" < ${pbs_file}

done
```
Use Shapeit to convert .haps and .sample files to .vcf. Bcftools index and normalize variants calls. split into reference, and query files for GNOMIX. Run according to the documentation. 
* https://samtools.github.io/bcftools/bcftools.html
* https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html

Example script:

```
geno_path=/path/to/your/files/

for i in {1..22}

  # Write the PBS commands to the file
do
echo 'cd '${geno_path}'' > Launch_shapeit_${i}.pbs
echo 'module load shapeit/v2r900' >> Launch_shapeit_${i}.pbs
echo 'module load bcftools' >>   Launch_shapeit_${i}.pbs
echo 'gunzip phased_query_and_refs_chr_'${i}'.haps.gz' >> Launch_shapeit_${i}.pbs
echo 'shapeit -convert --input-haps phased_query_and_refs_chr_'${i}' --output-vcf phased_query_and_refs_chr_'${i}'.vcf' >> Launch_shapeit_${i}.pbs
echo 'bcftools norm --check-ref -s -f B38.primary_assembly.genome.fa phased_query_and_refs_chr_'${i}'.vcf -Oz >  phased_norm_query_and_refs_chr_'${i}'.vcf.gz' >>  Launch_shapeit_${i}.pbs

# split into GNOMIX files:
#African American query
echo 'bcftools view -S aa_query.bcftools.keep --force-samples phased_norm_query_and_refs_chr_'${i}'.vcf.gz -Ov -o aa_query_phased_chr_'${i}'.vcf' >> Launch_shapeit_${i}.pbs

#Hispanic Latino query
echo 'bcftools view -S HL_query.bcftools.keep --force-samples phased_norm_query_and_refs_chr_'${i}'.vcf.gz -Ov -o hl_query_phased_chr_'${i}'.vcf' >>  Launch_shapeit_${i}.pbs

#Hispanic Latino Reference
echo 'bcftools view -S HL_ref.bcftools.keep --force-samples phased_norm_query_and_refs_chr_'${i}'.vcf.gz -Ov -o hl_ref_phased_chr_'${i}'.vcf' >> Launch_shapeit_${i}.pbs

#African American Reference
echo 'bcftools view -S AA_ref.bcftools.keep --force-samples phased_norm_query_and_refs_chr_'${i}'.vcf.gz -Ov -o aa_ref_phased_chr_'${i}'.vcf' >>  Launch_shapeit_${i}.pbs


# job submission parameters
bsub -q  submission_parameter  -P submission_parameter   -W 12:00 -M 70000 -o Launch_shapeit_${i}.log < Launch_shapeit_${i}.pbs

done

```

Run GNOMIX according to the documentation (https://github.com/AI-sandbox/gnomix). I used genotype data for LAI and so chose the genotype config file from GNOMIX (https://github.com/AI-sandbox/gnomix/blob/main/configs/config_array.yaml) 

Example script:

```
# Define the input and output prefixes from command line arguments
IN_PRE=$1         # Input directory prefix (e.g., /path/to/your/files/)
OUT_PRE=$2        # Output file prefix (e.g., launch_GNOMIX_)

# Get the current date in the format DDMMYY
NOW=$(date +"%d%m%y")

# Define the prefix for configuration files
CONFIG_PRE=${OUT_PRE}_

# Loop over each chromosome  to create separate job submission files for each
for i in {1..22}
do
    # Define the address (file path) for the current configuration file
    CONFIG_ADDR=${CONFIG_PRE}${i}_subfile_

    # Generate a job submission script for each chromosome and save it to ${CONFIG_ADDR}
    cat <<EOF > ${CONFIG_ADDR}
#!/bin/bash

#some example job submission parameters
#BSUB -J GNOMIX_AA_SHORT_$i         # Job name 
#BSUB -P submission_parameter        
#BSUB -q submission_parameter      
#BSUB -n 3     
#BSUB -R "span[hosts=1]"     
#BSUB -R himem              
#BSUB -R rusage[mem=120000]         
#BSUB -W 15:00                 
#BSUB -o ${CONFIG_PRE}${NOW}_${i}_AA_SHORT.out  # job output file
#BSUB -e ${CONFIG_PRE}${NOW}_${i}_AA_SHORT.err  # job error file
#BSUB -L /bin/bash                  # Use bash shell for job execution

# Load necessary modules and activate the required conda environment
ml anaconda3
ml bcftools
source activate /sc/arion/projects/kennylab/roohy/conda/envs/igh_gnomix/

# Execute the gnomix.py script with specified input files and parameters
python3 gnomix.py ${IN_PRE}aa_query_phased_chr_${i}.vcf \              # Input VCF file for the chromosome
                  ${IN_PRE}gnomix_chr${i}_local_ancestry_aa \          # Output file prefix for local ancestry
                  ${i} \                                               # Chromosome number
                  False \                                              # intent to use phasing correction 
                  ${IN_PRE}genetic_map_files/chr${i}.gmap \            # Genetic map file 
                  ${IN_PRE}aa_ref_phased_chr_${i}.vcf \                # Reference file
                  ${IN_PRE}aa_ref_panel.smap \                         # sample map file
                  configs/config_array.yaml                            # Configuration file

EOF

    # Submit the job with the generated script
    bsub < ${CONFIG_ADDR}
done
```

This will generate a directory for each chromosome with local ancestry calls summarized in .msp, .fb and .lai files. I then use the custom R script make_VCF_file_from_GNOMIX_AA.R to convert the local ancestry calls to VCF style format for two-way local ancestry,  make_VCF_file_from_GNOMIX_HL.R does the same for three-way local ancestry calls. The script takes a .msp file as input, in the case of three way local ancestry it will output three different VCF files, one for each local ancestry background (i.e. AFR, EUR, NAT). For two-way local ancestry only one VCF file is output.

Additionally:
* Rscript get_sum_local_ancestry_per_sample_and_lai_genome_density.R:
  - calculates the sum of local ancestry components in bp to plot the correlation between global admixture components and local ancestry proportions for  QC purposes.
  - calculates density of local ancestry calls for each ancestry background, genome-wide for  QC purposes. Plot local ancestry call density for each ancestry background remove the HLA region and then remove  regions that have local ancestry call density +/-3 SD from the median. 
  - script works for both two and three way local ancestry calls

```
Rscript get_sum_local_ancestry_per_sample_and_lai_genome_density.R lai_results_chr21.msp 21

```

* Centromeric regions should also be removed from local ancestry calls, they can be downloaded at: https://genome.ucsc.edu/cgi-bin/hgTables;Group: Mapping and Sequencing, Track: Centromeres, Table Centromeres

Merge vcf files of local ancestry calls and remove outlier samples and regions identified in QC for local ancestry inference.

#### 2) Run admixture mapping: Two-way and three-way ####

Prune genotype data for calculation of GRM:
```
plink2
  --indep-pairwise 500 50 0.2
  --out query_HIS_Shapeit_normalized_Phased_GDA_all_chr_pruned
  --vcf query_HIS_Shapeit_normalized_Phased_GDA_all_chr.vcf.gz
```


Run SAIGE step1, calculating full GRM and fitting null model:

```
NOW=$(date +"%d%m%y")

# mkdir $CONFIG_PRE
while read p;
do
CONFIG_ADDR=subfiles/STEP1_phecode_${p}_subfile_
cat <<EOF > ${CONFIG_ADDR}
#!/bin/bash

#BSUB -J SAIGE_STEP1_HIS_${p}
#BSUB -P acc_kennylab
#BSUB -q premium
#BSUB -M 20000
#BSUB -W 06:00
#BSUB -n 4
#BSUB -o logfiles/SAIGE_STEP1_${p}_${NOW}.out
#BSUB -e logfiles/SAIGE_STEP1_${p}_${NOW}.err
#BSUB -L /bin/bash
ml saige
Rscript /hpc/packages/minerva-centos7/saige/1.1.6/SAIGE/extdata/step1_fitNULLGLMM.R 
--plinkFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/query_HIS_Shapeit_normalized_Phased_GDA_all_chr_sparse_input \
--phenoFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/SAIGE_input_files/his_covariates_plus_phecodes.txt --phenoCol=${p}  \
--covarColList=YOB,SEX,chip,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10  \
--qCovarColList=SEX,chip \
--sampleIDColinphenoFile=MASKED_MRN \
--traitType=binary --outputPrefix=output_STEP1/step1_HIS_SAIGE_pheno_${p}  \
--IsOverwriteVarianceRatioFile=TRUE  \
 --nThreads=4 --skipVarianceRatioEstimation=FALSE  \
 --IsOverwriteVarianceRatioFile=TRUE \
  --isCovariateOffset=TRUE \

EOF
bsub < ${CONFIG_ADDR}
done < /sc/arion/projects/kennylab/Sinead/AncestryWAS/SAIGE_input_files/all_gender_phecode_gtr30_HIS.txt

```

Run SAIGE step 2 for a given local ancestry:

```

NOW=$(date +"%d%m%y")

# mkdir $CONFIG_PRE
while read p;
do
for c in {1..22};
do
CONFIG_ADDR=subfiles/STEP2_anc0_chr_${c}_phecode_${p}_subfile_
cat <<EOF > ${CONFIG_ADDR}
#!/bin/bash

#BSUB -J SAIGE_STEP2_anc0_${p}_${c}
#BSUB -P acc_kennylab
#BSUB -q premium
#BSUB -M 2000
#BSUB -W 04:00
#BSUB -o /sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/logfiles/anc0_log/SAIGE_STEP2_anc0_${p}_${c}_${NOW}.out
#BSUB -e /sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/logfiles/anc0_log/SAIGE_STEP2_anc0_${p}_${c}_${NOW}.err
#BSUB -L /bin/bash

ml saige
cd /sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE
Rscript /hpc/packages/minerva-centos7/saige/1.1.6/SAIGE/extdata/step2_SPAtests.R \
--bedFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/step2_plinkfiles/HIS_ancestry_0_HLA_centromere_density_filtered_chr_${c}.bed \
--bimFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/step2_plinkfiles/HIS_ancestry_0_HLA_centromere_density_filtered_chr_${c}.bim  \
--famFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/step2_plinkfiles/HIS_ancestry_0_HLA_centromere_density_filtered_chr_${c}.fam  \
--AlleleOrder=alt-first  \
--SAIGEOutputFile=output_step2_anc0/HIS_anc0_GNOMIX_step2_no_vr_phecode_${p}_chr${c}.txt \
--chrom=${c} \
--minMAF=0  \
--minMAC=0.5 \
--GMMATmodelFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/output_STEP1/step1_HIS_SAIGE_pheno_${p}.rda  \
--varianceRatioFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/output_STEP1/step1_HIS_SAIGE_pheno_${p}.varianceRatio.txt  \
--LOCO=TRUE  \
--is_Firth_beta=TRUE  \
--pCutoffforFirth=0.05  \
--is_output_moreDetails=TRUE 

EOF
bsub < ${CONFIG_ADDR}
done
done < /sc/arion/projects/kennylab/Sinead/AncestryWAS/SAIGE_input_files/all_gender_phecode_gtr30_HIS.txt

```
* STEAM package for calculating admixture mapping significance thresholds: https://github.com/GrindeLab/STEAM
 - Significance threshold AA: p<1.60x10-05  
 - Significance threshold HL: p<4.88x10-06 

I consolidate results into a single file that has results for all phecodes tested with original SAIGE output columns (with an addional column added to the results file that has the phecode name). I filter down this results file to retain associations below a significance threshold of p<5x10-4 to scan for admixture mapping peaks.  The top position of each AM peak was defined by the first position with the most significant p-value for a given association, and start and stop positions for the peak interval are defined as the positions at which the p-value signal drops below one order of magnitude from the peak. 

This script filters this consolidated results file, identifies the number of unique GWS peaks per phecode per chromosome, defines peak boundaries, and adds some additional annotation. Following use of this script, peaks should be plotted and boundaries manually inspected some admixture mapping peaks can have complex structure that require manual adjustment.

Phecode maps for the phecode definiton used can be downloaded at: https://phewascatalog.org/

```

Rscript process_SAIGE_output.R consolidateded_saige_results_file.txt  phecode_definitions1.2.csv GWS_significant_peaks_annotated.txt

```
#### 3) Fine-map admixture mapping peaks using conditional association testing ####

We next finemap admixture mapping peaks by performing conditional association testing using genotype data of variants falling within the boundaries of the admixture mapping peak. 
- Formula: Phecode ~ Local Ancestry at index LA SNP + Age + Sex + Genotype Chip + PCs1-10 + TEST SNP GENOTYPE X + full GRM

To do this, I made an R script that takes the position range and extracts this region from an input VCF file that contains genotype calls of individuals that were in the discovery admixture mapping cohort. It then converts genotype calls in the interval to a numeric covariate and adds these columns to existing phenotye + covariate information to use as input for conditional association testing in SAIGE. It also outputs a summary df that details the number of snps in the interval and is used later to launch the conditional SAIGE jobs.

```
#./prep_conditional_analysis_input_files.sh <assoc_id> <vcf_filepath> <chr> <start_bp> <stop_bp> <cov_file> <phecode_col_name>
./prep_conditional_analysis_input_files.sh  288.1_HIS_anc2_1 imputed_TOPMed_data.vcf 1 158907131 159795459  his_covariates_plus_phecodes.txt X288.1

```
