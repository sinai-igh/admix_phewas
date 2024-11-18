#!/bin/bash

#BSUB -P acc_kennylab
#BSUB -n 1
#BSUB -R span[hosts=1]
#BSUB -W 4:00
#BSUB -q premium
#BSUB -R rusage[mem=50000]
#BSUB -J "saige[1-217]"
#BSUB -o /sc/arion/projects/igh/kennylab/christa/labwas/logs/saige.out.%J.%I
#BSUB -e /sc/arion/projects/igh/kennylab/christa/labwas/logs/saige.err.%J.%I

. /sc/arion/projects/igh/kennylab/christa/miniforge/etc/profile.d/conda.sh
conda activate saige

############ input files and parameters
imputed_genotypes_prefix="geno/GSA_hg38_QC_FINAL_autosomes"
phenotype_file="final_cov_pheno_ambulatory_normalized.csv"
grm="grm/sparse_imputed"

pheno=$(awk "NR==$LSB_JOBINDEX" final_cov_pheno_ambulatory_pheno_names.txt)
output="output/all_pop/ambulatory/normalized/ambulatory_"$pheno
plot_dir="output/all_pop/ambulatory/normalized/plots"

echo $pheno

#no sparse GRM fitted - fit full GRM for better pop structure correction

########### Step 1 of SAIGE
    Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --useSparseGRMtoFitNULL=FALSE    \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=$pheno"_median_age",GENDER,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --qCovarColList=GENDER  \
        --sampleIDColinphenoFile=IID \
        --invNormalize=FALSE     \
        --traitType=quantitative        \
        --outputPrefix=$output \
        --nThreads=24   \
        --IsOverwriteVarianceRatioFile=TRUE


########### Step 2 of SAIGE
step2_SPAtests.R \
     --bedFile=$imputed_genotypes_prefix".bed" \
     --bimFile=$imputed_genotypes_prefix".bim" \
     --famFile=$imputed_genotypes_prefix".fam"  \
     --AlleleOrder=ref-first \
     --SAIGEOutputFile=$output \
     --minMAF=0 \
     --minMAC=1 \
     --GMMATmodelFile=$output".rda" \
     --varianceRatioFile=$output".varianceRatio.txt" \
     --LOCO=FALSE \
     --is_output_moreDetails=TRUE

########### generate QQ and manhattan plots for SAIGAE output

conda activate /sc/arion/projects/igh/kennylab/christa/miniconda