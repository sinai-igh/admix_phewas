
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
Rscript /hpc/packages/minerva-centos7/saige/1.1.6/SAIGE/extdata/step1_fitNULLGLMM.R --plinkFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/query_HIS_Shapeit_normalized_Phased_GDA_all_chr_sparse_input --phenoFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/SAIGE_input_files/his_covariates_plus_phecodes.txt --phenoCol=${p} --covarColList=YOB,SEX,chip,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --qCovarColList=SEX,chip --sampleIDColinphenoFile=MASKED_MRN --traitType=binary --outputPrefix=output_STEP1/step1_HIS_SAIGE_pheno_${p} --IsOverwriteVarianceRatioFile=TRUE  --nThreads=4 --skipVarianceRatioEstimation=FALSE --IsOverwriteVarianceRatioFile=TRUE --isCovariateOffset=TRUE

EOF
bsub < ${CONFIG_ADDR}
done < /sc/arion/projects/kennylab/Sinead/AncestryWAS/SAIGE_input_files/all_gender_phecode_gtr30_HIS.txt
