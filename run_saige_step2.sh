NOW=$(date +"%d%m%y")

# mkdir $CONFIG_PRE
while read p;
do
for c in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
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
Rscript /hpc/packages/minerva-centos7/saige/1.1.6/SAIGE/extdata/step2_SPAtests.R --bedFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/step2_plinkfiles/HIS_ancestry_0_HLA_centromere_density_filtered_chr_${c}.bed --bimFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/step2_plinkfiles/HIS_ancestry_0_HLA_centromere_density_filtered_chr_${c}.bim --famFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/step2_plinkfiles/HIS_ancestry_0_HLA_centromere_density_filtered_chr_${c}.fam --AlleleOrder=alt-first  --SAIGEOutputFile=output_step2_anc0/HIS_anc0_GNOMIX_step2_no_vr_phecode_${p}_chr${c}.txt --chrom=${c} --minMAF=0 --minMAC=0.5 --GMMATmodelFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/output_STEP1/step1_HIS_SAIGE_pheno_${p}.rda --varianceRatioFile=/sc/arion/projects/kennylab/Sinead/AncestryWAS/HIS_imputed_SAIGE/output_STEP1/step1_HIS_SAIGE_pheno_${p}.varianceRatio.txt --LOCO=TRUE --is_Firth_beta=TRUE  --pCutoffforFirth=0.05  --is_output_moreDetails=TRUE

EOF
bsub < ${CONFIG_ADDR}
done
done < /sc/arion/projects/kennylab/Sinead/AncestryWAS/SAIGE_input_files/all_gender_phecode_gtr30_HIS.txt



