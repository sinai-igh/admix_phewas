IN_PRE=$1 #/sc/arion/projects/kennylab/Sinead/LAI_GDA_GSA/HIS_AA_imputed_merge/GNOMIX_LAI/
OUT_PRE=$2 #launch_GNOMIX_

NOW=$(date +"%d%m%y")

CONFIG_PRE=${OUT_PRE}_



#convert_gnomix_to_vcf.R 
mkdir $CONFIG_PRE
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

