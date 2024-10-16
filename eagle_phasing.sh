geno_path=/sc/arion/projects/kennylab/Sinead/LAI_GDA_GSA/HIS_AA_imputed_merge/PHASING/
for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}

do
echo 'cd '${geno_path}'' > Launch_Eagle_${i}.pbs
echo 'module load eagle/2.4' >> Launch_Eagle_${i}.pbs
echo 'eagle --bfile HIS_AFR_plus_Ref_no_2nd_deg_rel --geneticMapFile /hpc/packages/minerva-common/eagle/2.4/Eagle_v2.4/tables/genetic_map_hg38_withX.txt.gz --chrom '${i}' --outPrefix Phased_HIS_AFR_plus_Ref_no_2nd_deg_rel'${i}'' >> Launch_Eagle_${i}.pbs
## if premium doesn't work for you, you may want to try changing to -q alloc

bsub -q premium -P acc_kennylab -n 10 -R span[ptile=10] -R rusage[mem=1200] -W 60:00 -o Launch_Eagle_${i}.log < Launch_Eagle_${i}.pbs

done