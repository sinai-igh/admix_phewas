geno_path=/sc/arion/projects/kennylab/Sinead/LAI_GDA_GSA/HIS_AA_imputed_merge/PHASING/

for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}

do
echo 'cd '${geno_path}'' > Launch_shapeit_${i}.pbs
echo 'module load shapeit/v2r900' >> Launch_shapeit_${i}.pbs
echo 'gunzip Phased_HIS_AFR_plus_Ref_no_2nd_deg_rel'${i}'.haps.gz' >> Launch_shapeit_${i}.pbs
echo 'shapeit -convert --input-haps Phased_HIS_AFR_plus_Ref_no_2nd_deg_rel'${i}' --output-vcf Phased_HIS_AFR_plus_Ref_no_2nd_deg_rel'${i}'.vcf' >> Launch_shapeit_${i}.pbs


bsub -q premium -P acc_kennylab -W 12:00 -M 70000 -o Launch_shapeit_${i}.log < Launch_shapeit_${i}.pbs


done