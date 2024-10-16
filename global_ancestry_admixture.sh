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

