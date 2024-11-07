#!/bin/bash

# Usage:
# ./run_analysis.sh <assoc_id> <vcf_filepath> <chr> <start_bp> <stop_bp> <cov_file> <phecode_col_name>

# Check for required arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <assoc_id> <vcf_filepath> <chr> <start_bp> <stop_bp> <cov_file> <phecode_col_name>"
    exit 1
fi

# Assign arguments to variables
assoc_id=$1
vcf_filepath=$2
chr=$3
start_bp=$4
stop_bp=$5
cov_file=$6
phecode_col_name=$7

# Paths and filenames
output_vcf="${assoc_id}_conditional_assoc_interval"
output_file="${output_vcf}.vcf"

# Step 1: Run PLINK command to export VCF file
plink --bfile "${vcf_filepath}" \
      --chr "${chr}" \
      --export vcf \
      --from-bp "${start_bp}" \
      --to-bp "${stop_bp}" \
      --out "${output_vcf}"

# Check if PLINK command was successful
if [ $? -ne 0 ]; then
    echo "Error: PLINK command failed."
    exit 1
fi

# Step 2: Run the R script to process the VCF file
Rscript -e "
# Load necessary libraries
require(data.table)

# Define input arguments
invcf <- '${output_file}'
assoc_id <- '${assoc_id}'
covariate_file_SAIGE <- '${cov_file}'
in_phe <- '${phecode_col_name}' 

# Read and preprocess VCF file
my_file <- fread(invcf)
my_file <- my_file[, -c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')]

# Create SNP names and transpose data
my_file <- cbind(paste0('snp', 1:nrow(my_file)), my_file)
colnames(my_file)[1] <- 'snpname'
my_file <- as.data.frame(t(my_file))
colnames(my_file) <- my_file[1, ]
my_file <- my_file[-1, ]

# Convert genotypes to numeric values
my_file_conv1 <- as.data.frame(lapply(my_file, function(x) gsub('1/1', '2', x)))
my_file_conv1 <- as.data.frame(lapply(my_file_conv1, function(x) gsub('0/0', '0', x)))
my_file_conv1 <- as.data.frame(lapply(my_file_conv1, function(x) gsub('0/1', '1', x)))
my_file_conv1 <- as.data.frame(lapply(my_file_conv1, function(x) gsub('1/0', '1', x)))
my_file_conv1 <- as.data.frame(lapply(my_file_conv1, function(x) gsub('\\./\\.', 'NA', x)))

# Add MRN identifier
my_file_conv1 <- cbind(MRN = rownames(my_file), my_file_conv1)

# Load covariate and phenotype data
covariate_file <- fread(covariate_file_SAIGE)

# pull covariates and single phecode variable for conditional analysis
phecode <- as.character(in_phe)
covariate_file <- cbind(covariate_file[, 1:15], covariate_file[[phecode]])

#splitting of MRN IDs to match IDs in covariate file because in VCF file they were informat ID_ID (may not be applicable)
my_file_conv1\$MRN <- sapply(strsplit(as.character(my_file_conv1\$MRN), '_'), '[', 1)

# Create SNP metadata matrix that writes out to dataframe that has information on number of snps being tested in the interval
holder_matrix <- matrix(nrow = (ncol(my_file_conv1) - 1), ncol = 4)
holder_matrix[, 1] <- paste0('snp', 1:(ncol(my_file_conv1) - 1))
holder_matrix[, 2] <- assoc_id
holder_matrix[, 3] <- sapply(strsplit(assoc_id, '_'), '[', 1)
holder_matrix[, 4] <- sapply(strsplit(assoc_id, '_'), '[', 4)

# Write SNP metadata and covariate files
write.table(holder_matrix, paste0('nsnps_in_interval', assoc_id, '.txt'), quote = FALSE, sep = '\\t', row.names = FALSE, col.names = FALSE)

covariate_file\$MASKED_MRN <- as.character(covariate_file\$MASKED_MRN)
covariate_file <- unique(merge(covariate_file, my_file_conv1, by.x = 'MASKED_MRN', by.y = 'MRN'))
colnames(covariate_file)[5] <- 'FID'
write.table(covariate_file, paste0('cov_and_pheno_file_', assoc_id, '.txt'), quote = FALSE, sep = '\\t', row.names = FALSE)
"

# Check if R script was successful
if [ $? -ne 0 ]; then
    echo "Error: R script failed."
    exit 1
fi

echo "Analysis completed successfully."
