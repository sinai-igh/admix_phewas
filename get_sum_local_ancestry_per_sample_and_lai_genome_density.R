#author Sinead Cullina
#gets sum of bp called as each local ancestry component from GNOMIX .msp file
#works regardless of n ancestries called

# Load required packages
library(data.table)
library(dplyr)

# Parse command line arguments 
args <- commandArgs(TRUE)
chr <- args[1]
msp_filepath <-args[2]


# Load input file for the given chromosome
my_data <- fread(msp_filepath)

# Calculate base-pair length (n_BP) for each row
my_data$n_BP <- my_data$epos - my_data$spos

# Extract relevant column names for ancestry data
ancestry_cols <- colnames(my_data)[7:(ncol(my_data)-1)]
ind_b <- ancestry_cols[seq(2, length(ancestry_cols), by = 2)]
ind_a <- ancestry_cols[seq(1, length(ancestry_cols), by = 2)]

# Split data into two ancestry groups and combine
ancestry_a <- cbind(my_data$n_BP, my_data[, ind_a, with = FALSE])
ancestry_b <- cbind(my_data$n_BP, my_data[, ind_b, with = FALSE])

# Set column names of both dataframes to match
colnames(ancestry_b) <- colnames(ancestry_a)

# Combine the two dataframes by rows
combined_data <- rbind(ancestry_b, ancestry_a)
colnames(combined_data)[1] = "n_BP"


# Summarize ancestry proportions for each individual
for (i in 2:ncol(combined_data)) {
  # Group by ancestry information and calculate the sum of base-pair lengths
  summary_df <- combined_data %>%
    group_by(combined_data[, ..i]) %>%
    summarise(sum_BP = sum(n_BP, na.rm = TRUE))
  
  # Add metadata columns to the summary
  sample_id <- colnames(summary_df)[1]
  colnames(summary_df)[1] <- "ID"
  summary_df$sample_name <- sample_id
  summary_df$total_BP <- sum(combined_data$n_BP)  # Total base-pair length
  summary_df$sample_name = sapply(strsplit(summary_df$sample_name,"\\."), `[`, 1)

  # Write the summary to an output file
  output_file <- paste0("ANCESTRY_sum_check", chr, ".txt")
  write.table(
    summary_df, file = output_file, 
    quote = FALSE, sep = "\t", 
    row.names = FALSE, col.names = FALSE, append = TRUE
  )
}

my_data = my_data[,-c("n_BP")] #remove n_BP column
#write frequency of each local ancestry call per position - this should be used to plot the density of calls for each ancestry component genome

for (i in 1:nrow(my_data)){
summary_df = cbind( as.data.frame(table(t(my_data[i,7:ncol(my_data)]))),  my_data[i,1:5])
write.table(summary_df , paste("ancestry_density_" , chr, ".txt", sep = ""), quote= F, sep = "\t", row.names = F, col.names = F, append = T)
}

