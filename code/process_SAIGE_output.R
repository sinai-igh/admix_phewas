args <- commandArgs(TRUE)

# Load necessary libraries
require(data.table)  # For efficient data manipulation
require(dplyr)       # For data manipulation (especially useful for grouping and summarizing data)


association_results= args[1] #path to results file
phecode_mapper= args[2] # e.g. phecode_definitions1.2.csv
outfile = arg[3] #output filename

# Read in main association result file
my_file <- fread(association_results)  # Update file path if needed

# Name SAIGE columns (adjust if neccesary)
colnames(my_file) <- c("CHR", "POS", "MarkerID", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2", 
                       "MissingRate", "BETA", "SE", "Tstat", "var", "p.value", "p.value.NA", 
                       "Is.SPA", "AF_case", "AF_ctrl", "N_case", "N_ctrl", "N_case_hom", 
                       "N_case_het", "N_ctrl_hom", "N_ctrl_het", "phecode")

# Load phecode definitions for phenotype mapping
mapper <- fread(phecode_mapper, colClasses = "character")

# Merge association results with phenotype definitions
out1 <- merge(my_file, mapper, by.x = "phecode", by.y = "phecode")

#make sure POS and p.value variables are numeric
out1$POS <- as.numeric(out1$POS) #make suer POS variable is numeric
out1$p.value <- as.numeric(out1$p.value)

# Identify phecode chromosome pairs that have more than two GWS peaks
# Filter for peaks where p-value <= 1.60e-05
peaks_below_threshold <- out1[out1$p.value <= 1.60e-05  , ]

# Group peaks by contiguous regions based on p-value threshold
# Here we define contiguous peaks based on non-overlapping position intervals
# using `rle` (run-length encoding) to count separate regions in each chromosome

# Define a new column to indicate separate peak regions
peaks_below_threshold <- peaks_below_threshold %>%
  arrange(CHR, POS) %>%
  mutate(Region = cumsum(c(1, diff(POS) > 2e6)))  # Define separate regions by distance (e.g., 2 Mb gap) or other threshold

# Count the number of distinct regions (peaks) per phecode and chromosome
peak_counts <- peaks_below_threshold %>%
  group_by(phecode, CHR, Region) %>%
  summarize(n_peaks = n_distinct(Region), .groups = "drop") %>%
  group_by(phecode, CHR) %>%
  summarize(total_peaks = n_distinct(Region), .groups = "drop")

# Filter for entries with more than two peaks passing the threshold
multi_peak_phecodes <- peak_counts[peak_counts$total_peaks > 2, ]

# Output multi-peak entries for review
print("Phecodes and chromosomes with more than two peaks passing the 5 x 10^-4 threshold:")
print(multi_peak_phecodes)

### Our results only had max one peak per chromosome per phecode

# Identify minimum p-value by phecode and chromosome
out1 <- out1 %>%
  group_by(phecode, CHR) %>%
  mutate(Minp.val = min(p.value, na.rm = TRUE))

# Keep records with p-values below one order of magnitude below peak
check2 <- out1[out1$p.value <= (out1$Minp.val * 10),]

# Calculate peak boundaries (min and max position per chromosome-phecode combination) #this is contingent on there being only one peak per phecode per chromosome
out_sig <- check2 %>%
  group_by(phecode, CHR) %>%
  mutate(MaxPOSbyChrPhe = max(POS, na.rm = TRUE),
         MinPOSByChrPhe = min(POS, na.rm = TRUE)) %>%
  as.data.frame()

# Retain only minimum p-value rows per chromosome-phecode group
out_sig <- out_sig[out_sig$p.value == out_sig$Minp.val,]

# Identify points where p-value changes to have one line per result 
switch_points <- which(out_sig$p.value != dplyr::lag(out_sig$p.value))
switch_points1 <- c(1, switch_points, nrow(out_sig))
out_sig <- out_sig[switch_points1,]

out_sig = out_sig[out_sig$Minp.val <= 1.60e-05  ,]
# Add odds ratio (OR) and confidence intervals (CI) for interpretation
out_sig$OR <- round(exp(out_sig$BETA), 2)
out_sig$Upper_CI <- round(exp(as.numeric(out_sig$BETA) + (1.96 * as.numeric(out_sig$SE))), 2)
out_sig$Lower_CI <- round(exp(as.numeric(out_sig$BETA) - (1.96 * as.numeric(out_sig$SE))), 2)

# Order by minimum p-value and create case-control ratio and genomic region information
out_sig <- out_sig[order(out_sig$Minp.val),]
out_sig$ratio_case_cont <- paste(out_sig$N_case, out_sig$N_ctrl, sep = ":")
out_sig$min_pos_collapsed <- paste(out_sig$CHR, ":", round(out_sig$MinPOSByChrPhe / 1e6, 3), "-", 
                                    round(out_sig$MaxPOSbyChrPhe / 1e6, 3), sep = "")
out_sig$OR_CI <- paste(out_sig$OR, " (", out_sig$Lower_CI, "-", out_sig$Upper_CI, ")", sep = "")
out_sig <- unique(out_sig)

out_sig$id <- paste(out_sig$phecode, "AA", out_sig_cond$CHR, sep = "_") #add id for association

# Write final processed table to output file
write.table(out_sig, outfile, quote = FALSE, sep = "\t", row.names = FALSE)

