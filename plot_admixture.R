# Accept parameters
working_dir <- args[1]
ordered_admixture <- args[2]  # e.g., /sc/private/regen/data/GSA_GDA/imputed_TOPMED_V2/addition/PC/GSA_GDA_PCA_V2.txt
output_plot <- args[3]
admixture_type <- args[4]  # Accept either "K3" or "K4"

library(data.table)

# Set working directory
setwd(working_dir)

# Load data
results_ordered <- fread(ordered_admixture)

# List of unique ethnic group codes to avoid redundancy
pop_groups <- unique(results_ordered$population)

# Determine the column indices based on admixture type (K3 or K4)
if (admixture_type == "K3") {
  col_start <- 3
  col_end <- 5  # For K3, we use 3 columns (V3 to V5)
  color_scheme <- c('#335C67', "#7EA16B", "#9E2A2B")  # Adjusted for K3 (3 colors)
  y_label <- "K3"
} else if (admixture_type == "K4") {
  col_start <- 3
  col_end <- 6  # For K4, we use 4 columns (V3 to V6)
  color_scheme <- c('#335C67', "#7EA16B", "#9E2A2B", '#E09F3E')  # For K4 (4 colors)
  y_label <- "K4"
} else {
  stop("Please provide a valid admixture type: 'K3' or 'K4'")
}

# Create a named list of matrices
group_matrices <- lapply(pop_groups, function(group) {
  as.matrix(results_ordered[results_ordered$V2 == group, col_start:col_end])
})

# Function to generate a barplot and save as TIFF
generate_barplot <- function(data, filename, color_scheme, y_label) {
  tiff(filename, width = 1800, height = 300)
  barplot(t(data), border = NA, space = 0, col = color_scheme, ylab = y_label, xaxt = 'n', cex.lab = 1.5)
  dev.off()
}

# Placeholder matrix to create gaps between samples
x <- matrix(data = 0, nrow = 10, ncol = (col_end - col_start + 1))  # Adjust size based on K3/K4 

# Combine matrices for the plot (example uses Thousand Genome populations, you can modify as needed)
plot_population_order <- rbind(
  group_matrices$ACB, x, x, group_matrices$ASW, x, x, group_matrices$YRI, x, x,
  group_matrices$ESN, x, x, group_matrices$MSL, x, x, group_matrices$GWD, x, x,
  group_matrices$LWK, x, x, x, x, group_matrices$BEB, x, x, group_matrices$STU, x, x,
  group_matrices$PJL, x, x, group_matrices$ITU, x, x, group_matrices$GIH, x, x, x, x,
  group_matrices$CDX, x, x, group_matrices$CHB, x, x, group_matrices$CHS, x, x,
  group_matrices$JPT, x, x, group_matrices$KHV, x, x, x, x, group_matrices$CEU,
  x, x, group_matrices$TSI, x, x, group_matrices$FIN, x, x, group_matrices$GBR,
  x, x, group_matrices$IBS, x, x, x, x, group_matrices$CLM, x, x, group_matrices$PUR,
  x, x, group_matrices$PEL, x, x, group_matrices$MXL
)

# Generate the plot
generate_barplot(plot_population_order, as.character(output_plot), color_scheme, y_label)
