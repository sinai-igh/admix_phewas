admixture_values= args[1] #ath to Admixture .Q file with admixture components
plink_fam= args[2] #Path to .fam file to map IDs to the admixture components (in my case column 2 had the ID i chose hence why ID2 is used for merging later)
mapper_file= args[3] #mapper with 3 colums. col1= ID, col2= population name , col3= cohort name e.g. "HGDP" or "Query"
output_file = args[4] #Path to output file

#load libraries
library(data.table)
library(dplyr)

#load files
Admixture_res = fread(admixture_values)
fam_file = fread(plink_fam)

# Combine ID columns from .fam file with the admixture results to map IDs to admixture components
fam_file = fam_file[,1:2]
Admixture_res = cbind(fam_file, Admixture_res)
colnames(Admixture_res)[1:2] = c("ID1", "ID2")

# Merge the mapper data with the admixture data based on matching IDs (ID2 in admixture and ID in mapper)
mapper = fread(mapper_file)
Admixture_res = merge(Admixture_res, mapper, by.x= "ID2" , by.y = "ID", all.x = T)

##plotting admixture results will help you to understand which admixture components in your query data are most similar to the reference panels. To plot neatly we need to order the file. 
## this next section orders the samples according to the major component within each population group.

# Filter out any rows where the "pop_name" column has missing data, as they are not needed for ordering
Admixture_res = Admixture_res[complete.cases(Admixture_res$pop_name),]

# Retrieve unique pop names from the "pop_name" column to define the population order
order= unique(Admixture_res$pop_name)

# Initialize an empty list for master
master <- list()

# Iterate through each row in Admixture_res to populate the master list by population
for (i in 1:nrow(Admixture_res)) {
  pop_id <- as.character(Admixture_res[i, 6]) #column 6 refers to the population label column within which we want to sort
  
  # If population ID doesn't already exist in master, initialize it with empty lists for each component
  if (!pop_id %in% names(master)) {
    master[[pop_id]] <- list(c1 = c(), c2 = c(), c3 = c())
  }
  
  # Append the admixture component values to the lists for this population
  master[[pop_id]]$c1 <- c(master[[pop_id]]$c1, as.numeric(Admixture_res[i, 3]))
  master[[pop_id]]$c2 <- c(master[[pop_id]]$c2, as.numeric(Admixture_res[i, 4]))
  master[[pop_id]]$c3 <- c(master[[pop_id]]$c3, as.numeric(Admixture_res[i, 5]))
}

# Initialize a list to store the dominant admixture component for each population
populations <- list()

# Calculate the mean admixture component values and identify the dominant component for each population
for (pop_id in names(master)) {
  temp <- c(
    c1 = mean(as.numeric(master[[pop_id]]$c1)),
    c2 = mean(as.numeric(master[[pop_id]]$c2)),
    c3 = mean(as.numeric(master[[pop_id]]$c3))
  )
  
  # Find the dominant component (the one with the maximum mean value)
  top_component <- names(temp)[which.max(temp)]
  populations[[pop_id]] <- top_component # Store the dominant component for the population
}

# Open the output file for writing the ordered results
output_file <- file(output_file, "w")

Admixture_res = as.data.frame(Admixture_res)# Convert Admixture_res to a data frame for simpler standard R indexing in the following code
# Iterate over the population order and write ordered results
for (pop in order) {
  print(pop)  # Display current population being processed
  
  temp <- list() # Temporary list to store lines to be ordered
  
  # Loop over each row in Admixture_res to find rows belonging to the current population
  for (i in 1:nrow(Admixture_res)) {
    pop_id <- as.character(Admixture_res[i, 6]) # Population label in current row
   
     # If the row's population matches the current population in the loop
    if (pop_id == pop) {
      # Check if a dominant component exists for the current population
      if (pop_id %in% names(populations)) {
        top <- populations[[pop_id]]  # Get the dominant component for this population
        # Extract the column number corresponding to the dominant component
        w <- switch(top,
                    "c1" = 3,
                    "c2" = 4,
                    "c3" = 5)
        # Store the line and value of the dominant component for sorting
        temp[[paste(Admixture_res[i, ], collapse = "\t")]] <- as.numeric(Admixture_res[i, w])
      } else {
        print(paste(pop_id, 'why not?'))
      }
    }
  }
  
  # Sort the temp list by the values of the dominant component for the current population
  sorted_temp <- temp[order(unlist(temp))]
  
  # Write the sorted data to the output file
  for (line in names(sorted_temp)) {
    writeLines(line, output_file)
  }
}

# Close the output file
close(output_file)
