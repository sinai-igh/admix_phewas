admixture_values= args[1] #Admixture. Q file
plink_fam= args[2] #admixture input fam file to map IDs to admixture components .Q file
mapper_file= args[3] #mapper with 3 colums. col1= ID, col2= population name , col3= cohort name e.g. "HGDP" or "Query"
output_file = args[4]

library(data.table)
library(dplyr)
setwd("~/Dropbox/Data_ADMIXTURE-PheWAS/Race_Ethnicity/Admixture_AA_and_HIS_GDA_GSA_Impute_TGP_HGDP_PAGE/")

fread(admixture_values)
fam_file = fread(plink_fam)
fam_file = fam_file[,1:2]
Admixture_res = cbind(fam_file, Admixture_res)
colnames(Admixture_res)[1:2] = c("ID1", "ID2")

## merge mapper file with admixture components
mapper = fread(mapper_file)
Admixture_res = merge(Admixture_res, mapper, by.x= "ID2" , by.y = "ID", all.x = T)

##plotting admixture results will help you to understand which admixture components in your query data are most similar to the reference panels. To plot neatly we need to order the file. 
## this next section orders the samples according to the major component within each population group.

# Process the data and populate the master list
for (i in 1:nrow(Admixture_res)) {
  pop_id <- as.character(Admixture_res[i, 7]) #column 7 refers to the population label column within which we want to sort
  
  if (!pop_id %in% names(master)) {
    master[[pop_id]] <- list(c1 = c(), c2 = c(), c3 = c(), c4 = c())
  }
  
  master[[pop_id]]$c1 <- c(master[[pop_id]]$c1, as.numeric(Admixture_res[i, 3]))
  master[[pop_id]]$c2 <- c(master[[pop_id]]$c2, as.numeric(Admixture_res[i, 4]))
  master[[pop_id]]$c3 <- c(master[[pop_id]]$c3, as.numeric(Admixture_res[i, 5]))
  master[[pop_id]]$c4 <- c(master[[pop_id]]$c4, as.numeric(Admixture_res[i, 6]))
}

# Calculate mean proportions and determine the dominant component for each population
populations <- list()

for (pop_id in names(master)) {
  temp <- c(
    c1 = mean(as.numeric(master[[pop_id]]$c1)),
    c2 = mean(as.numeric(master[[pop_id]]$c2)),
    c3 = mean(as.numeric(master[[pop_id]]$c3)),
    c4 = mean(as.numeric(master[[pop_id]]$c4))
  )
  
  # Find the dominant component (the one with the maximum mean value)
  top_component <- names(temp)[which.max(temp)]
  populations[[pop_id]] <- top_component
}

# Open the output file for writing the ordered results
output_file <- file(output_file, "w")

# Iterate over the population order and write ordered results
for (pop in order) {
  print(pop)
  
  temp <- list()
  
  for (i in 1:nrow(Admixture_res)) {
    pop_id <- as.character(Admixture_res[i, 7])
    
    if (pop_id == pop) {
      if (pop_id %in% names(populations)) {
        top <- populations[[pop_id]]
        # Extract the column number corresponding to the dominant component
        w <- switch(top,
                    "c1" = 3,
                    "c2" = 4,
                    "c3" = 5,
                    "c4" = 6)
        
        temp[[paste(admixture_data[i, ], collapse = "\t")]] <- as.numeric(admixture_data[i, w])
      } else {
        print(paste(pop_id, 'why not?'))
      }
    }
  }
  
  # Sort the temp list by the values of the dominant component
  sorted_temp <- temp[order(unlist(temp))]
  
  # Write the sorted data to the output file
  for (line in names(sorted_temp)) {
    writeLines(line, output_file)
  }
}

# Close the output file
close(output_file)
