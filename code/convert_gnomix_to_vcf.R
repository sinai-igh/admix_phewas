##args <- commandArgs(TRUE)
inner=args[1]
position_info= args[2]
chr= args[3]

require(data.table)
my_file =fread(inner)
positions_file = fread(position_info)
positions_file$FORMAT = "GT"

output_df <- data.frame()

i = 1
for (entry in 1:nrow(positions_file)){
pos = as.numeric(positions_file[entry,][,2])
pos_output= my_file[(my_file$spos <= pos) & (my_file$epos >= pos),]
output_df  <- rbind(output_df, pos_output)
i = i+1
}

output_df=output_df[,-c(1:6)]

sample_order = unique(sapply(strsplit(colnames(output_df),"\\."), `[`, 1))

matrix_holder = matrix(nrow = nrow(output_df), ncol = (ncol(output_df)/2))
i = 1


for (j in seq(2, ncol(output_df),2)){
 col1 = colnames(output_df)[j-1]
 col2 = colnames(output_df)[j]
 matrix_holder[,i] = paste(output_df[[col1]] ,"|", output_df[[col2]], sep= "")
 i = i+1
 }

matrix_holder = as.data.frame(matrix_holder)
#Ancestry0_conversion

matrix_holder_anc0 <- data.frame(lapply(matrix_holder, function(x) {
gsub("0", "3", x)
}))

matrix_holder_anc0 <- data.frame(lapply(matrix_holder_anc0, function(x) {
gsub("2", "0", x)
}))

matrix_holder_anc0 <- data.frame(lapply(matrix_holder_anc0, function(x) {
gsub("1", "0", x)
}))

matrix_holder_anc0 <- data.frame(lapply(matrix_holder_anc0, function(x) {
gsub("3", "1", x)
}))

colnames(matrix_holder_anc0) = sample_order
matrix_holder_anc0 = cbind(positions_file, matrix_holder_anc0)
write.table(matrix_holder_anc0, paste("chr", chr, "_ancestry_0.vcf", sep= ""), quote= F, sep= "\t", row.names = F)


#Ancestry1_conversion
matrix_holder_anc1 <- data.frame(lapply(matrix_holder, function(x) {
gsub("2", "0", x)
}))

colnames(matrix_holder_anc1) = sample_order
matrix_holder_anc1 = cbind(positions_file, matrix_holder_anc1)
write.table(matrix_holder_anc1, paste("chr", chr, "_ancestry_1.vcf", sep= ""), quote= F, sep= "\t", row.names = F)


#Ancestry2_conversion

matrix_holder_anc2 <- data.frame(lapply(matrix_holder, function(x) {
gsub("1", "0", x)
}))

matrix_holder_anc2 <- data.frame(lapply(matrix_holder_anc2, function(x) {
gsub("2", "1", x)
}))

colnames(matrix_holder_anc2) = sample_order
matrix_holder_anc2 = cbind(positions_file, matrix_holder_anc2)

write.table(matrix_holder_anc2, paste("chr", chr, "_ancestry_2.vcf", sep= ""), quote= F, sep= "\t", row.names = F)
