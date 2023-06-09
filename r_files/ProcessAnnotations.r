args = commandArgs(trailingOnly=TRUE)
cell_type <- args[1]
annotations_directory <- args[2]
output_directory <- args[3]

# get a list of all the module files for a specific resolution
intermediate_directory <- paste0(annotations_directory, "/", cell_type, "/")
output_directory <- paste0(args[3], "/", cell_type, "/")
list_of_annotation_files = list.files(path = intermediate_directory, pattern = paste0("full.txt"), full.names = TRUE)

# make raw + filtered files
raw_annotations <- data.frame()
filtered_annotations <- data.frame()
for(i in 1:(length(list_of_annotation_files))) {
    
    #add res
    enriched <- read.table(list_of_annotation_files[i], header=TRUE)
    
    #save raw annotations
    raw_annotations <- rbind(raw_annotations, enriched)
    
#     #save filtered annotations
#     enriched <- enriched[(enriched$FDR < 0.05),]
#     enriched <- enriched[(enriched$nOverlap > 4),]
#     filtered_annotations <- rbind(filtered_annotations, enriched) 
}

# raw_annotations$FDR <- p.adjust(raw_annotations$FDR, method='BH')
# raw_annotations$negLog10FDR <- -log10(raw_annotations$FDR)
filtered_annotations <- raw_annotations[raw_annotations$FDR<=0.05 & raw_annotations$nOverlap>=5,]
filtered_annotations <- filtered_annotations[order(filtered_annotations$FDR),]

#save the raw and filtered annotations
dir.create(output_directory, showWarnings=F)
write.table(raw_annotations, paste0(output_directory, "/RAW.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(filtered_annotations, paste0(output_directory, "/FILTERED.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

# #remove all module + .dat + heat map files
# unlink(output_directory, recursive=TRUE)

