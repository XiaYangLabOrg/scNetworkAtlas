#load lab functions
  source('r_files/pathway_enrichment_functions.R')

args = commandArgs(trailingOnly=TRUE)
#set user input parameters
input_dir <- args[1]
cell_type <- args[2]
q1 <- args[3]
dbs <- args[4]
output_directory <- args[5]
convertToHuman <- ifelse(args[6]=='True',T,F)
num_cores <- as.numeric(args[7])
# final_directory <- args[6]

print("processed all parameters!")
print(paste0("parameter 1 is ", input_dir))
print(paste0("parameter 2 is ", cell_type))
print(paste0("parameter 3 is ", q1))
print(paste0("parameter 4 is ", dbs))
print(paste0("parameter 5 is ", output_directory))
print(paste0("parameter 6 is ", num_cores))

# print(paste0("parameter 6 is ", final_directory))

cell_type = sapply(strsplit(cell_type, "/"), tail, 1)
print(paste0("cell type is now ", cell_type))

output_directory <- paste0(output_directory, "/", cell_type, "/")
dir.create(output_directory,showWarnings = F)

# modules <- data.frame(matrix(nrow=0,ncol=2))
modules <- read.csv(paste0(input_dir,cell_type,'.gene_membership.',q1,'.csv.gz'), header=TRUE)
modules <- modules[c("cluster_membership", "genes")]
# change the module names
modules$cluster_membership <- paste0(q1,'_',modules$cluster_membership)
# keep modules with >10 genes and <300 genes
module_sizes <- table(modules$cluster_membership)
filtered_modules <- names(module_sizes[which(module_sizes>10 & module_sizes < 300)])
modules <- modules[which(modules$cluster_membership %in% filtered_modules),]

# create pathway enrichment dataframe
t <- system.time({
  makePathwayEnrichmentDf(modules, dbs, output_Dir=output_directory, convertToHuman=convertToHuman, num_cores=num_cores, save_temp_dat=F)
})
cat('time elapsed: \n')
print(t)

# edit pathway enrichment dataframes
pathway_df = data.frame(stringsAsFactors = FALSE)
counter = 0
for(m in unique(modules$cluster_membership)){
  file_name <- paste0(output_directory,m, ".txt")
  temp <- read.delim(file_name, stringsAsFactors = FALSE, header = TRUE)
  pathway_df <- rbind(pathway_df, temp)
}
pathway_df = pathway_df[order(pathway_df$nOverlap, decreasing = TRUE),]
pathway_df$`negLog10FDR` = -log10(pathway_df$FDR)
pathway_df$q1 <- q1

outfile <- paste0(output_directory,q1,"_full.txt")
write.table(pathway_df, outfile, quote=F, sep='\t', row.names=F)

# remove temp files
delete_files <- list.files(output_directory, pattern=paste0(q1,'_'), full.names=T)
delete_files <- delete_files[!grepl('full.txt$',delete_files)]
unlink(delete_files)
