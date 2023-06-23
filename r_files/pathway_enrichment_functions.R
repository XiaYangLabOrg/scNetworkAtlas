# not good commenting right now and flexibility. will update later

library(ggplot2)
library(foreach)
library(doParallel)
# library(biomaRt)

tool.read <- function(file, vars=NULL) {
  if(is.null(file)) return(data.frame())
  if(file == "") return(data.frame())
  dat <- read.delim(file=file, header=TRUE,
                    na.strings=c("NA", "NULL", "null", ""),
                    colClasses="character", comment.char="",
                    stringsAsFactors=FALSE)
  if(is.null(vars) == FALSE) dat <- dat[,vars]
  dat <- na.omit(dat)
  return(dat)
}

# must have "GENE" column
# single-cell RNA-seq gene name convention applies 
#       -so mitochondrial genes are "mt-"/"MT-"
# if converting mouse genes for pathway enrichment, set forPathway = TRUE
#       This is because pathway databases will represent MT-ND1 as ND1, for example
convertDfGeneColumnMouseHuman <- function(df, toSpecies="human", forPathway=FALSE){
  MT_gene_names = c("MT-ND1"="ND1",  "MT-ND2"="ND2",  "MT-CO1"="COX1",  "MT-CO2"="COX2",  "MT-ATP8"="ATP8", 
                    "MT-ATP6"="ATP6", "MT-CO3"="COX3","MT-ND3" ="ND3", "MT-ND4L"="ND4L", "MT-ND4"="ND4",
                    "MT-ND5"="ND5",  "MT-ND6"="ND6",  "MT-CYB"="CYB")
  if(toSpecies=="mouse"){
    convertedtoMouse <- convertMouseGeneList(df$GENE)
    df$MOUSE = convertedtoMouse$MGI.symbol[match(df$GENE, convertedtoMouse$HGNC.symbol)]
    df$MOUSE[which(is.na(df$MOUSE))] <- tolower(df$GENE[which(is.na(df$MOUSE))])
    df$GENE <- df$MOUSE
    df$GENE <- paste0(toupper(x=substr(df$GENE, start = 1, stop = 1)),tolower(substring(df$GENE, first = 2)))
    # correct for mt- genes
    MT_mouse_gene_names = MT_gene_names
    new_names = c()
    for(name in names(MT_gene_names)){
      new_names = append(new_names, 
                         paste0(tolower(unlist(strsplit(name, split = "-"))[1]), "-",
                                paste0(toupper(x=substr(unlist(strsplit(name, split = "-"))[2], start = 1, stop = 1)),
                                       tolower(substring(unlist(strsplit(name, split = "-"))[2], first = 2)))))
    }
    names(MT_mouse_gene_names) = new_names
    
    for(gene in 1:nrow(df)){
      if(sum(MT_mouse_gene_names==df$gene[gene])>0){
        df$GENE[gene] = names(MT_mouse_gene_names)[MT_mouse_gene_names==df$gene[gene]]
      }
      else{
        next
      }
    }
    return(df)
  }
  else if(toSpecies=="human"){
    convertedToHuman <- convertMouseGeneList(df$GENE)
    df$HUMAN <- convertedToHuman$HGNC.symbol[match(df$GENE, convertedToHuman$MGI.symbol)]
    df$HUMAN[which(is.na(df$HUMAN))] <- toupper(df$GENE[which(is.na(df$HUMAN))])
    df$GENE <- df$HUMAN
    df$HUMAN <- NULL
    if(sum(is.na(df$GENE))>0) cat("Warning: NAs created.\n")
    if(forPathway){ # convert MT- genes
      for(gene in 1:nrow(df)){
        if(grepl("^MT-",df$GENE[gene])){
          df$GENE[gene] = MT_gene_names[df$GENE[gene]]
        }
        else{
          next
        }
      }
    }
    return(df)
  }
  else cat("Put a valid toSpecies value - 'human' or 'mouse'\n")
}

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  #Added archived database as host. DAN - 08.02.2022 
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , 
                   mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  # Print the first 6 genes found to the screen
  print(head(genesV2))
  return(genesV2)
}

# there may be bugs here
pathway_df_add_logFC <- function(pathway_df = pathway_df, 
                                 logFC_df, 
                                 orderBy = "nOverlap", cluster_meta = "Cell_type",
                                 calculateOnlyOverlapLogFC = FALSE, 
                                 compiled_pathways_df,
                                 onlySig=FALSE){
  if(calculateOnlyOverlapLogFC){
    avg_logFC = c()
    for(row in 1:nrow(pathway_df)){
      if(pathway_df$nOverlap[row]>0){
        ct_logFC_all = logFC_df[logFC_df[,cluster_meta]==pathway_df$Module[row],]
        genes = unlist(strsplit(pathway_df$Overlap[row], split = ","))
        for(gene in 1:length(genes)){
          #cat(genes[gene], "\n")
          if(sum(MT_gene_names==genes[gene])>0){
            #genes[gene] = names(MT_gene_names)[grep(genes[gene], MT_gene_names, fixed = TRUE)]
            tryCatch({genes[gene] = names(MT_gene_names)[grep(paste0("\\b",genes[gene],"\\b"), MT_gene_names)]}, 
                     warning=function(w) print(row))
            if(length(names(MT_gene_names)[grep(genes[gene], MT_gene_names, fixed = TRUE)]>1)){
              cat(row, "\n")
            }
          }
          else{
            next
          }
        }
        logFCs = ct_logFC_all$avg_logFC[match(x = genes, table = ct_logFC_all$GENE)]
        avg_logFC[row] = sum(logFCs[!is.na(logFCs)])
        if(length(logFCs)!=length(genes)){
          cat("Did not find equal number of logFC values in row ",row,"\n")
        }
      }
      else{
        avg_logFC[row] = 0
      }
    }
    pathway_df$avg_logFC = avg_logFC
  }
  else if(!calculateOnlyOverlapLogFC & onlySig){
    avg_logFC = c()
    #compiled_pathways_df = read.delim(compiled_pathways, header = TRUE, stringsAsFactors = FALSE)
    for(row in 1:nrow(pathway_df)){
      pathway_genes = compiled_pathways_df$GENE[compiled_pathways_df$MODULE==pathway_df$Pathway[row]]
      ct_logFC_all = logFC_df[logFC_df[,cluster_meta]==pathway_df$Module[row] & logFC_df$p_val_adj<0.05,] # take only those that are significant
      logFCs = ct_logFC_all$avg_logFC[match(x = pathway_genes, table = ct_logFC_all$GENE)]
      avg_logFC[row] = sum(logFCs[!is.na(logFCs)])
    }
    pathway_df$avg_logFC = avg_logFC
    return(pathway_df)
  }
  else{
    avg_logFC = c()
    #compiled_pathways_df = read.delim(compiled_pathways, header = TRUE, stringsAsFactors = FALSE)
    for(row in 1:nrow(pathway_df)){
      pathway_genes = compiled_pathways_df$GENE[compiled_pathways_df$MODULE==pathway_df$Pathway[row]]
      ct_logFC_all = logFC_df[logFC_df[,cluster_meta]==pathway_df$Module[row],]
      logFCs = ct_logFC_all$avg_logFC[match(x = pathway_genes, table = ct_logFC_all$GENE_human)]
      avg_logFC[row] = sum(logFCs[!is.na(logFCs)])
    }
    pathway_df$avg_logFC = avg_logFC
    return(pathway_df)
  }
}

# Apoptosis/Oxidative Phosphorylation/Immune system and maybe others will be duplicated (different databases, same pathway name)
makePathwayPretty <- function(vector){
  vector <- gsub("HALLMARK_","", vector)
  vector <- gsub("KEGG_","", vector)
  vector <- gsub("REACTOME_","", vector)
  vector <- gsub("BIOCARTA_","", vector)
  vector <- gsub("MHC_CLASS_II_ANTIGEN_PRESENTATION","MHC II antigens", vector)
  vector <- gsub("ADIPOCYTOKINE_SIGNALING_PATHWAY","Adipocytokine signaling", vector)
  vector <- gsub("XENOBIOTIC_METABOLISM","Xenobiotic metabolism", vector)
  vector <- gsub("INTERFERON_ALPHA_RESPONSE","Interferon alpha response", vector)
  vector <- gsub("MITOCHONDRIAL_PROTEIN_IMPORT","Mitochondrial protein import", vector)
  vector <- gsub("TRANS_GOLGI_NETWORK_VESICLE_BUDDING","Golgi vesicle budding", vector)
  vector <- gsub("AMYLOIDS","Amyloids", vector)
  vector <- gsub("COMPLEMENT","Complement", vector)
  vector <- gsub("EPITHELIAL_MESENCHYMAL_TRANSITION","Epithelial mesenchymal transition", vector)
  vector <- gsub("INFLAMMATORY_RESPONSE","Inflammatory Response", vector)
  vector <- gsub("IL6_JAK_STAT3_SIGNALING","IL6 JAK STAT3 Signaling", vector)
  vector <- gsub("P53HYPOXIA_PATHWAY","P53 Hypoxia pathway", vector)
  vector <- gsub("HYPOXIA","Hypoxia", vector)
  vector <- gsub("TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT","TCA cycle and ETC", vector)
  vector <- gsub("METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES","Amino acid metabolism", vector)
  vector <- gsub("GLYCOSPHINGOLIPID_BIOSYNTHESIS_GLOBO_SERIES","Glycosphingolipid biosynthesis", vector)
  vector <- gsub("GLUCONEOGENESIS","Gluconeogenesis", vector)
  vector <- gsub("NDKDYNAMIN_PATHWAY","NDK dyamin pathway", vector)
  vector <- gsub("LYSOSOME","Lysosome", vector)
  vector <- gsub("IRON_UPTAKE_AND_TRANSPORT","Iron uptake and transport", vector)
  vector <- gsub("APOPTOSIS","Apoptosis", vector)
  vector <- gsub("PROTEIN_FOLDING","Protein folding", vector)
  vector <- gsub("NITROGEN_METABOLISM","Nitrogen Metabolism", vector)
  vector <- gsub("REACTIVE_OXIGEN_SPECIES_PATHWAY","Reactive oxygen species", vector)
  vector <- gsub("ESTROGEN_RESPONSE_LATE","Late estrogen response", vector)
  vector <- gsub("INTERFERON_GAMMA_RESPONSE","Interferon gamma response", vector)
  vector <- gsub("ADIPOGENESIS","Adipogenesis", vector)
  vector <- gsub("TNFA_SIGNALING_VIA_NFKB","TNFa signaling via NFkB", vector)
  vector <- gsub("CHOLESTEROL_HOMEOSTASIS","Cholesterol homeostasis", vector)
  vector <- gsub("GLYCOLYSIS","Glycolysis", vector)
  vector <- gsub("AXON_GUIDANCE","Axon guidance", vector)
  vector <- gsub("ION_CHANNEL_TRANSPORT","Ion channel transport", vector)
  vector <- gsub("ESTROGEN_RESPONSE_EARLY","Early estrogen response", vector)
  vector <- gsub("PPAR_SIGNALING_PATHWAY","PPAR signaling", vector)
  vector <- gsub("STEROID_BIOSYNTHESIS","Steroid biosynthesis", vector)
  vector <- gsub("OXIDATIVE_PHOSPHORYLATION","Oxidative phosphorylation", vector)
  vector <- gsub("METABOLISM_OF_LIPIDS_AND_LIPOPROTEINS","Lipid and lipoprotein metabolism", vector)
  vector <- gsub("TRANSMEMBRANE_TRANSPORT_OF_SMALL_MOLECULES","Small molecule transport", vector)
  vector <- gsub("NEURONAL_SYSTEM","Neuronal system", vector)
  vector <- gsub("HEMOSTASIS","Hemostasis", vector)
  vector <- gsub("TIGHT_JUNCTION","Tight junction", vector)
  vector <- gsub("REGULATION_OF_ORNITHINE_DECARBOXYLASE_ODC","Regulation of ODC", vector)
  vector <- gsub("PROTEASOME","Proteasome", vector)
  vector <- gsub("MRNA_SPLICING","mRNA splicing", vector)
  vector <- gsub("SPLICEOSOME","Spliceosome", vector)
  vector <- gsub("SIGNALING_BY_WNT","Wnt signaling", vector)
  vector <- gsub("IMMUNE_SYSTEM","Immune System", vector)
  vector <- gsub("PGC1A_PATHWAY","PGC-1a pathway", vector)
  vector <- gsub("BRANCHED_CHAIN_AMINO_ACID_CATABOLISM","Branched chain amino acid catabolism", vector)
  vector <- gsub("ZINC_TRANSPORTERS","Zinc transporters", vector)
  vector <- gsub("ABC_TRANSPORTERS","ABC transporters", vector)
  vector <- gsub("GLYCOSPHINGOLIPID_METABOLISM","Glycosphingolipid metabolism", vector)
  vector <- gsub("FORMATION_OF_TUBULIN_FOLDING_INTERMEDIATES_BY_CCT_TRIC","Formation of tubulin folding intermediates", vector)
  vector <- gsub("MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION","Mitochondrial fatty acid beta oxidation", vector)
  vector <- gsub("CHOLESTEROL_BIOSYNTHESIS","Cholesterol biosynthesis", vector)
  vector <- gsub("UNFOLDED_PROTEIN_RESPONSE","Unfolded Protein Response", vector)
  vector <- gsub("PEROXISOME","Peroxisome", vector)
  vector <- gsub("MEMBRANE_TRAFFICKING","Membrane Trafficking", vector)
  vector <- gsub("TRANSCRIPTION","Transcription", vector)
  vector <- gsub("SYNTHESIS_SECRETION_AND_DEACYLATION_OF_GHRELIN","Ghrelin secretion", vector)
  vector <- gsub("DNA_REPAIR","DNA repair", vector)
  vector <- gsub("GAP_JUNCTION","Gap junction", vector)
  vector <- gsub("CELL_CYCLE","Cell cycle", vector)
  vector <- gsub("INSULIN_SIGNALING_PATHWAY","Insulin signaling", vector)
  vector <- gsub("SMOOTH_MUSCLE_CONTRACTION","Smooth Muscle Contraction", vector)
  vector <- gsub("CELL_CYCLE","Cell cycle", vector)
  vector <- gsub("CHEMOKINE_SIGNALING_PATHWAY","Chemokine signaling", vector)
  vector <- gsub("Complement_AND_COAGULATION_CASCADES","Complement and coagulation", vector)
  vector <- gsub("ANTIGEN_PROCESSING_AND_PRESENTATION","Antigen processing", vector)
  vector <- gsub("MET_PATHWAY","Met pathway", vector)
  vector <- gsub("PATHOGENIC_ESCHERICHIA_COLI_INFECTION","Infection", vector)
  vector <- gsub("APICAL_JUNCTION","Apical Junction", vector)
  vector <- gsub("REGULATION_OF_ACTIN_CYTOSKELETON","Actin cytoskeleton regulation", vector)
  vector <- gsub("P53_PATHWAY","p53 Pathway", vector)
  vector <- gsub("MAPK_SIGNALING_PATHWAY","MAPK Signaling", vector)
  vector <- gsub("IL2_STAT5_SIGNALING","IL2 STAT5 Signaling", vector)
  vector <- gsub("TRANSLATION","Translation", vector)
  vector <- gsub("PARKINSONS_DISEASE","Parkinson's disease", vector)
  vector <- gsub("ALZHEIMERS_DISEASE","Alzheimer's disease", vector)
  vector <- gsub("ANDROGEN_RESPONSE","Androgen response", vector)
  vector <- gsub("FATTY_ACID_TRIACYLGLYCEROL_AND_KETONE_BODY_METABOLISM","Fatty acid metabolism", vector)
  vector <- gsub("MTORC1_SIGNALING","mTORC1 signaling", vector)
  vector <- gsub("KRAS_SIGNALING_UP","K-Ras signaling", vector)
  vector <- gsub("INSULIN_PATHWAY","Insulin pathway", vector)
  
  vector <- gsub("TRANSPORT_OF_GLUCOSE_AND_OTHER_SUGARS_BILE_SALTS_AND_ORGANIC_ACIDS_METAL_IONS_AND_AMINE_COMPOUNDS",
                 "Transport of glucose, bile salts, organic acids, metal ions, amine compounds", vector)
  
  return(vector)
}

# takes in melted data frame and converts to a matrix
# First variable (column) are the rows
# Second variable are the columns
# Third variable is the values
convertToMatrix <- function(df){
  DEG_mat <- matrix(nrow = length(unique(df[,1])),
                    ncol = length(unique(df[,2])))
  rownames(DEG_mat) <- unique(df[,1])
  colnames(DEG_mat) <- unique(df[,2])
  for(col in colnames(DEG_mat)){
    for(row in rownames(DEG_mat)){
      if(length(df[,3][df[,1]==row & df[,2]==col])==0){
        DEG_mat[row,col] = 0 
      }
      else DEG_mat[row,col] = df[,3][df[,1]==row & df[,2]==col]
    }
  }
  return(DEG_mat)
}

# @param DEG_df - can be either 'MODULE' 'GENE' file or DEG output from Seurat with genes in 'GENE' column
# @param MODULE_column - if 'Cell_type' is not the desired "module" to run the analysis, set the desired module column here
# @param resources_path - path to resources directory with txt files of resources in 'module' 'gene' format
# @param num_cores - number of cores for parallel computing (default 1, no parallel computing)
# @param heatmap - saves heatmap if TRUE
# @param return_nonconcat - concatenates all module enrichments if TRUE
# @param remove_dat - remove .dat intermediate file after use if TRUE
# outputs the .txt reports for each "module"
# outputs an example heatmap of the top 50 consistent pathways
# returns a dataframe detailing each pathway enrichment result
# a note about logFC threshold:
#   i've seen opposite directions of pathway logFC when using 0.25 and 0.1 logFC thresholds, expectedly. something to consider when choosing a threshold
#   i.e. whether to be inclusive or only focus on DEGs with over 0.40546 logFC (this is a 1.5 fold change) (which you may claim are "valid" changes)
#   i'd prefer to be inclusive when discussing things on a pathway level
makePathwayEnrichmentDf <- function(DEG_df, 
                                    resources_path, 
                                    output_Dir="./PathwayEnrichmentResults", 
                                    convertToHuman=TRUE, 
                                    addlogFC=FALSE, 
                                    MODULE_column=NULL, 
                                    FDR_threshold=0.05, 
                                    logFC_threshold=0.1, 
                                    min_max=NULL,
                                    num_cores=1,
                                    heatmap=FALSE,
                                    return_nonconcat=TRUE,
                                    remove_dat=TRUE){
  
  # trim DEG output if necessary 
  if(length(colnames(DEG_df))>2){
    if(addlogFC){
      DEG_df_logFC_info = DEG_df # this is to add the logFC for the entire pathway
    }
    DEG_df = DEG_df[abs(DEG_df$avg_logFC)>logFC_threshold &
                      DEG_df$p_val_adj<FDR_threshold,]
    if(!is.null(MODULE_column)){
      DEG_df <- DEG_df[,c(MODULE_column,"GENE")]
    }
    else{
      DEG_df <- DEG_df[,c("Cell_type","GENE")]
    }
    colnames(DEG_df) <- c("MODULE","GENE")
  }else if(length(colnames(DEG_df))==2){
    colnames(DEG_df) <- c("MODULE","GENE")
  }
  else{
    cat("At least two columns are required, for example the geneset name in the 'MODULE' column and genes
        in the 'GENE' column.\n")
  }
  
  if(convertToHuman & addlogFC){
    temp <- convertDfGeneColumnMouseHuman(df = DEG_df_logFC_info, toSpecies = "human", forPathway = TRUE)
    DEG_df_logFC_info$HUMAN <- temp$GENE
    rm(temp)
    reference <- DEG_df_logFC_info[!duplicated(DEG_df_logFC_info[,c("GENE","HUMAN")]),]
    DEG_df$GENE <- DEG_df_logFC_info$HUMAN[match(DEG_df$GENE, DEG_df_logFC_info$GENE)]
    DEG_df_logFC_info$GENE <- DEG_df_logFC_info$HUMAN
  }
  if(convertToHuman & !addlogFC){
    DEG_df <- convertDfGeneColumnMouseHuman(df = DEG_df, toSpecies = "human", forPathway = TRUE)
  }
  
  # enrichment doesn't work if you have " " in the module name
  DEG_df$MODULE <- gsub(" ","_", DEG_df$MODULE)
  
  # run pathway enrichment ------------------------------------------------------------------------------------------
  # register local cluster of parallel computing cores
  registerDoParallel(num_cores)
  list_of_pathway_databases = list.files(path = resources_path, pattern = "*.txt", full.names = TRUE)
  ifelse(!dir.exists(file.path(output_Dir)), dir.create(file.path(output_Dir),recursive = T), FALSE)
  result = list()
  all_result = list()
  annotations = c()
  foreach(module=unique(DEG_df$MODULE), .combine=c) %dopar% {
    if (!file.exists(paste0(output_Dir,"/",module, ".txt"))){
      cat(module, "\n")
      deg_list = DEG_df[DEG_df$MODULE==module,] # why is this creating NAs?
      deg_list = deg_list[!is.na(deg_list$MODULE),]
      deg_list[["gene"]] <- deg_list$GENE
      deg_list[["module"]] <- deg_list$MODULE
      deg_list <- deg_list[,c("module","gene")]
      unique_module1 <- unique(deg_list$module)
      module1_len <- length(unique(deg_list$module))
      x <- list()
      for(z in 1:length(list_of_pathway_databases)){
        pathway_database <- list_of_pathway_databases[z]
        database_name = unlist(strsplit(pathway_database,"/"))
        database_name = database_name[length(database_name)]
        database_name = unlist(strsplit(database_name, ".txt"))[1]
        print(database_name)
        
        Module2 <- tool.read(pathway_database)
        Unique_module2 <- unique(Module2$module)
        
        Module2_len <- length(unique(Module2$module))
        
        data_matrix_for_enrichment <- data.frame()
        List_initial <- 1
        # go through the different databases
        for(k in 1:Module2_len){
          data_matrix_for_enrichment[List_initial,1] <- unique_module1
          data_matrix_for_enrichment[List_initial,2] <- length(deg_list$gene[which(match(deg_list$module, unique_module1)>0)])
          
          Overlapped_genes <- intersect(deg_list$gene[which(match(deg_list$module, unique_module1)>0)],
                                        Module2$gene[which(match(Module2$module, Unique_module2[k])>0)])
          data_matrix_for_enrichment[List_initial,3] <- length(Overlapped_genes)
          data_matrix_for_enrichment[List_initial,4] <- length(Module2$gene[which(match(Module2$module, Unique_module2[k])>0)])
          
          data_matrix_for_enrichment[List_initial, 5] <- 20000
          data_matrix_for_enrichment[List_initial,6] <- Unique_module2[k]
          
          if(length(Overlapped_genes)){
            data_matrix_for_enrichment[List_initial,7] <- paste(Overlapped_genes, collapse = ",")
          } else{
            data_matrix_for_enrichment[List_initial,7] <- c("NULL")
          }
          data_matrix_for_enrichment[List_initial,8] <- database_name
          
          List_initial=List_initial + 1
        }
        
        write.table(data_matrix_for_enrichment, file = paste0(output_Dir,"/",module,".dat"), quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = FALSE)
        # record_mat <- read.table(paste0("Csvs/pathway_enrichment/temp1.",module,".dat"))
        record_mat <- read.table(paste0(output_Dir,"/",module,".dat"))
        record_length <- dim(record_mat)
        enrichment_score <- data.frame()
        
        for(i in 1:record_length[1]){
          enrichment_score[i,1] <- record_mat[i,1]
          enrichment_score[i,2] <- phyper(record_mat[i,3], record_mat[i,4], record_mat[i,5]-record_mat[i,4], record_mat[i,2], lower.tail = FALSE)
          enrichment_score[i,3] <- record_mat[i,3]/record_mat[i,2]*record_mat[i,5]/record_mat[i,4]
          enrichment_score[i,4] <- 0
          enrichment_score[i,5]<-record_mat[i,8]
          enrichment_score[i,6]<-record_mat[i,6]
          enrichment_score[i,7]<-record_mat[i,3]
          enrichment_score[i,8]<-record_mat[i,7]
        }
        enrichment_score[,4]<-p.adjust(enrichment_score[,2], 'BH')
        colnames(enrichment_score) <- c("Module","Pval","Enrichment","FDR","PathwaySource","Pathway","nOverlap","Overlap")
        enrichment_score$ModuleGeneCount = length(deg_list$gene)
        enrichment_score <- enrichment_score[order(enrichment_score$FDR),]
        x[[database_name]] <- data.frame(enrichment_score)
        
        if(z==1){
          all_pathways_df <- enrichment_score
        }else{
          all_pathways_df <- rbind(all_pathways_df, enrichment_score)
        }
  
      }
      
      total = data.frame(stringsAsFactors = FALSE)
      for(iter in 1:length(x)){
        total = rbind(total, x[[iter]])
      }
      total = total[order(total$nOverlap, decreasing = TRUE),]
      write.table(total, paste0(output_Dir,"/",module, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
      if (remove_dat){
        file.remove(paste0(output_Dir,"/",module, ".dat"))
      }
    }else{
      cat(paste0(output_Dir,"/",module, ".txt already exists ... skipping \n"))
    }
  }

    # -------------------- for loop ends here ----------------------------------------------------------
  stopImplicitCluster() # clean up cluster
  if (return_nonconcat){
    return()
  }
  pathway_files = list.files(output_Dir)[grep(".txt", list.files(output_Dir))]
  pathway_df = data.frame(stringsAsFactors = FALSE)
  for(file in pathway_files){
    temp = read.delim(paste0(output_Dir,"/",file), stringsAsFactors = FALSE, header = TRUE)
    pathway_df = rbind(pathway_df, temp)
  }
  pathway_df = pathway_df[order(pathway_df$nOverlap, decreasing = TRUE),]
  pathway_df$`negLog10FDR` = -log10(pathway_df$FDR)
  
  # add summed logFC ------------------------------------------------------------------------------------------
  # can be very slow
  if(addlogFC){
    list_of_pathway_databases = list.files(path = resources_path, pattern = "*.txt", full.names = TRUE)
    compiled_pathways_df = data.frame(stringsAsFactors = FALSE)
    for(file in list_of_pathway_databases){
      compiled_pathways_df <- rbind(compiled_pathways_df, read.delim(file, stringsAsFactors = FALSE))
    }
    colnames(compiled_pathways_df) <- c("MODULE","GENE")
    
    if(!is.null(MODULE_column)){
      pathway_df <- pathway_df_add_logFC(pathway_df = pathway_df, 
                                         logFC_df = DEG_df_logFC_info, 
                                         cluster_meta = MODULE_column, 
                                         calculateOnlyOverlapLogFC = FALSE, onlySig = TRUE, 
                                         compiled_pathways_df = compiled_pathways_df)
    }
    else{
      pathway_df <- pathway_df_add_logFC(pathway_df = pathway_df, 
                                         logFC_df = DEG_df_logFC_info, 
                                         cluster_meta = "Cell_type",
                                         calculateOnlyOverlapLogFC = FALSE, onlySig = TRUE, 
                                         compiled_pathways_df = compiled_pathways_df)
    }
  }
  
  # make simple heatmap --------------------------------------------------------------------------------------------------
  # pick top consistent pathways
  pathway_df$numAppearances <- sapply(pathway_df$Pathway, function(x){return(nrow(pathway_df[pathway_df$Pathway==x,]))})
  pathway_df = pathway_df[order(pathway_df$numAppearances, decreasing = TRUE),]
  if(heatmap){
    topPathways <- unique(pathway_df$Pathway)[1:50]
    top_pathway_df = pathway_df[pathway_df$Pathway %in% topPathways,]
    
    # change pathway names to look better
    top_pathway_df$Pathway <- makePathwayPretty(top_pathway_df$Pathway)
    top_pathway_df <- top_pathway_df[!duplicated(top_pathway_df[,c("Module","Pathway")]),] # might be losing a "better" pathway here
    
    if(addlogFC){
      # hclust pathways by logFC
      pathway_heat_mat <- convertToMatrix(top_pathway_df[,c("Pathway","Module","avg_logFC")])
      pathway.dendro <- as.dendrogram(hclust(d = dist(x = pathway_heat_mat)))
      # get order
      pathway_order <- order.dendrogram(pathway.dendro)
      pathway_levels <- rownames(pathway_heat_mat)[pathway_order]
      
      # hclust pathways by Module
      if(length(unique(top_pathway_df$Module))>=2){
        pathway_heat_mat <- convertToMatrix(top_pathway_df[,c("Module","Pathway","nOverlap")])
        pathway.dendro <- as.dendrogram(hclust(d = dist(x = pathway_heat_mat)))
        ct_order <- order.dendrogram(pathway.dendro)
        ct_levels <- rownames(pathway_heat_mat)[ct_order]
        top_pathway_df$Module <- factor(top_pathway_df$Module, levels = ct_levels)
      }
    } else{
      # hclust pathways by nOverlap
      pathway_heat_mat <- convertToMatrix(top_pathway_df[,c("Pathway","Module","nOverlap")]) # can change to Enrichment or FDR
      pathway.dendro <- as.dendrogram(hclust(d = dist(x = pathway_heat_mat)))
      # get order
      pathway_order <- order.dendrogram(pathway.dendro)
      pathway_levels <- rownames(pathway_heat_mat)[pathway_order]
      
      # hclust pathways by Module
      if(length(unique(top_pathway_df$Module))>=2){
        pathway_heat_mat <- convertToMatrix(top_pathway_df[,c("Module","Pathway","nOverlap")])
        pathway.dendro <- as.dendrogram(hclust(d = dist(x = pathway_heat_mat)))
        ct_order <- order.dendrogram(pathway.dendro)
        ct_levels <- rownames(pathway_heat_mat)[ct_order]
        top_pathway_df$Module <- factor(top_pathway_df$Module, levels = ct_levels)
      }
    }
    
    top_pathway_df$Pathway <- factor(top_pathway_df$Pathway, levels = pathway_levels)
    
    # plot heat map
    if(addlogFC){
      if(!is.null(min_max)){
        top_pathway_df$avg_logFC <- MinMax(top_pathway_df$avg_logFC, min = min_max[1], max = min_max[2]) # need Seurat
      }
      heat <- ggplot(top_pathway_df, aes(x=Module, y=Pathway, fill=avg_logFC)) +#, label=nOverlap, )) +
        #geom_tile(aes(colour="gray"), size=.7) +
        geom_tile(colour="grey", size=.7) +
        #geom_text(size=2.7, aes(colour=FDR<0.05)) +
        theme_bw() +
        #facet_grid(.~Tissue, scales = "free", space = "free") +
        #scale_fill_gradient2(low="blue",high="red", mid = "white", breaks = c(0), labels=c("")) + 
        scale_fill_gradient2(low="blue",high="red", mid = "white") + 
        theme(#strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          # strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
          #                           hjust = .5, size=14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 11.5),
          axis.text.y = element_text(size = 11),
          panel.grid = element_blank(), axis.line = element_blank(), 
          legend.title = element_text(margin=margin(t=0,b=3,r=0,l=0)))+
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        #scale_color_manual(values = c("black", "darkgray"))+
        labs(fill="logFC")
      pdf(paste0(output_Dir,"/heatmap.pdf"), height = 15, width = 20)
      print(heat)
      dev.off()
    } else{
      if(!is.null(min_max)){
        top_pathway_df$nOverlap <- MinMax(top_pathway_df$nOverlap, min = min_max[1], max = min_max[2]) # need Seurat
      }
      heat <- ggplot(top_pathway_df, aes(x=Module, y=Pathway, fill=nOverlap)) +#, label=nOverlap, )) +
        #geom_tile(aes(colour="gray"), size=.7) +
        geom_tile(colour="grey", size=.7) +
        #geom_text(size=2.7, aes(colour=FDR<0.05)) +
        theme_bw() +
        #facet_grid(.~Tissue, scales = "free", space = "free") +
        #scale_fill_gradient2(low="blue",high="red", mid = "white", breaks = c(0), labels=c("")) + 
        scale_fill_gradient2(low="blue",high="red", mid = "white") + 
        theme(#strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          # strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
          #                           hjust = .5, size=14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 11.5),
          axis.text.y = element_text(size = 11),
          panel.grid = element_blank(), axis.line = element_blank(), 
          legend.title = element_text(margin=margin(t=0,b=3,r=0,l=0)))+
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        #scale_color_manual(values = c("black", "darkgray"))+
        labs(fill="nOverlap")
      pdf(paste0(output_Dir,"/heatmap.pdf"), height = 15, width = 20)
      print(heat)
      dev.off()
    }
  }
  
  return(pathway_df)
  
}

