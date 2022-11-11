dsH2GKB_generation <- function(disease, filter_pct = .85) {
  require(readxl)
  require(data.table)
  require(stringr)
  require(readr)
  
  # Download H2GKB
  # url_H2GKB <- "https://github.com/WGLab/Phen2Gene/releases/download/1.1.0/H2GKBs.zip"
  # download.file(url_H2GKB, destfile = "Data/H2GKBs.zip")
  # unzip("Data/H2GKBs.zip", exdir = "Data/H2GKBs/")
  # unlink("Data/H2GKBs.zip")
  
  
  # Store H2GKB data
  dir <- "Data/H2GKBs/" # Directory to store H2GKBs files
  files_skewness <- list.files(paste0(dir,"skewness"), full.names = TRUE)
  
  # Filter H2GKB using disease-specific HPO ids (figure out a weighting?)
  HPO_disease <- read_xlsx("Data/Disease-HPO-ids.xlsx", sheet = disease, col_names = TRUE)
  HPO_disease_ids <- unique(HPO_disease$HPO_ID)
  HPO_disease_ids <- gsub(":","_", HPO_disease_ids)
  
  H2GKB_ids <- unlist(lapply(files_skewness, function(x) {
    str_sub(x, -10, -1)
  }))
  H2GKB_ids <- intersect(H2GKB_ids, HPO_disease_ids)
  
  files_knowledgebase <- unlist(lapply(H2GKB_ids, function(x) {
    paste0(dir, "Knowledgebase/", x, ".candidate_gene_list")
  }))
  files_skewness <- unlist(lapply(H2GKB_ids, function(x) {
    paste0(dir, "skewness/", x)
  }))
  
  HPO_data <- lapply(files_knowledgebase, fread, header = TRUE)
  skew_data <- lapply(files_skewness, fread, header = FALSE)
  skew_data <- unlist(skew_data)
  
  # Create a gene list -- unfiltered
  gene_H2GKB <- lapply(HPO_data, function(x) {
    x$Gene
  })
  gene_H2GKB <- unique(unlist(gene_H2GKB))
  # write_rds(gene_H2GKB, "Data/gene_H2GKB.rds")
  # gene_STRING <- read_rds("Data/gene_STRING.rds")
  # gene_both <- intersect(gene_H2GKB, gene_STRING)
  gene_both <- gene_H2GKB
  
  # Weight HPO scores using skewness (recommended by Phen2Gene)
  HPO_weighted <- lapply(1:length(HPO_data), function(x) {
    HPO_data[[x]]$Score <- HPO_data[[x]]$Score * skew_data[x]
    return(HPO_data[[x]])
  })
  
  ## Keep 'Gene' and 'Score' columns
  HPO_weighted_sub <- lapply(HPO_weighted, function(x) {
    x[,.(Gene, Score)]
  })
  
  ## Make sure each HPO ID has all genes
  gene_both <- as.data.table(gene_both)
  HPO_allgenes <- lapply(HPO_weighted_sub, function(x) {
    HPO_merge <- merge.data.table(x, gene_both, by.x = "Gene", by.y = "gene_both", all.y = TRUE)
    HPO_merge$Score[is.na(HPO_merge$Score)] <- 0
    return(HPO_merge)
  })
  
  ## Generate a gene score matrix (n HPO IDs x 20 967 Genes)
  HPO_scores <- lapply(HPO_allgenes, function(x) {
    x[, .(Score)]
  })
  
  gene_score_matrix <- matrix(unlist(HPO_scores), nrow = length(HPO_scores), byrow = TRUE)
  colnames(gene_score_matrix) <- HPO_allgenes[[1]]$Gene
  rownames(gene_score_matrix) <- H2GKB_ids
  
  # Keep top scores
  cutoff <- quantile(gene_score_matrix, filter_pct)
  gene_score_matrix[gene_score_matrix<cutoff] <- 0
  
  # Clear null columns
  nullcols <- colSums(gene_score_matrix)==0
  gene_score_matrix <- gene_score_matrix[, !nullcols]
  
  gene_score_matrix
}