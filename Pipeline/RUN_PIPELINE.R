dnm_net_enhance <- function(disease) {
  require(data.table)
  require(readr)
  
  # pre-filter general human network (400, 700, top 50%, no filter, etc.)
  # source("Script/Pipeline/ppi_generation.R")
  # ppi <- ppi_generation(version = "11.5", min_score = 700)
  ppi <- read_rds("Data/string_ppi_700_weighted.rds")
  
  # retrieve disease-specific data from H2GKB
  source("Script/Pipeline/dsH2GKB_generation.R")
  dsH2GKB <- dsH2GKB_generation(disease = "CHD", filter_pct = .85)
  
  # map genes to ENSEMBL to increase overlap?
  
  # create a gene similarity matrix (and amplify it?)
  source("Script/Pipeline/similarity_matrix_generation.R")
  sim_mat <- similarity_matrix_generation(score_mat = dsH2GKB, type = "cosine")
  
  # enhance the general human network through random walk with the gene similarity matrix as an adjacency matrix
  source("Script/Pipeline/RWR_standard.R")
  net_enhanced <- RWR_standard(ppi, sim_mat, alpha = 0.75, difference = 1e-6)
  write_rds(net_enhanced, "Data/net_enhanced_CHD.rds")
  
  # denoise the network (NE)
  source("Script/Pipeline/denoise_NE.R")
  net_enhanced_dn <- denoise_NE(net_enhanced)
  colnames(net_enhanced_dn) <- colnames(net_enhanced)
  rownames(net_enhanced_dn) <- rownames(net_enhanced)
  write_rds(net_enhanced_dn, "Data/net_enhanced_dn_NE.rds")
  
  # Trim network to only genes with mut > 0
  gene_dnm <- fread("Data/Multi_NGS_Allgenes_0827_v3.csv")
  gene_dnm <- gene_dnm[mut > 0]
  gene_dnm <- gene_dnm$Gene
  gene_dnm <- intersect(gene_dnm, colnames(net_enhanced_dn))
  
  net_enhanced_dn <- net_enhanced_dn[gene_dnm, gene_dnm]
  write_rds(net_enhanced_dn, "Data/net_enhanced_dn_NE_filtered.rds")
  
  # Conduct another random walk
  ## input is known genes evenly split to total a probability of 1
  ## adjacency matrix is denoised enhanced network
  gene_known <- fread("Data/KnownCHD_n254.txt", header = FALSE, col.names = "Gene")
  gene_known <- gene_known$Gene
  
  source("Script/Pipeline/RWR_vector.R")
  results <- RWR_vector(gene_known, net_enhanced_dn, alpha = 0.25, difference = 1e-6)
  
  results_ranked <- results[order(-results$score),]
  
  return(results_ranked)
}
