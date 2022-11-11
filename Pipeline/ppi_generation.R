ppi_generation <- function(version = "11.5", min_score = 700) {
  require(STRINGdb)
  require(igraph)
  require(biomaRt)
  require(readr)
  
  string_db <- STRINGdb$new(version = version, species=9606, score_threshold = min_score)
  human_graph <- string_db$get_graph()
  adj_matrix <- as_adjacency_matrix(human_graph, attr = "combined_score")
  
  adj_ids <- unlist(rownames(adj_matrix))
  dict_names <- string_db$get_proteins()
  dict_names <- dict_names[c("protein_external_id", "preferred_name")]
  colnames(dict_names) <- c("STRING_id", "alias")
  dict_names <- dict_names[dict_names$STRING_id%in%adj_ids,]
  newnames <- dict_names$alias
  rownames(adj_matrix) <- newnames
  colnames(adj_matrix) <- newnames
  
  ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
  nullrows <- Matrix::rowSums(ppi)==0
  ppi <- ppi[!nullrows,!nullrows]
  
  # write_rds(ppi, "Data/string_ppi_700_weighted.rds")
  
  ppi
}
