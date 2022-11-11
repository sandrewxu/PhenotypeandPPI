# Initial Setup
# Libraries
library(MultiRNG) # Distributions
library(DiSNEP)

# Inputs
## STRING PPI network
STRING_full <- read_rds("Data/string_ppi_full_weighted.rds")

## H2GKB genes
# H2GKB_names <- read_rds("Data/H2GKB_names.rds")

## Mutability Table (EncoreDNM)
# mut_full <- readxl::read_xlsx("Data/elife-75551-supp1-v2.xlsx",
#                               sheet = "STable 20",
#                               range = "B2:B18456")

## Intersect genes
gene_all <- row.names(STRING_full)
# gene_all <- intersect(row.names(STRING_full), H2GKB_names)
# gene_all <- intersect(gene_all, mut_full$gene)


## Known Genes
gene_known <-
  unlist(fread("Data/KnownCHD_n254.txt", header = FALSE))

gene_known <- gene_all[which(gene_all %in% gene_known)]

simulation_KS <-
  function(num_known = 50,
           num_total = 2000,
           all_genes = gene_all,
           known_genes = gene_known,
           ppi_full = STRING_full,
           pheno_known = 200,
           pheno_ctrl = 100) {
    # Step 1: Generate Known Genes and Noise Genes
    true_known <- sample(known_genes, num_known)
    true_noise <-
      sample(setdiff(all_genes, known_genes), num_total - num_known)
    true_all <- c(true_known, true_noise)
    
    # Step 2: Simulate Phenotype Information using a Multivariate Beta Distribution
    HPO_sim <-
      draw.dirichlet(pheno_known + pheno_ctrl, num_total, rep(1, num_total), 10)
    HPO_sim <- HPO_sim * 8 / mean(HPO_sim)
    
    ## Associated
    high_asso_sim <-
      draw.dirichlet(pheno_known, num_known/2, rep(1, num_known/2), 10)
    high_asso_sim <- high_asso_sim * 10 / mean(high_asso_sim)
    
    med_asso_sim <-
      draw.dirichlet(pheno_known, num_known/2, rep(1, num_known/2), 10)
    med_asso_sim <- med_asso_sim * 9 / mean(med_asso_sim)
    
    HPO_sim[1:pheno_known, 1:(num_known/2)] <- high_asso_sim
    HPO_sim[1:pheno_known, (num_known/2+1):num_known] <- med_asso_sim
    
    # Step 3: Use the KS Test to generate initial p-values
    KS_pval <- apply(HPO_sim, 2, function (x) {
      ind <- pheno_known
      last <- length(x)
      res <- ks.test(x[1:ind], x[ind + 1:last])
      return(res$p.value)
    })
    
    KS_res <- data.frame(gene = true_all, pval = KS_pval)
    KS_res <- KS_res[order(KS_res$pval), ]
    
    # Step 4: Enhance using a PPI network (GeneWanderer)
    ppi <- ppi_full[true_all, true_all]
    # ppi_log <- log(ppi + 1 + 1e-16)
    # ppi_scale <- 999 / max(ppi_log)
    # ppi_log <- ppi_scale * ppi_log
    
    source("Script/Simulation/REDUCE.R")
    ppi_dn <- REDUCE(ppi, K = 3, alpha = 0.9)
    ppi_gen <- ppi_dn
    ppi_cut <- quantile(ppi_gen, 0.9)
    ppi_gen[ppi_gen < ppi_cut] <- 0
    ppi_gen[ppi_gen > 0] <- 1
    
    ppi_res <-
      DiSNEP::diffus_vec(KS_res,
                         ppi_gen,
                         type = "pvalue",
                         iter = 50,
                         top = 500)
    
    # Step 5: Enhance using the integrated network
    # gene_score_matrix <- HPO_sim[1:pheno_known,]
    # gene_similarity <-
    #   SMUT::eigenMapMatMult(t(gene_score_matrix), gene_score_matrix) / (sqrt(tcrossprod(colSums((gene_score_matrix) ^
    #                                                                                               2
    #   ))))
    gene_similarity <- matrix(rnorm(num_total*num_total, mean = 0.05, sd = 0.01), nrow = num_total, ncol = num_total)
    gene_similarity[1:num_known, 1:num_known] <- matrix(rnorm(num_known*num_known, mean = 0.1, sd = 0.02), nrow = num_known, ncol = num_known)
    gene_similarity[1:(num_known/2), 1:(num_known/2)] <- matrix(rnorm((num_known/2)*(num_known/2), mean = 0.25, sd = 0.05), nrow = num_known/2, ncol = num_known/2)
    
    colnames(gene_similarity) <- true_all
    rownames(gene_similarity) <- true_all
    
    source("Script/Simulation/INTEGRATE.R")
    net_integrated <- INTEGRATE(ppi, gene_similarity, r = 0.25)
    net_integrated <- REDUCE(net_integrated, K = 3, alpha = 0.9)
    novel_res <-
      DiSNEP::diffus_vec(
        KS_res,
        net_integrated,
        type = "pvalue",
        iter = 50,
        top = 500
      )
    
    res <- list(
      "KS" = KS_res,
      "PPI" = ppi_res,
      "Novel" = novel_res,
      "Known" = true_known
    )
  }


-
  # Run the Simulation
  N <- 5
allresults <- vector("list", N)

for (i in 1:N) {
  cat(paste("Simulation", i, "\n"))
  allresults[[i]] <- simulation_KS(num_known = 50, num_total = 2000)
}

write_rds(allresults, "Data/simulationresults102222_5.rds")
