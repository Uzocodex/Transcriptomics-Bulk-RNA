####Introduction####

# title: "Batch Effect Removal"
# author: "Uzochukwu Gospel Ukachukwu"
# date: "2026-03-15"
# R_version: "4.4.3"
# output: R_script 

#### End ####

#### Package setup ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MBECS")
install.packages("devtools")
install.packages("doParallel")
devtools::install_github("wdl2459/ConQuR")
install_github("cran/PERMANOVA")

# Load necessary libraries
library(phyloseq)
library(MBECS)    # For MBECS batch correction
library(doParallel)# ConQuR depends on it
library(ConQuR)   # For ConQuR batch correction
library(foreach)

#1. Extract OTU count matrix from phyloseq object (samples x taxa)
otu_mat <- as(otu_table(physeq), "matrix")
          if(taxa_are_rows(physeq)) {
          otu_mat <- t(otu_mat)
          }

#2. Prepare metadata - must include batch and biological group columns
# For example, batch in "batch" column, biological group in "group" column
meta <- data.frame(sample_data(physeq))

####Batch Correction with MBECS####
# Initialize MBEC object from phyloseq
mbec_obj <- mbecProcessInput(physeq, required.col = c("Batch", "Group"))

# Run batch correction with default method 'bat' (ComBat on CLR-transformed data)
mbec_obj <- mbecRunCorrections(mbec_obj,
                               model.vars = c("Batch", "Group"),
                               method = "bat",
                               type = "clr")

# Retrieve corrected phyloseq object (CLR corrected counts)
physeq_mbec_corrected <- mbecGetPhyloseq(mbec_obj, type = "clr")

####End####

####Batch Correction with ConQuR####
# Run ConQuR correction (default) - specify batch and condition columns accordingly

otu_mat1 <- t(otu_mat)
meta1 <- as.matrix(meta)

batchid <- as.factor(meta1[,"Batch"])
summary(batchid)
batchid1 <- relevel(batchid, ref = "Batch_1")
batchid2 <- relevel(batchid, ref = "Batch_2")
batchid3 <- relevel(batchid, ref = "Batch_3")


covariate <- as.factor(meta1[, "Group"])
covariate

conqr_corrected_mat3 <- ConQuR(  tax_tab = otu_mat1,
                                batchid = batchid3,
                                covariates = covariate,
                                batch_ref = "Batch_3",
                                num_core = 2)
      

conqr_corrected_mat1 
conqr_corrected_mat2
conqr_corrected_mat3

# Retrieve corrected phyloseq object 

otu_conqr_corrected1_1 <- otu_table(conqr_corrected_mat1, taxa_are_rows = FALSE)
physeq_conqr_corrected1_1 <- physeq
otu_table(physeq_conqr_corrected1_1) <- otu_conqr_corrected1_1
otu_conqr_corrected_df1_1 <- as.data.frame(t(otu_conqr_corrected1_1))


otu_conqr_corrected2_1 <- otu_table(conqr_corrected_mat2, taxa_are_rows = FALSE)
physeq_conqr_corrected2_1 <- physeq
otu_table(physeq_conqr_corrected2_1) <- otu_conqr_corrected2_1
otu_conqr_corrected_df2_1 <- as.data.frame(t(otu_conqr_corrected2_1))

otu_conqr_corrected3_1 <- otu_table(conqr_corrected_mat3, taxa_are_rows = FALSE)
physeq_conqr_corrected3_1 <- physeq
otu_table(physeq_conqr_corrected3_1) <- otu_conqr_corrected3_1
otu_conqr_corrected_df3_1 <- as.data.frame(t(otu_conqr_corrected3_1))

#### End ####

## Run ConQuR correction (Penalized) - specify batch and condition columns accordingly

conqrait_corr_mat2 <- ConQuR(  tax_tab = otu_mat1,
                               batchid = batchid2,
                               covariates = covariate,
                               batch_ref = "Batch_2",
                               logistic_lasso = TRUE,
                               quantile_type = "lasso",
                               lambda_quantile = "2p/logn",
                               interplt = TRUE,
                               taus = seq(0.005, 0.995, by = 0.005), # Default high resolution is usually good
                               num_core = 4) # Use multiple cores to speed up computation

conqrait_corr_mat1
conqrait_corr_mat2
conqrait_corr_mat3


otu_conqr_corrected1_2 <- otu_table(conqrait_corr_mat1, taxa_are_rows = FALSE)
physeq_conqr_corrected1_2 <- physeq
otu_table(physeq_conqr_corrected1_2 ) <- otu_conqr_corrected1_2
otu_conqr_corrected_df1_2 <- as.data.frame(t(otu_conqr_corrected1_2))

otu_conqr_corrected2_2 <- otu_table(conqrait_corr_mat2, taxa_are_rows = FALSE)
physeq_conqr_corrected2_2 <- physeq
otu_table(physeq_conqr_corrected2_2 ) <- otu_conqr_corrected2_2
otu_conqr_corrected_df2_2 <- as.data.frame(t(otu_conqr_corrected2_2))

otu_conqr_corrected3_2 <- otu_table(conqrait_corr_mat3, taxa_are_rows = FALSE)
physeq_conqr_corrected3_2 <- physeq
otu_table(physeq_conqr_corrected3_2 ) <- otu_conqr_corrected3_2
otu_conqr_corrected_df3_2 <- as.data.frame(t(otu_conqr_corrected3_2))


#### Evaluate if ConQuR worked with Bray-Curtis ####

par(mfrow=c(3, 3))  #display multiple plots on a single page

Plot_PCoA(TAX=otu_mat1, factor=batchid1, main="Before Correction, Bray-Curtis")
Plot_PCoA(TAX=otu_conqr_corrected1_1, factor=batchid1, main="ConQuR (Default), Bray-Curtis")
Plot_PCoA(TAX=otu_conqr_corrected1_2, factor=batchid1, main="ConQuR (Penalized), Bray-Curtis")

Plot_PCoA(TAX=otu_mat1, factor=batchid2, main="Before Correction, Bray-Curtis")
Plot_PCoA(TAX=otu_conqr_corrected2_1, factor=batchid2, main="ConQuR (Default), Bray-Curtis")
Plot_PCoA(TAX=otu_conqr_corrected2_2, factor=batchid2, main="ConQuR (Penalized), Bray-Curtis")

Plot_PCoA(TAX=otu_mat1, factor=batchid3, main="Before Correction, Bray-Curtis")
Plot_PCoA(TAX=otu_conqr_corrected3_1, factor=batchid3, main="ConQuR (Default), Bray-Curtis")
Plot_PCoA(TAX=otu_conqr_corrected3_2, factor=batchid3, main="ConQuR (Penalized), Bray-Curtis")

#### Evaluate if ConQuR worked with Aitchison ####

par(mfrow=c(3, 3))  #display multiple plots on a single page

Plot_PCoA(TAX=otu_mat1, factor=batchid1, dissimilarity="Aitch", main="Before Correction, Aitchison")
Plot_PCoA(TAX=otu_conqr_corrected1_1, factor=batchid1, dissimilarity="Aitch", main="ConQuR (Default), Aitchison")
Plot_PCoA(TAX=otu_conqr_corrected1_2, factor=batchid1, dissimilarity="Aitch", main="ConQuR (Penalized), Aitchison")

Plot_PCoA(TAX=otu_mat1, factor=batchid2, dissimilarity="Aitch", main="Before Correction, Aitchison")
Plot_PCoA(TAX=otu_conqr_corrected2_1, factor=batchid2, dissimilarity="Aitch", main="ConQuR (Default), Aitchison")
Plot_PCoA(TAX=otu_conqr_corrected2_2, factor=batchid2, dissimilarity="Aitch", main="ConQuR (Penalized), Aitchison")

Plot_PCoA(TAX=otu_mat1, factor=batchid3, main="Before Correction, Aitchison")
Plot_PCoA(TAX=otu_conqr_corrected3_1, factor=batchid3,dissimilarity="Aitch", main="ConQuR (Default), Aitchison")
Plot_PCoA(TAX=otu_conqr_corrected3_2, factor=batchid3,dissimilarity="Aitch", main="ConQuR (Penalized), Aitchison")



#### PERMANOVA EVALUATION ####
#### Edited from ConQuR::PERMANOVA_R2 by Perplexity ####
PERMANOVA_R2 <- function(TAX, batchid, covariates, key_index) {
  # Construct metadata with batchid and covariates
  metadata <- data.frame(batchid = batchid)
  if (!missing(covariates)) {
    metadata <- cbind(metadata, covariates)
  }
  
  tab_count <- tab_rel <- matrix(ncol = 3, nrow = 2)
  colnames(tab_count) <- colnames(tab_rel) <- c("standard", "sqrt.dist=T", "add=T")
  rownames(tab_count) <- rownames(tab_rel) <- c("batch", "key")
  
  # Create dynamic formulas
  formula_batch <- reformulate("batchid", response = "TAX")
  key_var_name <- colnames(metadata)[key_index]  # +1 due to batchid column at 1
  formula_key <- reformulate(key_var_name, response = "TAX")
  
  # Standard Bray-Curtis distance PERMANOVA
  res <- adonis2(formula_batch, data = metadata, permutations = 999)
  tab_count[1, 1] <- res$R2[1]
  res <- adonis2(formula_batch, data = metadata, permutations = 999, sqrt.dist = TRUE)
  tab_count[1, 2] <- res$R2[1]
  res <- adonis2(formula_batch, data = metadata, permutations = 999, add = TRUE)
  tab_count[1, 3] <- res$R2[1]
  
  res <- adonis2(formula_key, data = metadata, permutations = 999)
  tab_count[2, 1] <- res$R2[1]
  res <- adonis2(formula_key, data = metadata, permutations = 999, sqrt.dist = TRUE)
  tab_count[2, 2] <- res$R2[1]
  res <- adonis2(formula_key, data = metadata, permutations = 999, add = TRUE)
  tab_count[2, 3] <- res$R2[1]
  
  # Aitchison distance (Euclidean on clr transformed)
  Z <- coda.base::dist(TAX + 0.5, method = "aitchison")
  
  res <- adonis2(Z ~ batchid, data = metadata, permutations = 999, method = "euclidean")
  tab_rel[1, 1] <- res$R2[1]
  res <- adonis2(Z ~ batchid, data = metadata, permutations = 999, method = "euclidean", sqrt.dist = TRUE)
  tab_rel[1, 2] <- res$R2[1]
  res <- adonis2(Z ~ batchid, data = metadata, permutations = 999, method = "euclidean", add = TRUE)
  tab_rel[1, 3] <- res$R2[1]
  
  # For 'key' variable, create formula dynamically (no get())
  formula_key_2 <- reformulate(key_var_name, response = "Z")
  
  res <- adonis2(formula_key_2, data = metadata, permutations = 999, method = "euclidean")
  tab_rel[2, 1] <- res$R2[1]
  res <- adonis2(formula_key_2, data = metadata, permutations = 999, method = "euclidean", sqrt.dist = TRUE)
  tab_rel[2, 2] <- res$R2[1]
  res <- adonis2(formula_key_2, data = metadata, permutations = 999, method = "euclidean", add = TRUE)
  tab_rel[2, 3] <- res$R2[1]
  
  return(list(tab_count = tab_count, tab_rel = tab_rel))
}

#### Apply function on corrected tables accordingly ####
library(dplyr)
Original1 <- PERMANOVA_R2(TAX=otu_mat1, batchid=batchid1, 
                          covariates=covar %>% dplyr::select(where(~ !any(is.na(.)))), key_index=1)

ConQuR_default1 <- PERMANOVA_R2(TAX=otu_conqr_corrected1_1, batchid=batchid1, 
                                covariates=covar %>% dplyr::select(where(~ !any(is.na(.)))), key_index=1)

ConQuR_Penalized1 <- PERMANOVA_R2(TAX=otu_conqr_corrected1_2, batchid=batchid1, 
                                  covariates=covar %>% dplyr::select(where(~ !any(is.na(.)))), key_index=1)
#--------------------------------------------------------------------------------------------------------

Original2 <- PERMANOVA_R2(TAX=otu_mat1, batchid=batchid2, 
                          covariates=covar %>% dplyr::select(where(~ !any(is.na(.)))), key_index=1)

ConQuR_default2 <- PERMANOVA_R2(TAX=otu_conqr_corrected2_1, batchid=batchid2, 
                                covariates=covar %>% dplyr::select(where(~ !any(is.na(.)))), key_index=1)

ConQuR_Penalized2 <- PERMANOVA_R2(TAX=otu_conqr_corrected2_2, batchid=batchid2, 
                                  covariates=covar %>% dplyr::select(where(~ !any(is.na(.)))), key_index=1)
#-----------------------------------------------------------------------------------------------------------

Original3 <- PERMANOVA_R2(TAX=otu_mat1, batchid=batchid3, 
             covariates=covar %>% dplyr::select(where(~ !any(is.na(.)))), key_index=1)

ConQuR_default3 <- PERMANOVA_R2(TAX=otu_conqr_corrected3_1, batchid=batchid3, 
             covariates=covar %>% dplyr::select(where(~ !any(is.na(.)))), key_index=1)

ConQuR_Penalized3 <- PERMANOVA_R2(TAX=otu_conqr_corrected3_2, batchid=batchid3, 
             covariates=covar %>% dplyr::select(where(~ !any(is.na(.)))), key_index=1)

#### End ####