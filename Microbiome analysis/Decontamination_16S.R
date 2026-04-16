####Introduction####

# title: "Decontamination Steps"
# author: "Uzochukwu Gospel Ukachukwu"
# date: "2026-03-15"
# R_version: "4.4.3"
# output: R_script 



# Install packages if not already installed
# install.packages("BiocManager")
# BiocManager::install("decontam")
# BiocManager::install("phyloseq")

library(dada2)
library(phyloseq)
library(decontam)

# Example: Assume you already have:
# 1. An ASV table (sequence table) from DADA2, called seqtab.nochim
# 2. Sample metadata including information about negative controls

# ---------------------
# 1. Create phyloseq object
# ---------------------

# Example: load your metadata (make sure sample names match the ASV table column names)
       # metadata.csv should include a column indicating if a sample is a 'negative control'
sample_16smetadata_df1 <- read.csv("sample_metadata/16S_metadata1.csv", row.names = 1) 
sample_16smetadata_df2 <- read.csv("sample_metadata/16S_metadata2.csv", row.names = 1)
sample_16smetadata_df3 <- read.csv("sample_metadata/16S_metadata3.csv", row.names = 1)

# Step 1: Extract sample names from metadata
meta_names1 <- rownames(sample_16smetadata_df1)
meta_names2 <- rownames(sample_16smetadata_df2)
meta_names3 <- rownames(sample_16smetadata_df3)

# Step 2: Check if sample names are in rows or columns of seqtab.nochim
# If samples are in rownames

if (all(meta_names1 %in% rownames(seqtab.nochim1))) {
  seqtab.nochim1 <- seqtab.nochim1[meta_names1, , drop = FALSE]
}

if (all(meta_names2 %in% rownames(seqtab.nochim2))) {
  seqtab.nochim2 <- seqtab.nochim2[meta_names2, , drop = FALSE]
}

if (all(meta_names3 %in% rownames(seqtab.nochim3))) {
  seqtab.nochim3 <- seqtab.nochim3[meta_names3, , drop = FALSE]
}

# Step 3: Verify if the order now matches
identical(rownames(seqtab.nochim1), rownames(sample_16smetadata_df1))

identical(rownames(seqtab.nochim2), rownames(sample_16smetadata_df2))

identical(rownames(seqtab.nochim3), rownames(sample_16smetadata_df3))


# Load taxonomy if already assigned (optional at this stage)
# Here we assume tax is your taxonomy table; if not yet done, skip
# tax <- assignTaxonomy(seqtab.nochim, "path/to/reference.fasta", multithread=TRUE)

# Create phyloseq object (without taxonomy for now if not available)
ps_obj1 <- phyloseq(otu_table(seqtab.nochim1, taxa_are_rows=FALSE),
               sample_data(sample_16smetadata_df1)
               # ,tax_table(tax) # Optional if taxonomy done
)

ps_obj2 <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows=FALSE),
                    sample_data(sample_16smetadata_df2)
                    # ,tax_table(tax) # Optional if taxonomy done
)

ps_obj3 <- phyloseq(otu_table(seqtab.nochim3, taxa_are_rows=FALSE),
                    sample_data(sample_16smetadata_df3)
                    # ,tax_table(tax) # Optional if taxonomy done
)


# ---------------------
# 2. Identify contaminants
# ---------------------

# Method 1: Prevalence-based method (requires negative controls info)
# Your metadata should have a column like "Sample_or_Control" with values "Sample" or "Control"
ps.decontam1 <- isContaminant(ps_obj1, method="prevalence", 
                             neg = sample_data(ps_obj1)$Group == "N-Control",
                             threshold=0.5)    # threshold can be adjusted depending on stringency

ps.decontam2 <- isContaminant(ps_obj2, method="prevalence", 
                              neg = sample_data(ps_obj2)$Group == "N-Control",
                              threshold=0.5)

ps.decontam3 <- isContaminant(ps_obj3, method="prevalence", 
                              neg = sample_data(ps_obj3)$Group == "N-Control",
                              threshold=0.5)


table(ps.decontam1$contaminant) # How many were flagged
table(ps.decontam2$contaminant)
table(ps.decontam3$contaminant)


# ---------------------
# 3. Remove contaminants
# ---------------------
contaminants1 <- rownames(subset(ps.decontam1, contaminant == TRUE))
ps_obj1.cleaned <- prune_taxa(!taxa_names(ps_obj1) %in% contaminants1, ps_obj1)

contaminants2 <- rownames(subset(ps.decontam2, contaminant == TRUE))
ps_obj2.cleaned <- prune_taxa(!taxa_names(ps_obj2) %in% contaminants2, ps_obj2)

contaminants3 <- rownames(subset(ps.decontam3, contaminant == TRUE))
ps_obj3.cleaned <- prune_taxa(!taxa_names(ps_obj3) %in% contaminants1, ps_obj3)


# ---------------------
# 4. Save cleaned ASV table
# ---------------------
dir.create("Results/ASV_Tables", recursive = T)

asv_table_clean1 <- as(otu_table(ps_obj1.cleaned), "matrix")
write.csv(asv_table_clean1, "Results/ASV_Tables/ASV_table1.csv")

asv_table_clean2 <- as(otu_table(ps_obj2.cleaned), "matrix")
write.csv(asv_table_clean2, "Results/ASV_Tables/ASV_table2.csv")

asv_table_clean3 <- as(otu_table(ps_obj3.cleaned), "matrix")
write.csv(asv_table_clean3, "Results/ASV_Tables/ASV_table3.csv")
