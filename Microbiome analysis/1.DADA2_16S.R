####Introduction####

# title: "16S_rDNA metagenomics analysis - Stool DNA from Humans"
# author: "Uzochukwu Gospel Ukachukwu"
# date: "2026-03-15"
# R_version: "4.4.3"
# output: R_script 

####End ####

#### Package setup ####

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2")
#BiocManager::install("phyloseq")
#BiocManager::install("msa")
#BiocManager::install(c("microbiome","ComplexHeatmap"))
install.packages("remotes")
#remotes::install_github("jfq3/ggordiplots")
#remotes::install_github("kstagaman/phyloseqCompanion")
#remotes::install_github("igraph/rigraph")
remotes::install_github("david-barnett/microViz")

library(dada2); packageVersion("dada2")
library(vegan)
library(phyloseq)
library(decontam)
library(msa)
library(phangorn)
library(ggplot2)
library(gridExtra)
library(microViz)

#### End ####


#### Set Directory ####
setwd("C:/Users/Documents/Microbiome Analysis/16s_analysis")
getwd()
path1 <- "16s_batch1_trimmed_data/"   # Change to the directory containing the fastq files
path2 <- "16s_batch2_trimmed_data/" # Change to the directory containing the fastq files
path3 <- "16s_batch3_trimmed_data/" # Change to the directory containing the fastq files

list.files(path1,pattern = "\\.fastq?$", full.names = TRUE, recursive = TRUE)
list.files(path2,pattern = "\\.fastq?$", full.names = TRUE, recursive = TRUE)
list.files(path3,pattern = "\\.fastq?$", full.names = TRUE, recursive = TRUE)

#Default plotting theme for all subsequent ggplot2 plots in your current R session to the black-and-white "theme_bw" style
theme_set(theme_bw())

# Create directory to save results of the analysis and subdirectory for quality control

dir.create("Results/1.QC", recursive = T)

#### End ####


#### Quality Control and Filtering####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
#### For Batch1
fnFs1 <- sort(list.files(path1, pattern="_R1_001.fastq", full.names = TRUE, recursive = T))
fnFs1
fnRs1 <- sort(list.files(path1, pattern="_R2_001.fastq", full.names = TRUE,recursive = T))
fnRs1

# For Batch2
fnFs2 <- sort(list.files(path2, pattern="_R1.fastq", full.names = TRUE, recursive = T))
fnFs2
fnRs2 <- sort(list.files(path2, pattern="_R2.fastq", full.names = TRUE, recursive = T))
fnRs2

# For Batch3
fnFs3 <- sort(list.files(path3, pattern="_R1.fastq", full.names = TRUE, recursive = T))
fnFs3
fnRs3 <- sort(list.files(path3, pattern="_R2.fastq", full.names = TRUE, recursive = T))
fnRs3

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names1 <- sapply(strsplit(basename(fnFs1), "_"), `[`, 1)
sample.names1

sample.names1 <- sapply(strsplit(basename(fnRs1), "_"), `[`, 1)
sample.names1

sample.names2 <- sapply(strsplit(basename(fnFs2), "_"), `[`, 1)
sample.names2

sample.names2 <- sapply(strsplit(basename(fnRs2), "_"), `[`, 1)
sample.names2

sample.names3 <- sapply(strsplit(basename(fnFs3), "_"), `[`, 1)
sample.names3

sample.names3 <- sapply(strsplit(basename(fnRs3), "_"), `[`, 1)
sample.names3


#Inspect read quality profiles and Export plots for QC of unfiltered reads in Results/sub-directory

dir.create("Results/unfiltered.QCplots", recursive = T)

pdf("Results/unfiltered.QCplots/unfilt_f1.pdf", width=20, height=12)
plotQualityProfile(fnFs1)
dev.off()

pdf("Results/unfiltered.QCplots/unfilt_R1.pdf", width=20, height=12)
plotQualityProfile(fnRs1)
dev.off()

pdf("Results/unfiltered.QCplots/unfilt_f2.pdf", width=20, height=12)
plotQualityProfile(fnFs2)
dev.off()

pdf("Results/unfiltered.QCplots/unfilt_R2.pdf", width=20, height=12)
plotQualityProfile(fnRs2)
dev.off()

pdf("Results/unfiltered.QCplots/unfilt_f3.pdf", width=20, height=12)
plotQualityProfile(fnFs3)
dev.off()

pdf("Results/unfiltered.QCplots/unfilt_R3.pdf", width=20, height=12)
plotQualityProfile(fnRs3)
dev.off()


# Place filtered files in filtered/ subdirectory
filtFs1 <- file.path(path1, "trimmed1.QC", paste0(sample.names1, "_F1_filt.fastq"))
filtFs1
filtRs1 <- file.path(path1, "trimmed1.QC", paste0(sample.names1, "_R1_filt.fastq"))
filtRs1
filtFs2 <- file.path(path2, "trimmed2.QC", paste0(sample.names2, "_F2_filt.fastq"))
filtFs2
filtRs2 <- file.path(path2, "trimmed2.QC", paste0(sample.names2, "_R2_filt.fastq"))
filtRs2
filtFs3 <- file.path(path3, "trimmed3.QC", paste0(sample.names3, "_F3_filt.fastq"))
filtFs3
filtRs3 <- file.path(path3, "trimmed3.QC", paste0(sample.names3, "_R3_filt.fastq"))
filtRs3

names(filtFs1) <- sample.names1
names(filtRs1) <- sample.names1
names(filtFs2) <- sample.names2
names(filtRs2) <- sample.names2
names(filtFs3) <- sample.names3
names(filtRs3) <- sample.names3


# Filter and trim 

# We'll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE
# and  maxEE=2. The maxEE parameter sets the maximum number of "expected errors" allowed in a read,
# which is a better filter than simply averaging quality scores.
# Watch out with trunclen, reads have to overlap at the end, you have to try out.
# maxEE can be eased maxEE=c(2,5) if too many read are lost because of low quality.

trimmed_out1 <- filterAndTrim(fnFs1, filtFs1, fnRs1, filtRs1, truncLen=c(230,230),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
trimmed_out1
head(trimmed_out1)
df_trimmed_out1 <- as.data.frame(trimmed_out1)

trimmed_out2 <- filterAndTrim(fnFs2, filtFs2, fnRs2, filtRs2, truncLen=c(230,230),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
trimmed_out2
df_trimmed_out2 <- as.data.frame(trimmed_out2)

trimmed_out3 <- filterAndTrim(fnFs3, filtFs3, fnRs3, filtRs3, truncLen=c(230,230),
                              maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                              compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
trimmed_out3
df_trimmed_out3 <- as.data.frame(trimmed_out3)


# Check quality after filtering, use the new filtered files and output in a pdf file.

pdf("Results/1.QC/filt_f1.pdf", width=20, height=12)
plotQualityProfile(filtFs1)
dev.off()

pdf("Results/1.QC/filt_R1.pdf", width=20, height=12)
plotQualityProfile(filtRs1)
dev.off()

pdf("Results/1.QC/filt_f2.pdf", width=20, height=12)
plotQualityProfile(filtFs2)
dev.off()

pdf("Results/1.QC/filt_R2.pdf", width=20, height=12)
plotQualityProfile(filtRs2)
dev.off()

pdf("Results/1.QC/filt_f3.pdf", width=20, height=12)
plotQualityProfile(filtFs3)
dev.off()

pdf("Results/1.QC/filt_R3.pdf", width=20, height=12)
plotQualityProfile(filtRs3)
dev.off()

#### End ####


####Learn the Error Rates####

errF1 <- learnErrors(filtFs1, multithread=3) #set to FALSE after running command

errR1 <- learnErrors(filtRs1, multithread=3) #set to FALSE after running command

errF2 <- learnErrors(filtFs2, multithread=3) #set to FALSE after running command

errR2 <- learnErrors(filtRs2, multithread=3) #set to FALSE after running command

errF3 <- learnErrors(filtFs3, multithread=3) #set to FALSE after running command

errR3 <- learnErrors(filtRs3, multithread=3) #set to FALSE after running command

##Plot Error Rates
  #plotErrors(errF, nominalQ=TRUE)#
  #plotErrors(errR, nominalQ=TRUE)#

pdf("Results/error.rates/errF1.pdf", width=20, height=12)
plotErrors(errF1,nominalQ=TRUE )
dev.off()

pdf("Results/error.rates/errR1.pdf", width=20, height=12)
plotErrors(errR1, nominalQ=TRUE)
dev.off()

pdf("Results/error.rates/errF2.pdf", width=20, height=12)
plotErrors(errF2, nominalQ=TRUE)
dev.off()

pdf("Results/error.rates/errR2.pdf", width=20, height=12)
plotErrors(errR2, nominalQ=TRUE)
dev.off()

pdf("Results/error.rates/errF3.pdf", width=20, height=12)
plotErrors(errF3, nominalQ=TRUE)
dev.off()

pdf("Results/error.rates/errR3.pdf", width=20, height=12)
plotErrors(errR3, nominalQ=TRUE)
dev.off()

#### End ####


#### Dada and Merge####

# Apply the core sample inference algorithm to the the filtered and trimmed sequence data.
dadaFs1 <- dada(filtFs1, err=errF1, multithread=3)

dadaRs1 <- dada(filtRs1, err=errR1, multithread=3)

dadaFs2 <- dada(filtFs2, err=errF2, multithread=3)

dadaRs2 <- dada(filtRs2, err=errR2, multithread=3)

dadaFs3 <- dada(filtFs3, err=errF3, multithread=6)

dadaRs3 <- dada(filtRs3, err=errR3, multithread=6)

# Inspecting the returned dada-class object

dadaFs1[[8]]
dadaRs1[[8]]
dadaFs2[[4]]
dadaRs2[[4]]
dadaFs3[[4]]
dadaRs3[[4]]

#Merge paired reads

mergers1 <- mergePairs(dadaFs1, filtFs1, dadaRs1, filtRs1, verbose=TRUE) # min overlap is 12 as default, but can be adjusted
mergers2 <- mergePairs(dadaFs2, filtFs2, dadaRs2, filtRs2, verbose=TRUE) # min overlap is 12 as default, but can be adjusted
mergers3 <- mergePairs(dadaFs3, filtFs3, dadaRs3, filtRs3, verbose=TRUE) # min overlap is 12 as default, but can be adjusted

# Inspect the merger data.frame from the first sample

head(mergers1[[8]])
head(mergers2[[4]])
head(mergers3[[4]])

#Most of your reads should successfully merge. If that is not the case upstream parameters may need
#to be revisited: Did you trim away the overlap between your reads?


#Construct sequence table
seqtab1 <- makeSequenceTable(mergers1)  
seqtab2 <- makeSequenceTable(mergers2) 
seqtab3 <- makeSequenceTable(mergers3)

dim(seqtab1)
dim(seqtab2)
dim(seqtab3)

#If you wish to filter very short or long merged reads:-
seqtab1_1 <- seqtab1[,nchar(colnames(seqtab1)) %in% 300:387]
seqtab2_1 <- seqtab2[,nchar(colnames(seqtab2)) %in% 300:387]
seqtab3_1 <- seqtab3[,nchar(colnames(seqtab3)) %in% 300:387]

dim(seqtab1_1)
dim(seqtab2_1)
dim(seqtab3_1)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab1_1)))
table(nchar(getSequences(seqtab2_1)))
table(nchar(getSequences(seqtab3_1)))

# Get ASV lengths and counts
asv_lengths1 <- nchar(colnames(seqtab1_1))
length_df1 <- as.data.frame(table(Length = asv_lengths1))

asv_lengths2 <- nchar(colnames(seqtab2_1))
length_df2 <- as.data.frame(table(Length = asv_lengths2))

asv_lengths3 <- nchar(colnames(seqtab3_1))
length_df3 <- as.data.frame(table(Length = asv_lengths3))


# Basic ggplot of ASV lengths
ASV_LenDis3 <- ggplot(length_df3, aes(x = as.numeric(as.character(Length)), y = Freq)) +
  geom_col(fill = "skyblue") +
  geom_vline(xintercept = c(310, 340), colour = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Length Distribution of ASVs (300–387 bp)",
    x = "ASV length (bp)",
    y = "Number of ASVs"
  ) +
  theme_minimal()

# Print to RStudio viewer
print(ASV_LenDis3)

#create sub-directory to save plots
dir.create("Results/ASV_length_distribution", recursive = T)

# Save to file
ggsave("Results/ASV_length_distribution/ASV_length_distribution3.svg", plot = ASV_LenDis3, width = 10, height = 8, dpi = 300)


# Remove chimeras

seqtab.nochim1 <- removeBimeraDenovo(seqtab1, method="consensus", multithread=3, verbose=TRUE)
seqtab.nochim2 <- removeBimeraDenovo(seqtab2, method="consensus", multithread=3, verbose=TRUE)
seqtab.nochim3 <- removeBimeraDenovo(seqtab3, method="consensus", multithread=6, verbose=TRUE)

dim(seqtab.nochim1)
dim(seqtab.nochim2)
dim(seqtab.nochim3)

# Percentage of not chimeric reads

sum(seqtab.nochim1)/sum(seqtab1)
sum(seqtab.nochim2)/sum(seqtab2)
sum(seqtab.nochim3)/sum(seqtab3)

#Make a Plot of read counts before and after chimera removal using the R_script "Chimera_plot.R"

#### End ####

#### Track Pipeline ####
#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))

track1 <- cbind(trimmed_out1, sapply(dadaFs1, getN), sapply(dadaRs1, getN), sapply(mergers1, getN), rowSums(seqtab.nochim1))
track2 <- cbind(trimmed_out2, sapply(dadaFs2, getN), sapply(dadaRs2, getN), sapply(mergers2, getN), rowSums(seqtab.nochim2))
track3 <- cbind(trimmed_out3, sapply(dadaFs3, getN), sapply(dadaRs3, getN), sapply(mergers3, getN), rowSums(seqtab.nochim3))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
colnames(track3) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track1) <- sample.names1
rownames(track2) <- sample.names2
rownames(track3) <- sample.names3

head(track3)

# Save the track
write.csv(track1, file = "Results/1.QC/Pipeline_Track1_final.csv")
write.csv(track2, file = "Results/1.QC/Pipeline_Track2_final.csv")
write.csv(track3, file = "Results/1.QC/Pipeline_Track3_final.csv")

write.csv(sample.names1, file = "Results/1.QC/sample.names1.csv")
write.csv(sample.names2, file = "Results/1.QC/sample.names2.csv")
write.csv(sample.names3, file = "Results/1.QC/sample.names3.csv")

#### End ####

#### Proceed to the Decontamination step using the "Decontamination_16S.R script" ####

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

#### End ####

#### Combine the cleaned asv_tables ####

#Extract the unique ASV_column names
all_asvs_column <- unique(c(colnames(asv_table_clean1), colnames(asv_table_clean2), colnames(asv_table_clean3)))

#Make a function that finds missing ASV in each asv_table, compares them to all_asv_columns, and reorders them accordingly

pad_table <- function(asv_tab, all_asvs_column) {
  asv_tab <- as.data.frame(asv_tab)
  missing_cols <- setdiff(all_asvs_column, colnames(asv_tab))
  for (col in missing_cols) {
    asv_tab[, col] <- 0   # add new columns filled with zero
  }
  asv_tab[, all_asvs_column, drop=FALSE]    # reorder columns to match all_colnames
}

#call the "pad_table" function
# Standardize columns (ASVs) across all tables
asvtab1_padded <- pad_table(asv_table_clean1, all_asvs_column)
asvtab2_padded <- pad_table(asv_table_clean2, all_asvs_column)
asvtab3_padded <- pad_table(asv_table_clean3, all_asvs_column)


# Combine all samples (rows), keeping sample names distinct
combined_asvtab <- rbind(asvtab1_padded,asvtab2_padded, asvtab3_padded)

####End####

#### Assign taxonomy for Bacteria with GTDB####
ref_fasta1 = "C:/Users/Documents/Microbiome Analysis/16s_analysis/16S_Taxonomy_Databases/GTDB/GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz"
ref_fasta2 = "C:/Users/Documents/Microbiome Analysis/16s_analysis/16S_Taxonomy_Databases/GTDB/GTDB_bac120_arc53_ssu_r220_species.fa.gz"

#tryRC = TRUE -> reverse-complement orientation
set.seed(123456789)
combtaxa_GTDB <- assignTaxonomy(as.matrix(combined_asvtab), ref_fasta1, tryRC = TRUE, multithread = 8 )
unname(combtaxa_GTDB)

taxa.print_GTDB  <- combtaxa_GTDB # Removing sequence rownames for display only
rownames(taxa.print_GTDB) <- NULL
head(taxa.print_GTDB)


combtaxa_species_GTDB = addSpecies(combtaxa_GTDB,ref_fasta2, allowMultiple=TRUE)
taxa.print_spp_GTDB  <- combtaxa_species_GTDB  # Removing sequence rownames for display only
rownames(taxa.print_spp_GTDB) <- NULL
head(taxa.print_spp_GTDB)

####End####

#### Assign taxonomy for Bacteria with SILVA####
ref_fasta3 = "C:/Users/Documents/Microbiome Analysis/16s_analysis/16S_Taxonomy_Databases/SILVA/silva_nr99_v138.2_toGenus_trainset.fa.gz"
ref_fasta4 = "C:/Users/Microbiome Analysis/16s_analysis/16S_Taxonomy_Databases/SILVA/silva_v138.2_assignSpecies.fa.gz"

set.seed(123456789)
combtaxa_SILVA <- assignTaxonomy(as.matrix(combined_asvtab), ref_fasta3, tryRC = TRUE, multithread = 8 )
unname(combtaxa_SILVA)

taxa.print_SILVA <- combtaxa_SILVA # Removing sequence rownames for display only
rownames(taxa.print_SILVA) <- NULL
head(taxa.print_SILVA)


combtaxa_species_SILVA = addSpecies(combtaxa_SILVA,ref_fasta4, allowMultiple=TRUE)
taxa.print_spp_SILVA  <- combtaxa_species_SILVA  # Removing sequence rownames for display only
rownames(taxa.print_spp_SILVA) <- NULL
head(taxa.print_spp_SILVA)

####End####

#### Make Count and Taxa Tables for GTDB and SILVA derived taxa_table ####
# Giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs <- colnames(as.matrix(combined_asvtab))
asv_headers <- vector(dim(as.matrix(combined_asvtab))[2], mode="character")

for (i in 1:dim(as.matrix(combined_asvtab))[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}


# Making fasta of our final ASV seqs:

fasta_tab <- c(rbind(asv_headers, asv_seqs))

dir.create("Results/Fasta_tab", recursive = T)

write(fasta_tab, "Results/Fasta_tab/fasta_tab.fa")


# count table:

count_tab <- t(as.matrix(combined_asvtab))
row.names(count_tab) <- sub(">", "", asv_headers)
write.table(count_tab, "Results/ASV_tables/count_tab.tsv", sep="\t", quote=F, col.names=NA)


# tax table:

tax_tab_GTDB <- combtaxa_species_GTDB
row.names(tax_tab_GTDB) <- sub(">", "", asv_headers)
write.table(tax_tab_GTDB, "Results/Assigned_Taxonomy_GTDB/tax_tab_GTDB.tsv", sep="\t", quote=F, col.names=NA)


tax_tab_SILVA <- combtaxa_species_SILVA
row.names(tax_tab_SILVA) <- sub(">", "", asv_headers)
write.table(tax_tab_SILVA, "Results/Assigned_Taxonomy_Silva/tax_tab_SILVA.tsv", sep="\t", quote=F, col.names=NA)


#### End ####

#### Add samples metadata and match to count table####

#sample_info_tab <- read.table("Data/Sample_Info.txt",header=T, 
#                              row.names = 1, check.names=F)
library(openxlsx)
sample_info_tab <- read.xlsx("sample_metadata/16S_metadata.xlsx", rowNames = TRUE)

# Examine all tables
head(count_tab)
head(tax_tab_GTDB)
head(tax_tab_SILVA)
head(fasta_tab)
head(sample_info_tab)

# Name samples equally in all tables
colnames(count_tab) <- gsub("-", "_", colnames(count_tab)) 
rownames(sample_info_tab) <- gsub("-", "_", rownames(sample_info_tab))
sample_info_tab$Group <- gsub("-", "_", sample_info_tab$Group)

## Keep only controls
count_tab_samp <- count_tab[, !(colnames(count_tab) %in% c("NC_DNA_1", "NC_DNA_3", "NC_DNA_4", "NC_PCR", "PC_Fecal",
                                                    "PC_Gut", "PC_Mock", "PC_Mock1","NC_PCR1", "NC_DNA", "PC_Mock2","NC_DNA1"))]
head(count_tab_samp)
# head(count_tab)

## Examine rows with empty ASV
# Logical index for rows with all zeros (to be deleted)
deleted_rows_mask <- rowSums(count_tab_samp) == 0

# Display the actual rows to be deleted
deleted_rows <- count_tab_samp[deleted_rows_mask, ]
deleted_rows
deleted_rows <- as.data.frame(deleted_rows)

# # Delete rows with empty ASV
count_tab_samp1 <- count_tab_samp[!(rowSums(count_tab_samp) == 0),]
tax_tab_GTDB1 <- tax_tab_GTDB[rownames(tax_tab_GTDB) %in% rownames(count_tab_samp1),]
tax_tab_GTDB1 <- as.data.frame(tax_tab_GTDB1)
tax_tab_GTDB1$Species <- NULL
tax_tab_GTDB1<- as.matrix(tax_tab_GTDB1)


tax_tab_SILVA1 <- tax_tab_SILVA[rownames(tax_tab_SILVA) %in% rownames(count_tab_samp1),]


# Examine consistency in order between count_tab colnames and coldata rownames 
# (They have to be in the same order or Deseq2 won't work out)

all(rownames(sample_info_tab) %in% colnames(count_tab_samp1))
all(rownames(sample_info_tab) == colnames(count_tab_samp1)) # Mock community is not in sample_info_tab


all(rownames(tax_tab_GTDB1) %in% rownames(count_tab_samp1))
all(rownames(tax_tab_GTDB1) == rownames(count_tab_samp1))

all(rownames(tax_tab_SILVA1) %in% rownames(count_tab_samp1))
all(rownames(tax_tab_SILVA1) == rownames(count_tab_samp1))

gplots::venn(list(taxonomy=rownames(tax_tab_GTDB1), featuretable=rownames(count_tab_samp1)))
gplots::venn(list(taxonomy=rownames(tax_tab_SILVA1), featuretable=rownames(count_tab_samp1)))

#### End ####

#### Rarefaction curves ####

col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")

raremax <- min(rowSums(t(count_tab_samp1))) #minimum sequencing depth, you can tweak it as you wish with "mean" or "max"
raremax
rarefy <- rarefy(t(count_tab_samp1), raremax)
rarefy
# Option 1
rarecurve(t(count_tab_samp1), step = 100,  col = col, lwd=2, lty = lty, ylab = "ASVs", label = T)
abline(v=(raremax))

# Option 2
rarecurve(t(count_tab_samp1), step = 100, sample = raremax, col = col, lwd=2, lty = lty, ylab = "ASVs", label = T)

average_depth <- min(colSums(count_tab_samp1))
average_depth1 <- mean(rowSums(count_tab_samp1))

#### End ####



#### Phylogenetic tree ####
# Prepare the sequence table
phylog_count_tab <- t(as.matrix(combined_asvtab))
colnames(phylog_count_tab) <- gsub("-", "_", colnames(phylog_count_tab)) 
phylog_count_samp <- phylog_count_tab[, !(colnames(phylog_count_tab) %in% c("NC_DNA_1", "NC_DNA_3", "NC_DNA_4", "NC_PCR", "PC_Fecal",
                                                           "PC_Gut", "PC_Mock", "PC_Mock1","NC_PCR1", "NC_DNA", "PC_Mock2","NC_DNA1"))]

phylog_count_samp1 <- phylog_count_samp[!(rowSums(phylog_count_samp) == 0),]
phylog_count_samp2 <- t(as.matrix(phylog_count_samp1))


seqs <- getSequences(phylog_count_samp2)
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalOmega", type="dna", order="input")

aligned_asvs <- as(mult, "DNAStringSet")
names(aligned_asvs) <- seqs
writeXStringSet(aligned_asvs, filepath = file.path("aligned_asvs.fasta"), format = "fasta")

#-----> Fasttree in ubuntu
#tree <- read.tree(file.path(res.dir, "tree.nwk"))
tree <- read.tree("Results/phylog_tree/tree.nwk")

# The phangorn package is then used to construct a phylogenetic tree. 
# Here we first construct a neighbor-joining tree, and then fit a GTR+G+I maximum
# likelihood tree using the neighbor-joining tree as a starting point.

# phang.align <- as.phyDat(mult, type="dna", names=getSequence(phylog_count_samp2))
# dm <- dist.ml(phang.align)
# treeNJ <- NJ(dm) # Note, tip order != sequence order
# fit = pml(treeNJ, data=phang.align)

# Negative edges length changed to 0!
#fitGTR <- update(fit, k=4, inv=0.2)
#fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                   rearrangement = "stochastic", control = pml.control(trace = 0))

# will get back to you
# detach("package:phangorn", unload=TRUE)

#### End ####

#### Making our phyloseq object ####
# Examine all tables

head(count_tab_samp1)
head(tax_tab_GTDB1)
head(fasta_tab)
head(sample_info_tab)

physeq <- phyloseq(otu_table(count_tab_samp1, taxa_are_rows = T), #taxa_are_rows=F (if your taxa names on the column not the rows)
                   sample_data(sample_info_tab), 
                   tax_table(tax_tab_GTDB1))

## Adding ASV Fasta sequences and Phylogenetic tree to phyloseq object

asv_seqs <- colnames(as.matrix(phylog_count_samp2)) # use the latest modified count_table
dna <- Biostrings::DNAStringSet(asv_seqs)  # Making ASV Fasta sequences 
names(dna) <- taxa_names(physeq)

ph_tree = phy_tree(tree)            # Making Phylogenetic tree
taxa_names(ph_tree) = taxa_names(dna)

physeq_dna_tree <- merge_phyloseq(physeq, ph_tree, dna) #Merging  ASV Fasta sequences and Phylogenetic tree to phyloseq object

taxa_names(physeq_dna_tree) <- paste0("ASV", seq(ntaxa(physeq_dna_tree)))

#ph_tree = phy_tree(fitGTR$tree)            # Making Phylogenetic tree
#taxa_names(ph_tree) = taxa_names(dna)

#physeq_dna_tree <- merge_phyloseq(physeq, ph_tree, dna) #Merging  ASV Fasta sequences and Phylogenetic tree to phyloseq object

#taxa_names(physeq_dna_tree) <- paste0("ASV", seq(ntaxa(physeq_dna_tree)))

physeq
physeq_dna_tree
physeq = physeq_dna_tree

#### End ####

## Phyloseq Object Filtering ## First filtering step for low count samples and NA phyla ##

#### Remove samples with less than 100 total reads ####
sample_sums(physeq) #Nr of reads per Sample
min(sample_sums(physeq))
#No samples below 100 reads, thus no filtering necessary

#physeq_above100reads = prune_samples(sample_sums(physeq)>=100, physeq)

#physeq = physeq_above100reads

#### End ####

#### Plot Count matrix ####
# Extract total reads per sample
phy_count <- data.frame(Groups = physeq@sam_data$Group,
                        TotalReads = sample_sums(physeq))

# Summarize total reads by group (replace 'GroupVariable' with your variable name)
group_counts <- sample_counts %>%
  group_by(GroupVariable) %>%
  summarise(TotalReads = sum(TotalReads))

# Plot total reads per group
count_matrix <- ggplot(phy_count, aes(x = Groups, y = TotalReads)) +
  geom_boxplot(width = 0.4, fill = "steelblue") +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.3))) +
  labs(title = "Count matrix") +
  ##xlab("Groups") +
  ylab("Total Read Counts") +
  theme_minimal() +
  theme(strip.background = element_blank(),
        text = element_text(size=14),
        axis.title.y = element_text(size=14, margin = margin(r = 15),face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.text.x.bottom = element_text(angle = 45, size = 14,hjust = 1, vjust = 1, face = "bold"),
        panel.border = element_rect(linewidth = 2),
        strip.text = element_text(size=13, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")
  ) 
count_matrix

ggsave("Results/count_matrix.svg", plot = count_matrix, width = 5, height = 4, units = "in")

#### Remove controls and save controls in separate phyloseq object ####

# Phyloseq object with controls

physeq_K <- prune_samples(sample_names(physeq) %in% c("NC_DNA_1", "NC_DNA_3", "NC_DNA_4", "NC_PCR", "PC_Fecal",
                                                      "PC_Gut", "PC_Mock"), physeq)
physeq_K <- prune_taxa(taxa_sums(physeq_K) != 0, physeq_K)

# Remove controls from physeq

physeq_noK <- prune_samples(!(sample_names(physeq) %in% c("NC_DNA_1", "NC_DNA_3", "NC_DNA_4", "NC_PCR", "PC_Fecal",
                                                          "PC_Gut", "PC_Mock")), physeq)
physeq_noK <- prune_taxa(taxa_sums(physeq_noK) != 0, physeq_noK)

physeq = physeq_noK

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

dir.create("Results/RData")

save.image("Results/RData/1.physeq.original.RData")

#### End ####


## Remove Batch-Effect from physeq object(physeq)...

## Filter low Abundance Taxa and count table normalization ##

#### Remove Phyla with NA####

rank_names(physeq)

table(tax_table(physeq)[, "Phylum"], exclude = NULL)

physeq_oNA <- subset_taxa(physeq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#physeq = physeq_oNA

table(tax_table(physeq_oNA)[, "Phylum"], exclude = NULL)


#### End ####

####Define prevalence of each taxa (in how many samples did each taxa appear at least once)####

prev0 = apply(X = otu_table(physeq_oNA),
              MARGIN = ifelse(taxa_are_rows(physeq_oNA), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(physeq_oNA),
                    tax_table(physeq_oNA))



#save ASV Prevalence and Abundance table before filtering

write.table(prevdf, "Results/asv_prevalancedf.tsv", sep="\t", quote=F, col.names=NA)

#Plot Taxa prevalence v. total counts. Each point is a different taxa. 

ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(physeq_oNA),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#ggsave("Results/ASV_length_distribution/ASV_length_distribution3.svg", plot = ASV_LenDis3, width = 10, height = 8, dpi = 300)
ggsave("Results/Prev_Totalabundance.tiff", width = 12, height = 8, dpi = 600 )

#### End ####

#### Remove taxa not seen more than 3 times in at least 5% of the samples #### 
# This protects against an OTU with small mean & trivially large C.V.
# Setting filter parameters :

countperphyla = 3
Samplepercentage = 0.05

physeq_filtered = filter_taxa(physeq_oNA, function(x) sum(x > countperphyla) > (Samplepercentage*length(x)), TRUE)
physeq_oNA
physeq_filtered
#physeq = physeq_filtered

#### End ####

#### Normalize number of reads in each sample using median sequencing depth.####

total = median(sample_sums(physeq_filtered))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq_mednorm = transform_sample_counts(physeq_filtered, standf)
physeq_mednorm

# Transform to relative abundance. Save as new object.
physeq_re = transform_sample_counts(physeq_mednorm, function(x){x / sum(x)})
physeq_re

#### End ####

#### Exploratory plots after filtering and normalization ####
# Check individual phylum Abundance
# Abundance value transformation function

library(rlang)

plot_abundance = function(physeq_filtered, ylabn = "",
                          Facet = "Phylum",
                          Color = "Phylum",
                          n = NULL){
  mphyseq = psmelt(physeq_filtered)
  mphyseq <- subset(mphyseq, Abundance > 0)
  
  ggplot(data = mphyseq,
         mapping = aes(x = Group,
                       y = Abundance,
                       color = !!sym(Color),
                       fill = !!sym(Color))) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    geom_boxplot(color = "black", size = 0.5) +
    facet_wrap(vars(!!sym(Facet)), nrow = n) +
    ylab(ylabn) +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
#plot the abundance values before and after transformation

pl_ab_original  = plot_abundance(physeq_filtered,"Original Abundances")
print(pl_ab_original)
ggsave("Results/Abundance_plots/Phylum_Original_abundances.tiff", width = 12, height = 8, dpi = 600 )

pl_ab_original_norm  =plot_abundance(physeq_mednorm,"Normalized to squencing depth Abundances")
print(pl_ab_original_norm)
ggsave("Results/Abundance_plots/Phylum_Norm_to_seq_depth_Abundances.tiff", width = 12, height = 8, dpi = 600)

pl_ab_original_norm_re  =plot_abundance(physeq_re,"Normalized Relative Abundances")
print(pl_ab_original_norm_re)
ggsave("Results/Abundance_plots/Phylum_Norm_Relative_Abundances.svg", width = 12, height = 8, dpi = 600)

grid.arrange(pl_ab_original, pl_ab_original_norm, pl_ab_original_norm_re)

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image("Results/RData/2.phyloseq.filtered.RData")

#### End ####
