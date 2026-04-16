####Introduction####

# title: "Combining all asv_tables"
# author: "Uzochukwu Gospel Ukachukwu"
# date: "2026-03-15"
# R_version: "4.4.3"
# output: R_script 

#### End ####

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