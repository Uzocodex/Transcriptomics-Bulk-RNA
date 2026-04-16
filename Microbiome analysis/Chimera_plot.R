####Introduction####

# title: "Chimera_plot"
# author: "Uzochukwu Gospel Ukachukwu"
# date: "2026-03-15"
# R_version: "4.4.3"
# output: R_script 

#### End ####

library(ggplot2)
library(reshape2)

# Create data frame with before/after read counts per sample (replace with your own seqtab files)
reads_df3 <- data.frame(
  Sample   = rownames(seqtab3),
  Before   = rowSums(seqtab3),
  After    = rowSums(seqtab.nochim3)
)

# Reshape to long format for ggplot
reads_long3 <- melt(reads_df3, id.vars = "Sample", 
                   variable.name = "Stage", value.name = "ReadCount")

# Plot
remov.chimplot3 <- ggplot(reads_long3, aes(x = Sample, y = ReadCount, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Read Counts per Sample Before & After Chimera Removal",
       x = "Sample", y = "Read Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Show in viewer
print(remov.chimplot3)

# Save to file
ggsave("Results/chimera_remov_plot/read_counts_before_after_chimera3.png", plot = remov.chimplot3,
       width = 12, height = 6, dpi = 300)

#### End ####
