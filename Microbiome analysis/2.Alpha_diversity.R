####Introduction####

# title: "Alpha diversity"
# author: "Uzochukwu Gospel Ukachukwu"
# date: "2026-03-15"
# R_version: "4.4.3"
# output: R_script 

 
#### Package setup ####

library(phyloseq)
library(ggplot2)
library(ggpubr)
library(grid)
library(openxlsx)

#### End ####

#### Load output from DADA2 pipeline with original phyloseq object ####

#load("Results/RData/1.physeq.original.RData")

# Create directories for results
dir.create("Results/2.Alpha_stats")
dir.create("Results/3.Alpha_plots")

#### End ####

#### Alpha-diversity (Richness and diversity estimates) ####
# Should be done in non filtered physeq
sample_variables(physeq_oNA)

# Visualize alpha-diversity on unfiltered phyloseq object
#Plot richness
P1 = plot_richness(physeq_oNA, x="Group", color = "Group", title = "Alpha Diversity", measures=c("Observed", "Chao1", "Shannon", "InvSimpson"), nrow = 2)
print(P1)


P1.1 = P1+
  geom_boxplot(alpha = 0.5, width = 0.5, linewidth = 1) +
  theme_bw() + xlab(NULL) +
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

P1.1

ggsave(path = "Results/3.Alpha_plots", filename = "Alpha Diversity.tiff", 
       plot = P1.1, width = 10, height = 6, units = "in", dpi = 600)

# Violin plot

P1.2 <- P1 + geom_violin(alpha = 0.5) +
  theme_bw() +
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold")
        ) + 
  geom_point(position = position_dodge(width = 0.75))
P1.2

ggsave(path = "Results/3.Alpha_plots", filename = "Alpha_Diversity_violin.svg", 
       plot = P1.2, width = 10, height = 6, units = "in", dpi = 600)


#### End ####

#### Richness statistics####
# Has to be counts not relative abundances

richness_data <- estimate_richness(physeq_oNA, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
#rownames(rich) <- gsub("X", "", rownames(rich)) %>% gsub("\\.", "-", .)

# Merge with sample_data for easy plotting
#richness_data$Group <- sample_data(physeq_oNA)$Group
richness_data$original_order <- seq_len(nrow(richness_data)) #preserve order of rownames
richness_data <- merge(richness_data, sample_data(physeq_oNA), by = 0)
richness_data<- richness_data[order(richness_data$original_order), ]
richness_data$original_order <- NULL  # remove helper column if desired
# Save table with alpha diversity

write.xlsx(richness_data, "Results/2.Alpha_stats/richness.xlsx")

#Check if the indices are normally distributed with "Normality_check.R" script and come back to this one

#Shannon and InvSimpson are not normally distributed and unpaired -> Kruskal-Wallis / U-test
#Observed, Chao1  are normally distributed and unpaired -> one way Anova / unpaired t test
#However, we performed all tests in all indexes

rich <- richness_data

aov <- list()
for(v in c("Group", "Bristol.Scale", "Medication", "Diet", "Alcohol")){
  for(i in c("Observed", "Chao1", "Shannon", "InvSimpson")){
     aov_test <- aov(rich[,i] ~ rich[,v], data = rich)
     aov_sum <- summary(aov_test)
     tab <- c(i, v, aov_sum[[1]]$`Sum Sq`, aov_sum[[1]]$`Mean Sq`, aov_sum[[1]]$`F value`, aov_sum[[1]]$`Pr(>F)`)
     aov[[paste0(i, "_", v)]] <- tab
     }}
aov <- as.data.frame(do.call(rbind, aov)) ; aov <- aov[ , colSums(is.na(aov)) == 0]
colnames(aov) <- c("Alpha_index","Variable", "Sum Sq [i]", "Sum Sq Residuals", "Mean Sq [i]", "Mean Sq Residuals", "F value", "Pr(>F)")

write.table(aov, "Results/2.Alpha_stats/One_way_Anova.txt")
write.table(aov, "Results/2.Alpha_stats/One_way_Anova.csv")

# post-hoc test (Tukey Honest Significant Difference test)
library("agricolae")
aov_test <- aov(Shannon ~ Alcohol, data = rich)
hsd_test <- TukeyHSD(aov_test) # require the agricolae package  
hsd_res <- HSD.test(aov_test, "Alcohol", group=T)$groups  

##Continue from here###
#### Paired t-test ####
t.test(rich[,"Observed"] ~ rich[,"Sex"])  
#, paired = F) #independent samples
colnames(rich)[13] <- "Pre.existing.conditions" #renamed the 13 column

rich[,"Pre.existing.conditions"] = gsub(" ", "", rich[,"Pre.existing.conditions"])
rich[,"Smoke"] = gsub(" ", "", rich[,"Smoke"])
t <- list()
for(v in c("Sex", "Pre.existing.conditions", "Smoke")){    #removed "Condition"
   for(n in c("Observed", "Chao1", "Shannon", "InvSimpson")){
       t_test <- t.test(rich[!is.na(rich[,v]),i] ~ rich[!is.na(rich[,v]),v]) #, paired = F)
       tab <- c(t_test$method, t_test$statistic, t_test$parameter, t_test$estimate, t_test$p.value)
       t[[paste0(n, "_", v)]] <- tab
     }}
t <- as.data.frame(do.call(rbind, t)) 
colnames(t) <- c("Method", colnames(t)[2:3], "mean in group 1", "mean in group 2", "p-value")
write.xlsx(t, "Results/2.Alpha_stats/Paired_t_test.xlsx", rowNames = T)



#### ####

# Mann-Whitney-U-Test 

mwut <- list()

for (n in c("Observed", "Chao1", "Shannon", "InvSimpson")) {
  U_test <- wilcox.test(rich[, n] ~ rich[, "Group"], data = rich, paired = F, exact = T)
  z <- abs(qnorm(U_test$p.value/2))
  r <- z/sqrt(nrow(rich))
  
  tab <- c(U_test$method, paste("data:", n, "and Group"), U_test$statistic, U_test$p.value, r)
  mwut[[paste0(n)]] <- tab
  
}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variables", "Statistic", "p-value", "Effect size")

mwut$p.adj <- p.adjust(mwut$`p-value`, method = "BH")

write.xlsx(mwut, "Results/2.Alpha_stats/Mann_Whitney_u_test.xlsx", rowNames = T)

# Wilcoxon-Test 

wt <- list()

for (n in c("Observed", "Chao1", "Shannon", "InvSimpson")) {
  W_test <- wilcox.test(rich[, n] ~ rich[, "Group"], data = rich, paired = T, exact = T)
  z <- abs(qnorm(W_test$p.value/2))
  r <- z/sqrt(nrow(rich))
  
  tab <- c(W_test$method, paste("data:", n, "and Group"), W_test$statistic, W_test$p.value, r)
  wt[[paste0(n)]] <- tab
  
}

wt <- as.data.frame(do.call(rbind, wt))
colnames(wt) <- c("Test", "Variables", "Statistic", "p-value", "Effect size")

wt$p.adj <- p.adjust(wt$`p-value`, method = "BH")
wt$p.adj_notnorm <- c(NA, NA, p.adjust(wt$`p-value`[3:4], method = "BH"))

write.xlsx(wt, "Results/2.Alpha_stats/Wilcoxon_test.xlsx", rowNames = T)

# One-Way Anova

#https://scienceparkstudygroup.github.io/microbiome-lesson/aio/index.html

aov <- list()
for(v in c("Group", ""))
for(n in c("Observed", "Chao1", "Shannon", "InvSimpson")){
  aov_test <- aov(rich[,n] ~ rich[, "Group"], data = rich)
  aov_sum <- summary(aov_test)
  aov_test$df.residual
  tab <- c(n, aov_sum[[1]]$`Sum Sq`, aov_sum[[1]]$`Mean Sq`, aov_sum[[1]]$`F value`, aov_sum[[1]]$`Pr(>F)`)
  aov[[n]] <- tab
}
aov <- as.data.frame(do.call(rbind, aov)) ; aov <- aov[ , colSums(is.na(aov)) == 0]
colnames(aov) <- c("Variable", "Sum Sq [i]", "Sum Sq Residuals", "Mean Sq [i]", "Mean Sq Residuals", "F value", "Pr(>F)")

aov$p.adj <- p.adjust(aov$`Pr(>F)`, method = "BH")

write.xlsx(aov, "Results/2.Alpha_stats/One_way_Anova.xlsx")

# Paired t-test

t <- list()
for(n in c("Observed", "Chao1", "Shannon", "InvSimpson")){
  t_test <- t.test(rich[rich$Group == "HCB", n], rich[rich$Group == "MPYG", n], paired = T)
  tab <- c(t_test$method, t_test$statistic, t_test$parameter, t_test$estimate, t_test$p.value)
  t[[n]] <- tab
}

t <- as.data.frame(do.call(rbind, t)) 
colnames(t) <- c("Method", colnames(t)[2:4], "p-value")

t$p.adj <- p.adjust(t$`p-value`, method = "BH")
t$p.adj_norm <- c(p.adjust(t$`p-value`[1:2], method = "BH"), NA, NA)

write.xlsx(t, "Results/2.Alpha_stats/Paired_t_test.xlsx", rowNames = T)

#### End ####




#### Alpha-diversity (Richness and diversity estimates) ####

# Create plot with the right statistics

# Visualize alpha-diversity on unfiltered phyloseq object

P2 = plot_richness(physeq_oNA, x="Group", color = "Group", measures=c("Shannon", "InvSimpson"),nrow = 2)
#P2 = plot_richness(physeq_oNA, x="Group", color = "Group", measures=c("InvSimpson"),nrow = 1)

# Add a boxplot layer to the same data used in plot_richness (available as P2$data)
P2 <- P2+
  geom_boxplot(data = P2$data, aes(x = Group, y = value, color = Group),
               alpha = 0.5, width = 0.5, linewidth = 1) +
  theme_bw() + xlab(NULL) +
  theme(strip.background = element_blank(),
        text = element_text(size=14),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.text.x.bottom = element_text(angle = 45, size = 14,hjust = 1, vjust = 1, face = "bold"),
        panel.border = element_rect(linewidth = 2),
        strip.text = element_text(size=13, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  stat_kruskal_test(label.y = 1.8)
 
# stat_compare_means(method = "t.test", paired = T, label.y = 1.8) 
  
P2  

P3 = plot_richness(physeq_oNA, x="Group", color = "Group", measures=c("Observed", "Chao1"), nrow = 2)
#P3 = plot_richness(physeq_oNA, x="Group", color = "Group", measures=c("Chao1"), nrow = 1)

# Add a boxplot layer to the same data used in plot_richness (available as P2$data)
P3 <- P3+
  geom_boxplot(data = P3$data, aes(x = Group, y = value, color = Group),
               alpha = 0.5, width = 0.5, linewidth = 1) +
  theme_bw() + xlab(NULL) +
  theme(strip.background = element_blank(),
        text = element_text(size=14),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.text.x.bottom = element_text(angle = 45, size = 14,hjust = 1, vjust = 1, face = "bold"),
        panel.border = element_rect(linewidth = 2),
        strip.text = element_text(size=13, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  stat_compare_means(method = "anova", label.y = 1.8)
  
  
P3

# Combine all plots
combined_plot <- ggarrange(P2, P3,
                           #ncol = 2,      # Number of columns
                           #nrow = 2,      # Number of rows
                           #labels = c("A", "B"),   # Optional labels for plots
                           common.legend = TRUE,   # If you want a common legend
                           legend = "right")         # Position of the common legend

#combined_plot <- annotate_figure(combined_plot,
                                 #top = text_grob("Alpha Diversity", 
                                                 #face = "bold", size = 16),
                                 #left = text_grob("Alpha Diversity Measure", 
                                                  #rot = 90, size = 14, face = "bold"))
# Print the combined plot
print(combined_plot)  

ggsave(path = "Results/3.Alpha_plots", filename = "Alpha Diversity-stat.svg", 
       plot = combined_plot, width = 10, height = 6, units = "in", dpi = 600)

ggsave(path = "Results/3.Alpha_plots", filename = "Alpha Diversity-stat_presentation.tiff", 
       plot = combined_plot, width = 10, height = 6, units = "in", dpi = 600)

#### End ####