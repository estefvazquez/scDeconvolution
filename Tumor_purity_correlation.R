# Make a simple stacked barplot for visualization

stacked_df <- rbind(bisque_estimates, bp_estimates) # Patients as rows and fractions as columns
stacked_df <- stacked_df %>%
  pivot_longer(cols = c("Malignant","Macrophages","CD8_T","CD4_T","Cycling_T","Neutrophils","NK","CAF",
                        "Monocytes","Dendritic","Mast cells","Vascular","Endothelial","Endocrine",
                        "B","Plasma","Stellate","Schwann","Treg","Ductal","Acinar"), 
               names_to = "cell_type", 
               values_to = "fraction")

stacked_df <- stacked_df %>%
  group_by(Method, cell_type) %>%
  summarise(mean_fraction = mean(fraction), .groups = 'drop')

ggplot(stacked_df, aes(x = Method, y = mean_fraction, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell Proportions by Method",
       x = "Method",
       y = "Mean Proportion",
       fill = "Cell Type") +
  theme_minimal()+scale_fill_manual(values = c("Acinar" = "royalblue",
                                               "B" = "#F72585",
                                               "Stellate" = "darkorange",
                                               "Vascular" = "darkred",
                                               "CAF" = "#986960",
                                               "Malignant" = "#7209B7",
                                               "Dendritic" = "dimgray",
                                               "Macrophages" = "green",
                                               "Endothelial" = "red",
                                               "Neutrophils" = "thistle",
                                               "Monocytes" = "darkseagreen",
                                               "CD4_T" = "cornflowerblue",
                                               "CD8_T" = "darkgreen",
                                               "NK" = "darkgoldenrod",
                                               "Cycling_T" = "black",
                                               "Treg" = "navy",
                                               "Ductal" = "lightcoral",
                                               "Plasma" = "orchid",
                                               "Mast cells" = "limegreen",
                                               "Schwann" = "darkturquoise",
                                               "Endocrine" = "rosybrown"))

# Load clinical data --> This clinical metadata contains information for CNV_status (HIGH or LOW); ABSOLUTE purity score was available in the clinical metadata downloaded from TCGAbiolinks
# ABSOLUTE infers tumor purity and malignant cell ploidy directly from analysis of somatic DNA alterations (for more details check https://pmc.ncbi.nlm.nih.gov/articles/PMC4383288/)
# stromal and immune columns represent scores from the ESTIMATE package that infers the levels of stromal and immune cells presents in bulk RNA data
# ESTIMATE_purity score represents the levels of malignant cell fraction (similar to the ABSOLUTE algorithm, but this method was developed for bulk RNA)
# As expected, immune score negatively correlates with ESTIMATE_purity score

clindata <- readRDS('clinical_bulk_2025.rds')
                                                
# Remove additional suffixes to patient barcode

rownames(bisque_estimates) <- substr(rownames(bisque_estimates), 1, 12)
rownames(bp_estimates) <- substr(rownames(bp_estimates), 1, 12)

# Merge clinical data with cell fractions to correlate tumor purity with malignant cell estimation from deconvolution

bisque_estimates$Patient <- rownames(bisque_estimates)
bp_estimates$Patient <- rownames(bp_estimates)

purity_bisque <- merge(bisque_estimates, dwls1, by = "Patient")
purity_bayes <- merge(bp_estimates, dwls2, by = "Patient")

purity_bisque <- na.omit(purity_bisque) # Remove one patient that has CNV_status as 'NA'
purity_bayes <- na.omit(purity_bayes) # Remove one patient that has CNV_status as 'NA'

# Make scatter plot for tumor purity and malignant estimation

gg1 <- ggplot(purity_bisque, aes(x = ABSOLUTE_purity, y = Malignant)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add a linear regression line
  labs(title = "Correlation between Tumor Purity and Malignant Cell Fraction - Bisque",
       x = "Tumor Purity",
       y = "Malignant Cell Fraction") +
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 0.05, label.y = 0.95, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")))

gg2 <- ggplot(purity_bayes, aes(x = ABSOLUTE_purity, y = Malignant)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add a linear regression line
  labs(title = "Correlation between Tumor Purity and Malignant Cell Fraction - BayesPrism",
       x = "Tumor Purity",
       y = "Malignant Cell Fraction") +
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 0.05, label.y = 0.95, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")))

gg1+gg2 # Side by side plot

# Which method showed a stronger correlation with purity scores ?
# What else could we do with this data ?
# For more details keep an eye on the github repo for any updates on downstream analysis and benchmarking your own deconvolution