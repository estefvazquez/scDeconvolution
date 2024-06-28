#Subset matrix to contain top25 marker genes for each cluster to run deconvolution
#Some methods benefit more from a signature matrix than using the full matrix
#Load markers
#Running Bayes first
markers = read.csv('D:/Scanpy/markers_for_decon.csv')
markers$X <- NULL
df1 <- markers[order(markers$scores,decreasing=TRUE),]
ref_markers <- df1 %>%
  group_by(group) %>%
  dplyr::filter(logfoldchanges > 3) %>%
  slice_head(n = 25) %>%
  ungroup()

genes <- make.unique(ref_markers$names)

#Load reference data
sobj <- readRDS('C:/Users/Gabriel/Desktop/sc_ref_pdac_seurat.rds')
metadata <- read.csv('C:/Users/Gabriel/Desktop/metadata_complete.csv')
sobj@meta.data <- metadata

#Subset seurat
sobj_subset <- subset(sobj, features = genes)
dim(sobj_subset)

#Compare Bayes and Bisque results full matrix
bisque_res <- readRDS('C:/Users/Gabriel/Desktop/scdeconv_bisque.rds')
bp_res <- readRDS('C:/Users/Gabriel/Desktop/sc_deconv_bp.rds')
bisque_estimates <- bisque_res$bulk.props
bisque_estimates <- as.data.frame(t(bisque_estimates))
bp_estimates <- get.fraction(bp_res, which.theta = 'final', state.or.type = 'state')
bp_estimates <- as.data.frame(bp_estimates)


bisque_estimates$Method <- "Bisque"
bp_estimates$Method <- "BayesPrism"
bisque_estimates$Patient <- rownames(bisque_estimates)
bp_estimates$Patient <- rownames(bp_estimates)
df_combined <- rbind(bp_estimates, bisque_estimates)


# Melt the dataframe to long format
df_long <- melt(df_combined, id.vars = c("Patient", "Method"), 
                variable.name = "CellType", value.name = "Fraction")


# Create the paired barchart
plot<-ggplot(df_long, aes(x = Patient, y = Fraction, fill = Method)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ CellType, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
  labs(title = "Comparison of Cell-Type Fractions between Methods",
       x = "Patient", y = "Fraction") +
  scale_fill_manual(values = c("Bisque" = "blue", "BayesPrism" = "red"))

print(plot)


#Create a stacked barplor for better visualization
stacked_df <- rbind(bisque_estimates, bp_estimates)
stacked_df <- stacked_df %>%
  pivot_longer(cols = c("Malignant","Macrophages","CD8_T","CD4_T","Cycling_T","Neutrophils","NK","CAF",
                        "Monocytes","Dendritic","Mast cells","Vascular","Endothelial","Endocrine",
                        "B","Plasma","Stellate","Schwann","Treg","Ductal","Acinar"), 
               names_to = "cell_type", 
               values_to = "fraction")

stacked_df <- stacked_df %>%
  group_by(Method, cell_type) %>%
  summarise(mean_fraction = mean(fraction), .groups = 'drop')

p <- ggplot(stacked_df, aes(x = Method, y = mean_fraction, fill = cell_type)) +
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

#Create a stacked barplor for better visualization (top25 markers results)
stacked_df25 <- rbind(bisque_estimates25, bp_estimates25)
stacked_df25 <- stacked_df25 %>%
  pivot_longer(cols = c("Malignant","Macrophages","CD8_T","CD4_T","Cycling_T","Neutrophils","NK","CAF",
                        "Monocytes","Dendritic","Mast cells","Vascular","Endothelial","Endocrine",
                        "B","Plasma","Stellate","Schwann","Treg","Ductal","Acinar"), 
               names_to = "cell_type", 
               values_to = "fraction")

stacked_df25 <- stacked_df25 %>%
  group_by(Method, cell_type) %>%
  summarise(mean_fraction = mean(fraction), .groups = 'drop')

p25 <- ggplot(stacked_df25, aes(x = Method, y = mean_fraction, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell Proportions by Method (Top25 markers)",
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


# Display the plot
print(p25)

#Compare Bayes and Bisque results top25 matrix
bisque_res25 <- readRDS('D:/Scanpy/bisque_top25.rds')
bp_res25 <- readRDS('D:/Scanpy/bayes_top25.rds')
bisque_estimates25 <- bisque_res25$bulk.props
bisque_estimates25 <- as.data.frame(t(bisque_estimates25))
bp_estimates25 <- bp_res25@posterior.initial.cellState@theta
bp_estimates25 <- as.data.frame(bp_estimates25)


bisque_estimates25$Method <- "Bisque"
bp_estimates25$Method <- "BayesPrism"
bisque_estimates25$Patient <- rownames(bisque_estimates25)
bp_estimates25$Patient <- rownames(bp_estimates25)
df_combined25 <- rbind(bp_estimates25, bisque_estimates25)


# Melt the dataframe to long format
df_long25 <- melt(df_combined25, id.vars = c("Patient", "Method"), 
                  variable.name = "CellType", value.name = "Fraction")


# Create the paired barchart
plot25<-ggplot(df_long25, aes(x = Patient, y = Fraction, fill = Method)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ CellType, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
  labs(title = "Comparison of Cell-Type Fractions between Methods (Top 25 markers)",
       x = "Patient", y = "Fraction") +
  scale_fill_manual(values = c("Bisque" = "blue", "BayesPrism" = "red"))

print(plot25)

bayes_purity25_malignant <- merge(clinical_purity25, bp_estimates25, by = "Patient")
bisque_purity25_malignant <- merge(clinical_purity25, bisque_estimates25, by = "Patient")
bayes_purity25_malignant$Tumor_Ploidy_score <- NULL
bayes_purity25_malignant$TumorPurity <- NULL
bayes_purity25_malignant$Pathologist_tumor_cellularity <- NULL

bisque_purity25_malignant$Tumor_Ploidy_score <- NULL
bisque_purity25_malignant$TumorPurity <- NULL
bisque_purity25_malignant$Pathologist_tumor_cellularity <- NULL

bisque_purity25_malignant <- bisque_purity25_malignant %>% filter(bisque_purity25_malignant$CNV_classification == "High" |
                                                                    bisque_purity25_malignant$CNV_classification == "Low")

bayes_purity25_malignant <- bayes_purity25_malignant %>% filter(bayes_purity25_malignant$CNV_classification == "High" |
                                                                  bayes_purity25_malignant$CNV_classification == "Low")


bayes_stacked_df25 <- bayes_purity25_malignant %>%
  pivot_longer(cols = c("Malignant_low_cnv","Malignant_high_cnv","Macrophages","CD8_T","CD4_T","Cycling_T","Neutrophils","NK","CAF",
                        "Monocytes","Dendritic","Mast cells","Vascular","Endothelial","Endocrine",
                        "B","Plasma","Stellate","Schwann","Treg","Ductal","Acinar"), 
               names_to = "cell_type", 
               values_to = "fraction")

bayes_stacked_df25 <- bayes_stacked_df25 %>%
  group_by(CNV_classification, cell_type) %>%
  summarise(mean_fraction = mean(fraction), .groups = 'drop')

bayesp25 <- ggplot(bayes_stacked_df25, aes(x = CNV_classification, y = mean_fraction, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell Proportions by CNV class (Top25 markers) - BayesPrism",
       x = "CNV",
       y = "Mean Proportion",
       fill = "Cell Type") +
  theme_minimal()+scale_fill_manual(values = c("Acinar" = "royalblue",
                                               "B" = "#F72585",
                                               "Stellate" = "darkorange",
                                               "Vascular" = "darkred",
                                               "CAF" = "#986960",
                                               "Malignant_high_cnv" = "#7209B7",
                                               "Malignant_low_cnv" = "gold",
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


# Display the plot
print(bayesp25)


bisque_stacked_df25 <- bisque_purity25_malignant %>%
  pivot_longer(cols = c("Malignant_low_cnv","Malignant_high_cnv","Macrophages","CD8_T","CD4_T","Cycling_T","Neutrophils","NK","CAF",
                        "Monocytes","Dendritic","Mast cells","Vascular","Endothelial","Endocrine",
                        "B","Plasma","Stellate","Schwann","Treg","Ductal","Acinar"), 
               names_to = "cell_type", 
               values_to = "fraction")

bisque_stacked_df25 <- bisque_stacked_df25 %>%
  group_by(CNV_classification, cell_type) %>%
  summarise(mean_fraction = mean(fraction), .groups = 'drop')

bisquep25 <- ggplot(bisque_stacked_df25, aes(x = CNV_classification, y = mean_fraction, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell Proportions by CNV class (Top25 markers) - BisqueRNA",
       x = "CNV",
       y = "Mean Proportion",
       fill = "Cell Type") +
  theme_minimal()+scale_fill_manual(values = c("Acinar" = "royalblue",
                                               "B" = "#F72585",
                                               "Stellate" = "darkorange",
                                               "Vascular" = "darkred",
                                               "CAF" = "#986960",
                                               "Malignant_high_cnv" = "#7209B7",
                                               "Malignant_low_cnv" = "gold",
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


# Display the plot
print(bisquep25)


#Get clinical data with ABSOLUTE purity scores for each patient and correlate with Malignant estimates from Bayes and Bisque
clinical <- readRDS('C:/Users/Gabriel/Desktop/tcga_clindata_downsampled.rds')
clinical <- clinical[,c("paper_ABSOLUTE.Purity","barcode","vital_status","days_to_death","days_to_last_follow_up")]
clinical_purity <- clinical[,c("paper_ABSOLUTE.Purity","barcode")]
clinical_purity$barcode <- NULL
bisque_estimates_malignant <- as.data.frame(bisque_estimates[,"Malignant"])
bp_estimates_malignant <- as.data.frame(bp_estimates[,"Malignant"])

clinical_purity$Patient <- rownames(clinical_purity)

ids <- rownames(clinical_purity)
rownames(bp_estimates_malignant) <- ids
rownames(bisque_estimates_malignant) <- ids
bp_estimates_malignant$Patient <- rownames(bp_estimates_malignant)
bisque_estimates_malignant$Patient <- rownames(bisque_estimates_malignant)

colnames(clinical_purity) <- c("TumorPurity", "Patient")
colnames(bp_estimates_malignant) <- c("MalignantFraction", "Patient")
colnames(bisque_estimates_malignant) <- c("MalignantFraction", "Patient")

clinical_purity <- clinical_purity %>% replace_na(list(TumorPurity = 0))

bayes_purity <- merge(clinical_purity, bp_estimates_malignant, by = "Patient")
bisque_purity <- merge(clinical_purity, bisque_estimates_malignant, by = "Patient")

bayes_purity <- bayes_purity %>% filter(TumorPurity != 0)
bisque_purity <- bisque_purity %>% filter(TumorPurity != 0)

p10<-ggplot(bayes_purity, aes(x = TumorPurity, y = MalignantFraction)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add a linear regression line
  labs(title = "Correlation between Tumor Purity and Malignant Cell Fraction - BayesPrism",
       x = "Tumor Purity",
       y = "Malignant Cell Fraction") +
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 0.05, label.y = 0.95, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")))

p11<-ggplot(bisque_purity, aes(x = TumorPurity, y = MalignantFraction)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add a linear regression line
  labs(title = "Correlation between Tumor Purity and Malignant Cell Fraction - BisqueRNA",
       x = "Tumor Purity",
       y = "Malignant Cell Fraction") +
  theme_minimal() +
  stat_cor(method = "pearson", label.x = 0.05, label.y = 0.95, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")))


#Comparison with clinical and other molecular features present in metadata. Now for the results from the top 25 markers

clinical25 <- clinical_full[,c("paper_ABSOLUTE.Purity","barcode","vital_status","days_to_death","days_to_last_follow_up","paper_Ploidy",
                               "paper_Copy.Number.Clusters..All.150.Samples.","paper_Pathologist.Reviewed.Tumor.Cellularity")]
clinical_purity25 <- clinical25[,c("paper_ABSOLUTE.Purity","barcode","paper_Ploidy","paper_Copy.Number.Clusters..All.150.Samples.","paper_Pathologist.Reviewed.Tumor.Cellularity")]
clinical_purity25$barcode <- NULL
# Sum the proportions for the new classes and remove the original columns
bisque_estimates25$Malignant <- bisque_estimates25$Malignant_high_cnv + bisque_estimates25$Malignant_low_cnv
bisque_estimates25 <- bisque_estimates25[, !(names(bisque_estimates25) %in% c("Malignant_high_cnv", "Malignant_low_cnv"))]

bp_estimates25$Malignant <- bp_estimates25$Malignant_high_cnv + bp_estimates25$Malignant_low_cnv
bp_estimates25 <- bp_estimates25[, !(names(bp_estimates25) %in% c("Malignant_high_cnv", "Malignant_low_cnv"))]

bisque_estimates_malignant25 <- as.data.frame(bisque_estimates25[,"Malignant"])
bp_estimates_malignant25 <- as.data.frame(bp_estimates25[,"Malignant"])

clinical_purity25$Patient <- rownames(clinical_purity25)

ids <- rownames(clinical_purity25)
rownames(bp_estimates_malignant25) <- ids
rownames(bisque_estimates_malignant25) <- ids
bp_estimates_malignant25$Patient <- rownames(bp_estimates_malignant25)
bisque_estimates_malignant25$Patient <- rownames(bisque_estimates_malignant25)

colnames(clinical_purity25) <- c("TumorPurity", "Tumor_Ploidy_score","CNV_classification","Pathologist_tumor_cellularity","Patient")
colnames(bp_estimates_malignant25) <- c("MalignantFraction","Patient")
colnames(bisque_estimates_malignant25) <- c("MalignantFraction","Patient")

clinical_purity25 <- clinical_purity25 %>% replace_na(list(TumorPurity = 0))

bayes_purity25 <- merge(clinical_purity25, bp_estimates_malignant25, by = "Patient")
bisque_purity25 <- merge(clinical_purity25, bisque_estimates_malignant25, by = "Patient")

bisque_purity25 <- bisque_purity25 %>% filter(TumorPurity != 0)
bayes_purity25 <- bayes_purity25 %>% filter(TumorPurity != 0)

p12<-ggplot(bayes_purity25, aes(x = TumorPurity, y = MalignantFraction)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add a linear regression line
  labs(title = "Correlation between Tumor Purity and Malignant Cell Fraction - BayesPrism (Top25 markers)",
       x = "Tumor Purity",
       y = "Malignant Cell Fraction") +
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 0.05, label.y = 0.95, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")))

p13<-ggplot(bisque_purity25, aes(x = TumorPurity, y = MalignantFraction)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add a linear regression line
  labs(title = "Correlation between Tumor Purity and Malignant Cell Fraction - BisqueRNA (Top25 markers)",
       x = "Tumor Purity",
       y = "Malignant Cell Fraction") +
  theme_minimal() +
  stat_cor(method = "pearson", label.x = 0.05, label.y = 0.95, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")))

#Compare High and Low CNV groups (Bayes and Prism)
bayes_malignant <- bp_estimates25
df_long <- df %>%
  pivot_longer(cols = c(Malignant_high_cnv, Malignant_low_cnv), 
               names_to = "cnv_type", 
               values_to = "fraction")

p <- ggplot(df_long, aes(x = group, y = fraction, fill = cnv_type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Proportion of Malignant CNVs by Group",
       x = "Group",
       y = "Proportion",
       fill = "CNV Type") +
  theme_minimal()

# Display the plot
print(p)