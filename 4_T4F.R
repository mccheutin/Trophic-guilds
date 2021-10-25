## T4F -----
setwd( "/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Carnivory")
dir.create("./Functions", recursive = T)
setwd("./Functions")

# Setting T4F ----
library(Tax4Fun2)

# Build the reference data and dependencies
# library(Tax4Fun2)
buildReferenceData(path_to_working_directory = ".", install_suggested_packages = T)
buildDependencies(path_to_reference_data = "Tax4Fun2_ReferenceData_v2")

# Predict functions for the gut core ----
load("../physeq_object/gut_core.RData")
otu_table <-as.data.frame(otu_table(gut_core))
write.table(otu_table, sep = '\t', file = 'otu_table.txt', quote = F, row.names = T, col.names = T )
# Add 'ID' as colnames for the ASVs rows
tax_table <- as.data.frame(tax_table(gut_core))
write.table(tax_table, sep = '\t', 'tax_table.txt', quote = F)
Biostrings::writeXStringSet(gut_core@refseq, file = "gut_core.fasta")

runRefBlast(path_to_otus = "gut_core.fasta", path_to_reference_data = "Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "core_Ref99NR", database_mode = "Ref99NR", use_force = T, num_threads = 1)
makeFunctionalPrediction(path_to_otu_table = "otu_table.txt", 
                         path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
                         path_to_temp_folder = "core_Ref99NR", 
                         database_mode = "Ref99NR", 
                         normalize_by_copy_number = TRUE, 
                         min_identity_to_reference = 0.90, #min_identity_to_reference = round(min(read.table("./core_Ref99NR/ref_blast.txt")[,3]),0),
                         normalize_pathways = T)

calculateFunctionalRedundancy(path_to_otu_table = "otu_table.txt",
                              path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
                              path_to_temp_folder = "core_Ref99NR",
                              database_mode = "Ref99NR",
                              min_identity_to_reference = 0.97)



prediction_core <- as.data.frame(read_delim("core_Ref99NR/functional_prediction.txt", "\t", escape_double = FALSE, trim_ws = TRUE))
colnames(prediction_core)[1] <- "X.OTU_IDs"
ko_pathway_sort <- read_delim("ko_pathway_sort.txt","\t", escape_double = FALSE, trim_ws = TRUE)
prediction_core$X.OTU_IDs <- str_replace_all(prediction_core$X.OTU_IDs, "K", "ko")
write.table(prediction_core , file = "prediction_core.txt", sep = "\t", quote = F, row.names = F)
prediction_core <- read.table(file = "prediction_core.txt", sep = "\t" , header = T)

prediction_core2 <- left_join(prediction_core, ko_pathway_sort, by = "X.OTU_IDs")

# How many KO redudant?
length(which(table(prediction_core2$X.OTU_IDs)>1, T)) # 1771
redundant_ko <- rownames(which(table(prediction_core2$X.OTU_IDs)>1, T))

prediction_core3 <- prediction_core2[-which(prediction_core2$X.OTU_IDs %in% redundant_ko), ]
prediction_rel <- prediction_core3
prediction_rel[,2:379] <- prediction_rel[,2:379]/sum(prediction_rel[,2:379])
save(prediction_core, prediction_core2, prediction_core3, prediction_rel, file = "prediction_files.RData")

db_core <- as.data.frame(sample_data(gut_core))
rownames(db_core) <- str_replace_all(rownames(db_core), "-", ".")
save(db_core, file = "db_core.RData")
######## ***********************************####################
######## *****  LEFSE analysis ****** #####
######## ***********************************####################
#  LEFSE on the three dietary levels -----
dir.create("./lefse", recursive = T)
setwd("./lefse")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
### __ For KO ----#########################
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
#
lefse <- data.frame(db_core$diet3, stringsAsFactors = FALSE)
ko_core = prediction_rel
rownames(ko_core) <- ko_core$X.OTU_IDs
ko_core <- ko_core[,-1]
ko_core <- t(ko_core)

lefse <- cbind(lefse, ko_core[1:378,])
KO_table_core <- prediction_rel$X.OTU_IDs
tax_name_core <- KO_table_core
tax_name_core <- c("Diet",tax_name_core)
lefse2 <- cbind(tax_name_core,t(lefse))
write.table(lefse2, file="lefse_diet_pred_KO.txt",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


# Run your file in LEfSE software on Galaxy, download the results, move it on your active folder "dir_data_cleaning" and rename the file "LDA_results".
LDA_Effect_Size <- read.delim("LDA_results_KO", header=FALSE)
LDA_carn <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Carnivorous")
LDA_herb <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Herbivorous")
LDA_omni <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Omnivorous")

LDA_all <- rbind(LDA_carn,LDA_herb,LDA_omni)

colnames(LDA_all) <- c("X.OTU_IDs" , "LDA_res" , "Diet" , "LDA_res2" , "sig")
write.table(LDA_all, file = "LDA_all_ko.txt", sep = "\t", quote = T,row.names = F )

ko_lda <- LDA_all$X.OTU_IDs
prediction_lda <- prediction_rel[which(prediction_rel$X.OTU_IDs %in% ko_lda,T),]
write.table(prediction_lda , file = "prediction_table.txt", sep = "\t", row.names = F)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
### __ For KEGG description ----#########################
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 

lefse <- data.frame(db_core$diet3, stringsAsFactors = FALSE)
ko_core = aggregate(prediction_rel[,2:379], by=list(description = prediction_rel$description), FUN=sum)
rownames(ko_core) <- ko_core$description
ko_core <- ko_core[,-1]
ko_core <- t(ko_core)

lefse <- cbind(lefse, ko_core)
description_table_core <- colnames(ko_core)
tax_name_core <- description_table_core
tax_name_core <- c("Diet",tax_name_core)
lefse2 <- cbind(tax_name_core,t(lefse))
write.table(lefse2, file="lefse_diet_pred_description.txt",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


# Run your file in LEfSE software on Galaxy, download the results, move it on your active folder "dir_data_cleaning" and rename the file "LDA_results".
LDA_Effect_Size <- read.delim("LDA_results_description", header=FALSE)
LDA_carn <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Carnivorous")
LDA_herb <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Herbivorous")
LDA_omni <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Omnivorous")

LDA_all <- rbind(LDA_carn,LDA_herb,LDA_omni)

colnames(LDA_all) <- c("Kegg" , "LDA_res" , "Diet" , "LDA_res2" , "sig")
write.table(LDA_all, file = "LDA_all_description.txt", sep = "\t", quote = T,row.names = F )

kegg_lda <- LDA_all$Kegg
write.table(kegg_lda, file = "kegg_lda.txt", sep = "\t", row.names= F) # To edit the name of the descriptions which have been modificated during Lefse import
kegg_lda <- read.table(file = "kegg_lda.txt",  header = T, sep = "\t")

prediction_lda_kegg <- prediction_rel[which(prediction_rel$description %in% kegg_lda$description,T),]
write.table(prediction_lda_kegg , file = "prediction_table_kegg.txt", sep = "\t", row.names = F)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
# __ For level C description----
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
#
lefse <- data.frame(db_core$diet3, stringsAsFactors = FALSE)
ko_core = aggregate(prediction_rel[,2:379], by=list(level_C = prediction_rel$C), FUN=sum)
rownames(ko_core) <- ko_core$level_C
ko_core <- ko_core[,-1]
ko_core <- t(ko_core)

lefse <- cbind(lefse, ko_core)
C_table_core <- colnames(ko_core)
tax_name_core <- C_table_core
tax_name_core <- c("Diet",tax_name_core)
lefse2 <- cbind(tax_name_core,t(lefse))
write.table(lefse2, file="lefse_diet_pred_C.txt",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


# Run your file in LEfSE software on Galaxy, download the results, move it on your active folder "dir_data_cleaning" and rename the file "LDA_results".
LDA_Effect_Size <- read.delim("LDA_results_C", header=FALSE)
LDA_carn <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Carnivorous")
LDA_herb <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Herbivorous")
LDA_omni <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Omnivorous")

LDA_all <- rbind(LDA_carn,LDA_herb,LDA_omni)

colnames(LDA_all) <- c("Level_C" , "LDA_res" , "Diet" , "LDA_res2" , "sig")
write.table(LDA_all, file = "LDA_all_C.txt", sep = "\t", quote = T,row.names = F )

C_lda <- LDA_all$Level_C
write.table(C_lda, file = "C_lda.txt", sep = "\t", row.names= F) # To edit the name of the descriptions which have been modificated during Lefse import
C_lda <- read.table(file = "C_lda.txt",  header = T, sep = "\t")

prediction_lda_C <- prediction_rel[which(prediction_rel$C %in% C_lda$description,T),]
write.table(prediction_lda_C , file = "prediction_table_C.txt", sep = "\t", row.names = F)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
# __ For level B description----
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
lefse <- data.frame(db_core$diet3, stringsAsFactors = FALSE)
ko_core = aggregate(prediction_rel[,2:379], by=list(level_B = prediction_rel$B), FUN=sum)
rownames(ko_core) <- ko_core$level_B
ko_core <- ko_core[,-1]
ko_core <- t(ko_core)

lefse <- cbind(lefse, ko_core)
B_table_core <- colnames(ko_core)
tax_name_core <- B_table_core
tax_name_core <- c("Diet",tax_name_core)
lefse2 <- cbind(tax_name_core,t(lefse))
write.table(lefse2, file="lefse_diet_pred_B.txt",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


# Run your file in LEfSE software on Galaxy, download the results, move it on your active folder "dir_data_cleaning" and rename the file "LDA_results".
LDA_Effect_Size <- read.delim("LDA_results_B", header=FALSE)
LDA_carn <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Carnivorous")
LDA_herb <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Herbivorous")
LDA_omni <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Omnivorous")

LDA_all <- rbind(LDA_carn,LDA_herb,LDA_omni)

colnames(LDA_all) <- c("Level_B" , "LDA_res" , "Diet" , "LDA_res2" , "sig")
write.table(LDA_all, file = "LDA_all_B.txt", sep = "\t", quote = T,row.names = F )
B_lda <- LDA_all$Level_B
write.table(B_lda , file = "B_lda.txt", sep = "\t", quote = F, row.names = F)
B_lda <- read.table(file = "B_lda.txt", sep = "\t", header = T)

prediction_lda_B <- prediction_rel[which(prediction_rel$B %in% B_lda$description,T),]
write.table(prediction_lda_B , file = "prediction_table_B.txt", sep = "\t", row.names = F)

## Plot the Lefse results in boxplot ----
## __ For the KO ----
# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(prediction_lda[,2:379])))
# 2 - Create the database with the Ko description
ko_db <- cbind("X.OTU_IDs" = rep(prediction_lda$X.OTU_IDs, each = dim(prediction_lda[,2:379])[2]),
                 "description" = rep(prediction_lda$description, each = dim(prediction_lda[,2:379])[2]),
                 "A" = rep(prediction_lda$A, each = dim(prediction_lda[,2:379])[2]),
                 "B" = rep(prediction_lda$B, each = dim(prediction_lda[,2:379])[2]),
                 "C" = rep(prediction_lda$C, each = dim(prediction_lda[,2:379])[2]),
                 "kegg" = rep(prediction_lda$kegg, each = dim(prediction_lda[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(prediction_lda[2:379]) , dim(prediction_lda[,2:379])[1])
# 4 - Create the dataframe
ko_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, ko_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(prediction_lda[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(prediction_lda[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(prediction_lda[,2:379])[1]))))
# 6 - Merge the database
ko_df_2 <- as.data.frame(cbind(ko_df , diet_df[,c(2,3)]))
ko_df_2$value <- as.numeric(as.character(ko_df_2$value))
save(ko_df_2, file = "ko_df_2.RData")

# By diet predicted in Parraviccini et al., 2020
ko_plot_diet_pred  <- ggplot(ko_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ X.OTU_IDs,  scales = "free")

ko_plot_diet3  <- ggplot(ko_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ X.OTU_IDs,  scales = "free")

pdf(file = "KO_plots.pdf" , he = 7 , wi= 12)
ko_plot_diet_pred
ko_plot_diet3
dev.off()

## __ For the Kegg pathways ----

# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(prediction_lda_kegg[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(prediction_lda_kegg$X.OTU_IDs, each = dim(prediction_lda_kegg[,2:379])[2]),
               "description" = rep(prediction_lda_kegg$description, each = dim(prediction_lda_kegg[,2:379])[2]),
               "A" = rep(prediction_lda_kegg$A, each = dim(prediction_lda_kegg[,2:379])[2]),
               "B" = rep(prediction_lda_kegg$B, each = dim(prediction_lda_kegg[,2:379])[2]),
               "C" = rep(prediction_lda_kegg$C, each = dim(prediction_lda_kegg[,2:379])[2]),
               "kegg" = rep(prediction_lda_kegg$kegg, each = dim(prediction_lda_kegg[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(prediction_lda_kegg[2:379]) , dim(prediction_lda_kegg[,2:379])[1])
# 4 - Create the dataframe
kegg_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(prediction_lda_kegg[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(prediction_lda_kegg[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(prediction_lda_kegg[,2:379])[1]))))
# 6 - Merge the database
kegg_df_2 <- as.data.frame(cbind(kegg_df , diet_df[,c(2,3)]))
kegg_df_2$value <- as.numeric(as.character(kegg_df_2$value))
save(kegg_df_2, file = "kegg_df_2.RData")

# By diet predicted in Parraviccini et al., 2020
kegg_plot_diet_pred  <- ggplot(kegg_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ description,  scales = "free")

kegg_plot_diet3  <- ggplot(kegg_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ description,  scales = "free")

pdf(file = "kegg_plots.pdf" , he = 7 , wi= 15)
kegg_plot_diet_pred
kegg_plot_diet3
dev.off()

## __ For the C level ----

# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(prediction_lda_C[,2:379])))
# 2 - Create the database with the Ko description
C_db <- cbind("X.OTU_IDs" = rep(prediction_lda_C$X.OTU_IDs, each = dim(prediction_lda_C[,2:379])[2]),
                 "description" = rep(prediction_lda_C$description, each = dim(prediction_lda_C[,2:379])[2]),
                 "A" = rep(prediction_lda_C$A, each = dim(prediction_lda_C[,2:379])[2]),
                 "B" = rep(prediction_lda_C$B, each = dim(prediction_lda_C[,2:379])[2]),
                 "C" = rep(prediction_lda_C$C, each = dim(prediction_lda_C[,2:379])[2]),
                 "kegg" = rep(prediction_lda_C$kegg, each = dim(prediction_lda_C[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(prediction_lda_C[2:379]) , dim(prediction_lda_C[,2:379])[1])
# 4 - Create the dataframe
C_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, C_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(prediction_lda_C[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(prediction_lda_C[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(prediction_lda_C[,2:379])[1]))))
# 6 - Merge the database
C_df_2 <- as.data.frame(cbind(C_df , diet_df[,c(2,3)]))
C_df_2$value <- as.numeric(as.character(C_df_2$value))
save(C_df_2, file = "C_df_2.RData")
# By diet predicted in Parraviccini et al., 2020
C_plot_diet_pred  <- ggplot(C_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ C,  scales = "free")

C_plot_diet3  <- ggplot(C_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ C,  scales = "free")

pdf(file = "C_plots.pdf" , he = 7 , wi= 10)
C_plot_diet_pred
C_plot_diet3
dev.off()

## __ For the B level ----

# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(prediction_lda_B[,2:379])))
# 2 - Create the database with the Ko description
B_db <- cbind("X.OTU_IDs" = rep(prediction_lda_B$X.OTU_IDs, each = dim(prediction_lda_B[,2:379])[2]),
                 "description" = rep(prediction_lda_B$description, each = dim(prediction_lda_B[,2:379])[2]),
                 "A" = rep(prediction_lda_B$A, each = dim(prediction_lda_B[,2:379])[2]),
                 "B" = rep(prediction_lda_B$B, each = dim(prediction_lda_B[,2:379])[2]),
                 "C" = rep(prediction_lda_B$C, each = dim(prediction_lda_B[,2:379])[2]),
                 "kegg" = rep(prediction_lda_B$kegg, each = dim(prediction_lda_B[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(prediction_lda_B[2:379]) , dim(prediction_lda_B[,2:379])[1])
# 4 - Create the dataframe
B_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, B_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(prediction_lda_B[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(prediction_lda_B[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(prediction_lda_B[,2:379])[1]))))
# 6 - Merge the database
B_df_2 <- as.data.frame(cbind(B_df , diet_df[,c(2,3)]))
B_df_2$value <- as.numeric(as.character(B_df_2$value))
save(B_df_2, file = "B_df_2.RData")

# By diet predicted in Parraviccini et al., 2020
B_plot_diet_pred  <- ggplot(B_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ B,  scales = "free")

B_plot_diet3  <- ggplot(B_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ B,  scales = "free")

pdf(file = "B_plots.pdf" , he = 7 , wi= 10)
B_plot_diet_pred
B_plot_diet3
dev.off()

######## ***********************************####################
########  *** ALL FUNCTIONS ***** -----
######## ***********************************####################
setwd("../")

## __ Urate associated functions ----
urate_gene <- prediction_rel[rownames(prediction_rel)[which(str_detect(prediction_rel$description , "urate")
                                                                | str_detect(prediction_rel$description, "ureid")
                                                                | str_detect(prediction_rel$description, "isourate")
                                                                | str_detect(prediction_rel$description, "allanto")
                                                                | str_detect(prediction_rel$kegg , "urate")
                                                                | str_detect(prediction_rel$kegg, "ureid")
                                                                | str_detect(prediction_rel$kegg, "isourate")
                                                                | str_detect(prediction_rel$kegg, "allanto")
)],]$description

urate_db <- prediction_rel %>% filter(description %in% urate_gene)
# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(urate_db[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(urate_db$X.OTU_IDs, each = dim(urate_db[,2:379])[2]),
                 "description" = rep(urate_db$description, each = dim(urate_db[,2:379])[2]),
                 "A" = rep(urate_db$A, each = dim(urate_db[,2:379])[2]),
                 "B" = rep(urate_db$B, each = dim(urate_db[,2:379])[2]),
                 "C" = rep(urate_db$C, each = dim(urate_db[,2:379])[2]),
                 "kegg" = rep(urate_db$kegg, each = dim(urate_db[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(urate_db[2:379]) , dim(urate_db[,2:379])[1])
# 4 - Create the dataframe
urate_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(urate_db[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(urate_db[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(urate_db[,2:379])[1]))))
# 6 - Merge the database
urate_df_2 <- as.data.frame(cbind(urate_df , diet_df[,c(2,3)]))
urate_df_2$value <- as.numeric(as.character(urate_df_2$value))
save(urate_df_2 , file = "urate_df_2.RData")

# __ * Plot ----
ureic_plot_diet_pred  <- ggplot(urate_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ description,  scales = "free")

ureic_plot_diet3  <- ggplot(urate_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  geom_signif(comparisons = list(c("Carnivorous", "Herbivorous"),c("Carnivorous","Omnivorous"),c("Herbivorous", "Omnivorous")), 
              map_signif_level = TRUE, textsize=4, color="black", family = "serif", vjust = 0.5)+
  facet_wrap( ~ description,  scales = "free")

pdf(file = "./plots/ureic_plot.pdf" , he = 7 , wi= 15)
ureic_plot_diet_pred
ureic_plot_diet3
dev.off()

allantoinase <- urate_df_2 %>% 
  filter(description == "allantoinase [EC:3.5.2.5]") %>%
  group_by(diet3) %>%
  summarize(Mean = mean(value, na.rm=TRUE)*100*10^6)


## __ Chitin associated functions ----
chitinase_gene <- prediction_rel[rownames(prediction_rel)[which(str_detect(prediction_rel$description , "chitin")
                                                                | str_detect(prediction_rel$description, "endo-beta-N-acetylglucosaminidase"))],]$description

chitin_db <- prediction_rel %>% filter(description %in% chitinase_gene)

# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(chitin_db[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(chitin_db$X.OTU_IDs, each = dim(chitin_db[,2:379])[2]),
                 "description" = rep(chitin_db$description, each = dim(chitin_db[,2:379])[2]),
                 "A" = rep(chitin_db$A, each = dim(chitin_db[,2:379])[2]),
                 "B" = rep(chitin_db$B, each = dim(chitin_db[,2:379])[2]),
                 "C" = rep(chitin_db$C, each = dim(chitin_db[,2:379])[2]),
                 "kegg" = rep(chitin_db$kegg, each = dim(chitin_db[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(chitin_db[2:379]) , dim(chitin_db[,2:379])[1])
# 4 - Create the dataframe
chitin_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(chitin_db[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(chitin_db[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(chitin_db[,2:379])[1]))))
# 6 - Merge the database
chitin_df_2 <- as.data.frame(cbind(chitin_df , diet_df[,c(2,3)]))
chitin_df_2$value <- as.numeric(as.character(chitin_df_2$value))
save(chitin_df_2 , file = "chitin_df_2.RData")

# __ * Plot ----
chitin_plot_diet_pred  <- ggplot(chitin_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ description,  scales = "free")

chitin_plot_diet3  <- ggplot(chitin_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  geom_signif(comparisons = list(c("Carnivorous", "Herbivorous"),c("Carnivorous","Omnivorous"),c("Herbivorous", "Omnivorous")), 
              map_signif_level = TRUE, textsize=4, color="black", family = "serif", vjust = 0.5)+
  facet_wrap( ~ description,  scales = "free")

pdf(file = "chitin_plot.pdf" , he = 7 , wi= 15)
chitin_plot_diet_pred
chitin_plot_diet3
dev.off()

## __ Lipase associated functions ----
lipase_gene <- prediction_rel[rownames(prediction_rel)[which(str_detect(prediction_rel$description , "lipase")
                                                               # | str_detect(prediction_rel$description, "endo-beta-N-acetylglucosaminidase")
                                                               )],]$description

lipase_db <- prediction_rel %>% filter(description %in% lipase_gene)

# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(lipase_db[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(lipase_db$X.OTU_IDs, each = dim(lipase_db[,2:379])[2]),
                 "description" = rep(lipase_db$description, each = dim(lipase_db[,2:379])[2]),
                 "A" = rep(lipase_db$A, each = dim(lipase_db[,2:379])[2]),
                 "B" = rep(lipase_db$B, each = dim(lipase_db[,2:379])[2]),
                 "C" = rep(lipase_db$C, each = dim(lipase_db[,2:379])[2]),
                 "kegg" = rep(lipase_db$kegg, each = dim(lipase_db[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(lipase_db[2:379]) , dim(lipase_db[,2:379])[1])
# 4 - Create the dataframe
lipase_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(lipase_db[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(lipase_db[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(lipase_db[,2:379])[1]))))
# 6 - Merge the database
lipase_df_2 <- as.data.frame(cbind(lipase_df , diet_df[,c(2,3)]))
lipase_df_2$value <- as.numeric(as.character(lipase_df_2$value))
save(lipase_df_2 , file = "lipase_df_2.RData")

# __ * Plot ----
lipase_plot_diet_pred  <- ggplot(lipase_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ description,  scales = "free")

lipase_plot_diet3  <- ggplot(lipase_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  geom_signif(comparisons = list(c("Carnivorous", "Herbivorous"),c("Carnivorous","Omnivorous"),c("Herbivorous", "Omnivorous")), 
              map_signif_level = TRUE, textsize=4, color="black", family = "serif", vjust = 0.5)+
  facet_wrap( ~ description,  scales = "free")

pdf(file = "./plots/lipase_plot.pdf" , he = 7 , wi= 15)
lipase_plot_diet3
lipase_plot_diet_pred
dev.off()

## __ Carnitine associated functions ----
carnitine_gene <- prediction_rel[rownames(prediction_rel)[which(str_detect(prediction_rel$description , "carnitin")
                                                             # | str_detect(prediction_rel$description, "endo-beta-N-acetylglucosaminidase")
)],]$description

carnitine_db <- prediction_rel %>% filter(description %in% carnitine_gene)

# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(carnitine_db[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(carnitine_db$X.OTU_IDs, each = dim(carnitine_db[,2:379])[2]),
                 "description" = rep(carnitine_db$description, each = dim(carnitine_db[,2:379])[2]),
                 "A" = rep(carnitine_db$A, each = dim(carnitine_db[,2:379])[2]),
                 "B" = rep(carnitine_db$B, each = dim(carnitine_db[,2:379])[2]),
                 "C" = rep(carnitine_db$C, each = dim(carnitine_db[,2:379])[2]),
                 "kegg" = rep(carnitine_db$kegg, each = dim(carnitine_db[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(carnitine_db[2:379]) , dim(carnitine_db[,2:379])[1])
# 4 - Create the dataframe
carnitine_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(carnitine_db[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(carnitine_db[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(carnitine_db[,2:379])[1]))))
# 6 - Merge the database
carnitine_df_2 <- as.data.frame(cbind(carnitine_df , diet_df[,c(2,3)]))
carnitine_df_2$value <- as.numeric(as.character(carnitine_df_2$value))
save(carnitine_df_2 , file = "carnitine_df_2.RData")

# __ * Plot ----
carnitine_plot_diet_pred  <- ggplot(carnitine_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ description,  scales = "free")

carnitine_plot_diet3  <- ggplot(carnitine_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  geom_signif(comparisons = list(c("Carnivorous", "Herbivorous"),c("Carnivorous","Omnivorous"),c("Herbivorous", "Omnivorous")), 
              map_signif_level = TRUE, textsize=4, color="black", family = "serif", vjust = 0.5)+
  facet_wrap( ~ description,  scales = "free")

pdf(file = "carnitine_plot.pdf" , he = 7 , wi= 15)
carnitine_plot_diet_pred
carnitine_plot_diet3
dev.off()

## __ Cellulase associated functions ----
cellulase_gene <- prediction_rel[rownames(prediction_rel)[which(str_detect(prediction_rel$description , "cellulose")
                                                                | str_detect(prediction_rel$description , "cellobios")
                                                                |str_detect(prediction_rel$description ,  "beta-glucosidase")
                                                                |str_detect(prediction_rel$description ,  "arabinofuranosidase")
                                                                |str_detect(prediction_rel$description ,  "pect")
                                                                |str_detect(prediction_rel$description ,  "endoglucanase")
)],]$description

cellulase_db <- prediction_rel %>% filter(description %in% cellulase_gene)

# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(cellulase_db[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(cellulase_db$X.OTU_IDs, each = dim(cellulase_db[,2:379])[2]),
                 "description" = rep(cellulase_db$description, each = dim(cellulase_db[,2:379])[2]),
                 "A" = rep(cellulase_db$A, each = dim(cellulase_db[,2:379])[2]),
                 "B" = rep(cellulase_db$B, each = dim(cellulase_db[,2:379])[2]),
                 "C" = rep(cellulase_db$C, each = dim(cellulase_db[,2:379])[2]),
                 "kegg" = rep(cellulase_db$kegg, each = dim(cellulase_db[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(cellulase_db[2:379]) , dim(cellulase_db[,2:379])[1])
# 4 - Create the dataframe
cellulase_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(cellulase_db[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(cellulase_db[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(cellulase_db[,2:379])[1]))))
# 6 - Merge the database
cellulase_df_2 <- as.data.frame(cbind(cellulase_df , diet_df[,c(2,3)]))
cellulase_df_2$value <- as.numeric(as.character(cellulase_df_2$value))
save(cellulase_df_2 , file = "cellulase_df_2.RData")

# __ * Plot ----
cellulase_plot_diet_pred  <- ggplot(cellulase_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ description,  scales = "free")

cellulase_plot_diet3  <- ggplot(cellulase_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  geom_signif(comparisons = list(c("Carnivorous", "Herbivorous"),c("Carnivorous","Omnivorous"),c("Herbivorous", "Omnivorous")), 
              map_signif_level = TRUE, textsize=4, color="black", family = "serif", vjust = 0.5)+
  facet_wrap( ~ description,  scales = "free")

pdf(file = "cellulase_plot.pdf" , he = 7 , wi= 15)
cellulase_plot_diet_pred
cellulase_plot_diet3
dev.off()

## __ Amylase associated functions ----
amylase_gene <- prediction_rel[rownames(prediction_rel)[which(str_detect(prediction_rel$description , "amyla")
)],]$description

amylase_db <- prediction_rel %>% filter(description %in% amylase_gene)

# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(amylase_db[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(amylase_db$X.OTU_IDs, each = dim(amylase_db[,2:379])[2]),
                 "description" = rep(amylase_db$description, each = dim(amylase_db[,2:379])[2]),
                 "A" = rep(amylase_db$A, each = dim(amylase_db[,2:379])[2]),
                 "B" = rep(amylase_db$B, each = dim(amylase_db[,2:379])[2]),
                 "C" = rep(amylase_db$C, each = dim(amylase_db[,2:379])[2]),
                 "kegg" = rep(amylase_db$kegg, each = dim(amylase_db[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(amylase_db[2:379]) , dim(amylase_db[,2:379])[1])
# 4 - Create the dataframe
amylase_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(amylase_db[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(amylase_db[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(amylase_db[,2:379])[1]))))
# 6 - Merge the database
amylase_df_2 <- as.data.frame(cbind(amylase_df , diet_df[,c(2,3)]))
amylase_df_2$value <- as.numeric(as.character(amylase_df_2$value))
save(amylase_df_2 , file = "amylase_df_2.RData")

# __ * Plot ----
amylase_plot_diet_pred  <- ggplot(amylase_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ description,  scales = "free")

amylase_plot_diet3  <- ggplot(amylase_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  geom_signif(comparisons = list(c("Carnivorous", "Herbivorous"),c("Carnivorous","Omnivorous"),c("Herbivorous", "Omnivorous")), 
              map_signif_level = TRUE, textsize=4, color="black", family = "serif", vjust = 0.5)+
  facet_wrap( ~ description,  scales = "free")

pdf(file = "amylase_plot.pdf" , he = 7 , wi= 15)
amylase_plot_diet_pred
amylase_plot_diet3
dev.off()

## __ Glutamine/Glutamate metabolism associated functions ----
glut_gene <- prediction_rel[rownames(prediction_rel)[which(str_detect(prediction_rel$description , "glutamin")
                                                               |str_detect(prediction_rel$description , "glutamat")
)],]$description

glut_db <- prediction_rel %>% filter(description %in% glut_gene)
# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(glut_db[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(glut_db$X.OTU_IDs, each = dim(glut_db[,2:379])[2]),
                 "description" = rep(glut_db$description, each = dim(glut_db[,2:379])[2]),
                 "A" = rep(glut_db$A, each = dim(glut_db[,2:379])[2]),
                 "B" = rep(glut_db$B, each = dim(glut_db[,2:379])[2]),
                 "C" = rep(glut_db$C, each = dim(glut_db[,2:379])[2]),
                 "kegg" = rep(glut_db$kegg, each = dim(glut_db[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(glut_db[2:379]) , dim(glut_db[,2:379])[1])
# 4 - Create the dataframe
glut_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(glut_db[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(glut_db[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(glut_db[,2:379])[1]))))
# 6 - Merge the database
glut_df_2 <- as.data.frame(cbind(glut_df , diet_df[,c(2,3)]))
glut_df_2$value <- as.numeric(as.character(glut_df_2$value))
save(glut_df_2 , file = "glut_df_2.RData")

# __ * Plot ----
glut_plot_diet_pred  <- ggplot(glut_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ C,  scales = "free")

glut_plot_diet3  <- ggplot(glut_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  geom_signif(comparisons = list(c("Carnivorous", "Herbivorous"),c("Carnivorous","Omnivorous"),c("Herbivorous", "Omnivorous")), 
              map_signif_level = TRUE, textsize=4, color="black", family = "serif", vjust = 0.5)+
  facet_wrap( ~ C,  scales = "free")

pdf(file = "glut_plot.pdf" , he = 7 , wi= 15)
glut_plot_diet_pred
glut_plot_diet3
dev.off()

## __ GS, GDH, ALT metabolism associated functions ----
gs_gdh_alt_gene <- prediction_rel[rownames(prediction_rel)[which(str_detect(prediction_rel$description , "glutamine synthetase")
                                             |str_detect(prediction_rel$description , "glutamate dehydrogenase")
                                             |str_detect(prediction_rel$description , "alanine-synthesizing transaminase")
)],]$description

gs_gdh_alt_db <- prediction_rel %>% filter(description %in% gs_gdh_alt_gene)
# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(gs_gdh_alt_db[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(gs_gdh_alt_db$X.OTU_IDs, each = dim(gs_gdh_alt_db[,2:379])[2]),
                 "description" = rep(gs_gdh_alt_db$description, each = dim(gs_gdh_alt_db[,2:379])[2]),
                 "A" = rep(gs_gdh_alt_db$A, each = dim(gs_gdh_alt_db[,2:379])[2]),
                 "B" = rep(gs_gdh_alt_db$B, each = dim(gs_gdh_alt_db[,2:379])[2]),
                 "C" = rep(gs_gdh_alt_db$C, each = dim(gs_gdh_alt_db[,2:379])[2]),
                 "kegg" = rep(gs_gdh_alt_db$kegg, each = dim(gs_gdh_alt_db[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(gs_gdh_alt_db[2:379]) , dim(gs_gdh_alt_db[,2:379])[1])
# 4 - Create the dataframe
gs_gdh_alt_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(gs_gdh_alt_db[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(gs_gdh_alt_db[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(gs_gdh_alt_db[,2:379])[1]))))
# 6 - Merge the database
gs_gdh_alt_df_2 <- as.data.frame(cbind(gs_gdh_alt_df , diet_df[,c(2,3)]))
gs_gdh_alt_df_2$value <- as.numeric(as.character(gs_gdh_alt_df_2$value))
save(gs_gdh_alt_df_2 , file = "gs_gdh_alt_df_2.RData")

# __ * Plot ----
gs_gdh_alt_plot_diet_pred  <- ggplot(gs_gdh_alt_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ description,  scales = "free")

gs_gdh_alt_plot_diet3  <- ggplot(gs_gdh_alt_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  geom_signif(comparisons = list(c("Carnivorous", "Herbivorous"),c("Carnivorous","Omnivorous"),c("Herbivorous", "Omnivorous")), 
              map_signif_level = TRUE, textsize=4, color="black", family = "serif", vjust = 0.5)+
  facet_wrap( ~ description,  scales = "free")

pdf(file = "gs_gdh_alt_plot.pdf" , he = 7 , wi= 15)
gs_gdh_alt_plot_diet_pred
gs_gdh_alt_plot_diet3
dev.off()

## __ Carbohydrate metabolism associated functions ----
carbohydrate_gene <- prediction_rel[rownames(prediction_rel)[which(str_detect(prediction_rel$description , "carbohydra")
)],]$description

carbohydrate_db <- prediction_rel %>% filter(description %in% carbohydrate_gene)
# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(carbohydrate_db[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(carbohydrate_db$X.OTU_IDs, each = dim(carbohydrate_db[,2:379])[2]),
                 "description" = rep(carbohydrate_db$description, each = dim(carbohydrate_db[,2:379])[2]),
                 "A" = rep(carbohydrate_db$A, each = dim(carbohydrate_db[,2:379])[2]),
                 "B" = rep(carbohydrate_db$B, each = dim(carbohydrate_db[,2:379])[2]),
                 "C" = rep(carbohydrate_db$C, each = dim(carbohydrate_db[,2:379])[2]),
                 "kegg" = rep(carbohydrate_db$kegg, each = dim(carbohydrate_db[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(carbohydrate_db[2:379]) , dim(carbohydrate_db[,2:379])[1])
# 4 - Create the dataframe
carbohydrate_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(carbohydrate_db[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(carbohydrate_db[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(carbohydrate_db[,2:379])[1]))))
# 6 - Merge the database
carbohydrate_df_2 <- as.data.frame(cbind(carbohydrate_df , diet_df[,c(2,3)]))
carbohydrate_df_2$value <- as.numeric(as.character(carbohydrate_df_2$value))
save(carbohydrate_df_2 , file = "carbohydrate_df_2.RData")

# __ * Plot ----
carbohydrate_plot_diet_pred  <- ggplot(carbohydrate_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ description,  scales = "free")

carbohydrate_plot_diet3  <- ggplot(carbohydrate_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  geom_signif(comparisons = list(c("Carnivorous", "Herbivorous"),c("Carnivorous","Omnivorous"),c("Herbivorous", "Omnivorous")), 
              map_signif_level = TRUE, textsize=4, color="black", family = "serif", vjust = 0.5)+
  facet_wrap( ~ description,  scales = "free")

pdf(file = "carbohydrate_plot.pdf" , he = 7 , wi= 15)
carbohydrate_plot_diet_pred
carbohydrate_plot_diet3
dev.off()

## __ Purine metabolism associated functions ----
purine_db <- prediction_rel %>% filter(C == "Purine metabolism")

# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(purine_db[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(purine_db$X.OTU_IDs, each = dim(purine_db[,2:379])[2]),
                 "description" = rep(purine_db$description, each = dim(purine_db[,2:379])[2]),
                 "A" = rep(purine_db$A, each = dim(purine_db[,2:379])[2]),
                 "B" = rep(purine_db$B, each = dim(purine_db[,2:379])[2]),
                 "C" = rep(purine_db$C, each = dim(purine_db[,2:379])[2]),
                 "kegg" = rep(purine_db$kegg, each = dim(purine_db[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(purine_db[2:379]) , dim(purine_db[,2:379])[1])
# 4 - Create the dataframe
purine_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(purine_db[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(purine_db[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(purine_db[,2:379])[1]))))
# 6 - Merge the database
purine_df_2 <- as.data.frame(cbind(purine_df , diet_df[,c(2,3)]))
purine_df_2$value <- as.numeric(as.character(purine_df_2$value))
save(purine_df_2 , file = "purine_df_2.RData")

save(urate_gene, chitinase_gene, lipase_gene,carnitine_gene,cellulase_gene, amylase_gene,glut_gene, carbohydrate_gene, gs_gdh_alt_gene,file = "selected_genes.RData")
# __ * Plot ----
purine_plot_diet_pred  <- ggplot(purine_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ C,  scales = "free")

purine_plot_diet3  <- ggplot(purine_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  geom_signif(comparisons = list(c("Carnivorous", "Herbivorous"),c("Carnivorous","Omnivorous"),c("Herbivorous", "Omnivorous")), 
              map_signif_level = TRUE, textsize=4, color="black", family = "serif", vjust = 0.5)+
  facet_wrap( ~ C,  scales = "free")

pdf(file = "purine_plot.pdf" , he = 7 , wi= 7)
purine_plot_diet_pred
purine_plot_diet3
dev.off()

## __ All  functions ----

# Transform the db to plot boxplot per sample 
# 1- Transform the matrix into a vector of "Value"
value <-  as.vector(t(as.matrix(prediction_rel[,2:379])))
# 2 - Create the database with the Ko description
kegg_db <- cbind("X.OTU_IDs" = rep(prediction_rel$X.OTU_IDs, each = dim(prediction_rel[,2:379])[2]),
                 "description" = rep(prediction_rel$description, each = dim(prediction_rel[,2:379])[2]),
                 "A" = rep(prediction_rel$A, each = dim(prediction_rel[,2:379])[2]),
                 "B" = rep(prediction_rel$B, each = dim(prediction_rel[,2:379])[2]),
                 "C" = rep(prediction_rel$C, each = dim(prediction_rel[,2:379])[2]),
                 "kegg" = rep(prediction_rel$kegg, each = dim(prediction_rel[,2:379])[2]))
# 3 - Create the "Sample" vector
Sample <- rep(colnames(prediction_rel[2:379]) , dim(prediction_rel[,2:379])[1])
# 4 - Create the dataframe
prediction_df <- as.data.frame(cbind("Sample" = Sample, "value" = value, kegg_db))
# 5 - Join the information concerning the diet
diet_info <- db_core[,c(1,13,14)]
diet_df <- as.data.frame(cbind("Sample" = as.character(rep(diet_info$Sample , dim(prediction_rel[,2:379])[1])),
                               "diet3" = as.character(rep(diet_info$diet3, dim(prediction_rel[,2:379])[1])),
                               "diet_pred" = as.character(rep(diet_info$diet_pred, dim(prediction_rel[,2:379])[1]))))
# 6 - Merge the database
prediction_df_2 <- as.data.frame(cbind(prediction_df , diet_df[,c(2,3)]))
prediction_df_2$value <- as.numeric(as.character(prediction_df_2$value))
save(prediction_df_2 , file = "prediction_df_2.RData")

# __ * Plot ----
plot_diet_pred  <- ggplot(prediction_df_2, aes(x = diet_pred, y = value, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  facet_wrap( ~ A,  scales = "free")

amylase_plot_diet3  <- ggplot(prediction_df_2, aes(x = diet3, y = value, color = diet3)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 0, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, family = "serif"),
        legend.position="none")+
  geom_signif(comparisons = list(c("Carnivorous", "Herbivorous"),c("Carnivorous","Omnivorous"),c("Herbivorous", "Omnivorous")), 
              map_signif_level = TRUE, textsize=4, color="black", family = "serif", vjust = 0.5)+
  facet_wrap( ~ description,  scales = "free")

pdf(file = "amylase_plot.pdf" , he = 7 , wi= 15)
amylase_plot_diet_pred
amylase_plot_diet3
dev.off()

### ______ On level B----
prediction_df_3 <- prediction_df_2
prediction_df_3$B <- as.factor(str_to_lower(prediction_df_3$B))

prediction_df_3$B <- str_replace_all(prediction_df_3$B,c("replication, recombination and repair proteins" =  "replication and repair", 
                                                         "translation proteins" = "translation",
                                                         "transcription related proteins" = "transcription",
                                                         "biosynthesis of other secondary metabolites" = "biosynthesis and biodegradation of secondary metabolites",
                                                         "signal transduction"= "Cellular processes and Signaling mechanisms",
                                                         "signal transduction mechanisms" = "Cellular processes and Signaling mechanisms",
                                                         "signaling molecules and interaction" = "Cellular processes and Signaling mechanisms", 
                                                         "Cellular processes and Signaling and secretion" = "Cellular processes and Signaling mechanisms",
                                                         "cell growth and death" = "Cellular processes and Signaling mechanisms",
                                                         "Cellular processes and Signaling mechanisms mechanisms" = "Cellular processes and Signaling mechanisms",
                                                         "Cellular processes and Signaling mechanisms and secretion" = "Cellular processes and Signaling mechanisms",
                                                         "cell division" = "Cellular processes and Signaling mechanisms",
                                                         "cell motility" = "Cellular processes and Signaling mechanisms",
                                                         "membrane and intracellular structural molecules" = "Cellular processes and Signaling mechanisms",
                                                         "membrane transport" = "transport and catabolism",
                                                         "other ion-coupled transporters" = "transport and catabolism",
                                                         "inorganic ion transport and metabolism" =  "transport and catabolism",
                                                         "other transporters" = "transport and catabolism",
                                                         "pores ion channels"="transport and catabolism",
                                                         "electron transfer carriers" =  "transport and catabolism",
                                                         "protein folding and associated processing" = "folding, sorting and degradation",
                                                         "metabolism of other amino acids" = "amino acid metabolism",
                                                         "germination" = "sporulation"))


B_ab <- prediction_df_3 %>% 
  group_by(B) %>%
  summarize(Sum = sum(value, na.rm=TRUE))

B_ab <- B_ab[order(-B_ab$Sum),]
B_ab$B
B_selected <- B_ab %>% filter(!B %in% c(NA,"general function prediction only", "function unknown", "enzyme families", "others"))

B_12 <- as.character(B_selected[1:12,]$B)
B_10 <- as.character(B_selected[1:10,]$B)

B_other <- B_selected %>% filter(!B %in% B_10) %>% select(B)
B_other <- as.character(B_other$B)

prediction_df_3$B = as.character(prediction_df_3$B) # Avoid error message with factor for next step
save(prediction_df_3, file = "prediction_df_3.RData")
load("prediction_df_3.RData")
prediction_df_3[!prediction_df_3$B %in% B_10,
                  which(names(prediction_df_3) == "B", T)] <- "z-Other" #Change phylum to other for those not included in the list

# Transform the data frame in sub dataframe in order to calculate the % relatives in each dietary classes.
cora_pred <- prediction_df_3 %>% filter(diet_pred == "Corallivores")
cora_pred$value <- cora_pred$value/sum(cora_pred$value)
macroinv_pred <- prediction_df_3 %>% filter(diet_pred == "Macroinvertivores")
macroinv_pred$value <- macroinv_pred$value/sum(macroinv_pred$value)
plank_pred <- prediction_df_3 %>% filter(diet_pred == "Planktivores")
plank_pred$value <- plank_pred$value/sum(plank_pred$value)
crust_pred <- prediction_df_3 %>% filter(diet_pred == "Crustacivores")
crust_pred$value <- crust_pred$value/sum(crust_pred$value)
microinv_pred <- prediction_df_3 %>% filter(diet_pred == "Microinvertivores")
microinv_pred$value <- microinv_pred$value/sum(microinv_pred$value)
si_pred <- prediction_df_3 %>% filter(diet_pred == "sessile invertivores")
si_pred$value <- si_pred$value/sum(si_pred$value)
pisci_pred <- prediction_df_3 %>% filter(diet_pred == "Piscivores")
pisci_pred$value <- pisci_pred$value/sum(pisci_pred$value)
hmd_pred <- prediction_df_3 %>% filter(diet_pred == "Herbivores Microvores Detritivores")
hmd_pred$value <- hmd_pred$value/sum(hmd_pred$value)

prediction_df_rel_diet_pred <- rbind(cora_pred,macroinv_pred,plank_pred,crust_pred,microinv_pred,si_pred,
                                     pisci_pred,hmd_pred)

barplot <- prediction_df_rel_diet_pred %>% filter(!B %in% c(NA , "z-Other"))
barplot$B <- factor(barplot$B, levels = B_selected$B)

test <- barplot %>% group_by(B,diet_pred) %>% summarise(Sum = sum(value))

pdf(file = "test.pdf" , he = 10 ,wi = 7)
ggplot(test) +
  geom_bar(aes(x= Sum , y = diet_pred, fill = diet_pred), stat="identity",alpha=0.7) +
  theme_bw() +
  scale_fill_manual(values= col_diet_pred) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(size=0, family = "serif"),
        axis.text = element_text(family = "serif",size = 0),
        axis.text.x = element_text(family = "serif",size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=0, family = "serif"),
        strip.background = element_rect(colour = "black", fill = "white"),
        legend.position="none")+
  facet_grid(B ~. , scales = "free")
dev.off()



