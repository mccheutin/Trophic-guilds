#### Analysis of functional diversity -------

setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Carnivory/Functions")
load("prediction_files.RData") ; load("db_core.RData")

## Prepare the data ---- 
library(dplyr)
library(phyloseq)

otu_table_fun <- prediction_rel %>%
  select(!c(X.OTU_IDs, description, A, B, C, kegg))
rownames(otu_table_fun) <-  prediction_rel$X.OTU_IDs
otu_table_fun <- otu_table(t(otu_table_fun), taxa_are_rows = FALSE)

tax_table_fun <- prediction_rel %>%
  select(c(X.OTU_IDs, description, A, B, C, kegg))
tax_table_fun <- tax_table(tax_table_fun)
rownames(tax_table_fun) <- colnames(otu_table_fun)

ps_fun <- phyloseq(otu_table_fun,tax_table_fun)
DAT <- sample_data(db_core)
sort(sample_names(DAT)) == sort(sample_names(ps_fun))
ps_fun <- merge_phyloseq(ps_fun, DAT)
save(ps_fun, file = "ps_fun.RData")


########  *** PCOA ***** -----
library(vegan)
library(ape)
fun.bc <- vegdist(ps_fun@otu_table, method = "bray")
pcoa.sub.bc.fun <- pcoa(fun.bc)
pcoa_coord.bc.fun <- pcoa.sub.bc.fun$vectors[,1:3]

fun.eu <- vegdist(ps_fun@otu_table, method= "euclidian")
pcoa.sub.eu.fun <- pcoa(fun.eu)
pcoa_coord.eu.fun <- pcoa.sub.eu.fun$vectors[,1:3]

# Contruction of the table for graphic 
library(stringr)
hull.bc.fun <- cbind(pcoa_coord.bc.fun, db_core)
hull.eu.fun <- cbind(pcoa_coord.eu.fun, db_core)
save(hull.bc.fun,hull.eu.fun, file = "hulls.fun.RData")

# What is the percentage of the explicative variance? 
library(scales)
paste("Axis 1 :",percent(pcoa.sub.bc.fun$values$Relative_eig[1])) # 36 %
paste("Axis 2 :",percent(pcoa.sub.bc.fun$values$Relative_eig[2])) # 11 %
paste("Axis 3 :",percent(pcoa.sub.bc.fun$values$Relative_eig[3])) # 10 %
paste("Axis 4 :",percent(pcoa.sub.bc.fun$values$Relative_eig[4])) # 6 %
save(fun.bc, hull.bc.fun, file = "beta.bc.fun.RData")

paste("Axis 1 :",percent(pcoa.sub.eu.fun$values$Relative_eig[1])) # 28 %
paste("Axis 2 :",percent(pcoa.sub.eu.fun$values$Relative_eig[2])) # 14 %
paste("Axis 3 :",percent(pcoa.sub.eu.fun$values$Relative_eig[3])) # 9 %
paste("Axis 4 :",percent(pcoa.sub.eu.fun$values$Relative_eig[4])) # 9 %
save(fun.eu, hull.eu.fun, file = "beta.eu.fun.RData")

## __ Plot ----
col_diet_pred <- c('Corallivores' =  'indianred1',
                   'Crustacivores' = 'orange',
                   'Herbivores Microvores Detritivores' = 'darkgreen',
                   'Macroinvertivores' = 'darkslategray',
                   'Microinvertivores' = 'bisque3',
                   'Piscivores' = 'darkblue',
                   'Planktivores' = 'darkred',
                   'sessile invertivores' = 'coral4',
                   'Omnivorous' = 'gray',
                   'Carnivorous' = 'darkslateblue')


pcoa.diet.pred.eu <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = col_diet_pred) +
  geom_point(data = hull.eu.fun, aes(x=Axis.1, y=Axis.2, color = diet_pred), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.eu.fun$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.eu.fun$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  #stat_ellipse(data = hull.bc,aes(x=Axis.1, y=Axis.2,color=diet_pred),type = "norm")+
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=14, family = "serif"),
        axis.text.y = element_text(size=14, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.position = "none",
        legend.text = element_text(size = 14, family = "serif"),
        legend.title = element_text(size = 16,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))+ labs(color = "Diet prediction") #+facet_wrap(~ region)
        
pcoa.diet.pred.eu

pcoa.eu.diet3 <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  geom_point(data = hull.eu.fun, aes(x=Axis.1, y=Axis.2, color = diet3), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.eu.fun$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.eu.fun$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  #stat_ellipse(data = hull.bc,aes(x=Axis.1, y=Axis.2,color=diet_pred),type = "norm")+
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=14, family = "serif"),
        axis.text.y = element_text(size=14, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.position = "none",
        legend.text = element_text(size = 14, family = "serif"),
        legend.title = element_text(size = 16,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))+ labs(color = "Trophic guilds") #+facet_wrap(~ region)

pcoa.eu.diet3

# Save the plots 
library(gridExtra)
A = pcoa.diet.pred.bc + xlim(-0.25,0.4) +ylim(-0.2,0.35)
B = pcoa.diet3 + xlim(-0.25,0.4) +ylim(-0.2,0.35)

pdf(file = "function_bc_pcoa.pdf", he = 7 , wi = 7)
plot_bc_fun <- grid.arrange(A, B, ncol=2, nrow = 1)
dev.off()

pdf(file = "function_eu_pcoa.pdf", he = 7 , wi = 7)
plot_eu_fun <- grid.arrange(pcoa.diet.pred.eu, pcoa.eu.diet3, ncol=2, nrow = 1)
dev.off()

### __ PERMANOVA ------
library(vegan)
a.bc = adonis(fun.bc ~ diet_pred, data = db_core)
b.bc = adonis(fun.bc ~ diet3, data = db_core)

c.eu = adonis(fun.eu ~ diet_pred, data = db_core)
d.eu = adonis(fun.eu ~ diet3, data = db_core)

########  *** Heatmap ***** -----
ko_hm_kegg <- aggregate(prediction_rel[,2:379], by=list(description= prediction_rel$description), FUN=sum)
rownames(ko_hm_kegg) <- ko_hm_kegg$description
ko_hm_kegg <- ko_hm_kegg %>% select(!description)
# We transform the matrices in %
ko_hm_100 <- ko_hm_kegg * 100 # ko_hm_kegg is already relative

# Change the names of the samples
fish_tax <- cbind(ID = rownames(sample_data(ps_fun)) , sample_data(ps_fun)[,4],sample_data(ps_fun)[,13])
fish_id <- as.data.frame(cbind(colnames(ko_hm_100) ,rep("bla", length(colnames(ko_hm_100))),rep("blo", length(colnames(ko_hm_100)))))
colnames(fish_id) <- c("ID" , "tax","diet3")
fish_id <- left_join(fish_id, fish_tax, by = "ID")
colnames(ko_hm_100) <- fish_tax$diet


my_group <- c(rep(1,length(rownames(db_core %>% filter(diet3 == "Carnivorous")))),
              rep(2, length(rownames(db_core %>% filter(diet3 == "Herbivorous")))),
              rep(3, length(rownames(db_core %>% filter(diet3 == "Omnivorous")))))

colSide <- c("darkred", "darkgreen","darkblue")[my_group]


#selected_functions <- unique(c(urate_gene, chitinase_gene, lipase_gene,carnitine_gene,cellulase_gene, amylase_gene,glut_gene, carbohydrate_gene, gs_gdh_alt_gene))
selected_functions <- unique(c("carbohydrate diacid regulator","carbohydrate-specific outer membrane porin","RpiR family transcriptional regulator, carbohydrate utilization regulator",
                               "chitinase [EC:3.2.1.14]","chitin-binding protein","chitin disaccharide deacetylase [EC:3.5.1.105]",
                               "allantoinase [EC:3.5.2.5]","allantoin racemase [EC:5.1.99.3]","FAD-dependent urate hydroxylase [EC:1.14.13.113]", "ureidoglycolate lyase [EC:4.3.2.3]","ureidomalonase [EC:3.5.1.95]",
                               "pectate lyase [EC:4.2.2.2]","alpha-N-arabinofuranosidase [EC:3.2.1.55]" ,"cellobiose phosphorylase [EC:2.4.1.20]","cellulose synthase operon protein B",
                               "MerR family transcriptional regulator, glutamine synthetase repressor","N-[(2S)-2-amino-2-carboxyethyl]-L-glutamate dehydrogenase [EC:1.5.1.51]","methylglutamate dehydrogenase subunit C [EC:1.5.99.5]",
                               "alanine-synthesizing transaminase [EC:2.6.1.-]","methylglutamate dehydrogenase subunit B [EC:1.5.99.5]","methylglutamate dehydrogenase subunit D [EC:1.5.99.5]"
                               ))

ko_hm_100_sel <- ko_hm_100[which(rownames(ko_hm_100) %in% selected_functions),]


j

pdf(file = "./plots/heatmap_fun.pdf", he = 10 , wi = 10)
heatmap.2(as.matrix(ko_hm_100_sel), Rowv = F,Colv=F,scale = "row", col = greenred(100),
          trace = "none", density.info = "none", 
          ColSideColors=colSide, dendrogram = "none",
          key.title = "Functional Contribution")

dev.off()


########
########
p <- ggplot(ko_hm_100_sel, aes(x = ASV, y = Host_family)) + 
  geom_tile(aes(fill = as.factor(value)),colour = "white") + 
  scale_fill_manual(values=c("gray","red","black"))+
  theme(axis.line = element_line(size = 0), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 3, family = "serif", angle = 90,vjust = 0.65), 
        axis.text.y = element_text(family = "serif", size = 12,face = "italic"))

library(RColorBrewer)
heatmap(as.matrix(ko_hm_100_sel),ColSideColors=colSide, scale="column",Colv = NA, Rowv = NA, col= colorRampPalette(brewer.pal(8, "Blues"))(25))



