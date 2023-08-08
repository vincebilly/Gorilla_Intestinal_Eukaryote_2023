################################################################################################################
################################################################################################################
################################################################################################################
###                                                                                                          ###
###                   Analysis of 18S V4 Gorilla - interaction_btw_euk_gorilla_features.R                    ###
###                                                                                                          ###
################################################################################################################
################################################################################################################
################################################################################################################

#Clear environment
rm(list=ls())

########################################################
###                 Import Libraries                 ###
########################################################

library("ggridges")    #
library("treeio")      #
library("tidyverse")   #
library("reshape2")    #Load the reshape2 packAge for converting between long and wide format data
library("stringr")     #Load the stringr packAge for improved filtering of text
library("ape")         #Load the ape packAge for reading and modifying phylogenetic trees
library("phyloseq")    #Load the phyloseq packAge for microbial community analysis
library("data.table")  #Load the data.table packAge for better metadata manipulation
library("viridis")     #Load the viridis packAge for colour palettes for continuous data
library("qualpalr")    #Load the qualpalr packAge for colour palettes for qualitative data
library("ggplot2")     #load the ggplot2 packAge for visualization of data
library("vegan")       # 
#library("gplots")      #
library("ggpubr")      #
library("ggforce")     #
library("plyr")        #
library("lattice")     #
library("Rmisc")       #
library("funrar")      #
library("permute")     #
library("vegan")       #
library("readr")       #
library("forcats")     #
library("RColorBrewer")#
library("ade4")        #
library("cluster")     #
library("clue")        #
library("dada2")       #
library("permute")     #
library("broom")       #
library("viridisLite") #
library("seqinr")      #
library("ShortRead")   #
library("Biostrings")  #
library("VennDiagram") #
library("DESeq2")      #
library("ggtree")      #
library("aplot")       #
library("tidyr")       #
library("ggimAge")     #
library("lubridate")   #Add vertical line in time serie plot
library("onewaytests") #
library("gridExtra")   #
library("ggeffects")   #
library("cowplot")
library("dplyr")
library("ggstance")
library("tidytree")
library("ggalt")
library("tibble")


########################################################

########################################################
###              Importation of data                 ###
########################################################

##Importation of OTU table, which also include taxonomy
asv_rnk_spl_neg_df <- fread("/Users/vincebilly/Desktop/vince/phd/assignment/gorilla/sequence_table.18s.R1.uniqueIDs_modified3.csv", sep="\t", header = TRUE, data.table=FALSE);dim(asv_rnk_spl_neg_df)
rownames(asv_rnk_spl_neg_df) <- asv_rnk_spl_neg_df$SeqID                             #Rename row.name with ASV numbers
asv_rnk_spl_neg_df2 <- asv_rnk_spl_neg_df[,-1]                                       #New table without column with ASV IDs
head(asv_rnk_spl_neg_df2);dim(asv_rnk_spl_neg_df2);class(asv_rnk_spl_neg_df2)        #Visualisation for sanity check

##Importation of metadata
metadata_spl_ctl_df <- fread("/Users/vincebilly/Desktop/vince/phd/assignment/gorilla/mapping_file_Gor_Ren_18s_2018-5.csv", sep="\t", header = TRUE, data.table=FALSE)
rownames(metadata_spl_ctl_df) <- metadata_spl_ctl_df$SampleID                         #Rename row.name with sample IDS
head(metadata_spl_ctl_df); str(metadata_spl_ctl_df) ; dim(metadata_spl_ctl_df)        #Vizualisation for sanity check
#NOTES: 96 rows (samples and control) and 16 columns (metadata)

##Import trees
tree_blastocystis <- read.tree("/Users/vincebilly/Desktop/vince/phd/phylogeny/18s/protist/phylo_blastocystis/constrained_tree/tree/gorilla/70/reroot.tre")
tree_diplomonad <- read.tree("/Users/vincebilly/Desktop/vince/phd/phylogeny/18s/protist/phylo_dipent/tree/epa/gorilla/80/reroot.tre")
tree_trichostomatia <- read.tree("/Users/vincebilly/Desktop/vince/phd/phylogeny/18s/protist/phylo_trichostomatia/tree/epa/epa_gorilla/90/tree_trichostomatia.tre")
tree_trichomonad <- read.tree("/Users/vincebilly/Desktop/vince/phd/phylogeny/18s/protist/phylo_trichomonad/tree/epa/epa_gorilla/96/tree_trichomonad.tre")
tree_entamoeba <- read.tree("/Users/vincebilly/Desktop/vince/phd/phylogeny/18s/protist/phylo_entamoeba/tree/epa_gorilla/83/tree_entamoeba.tre")
tree_nematoda <- read.tree("/Users/vincebilly/Desktop/vince/phd/phylogeny/18s/worm/phylo_nematoda/tree/epa/epa_gorilla/90/tree_nematoda.tre")
#Combine trees
new_tree <- bind.tree(tree_blastocystis, tree_diplomonad, where = "root", position = 0, interactive = FALSE)
new_tree <- bind.tree(new_tree,tree_trichostomatia)
new_tree <- bind.tree(new_tree,tree_trichomonad)
new_tree <- bind.tree(new_tree,tree_entamoeba)
new_tree <- bind.tree(new_tree,tree_nematoda)
#test plot
plot(new_tree)

########################################################

################################################################################################################
################################################################################################################
################################################################################################################
###                                                                                                          ###
###                                 ANALYSIS OF GORILLAS METADATA                                            ###
###                                                                                                          ###
################################################################################################################
################################################################################################################
################################################################################################################

########################################################
### Co-association between factors for fecal samples ###
########################################################

metadata_spl_df <- subset(metadata_spl_ctl_df, Gorilla_name %in% setdiff(metadata_spl_ctl_df$Gorilla_name,c("Neg1","Neg2","Neg3","Neg4")))
metadata_spl_df <- metadata_spl_df[,c("Site","Season","Yaws","Social_Status","Age","Cluster_Duplicate")]
metadata_spl_df <- unique(metadata_spl_df)


#Site and Season
table <- table(metadata_spl_df$Site, metadata_spl_df$Season);table
fisher.test(table)

#Site and Season without Maya Moba
metadata_spl_df_rom_lok <- subset(metadata_spl_df, Site %in% c("Romani","Lokoue"))
table <- table(metadata_spl_df_rom_lok$Site, metadata_spl_df_rom_lok$Season);table
fisher.test(table)

#Site and Age
table <- table(metadata_spl_df_rom_lok$Site, metadata_spl_df_rom_lok$Age);table
fisher.test(table)

#Site and Age without subadult
metadata_spl_df_rom_lok_bb_sb <- subset(metadata_spl_df, Age %in% c("BB","SB"))
table <- table(metadata_spl_df_rom_lok_bb_sb$Site, metadata_spl_df_rom_lok_bb_sb$Age);table
fisher.test(table)

#Site and Social status
table <- table(metadata_spl_df_rom_lok$Site, metadata_spl_df_rom_lok$Social_Status);table
fisher.test(table)

#Site and Yaws
table <- table(metadata_spl_df_rom_lok$Site, metadata_spl_df_rom_lok$Yaws);table
fisher.test(table)

#Yaws and Age
table <- table(metadata_spl_df_rom_lok$Yaws, metadata_spl_df_rom_lok$Age);table
fisher.test(table)

#Yaws and Age without subadult
table <- table(metadata_spl_df_rom_lok_bb_sb$Yaws, metadata_spl_df_rom_lok_bb_sb$Age);table
fisher.test(table)

#Yaws and social status
table <- table(metadata_spl_df_rom_lok$Yaws, metadata_spl_df_rom_lok$Social_Status);table
fisher.test(table)

#Yaws and Season
table <- table(metadata_spl_df_rom_lok$Yaws, metadata_spl_df_rom_lok$Season);table
fisher.test(table)

#Social status and Age
table <- table(metadata_spl_df_rom_lok$Social_Status, metadata_spl_df_rom_lok$Age);table
fisher.test(table)

#Social status and Age without subadult
table <- table(metadata_spl_df_rom_lok_bb_sb$Social_Status, metadata_spl_df_rom_lok_bb_sb$Age);table
fisher.test(table)

#Social status and Season
table <- table(metadata_spl_df_rom_lok$Social_Status, metadata_spl_df_rom_lok$Season);table
fisher.test(table)

#Season and Age 
table <- table(metadata_spl_df_rom_lok$Season, metadata_spl_df_rom_lok$Age);table
fisher.test(table)

#Season and Age without subadult
table <- table(metadata_spl_df_rom_lok_bb_sb$Season, metadata_spl_df_rom_lok_bb_sb$Age);table
fisher.test(table)

########################################################

########################################################
### Co-association between factors for individuals   ###
########################################################

metadata_spl_df <- subset(metadata_spl_ctl_df, gorilla_name %in% setdiff(metadata_spl_ctl_df$gorilla_name,c("Neg1","Neg2","Neg3","Neg4")))
metadata_spl_df <- metadata_spl_df[,c("gorilla_name","Site","Yaws","Social_Status","Age")]

metadata_spl_df <- unique(metadata_spl_df)
metadata_spl_df_rom_lok <- subset(metadata_spl_df, Site %in% c("Romani","Lokoue"))


#Site and Season without Maya Moba
#table <- table(metadata_spl_df_rom_lok$Site, metadata_spl_df_rom_lok$Season);table
#fisher.test(table)

#Site and Yaws
table <- table(metadata_spl_df_rom_lok$Site, metadata_spl_df_rom_lok$Yaws);table
fisher.test(table)

#Site and Social status
table <- table(metadata_spl_df_rom_lok$Site, metadata_spl_df_rom_lok$Social_Status);table
fisher.test(table)

#Site and Age
table <- table(metadata_spl_df_rom_lok$Site, metadata_spl_df_rom_lok$Age);table
fisher.test(table)

#Site and Age without subadult
metadata_spl_df_rom_lok_bb_sb <- subset(metadata_spl_df, Age %in% c("BB","SB"))
table <- table(metadata_spl_df_rom_lok_bb_sb$Site, metadata_spl_df_rom_lok_bb_sb$Age);table
fisher.test(table)

#Yaws and Age
table <- table(metadata_spl_df_rom_lok$Yaws, metadata_spl_df_rom_lok$Age);table
fisher.test(table)

#Yaws and Age without subadult
table <- table(metadata_spl_df_rom_lok_bb_sb$Yaws, metadata_spl_df_rom_lok_bb_sb$Age);table
fisher.test(table)

#Yaws and social status
table <- table(metadata_spl_df_rom_lok$Yaws, metadata_spl_df_rom_lok$Social_Status);table
fisher.test(table)

#Social status and Age
table <- table(metadata_spl_df_rom_lok$Social_Status, metadata_spl_df_rom_lok$Age);table
fisher.test(table)

#Social status and Age without subadult
table <- table(metadata_spl_df_rom_lok_bb_sb$Social_Status, metadata_spl_df_rom_lok_bb_sb$Age);table
fisher.test(table)

########################################################

################################################################################################################
###                                  Sequencing depth across factors                 ###                     ###
################################################################################################################

########################################################
## Extract taxonomy, ASV and control tables from data ##
########################################################

### Extract taxonomy table ###
names_rnk <- c("Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Rank8","Accession_number")
asv_tax.df <- asv_rnk_spl_neg_df2[,names_rnk]
#Because of presence of NA, need to rename unknwown taxa at each rank (once printing makes easier selection of taxa to create table)
asvUn_tax.df <- asv_tax.df                          #Create new table  
asvUn_tax.df[is.na(asvUn_tax.df)] <- "Un"
#create a function to paste "Un" to the name of the previous rank
mypaste = function(s) {
  s = gsub("_Un","",s)
  return(paste0(s,"_Un"))
}
#apply the function to all ranks
for (i in 2:9){
  
  if (i == 9)
  {
    colName = "Accession_number"
    colNamePrevious = "Rank8"
  } else
  {
    colName = paste0("Rank", i)
    colNamePrevious = paste0("Rank", i-1)
  }
  w = which(asvUn_tax.df[[colName]] == "Un")
  asvUn_tax.df[[colName]][w] = mypaste(asvUn_tax.df[[colNamePrevious]][w])
}
asvUn_tax.mtx <- as.matrix(asvUn_tax.df);dim(asvUn_tax.mtx)


### Extract ASV-Samples table ###
names_spl_only <- colnames(asv_rnk_spl_neg_df2)
asv_spl.df <- asv_rnk_spl_neg_df2[,names_spl_only];str(asv_spl.df)
asv_spl.df$G0009 <- as.numeric(asv_spl.df$G0009)
asv_spl.df <- asv_spl.df %>% mutate_all(as.numeric)

########################################################

########################################################
x <- asv_rnk_spl_neg_df2[,grep("G",colnames(asv_rnk_spl_neg_df2))]
x <- x %>% mutate_all(as.numeric)
x[is.na(x)] <- 0

Seq_depth <- as.data.frame(colSums(x))

x2 <- merge(Seq_depth,metadata_spl_ctl_df,by="row.names",all=FALSE)
colnames(x2)[2] <- "Seq_depth"

ggplot(x2, aes(x=Site, y=Seq_depth)) +
  geom_boxplot() +
  geom_jitter() + 
  stat_compare_means(paired = TRUE)

x3 <- subset(x2, Site !="Moba")

ggplot(x3, aes(x=Site, y=Seq_depth)) +
  geom_boxplot() +
  geom_jitter() + 
  stat_compare_means()
########################################################

################################################################################################################
###                         Comparison between gut ASV and non-gut ASV (plants, fungi , ...)                 ###                     ###
################################################################################################################

########################################################
## Extract taxonomy, ASV and control tables from data ##
########################################################

### Extract taxonomy table ###
names_rnk <- c("Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Rank8","Accession_number")
asv_tax.df <- asv_rnk_spl_neg_df2[,names_rnk]
#Because of presence of NA, need to rename unknwown taxa at each rank (once printing makes easier selection of taxa to create table)
asvUn_tax.df <- asv_tax.df                          #Create new table  
asvUn_tax.df[is.na(asvUn_tax.df)] <- "Un"
#create a function to paste "Un" to the name of the previous rank
mypaste = function(s) {
  s = gsub("_Un","",s)
  return(paste0(s,"_Un"))
}
#apply the function to all ranks
for (i in 2:9){
  
  if (i == 9)
  {
    colName = "Accession_number"
    colNamePrevious = "Rank8"
  } else
  {
    colName = paste0("Rank", i)
    colNamePrevious = paste0("Rank", i-1)
  }
  w = which(asvUn_tax.df[[colName]] == "Un")
  asvUn_tax.df[[colName]][w] = mypaste(asvUn_tax.df[[colNamePrevious]][w])
}
asvUn_tax.mtx <- as.matrix(asvUn_tax.df);dim(asvUn_tax.mtx)


### Extract ASV-Samples table ###
names_spl_only <- colnames(asv_rnk_spl_neg_df2)
asv_spl.df <- asv_rnk_spl_neg_df2[,names_spl_only];str(asv_spl.df)
asv_spl.df$G0009 <- as.numeric(asv_spl.df$G0009)
asv_spl.df <- asv_spl.df %>% mutate_all(as.numeric)

########################################################

########################################################
###       Analysis with Gut and Non-Gut              ###
########################################################

asv_rnk_spl_neg_df2$Origin <- asv_rnk_spl_neg_df2$Symbiont
asv_rnk_spl_neg_df2$Origin <- ifelse(asv_rnk_spl_neg_df2$Origin %in% c("gut_worm","gut_protist"),"gut","non-gut")

to_keep <- colnames(asv_rnk_spl_neg_df2)[grep("G0|G1|neg",colnames(asv_rnk_spl_neg_df2))]
x <- asv_rnk_spl_neg_df2[,c("Origin",to_keep)]
x$G0009 <- as.numeric(x$G0009)
x[,-1] <- x[,-1] %>% mutate_all(as.numeric)

x1 <- as.data.frame(x %>% group_by(Origin) %>% 
  summarise_each(funs(sum)))
x1[is.na(x1)] <- 0

#Pie chart proportion number of reads from gut and non gut
total_gut_nongut <- x1
row.names(total_gut_nongut) <- total_gut_nongut$Origin
total_gut_nongut$Origin <- NULL
total_reads_origin <- rowSums(total_gut_nongut)
count <- c(total_reads_origin[1], total_reads_origin[2])
colours = c("grey40","gold")
pie(count, col = colours)
count[1]/(count[1]+count[2])*100


#Barplot with total reads and total gut and non gut reads
Total_reads <- as.data.frame(colSums(x1[,-1]))
x2 <- melt(x1)
x3 <- merge(x2,Total_reads, by.x="variable",by.y="row.names")
colnames(x3)[4] <- "Total_reads"
x4 <- merge(x3,metadata_spl_ctl_df, by.x="variable",by.y="row.names", all.y=FALSE)
x4$Origin <- factor(x4$Origin, levels = c("non-gut","gut"))
colours = c("gold","grey40")
x4$Site <- gsub("neg.*","Neg",x4$Site)
x4$Site <- factor(x4$Site, level=c("Lokoue","Romani","Maya","Moba","Neg"))
x4$gorilla_name_sampleID <- paste0(x4$gorilla_name, "_",x4$variable)
plot1 <- ggplot(x4, aes(x=reorder(x4$gorilla_name_sampleID,x4$gorilla_name_sampleID), y=value, fill=Origin)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid", size=0.01) + #height of the column equal the value, stacked
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=15),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=15),
        legend.position = "left",
        legend.spacing.y = unit(0.2, 'cm'),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        strip.text.x = element_text(size = 20)) +
  guides(fill=guide_legend(ncol=1,byrow = TRUE))
print(plot1)

#distribution of reads per samples
summary(x4$Total_reads)
x5 <- subset(x4, Origin =="gut")
write_csv2(x5, "/Users/vincebilly/Desktop/Sup_table1.csv")

#relative abundance
x1 <- as.data.frame(x %>% group_by(Origin) %>% 
  summarise_each(funs(sum)))

row.names(x1) <- x1$Origin
x1$Origin <- NULL
x1[is.na(x1)] <- 0
x1 <- t(x1)
x1 <- x1/rowSums(x1)
x1 <- t(x1)

x2 <- melt(x1)
x3 <- merge(x2,Total_reads, by.x="Var2",by.y="row.names")
colnames(x3)[4] <- "Total_reads"
x4 <- merge(x3,metadata_spl_ctl_df, by.x="Var2",by.y="row.names", all.y=FALSE)
x4$Var1 <- factor(x4$Var1, levels = c("non-gut","gut"))
colours = c("gold","grey40")
x4$Site <- gsub("neg.*","Neg",x4$Site)
x4$Site <- factor(x4$Site, level=c("Lokoue","Romani","Maya","Moba","Neg"))
x4$gorilla_name_sampleID <- paste0(x4$gorilla_name, "_",x4$Var2)
plot2 <- ggplot(x4, aes(x=reorder(Var2,Total_reads), y=value, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid", size=0.01) + #height of the column equal the value, stacked
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=15),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        legend.position = "left",
        legend.text=element_blank(),
        legend.title=element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank()) +
  guides(fill=guide_legend(ncol=1,byrow = TRUE))
print(plot2)
dev.off()

#Combined plot together
combined_plot <- ggarrange(plot1, plot2, ncol = 2, widths = c(1.2,0.8))
ggsave(filename = "Sup_Figure_3_fecal_bytotalreads.pdf", path = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/", plot= combined_plot, width = 17, height = 20, device = "pdf" )
########################################################

########################################################
###       Analysis with type of Eukaryotes           ###
########################################################

to_keep <- colnames(asv_rnk_spl_neg_df2)[grep("G0|G1|neg",colnames(asv_rnk_spl_neg_df2))]
x <- asv_rnk_spl_neg_df2[,c("Symbiont",to_keep)]
x$G0009 <- as.numeric(x$G0009)
x[,-1] <- x[,-1] %>% mutate_all(as.numeric)

#Pie chart for number of ASV per Origin
Rabk8_ASV.tbl <- x
Rabk8_ASV.tbl$ASV <- row.names(Rabk8_ASV.tbl)
Rabk8_ASV.tbl <- Rabk8_ASV.tbl[,c("Symbiont","ASV")]
Rabk8_ASV_count <- as.data.frame(table(Rabk8_ASV.tbl$Symbiont))
colnames(Rabk8_ASV_count) <- c("Symbiont", "Amount")
Rabk8_ASV_count$Relative_Abundance <- Rabk8_ASV_count$Amount/sum(Rabk8_ASV_count$Amount)*100
Rabk8_ASV_count$Symbiont <- factor(Rabk8_ASV_count$Symbiont, level = c("Arthropoda","Chimeras","Environmental micro-eukaryotes","Fungi","Gut_protist",
                                                                       "Gut_worm","Insect_parasite","Mammalia","Non-gut_parasite",
                                                                       "Occasional_parasite","Plant","Unknown"))


pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/filtering_taxa/pie_chart_ASV.pdf", height =  10, width = 15)
plot_ASV <- ggplot(Rabk8_ASV_count, aes(x="", y=Amount, fill=Symbiont)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = Amount), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom") +
  guides(fill=guide_legend(nrow=2,byrow = TRUE)) +
  scale_fill_manual(values=c("purple","orange","dodgerblue3","chocolate3","grey40","grey70","plum","gold","red","tomato","forestgreen","wheat"))
print(plot_ASV)
dev.off()

#Pie chart proportion number of reads from gut and non gut
x1 <- as.data.frame(x %>% group_by(Symbiont) %>% 
                      summarise_each(funs(sum)))
x1[is.na(x1)] <- 0
Rank8_reads.tbl <- x1
row.names(Rank8_reads.tbl) <- Rank8_reads.tbl$Symbiont
Rank8_reads.tbl$Symbiont <- NULL
Rabk8_reads_count <- as.data.frame(rowSums(Rank8_reads.tbl))
Rabk8_reads_count$Symbiont <- row.names(Rabk8_reads_count)
colnames(Rabk8_reads_count) <- c("Amount","Symbiont")
Rabk8_reads_count$Relative_Abundance <- Rabk8_reads_count$Amount/sum(Rabk8_reads_count$Amount)*100
Rabk8_reads_count$Symbiont <- factor(Rabk8_reads_count$Symbiont, level = c("Arthropoda","Chimeras","Environmental micro-eukaryotes","Fungi","Gut_protist",
                                                                           "Gut_worm","Insect_parasite","Mammalia","Non-gut_parasite",
                                                                           "Occasional_parasite","Plant","Unknown"))

plot_reads <- ggplot(Rabk8_reads_count, aes(x="", y=Amount, fill=Symbiont)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = Amount), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")  +
  guides(fill=guide_legend(nrow=2,byrow = TRUE)) +
  scale_fill_manual(values=c("purple","orange","dodgerblue3","chocolate3","grey40","grey70","plum","gold","red","tomato","forestgreen","wheat"))
print(plot_reads)

#Combined plot together
combined_plot <- ggarrange(plot_ASV, plot_reads, ncol = 2, widths = c(1,1))
ggsave(filename = "Pie_chart_ASV_reads.pdf", path = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/filtering_taxa/", plot= combined_plot, width = 15, height = 7, device = "pdf" )

#Barplot
x1 <- as.data.frame(x %>% group_by(Symbiont) %>% 
                      summarise_each(funs(sum)))
x1[is.na(x1)] <- 0

Total_reads <- as.data.frame(colSums(x1[,-1]))
x2 <- melt(x1)
x3 <- merge(x2,Total_reads, by.x="variable",by.y="row.names")
colnames(x3)[4] <- "Total_reads"
x4 <- merge(x3,metadata_spl_ctl_df, by.x="variable",by.y="row.names", all.y=FALSE)
x4$Symbiont <- factor(x4$Symbiont, level = c("Arthropoda","Chimeras","Environmental micro-eukaryotes","Fungi","Gut_protist",
                                             "Gut_worm","Insect_parasite","Mammalia","Non-gut_parasite",
                                             "Occasional_parasite","Plant","Unknown"))

colours =c("purple","orange","dodgerblue3","chocolate3","grey40","grey70","plum","gold","red","tomato","forestgreen","wheat")
x4$Site <- gsub("neg.*","Neg",x4$Site)
x4$Site <- factor(x4$Site, level=c("Lokoue","Romani","Maya","Moba","Neg"))
x4$Gorilla_name_SampleID <- paste0(x4$Gorilla_name," ", x4$SampleID)
plot1 <- ggplot(x4, aes(x=reorder(Gorilla_name_SampleID,Total_reads), y=value, fill=Symbiont)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid", size=0.01) + #height of the column equal the value, stacked
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=18),
        legend.text=element_blank(),
        legend.title=element_blank(),
        legend.position = "left",
        legend.spacing.y = unit(0.2, 'cm'),
        strip.text.x = element_text(size = 20)) +
  guides(fill=guide_legend(ncol=1,byrow = TRUE))
print(plot1)

#relative abundance
x1 <- as.data.frame(x %>% group_by(Symbiont) %>% 
                      summarise_each(funs(sum)))

row.names(x1) <- x1$Symbiont
x1$Symbiont <- NULL
x1[is.na(x1)] <- 0
x1 <- t(x1)
x1 <- x1/rowSums(x1)
x1 <- t(x1)

x2 <- melt(x1)
x3 <- merge(x2,Total_reads, by.x="Var2",by.y="row.names")
colnames(x3)[4] <- "Total_reads"
x4 <- merge(x3,metadata_spl_ctl_df, by.x="Var2",by.y="row.names")
x4$Var1 <- factor(x4$Var1, level = c("Arthropoda","Chimeras","Environmental micro-eukaryotes","Fungi","Gut_protist",
                                     "Gut_worm","Insect_parasite","Mammalia","Non-gut_parasite",
                                     "Occasional_parasite","Plant","Unknown"))

colours =c("purple","orange","dodgerblue3","chocolate3","grey40","grey70","plum","gold","red","tomato","forestgreen","wheat")
x4$Site <- gsub("neg.*","Neg",x4$Site)
x4$Site <- factor(x4$Site, level=c("Lokoue","Romani","Maya","Moba","Neg"))
x4$Gorilla_name_SampleID <- paste0(x4$Gorilla_name," ", x4$SampleID)
plot2 <- ggplot(x4, aes(x=reorder(Gorilla_name_SampleID,Total_reads), y=value, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid", size=0.01) + #height of the column equal the value, stacked
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_blank(),
        legend.title=element_blank(),
        strip.text.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  guides(fill=guide_legend(ncol=1,byrow = TRUE))
print(plot2)
dev.off()

#Combined plot together
combined_plot <- ggarrange(plot1, plot2, ncol = 2, widths = c(1.4,0.8))
ggsave(filename = "Sup_Figure_3_fecal_rel2.pdf", path = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/", plot= combined_plot, width = 17, height = 20, device = "pdf" )


########################################################

#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
###                                                                                                                                                                       ###
###                                                               ANALYSIS WITHOUT COMBINING FECAL SAMPLES                                                                ###
###                                                                                                                                                                       ###
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################

########################################################
## Extract taxonomy, ASV and control tables from data ##
########################################################

#Select only gut eukaryote
asv_gut_rnk_spl_neg_df2 <- subset(asv_rnk_spl_neg_df2, Symbiont %in% c("Gut_worm","Gut_protist"))

### Extract taxonomy table ###
names_rnk <- c("Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Rank8","Accession_number")
asv_tax.df <- asv_gut_rnk_spl_neg_df2[,names_rnk]
#Because of presence of NA, need to rename unknwown taxa at each rank (once printing makes easier selection of taxa to create table)
asvUn_tax.df <- asv_tax.df                          #Create new table  
asvUn_tax.df[is.na(asvUn_tax.df)] <- "Un"
#create a function to paste "Un" to the name of the previous rank
mypaste = function(s) {
  s = gsub("_Un","",s)
  return(paste0(s,"_Un"))
}
#apply the function to all ranks
for (i in 2:9){
  
  if (i == 9)
  {
    colName = "Accession_number"
    colNamePrevious = "Rank8"
  } else
  {
    colName = paste0("Rank", i)
    colNamePrevious = paste0("Rank", i-1)
  }
  w = which(asvUn_tax.df[[colName]] == "Un")
  asvUn_tax.df[[colName]][w] = mypaste(asvUn_tax.df[[colNamePrevious]][w])
}
asvUn_tax.mtx <- as.matrix(asvUn_tax.df);dim(asvUn_tax.mtx)

### Extract negative controls tables ###
names_neg <- c("neg1","neg2","neg3","neg4")
asv_neg.df <- asv_gut_rnk_spl_neg_df2[,names_neg]
class(asv_neg.df);head(asv_neg.df);dim(asv_neg.df);str(asv_neg.df)

### Extract ASV-Samples table ###
names_spl_only <- colnames(asv_gut_rnk_spl_neg_df2)[grep("G",colnames(asv_gut_rnk_spl_neg_df2))]
names_spl <- c("ASV_cluster",names_spl_only)
asv_spl.df <- asv_gut_rnk_spl_neg_df2[,names_spl]
asv_spl.df[,-1] <- asv_spl.df[,-1] %>% mutate_all(as.numeric)
asv_spl.df <- as.data.frame(asv_spl.df %>% group_by(ASV_cluster) %>% summarise_each(funs(sum)))
row.names(asv_spl.df) <- asv_spl.df$ASV_cluster
asv_spl.df$ASV_cluster <- NULL
head(asv_spl.df); dim(asv_spl.df) ; str(asv_spl.df)

########################################################

########################################################
###               Filtering steps                    ###
########################################################

#Create phyloseq object
metadata <- sample_data(metadata_spl_ctl_df)
otu <- otu_table(asv_spl.df, taxa_are_rows=TRUE)
tax <- tax_table(asvUn_tax.mtx)
project_data <- phyloseq(metadata,otu,tax);project_data

project_data <- subset_samples(project_data, SampleID != "G0295") #Select experiment_title

########################################################

################################################################################################################
###                              Prevalence and Abundance of intestinal eukaryotes                           ###
################################################################################################################

########################################################
###   Calculate Prevalence, and Abundance of reads   ###
########################################################

#Calculate total number of reads per sample (reads coverage)
asv_spl_df2 <- asv_rnk_spl_neg_df2[,grep("G",colnames(asv_rnk_spl_neg_df2))]
asv_spl_df2$G0009 <- as.numeric(asv_spl_df2$G0009)
asv_spl_df2[is.na(asv_spl_df2)] <- 0
asv_spl_df2 <- asv_spl_df2 %>% mutate_all(as.numeric)
Total_reads_MiSeq_run <- sum(asv_spl_df2)

#---- Create ASV table for ----#
asv_spl.df <- as.data.frame(otu_table(project_data) )#create ASV table
asv_tax.df <- as.data.frame(tax_table(project_data))
rownames(asv_spl.df)==rownames(asv_tax.df)           #Sanity check to see if ASV from both table are similar and same order
asv_tax_spl.df <- cbind(asv_tax.df,asv_spl.df)       #Bind both table to have ASV and TAX table together

#---- Build prevalence table at Rank8 ----#
asv_tax_spl_rnk8 <- asv_tax_spl.df[,setdiff(colnames(asv_tax_spl.df),c(names_rnk[-8]))]                 #Select table with samples and rank8 only
row.names(asv_tax_spl_rnk8) <- asv_tax_spl_rnk8$Rank8
asv_tax_spl_rnk8$Rank8 <- NULL
Total_Absolute_abundance <- rowSums(asv_tax_spl_rnk8)
#Calculate prevalence
asv_tax_spl_rnk8_bin <- asv_tax_spl_rnk8
asv_tax_spl_rnk8_bin[] <- +(asv_tax_spl_rnk8_bin >= 1) #Convert table into binary 1:0 (Presence absence)
Prevalence <- round(rowSums(asv_tax_spl_rnk8_bin)/length(colnames(asv_tax_spl_rnk8_bin))*100,2)
#Combine Prevalance, Total_abundance, Total_abundance_log to tax table
d1 <- data.frame(Rank8 = row.names(asv_tax_spl_rnk8), 
                 Prevalence = Prevalence,
                 Absolute_abundance = Total_Absolute_abundance)
d1$Rank8 <- as.character(d1$Rank8)
d1$Relative_GutEuk <- (d1$Absolute_abundance)/sum(d1$Prevalence)
d1$Relative_AllEuk <- (d1$Absolute_abundance)/Total_reads_MiSeq_run

ggplot(d1, aes(x=Prevalence))+
  geom_histogram() +
  scale_x_continuous(trans='log10')

gut_residents_categories <- data.frame(Rank8 = c("Strongylida_ASV1","Strongylida_ASV2","Tetratrichomonas_buttreyi","Cycloposthium_bipalmatum","Nematoda_sp._ASV29","Blastocystis_ST2","Parentodinium_sp._ASV76","Blastocystis_Unknown_ASV28","Parentodinium_sp._ASV132","Enteromonad_Unknown_ASV73","Metastrongyloidea_ASV152","Blastocystis_Unknown_ASV45","Blastocystis_ST3","Tetratrichomonas_prowazecki","Entamoeba_hartmanni","Strongyloides_fuelleborni","Spiruroidea_ASV58","Enteromonas_hominis","Strongylida_ASV154","Strongylida_ASV177","Enteromonad_Unknown_ASV663","Blastocystis_ST1","Entamoeba_coli_ST2","Charonina_ventriculi","Ascaridia_ASV389","Pseudoentodinium_elephantis","Metastrongyloidea_ASV482","Hemiprorodon_gymnoprosthium","Latteuria_media","Retortamonas_sp._CladeA","Blastocystis_ST11","Metastrongyloidea_ASV204","Entamoeba_coli_ST1","Trepomonas_Unknown_ASV294","Spiruroidea_ASV311","Abbreviata_ASV495","Triplumaria_fulgora","Helicozoster_indicus","Diplomonad_Unknown","Trepomonas_agilis"),
                                       Resident_status = c("Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Common gut residents","Occasional gut residents","Common gut residents","Common gut residents","Common gut residents","Common gut residents","Common gut residents","Common gut residents","Common gut residents","Occasional gut residents","Occasional gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents"))

d2 <- merge(d1,gut_residents_categories,by="Rank8")
d2 <- d2[order(-d2$Prevalence),]
d2$Resident_status <- factor(d2$Resident_status, level =c("Dominant gut residents","Common gut residents","Occasional gut residents","Non-gut residents"))

d2_melt <- melt(d2, id.vars = c("Rank8","Prevalence","Resident_status"), variable = "Abundance", na.rm=TRUE)

colours <- c("darkgreen","blue","red","violet")


pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/prevalence_abundance/scatter_prev_abun_nametaxa_logboth.pdf", height =  15, width = 15)
plot <- ggplot(d2_melt, aes(x=Prevalence, y=value, color=Resident_status)) +
  geom_point() +
  geom_text(data=d2_melt, aes(label=Rank8, vjust = 2), size=3) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_vline(xintercept=50, linetype='dashed', color = "orange") +
  geom_vline(xintercept=10, linetype='dashed', color = "chocolate") +
  scale_color_manual(values=colours) +
  theme_light() +
  facet_grid(rows = vars(Abundance), scales = "free_y", switch = "y", space = "free_y") + 
  labs(y="Absolute and relative abundance (number of reads)", x="Prevalence of intestinal eukaryotes (%)") +
  theme(axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.position = 'bottom',
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        strip.text.y = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
print(plot)
dev.off()



d12 <- subset(d1, Rank8 !="Strongylida_ASV1" ) %>%
  subset(Rank8 !="Strongylida_ASV2" )

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/prevalence_abundance/scatter_prev_abun_nametaxa_noASV1andASV2.pdf", height =  10, width = 10)
plot <- ggplot(d12, aes(x=Prevalence, y=Total_abundance)) +
  geom_point() +
  geom_text(data=d12, aes(label=row.names(d12), vjust = 2))
print(plot)
dev.off()

d11 <- merge(d1,total_nb_per_sample,)


########################################################

########################################################
###           Calculate Relative abundance           ###
########################################################

#---- Create ASV table for ----#
asv_spl_df2 <- asv_rnk_spl_neg_df2[,grep("G",colnames(asv_rnk_spl_neg_df2))]
asv_spl_df2$G0009 <- as.numeric(asv_spl_df2$G0009)
asv_spl_df2[is.na(asv_spl_df2)] <- 0
asv_spl_df2 <- asv_spl_df2 %>% mutate_all(as.numeric)
asv_spl_df2_rel <- (asv_spl_df2)/(colSums(asv_spl_df2))
colSums(asv_spl_df2_rel)
asv_tax_rank8 <- asv_rnk_spl_neg_df2[,c("Symbiont","Rank8")]
asv_rank8_rel_all <- merge(asv_tax_rank8,asv_spl_df2_rel,by="row.names")
asv_rank8_rel_gut <- subset(asv_rank8_rel_all, Symbiont %in% c("Gut_protist","Gut_worm"))
asv_rank8_rel_gut[,1:2] <- NULL
otu_tax_agg.df <- aggregate(asv_rank8_rel_gut, by = list(asv_rank8_rel_gut$Rank8), FUN = max)[,-1]
Relative_abundance_allEuk_agg_melt <- melt(otu_tax_agg.df,id.vars=c("Rank8"))
Relative_abundance_allEuk_agg_melt$Abundance <- "Relative_AllEuk"


#---- Create ASV table for ----#
asv_spl.df <- as.data.frame(otu_table(project_data) )#create ASV table
asv_tax.df <- as.data.frame(tax_table(project_data))
rownames(asv_spl.df)==rownames(asv_tax.df)           #Sanity check to see if ASV from both table are similar and same order
asv_tax_spl.df <- cbind(asv_tax.df,asv_spl.df)       #Bind both table to have ASV and TAX table together
asv_tax_spl_rnk8 <- asv_tax_spl.df[,setdiff(colnames(asv_tax_spl.df),c(names_rnk[-8]))]                 #Select table with samples and rank8 only
Absolute_abundance_agg_melt <- melt(asv_tax_spl_rnk8,id.vars=c("Rank8"))
Absolute_abundance_agg_melt$Abundance <- "Absolute_abundance"

#---- Calculate relative abundance
asv_tax_spl_rnk8 <- asv_tax_spl.df[,setdiff(colnames(asv_tax_spl.df),c(names_rnk[-8]))]                 #Select table with samples and rank8 only
row.names(asv_tax_spl_rnk8) <- asv_tax_spl_rnk8$Rank8
asv_tax_spl_rnk8$Rank8 <- NULL
asv_tax_spl_rnk8 <- as.data.frame(t(asv_tax_spl_rnk8))
asv_tax_spl_rnk8 <- asv_tax_spl_rnk8 %>% mutate_all(as.numeric)
asv_tax_spl_rnk8_rel <- asv_tax_spl_rnk8/rowSums(asv_tax_spl_rnk8)
asv_tax_spl_rnk8_rel_melt <- as.data.frame(t(asv_tax_spl_rnk8_rel))
asv_tax_spl_rnk8_rel_melt$Rank8 <- row.names(asv_tax_spl_rnk8_rel_melt) 
Relative_abundance_GutEuk_agg_melt <- melt(asv_tax_spl_rnk8_rel_melt,id.vars=c("Rank8"))
Relative_abundance_GutEuk_agg_melt$Abundance <- "Avsolute_GutEuk"

#Merge with prevalence data
d4 <- rbind(Relative_abundance_allEuk_agg_melt,Absolute_abundance_agg_melt,Relative_abundance_GutEuk_agg_melt)
d5 <- merge(d4,d1, by="Rank8")
gut_residents_categories <- data.frame(Rank8 = c("Strongylida_ASV1","Strongylida_ASV2","Tetratrichomonas_buttreyi","Cycloposthium_bipalmatum","Nematoda_sp._ASV29","Blastocystis_ST2","Parentodinium_sp._ASV76","Blastocystis_Unknown_ASV28","Parentodinium_sp._ASV132","Enteromonad_Unknown_ASV73","Metastrongyloidea_ASV152","Blastocystis_Unknown_ASV45","Blastocystis_ST3","Tetratrichomonas_prowazecki","Entamoeba_hartmanni","Strongyloides_fuelleborni","Spiruroidea_ASV58","Enteromonas_hominis","Strongylida_ASV154","Strongylida_ASV177","Enteromonad_Unknown_ASV663","Blastocystis_ST1","Entamoeba_coli_ST2","Charonina_ventriculi","Ascaridia_ASV389","Pseudoentodinium_elephantis","Metastrongyloidea_ASV482","Hemiprorodon_gymnoprosthium","Latteuria_media","Retortamonas_sp._CladeA","Blastocystis_ST11","Metastrongyloidea_ASV204","Entamoeba_coli_ST1","Trepomonas_Unknown_ASV294","Spiruroidea_ASV311","Abbreviata_ASV495","Triplumaria_fulgora","Helicozoster_indicus","Diplomonad_Unknown","Trepomonas_agilis"),
                                       Resident_status = c("Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Common gut residents","Occasional gut residents","Common gut residents","Common gut residents","Common gut residents","Common gut residents","Common gut residents","Common gut residents","Occasional gut residents","Occasional gut residents","Common gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents"))

d6 <- merge(d5,gut_residents_categories,by="Rank8")
d6 <- d6[order(-d6$Prevalence),]
d6$Resident_status <- factor(d6$Resident_status, level =c("Dominant gut residents","Common gut residents","Occasional gut residents","Non-gut residents"))
colours <- c("darkgreen","blue","red","violet")

d6$Rank8 <- gsub("_", " ", d6$Rank8)

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/prevalence_abundance/JitterBoxplot_relabun_prev_nametaxa_logboth.pdf", height =  15, width = 19)
plot <- ggplot(d6, aes(x=reorder(Rank8,Prevalence), y=value, color=Resident_status)) +
  geom_jitter(position=position_dodge(0.8), size = 1) +
  scale_y_continuous(trans='log10') +
  geom_vline(xintercept="Enteromonas_hominis", linetype='dashed', color = "brown") +
  geom_vline(xintercept="Parentodinium_sp._ASV76", linetype='dashed', color = "orange") +
  scale_color_manual(values=colours) +
  theme_light() +
  facet_grid(rows = vars(Abundance), scales = "free_y", switch = "y", space = "free_y") + 
  labs(y="Absolute and relative abundance (number of reads)",x="Intestinal eukaryote taxa") +
  theme(axis.text.x = element_text(size = 19, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.position = 'bottom',
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        strip.text.y = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
print(plot)
dev.off()

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/prevalence_abundance/JitterBoxplot_relabun_prev_nametaxa.pdf", height =  15, width = 19)
plot <- ggplot(d6, aes(x=reorder(Rank8,Prevalence), y=value, color=Resident_status)) +
  geom_jitter(position=position_dodge(0.8), size = 1) +
  geom_vline(xintercept="Enteromonas_hominis", linetype='dashed', color = "brown") +
  geom_vline(xintercept="Parentodinium_sp._ASV76", linetype='dashed', color = "orange") +
  scale_color_manual(values=colours) +
  theme_light() +
  facet_grid(rows = vars(Abundance), scales = "free_y", switch = "y", space = "free_y") + 
  labs(y="Absolute and relative abundance (number of reads)",x="Intestinal eukaryote taxa") +
  theme(axis.text.x = element_text(size = 19, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.position = 'bottom',
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        strip.text.y = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
print(plot)
dev.off()

########################################################

########################################################
###     Filter tree using only ASV from ASV table    ###
########################################################

asv_to_drop <- setdiff(new_tree$tip.label,intersect(row.names(asv_tax.df),new_tree$tip.label))
new_tree <- drop.tip(new_tree, asv_to_drop)
plot(new_tree)
asv_tax.df$ASV_cluster <- row.names(asv_tax.df)
asv_tax.df <- asv_tax.df %>% arrange(factor(ASV_cluster, levels = new_tree$tip.label))
new_tree$tip.label == asv_tax.df$ASV_cluster
new_tree$tip.label <- asv_tax.df$Rank8
plot(new_tree)

########################################################

########################################################
###                  Create figure                   ###
########################################################

#Create tree
g <- ggtree(new_tree, branch.length="none")+ 
  geom_tiplab(size=6)+
  ggplot2::xlim(0, 20)

#Create Prevalence barplot
p1 <- ggplot(d1, aes(x=Rank8, y=Prevalence)) + 
  geom_col(width = 0.5) + 
  coord_flip() + 
  theme_half_open() +
  background_grid() +
  scale_y_continuous(limits=c(0, 100), breaks=c(0, 25, 50, 75, 100)) +
  labs(y="Prevalence  (% of fecal samples)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

#Create boxplot for abundance distribution log-transformed
p2b_log <- ggplot(asv_tax_spl_rnk8_rel_melt, aes(x=Rank8, y=value)) + 
  geom_boxplot(size = 0.5, width = 0.5)+
  scale_y_continuous(trans='log10') +
  geom_jitter(position=position_dodge(0.8), size = 1) +
  coord_flip() +
  theme_half_open() +
  background_grid() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.text.y = element_blank())

#Create boxplot for abundance distribution log-transformed
p2b <- ggplot(asv_tax_spl_rnk8_rel_melt, aes(x=Rank8, y=value)) + 
  geom_boxplot(size = 0.4, width = 0.5) +
  geom_jitter(position=position_dodge(0.8), size = 1) +
  coord_flip() +
  theme_half_open() +
  background_grid() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.text.y = element_blank())

#Save figure in directory
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/prevalence_abundance/phylo_all_prevalence_fecal_samples_not_combined.pdf", height =  10, width = 10)
plot_phylogeny_prevalence <- p1 %>% insert_left(g, width = 2) 
print(plot_phylogeny_prevalence)
dev.off()

#Save figure in directory
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/prevalence_abundance/phylo_all_rel_abun_logscale_fecal_samples_not_combined.pdf", height =  10, width = 13)
plot_phylogeny_rel_abun_logscale <- p2b_log %>% insert_left(g, width = 2) 
print(plot_phylogeny_rel_abun_logscale)
dev.off()

#Save figure in directory
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/prevalence_abundance/phylo_all_rel_abun_fecal_samples_not_combined.pdf", height =  10, width = 13)
plot_phylogeny_rel_abun <- p2b %>% insert_left(g, width = 2) 
print(plot_phylogeny_rel_abun)
dev.off()

plot_phylo_prev_rel_abun <- ggarrange(plot_phylogeny_prevalence,plot_phylogeny_rel_abun, ncol=2)
ggsave(filename = "Figure1_fecal_V1.pdf", path = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/prevalence_abundance/", plot= plot_phylo_prev_rel_abun, width = 30, height = 15, device = "pdf" )

########################################################

########################################################
###          Figure for committee meeting            ###
########################################################


gut_residents_categories <- data.frame(Rank8 = c("Strongylida_ASV1","Strongylida_ASV2","Tetratrichomonas_buttreyi","Cycloposthium_bipalmatum","Nematoda_sp._ASV29","Blastocystis_ST2","Parentodinium_sp._ASV76","Blastocystis_Unknown_ASV28","Parentodinium_sp._ASV132","Enteromonad_Unknown_ASV73","Metastrongyloidea_ASV152","Blastocystis_Unknown_ASV45","Blastocystis_ST3","Tetratrichomonas_prowazecki","Entamoeba_hartmanni","Strongyloides_fuelleborni","Spiruroidea_ASV58","Enteromonas_hominis","Strongylida_ASV154","Strongylida_ASV177","Enteromonad_Unknown_ASV663","Blastocystis_ST1","Entamoeba_coli_ST2","Charonina_ventriculi","Ascaridia_ASV389","Pseudoentodinium_elephantis","Metastrongyloidea_ASV482","Hemiprorodon_gymnoprosthium","Latteuria_media","Retortamonas_sp._CladeA","Blastocystis_ST11","Metastrongyloidea_ASV204","Entamoeba_coli_ST1","Trepomonas_Unknown_ASV294","Spiruroidea_ASV311","Abbreviata_ASV495","Triplumaria_fulgora","Helicozoster_indicus","Diplomonad_Unknown","Trepomonas_agilis"),
                                       Resident_status = c("Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Dominant gut residents","Common gut residents","Occasional gut residents","Common gut residents","Common gut residents","Common gut residents","Common gut residents","Common gut residents","Common gut residents","Common gut residents","Occasional gut residents","Occasional gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents","Non-gut residents"))

d2 <- merge(d1,gut_residents_categories,by="Rank8")
d2 <- d2[order(-d2$Prevalence),]
d2$Resident_status <- factor(d2$Resident_status, level =c("Dominant gut residents","Common gut residents","Occasional gut residents","Non-gut residents"))

colours <- c("darkgreen","blue","red","violet")

d2$Rank8 <- gsub("_"," ",d2$Rank8)

#Create Prevalence barplot
p1 <- ggplot(d2, aes(x=reorder(Rank8,Prevalence), y=Prevalence)) + 
  geom_col(width = 0.5) +
  geom_hline(yintercept=50, linetype='dashed', color = "orange") +
  geom_hline(yintercept=10, linetype='dashed', color = "chocolate") +
  coord_flip() + 
  theme_half_open() +
  background_grid() +
  scale_y_continuous(limits=c(0, 100), breaks=c(0, 25, 50, 75, 100)) +
  labs(y="Prevalence  (% of fecal samples)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank())


d7 <- Relative_abundance_GutEuk_agg_melt
d8 <- merge(d7,d1,by="Rank8")
p3 <- ggplot(d8, aes(x=value, y=reorder(Rank8,Prevalence))) +
  geom_boxplot() + 
  theme_half_open() +
  background_grid() +
  labs(y="Abundance (Number of reads)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
  
 


#Build sample table/Rank8
otu <- otu_table(project_data)
tax <- as.data.frame(tax_table(project_data))
otu_tax.df <- merge(otu,tax, by= "row.names")[,c("Rank8", colnames(otu))]
row.names(otu_tax.df) <- otu_tax.df$Rank8
otu_tax.df$Rank8 <- NULL
otu_tax.df <- t(otu_tax.df)
colnames(otu_tax.df)

#Add metadata with gorilla names, sampleID and dates
data <- as.data.frame(as.matrix(sample_data(project_data)))
data <- subset(data, Gorilla_name %in% c("Arthur","Kokima","Ribambel","Tino","Bangwetu","Diyo_L15","Minimee_L61","Murphy","Papillon","Riba","Robin_GS14","Vidole","Yugos"))

#Calculate number of day between collection of fecal samples per gorillas
metadata <- data[,c("Gorilla_name","Date")]
metadata$Date <- as.Date(metadata$Date, "%d-%m-%Y")
go_name <- unique(metadata$Gorilla_name)
stat_list <- list()
for (j in 1:length(go_name)) {
  table <- metadata[metadata$Gorilla_name==go_name[j],]
  date_recent <- sort(table$Date,decreasing = TRUE)[1]
  date_old <- sort(table$Date,decreasing = FALSE)[1]
  date_diff <- time_length(difftime(date_recent,date_old), "days")
  table_taxa <- data.frame(Gorilla_name = go_name[j],
                           Days_diff = date_diff)
  stat_list[[length(stat_list)+1]] <- table_taxa
}

Number_days_per_gorilla <- data.frame()
for (i in 1:length(stat_list)){
  Number_days_per_gorilla <- rbind(Number_days_per_gorilla,stat_list[[i]])
}


otu_tax_agg_meta.df <- merge(otu_tax.df, data, by="row.names", all.X=FALSE)
otu_tax_agg_meta.df <- otu_tax_agg_meta.df[,c(colnames(otu_tax.df),"Gorilla_name","SampleID","Date")]
otu_tax_agg_meta.df <- aggregate(otu_tax_agg_meta.df, by = list(otu_tax_agg_meta.df$Gorilla_name, otu_tax_agg_meta.df$Date), FUN=max)[,-c(1,2)]
otu_tax_agg_meta_melt.df <- melt(otu_tax_agg_meta.df, id.vars=c("Gorilla_name","SampleID","Date"), value.name = "Taxa")
otu_tax_agg_meta_melt.df$PresenceAbsence <- ifelse(otu_tax_agg_meta_melt.df$Taxa >0, "Presence","Absence")
otu_tax_agg_meta_melt.df$Date <- as.Date(otu_tax_agg_meta_melt.df$Date, "%d-%m-%Y")
otu_tax_agg_meta_melt.df$Connect <- paste0(otu_tax_agg_meta_melt.df$Gorilla_name,"_",otu_tax_agg_meta_melt.df$variable)

##Calculate the number of days between each collected fecal sample for each eukaryote and each gorilla
go_name <- unique(otu_tax_agg_meta_melt.df$Gorilla_name)
tax <- unique(otu_tax_agg_meta_melt.df$variable)
stat_list <- list()
for (i in 1:length(tax)){
  for (j in 1:length(go_name)) {
    table <- otu_tax_agg_meta_melt.df[otu_tax_agg_meta_melt.df$variable==tax[i],]
    table <- table[table$Gorilla_name==go_name[j],]
    x <- length(unique(table$PresenceAbsence))
    table_taxa <- data.frame(Taxa = tax[i],
                             Gorilla_name = go_name[j],
                             date_diff = x )
    stat_list[[length(stat_list)+1]] <- table_taxa
  }
}
summary_stat <- data.frame()
for (i in 1:length(stat_list)){
  summary_stat <- rbind(summary_stat,stat_list[[i]])
}
summary_stat$Connect <- paste0(summary_stat$Gorilla_name,"_",summary_stat$Taxa)

otu_tax_agg_meta_melt.df <- otu_tax_agg_meta_melt.df[,c("PresenceAbsence","Connect")]
test <- merge(summary_stat,otu_tax_agg_meta_melt.df,by="Connect")
test <- unique(test)
test <- merge(test,Number_days_per_gorilla,by="Gorilla_name")
test <- test %>% mutate(Status = case_when(
  date_diff == 2 ~ "Not_stay",
  date_diff == 1  & PresenceAbsence == "Absence" ~ "Absence",
  date_diff == 1  & PresenceAbsence == "Presence" ~ "Presence"
))


test <- subset(test, Status !="Absence")


# The palette with grey:
colours <- c("mediumorchid1",
             "navyblue",
             "yellow",
             "cyan",
             "gray33",
             "sienna1",
             "red",
             "forestgreen",
             "firebrick",
             "chartreuse2",
             "plum",
             "blue",
             "chocolate")

summary_stat_max_2 <- test

#diff_rank8 <- setdiff(d1$Rank8,summary_stat_max_2$Taxa)
#tablez <- as.data.frame(matrix(ncol = length(colnames(summary_stat_max_2)), nrow = length(diff_rank8)))
#colnames(tablez) <- colnames(summary_stat_max_2)
#tablez$Taxa <- diff_rank8

#taleg <- rbind(summary_stat_max_2,tablez)

#to_keep <- unique(summary_stat_max_2$Taxa)
#prevalence_taxa <- d1[d1$Rank8 %in% to_keep,]
tablex <-merge(summary_stat_max_2, d1, by.x = "Taxa",by.y="Rank8", all.y=TRUE)

d2 <- arrange(d2, -Prevalence)
tablex$Taxa <- gsub("_"," ",tablex$Taxa)
tablex$Taxa <- factor(tablex$Taxa, level = c(d2$Rank8))



pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/temporal_analysis/stay_not_stay.pdf", height =  8, width = 15)
p2 <- ggplot(tablex, aes(y = reorder(Taxa,Prevalence), x = Days_diff, shape=Gorilla_name)) + 
  geom_point(aes(color = as.factor(Status)), size = 5) +
  scale_colour_manual(values=colours) +
  scale_shape_manual(values=c(1:20)) +
  scale_x_continuous(trans='log10') +
  theme_half_open(line_size = 0.5) +
  background_grid() +
  labs(y="Days between samples collection") +
  theme(axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
print(p2)
dev.off()

plot_ggarrange <- ggarrange(p1,p3)
plot_ggarrange <- ggarrange(p1,p3,p2, ncol=3, widths = c(0.7,0.5,0.8)) #vince/phd/statistic/gorilla/figure/community_structure
ggsave(filename = "committee_meeting.pdf", path = "/Users/vincebilly/Desktop/", plot= plot_ggarrange, width = 14, height = 8, device = "pdf" )

plot_phylogeny_prevalence <- p1 %>% insert_left(p2, width = 2) 



########################################################


################################################################################################################
###                      Relative abundance for protist, eukaryote, and protist/eukaryote                    ###
################################################################################################################

########################################################
### Relative abundance within Protists Rank8
########################################################

sort(sample_sums(project_data))
#project_data_rarefied <- rarefy_even_depth(project_data, sample.size = 14418);project_data_rarefied
project_data_rarefied <- project_data
#Select ASV that belong to intestinal protist 
asv_protist.df <- subset(asv_rnk_spl_neg_df2, Symbiont %in% c("Gut_protist"))
asv_gut_protist <- row.names(asv_protist.df)
project_data_protist <- prune_taxa(asv_gut_protist, project_data_rarefied);project_data_protist

#Extract otu table and calculate relative abundance
otu <- t(otu_table(project_data_protist))
otu_rel <- otu/rowSums(otu)
otu_rel <- as.data.frame(otu_rel[,order(colSums(otu_rel), decreasing = TRUE)])
rowSums(otu_rel)  #sanity check, all value should be equal to 1
colSums(otu_rel)  #should be sorted from the highest to the lowest
otu_rel$sample_id <- row.names(otu_rel)
otu_rel_melt = melt(otu_rel, id = c("sample_id"))
otu_abs <- as.data.frame(otu)
otu_abs$sample_id <- row.names(otu_abs)
otu_abs_melt = melt(otu_abs, id = c("sample_id"))
otu_abs_rel_melt <-  merge(otu_rel_melt,otu_abs_melt,by=c("sample_id", "variable"))

#Extract taxonomy table
tax <- as.data.frame(tax_table(project_data_protist))
#Extract metadata table
data <- as.data.frame(as.matrix(sample_data(project_data_protist)))

#Merge tables (otu, taxand data)
otu_abs_rel_melt_tax <- merge(otu_abs_rel_melt,tax, by.x= "variable", by.y="row.names",all.x=TRUE, all.y=FALSE)
otu_abs_rel_melt_tax_data <- merge(otu_abs_rel_melt_tax,data,data, by.x="sample_id",by.y="row.names")

#Create column with presence absence of ASV
otu_abs_rel_melt_tax_data$value_bin <- ifelse(otu_abs_rel_melt_tax_data$value.x >0 , 1, 0) 

otu_abs_rel_melt_tax_data$gorilla_name_sampleID <- paste0(otu_abs_rel_melt_tax_data$Gorilla_name, "_", otu_abs_rel_melt_tax_data$SampleID)
otu_abs_rel_melt_tax_data$Rank8 <- gsub("_"," ",otu_abs_rel_melt_tax_data$Rank8 )
otu_abs_rel_melt_tax_data$Rank8 <- factor(otu_abs_rel_melt_tax_data$Rank8, 
                                          level= c("Blastocystis ST1","Blastocystis ST11","Blastocystis ST2","Blastocystis ST3","Blastocystis Unknown ASV28","Blastocystis Unknown ASV45",
                                                   "Diplomonad Unknown","Enteromonad Unknown ASV663","Enteromonad Unknown ASV73","Enteromonas hominis","Trepomonas agilis","Trepomonas Unknown ASV294","Retortamonas sp. CladeA",
                                                   "Charonina ventriculi","Cycloposthium bipalmatum","Helicozoster indicus","Hemiprorodon gymnoprosthium","Latteuria media","Parentodinium sp. ASV132","Parentodinium sp. ASV76","Pseudoentodinium elephantis","Triplumaria fulgora",
                                                   "Entamoeba coli ST1","Entamoeba coli ST2","Entamoeba hartmanni",
                                                   "Tetratrichomonas buttreyi","Tetratrichomonas prowazecki"))

otu_abs_rel_melt_tax_data$gorilla_name_sampleID <- paste0(otu_abs_rel_melt_tax_data$Gorilla_name,"_",otu_abs_rel_melt_tax_data$SampleID)
otu_abs_rel_melt_tax_data$Site <- factor(otu_abs_rel_melt_tax_data$Site, 
                                         level= c("Lokoue","Romani","Moba","Maya"))

colours <- c("purple","purple4","hotpink","pink","violet","palevioletred1","navy","dodgerblue","cyan","blue","lightcyan1","blue4","dodgerblue","olivedrab4","greenyellow","khaki1","aquamarine4","gold3","darkseagreen1","green4","seagreen2","yellow","grey10","grey50","grey60","coral4","peru")

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/Taxa_abundance_rank8_protist_fecal_sample.pdf", width = 15, height = 20)
Taxa_abundance_rank8_protist <- ggplot(otu_abs_rel_melt_tax_data, aes(x = reorder(gorilla_name_sampleID,Gorilla_name), y= value.x, fill=Rank8, color = colours)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid", size=0.01) + #height of the column equal the value, stacked
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        legend.position = "right",
        legend.spacing.y = unit(0.2, 'cm'),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  guides(fill=guide_legend(ncol=1,byrow = TRUE))
print(Taxa_abundance_rank8_protist)
dev.off()

########################################################

########################################################
### Nematoda community and Protist
########################################################

#Extract otu table and calculate relative abundance
otu <- t(otu_table(project_data_rarefied))
otu_rel <- otu/rowSums(otu)
otu_rel <- as.data.frame(otu_rel[,order(colSums(otu_rel), decreasing = TRUE)])
rowSums(otu_rel)  #sanity check, all value should be equal to 1
colSums(otu_rel)  #should be sorted from the highest to the lowest
otu_rel$sample_id <- row.names(otu_rel)
otu_rel_melt = melt(otu_rel, id = c("sample_id"))
otu_abs <- as.data.frame(otu)
otu_abs$sample_id <- row.names(otu_abs)
otu_abs_melt = melt(otu_abs, id = c("sample_id"))
otu_abs_rel_melt <-  merge(otu_rel_melt,otu_abs_melt,by=c("sample_id", "variable"))

#Extract taxonomy table
tax <- as.data.frame(tax_table(project_data_rarefied))
#Extract metadata table
data <- as.data.frame(as.matrix(sample_data(project_data_rarefied)))

#Merge tables (otu, taxand data)
otu_abs_rel_melt_tax <- merge(otu_abs_rel_melt,tax, by.x= "variable", by.y="row.names",all.x=TRUE, all.y=FALSE)
otu_abs_rel_melt_tax_data <- merge(otu_abs_rel_melt_tax,data,data, by.x="sample_id",by.y="row.names")

#Create column with presence absence of ASV
otu_abs_rel_melt_tax_data$value_bin <- ifelse(otu_abs_rel_melt_tax_data$value.x >0 , 1, 0) 

#Convert all protist lineAges into "protist" category
Protist_lineAges <- c("Blastocystis_ST1|Blastocystis_ST11|Blastocystis_ST2|Blastocystis_ST3|Blastocystis_Unknown_ASV28|Blastocystis_Unknown_ASV45|Diplomonad_Unknown|Enteromonad_Unknown_ASV663|Enteromonad_Unknown_ASV73|Enteromonas_hominis|Trepomonas_agilis|Trepomonas_Unknown_ASV294|Retortamonas_sp._CladeA|Charonina_ventriculi|Cycloposthium_bipalmatum|Helicozoster_indicus|Hemiprorodon_gymnoprosthium|Latteuria_media|Parentodinium_sp._ASV132|Parentodinium_sp._ASV76|Pseudoentodinium_elephantis|Triplumaria_fulgora|Entamoeba_coli_ST1|Entamoeba_coli_ST2|Entamoeba_hartmanni|Tetratrichomonas_buttreyi|Tetratrichomonas_prowazecki")
otu_abs_rel_melt_tax_data$Rank8 <- gsub(Protist_lineAges,"Protist",otu_abs_rel_melt_tax_data$Rank8)

otu_abs_rel_melt_tax_data$gorilla_name_sampleID <- paste0(otu_abs_rel_melt_tax_data$Gorilla_name, "_", otu_abs_rel_melt_tax_data$SampleID)
otu_abs_rel_melt_tax_data$Rank8 <- gsub("_"," ",otu_abs_rel_melt_tax_data$Rank8)
otu_abs_rel_melt_tax_data$Rank8 <- factor(otu_abs_rel_melt_tax_data$Rank8, 
                                          level= c("Protist","Strongylida ASV2","Strongylida ASV1","Strongylida ASV177","Strongylida ASV154",
                                                   "Metastrongyloidea ASV482","Metastrongyloidea ASV152","Metastrongyloidea ASV204",
                                                   "Strongyloides fuelleborni","Spiruroidea ASV58","Spiruroidea ASV311","Abbreviata ASV495","Ascaridia ASV389","Nematoda sp. ASV29"))


otu_abs_rel_melt_tax_data$Site <- factor(otu_abs_rel_melt_tax_data$Site, 
                                         level= c("Lokoue","Romani","Moba","Maya"))
otu_abs_rel_melt_tax_data$gorilla_name_sampleID <- gsub("_"," ", otu_abs_rel_melt_tax_data$gorilla_name_sampleID)

colours <- c("grey70","orange","red","tomato","red4","deepskyblue4","lightgoldenrod4","black","sienna1","aquamarine4","purple","forestgreen","white","dodgerblue")
otu_abs_rel_melt_tax_data <- otu_abs_rel_melt_tax_data %>% select(c("Rank8","gorilla_name_sampleID","Gorilla_name","Site","value.x"))
x <- otu_abs_rel_melt_tax_data %>% 
  group_by(Rank8,gorilla_name_sampleID,Gorilla_name,Site)  %>%
  summarise(value.x = sum(value.x))

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/Taxa_abundance_rank8_protist_nematode_fecal.pdf", width = 15, height = 20)
Taxa_abundance_rank8_protist_nematode <- ggplot(x, aes(x = reorder(gorilla_name_sampleID,Gorilla_name), y= value.x, fill=Rank8, color = colours)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid",size=0.1) + #height of the column equal the value, stacked
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=20),
        legend.spacing.y = unit(0.2, 'cm'),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position = "left",
        strip.text.x = element_text(size = 20)) +
  guides(fill = guide_legend(byrow = TRUE))
print(Taxa_abundance_rank8_protist_nematode)
dev.off()

########################################################

########################################################
### Nematoda and Protist
########################################################

sort(sample_sums(project_data))
project_data_rarefied <- rarefy_even_depth(project_data, sample.size = 14418);project_data

otu <- t(otu_table(project_data_rarefied))
otu_rel <- otu/rowSums(otu)
otu_rel <- as.data.frame(otu_rel[,order(colSums(otu_rel), decreasing = TRUE)])
rowSums(otu_rel)  #sanity check, all value should be equal to 1
colSums(otu_rel)  #should be sorted from the highest to the lowest
otu_rel$sample_id <- row.names(otu_rel)
otu_rel_melt = melt(otu_rel, id = c("sample_id"))
otu_abs <- as.data.frame(otu)
otu_abs$sample_id <- row.names(otu_abs)
otu_abs_melt = melt(otu_abs, id = c("sample_id"))
otu_abs_rel_melt <-  merge(otu_rel_melt,otu_abs_melt,by=c("sample_id", "variable"))

#Extract taxonomy table
tax <- as.data.frame(tax_table(project_data_rarefied))
#Extract metadata table
data <- as.data.frame(as.matrix(sample_data(project_data_rarefied)))

#Merge tables (otu, taxand data)
otu_abs_rel_melt_tax <- merge(otu_abs_rel_melt,tax, by.x= "variable", by.y="row.names",all.x=TRUE, all.y=FALSE)
otu_abs_rel_melt_tax_data <- merge(otu_abs_rel_melt_tax,data,data, by.x="sample_id",by.y="row.names")

#Create column with presence absence of ASV
otu_abs_rel_melt_tax_data$value_bin <- ifelse(otu_abs_rel_melt_tax_data$value.x >0 , 1, 0) 

otu_abs_rel_melt_tax_data$Rank4 <- gsub("Entamoebida|Blastocystis|Fornicata|Litostomatea|Parabasalia","Protist",otu_abs_rel_melt_tax_data$Rank4)
otu_abs_rel_melt_tax_data$gorilla_name_sampleID <- paste0(otu_abs_rel_melt_tax_data$gorilla_name,"_",otu_abs_rel_melt_tax_data$sampleID)
otu_abs_rel_melt_tax_data$Site <- factor(otu_abs_rel_melt_tax_data$Site, 
                                         level= c("Lokoue","Romani","Moba","Maya"))

colours  = c("grey69","grey28")

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/Taxa_abundance_rank4_protist_nematode_fecal_sample.pdf", width = 10, height = 10)
Taxa_abundance_rank4_protist_nematode <- ggplot(otu_abs_rel_melt_tax_data, aes(x = reorder(gorilla_name_sampleID,gorilla_name), y= value.x, fill=Rank4, color = colours)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid",size=0.1) + #height of the column equal the value, stacked
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=15),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  theme(legend.position="none")
print(Taxa_abundance_rank4_protist_nematode)
dev.off()

########################################################

########################################################
### Relative abundance within Nematodes Rank8
########################################################

# Select only worm
asv_worm.df <- subset(asv_rnk_spl_neg_df2, Symbiont %in% c("gut_worm"))
asv_gut_worm <- row.names(asv_worm.df)
project_data_worm <- prune_taxa(asv_gut_worm, project_data_rarefied);project_data_worm

# Buil relative abudance table
otu <- t(otu_table(project_data_worm))
otu_rel <- otu/rowSums(otu)
otu_rel <- as.data.frame(otu_rel[,order(colSums(otu_rel), decreasing = TRUE)])
rowSums(otu_rel)  #sanity check, all value should be equal to 1
colSums(otu_rel)  #should be sorted from the highest to the lowest
otu_rel$sample_id <- row.names(otu_rel)
otu_rel_melt = melt(otu_rel, id = c("sample_id"))
otu_abs <- as.data.frame(otu)
otu_abs$sample_id <- row.names(otu_abs)
otu_abs_melt = melt(otu_abs, id = c("sample_id"))
otu_abs_rel_melt <-  merge(otu_rel_melt,otu_abs_melt,by=c("sample_id", "variable"))

#Extract taxonomy table
tax <- as.data.frame(tax_table(project_data_worm))
#Extract metadata table
data <- as.data.frame(as.matrix(sample_data(project_data_worm)))

#Merge tables (otu, taxand data)
otu_abs_rel_melt_tax <- merge(otu_abs_rel_melt,tax, by.x= "variable", by.y="row.names",all.x=TRUE, all.y=FALSE)
otu_abs_rel_melt_tax_data <- merge(otu_abs_rel_melt_tax,data,data, by.x="sample_id",by.y="row.names")

#Create column with presence absence of ASV
otu_abs_rel_melt_tax_data$value_bin <- ifelse(otu_abs_rel_melt_tax_data$value.x >0 , 1, 0) 

otu_abs_rel_melt_tax_data$gorilla_name_sampleID <- paste0(otu_abs_rel_melt_tax_data$gorilla_name, "_", otu_abs_rel_melt_tax_data$sampleID)
otu_abs_rel_melt_tax_data$Rank8 <- gsub("_"," ",otu_abs_rel_melt_tax_data$Rank8)
otu_abs_rel_melt_tax_data$Rank8 <- factor(otu_abs_rel_melt_tax_data$Rank8, 
                                          level= c("Strongylida ASV2","Strongylida ASV1","Strongylida ASV177","Strongylida ASV154",
                                                   "Metastrongyloidea ASV482","Metastrongyloidea ASV152","Metastrongyloidea ASV204",
                                                   "Strongyloides fuelleborni","Spiruroidea ASV58","Spiruroidea ASV311","Abbreviata ASV495","Ascaridia ASV389","Nematoda sp. ASV29"))

otu_abs_rel_melt_tax_data$gorilla_name_sampleID <- paste0(otu_abs_rel_melt_tax_data$gorilla_name,"_",otu_abs_rel_melt_tax_data$sampleID)
otu_abs_rel_melt_tax_data$Site <- factor(otu_abs_rel_melt_tax_data$Site, 
                                         level= c("Lokoue","Romani","Moba","Maya"))

colours <- c("orange","red","tomato","red4","green","darkgreen","greenyellow","gold3","aquamarine4","yellow2","grey10","dodgerblue","grey70")

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/Taxa_abundance_rank8_worm_fecal_samples.pdf", width = 15, height = 12)
Taxa_abundance_rank8_worm <- ggplot(otu_abs_rel_melt_tax_data, aes(x = reorder(gorilla_name_sampleID,gorilla_name), y= value.x, fill=Rank8, color = colours)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid", size=0.01) + #height of the column equal the value, stacked
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=15),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.spacing.y = unit(0.2, 'cm'),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE))
print(Taxa_abundance_rank8_worm)
dev.off()

########################################################

########################################################
### Combine Plot into Figure 2
########################################################

plot_ggarrange <- ggarrange(Taxa_abundance_rank8_protist,Taxa_abundance_rank4_protist_nematode,Taxa_abundance_rank8_worm, ncol=3, widths = c(1.5,0.3,1.2))
ggsave(filename = "Figure2_fecal_V1.pdf", path = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/", plot= plot_ggarrange, width = 25, height = 15, device = "pdf" )

plot_ggarrange <- ggarrange(Taxa_abundance_rank8_protist_nematode,Taxa_abundance_rank8_protist, ncol=2, widths = c(1.3,1))
ggsave(filename = "Figure2_fecal_V2_2.pdf", path = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/", plot= plot_ggarrange, width = 25, height = 20, device = "pdf" )


########################################################

################################################################################################################
###                             RELATIONSHIP BETWEEN INTESINAL EUKARYOTES                                    ###
################################################################################################################

### Co-association (Presence/Absence)
########################################################

#Set working directory
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test")
sort(sample_sums(project_data))
#Create list for final table with
permute_list <- list()
permutation_vector <- rep(14418,20)
reads_threshold <- c(2,5,10,25,50)
min_sample_size <- c(5)
for(perm in 1:length(permutation_vector)){
  project_data_rarefied <- rarefy_even_depth(project_data, sample.size = permutation_vector[perm]);project_data
  #---- Create ASV table for ----#
  asv_spl.df <- as.data.frame(otu_table(project_data_rarefied) )#create ASV table
  asv_tax.df <- as.data.frame(tax_table(project_data_rarefied))
  rownames(asv_spl.df)==rownames(asv_tax.df)           #Sanity check to see if ASV from both table are similar and same order
  asv_tax_spl.df <- cbind(asv_tax.df,asv_spl.df)       #Bind both table to have ASV and TAX table together
  #---- Build prevalence table at Rank8 ----#
  for (nbreads in 1:length(reads_threshold)){
    asv_tax_spl_rnk8_before <- asv_tax_spl.df[,setdiff(colnames(asv_tax_spl.df),c(names_rnk[-8]))]                 #Select table with samples and rank8 only
    row.names(asv_tax_spl_rnk8_before) <- asv_tax_spl_rnk8_before$Rank8
    asv_tax_spl_rnk8_before$Rank8 <- NULL
    asv_tax_spl_rnk8_before <- as.data.frame(t(asv_tax_spl_rnk8_before))
    asv_tax_spl_rnk8_before[asv_tax_spl_rnk8_before<reads_threshold[nbreads]] <- 0
    asv_tax_spl_rnk8_before[asv_tax_spl_rnk8_before>=reads_threshold[nbreads]] <- 1
    #Only take taxa with at least 5 sample presence  and at least 5 sample absent 
    for (j in 1:length(min_sample_size)){
    to_keep <- which(colSums(asv_tax_spl_rnk8_before) < length(sample_names(project_data_rarefied))- min_sample_size[j] & colSums(asv_tax_spl_rnk8_before) > min_sample_size[j])
    asv_tax_spl_rnk8 <- asv_tax_spl_rnk8_before[,to_keep]
    asv_tax_spl_rnk8 <- ifelse(asv_tax_spl_rnk8 ==1, "Presence", "Absence")
    #Run fisher test for all pairwise comparison
    stat_list <- list()
    for (taxa1 in colnames(asv_tax_spl_rnk8)){
      for (taxa2 in colnames(asv_tax_spl_rnk8)) {
        table_fisher <- table(asv_tax_spl_rnk8[,taxa1],asv_tax_spl_rnk8[,taxa2])
        is_fisher_test <-  
          tryCatch(expr =  { 
            fisher_test <- fisher.test(table_fisher)
            TRUE},
            error=function(e){
              return(FALSE)}
          )
        if(is_fisher_test)
          c(summary_table <- data.frame(Reads_Threshold = reads_threshold[nbreads],
                                        Iterate =perm,
                                        Min_Sample_Size = min_sample_size[j], 
                                        Taxa1=taxa1,
                                        Spl_taxa1= length(which(asv_tax_spl_rnk8[,taxa1]=="Presence")),
                                        Taxa2=taxa2,
                                        Spl_taxa2= length(which(asv_tax_spl_rnk8[,taxa2]=="Presence")),
                                        Pvalue=fisher_test$p.value,
                                        Odds_ratio =round(fisher_test$estimate,3)
          )
          )
        if(!is_fisher_test)
          c(summary_table <- data.frame(Reads_Threshold = reads_threshold[nbreads],
                                        Iterate =perm,
                                        Min_Sample_Size = min_sample_size[j],
                                        Taxa1=taxa1,
                                        Spl_taxa1= length(which(asv_tax_spl_rnk8[,taxa1]=="Presence")),
                                        Taxa2=taxa2,
                                        Spl_taxa2= length(which(asv_tax_spl_rnk8[,taxa2]=="Presence")),
                                        Pvalue=NA,
                                        Odds_ratio =NA
          )
          )
        stat_list[[length(stat_list)+1]] <- summary_table
      }
    }
    #Generate table that summarize all alpha and beta statistic
    Fisher_test_summary <- data.frame()
    for (i in 1:length(stat_list)){
      Fisher_test_summary <- rbind(Fisher_test_summary,stat_list[[i]])
      }
    Fisher_test_summary$Pvalue_Adjust <- round(p.adjust(Fisher_test_summary$Pvalue, method = "BH", n = length(Fisher_test_summary$Pvalue)),3)
    Fisher_test_summary$duplicate <- Fisher_test_summary$Spl_taxa1==Fisher_test_summary$Spl_taxa2
    Fisher_test_summary <- subset(Fisher_test_summary, duplicate == FALSE)
    Fisher_test_summary$Odds_sup1 <- ifelse(Fisher_test_summary$Odds_ratio>=1, ">1 (Positive association)","<1 (Negative association)")
    Fisher_test_summary$Significance <- ifelse(Fisher_test_summary$Pvalue_Adjust < 0.05,"Significant","Non-significant")
    row.names(Fisher_test_summary) <- NULL
    temp <- Fisher_test_summary[,c("Taxa1","Taxa2")]
    newDf <- data.frame(t(apply(temp,1,sort)))
    Fisher_test_summary <- Fisher_test_summary[!duplicated(newDf),]
    #If no single significant taxa, do not add to list
    if(dim(Fisher_test_summary)[1]>0){
      permute_list[[length(permute_list)+1]] <- Fisher_test_summary}
    #Plot histogram
    hist_plot <- ggplot(Fisher_test_summary, aes(x=Odds_ratio)) +
      geom_histogram(binwidth = 0.2)+
      geom_vline(xintercept=1, linetype='dashed') +
      theme_light() +
      theme(axis.title.y = element_text(size=15),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=15),
            axis.text.x = element_text(size=15))
    #Plot Bar plot number of odd ratio below and above 1
    bar_plot <- ggplot(Fisher_test_summary, aes(x=Odds_sup1)) +
      geom_bar() +
      theme_light() +
      theme(axis.title.y = element_text(size=15),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=15),
            axis.text.x = element_text(size=9))
    #Scatter plot, Odds ratio and Pvalue
    bubble_plot <- ggplot(Fisher_test_summary,aes(x=Odds_ratio, y=Pvalue_Adjust, color=Significance)) +
      geom_point() +
      xlim(0,20) +
      ylim(1,0.0001) +
      theme_light() +
      geom_vline(xintercept=1, linetype='dashed') +
      geom_hline(yintercept=0.05, linetype='dashed') +
      theme(axis.title.y = element_text(size=15),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=15),
            axis.text.x = element_text(size=15))
    combined <- ggarrange(
     bubble_plot,                # First row with line plot
     widths = c(1, 0.2),
    # Second row with box and histogram
    ggarrange(bar_plot, hist_plot, ncol = 2, labels = c("B", "C"), widths = c(0.7, 1)), 
      nrow = 2, 
      labels = "A"       # Label of the line plot
  )
  ggsave(combined,filename = paste0(perm,"_",reads_threshold[i],"_","fisher_test_multipannel.pdf"), width = 10, height = 10)
  }
  }
}


#Convert list into table
Fisher_permute_summary <- data.frame()
for (i in 1:length(permute_list)){
  Fisher_permute_summary <- rbind(Fisher_permute_summary,permute_list[[i]])
}

#Add additional information to table
Fisher_permute_summary$nb_significant <- 1 #To later calculate number of time pairwise was significant
Fisher_permute_summary$pair <- paste0(Fisher_permute_summary$Taxa1,"_",Fisher_permute_summary$Taxa2)
Fisher_permute_summary$pairwise <- ifelse(Fisher_permute_summary$Significance=="Significant",paste0(Fisher_permute_summary$Taxa1,"_",Fisher_permute_summary$Taxa2),"NS")
to_keep <- unique(Fisher_permute_summary$pairwise)
to_keep <- to_keep[!to_keep %in% "NS"]
Fisher_permute_summary_sg <- subset(Fisher_permute_summary, pair %in% to_keep)
Fisher_permute_summary_sg$pairwise <- Fisher_permute_summary_sg$pair
to_keep2 <- setdiff(Fisher_permute_summary$pair, to_keep)
Fisher_permute_summary_nsg <- subset(Fisher_permute_summary, pair %in% to_keep2)
Fisher_permute_summary <- rbind(Fisher_permute_summary_sg,Fisher_permute_summary_nsg)

#Plot histogram
hist_plot <- ggplot(Fisher_permute_summary, aes(x=Odds_ratio)) +
  geom_histogram(binwidth = 0.2)+
  geom_vline(xintercept=1, linetype='dashed') +
  theme_light() +
  facet_grid(rows = vars(Reads_Threshold), scales = "free_y", switch = "y", space = "free_y")+ 
  theme(axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=15))
plot(hist_plot)

#Plot Bar plot number of odd ratio below and above 1
bar_plot <- ggplot(Fisher_permute_summary, aes(x=Odds_sup1)) +
  geom_bar() +
  theme_light() +
  theme(axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=9))
print(bar_plot)

#Bubble plot for proportion of time test is significant
Fisher_permute_summary_sig <- subset(Fisher_permute_summary, Pvalue_Adjust <0.05) # Select only signiifcant adjusted pvalue
Fisher_permute_summary_groupby <- Fisher_permute_summary_sig %>% 
  group_by(Taxa1,Taxa2) %>%
  summarise(avg_adpvl = mean(Pvalue_Adjust),    #Calculate mean of pvalue 
            Mean_OddRatio = mean(Odds_ratio),   #Calculate mean of odds ratio
            Permutation = sum(nb_significant))  #Calculate total number of time the pairwise is significant
Fisher_permute_summary_groupby$Permutation <- Fisher_permute_summary_groupby$Permutation/perm*100
setorder(Fisher_permute_summary_groupby, cols= - "Permutation") #Reoerder table by the number of time it was significant
order_legend_pairwise <- paste0(Fisher_permute_summary_groupby$Taxa1,"_",Fisher_permute_summary_groupby$Taxa2) #Create vector to reorder scatter plot color legend
Fisher_permute_summary_groupby$Taxa1 <- gsub("_"," ",Fisher_permute_summary_groupby$Taxa1) #Remove underscore for aesthetic
Fisher_permute_summary_groupby$Taxa2 <- gsub("_"," ",Fisher_permute_summary_groupby$Taxa2) #Remove underscore for aesthetic
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test/bubble_plot_fisher_test.pdf", height =  5, width = 10)
plot <- ggplot(Fisher_permute_summary_groupby, aes(x=Taxa1, y=Taxa2, color=Mean_OddRatio, size=Permutation)) +
  geom_point() +
  scale_size(range = c(0, 15)) +
  theme_light() +
  scale_fill_gradient(low = "royalblue4", high = "red") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=15, hjust = 1, angle = 20, ))
print(plot)
dev.off()

#Scatter plot, Odds ratio and Pvalue
colours <- c("grey75","dodgerblue","coral1","chartreuse4","mAgenta","grey14","seagreen1","darkorange","purple","forestgreen","firebrick","aquamarine","lightpink","grey14","seagreen1","darkorange","purple","forestgreen","firebrick","aquamarine","lightpink","purple","forestgreen","firebrick","aquamarine","lightpink","grey14","seagreen1","darkorange","purple","forestgreen","firebrick","aquamarine","lightpink","dodgerblue")
Fisher_permute_summary$pairwise <- factor(Fisher_permute_summary$pairwise, level=c("NS",order_legend_pairwise))
Fisher_permute_summary$Reads_Threshold <- paste0(Fisher_permute_summary$Reads_Threshold," reads")
Fisher_permute_summary$Reads_Threshold <- factor(Fisher_permute_summary$Reads_Threshold, level = c("2 reads","5 reads","10 reads","25 reads","50 reads"))
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test/Scatter_fisher_20perm_MultiThreshold_fecal_5spl_legend.pdf", height =  40, width = 60)
scatter_plot <- ggplot(Fisher_permute_summary,aes(x=Odds_ratio, y=Pvalue_Adjust, color=pairwise, shape=pairwise)) +
  geom_point(size=5) +
  scale_shape_manual(values=c(1:31))+
  xlim(0,20) +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values=colours) +
  theme_light() +
  geom_vline(xintercept=1, linetype='dashed', colour = "red") +
  geom_text(aes(x=2.3, label="\nOdds Ratio = 1", y=0.002), colour="red", angle=0, size=17) +
  geom_hline(yintercept=0.05, linetype='dashed', colour = "orange") +
  geom_text(aes(y=0.2, label="\nAdj.Pvalue = 0.05", x=17.5), colour="orange", angle=0, size=17) +
  facet_grid(rows = vars(Reads_Threshold), scales = "free_y", switch = "y", space = "free_y")+ 
  labs(x="Odds Ratio", y="Adjusted p-value") +
  theme(axis.title.y = element_text(size=40),
        axis.text.y = element_text(size=40),
        axis.title.x = element_text(size=40),
        axis.text.x = element_text(size=40),
        legend.text = element_text(size=28),
        strip.text.y = element_text(size=40),
        legend.position = 'bottom') +
  guides(col=guide_legend(nrow=10))
print(scatter_plot)
dev.off()

#Bubble plot significant pairwise across Reads_Threshold and Iterations
Fisher_permute_summary_remove_sign <- subset(Fisher_permute_summary, Significance == "Significant")
Fisher_permute_summary_remove_sign$Min_Sample_Size <- as.factor(Fisher_permute_summary_remove_sign$Min_Sample_Size)
Fisher_permute_summary_remove_sign$Min_Sample_Size_Reads_Threshold <- paste0(Fisher_permute_summary_remove_sign$Min_Sample_Size,"_",Fisher_permute_summary_remove_sign$Reads_Threshol)
Fisher_permute_summary_remove_sign$Min_Sample_Size_Reads_Threshold <- factor(Fisher_permute_summary_remove_sign$Min_Sample_Size_Reads_Threshold, level = c("5_1","5_2","5_5","5_10","5_25","5_50","5_100","5_200","10_500","5_500","10_1","10_2","10_5","10_10","10_25","10_50","10_100","10_200"))
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test/bubble_fisher_20perm_MultiThreshold_5-10spl_fecal.pdf", height =  15, width = 35)
Bubble_plot <- ggplot(Fisher_permute_summary_remove_sign,aes(x=Iterate, y=pairwise, color=Odds_sup1, size=Odds_ratio)) +
  geom_point() +
  theme_light() +
  facet_grid(cols = vars(Min_Sample_Size_Reads_Threshold), scales = "free_x", switch = "x", space = "free_x") + 
  theme(axis.title.y = element_text(size=20),
      axis.text.y = element_text(size=20),
      axis.title.x = element_text(size=20),
      axis.text.x = element_text(size=20),
      legend.text = element_text(size=14),
      legend.position = 'right')
plot(Bubble_plot)
dev.off()



combined <- ggarrange(
  bubble_plot,                # First row with line plot
  # Second row with box and histogram
  ggarrange(plot,bar_plot,hist_plot, widths = c(0.6,0.2,0.3), ncol=3),
  nrow = 2, 
  labels = "A",       # Label of the line plot
  heights = c(1,0.5),
  widths = c(1,1)
)
ggsave(combined,filename = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test/fisher_test_multipannel.pdf", width = 30, height = 20)

write.table(Fisher_permute_summary,"/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/bubble_plot_fisher_test.txt", sep="\t", row.names = FALSE)
write.table(Fisher_test_summary,"/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test/fisher_test.txt", sep ="\t" )
Fisher_permute_summary <- fread("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/bubble_plot_fisher_test.txt")


#Understand Infinite value
x <- subset(Fisher_permute_summary, Odds_ratio == "Inf")
x2 <- subset(x, Pvalue_Adjust<0.05 )
taxa1_taxa2_cont_table <- table(asv_tax_spl_rnk8[,"Parentodinium_sp._ASV132"],asv_tax_spl_rnk8[,"Cycloposthium_bipalmatum"])
fisher.test(taxa1_taxa2_cont_table)
dimnames(taxa1_taxa2_cont_table) <- list(Parentodinium_sp._ASV132 = c("Absence_Taxa1", "Presence_Taxa1"),
                                         Cycloposthium_bipalmatum = c("Absence_Taxa2", "Presence_Taxa2"))
table <- as.data.frame(taxa1_taxa2_cont_table)
colours <- c("orangered1","blue")
plot <- ggplot(table, aes(x=Cycloposthium_bipalmatum, y=Freq, fill=Parentodinium_sp._ASV132, color = colours )) +
  geom_bar(stat="identity", width=1, color="black", linetype = "solid") +  #height of the column equal the value, stacked
  scale_fill_manual(values=colours) +
  theme_bw()+
  labs(x ="Cycloposthium_bipalmatum" ) +
  #ggtitle(paste0("Total:",total_sample,", Sens:",sensitivity, ", Spec:",specificity,", Accu:",accuracy,", Pval(Fisher):",fisher_Pval))+
  theme(legend.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        plot.title=element_text(size=8,face="bold"))

#
x <- subset(Fisher_permute_summary, pairwise == "Blastocystis_Unknown_ASV45_Enteromonad_Unknown_ASV73")
x2 <- subset(x, Pvalue_Adjust>0.2 )

ggplot(x, aes(x=Pvalue_Adjust, y=Spl_taxa1)) +
  geom_point() +
  geom_vline(xintercept=0.05, linetype='dashed')

test <- cor.test(x$Pvalue_Adjust,x$Spl_taxa1)
summary(test)

taxa1_taxa2_cont_table <- table(asv_tax_spl_rnk8[,"Blastocystis_Unknown_ASV45"],asv_tax_spl_rnk8[,"Enteromonad_Unknown_ASV73"])
dimnames(taxa1_taxa2_cont_table) <- list(Taxa1 = c("Absence_Taxa1", "Presence_Taxa1"),
                                         Taxa2 = c("Absence_Taxa2", "Presence_Taxa2"))

colours <- c("orangered1","blue")
plot <- ggplot(taxa1_taxa2_cont_table, aes(x=Taxa1, y=Freq, fill=Taxa2, color = colours )) +
  geom_bar(stat="identity", width=1, color="black", linetype = "solid") +  #height of the column equal the value, stacked
  scale_fill_manual(values=colours) +
  theme_bw()+
  labs(y ="Taxa2" , x ="Taxa1" ) +
  #ggtitle(paste0("Total:",total_sample,", Sens:",sensitivity, ", Spec:",specificity,", Accu:",accuracy,", Pval(Fisher):",fisher_Pval))+
  theme(legend.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        plot.title=element_text(size=8,face="bold"))



########################################################

### Heatmap with FISHER
########################################################

#Set working directory
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test")
sort(sample_sums(project_data))
#Create list for final table with
permute_list <- list()
permutation_vector <- rep(14418,1)
reads_threshold <- c(10)
min_sample_size <- c(5)
for(perm in 1:length(permutation_vector)){
  project_data_rarefied <- rarefy_even_depth(project_data, sample.size = permutation_vector[perm]);project_data
  #---- Create ASV table for ----#
  asv_spl.df <- as.data.frame(otu_table(project_data_rarefied) )#create ASV table
  asv_tax.df <- as.data.frame(tax_table(project_data_rarefied))
  rownames(asv_spl.df)==rownames(asv_tax.df)           #Sanity check to see if ASV from both table are similar and same order
  asv_tax_spl.df <- cbind(asv_tax.df,asv_spl.df)       #Bind both table to have ASV and TAX table together
  #---- Build prevalence table at Rank8 ----#
  for (nbreads in 1:length(reads_threshold)){
    asv_tax_spl_rnk8_before <- asv_tax_spl.df[,setdiff(colnames(asv_tax_spl.df),c(names_rnk[-8]))]                 #Select table with samples and rank8 only
    row.names(asv_tax_spl_rnk8_before) <- asv_tax_spl_rnk8_before$Rank8
    asv_tax_spl_rnk8_before$Rank8 <- NULL
    asv_tax_spl_rnk8_before <- as.data.frame(t(asv_tax_spl_rnk8_before))
    asv_tax_spl_rnk8_before[asv_tax_spl_rnk8_before<reads_threshold[nbreads]] <- 0
    asv_tax_spl_rnk8_before[asv_tax_spl_rnk8_before>=reads_threshold[nbreads]] <- 1
    #Only take taxa with at least 5 sample presence  and at least 5 sample absent 
    for (j in 1:length(min_sample_size)){
      to_keep <- which(colSums(asv_tax_spl_rnk8_before) < length(sample_names(project_data_rarefied))- min_sample_size[j] & colSums(asv_tax_spl_rnk8_before) > min_sample_size[j])
      asv_tax_spl_rnk8 <- asv_tax_spl_rnk8_before[,to_keep]
      asv_tax_spl_rnk8 <- ifelse(asv_tax_spl_rnk8 ==1, "Presence", "Absence")
      #Run fisher test for all pairwise comparison
      stat_list <- list()
      for (taxa1 in colnames(asv_tax_spl_rnk8)){
        for (taxa2 in colnames(asv_tax_spl_rnk8)) {
          table_fisher <- table(asv_tax_spl_rnk8[,taxa1],asv_tax_spl_rnk8[,taxa2])
          is_fisher_test <-  
            tryCatch(expr =  { 
              fisher_test <- fisher.test(table_fisher)
              TRUE},
              error=function(e){
                return(FALSE)}
            )
          if(is_fisher_test)
            c(summary_table <- data.frame(Reads_Threshold = reads_threshold[nbreads],
                                          Iterate =perm,
                                          Min_Sample_Size = min_sample_size[j], 
                                          Taxa1=taxa1,
                                          Spl_taxa1= length(which(asv_tax_spl_rnk8[,taxa1]=="Presence")),
                                          Taxa2=taxa2,
                                          Spl_taxa2= length(which(asv_tax_spl_rnk8[,taxa2]=="Presence")),
                                          Pvalue=fisher_test$p.value,
                                          Odds_ratio =round(fisher_test$estimate,9)
            )
            )
          if(!is_fisher_test)
            c(summary_table <- data.frame(Reads_Threshold = reads_threshold[nbreads],
                                          Iterate =perm,
                                          Min_Sample_Size = min_sample_size[j],
                                          Taxa1=taxa1,
                                          Spl_taxa1= length(which(asv_tax_spl_rnk8[,taxa1]=="Presence")),
                                          Taxa2=taxa2,
                                          Spl_taxa2= length(which(asv_tax_spl_rnk8[,taxa2]=="Presence")),
                                          Pvalue=NA,
                                          Odds_ratio =NA
            )
            )
          stat_list[[length(stat_list)+1]] <- summary_table
        }
      }
      #Generate table that summarize all alpha and beta statistic
      Fisher_test_summary <- data.frame()
      for (i in 1:length(stat_list)){
        Fisher_test_summary <- rbind(Fisher_test_summary,stat_list[[i]])
      }
      Fisher_test_summary$Pvalue_Adjust <- round(p.adjust(Fisher_test_summary$Pvalue, method = "BH", n = length(Fisher_test_summary$Pvalue)),3)
      Fisher_test_summary$Odds_sup1 <- ifelse(Fisher_test_summary$Odds_ratio>=0, ">0 (Positive association)","<0 (Negative association)")
      Fisher_test_summary$Significance <- ifelse(Fisher_test_summary$Pvalue_Adjust < 0.05,"Significant","Non-significant")
      row.names(Fisher_test_summary) <- NULL
    }
  }
}

#Heatmapt at 10 reads
Fisher_permute_summary_heatmap <- Fisher_test_summary
Fisher_permute_summary_heatmap$Odds_ratio_transformed <- ifelse(Fisher_permute_summary_heatmap$Odds_ratio=="Inf",NA,Fisher_permute_summary_heatmap$Odds_ratio)
Fisher_permute_summary_heatmap$Odds_ratio_transformed <- ifelse(Fisher_permute_summary_heatmap$Odds_ratio_transformed==0,NA,Fisher_permute_summary_heatmap$Odds_ratio_transformed)
Fisher_permute_summary_heatmap$Odds_ratio_transformed <- ifelse(Fisher_permute_summary_heatmap$Odds_ratio_transformed<1,-1/Fisher_permute_summary_heatmap$Odds_ratio_transformed,Fisher_permute_summary_heatmap$Odds_ratio_transformed)
Fisher_permute_summary_heatmap$Odds_ratio_transformed <- ifelse(Fisher_permute_summary_heatmap$Taxa1==Fisher_permute_summary_heatmap$Taxa2,NA,Fisher_permute_summary_heatmap$Odds_ratio_transformed)

Fisher_permute_summary_heatmap$r_if_sig <- ifelse(Fisher_permute_summary_heatmap$Significance=="Significant",Fisher_permute_summary_heatmap$Odds_ratio_transformed,NA)
#Fisher_permute_summary_heatmap$r_if_sig <- ifelse(Fisher_permute_summary_heatmap$Taxa1==Fisher_permute_summary_heatmap$Taxa2,NA,Fisher_permute_summary_heatmap$Odds_ratio_transformed)

Fisher_permute_summary_heatmap$Taxa1 <- gsub("_"," ",Fisher_permute_summary_heatmap$Taxa1)
Fisher_permute_summary_heatmap$Taxa2 <- gsub("_"," ",Fisher_permute_summary_heatmap$Taxa2)

Fisher_permute_summary_heatmap <- Fisher_permute_summary_heatmap %>% arrange(Taxa2)


#---- Create heatmap plot with pearson signifiancnt value
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test/heatmap_correlation.pdf", height =  5, width = 8)
plot <- Fisher_permute_summary_heatmap %>% 
  ggplot(aes(Taxa1, Taxa2, fill=Odds_ratio_transformed, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, 
       fill = "Odds Ratio\n(Transformed*)") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#A63446",high="#0C6291") +
  geom_text(size=2.9) +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10))
print(plot)
dev.off()

temp <- Fisher_test_summary[,c("Taxa1","Taxa2")]
newDf <- data.frame(t(apply(temp,1,sort)))
Fisher_test_summary_half <- Fisher_test_summary[!duplicated(newDf),]
Fisher_test_summary_half <- Fisher_test_summary_half %>% arrange(Taxa2)


plot <- Fisher_test_summary_half %>% 
  ggplot(aes(Taxa1, Taxa2, fill=Odds_ratio, label=round(Odds_ratio,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, 
       fill = "Odds Ratio\n(Transformed*)") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#A63446",high="#0C6291") +
  geom_text(size=2.9) +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10))

########################################################

### Probabilistic approach 
########################################################


project_data_rarefied <- rarefy_even_depth(project_data, sample.size = 14418);project_data
#---- Create ASV table for ----#
asv_spl.df <- as.data.frame(otu_table(project_data_rarefied) )#create ASV table
asv_tax.df <- as.data.frame(tax_table(project_data_rarefied))
rownames(asv_spl.df)==rownames(asv_tax.df)           #Sanity check to see if ASV from both table are similar and same order
asv_tax_spl.df <- cbind(asv_tax.df,asv_spl.df)       #Bind both table to have ASV and TAX table together
#---- Build prevalence table at Rank8 ----#
for (nbreads in 1:length(reads_threshold)){
  asv_tax_spl_rnk8_before <- asv_tax_spl.df[,setdiff(colnames(asv_tax_spl.df),c(names_rnk[-8]))]                 #Select table with samples and rank8 only
  row.names(asv_tax_spl_rnk8_before) <- asv_tax_spl_rnk8_before$Rank8
  asv_tax_spl_rnk8_before$Rank8 <- NULL
  asv_tax_spl_rnk8_before <- as.data.frame(t(asv_tax_spl_rnk8_before))
  asv_tax_spl_rnk8_before[asv_tax_spl_rnk8_before<10] <- 0
  asv_tax_spl_rnk8_before[asv_tax_spl_rnk8_before>=10] <- 1
  #Only take taxa with at least 5 sample presence  and at least 5 sample absent 
  for (j in 1:length(min_sample_size)){
    to_keep <- which(colSums(asv_tax_spl_rnk8_before) < length(sample_names(project_data_rarefied))- min_sample_size[j] & colSums(asv_tax_spl_rnk8_before) > min_sample_size[j])
    asv_tax_spl_rnk8 <- asv_tax_spl_rnk8_before[,to_keep]
}


  
t <- t(asv_tax_spl_rnk8)
cooccur.finches <- cooccur(mat=t,
                           type="spp_site",
                           thresh=TRUE,
                           spp_names=TRUE)

summary(cooccur.finches) 
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/co-occur",)
plot(cooccur.finches)

summary_test <- cooccur.finches$results
summary_test$

  
d41_melt.df <- melt(summary_test, id.vars=c("sp1","sp2","sp1_inc","sp2_inc","obs_cooccur","prob_cooccur","exp_cooccur","sp1_name","sp2_name"),
                      measure.vars = c("p_lt","p_gt"))
  
cooccur.finches$co_occurrences

pair(mod = cooccur.finches, "Parentodinium_sp._ASV76")
pair.attributes(cooccur.finches)
pair.profile(cooccur.finches)


cooccur(mat = t, type = "spp_site", thresh = FALSE,
        spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE,
        eff_matrix = TRUE)

obs.v.exp(cooccur.finches)

library("gridExtra")
data("rodents")
data("beetles")
cooc_mod <- lapply(list(beetles, rodents),
                   FUN = function(x) cooccur(mat = x, thresh = FALSE))
cooccur.beetles <- cooc_mod[[1]]
cooccur.rodents <- cooc_mod[[2]]
grid.arrange(pair.profile(cooccur.beetles), obs_v_exp(cooccur.beetles),
             ncol = 2,main = textGrob("Beetles", gp = gpar(cex = 2), just = "top",
                                      vjust = 0.75))
grid.arrange(pair.profile(cooccur.rodents), obs_v_exp(cooccur.rodents),
                ncol = 2,main = textGrob("Rodents", gp = gpar(cex = 2), just = "top",
                                         vjust = 0.75))

########################################################

### Correlation (Abundance)
########################################################

#Set working directory
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/correlation")
sort(sample_sums(project_data))
#Create list for final table with
permute_list <- list()
permutation_vector <- rep(14418,20)
reads_threshold <- c(2,5,10,25,50)
min_sample_size <- c(5)
for(perm in 1:length(permutation_vector)){
  project_data_rarefied <- rarefy_even_depth(project_data, sample.size = permutation_vector[perm]);project_data
  #---- Create ASV table for ----#
  asv_spl.df <- as.data.frame(otu_table(project_data_rarefied) )#create ASV table
  asv_tax.df <- as.data.frame(tax_table(project_data_rarefied))
  rownames(asv_spl.df)==rownames(asv_tax.df)           #Sanity check to see if ASV from both table are similar and same order
  asv_tax_spl.df <- cbind(asv_tax.df,asv_spl.df)       #Bind both table to have ASV and TAX table together
  #---- Build prevalence table at Rank8 ----#
  for (nbreads in 1:length(reads_threshold)){
    asv_tax_spl_rnk8_before <- asv_tax_spl.df[,setdiff(colnames(asv_tax_spl.df),c(names_rnk[-8]))]                 #Select table with samples and rank8 only
    row.names(asv_tax_spl_rnk8_before) <- asv_tax_spl_rnk8_before$Rank8
    asv_tax_spl_rnk8_before$Rank8 <- NULL
    asv_tax_spl_rnk8_before <- as.data.frame(t(asv_tax_spl_rnk8_before))
    asv_tax_spl_rnk8_before[asv_tax_spl_rnk8_before<reads_threshold[nbreads]] <- 0
    #Only take taxa with at least 5 sample presence  and at least 5 sample absent 
    asv_tax_spl_rnk8_before_bin <- asv_tax_spl_rnk8_before
    asv_tax_spl_rnk8_before_bin[asv_tax_spl_rnk8_before_bin>0] <- 1
    for (j in 1:length(min_sample_size)){
      to_keep <- which(colSums(asv_tax_spl_rnk8_before_bin)>min_sample_size[j])
      asv_tax_spl_rnk8 <- asv_tax_spl_rnk8_before[,to_keep]
      #Run fisher test for all pairwise comparison
      stat_list <- list()
      for (taxa1 in colnames(asv_tax_spl_rnk8)){
        for (taxa2 in colnames(asv_tax_spl_rnk8)) {
          table_Spearman <- asv_tax_spl_rnk8[,c(taxa1,taxa2)]
          table_Spearman$double_zero <- table_Spearman[,1] + table_Spearman[,2]
          table_Spearman <- subset(table_Spearman, double_zero >0)
          is_Spearman_test <-  
            tryCatch(expr =  { 
              Spearman_test <- cor.test(table_Spearman[,1],table_Spearman[,2])
              TRUE},
              error=function(e){
                return(FALSE)}
            )
          if(is_Spearman_test)
            c(summary_table <- data.frame(Reads_Threshold = reads_threshold[nbreads],
                                          Iterate =perm,
                                          Min_Sample_Size = min_sample_size[j], 
                                          Taxa1=taxa1,
                                          Sample_Size= dim(table_Spearman)[1],
                                          Taxa2=taxa2,
                                          #Spl_taxa2= length(which(asv_tax_spl_rnk8[,taxa2]=="Presence")),
                                          Pvalue=Spearman_test$p.value,
                                          Cor =round(Spearman_test$estimate,3)
            )
            )
          if(!is_Spearman_test)
            c(summary_table <- data.frame(Reads_Threshold = reads_threshold[nbreads],
                                          Iterate =perm,
                                          Min_Sample_Size = min_sample_size[j],
                                          Taxa1=taxa1,
                                          #Spl_taxa1= length(which(asv_tax_spl_rnk8[,taxa1]=="Presence")),
                                          Taxa2=taxa2,
                                          #Spl_taxa2= length(which(asv_tax_spl_rnk8[,taxa2]=="Presence")),
                                          Pvalue=NA,
                                          Cor =NA
            )
            )
          stat_list[[length(stat_list)+1]] <- summary_table
        }
      }
      #Generate table that summarize all alpha and beta statistic
      Fisher_test_summary <- data.frame()
      for (i in 1:length(stat_list)){
        Fisher_test_summary <- rbind(Fisher_test_summary,stat_list[[i]])
      }
      Fisher_test_summary$Pvalue_Adjust <- round(p.adjust(Fisher_test_summary$Pvalue, method = "BH", n = length(Fisher_test_summary$Pvalue)),3)
      Fisher_test_summary$duplicate <- Fisher_test_summary$Taxa1==Fisher_test_summary$Taxa2
      Fisher_test_summary <- subset(Fisher_test_summary, duplicate == FALSE)
      Fisher_test_summary$Odds_sup1 <- ifelse(Fisher_test_summary$Cor>0, ">0 (Positive association)","<0 (Negative association)")
      Fisher_test_summary$Significance <- ifelse(Fisher_test_summary$Pvalue_Adjust < 0.05,"Significant","Non-significant")
      row.names(Fisher_test_summary) <- NULL
      temp <- Fisher_test_summary[,c("Taxa1","Taxa2")]
      newDf <- data.frame(t(apply(temp,1,sort)))
      Fisher_test_summary <- Fisher_test_summary[!duplicated(newDf),]
      if(dim(Fisher_test_summary)[1]>0){
        permute_list[[length(permute_list)+1]] <- Fisher_test_summary}
      #Plot histogram
      hist_plot <- ggplot(Fisher_test_summary, aes(x=Cor)) +
        geom_histogram(binwidth = 0.2)+
        geom_vline(xintercept=1, linetype='dashed') +
        theme_light() +
        theme(axis.title.y = element_text(size=15),
              axis.text.y = element_text(size=15),
              axis.title.x = element_text(size=15),
              axis.text.x = element_text(size=15))
      #Plot Bar plot number of odd ratio below and above 1
      bar_plot <- ggplot(Fisher_test_summary, aes(x=Odds_sup1)) +
        geom_bar() +
        theme_light() +
        theme(axis.title.y = element_text(size=15),
              axis.text.y = element_text(size=15),
              axis.title.x = element_text(size=15),
              axis.text.x = element_text(size=9))
      #Scatter plot, Odds ratio and Pvalue
      bubble_plot <- ggplot(Fisher_test_summary,aes(x=Cor, y=Pvalue_Adjust, color=Significance)) +
        geom_point() +
        xlim(0,20) +
        ylim(1,0.0001) +
        theme_light() +
        geom_vline(xintercept=1, linetype='dashed') +
        geom_hline(yintercept=0.05, linetype='dashed') +
        theme(axis.title.y = element_text(size=15),
              axis.text.y = element_text(size=15),
              axis.title.x = element_text(size=15),
              axis.text.x = element_text(size=15))
      combined <- ggarrange(
        bubble_plot,                # First row with line plot
        widths = c(1, 0.2),
        # Second row with box and histogram
        ggarrange(bar_plot, hist_plot, ncol = 2, labels = c("B", "C"), widths = c(0.7, 1)), 
        nrow = 2, 
        labels = "A"       # Label of the line plot
      )
      #ggsave(combined,filename = paste0(perm,"_",reads_threshold[i],"_","Spearman_test_multipannel.pdf"), width = 10, height = 10)
    }
  }
}


#Convert list into table
Fisher_permute_summary <- data.frame()
for (i in 1:length(permute_list)){
  Fisher_permute_summary <- rbind(Fisher_permute_summary,permute_list[[i]])
}


#Add additional information to table
Fisher_permute_summary$nb_significant <- 1 #To later calculate number of time pairwise was significant
Fisher_permute_summary$pair <- paste0(Fisher_permute_summary$Taxa1,"_",Fisher_permute_summary$Taxa2)
Fisher_permute_summary$pairwise <- ifelse(Fisher_permute_summary$Significance=="Significant",paste0(Fisher_permute_summary$Taxa1,"_",Fisher_permute_summary$Taxa2),"NS")
to_keep <- unique(Fisher_permute_summary$pairwise)
to_keep <- to_keep[!to_keep %in% "NS"]
Fisher_permute_summary_sg <- subset(Fisher_permute_summary, pair %in% to_keep)
Fisher_permute_summary_sg$pairwise <- Fisher_permute_summary_sg$pair
to_keep2 <- setdiff(Fisher_permute_summary$pair, to_keep)
Fisher_permute_summary_nsg <- subset(Fisher_permute_summary, pair %in% to_keep2)
Fisher_permute_summary <- rbind(Fisher_permute_summary_sg,Fisher_permute_summary_nsg)

#Plot histogram
hist_plot <- ggplot(Fisher_permute_summary, aes(x=Cor)) +
  geom_histogram(binwidth = 0.2)+
  geom_vline(xintercept=0, linetype='dashed') +
  theme_light() +
  facet_grid(rows = vars(Reads_Threshold), scales = "free_y", switch = "y", space = "free_y")+ 
  theme(axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=15))
plot(hist_plot)

#Plot Bar plot number of odd ratio below and above 1
bar_plot <- ggplot(Fisher_permute_summary, aes(x=Odds_sup1)) +
  geom_bar() +
  theme_light() +
  theme(axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=9))
print(bar_plot)

#Bubble plot for proportion of time test is significant
Fisher_permute_summary_sig <- subset(Fisher_permute_summary, Pvalue_Adjust <0.05) # Select only signiifcant adjusted pvalue
Fisher_permute_summary_groupby <- Fisher_permute_summary_sig %>% 
  group_by(Taxa1,Taxa2) %>%
  summarise(avg_adpvl = mean(Pvalue_Adjust),    #Calculate mean of pvalue 
            Mean_OddRatio = mean(Cor),   #Calculate mean of odds ratio
            Permutation = sum(nb_significant))  #Calculate total number of time the pairwise is significant
Fisher_permute_summary_groupby$Permutation <- Fisher_permute_summary_groupby$Permutation/perm*100
setorder(Fisher_permute_summary_groupby, cols= - "Permutation") #Reoerder table by the number of time it was significant
order_legend_pairwise <- paste0(Fisher_permute_summary_groupby$Taxa1,"_",Fisher_permute_summary_groupby$Taxa2) #Create vector to reorder scatter plot color legend
Fisher_permute_summary_groupby$Taxa1 <- gsub("_"," ",Fisher_permute_summary_groupby$Taxa1) #Remove underscore for aesthetic
Fisher_permute_summary_groupby$Taxa2 <- gsub("_"," ",Fisher_permute_summary_groupby$Taxa2) #Remove underscore for aesthetic
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/correlation/bubble_Spearman_test.pdf", height =  5, width = 10)
plot <- ggplot(Fisher_permute_summary_groupby, aes(x=Taxa1, y=Taxa2, color=Mean_OddRatio, size=Permutation)) +
  geom_point() +
  scale_size(range = c(0, 15)) +
  theme_light() +
  scale_fill_gradient(low = "royalblue4", high = "red") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=15, hjust = 1, angle = 20, ))
print(plot)
dev.off()

#Scatter plot, Odds ratio and Pvalue
colours <- c("grey75","dodgerblue","coral1","chartreuse4","mAgenta","grey14","seagreen1","darkorange","purple","forestgreen","firebrick","aquamarine","lightpink","grey14","seagreen1","darkorange","purple","forestgreen","firebrick","aquamarine","lightpink","purple","forestgreen","firebrick","aquamarine","lightpink","grey14","seagreen1","darkorange","purple","forestgreen","firebrick","aquamarine","lightpink","dodgerblue",
             "mAgenta","grey14","seagreen1","darkorange","purple","forestgreen","firebrick","aquamarine","lightpink","grey14","seagreen1","darkorange","purple","forestgreen","firebrick","aquamarine","lightpink","purple","forestgreen","firebrick","aquamarine","lightpink","grey14","seagreen1","darkorange","purple","forestgreen","firebrick","aquamarine","lightpink","dodgerblue")
Fisher_permute_summary$pairwise <- factor(Fisher_permute_summary$pairwise, level=c("NS",order_legend_pairwise))
Fisher_permute_summary$Reads_Threshold <- paste0(Fisher_permute_summary$Reads_Threshold," reads")
Fisher_permute_summary$Reads_Threshold <- factor(Fisher_permute_summary$Reads_Threshold, level = c("2 reads","5 reads","10 reads","25 reads","50 reads"))
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/correlation/Scatter_Spearman_20perm_MultiThreshold_fecal.pdf", height =  40, width = 60)
scatter_plot <- ggplot(Fisher_permute_summary,aes(x=Cor, y=Pvalue_Adjust, color=pairwise, shape=pairwise)) +
  geom_point(size=5) +
  scale_shape_manual(values=c(1:38))+
  xlim(-1,1) +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values=colours) +
  theme_light() +
  geom_vline(xintercept=0, linetype='dashed', colour = "red") +
  geom_text(aes(x=-0.1, label="\nCor = 1", y=0.035), colour="red", angle=0, size=17) +
  geom_hline(yintercept=0.05, linetype='dashed', colour = "orange") +
  geom_text(aes(y=0.15, label="\nAdj.Pvalue = 0.05", x=0.8), colour="orange", angle=0, size=17) +
  facet_grid(rows = vars(Reads_Threshold), scales = "free_y", switch = "y", space = "free_y")+ 
  labs(x="Coefficient correlation", y="Adjusted p-value") +
  theme(axis.title.y = element_text(size=40),
        axis.text.y = element_text(size=40),
        axis.title.x = element_text(size=40),
        axis.text.x = element_text(size=40),
        legend.text = element_text(size=28),
        strip.text.y = element_text(size=60),
        legend.position = 'bottom') +
  guides(col=guide_legend(ncol=2))
print(scatter_plot)
dev.off()

#Bubble plot significant pairwise across Reads_Threshold and Iterations
Fisher_permute_summary_remove_sign <- subset(Fisher_permute_summary, Significance == "Significant")
Fisher_permute_summary_remove_sign$Min_Sample_Size <- as.factor(Fisher_permute_summary_remove_sign$Min_Sample_Size)
Fisher_permute_summary_remove_sign$Min_Sample_Size_Reads_Threshold <- paste0(Fisher_permute_summary_remove_sign$Min_Sample_Size,"_",Fisher_permute_summary_remove_sign$Reads_Threshol)
Fisher_permute_summary_remove_sign$Min_Sample_Size_Reads_Threshold <- factor(Fisher_permute_summary_remove_sign$Min_Sample_Size_Reads_Threshold, level = c("5_1","5_2","5_5","5_10","5_25","5_50","5_100","5_200","10_500","5_500","10_1","10_2","10_5","10_10","10_25","10_50","10_100","10_200"))
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/correlation/bubble_Spearman_20perm_MultiThreshold_5-10spl_fecal.pdf", height =  15, width = 35)
Bubble_plot <- ggplot(Fisher_permute_summary_remove_sign,aes(x=Iterate, y=pairwise, color=Odds_sup1, size=Cor)) +
  geom_point() +
  theme_light() +
  facet_grid(cols = vars(Min_Sample_Size_Reads_Threshold), scales = "free_x", switch = "x", space = "free_x") + 
  theme(axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=20),
        legend.text = element_text(size=14),
        legend.position = 'right')
plot(Bubble_plot)
dev.off()



combined <- ggarrange(
  bubble_plot,                # First row with line plot
  # Second row with box and histogram
  ggarrange(plot,bar_plot,hist_plot, widths = c(0.6,0.2,0.3), ncol=3),
  nrow = 2, 
  labels = "A",       # Label of the line plot
  heights = c(1,0.5),
  widths = c(1,1)
)
ggsave(combined,filename = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test/fisher_test_multipannel.pdf", width = 30, height = 20)

write.table(Fisher_permute_summary,"/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/bubble_plot_fisher_test.txt", sep="\t", row.names = FALSE)
write.table(Fisher_test_summary,"/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test/fisher_test.txt", sep ="\t" )
Fisher_permute_summary <- fread("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/bubble_plot_fisher_test.txt")


#Understand Infinite value
x <- subset(Fisher_permute_summary, Odds_ratio == "Inf")
x2 <- subset(x, Pvalue_Adjust<0.05 )
taxa1_taxa2_cont_table <- table(asv_tax_spl_rnk8[,"Parentodinium_sp._ASV132"],asv_tax_spl_rnk8[,"Cycloposthium_bipalmatum"])
fisher.test(taxa1_taxa2_cont_table)
dimnames(taxa1_taxa2_cont_table) <- list(Parentodinium_sp._ASV132 = c("Absence_Taxa1", "Presence_Taxa1"),
                                         Cycloposthium_bipalmatum = c("Absence_Taxa2", "Presence_Taxa2"))
table <- as.data.frame(taxa1_taxa2_cont_table)
colours <- c("orangered1","blue")
plot <- ggplot(table, aes(x=Cycloposthium_bipalmatum, y=Freq, fill=Parentodinium_sp._ASV132, color = colours )) +
  geom_bar(stat="identity", width=1, color="black", linetype = "solid") +  #height of the column equal the value, stacked
  scale_fill_manual(values=colours) +
  theme_bw()+
  labs(x ="Cycloposthium_bipalmatum" ) +
  #ggtitle(paste0("Total:",total_sample,", Sens:",sensitivity, ", Spec:",specificity,", Accu:",accuracy,", Pval(Fisher):",fisher_Pval))+
  theme(legend.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        plot.title=element_text(size=8,face="bold"))

#
x <- subset(Fisher_permute_summary, pairwise == "Blastocystis_Unknown_ASV45_Enteromonad_Unknown_ASV73")
x2 <- subset(x, Pvalue_Adjust>0.2 )

ggplot(x, aes(x=Pvalue_Adjust, y=Spl_taxa1)) +
  geom_point() +
  geom_vline(xintercept=0.05, linetype='dashed')

test <- cor.test(x$Pvalue_Adjust,x$Spl_taxa1)
summary(test)

taxa1_taxa2_cont_table <- table(asv_tax_spl_rnk8[,"Blastocystis_Unknown_ASV45"],asv_tax_spl_rnk8[,"Enteromonad_Unknown_ASV73"])
dimnames(taxa1_taxa2_cont_table) <- list(Taxa1 = c("Absence_Taxa1", "Presence_Taxa1"),
                                         Taxa2 = c("Absence_Taxa2", "Presence_Taxa2"))

colours <- c("orangered1","blue")
plot <- ggplot(taxa1_taxa2_cont_table, aes(x=Taxa1, y=Freq, fill=Taxa2, color = colours )) +
  geom_bar(stat="identity", width=1, color="black", linetype = "solid") +  #height of the column equal the value, stacked
  scale_fill_manual(values=colours) +
  theme_bw()+
  labs(y ="Taxa2" , x ="Taxa1" ) +
  #ggtitle(paste0("Total:",total_sample,", Sens:",sensitivity, ", Spec:",specificity,", Accu:",accuracy,", Pval(Fisher):",fisher_Pval))+
  theme(legend.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        plot.title=element_text(size=8,face="bold"))



########################################################


### Heatmap with Pearson Correlation                 ###
########################################################

#Set working directory
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/correlation")
sort(sample_sums(project_data))
#Create list for final table with
permute_list <- list()
permutation_vector <- rep(14418,1)
reads_threshold <- c(10)#
min_sample_size <- c(5)
for(perm in 1:length(permutation_vector)){
  project_data_rarefied <- rarefy_even_depth(project_data, sample.size = permutation_vector[perm]);project_data
  #---- Create ASV table for ----#
  asv_spl.df <- as.data.frame(otu_table(project_data_rarefied) )#create ASV table
  asv_tax.df <- as.data.frame(tax_table(project_data_rarefied))
  rownames(asv_spl.df)==rownames(asv_tax.df)           #Sanity check to see if ASV from both table are similar and same order
  asv_tax_spl.df <- cbind(asv_tax.df,asv_spl.df)       #Bind both table to have ASV and TAX table together
  #---- Build prevalence table at Rank8 ----#
  for (nbreads in 1:length(reads_threshold)){
    asv_tax_spl_rnk8_before <- asv_tax_spl.df[,setdiff(colnames(asv_tax_spl.df),c(names_rnk[-8]))]                 #Select table with samples and rank8 only
    row.names(asv_tax_spl_rnk8_before) <- asv_tax_spl_rnk8_before$Rank8
    asv_tax_spl_rnk8_before$Rank8 <- NULL
    asv_tax_spl_rnk8_before <- as.data.frame(t(asv_tax_spl_rnk8_before))
    asv_tax_spl_rnk8_before[asv_tax_spl_rnk8_before<reads_threshold[nbreads]] <- 0
    #Only take taxa with at least 5 sample presence  and at least 5 sample absent 
    asv_tax_spl_rnk8_before_bin <- asv_tax_spl_rnk8_before
    asv_tax_spl_rnk8_before_bin[asv_tax_spl_rnk8_before_bin>0] <- 1
    for (j in 1:length(min_sample_size)){
      to_keep <- which(colSums(asv_tax_spl_rnk8_before_bin)>min_sample_size[j])
      asv_tax_spl_rnk8 <- asv_tax_spl_rnk8_before[,to_keep]
      #Run fisher test for all pairwise comparison
      stat_list <- list()
      for (taxa1 in colnames(asv_tax_spl_rnk8)){
        for (taxa2 in colnames(asv_tax_spl_rnk8)) {
          table_Spearman <- asv_tax_spl_rnk8[,c(taxa1,taxa2)]
          table_Spearman$double_zero <- table_Spearman[,1] + table_Spearman[,2]
          table_Spearman <- subset(table_Spearman, double_zero >0)
          is_Spearman_test <-  
            tryCatch(expr =  { 
              Spearman_test <- cor.test(table_Spearman[,1],table_Spearman[,2])
              TRUE},
              error=function(e){
                return(FALSE)}
            )
          if(is_Spearman_test)
            c(summary_table <- data.frame(Reads_Threshold = reads_threshold[nbreads],
                                          Iterate =perm,
                                          Min_Sample_Size = min_sample_size[j], 
                                          Taxa1=taxa1,
                                          Sample_Size= dim(table_Spearman)[1],
                                          Taxa2=taxa2,
                                          #Spl_taxa2= length(which(asv_tax_spl_rnk8[,taxa2]=="Presence")),
                                          Pvalue=Spearman_test$p.value,
                                          Cor =round(Spearman_test$estimate,3)
            )
            )
          if(!is_Spearman_test)
            c(summary_table <- data.frame(Reads_Threshold = reads_threshold[nbreads],
                                          Iterate =perm,
                                          Min_Sample_Size = min_sample_size[j],
                                          Taxa1=taxa1,
                                          #Spl_taxa1= length(which(asv_tax_spl_rnk8[,taxa1]=="Presence")),
                                          Taxa2=taxa2,
                                          #Spl_taxa2= length(which(asv_tax_spl_rnk8[,taxa2]=="Presence")),
                                          Pvalue=NA,
                                          Cor =NA
            )
            )
          stat_list[[length(stat_list)+1]] <- summary_table
        }
      }
      #Generate table that summarize all alpha and beta statistic
      Fisher_test_summary <- data.frame()
      for (i in 1:length(stat_list)){
        Fisher_test_summary <- rbind(Fisher_test_summary,stat_list[[i]])
      }
      Fisher_test_summary$Pvalue_Adjust <- round(p.adjust(Fisher_test_summary$Pvalue, method = "BH", n = length(Fisher_test_summary$Pvalue)),3)
      Fisher_test_summary$Odds_sup1 <- ifelse(Fisher_test_summary$Cor>0, ">0 (Positive association)","<0 (Negative association)")
      Fisher_test_summary$Significance <- ifelse(Fisher_test_summary$Pvalue_Adjust < 0.05,"Significant","Non-significant")
      row.names(Fisher_test_summary) <- NULL
      temp <- Fisher_test_summary[,c("Taxa1","Taxa2")]
      newDf <- data.frame(t(apply(temp,1,sort)))
      Fisher_test_summary <- Fisher_test_summary[!duplicated(newDf),]
      if(dim(Fisher_test_summary)[1]>0){
        permute_list[[length(permute_list)+1]] <- Fisher_test_summary}
    }
  }
}



#Heatmapt at 10 reads
Fisher_permute_summary_heatmap <- Fisher_test_summary
Fisher_permute_summary_heatmap$r_if_sig <- ifelse(Fisher_permute_summary_heatmap$Pvalue_Adjust<0.05,Fisher_permute_summary_heatmap$Cor,NA)
Fisher_permute_summary_heatmap <- Fisher_permute_summary_heatmap %>% arrange(Taxa2)
Fisher_permute_summary_heatmap$Taxa1 <- gsub("_"," ",Fisher_permute_summary_heatmap$Taxa1)
Fisher_permute_summary_heatmap$Taxa2 <- gsub("_"," ",Fisher_permute_summary_heatmap$Taxa2)
Fisher_permute_summary_heatmap$Cor <- ifelse(Fisher_permute_summary_heatmap$Taxa1==Fisher_permute_summary_heatmap$Taxa2,NA,Fisher_permute_summary_heatmap$Cor)
Fisher_permute_summary_heatmap$r_if_sig <- ifelse(Fisher_permute_summary_heatmap$Taxa1==Fisher_permute_summary_heatmap$Taxa2,NA,Fisher_permute_summary_heatmap$r_if_sig)


#---- Create heatmap plot with pearson signifiancnt value
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/correlation/heatmap_correlation.pdf", height =  5, width = 8)
plot <- Fisher_permute_summary_heatmap %>% 
  ggplot(aes(Taxa1, Taxa2, fill=Cor, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, 
       fill = "Pearson's\nCorrelation") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#A63446",high="#0C6291", limits=c(-1,1)) +
  geom_text(size=2.9) +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10))
print(plot)
dev.off()



########################################################


### Association between Presence and alpha diversity
########################################################

########################################################

### Association between Presence and beta diversity
########################################################

########################################################

################################################################################################################
###                  RELATIONSHIP BETWEEN INTESINAL EUKARYOTES and ENVIRONMENT AND HOST FACTORS              ###
################################################################################################################

########################################################
###                Beta diversity                    ###
########################################################

sort(sample_sums(project_data))
project_data_rarefied <- rarefy_even_depth(project_data, sample.size = 14418);project_data_rarefied

### Sites
########################################################

#All Site
sample_df <- data.frame(sample_data(project_data_rarefied))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Site, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Site, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x = NMDS.jaccard1, y = NMDS.jaccard2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Moba
project_data_rarefied_no_moba <- project_data_rarefied %>%
  subset_samples(Site != "Moba")

sample_df <- data.frame(sample_data(project_data_rarefied_no_moba))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_moba, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Site, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_no_moba, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_moba))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_no_moba))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_moba, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Site, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied_no_moba, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_moba))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Maya
project_data_rarefied_no_maya <- project_data_rarefied %>%
  subset_samples(Site != "Maya")

sample_df <- data.frame(sample_data(project_data_rarefied_no_maya))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_maya, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Site, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_no_maya, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_maya))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_no_maya))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_maya, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Site, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied_no_maya, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_maya))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Romani
project_data_rarefied_no_romani <- project_data_rarefied %>%
  subset_samples(Site != "Romani")

sample_df <- data.frame(sample_data(project_data_rarefied_no_romani))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_romani, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Site, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_no_romani, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_romani))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_no_romani))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_romani, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Site, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied_no_romani, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_romani))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Lokoue
project_data_rarefied_no_lokoue <- project_data_rarefied %>%
  subset_samples(Site != "Lokoue")

sample_df <- data.frame(sample_data(project_data_rarefied_no_lokoue))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_lokoue, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Site, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_no_lokoue, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_lokoue))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_no_lokoue))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_lokoue, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Site, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied_no_lokoue, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_lokoue))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()


########################################################

### Season
########################################################

#All Season
sample_df <- data.frame(sample_data(project_data_rarefied))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Season, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x = NMDS.jaccard1, y = NMDS.jaccard2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Long-Dry
project_data_rarefied_no_longdry <- project_data_rarefied %>%
  subset_samples(Season != "Long-Dry")

sample_df <- data.frame(sample_data(project_data_rarefied_no_longdry))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_longdry, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_no_longdry, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_longdry))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_no_longdry))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_longdry, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Season, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied_no_longdry, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_longdry))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Short dry
project_data_rarefied_no_maya <- project_data_rarefied %>%
  subset_samples(Season != "Short-Dry")

sample_df <- data.frame(sample_data(project_data_rarefied_no_maya))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_maya, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_no_maya, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_maya))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_no_maya))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_maya, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Season, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied_no_maya, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_maya))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Long rain
project_data_rarefied_no_romani <- project_data_rarefied %>%
  subset_samples(Season != "Long-Rain")

sample_df <- data.frame(sample_data(project_data_rarefied_no_romani))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_romani, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_no_romani, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_romani))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_no_romani))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_romani, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Season, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied_no_romani, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_romani))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Short rain
project_data_rarefied_no_lokoue <- project_data_rarefied %>%
  subset_samples(Season != "Short-Rain")

sample_df <- data.frame(sample_data(project_data_rarefied_no_lokoue))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_lokoue, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_no_lokoue, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_lokoue))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_no_lokoue))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_lokoue, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Season, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied_no_lokoue, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_lokoue))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Long Rain and Short dry
project_data_rarefied_no_moba <- project_data_rarefied %>%
  subset_samples(Season != "Long-Rain")%>%
  subset_samples(Season != "Short-Dry")

sample_df <- data.frame(sample_data(project_data_rarefied_no_moba))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_moba, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_no_moba, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_moba))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  geom_text(data=NMDS, aes(label=Site), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_no_moba))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_moba, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Season, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied_no_moba, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_moba))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Moba
project_data_rarefied_no_moba <- project_data_rarefied %>%
  subset_samples(Site != "Moba")

sample_df <- data.frame(sample_data(project_data_rarefied_no_moba))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_moba, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_no_moba, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_moba))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  geom_text(data=NMDS, aes(label=Site), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_no_moba))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_rarefied_no_moba, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Season, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_rarefied_no_moba, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_no_moba))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()



# No Moba no Maya
project_data_Season <- project_data_rarefied %>%
  subset_samples(Site != "Moba") %>%
  subset_samples(Site != "Maya")

sample_df <- data.frame(sample_data(project_data_Season))
project_bray.rarefied.bray <- phyloseq::distance(project_data_Season, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~Season*Site*Age*Yaws*Social_Status, data=sample_df, method="bray")
res.adonis.rarefied.bray

#Season
sample_df <- data.frame(sample_data(project_data_Season))
project_bray.rarefied.bray <- phyloseq::distance(project_data_Season, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_Season, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_Season))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Season_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Season
sample_df <- data.frame(sample_data(project_data_Season))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_Season, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Season, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.jaccard <- ordinate(physeq = project_data_Season, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_Season))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Season_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()


########################################################

### Maturity
########################################################

#Remove data from Moba and Maya Sites (no information on dates)
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Moba") %>%
  subset_samples(Site !="Maya")

#Age
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Age, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Age, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Age_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Age)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","blue2","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Age
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Age, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Age, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Age_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Age)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  scale_color_manual(values=c("darkgreen","blue2","orangered3")) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()


#Remove Maya, Moba and Romani
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Romani")%>%
  subset_samples(Site !="Maya")%>%
  subset_samples(Site !="Moba")

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Age, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Age)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  geom_text(data=NMDS, aes(label=Gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Age, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Age)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  geom_text(data=NMDS, aes(label=Gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Maya, Moba and Lokoue
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Lokoue")%>%
  subset_samples(Site !="Maya")%>%
  subset_samples(Site !="Moba")

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Age, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Age)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  geom_text(data=NMDS, aes(label=Gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","blue2","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Age, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Age)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  geom_text(data=NMDS, aes(label=Gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","blue2","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()


########################################################

### Social Status
########################################################

#Remove data from Moba and Maya Sites (no information on dates)
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Moba") %>%
  subset_samples(Site !="Maya")

#Social_Status
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Social_Status, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Social_Status, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Social_Status_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Social_Status)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Social_Status
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Social_Status, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Social_Status, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Social_Status_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Social_Status)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Maya, Moba and Romani
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Romani")%>%
  subset_samples(Site !="Maya")%>%
  subset_samples(Site !="Moba")

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Social_Status, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Social_Status)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  geom_text(data=NMDS, aes(label=Gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Yaws
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Social_Status, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Social_Status)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  geom_text(data=NMDS, aes(label=Gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Maya, Moba and Lokoue
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Lokoue")%>%
  subset_samples(Site !="Maya")%>%
  subset_samples(Site !="Moba")

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Social_Status, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Social_Status)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  geom_text(data=NMDS, aes(label=Gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Social_Status, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Social_Status)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  geom_text(data=NMDS, aes(label=Gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

########################################################

### Yaws
########################################################

#Remove data from Moba and Maya Sites (no information on dates)
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Moba") %>%
  subset_samples(Site !="Maya")

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Yaws, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Yaws)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Yaws
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Yaws, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Yaws)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Maya, Moba and Romani
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Romani")%>%
  subset_samples(Site !="Maya")%>%
  subset_samples(Site !="Moba")

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Yaws, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Yaws)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Yaws, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Yaws)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove Maya, Moba and Lokoue
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Lokoue")%>%
  subset_samples(Site !="Maya")%>%
  subset_samples(Site !="Moba")

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Yaws, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS$Gorilla_name_Date <- paste0(NMDS$Gorilla_name,"_",NMDS$Date)
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_with_gorilla_features/NMDS_bray_yaws_Romani.pdf", height =  5, width = 7)
NMDS_Yaws_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Yaws)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  geom_text(data=NMDS, aes(label=Gorilla_name_Date), size=1.5, vjust = 3) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic() +
  labs(y="NMDS AXE2", x="NMDS AXE1") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

print(NMDS_Yaws_bray)
dev.off()


sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Yaws, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Yaws)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  geom_text(data=NMDS, aes(label=Gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

########################################################

### Yaws, Maturity, Social Status, Site, Season 
########################################################

#Remove data from Moba and Maya Sites (no information on dates)
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Moba") %>%
  subset_samples(Site !="Maya")

sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Yaws*Age*Site*Social_Status*Season, data=sample_df, method="bray")
res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Age*Site*Social_Status*Season*Yaws, data=sample_df, method="bray")
res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Site*Social_Status*Season*Yaws*Age, data=sample_df, method="bray")
res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Social_Status*Season*Yaws*Age*Site, data=sample_df, method="bray")
res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season*Yaws*Age*Site*Social_Status, data=sample_df, method="bray")
res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season*Social_Status, data=sample_df, method="bray")
res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Social_Status*Season, data=sample_df, method="bray")
res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season*Social_Status*Site, data=sample_df, method="bray")
res.adonis.rarefied.bray

########################################################

### Gorilla name and temporal data
########################################################

#Remove data from Moba and Maya Sites (no information on dates)
project_data_exp <- project_data %>% subset_samples(temporal == "yes") 

#gorilla name
sample_df <- data.frame(sample_data(project_data_exp))
project_bray.rarefied.bray <- phyloseq::distance(project_data_exp, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ gorilla_name, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_exp, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_exp))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$gorilla_name, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
## Create Colour Palettes ###
n <- length(unique(NMDS$gorilla_name))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-16]
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_gorilla_name_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = gorilla_name)) + 
  geom_point(size=4) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=col_vector) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#gorilla_name
sample_df <- data.frame(sample_data(project_data_exp))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_exp, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ gorilla_name, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
NMDS.jaccard <- ordinate(physeq = project_data_exp, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_exp))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$gorilla_name, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_gorilla_name_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = gorilla_name)) + 
  geom_point(size=4) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=col_vector) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

########################################################

### Loop of all analysis
########################################################
#Season
project_data_rarefied <- project_data_rarefied %>%
  subset_samples(Site != "Moba")

#Site and Season
sample_df <- data.frame(sample_data(project_data_rarefied))
project_bray.rarefied.jaccard <- phyloseq::distance(project_data_rarefied, method = "jaccard")

index_beta_vector <- c("bray","jaccard")
factors_vector <- c("Site","Season","Diet","Season_combined")
NMDS_all_list <- list()
PERMANOVA_list <- list()
for (i in 1:length(factors_vector)){
  table_PERMANOVA_bray_jaccard <- data.frame(Factor = factors_vector[i])
  for (j in 1:length(index_beta_vector)){
    phyloseq_distance_index_beta <- phyloseq::distance(project_data_rarefied, method = index_beta_vector[j])
    form <- as.formula(paste("phyloseq_distance_index_beta", factors_vector[i], sep="~"))
    res.adonis.rarefied.bray <- adonis(form, data=sample_df, method=index_beta_vector[j])
    r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
    pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
    table_PERMANOVA <- data.frame(Pvalue = pval_treatment,
                                  Rsquared = r2_treatment)
    colnames(table_PERMANOVA) <- paste0(colnames(table_PERMANOVA),"_",index_beta_vector[j])
    table_PERMANOVA_bray_jaccard <- cbind(table_PERMANOVA_bray_jaccard,table_PERMANOVA)
    }
  PERMANOVA_list[[length(PERMANOVA_list)+1]] <- table_PERMANOVA_bray_jaccard
}
#Stat table
table_test <- data.frame()
for (i in 1:length(PERMANOVA_list)){
  table_test <- rbind(table_test,PERMANOVA_list[[i]])
}


    
    NMDS.bray <- ordinate(physeq = project_data_rarefied, method = "NMDS", distance = "bray")
    NMDS <- as.data.frame(sample_data(project_data_rarefied))
    bray <- as.data.frame(NMDS.bray$points)
    row.names(bray) == row.names(NMDS) #sanity check #tests as true
    NMDS$NMDS.bray1 <- bray$MDS1
    NMDS$NMDS.bray2 <- bray$MDS2
    ## Create Colour Palettes ###
    n <- length(unique(NMDS$Gorilla_name))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-16]
    #plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
    title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
    NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = get(factors_vector[i]))) + 
      geom_point(size=4) +
      stat_ellipse(level = 0.95, linetype = 2) +
      ggtitle(as.character(title)) +
      #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
      scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      theme_classic()
    # Add plot to list
  NMDS_all_list[[length(NMDS_all_list)+1]] <- NMDS_Site_bray}
  ##Run Jaccard analysis
  form <- as.formula(paste("project_bray.rarefied.jaccard", i, sep="~"))
  res.adonis.rarefied.bray <- adonis(form, data=sample_df, method="jaccard")
  r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
  pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
  title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
  NMDS.bray <- ordinate(physeq = project_data_rarefied, method = "NMDS", distance = "jaccard")
  NMDS <- as.data.frame(sample_data(project_data_rarefied))
  bray <- as.data.frame(NMDS.bray$points)
  row.names(bray) == row.names(NMDS) #sanity check #tests as true
  NMDS$NMDS.bray1 <- bray$MDS1
  NMDS$NMDS.bray2 <- bray$MDS2
  #plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
  NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = get(i))) + 
    geom_point(size=4) +
    stat_ellipse(level = 0.95, linetype = 2) +
    ggtitle(as.character(title)) +
    #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
    scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme_classic()
  # Add plot to list
  NMDS_list[[length(NMDS_list)+1]] <- NMDS_Site_jaccard}

  


########################################################

#control <- ggarrange(NMDS_Season_jaccard,NMDS_Season_bray, labels = c("A","B"))
#control2 <- ggarrange(control,alpha_Site,nrow = 2, labels = c("","C"))










### combine plot
########################################################
NMDS_fecal_samples <- ggarrange(NMDS_Site_bray,NMDS_Site_jaccard,NMDS_Season_bray,NMDS_Season_jaccard,NMDS_Age_bray,NMDS_Age_jaccard,NMDS_Social_Status_bray,NMDS_Social_Status_jaccard,NMDS_Yaws_bray,NMDS_Yaws_jaccard,NMDS_gorilla_name_bray,NMDS_gorilla_name_jaccard,
                                labels = c("A","B","C","D","E","F","G","H","I","J","K","L"), nrow = 6, ncol = 2)
ggsave(NMDS_fecal_samples,filename = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_with_gorilla_features/beta/interaction_with_features_NMDS_fecal.pdf", width = 15, height = 25)



NMDS_individuals_fecal_samples <- ggarrange(NMDS_fecal_samples,NMDS_individuals,
                                            nrow = 1, ncol = 2)
ggsave(NMDS_individuals_fecal_samples,filename = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_with_gorilla_features/beta/interaction_with_features_NMDS_all.pdf", width = 20, height = 30)

########################################################

########################################################
###                Alpha diversity                   ###
########################################################

#Site
########################################################
my_comparisons <- list(c("Romani", "Lokoue"),c("Romani", "Maya"),c("Romani", "Moba"),
                       c("Lokoue", "Maya"),c("Lokoue", "Moba"),
                       c("Maya", "Moba"))

sort(sample_sums(project_data))
project_data_rarefied <- rarefy_even_depth(project_data, sample.size = 14418);project_data_rarefied


alpha_Site <- plot_richness(project_data_rarefied, x = "Site", color = "Site", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot()+
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  stat_compare_means(comparisons=my_comparisons, aes(label = paste0("p =", ..p.format..)))
########################################################

#Season
########################################################
project_data_rarefied_Season <- project_data_rarefied %>%
  subset_samples(Season != "NA")

my_comparisons <- list(c("Short-Dry", "Short-Rain"),c("Short-Dry", "Long-Dry"),c("Short-Dry", "Long-Rain"),
                       c("Short-Rain", "Long-Rain"),c("Short-Rain", "Long-Dry"),
                       c("Long-Rain", "Long-Dry"))

alpha_Season <- plot_richness(project_data_rarefied_Season, x = "Season", color = "Season", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot()+
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  stat_compare_means(comparisons=my_comparisons, aes(label = paste0("p =", ..p.format..)))
########################################################

#Age
########################################################
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Moba") %>%
  subset_samples(Site !="Maya")

my_comparisons <- list(c("BB", "SB"),c("BB", "sub"),c("SB", "sub"))

alpha_Age <- plot_richness(project_data_rarefied_host_factor, x = "Age", color = "Age", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot()+
  scale_color_manual(values=c("darkgreen","blue2","orangered3")) +
  stat_compare_means(comparisons=my_comparisons, aes(label = paste0("p =", ..p.format..)))

########################################################

#Social status 
########################################################
alpha_Social_Status <- plot_richness(project_data_rarefied_host_factor, x = "Social_Status", color = "Social_Status", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot()+
  scale_color_manual(values=c("orangered3","darkgreen")) +
  stat_compare_means(aes(label = paste0("p =", ..p.format..)))
########################################################

#Yaws
########################################################

project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Moba") %>%
  subset_samples(Site !="Maya")

alpha_Yaws <- plot_richness(project_data_rarefied_host_factor, x = "Yaws", color = "Yaws", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot()+
  scale_color_manual(values=c("darkgreen","orangered3")) +
  stat_compare_means(aes(label = paste0("p =", ..p.format..)))

project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Moba") %>%
  subset_samples(Site !="Maya") %>%
  subset_samples(Site !="Lokoue")

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_with_gorilla_features/Boxplot_Observed_Richness_yaws_Romani.pdf", height =  5, width = 5)
alpha_Yaws <- plot_richness(project_data_rarefied_host_factor, x = "Yaws", color = "Yaws", measures = c("Observed", "Shannon")) +
  geom_point()+
  geom_boxplot() +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  stat_compare_means(aes(label = paste0("p =", ..p.format..))) +
  theme_light() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
print(alpha_Yaws)
dev.off()
  

project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Moba") %>%
  subset_samples(Site !="Maya") %>%
  subset_samples(Site !="Romani")

alpha_Yaws <- plot_richness(project_data_rarefied_host_factor, x = "Yaws", color = "Yaws", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot()+
  scale_color_manual(values=c("darkgreen","orangered3")) +
  stat_compare_means(aes(label = paste0("p =", ..p.format..)))

########################################################

#combine plot
########################################################
alpha <- ggarrange(alpha_Site,alpha_Age,alpha_Social_Status,alpha_Yaws,alpha_Season, labels = c("A","B","C","D","E"), nrow = 3, ncol = 2)
ggsave(alpha,filename = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_with_gorilla_features/alpha/alpha_individual_allMetrics.pdf", width = 30, height = 20)

#Yaws per Site
alpha_Yaws_Site <- plot_richness(project_data, x = "Site", color = "Yaws", measures = c("Shannon")) +
  geom_boxplot(position=position_dodge2(padding=0.5, width = 0.5)) + 
  stat_compare_means(aes(label = paste0("p =", ..p.format..))) +
  scale_color_manual(values=c("darkgreen","orangered3"))



sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
otu <- otu_table(t(project_data_rarefied_host_factor))
shannon <- diversity(otu, index = "shannon", MARGIN = 1, base = exp(1))
otu[otu >0] <- 1
number_ASVs <- rowSums(otu)
meta_shannon <- cbind(sample_df,shannon,number_ASVs)
meta_shannon <- meta_shannon %>% subset(Age !="sub")

model <- lm(shannon ~ Yaws * Social_Status * Age * Site * Season, data = meta_shannon)
summary(model)

my_comparisons <- list(c("no", "yes"))


alpha_shannon_treatment_without <- ggplot(meta_shannon, aes(x = Site, y =shannon, color=Yaws)) + 
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8))+
  xlab("Site") +
  ylab("Shannon") +
  theme_bw() +
  stat_compare_means(aes(label = paste0("p =", ..p.format..))) +
  labs(title="Shannon") +
  scale_color_manual(values=c("darkgreen","orangered3"))

lm <- lm(shannon ~ Yaws + Site + Social_Status + Age,  data = meta_shannon)
anova(lm)
TukeyHSD(aov(lm))

lm <- lm(shannon ~ Yaws + Site + Social_Status,  data = meta_shannon)
anova(lm)
TukeyHSD(aov(lm))

lm <- lm(shannon ~ Yaws * Site * Social_Status * Age,  data = meta_shannon)
anova(lm)
TukeyHSD(aov(lm))

lm <- lm(shannon ~ Yaws + Site + Social_Status + Age + sampleID,  data = meta_shannon)
anova(lm)
TukeyHSD(aov(lm))

lm <- lm(shannon ~ Yaws * Site * Social_Status * Age,  data = meta_shannon)
anova(lm)
TukeyHSD(aov(lm))

########################################################




################################################################################################################
###                                               TEMPORAL ANALYSIS                                          ###
################################################################################################################


########################################################
### Chronological timeline
########################################################

otu <- otu_table(project_data)
tax <- as.data.frame(tax_table(project_data))
otu_tax.df <- merge(otu,tax, by= "row.names")[,c("Rank8", colnames(otu))]
row.names(otu_tax.df) <- otu_tax.df$Rank8
otu_tax.df$Rank8 <- NULL
otu_tax.df <- t(otu_tax.df)
colnames(otu_tax.df)

metadata <- as.data.frame(as.matrix(sample_data(project_data)))
metadata$Month_Year <-  paste0(metadata$Month_nb,"-",metadata$Year)

data <- metadata[,c("Date","Gorilla_name","SampleID","Month_Year","Season")]
colnames(data) <- c("start_date","event","SampleID","Month_Year","Season")

#Select length for each moth to improve plot clarity
metadata_month_year <- metadata[,c("Month_Year","SampleID")]
x <- as.data.frame(table(metadata_month_year$Month_Year))
colnames(x) <- c("Month_Year","displ")


data <- merge(data, x, by="Month_Year")

data$start_date <- as.Date(data$start_date, "%d-%m-%Y")
data$displ <- data$displ/10


#Function to shift x-axis to 0 adapted from link shown above

shift_axis <- function(p, xmin, xmax, y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(y=y)
  ax <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  p + annotation_custom(grid::grobTree(ax, vp = grid::viewport(y=1, height=sum(ax$height))), 
                        ymax=y, ymin=y) +
    annotate("segment", y = 0, yend = 0, x = xmin, xend = xmax, 
             arrow = arrow(length = unit(0.1, "inches"))) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank())
  
}


#Conditionally set whether text will be above or below the point
vjust = ifelse(data$displ > 0, -1, 1.5)

#plot
p1 <- data %>% 
  ggplot(aes(start_date, displ)) +
  geom_lollipop(point.size = 1) +
  geom_text(aes(x = start_date, y = displ, label = event), data = data,
            hjust = 0, vjust = vjust, size = 10, angle = 45) +
  theme_light() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 40)) +
  expand_limits(x = c(ymd(20000101), ymd(20150310)), y = 1.2) +
  scale_x_date(breaks = scales::pretty_breaks(n = 20))

#and run the function from above
timeline <- shift_axis(p1, ymd(20000101), ymd(20150310))

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/temporal_analysis/chronology.pdf", height =  10, width = 40)
print(timeline)
dev.off()

########################################################


########################################################
###              Laura's suggestion plot             ###
########################################################

#Build sample table/Rank8
otu <- otu_table(project_data)
tax <- as.data.frame(tax_table(project_data))
otu_tax.df <- merge(otu,tax, by= "row.names")[,c("Rank8", colnames(otu))]
row.names(otu_tax.df) <- otu_tax.df$Rank8
otu_tax.df$Rank8 <- NULL
otu_tax.df <- t(otu_tax.df)
colnames(otu_tax.df)

#Add metadata with gorilla names, sampleID and dates
data <- as.data.frame(as.matrix(sample_data(project_data)))
data <- subset(data, Gorilla_name %in% c("Arthur","Kokima","Ribambel","Tino","Bangwetu","Diyo_L15","Minimee_L61","Murphy","Papillon","Riba","Robin_GS14","Vidole","Yugos"))
otu_tax_agg_meta.df <- merge(otu_tax.df, data, by="row.names", all.X=FALSE)
otu_tax_agg_meta.df <- otu_tax_agg_meta.df[,c(colnames(otu_tax.df),"Gorilla_name","SampleID","Date")]
otu_tax_agg_meta.df <- aggregate(otu_tax_agg_meta.df, by = list(otu_tax_agg_meta.df$Gorilla_name, otu_tax_agg_meta.df$Date), FUN=max)[,-c(1,2)]

otu_tax_agg_meta_melt.df <- melt(otu_tax_agg_meta.df, id.vars=c("Gorilla_name","SampleID","Date"), value.name = "Taxa")
otu_tax_agg_meta_melt.df <- otu_tax_agg_meta_melt.df[otu_tax_agg_meta_melt.df$Taxa>0,]

otu_tax_agg_meta_melt.df$Date <- as.Date(otu_tax_agg_meta_melt.df$Date, "%d-%m-%Y")
otu_tax_agg_meta_melt.df$connect <- paste0(otu_tax_agg_meta_melt.df$Gorilla_name,"_",otu_tax_agg_meta_melt.df$variable)


##Calculate the number of days between each collected fecal sample for each eukaryote and each gorilla
go_name <- unique(otu_tax_agg_meta_melt.df$Gorilla_name)
tax <- unique(otu_tax_agg_meta_melt.df$variable)
stat_list <- list()
for (i in 1:length(tax)){
  for (j in 1:length(go_name)) {
    table <- otu_tax_agg_meta_melt.df[otu_tax_agg_meta_melt.df$variable==tax[i],]
    table <- table[table$Gorilla_name==go_name[j],]
    date_recent <- sort(table$Date,decreasing = TRUE)[1]
    date_old <- sort(table$Date,decreasing = FALSE)[1]
    date_diff <- time_length(difftime(date_recent,date_old), "days")
    table_taxa <- data.frame(Taxa = tax[i],
                             gorilla_name = go_name[j],
                             date_diff = date_diff)
    stat_list[[length(stat_list)+1]] <- table_taxa
  }
}
summary_stat <- data.frame()
for (i in 1:length(stat_list)){
  summary_stat <- rbind(summary_stat,stat_list[[i]])
}
summary_stat[is.na(summary_stat)] <- 0
summary_stat <- subset(summary_stat, date_diff >0 )

#summary_stat_months <- subset(summary_stat, date_diff <365 )
#summary_stat_months$scale <- "Months"
#summary_stat_months$date_diff <- summary_stat_months$date_diff/30

#summary_stat_years <- subset(summary_stat, date_diff >365 )
#summary_stat_years$scale <- "Years"
#summary_stat_years$date_diff <- summary_stat_years$date_diff/365

#summary_stat <- rbind(summary_stat_months,summary_stat_years)


tax <- unique(summary_stat$Taxa)
stat_list <- list()
for (i in 1:length(tax)) {
  max_date <- sort(summary_stat[summary_stat$Taxa == tax[i],3], decreasing = TRUE)[1]
  table_taxa_max <- data.frame(Taxa = tax[i],
                               max_date = max_date)
  stat_list[[length(stat_list)+1]] <- table_taxa_max
}
summary_stat_max <- data.frame()
for (i in 1:length(stat_list)){
  summary_stat_max <- rbind(summary_stat_max,stat_list[[i]])
}

summary_stat_max_2 <- merge(summary_stat, summary_stat_max, by="Taxa")
colnames(summary_stat_max_2)[3] <- "Days"
#

# The palette with grey:
colours <- c("mediumorchid1",
             "navyblue",
             "grey70",
             "cyan",
             "gray33",
             "sienna1",
             "red",
             "forestgreen",
             "firebrick",
             "chartreuse2",
             "plum",
             "blue",
             "chocolate")

summary_stat_max_2$Taxa <- gsub("_"," ",summary_stat_max_2$Taxa)

p2 <- ggplot(summary_stat_max_2, aes(x = reorder(Taxa, max_date), Days)) + 
  geom_point(aes(color = gorilla_name)) +
  scale_colour_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  coord_flip() +
  labs(y="Days between samples collection") +
  theme(axis.text.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=15),
        legend.position="none")


summary_stat_max_3 <- merge(otu_tax_agg_meta_melt.df,summary_stat_max,by.x="variable", by.y="Taxa")

p1 <- ggplot(summary_stat_max_3, aes(x=reorder(variable, max_date), y=Date)) +
  geom_line(aes(color = Gorilla_name), 
            position=position_dodge(width=c(0.7)), size=0.8) +
  scale_colour_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  coord_flip() +
  labs(y="Time interval between collection of samples") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=15),
        legend.text = element_text(size=20),
        legend.title =  element_text(size=20))

plot_ggarrange <- plot_grid(p2,p1, align = "h", ncol = 2, rel_widths = c(2/5, 3/5))
ggsave(filename = "rank8.pdf", path = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/temporal_analysis", plot= plot_ggarrange, width = 15, height = 10, device = "pdf" )

########################################################

########################################################
###              Stay not Stay plot                  ###
########################################################

#Build sample table/Rank8
otu <- otu_table(project_data)
tax <- as.data.frame(tax_table(project_data))
otu_tax.df <- merge(otu,tax, by= "row.names")[,c("Rank8", colnames(otu))]
row.names(otu_tax.df) <- otu_tax.df$Rank8
otu_tax.df$Rank8 <- NULL
otu_tax.df <- t(otu_tax.df)
colnames(otu_tax.df)

#Add metadata with gorilla names, sampleID and dates
data <- as.data.frame(as.matrix(sample_data(project_data)))
data <- subset(data, Gorilla_name %in% c("Arthur","Kokima","Ribambel","Tino","Bangwetu","Diyo_L15","Minimee_L61","Murphy","Papillon","Riba","Robin_GS14","Vidole","Yugos"))

#Calculate number of day between collection of fecal samples per gorillas
metadata <- data[,c("Gorilla_name","Date")]
metadata$Date <- as.Date(metadata$Date, "%d-%m-%Y")
go_name <- unique(metadata$Gorilla_name)
stat_list <- list()
for (j in 1:length(go_name)) {
    table <- metadata[metadata$Gorilla_name==go_name[j],]
    date_recent <- sort(table$Date,decreasing = TRUE)[1]
    date_old <- sort(table$Date,decreasing = FALSE)[1]
    date_diff <- time_length(difftime(date_recent,date_old), "days")
    table_taxa <- data.frame(Gorilla_name = go_name[j],
                             Days_diff = date_diff)
    stat_list[[length(stat_list)+1]] <- table_taxa
}

Number_days_per_gorilla <- data.frame()
for (i in 1:length(stat_list)){
  Number_days_per_gorilla <- rbind(Number_days_per_gorilla,stat_list[[i]])
}


otu_tax_agg_meta.df <- merge(otu_tax.df, data, by="row.names", all.X=FALSE)
otu_tax_agg_meta.df <- otu_tax_agg_meta.df[,c(colnames(otu_tax.df),"Gorilla_name","SampleID","Date")]
otu_tax_agg_meta.df <- aggregate(otu_tax_agg_meta.df, by = list(otu_tax_agg_meta.df$Gorilla_name, otu_tax_agg_meta.df$Date), FUN=max)[,-c(1,2)]
otu_tax_agg_meta_melt.df <- melt(otu_tax_agg_meta.df, id.vars=c("Gorilla_name","SampleID","Date"), value.name = "Taxa")
otu_tax_agg_meta_melt.df$PresenceAbsence <- ifelse(otu_tax_agg_meta_melt.df$Taxa >0, "Presence","Absence")
otu_tax_agg_meta_melt.df$Date <- as.Date(otu_tax_agg_meta_melt.df$Date, "%d-%m-%Y")
otu_tax_agg_meta_melt.df$Connect <- paste0(otu_tax_agg_meta_melt.df$Gorilla_name,"_",otu_tax_agg_meta_melt.df$variable)

##Calculate the number of days between each collected fecal sample for each eukaryote and each gorilla
go_name <- unique(otu_tax_agg_meta_melt.df$Gorilla_name)
tax <- unique(otu_tax_agg_meta_melt.df$variable)
stat_list <- list()
for (i in 1:length(tax)){
  for (j in 1:length(go_name)) {
    table <- otu_tax_agg_meta_melt.df[otu_tax_agg_meta_melt.df$variable==tax[i],]
    table <- table[table$Gorilla_name==go_name[j],]
    x <- length(unique(table$PresenceAbsence))
    table_taxa <- data.frame(Taxa = tax[i],
                             Gorilla_name = go_name[j],
                             date_diff = x )
    stat_list[[length(stat_list)+1]] <- table_taxa
  }
}
summary_stat <- data.frame()
for (i in 1:length(stat_list)){
  summary_stat <- rbind(summary_stat,stat_list[[i]])
}
summary_stat$Connect <- paste0(summary_stat$Gorilla_name,"_",summary_stat$Taxa)

otu_tax_agg_meta_melt.df <- otu_tax_agg_meta_melt.df[,c("PresenceAbsence","Connect")]
test <- merge(summary_stat,otu_tax_agg_meta_melt.df,by="Connect")
test <- unique(test)
test <- merge(test,Number_days_per_gorilla,by="Gorilla_name")
test <- test %>% mutate(Status = case_when(
  date_diff == 2 ~ "Not_stay",
  date_diff == 1  & PresenceAbsence == "Absence" ~ "Absence",
  date_diff == 1  & PresenceAbsence == "Presence" ~ "Presence"
))


test <- subset(test, Status !="Absence")


# The palette with grey:
colours <- c("mediumorchid1",
             "navyblue",
             "yellow",
             "cyan",
             "gray33",
             "sienna1",
             "red",
             "forestgreen",
             "firebrick",
             "chartreuse2",
             "plum",
             "blue",
             "chocolate")

summary_stat_max_2 <- test 
to_keep <- unique(summary_stat_max_2$Taxa)
prevalence_taxa <- d1[d1$Rank8 %in% to_keep,]
summary_stat_max_2 <-merge(summary_stat_max_2, prevalence_taxa, by.x = "Taxa",by.y="Rank8")

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/temporal_analysis/stay_not_stay.pdf", height =  8, width = 15)
p2 <- ggplot(summary_stat_max_2, aes(x = reorder(Taxa,Prevalence) , y = Days_diff, shape=Gorilla_name)) + 
  geom_point(aes(color = as.factor(Status)), size = 5) +
  scale_colour_manual(values=colours) +
  scale_shape_manual(values=c(1:20))+
  scale_y_continuous(trans='log10') +
  theme_half_open(line_size = 0.5) +
  background_grid() +
  coord_flip() +
  labs(y="Days between samples collection") +
  theme(axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=15))
print(p2)
dev.off()


summary_stat_max_3 <- merge(otu_tax_agg_meta_melt.df,summary_stat_max,by.x="variable", by.y="Taxa")

p1 <- ggplot(summary_stat_max_3, aes(x=reorder(variable, max_date), y=Date)) +
  geom_line(aes(color = Gorilla_name), 
            position=position_dodge(width=c(0.7)), size=0.8) +
  scale_colour_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  coord_flip() +
  labs(y="Time interval between collection of samples") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=15),
        legend.text = element_text(size=20),
        legend.title =  element_text(size=20))

plot_ggarrange <- plot_grid(p2,p1, align = "h", ncol = 2, rel_widths = c(2/5, 3/5))
ggsave(filename = "rank8.pdf", path = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/temporal_analysis", plot= plot_ggarrange, width = 15, height = 10, device = "pdf" )

########################################################



########################################################
### Beta diversity
########################################################


########################################################


########################################################
###               Alpha diversity                    ###
########################################################

## Create Colour Palettes ###
n <- length(unique(sample_data(project_data_rarefied)$gorilla_name))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-16]

#Gorilla name
project_data_exp <- project_data_rarefied %>% subset_samples(temporal == "yes") 
alpha_Season <- plot_richness(project_data_exp, x = "gorilla_name", color = "gorilla_name", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_point() +
  scale_color_manual(values=col_vector)

project_data_exp <- project_data_rarefied %>% subset_samples(temporal == "yes") 
alpha_Season <- plot_richness(project_data_exp, x = "gorilla_name", color = "gorilla_name", measures = c("Shannon")) +
  geom_point() +
  scale_color_manual(values=col_vector)
########################################################




################################################################################################################
################################################################################################################
################################################################################################################
###                                                                                                          ###
###                      ANALYSIS by COMBINING SAMPLES FROM MULTIPLE TIME POINT                              ###
###                                                                                                          ###
################################################################################################################
################################################################################################################
################################################################################################################

########################################################
## Extract taxonomy, ASV and control tables from data ##
########################################################

#Select only gut eukaryotes
asv_gut_rnk_spl_neg_df2 <- subset(asv_rnk_spl_neg_df2, Symbiont %in% c("gut_worm","gut_protist"))

### Extract taxonomy table ###
names_rnk <- c("Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Rank8","Accession_number")
asv_tax.df <- asv_gut_rnk_spl_neg_df2[,names_rnk]
#Because of presence of NA, need to rename unknwown taxa at each rank (once printing makes easier selection of taxa to create table)
asvUn_tax.df <- asv_tax.df                          #Create new table  
asvUn_tax.df[is.na(asvUn_tax.df)] <- "Un"
#create a function to paste "Un" to the name of the previous rank
mypaste = function(s) {
  s = gsub("_Un","",s)
  return(paste0(s,"_Un"))
}
#apply the function to all ranks
for (i in 2:9){
  
  if (i == 9)
  {
    colName = "Accession_number"
    colNamePrevious = "Rank8"
  } else
  {
    colName = paste0("Rank", i)
    colNamePrevious = paste0("Rank", i-1)
  }
  w = which(asvUn_tax.df[[colName]] == "Un")
  asvUn_tax.df[[colName]][w] = mypaste(asvUn_tax.df[[colNamePrevious]][w])
}
asvUn_tax.mtx <- as.matrix(asvUn_tax.df);dim(asvUn_tax.mtx)

### Extract negative controls tables ###
names_neg <- c("neg1","neg2","neg3","neg4")
asv_neg.df <- asv_gut_rnk_spl_neg_df2[,names_neg]
class(asv_neg.df);head(asv_neg.df);dim(asv_neg.df);str(asv_neg.df)

### Extract ASV-Samples table ###
names_spl_only <- colnames(asv_gut_rnk_spl_neg_df2)[grep("G",colnames(asv_gut_rnk_spl_neg_df2))]
names_spl <- c("ASV_cluster",names_spl_only)
asv_spl.df <- asv_gut_rnk_spl_neg_df2[,names_spl]
asv_spl.df[,-1] <- asv_spl.df[,-1] %>% mutate_all(as.numeric)
asv_spl.df <- as.data.frame(asv_spl.df %>% group_by(ASV_cluster) %>% summarise_each(funs(sum)))
row.names(asv_spl.df) <- asv_spl.df$ASV_cluster
asv_spl.df$ASV_cluster <- NULL
asv_spl.df$Group.1 <- NULL
head(asv_spl.df); dim(asv_spl.df) ; str(asv_spl.df)

########################################################

########################################################
###               Filtering steps                    ###
########################################################

#Merge sample that belong to similar gorilla name (avoid pseudo replication)
metadata_spl_ctl_GorillaName <- metadata_spl_ctl_df[,c("gorilla_name","Cluster_sampleID")]
spl_asv.df <- as.data.frame(t(asv_spl.df))
spl_asv_GorillaName <- merge(metadata_spl_ctl_GorillaName,spl_asv.df, by="row.names")
spl_asv_GorillaName <- spl_asv_GorillaName[,-c(1,2)]
spl_asv_GorillaName_agg <- as.data.frame(spl_asv_GorillaName %>% group_by(Cluster_sampleID) %>% summarise_each(funs(sum)))
row.names(spl_asv_GorillaName_agg) <- spl_asv_GorillaName_agg$Cluster_sampleID
spl_asv_GorillaName_agg$Cluster_sampleID <- NULL
spl_asv_GorillaName_agg <- as.data.frame(t(spl_asv_GorillaName_agg))


#Create phyloseq object
metadata <- sample_data(metadata_spl_ctl_df)
otu <- otu_table(spl_asv_GorillaName_agg, taxa_are_rows=TRUE)
tax <- tax_table(asvUn_tax.mtx)
project_data <- phyloseq(metadata,otu,tax);project_data

project_data <- project_data <- project_data %>% 
  subset_samples(sampleID !="G0295" );project_data

########################################################

################################################################################################################
###                              Prevalence and Abundance of intestinal eukaryotes                           ###
################################################################################################################

########################################################
###   Calculate Prevalence, and Abundance of reads   ###
########################################################

#---- Create ASV table for ----#
asv_spl.df <- as.data.frame(otu_table(project_data) )#create ASV table
asv_tax.df <- as.data.frame(tax_table(project_data))
rownames(asv_spl.df)==rownames(asv_tax.df)           #Sanity check to see if ASV from both table are similar and same order
asv_tax_spl.df <- cbind(asv_tax.df,asv_spl.df)       #Bind both table to have ASV and TAX table together

#---- Build prevalence table at Rank8 ----#
asv_tax_spl_rnk8 <- asv_tax_spl.df[,setdiff(colnames(asv_tax_spl.df),c(names_rnk[-8]))]                 #Select table with samples and rank8 only
row.names(asv_tax_spl_rnk8) <- asv_tax_spl_rnk8$Rank8
asv_tax_spl_rnk8$Rank8 <- NULL
Total_abundance <- rowSums(asv_tax_spl_rnk8)
Total_abundance_log <- round(log(Total_abundance, base=10),3)
#Calculate prevalence
asv_tax_spl_rnk8_bin <- asv_tax_spl_rnk8
asv_tax_spl_rnk8_bin[] <- +(asv_tax_spl_rnk8_bin >= 1) #Convert table into binary 1:0 (Presence absence)
Prevalence <- round(rowSums(asv_tax_spl_rnk8_bin)/length(colnames(asv_tax_spl_rnk8_bin))*100,2)
#Combine Prevelance, Total_abundance, Total_abundance_log to tax table
d1 <- data.frame(Rank8 = row.names(asv_tax_spl_rnk8), 
                 Prevalence = Prevalence,
                 Total_abundance = Total_abundance,
                 Total_abundance_log = Total_abundance_log)
d1$Rank8 <- as.character(d1$Rank8)

########################################################

########################################################
###           Calculate Relative abundance           ###
########################################################

#---- Create ASV table for ----#
asv_spl.df <- as.data.frame(otu_table(project_data) )#create ASV table
asv_tax.df <- as.data.frame(tax_table(project_data))
rownames(asv_spl.df)==rownames(asv_tax.df)           #Sanity check to see if ASV from both table are similar and same order
asv_tax_spl.df <- cbind(asv_tax.df,asv_spl.df)       #Bind both table to have ASV and TAX table together

#---- Calculate relative abundance
asv_tax_spl_rnk8 <- asv_tax_spl.df[,setdiff(colnames(asv_tax_spl.df),c(names_rnk[-8]))]                 #Select table with samples and rank8 only
row.names(asv_tax_spl_rnk8) <- asv_tax_spl_rnk8$Rank8
asv_tax_spl_rnk8$Rank8 <- NULL
asv_tax_spl_rnk8 <- as.data.frame(t(asv_tax_spl_rnk8))
asv_tax_spl_rnk8 <- asv_tax_spl_rnk8 %>% mutate_all(as.numeric)
asv_tax_spl_rnk8_rel <- asv_tax_spl_rnk8/rowSums(asv_tax_spl_rnk8)
asv_tax_spl_rnk8_rel_melt <- as.data.frame(t(asv_tax_spl_rnk8_rel))
asv_tax_spl_rnk8_rel_melt$Rank8 <- row.names(asv_tax_spl_rnk8_rel_melt) 
asv_tax_spl_rnk8_rel_melt <- melt(asv_tax_spl_rnk8_rel_melt,id.vars=c("Rank8"))

########################################################

########################################################
###     Filter tree using only ASV from ASV table    ###
########################################################

asv_to_drop <- setdiff(new_tree$tip.label,intersect(row.names(asv_tax.df),new_tree$tip.label))
new_tree <- drop.tip(new_tree, asv_to_drop)
plot(new_tree)
asv_tax.df$ASV_cluster <- row.names(asv_tax.df)
asv_tax.df <- asv_tax.df %>% arrange(factor(ASV_cluster, levels = new_tree$tip.label))
new_tree$tip.label == asv_tax.df$ASV_cluster
new_tree$tip.label <- asv_tax.df$Rank8
plot(new_tree)

########################################################

########################################################
###                  Create figure                   ###
########################################################

#Create tree
g <- ggtree(new_tree, branch.length="none")+ 
  geom_tiplab(size=6)+
  ggplot2::xlim(0, 20)

#Create Prevalence barplot
p1 <- ggplot(d1, aes(x=Rank8, y=Prevalence)) + 
  geom_col(width = 0.5) + 
  coord_flip() + 
  theme_half_open() +
  background_grid() +
  scale_y_continuous(limits=c(0, 100), breaks=c(0, 25, 50, 75, 100)) +
  labs(y="Prevalence  (%)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

#Create boxplot for abundance distribution log-transformed
p2b_log <- ggplot(asv_tax_spl_rnk8_rel_melt, aes(x=Rank8, y=value)) + 
  geom_boxplot(size = 0.5, width = 0.5)+
  scale_y_continuous(trans='log10') +
  geom_jitter(position=position_dodge(0.8), size = 1) +
  coord_flip() +
  theme_half_open() +
  background_grid() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.text.y = element_blank())

#Save figure in directory
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/prevalence_abundance/phylo_all_prevalence_abundance3.pdf", height =  10, width = 10)
p1 %>% insert_left(g, width = 2) 
dev.off()

#Save figure in directory
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/prevalence_abundance/phylo_all_rel_abun_across_samples.pdf", height =  10, width = 13)
p2b_log %>% insert_left(g, width = 2) 
dev.off()

########################################################

################################################################################################################
###                      Relative abundance for protist, eukaryote, and protist/eukaryote                    ###
################################################################################################################


########################################################
### Relative abundance within Protists Rank8
########################################################

sort(sample_sums(project_data))
project_data_rarefied <- rarefy_even_depth(project_data, sample.size = 14418);project_data

asv_protist.df <- subset(asv_rnk_spl_neg_df2, Symbiont %in% c("gut_protist"))
asv_gut_protist <- row.names(asv_protist.df)
#asv_blastocystis <- c("ASV28","ASV45","ASV1146","ASV117","ASV163","ASV35")

project_data_protist <- prune_taxa(asv_gut_protist, project_data_rarefied);project_data_protist

otu <- t(otu_table(project_data_protist))
otu_rel <- otu/rowSums(otu)
otu_rel <- as.data.frame(otu_rel[,order(colSums(otu_rel), decreasing = TRUE)])
rowSums(otu_rel)  #sanity check, all value should be equal to 1
colSums(otu_rel)  #should be sorted from the highest to the lowest
otu_rel$sample_id <- row.names(otu_rel)
otu_rel_melt = melt(otu_rel, id = c("sample_id"))
otu_abs <- as.data.frame(otu)
otu_abs$sample_id <- row.names(otu_abs)
otu_abs_melt = melt(otu_abs, id = c("sample_id"))
otu_abs_rel_melt <-  merge(otu_rel_melt,otu_abs_melt,by=c("sample_id", "variable"))

#Extract taxonomy table
tax <- as.data.frame(tax_table(project_data_protist))
#Extract metadata table
data <- as.data.frame(as.matrix(sample_data(project_data_protist)))

#Merge tables (otu, taxand data)
otu_abs_rel_melt_tax <- merge(otu_abs_rel_melt,tax, by.x= "variable", by.y="row.names",all.x=TRUE, all.y=FALSE)
otu_abs_rel_melt_tax_data <- merge(otu_abs_rel_melt_tax,data,data, by.x="sample_id",by.y="row.names")

#Create column with presence absence of ASV
otu_abs_rel_melt_tax_data$value_bin <- ifelse(otu_abs_rel_melt_tax_data$value.x >0 , 1, 0) 

otu_abs_rel_melt_tax_data$gorilla_name_sampleID <- paste0(otu_abs_rel_melt_tax_data$gorilla_name, "_", otu_abs_rel_melt_tax_data$sampleID)
otu_abs_rel_melt_tax_data$Rank8 <- gsub("_"," ",otu_abs_rel_melt_tax_data$Rank8 )
otu_abs_rel_melt_tax_data$Rank8 <- factor(otu_abs_rel_melt_tax_data$Rank8, 
                                          level= c("Blastocystis ST1","Blastocystis ST11","Blastocystis ST2","Blastocystis ST3","Blastocystis Unknown ASV28","Blastocystis Unknown ASV45",
                                                   "Diplomonad Unknown","Enteromonad Unknown ASV663","Enteromonad Unknown ASV73","Enteromonas hominis","Trepomonas agilis","Trepomonas Unknown ASV294","Retortamonas sp. CladeA",
                                                   "Charonina ventriculi","Cycloposthium bipalmatum","Helicozoster indicus","Hemiprorodon gymnoprosthium","Latteuria media","Parentodinium sp. ASV132","Parentodinium sp. ASV76","Pseudoentodinium elephantis","Triplumaria fulgora",
                                                   "Entamoeba coli ST1","Entamoeba coli ST2","Entamoeba hartmanni",
                                                   "Tetratrichomonas buttreyi","Tetratrichomonas prowazecki"))

#otu_abs_rel_melt_tax_data$Rank8 <- factor(otu_abs_rel_melt_tax_data$Rank8, 
#                                          level= c("Blastocystis_ST1","Blastocystis_ST11","Blastocystis_ST2","Blastocystis_ST3","Blastocystis_Unknown_ASV28","Blastocystis_Unknown_ASV45"))

otu_abs_rel_melt_tax_data$gorilla_name <- gsub("_"," ",otu_abs_rel_melt_tax_data$gorilla_name)

otu_abs_rel_melt_tax_data$Site <- factor(otu_abs_rel_melt_tax_data$Site, 
                                         level= c("Lokoue","Romani","Moba","Maya"))

colours <- c("purple4","purple","hotpink","pink","violet","palevioletred1","navy","dodgerblue","cyan","blue","lightcyan1","blue4","dodgerblue","olivedrab4","greenyellow","khaki1","aquamarine4","gold3","darkseagreen1","green4","seagreen2","yellow","grey10","grey50","grey60","coral4","peru")

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/Taxa_abundance_rank8_protist_individual.pdf", width = 15, height = 12)
Taxa_abundance_rank8_protist <- ggplot(otu_abs_rel_melt_tax_data, aes(x = gorilla_name, y= value.x, fill=Rank8, color = colours)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid", size=0.01) + #height of the column equal the value, stacked
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=15),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=15),
        legend.position = "left",
        legend.spacing.y = unit(0.2, 'cm'),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        strip.text.x = element_text(size = 20)) +
  guides(fill=guide_legend(ncol=1,byrow = TRUE))
print(Taxa_abundance_rank8_protist)
dev.off()

########################################################

########################################################
### Nematoda community and Protist
########################################################

sort(sample_sums(project_data))
project_data_rarefied <- rarefy_even_depth(project_data, sample.size = 14418);project_data

otu <- t(otu_table(project_data_rarefied))
otu_rel <- otu/rowSums(otu)
otu_rel <- as.data.frame(otu_rel[,order(colSums(otu_rel), decreasing = TRUE)])
rowSums(otu_rel)  #sanity check, all value should be equal to 1
colSums(otu_rel)  #should be sorted from the highest to the lowest
otu_rel$sample_id <- row.names(otu_rel)
otu_rel_melt = melt(otu_rel, id = c("sample_id"))
otu_abs <- as.data.frame(otu)
otu_abs$sample_id <- row.names(otu_abs)
otu_abs_melt = melt(otu_abs, id = c("sample_id"))
otu_abs_rel_melt <-  merge(otu_rel_melt,otu_abs_melt,by=c("sample_id", "variable"))

#Extract taxonomy table
tax <- as.data.frame(tax_table(project_data_rarefied))
#Extract metadata table
data <- as.data.frame(as.matrix(sample_data(project_data_rarefied)))

#Merge tables (otu, taxand data)
otu_abs_rel_melt_tax <- merge(otu_abs_rel_melt,tax, by.x= "variable", by.y="row.names",all.x=TRUE, all.y=FALSE)
otu_abs_rel_melt_tax_data <- merge(otu_abs_rel_melt_tax,data,data, by.x="sample_id",by.y="row.names")

#Create column with presence absence of ASV
otu_abs_rel_melt_tax_data$value_bin <- ifelse(otu_abs_rel_melt_tax_data$value.x >0 , 1, 0) 

#Convert all protist lineAges into "protist" category
Protist_lineAges <- c("Blastocystis_ST1|Blastocystis_ST11|Blastocystis_ST2|Blastocystis_ST3|Blastocystis_Unknown_ASV28|Blastocystis_Unknown_ASV45|Diplomonad_Unknown|Enteromonad_Unknown_ASV663|Enteromonad_Unknown_ASV73|Enteromonas_hominis|Trepomonas_agilis|Trepomonas_Unknown_ASV294|Retortamonas_sp._CladeA|Charonina_ventriculi|Cycloposthium_bipalmatum|Helicozoster_indicus|Hemiprorodon_gymnoprosthium|Latteuria_media|Parentodinium_sp._ASV132|Parentodinium_sp._ASV76|Pseudoentodinium_elephantis|Triplumaria_fulgora|Entamoeba_coli_ST1|Entamoeba_coli_ST2|Entamoeba_hartmanni|Tetratrichomonas_buttreyi|Tetratrichomonas_prowazecki")
otu_abs_rel_melt_tax_data$Rank8 <- gsub(Protist_lineAges,"Protist",otu_abs_rel_melt_tax_data$Rank8)

otu_abs_rel_melt_tax_data$gorilla_name_sampleID <- paste0(otu_abs_rel_melt_tax_data$gorilla_name, "_", otu_abs_rel_melt_tax_data$sampleID)
otu_abs_rel_melt_tax_data$Rank8 <- gsub("_"," ",otu_abs_rel_melt_tax_data$Rank8)
otu_abs_rel_melt_tax_data$Rank8 <- factor(otu_abs_rel_melt_tax_data$Rank8, 
                                          level= c("Strongylida ASV2","Strongylida ASV1","Strongylida ASV177","Strongylida ASV154",
                                                   "Metastrongyloidea ASV482","Metastrongyloidea ASV152","Metastrongyloidea ASV204",
                                                   "Strongyloides fuelleborni","Spiruroidea ASV58","Spiruroidea ASV311","Abbreviata ASV495","Ascaridia ASV389","Nematoda sp. ASV29","Protist"))


otu_abs_rel_melt_tax_data$Site <- factor(otu_abs_rel_melt_tax_data$Site, 
                                         level= c("Lokoue","Romani","Moba","Maya"))

colours <- c("orange","red","tomato","red4","green","darkgreen","greenyellow","gold3","aquamarine4","yellow2","grey10","white","dodgerblue","grey70")

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/Taxa_abundance_rank8_protist_nematode_individual.pdf", width = 10, height = 10)
Taxa_abundance_rank8_protist_nematode <- ggplot(otu_abs_rel_melt_tax_data, aes(x = gorilla_name, y= value.x, fill=Rank8, color = colours)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid",size=0.1) + #height of the column equal the value, stacked
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=15),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.spacing.y = unit(0.2, 'cm'),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE))
print(Taxa_abundance_rank8_protist_nematode)
dev.off()

########################################################

########################################################
### Nematoda and Protist
########################################################

sort(sample_sums(project_data))
project_data_rarefied <- rarefy_even_depth(project_data, sample.size = 14418);project_data


otu <- t(otu_table(project_data_rarefied))
otu_rel <- otu/rowSums(otu)
otu_rel <- as.data.frame(otu_rel[,order(colSums(otu_rel), decreasing = TRUE)])
rowSums(otu_rel)  #sanity check, all value should be equal to 1
colSums(otu_rel)  #should be sorted from the highest to the lowest
otu_rel$sample_id <- row.names(otu_rel)
otu_rel_melt = melt(otu_rel, id = c("sample_id"))
otu_abs <- as.data.frame(otu)
otu_abs$sample_id <- row.names(otu_abs)
otu_abs_melt = melt(otu_abs, id = c("sample_id"))
otu_abs_rel_melt <-  merge(otu_rel_melt,otu_abs_melt,by=c("sample_id", "variable"))

#Extract taxonomy table
tax <- as.data.frame(tax_table(project_data_rarefied))
#Extract metadata table
data <- as.data.frame(as.matrix(sample_data(project_data_rarefied)))

#Merge tables (otu, taxand data)
otu_abs_rel_melt_tax <- merge(otu_abs_rel_melt,tax, by.x= "variable", by.y="row.names",all.x=TRUE, all.y=FALSE)
otu_abs_rel_melt_tax_data <- merge(otu_abs_rel_melt_tax,data,data, by.x="sample_id",by.y="row.names")

#Create column with presence absence of ASV
otu_abs_rel_melt_tax_data$value_bin <- ifelse(otu_abs_rel_melt_tax_data$value.x >0 , 1, 0) 

otu_abs_rel_melt_tax_data$Rank4 <- gsub("Entamoebida|Blastocystis|Fornicata|Litostomatea|Parabasalia","Protist",otu_abs_rel_melt_tax_data$Rank4)

otu_abs_rel_melt_tax_data$Site <- factor(otu_abs_rel_melt_tax_data$Site, 
                                         level= c("Lokoue","Romani","Moba","Maya"))

colours  = c("grey69","grey28")

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/Taxa_abundance_rank4_protist_nematode_individual.pdf", width = 10, height = 10)
Taxa_abundance_rank4_protist_nematode <- ggplot(otu_abs_rel_melt_tax_data, aes(x = gorilla_name, y= value.x, fill=Rank4, color = colours)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid",size=0.1) + #height of the column equal the value, stacked
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=15),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  theme(legend.position="none")
print(Taxa_abundance_rank4_protist_nematode)
dev.off()

########################################################

########################################################
### Relative abundance within Nematodes Rank8
########################################################

# Select only worm
asv_worm.df <- subset(asv_rnk_spl_neg_df2, Symbiont %in% c("gut_worm"))
asv_gut_worm <- row.names(asv_worm.df)
project_data_worm <- prune_taxa(asv_gut_worm, project_data_rarefied);project_data_worm

# Buil relative abudance table
otu <- t(otu_table(project_data_worm))
otu_rel <- otu/rowSums(otu)
otu_rel <- as.data.frame(otu_rel[,order(colSums(otu_rel), decreasing = TRUE)])
rowSums(otu_rel)  #sanity check, all value should be equal to 1
colSums(otu_rel)  #should be sorted from the highest to the lowest
otu_rel$sample_id <- row.names(otu_rel)
otu_rel_melt = melt(otu_rel, id = c("sample_id"))
otu_abs <- as.data.frame(otu)
otu_abs$sample_id <- row.names(otu_abs)
otu_abs_melt = melt(otu_abs, id = c("sample_id"))
otu_abs_rel_melt <-  merge(otu_rel_melt,otu_abs_melt,by=c("sample_id", "variable"))

#Extract taxonomy table
tax <- as.data.frame(tax_table(project_data_worm))
#Extract metadata table
data <- as.data.frame(as.matrix(sample_data(project_data_worm)))

#Merge tables (otu, taxand data)
otu_abs_rel_melt_tax <- merge(otu_abs_rel_melt,tax, by.x= "variable", by.y="row.names",all.x=TRUE, all.y=FALSE)
otu_abs_rel_melt_tax_data <- merge(otu_abs_rel_melt_tax,data,data, by.x="sample_id",by.y="row.names")

#Create column with presence absence of ASV
otu_abs_rel_melt_tax_data$value_bin <- ifelse(otu_abs_rel_melt_tax_data$value.x >0 , 1, 0) 

otu_abs_rel_melt_tax_data$gorilla_name_sampleID <- paste0(otu_abs_rel_melt_tax_data$gorilla_name, "_", otu_abs_rel_melt_tax_data$sampleID)
otu_abs_rel_melt_tax_data$Rank8 <- gsub("_"," ",otu_abs_rel_melt_tax_data$Rank8)
otu_abs_rel_melt_tax_data$Rank8 <- factor(otu_abs_rel_melt_tax_data$Rank8, 
                                          level= c("Strongylida ASV2","Strongylida ASV1","Strongylida ASV177","Strongylida ASV154",
                                                   "Metastrongyloidea ASV482","Metastrongyloidea ASV152","Metastrongyloidea ASV204",
                                                   "Strongyloides fuelleborni","Spiruroidea ASV58","Spiruroidea ASV311","Abbreviata ASV495","Ascaridia ASV389","Nematoda sp. ASV29"))

otu_abs_rel_melt_tax_data$Site <- factor(otu_abs_rel_melt_tax_data$Site, 
                                         level= c("Lokoue","Romani","Moba","Maya"))

colours <- c("orange","red","tomato","red4","green","darkgreen","greenyellow","gold3","aquamarine4","yellow2","grey10","dodgerblue","grey70")

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/Taxa_abundance_rank8_worm_individual.pdf", width = 15, height = 12)
Taxa_abundance_rank8_worm <- ggplot(otu_abs_rel_melt_tax_data, aes(x = gorilla_name, y= value.x, fill=Rank8, color = colours)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid", size=0.01) + #height of the column equal the value, stacked
  facet_grid(rows = vars(Site), scales = "free_y", switch = "y", space = "free_y")+ 
  coord_flip() +
  scale_fill_manual(values=colours) +
  theme_half_open() +
  background_grid() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=15),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.spacing.y = unit(0.2, 'cm'),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE))
print(Taxa_abundance_rank8_worm)
dev.off()

########################################################

########################################################
### Combine Plot into Figure 2
########################################################

plot_ggarrange <- ggarrange(Taxa_abundance_rank8_protist,Taxa_abundance_rank4_protist_nematode,Taxa_abundance_rank8_worm, ncol=3, widths = c(1.5,0.3,1.2))
ggsave(filename = "Figure2_individual.pdf", path = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/", plot= plot_ggarrange, width = 25, height = 12, device = "pdf" )

plot_ggarrange <- ggarrange(Taxa_abundance_rank8_protist,Taxa_abundance_rank8_protist_nematode, ncol=2, widths = c(1.3,1))
ggsave(filename = "Figure2_individual_V2.pdf", path = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/community_structure/", plot= plot_ggarrange, width = 25, height = 12, device = "pdf" )


########################################################


################################################################################################################
###                             RELATIONSHIP BETWEEN INTESINAL EUKARYOTES                                    ###
################################################################################################################

### Co-association (Presence/Absence)
########################################################

#Set workinf directory
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test")
sort(sample_sums(project_data))
#Create list for final table with
permute_list <- list()
permutation_vector <- rep(14418,20)
for(perm in 1:length(permutation_vector)){
  project_data_rarefied <- rarefy_even_depth(project_data, sample.size = permutation_vector[perm]);project_data
  #---- Create ASV table for ----#
  asv_spl.df <- as.data.frame(otu_table(project_data_rarefied) )#create ASV table
  asv_tax.df <- as.data.frame(tax_table(project_data_rarefied))
  rownames(asv_spl.df)==rownames(asv_tax.df)           #Sanity check to see if ASV from both table are similar and same order
  asv_tax_spl.df <- cbind(asv_tax.df,asv_spl.df)       #Bind both table to have ASV and TAX table together
  #---- Build prevalence table at Rank8 ----#
  asv_tax_spl_rnk8 <- asv_tax_spl.df[,setdiff(colnames(asv_tax_spl.df),c(names_rnk[-8]))]                 #Select table with samples and rank8 only
  row.names(asv_tax_spl_rnk8) <- asv_tax_spl_rnk8$Rank8
  asv_tax_spl_rnk8$Rank8 <- NULL
  asv_tax_spl_rnk8 <- as.data.frame(t(asv_tax_spl_rnk8))
  asv_tax_spl_rnk8[asv_tax_spl_rnk8<=400] <- 0
  asv_tax_spl_rnk8[asv_tax_spl_rnk8>400] <- 1
  #Only take taxa with at least 5 sample presence  and at least 5 sample absent 
  to_keep <- which(colSums(asv_tax_spl_rnk8) < length(sample_names(project_data_rarefied))-3 & colSums(asv_tax_spl_rnk8) >4)
  asv_tax_spl_rnk8 <- asv_tax_spl_rnk8[,to_keep]
  asv_tax_spl_rnk8 <- ifelse(asv_tax_spl_rnk8 ==1, "Presence", "Absence")
  #Run fisher test for all pairwise comparison
  stat_list <- list()
  for (taxa1 in colnames(asv_tax_spl_rnk8)){
    for (taxa2 in colnames(asv_tax_spl_rnk8)) {
      table_fisher <- table(asv_tax_spl_rnk8[,taxa1],asv_tax_spl_rnk8[,taxa2])
      is_fisher_test <-  
        tryCatch(expr =  { 
          fisher_test <- fisher.test(table_fisher)
          TRUE},
          error=function(e){
            return(FALSE)}
        )
      if(is_fisher_test)
        c(summary_table <- data.frame(
          Taxa1=taxa1,
          Spl_taxa1= length(which(asv_tax_spl_rnk8[,taxa1]=="Presence")),
          Taxa2=taxa2,
          Spl_taxa2= length(which(asv_tax_spl_rnk8[,taxa2]=="Presence")),
          Pvalue=fisher_test$p.value,
          Odds_ratio =round(fisher_test$estimate,3)
        )
        )
      if(!is_fisher_test)
        c(summary_table <- data.frame(
          Taxa1=taxa1,
          Spl_taxa1= length(which(asv_tax_spl_rnk8[,taxa1]=="Presence")),
          Taxa2=taxa2,
          Spl_taxa2= length(which(asv_tax_spl_rnk8[,taxa2]=="Presence")),
          Pvalue=NA,
          Odds_ratio =NA
        )
        )
      stat_list[[length(stat_list)+1]] <- summary_table
    }
  }
  #Generate table that summarize all alpha and beta statistic
  Fisher_test_summary <- data.frame()
  for (i in 1:length(stat_list)){
    Fisher_test_summary <- rbind(Fisher_test_summary,stat_list[[i]])
  }
  Fisher_test_summary$Pvalue_Adjust <- round(p.adjust(Fisher_test_summary$Pvalue, method = "BH", n = length(Fisher_test_summary$Pvalue)),3)
  Fisher_test_summary$duplicate <- Fisher_test_summary$Spl_taxa1==Fisher_test_summary$Spl_taxa2
  Fisher_test_summary <- subset(Fisher_test_summary, duplicate == FALSE)
  Fisher_test_summary$Odds_sup1 <- ifelse(Fisher_test_summary$Odds_ratio>=1, ">1 (Positive association)","<1 (Negative association)")
  Fisher_test_summary$Significance <- ifelse(Fisher_test_summary$Pvalue_Adjust < 0.05,"Significant","Non-significant")
  row.names(Fisher_test_summary) <- NULL
  temp <- Fisher_test_summary[,c("Taxa1","Taxa2")]
  newDf <- data.frame(t(apply(temp,1,sort)))
  Fisher_test_summary <- Fisher_test_summary[!duplicated(newDf),]
  if(dim(Fisher_test_summary)[1]>0){
    permute_list[[length(permute_list)+1]] <- Fisher_test_summary}
  #Plot histogram
  hist_plot <- ggplot(Fisher_test_summary, aes(x=Odds_ratio)) +
    geom_histogram(binwidth = 0.2)+
    geom_vline(xintercept=1, linetype='dashed') +
    theme_light() +
    theme(axis.title.y = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=15),
          axis.text.x = element_text(size=15))
  #Plot Bar plot number of odd ratio below and above 1
  bar_plot <- ggplot(Fisher_test_summary, aes(x=Odds_sup1)) +
    geom_bar() +
    theme_light() +
    theme(axis.title.y = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=15),
          axis.text.x = element_text(size=9))
  #Scatter plot, Odds ratio and Pvalue
  bubble_plot <- ggplot(Fisher_test_summary,aes(x=Odds_ratio, y=Pvalue_Adjust, color=Significance)) +
    geom_point() +
    xlim(0,20) +
    ylim(1,0.0001) +
    theme_light() +
    geom_vline(xintercept=1, linetype='dashed') +
    geom_hline(yintercept=0.05, linetype='dashed') +
    theme(axis.title.y = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=15),
          axis.text.x = element_text(size=15))
  combined <- ggarrange(
    bubble_plot,                # First row with line plot
    widths = c(1, 0.2),
    # Second row with box and histogram
    ggarrange(bar_plot, hist_plot, ncol = 2, labels = c("B", "C"), widths = c(0.7, 1)), 
    nrow = 2, 
    labels = "A"       # Label of the line plot
    )
  ggsave(combined,filename = paste0(perm,"_","fisher_test_multipannel.pdf"), width = 10, height = 10)
}

Fisher_permute_summary <- data.frame()
for (i in 1:length(permute_list)){
  Fisher_permute_summary <- rbind(Fisher_permute_summary,permute_list[[i]])
}

#Add additional information to table
Fisher_permute_summary$nb_significant <- 1 #To later calculate number of time pairwise was significant
Fisher_permute_summary$pair <- paste0(Fisher_permute_summary$Taxa1,"_",Fisher_permute_summary$Taxa2)
Fisher_permute_summary$pairwise <- ifelse(Fisher_permute_summary$Significance=="Significant",paste0(Fisher_permute_summary$Taxa1,"_",Fisher_permute_summary$Taxa2),"NS")
to_keep <- unique(Fisher_permute_summary$pairwise)
to_keep <- to_keep[!to_keep %in% "NS"]
Fisher_permute_summary_sg <- subset(Fisher_permute_summary, pair %in% to_keep)
Fisher_permute_summary_sg$pairwise <- Fisher_permute_summary_sg$pair
to_keep2 <- setdiff(Fisher_permute_summary$pair, to_keep)
Fisher_permute_summary_nsg <- subset(Fisher_permute_summary, pair %in% to_keep2)
Fisher_permute_summary <- rbind(Fisher_permute_summary_sg,Fisher_permute_summary_nsg)

#Plot histogram
hist_plot <- ggplot(Fisher_permute_summary, aes(x=Odds_ratio)) +
  geom_histogram(binwidth = 0.2)+
  geom_vline(xintercept=1, linetype='dashed') +
  theme_light() +
  theme(axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=15))

#Plot Bar plot number of odd ratio below and above 1
bar_plot <- ggplot(Fisher_permute_summary, aes(x=Odds_sup1)) +
  geom_bar() +
  theme_light() +
  theme(axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=9))

#Bubble plot for proportion of time test is significant
Fisher_permute_summary_sig <- subset(Fisher_permute_summary, Pvalue_Adjust <0.05) # Select only signiifcant adjusted pvalue
Fisher_permute_summary_groupby <- Fisher_permute_summary_sig %>% 
  group_by(Taxa1,Taxa2) %>%
  summarise(avg_adpvl = mean(Pvalue_Adjust),    #Calculate mean of pvalue 
            Mean_OddRatio = mean(Odds_ratio),   #Calculate mean of odds ratio
            Permutation = sum(nb_significant))  #Calculate total number of time the pairwise is significant
Fisher_permute_summary_groupby$Permutation <- Fisher_permute_summary_groupby$Permutation/perm*100
setorder(Fisher_permute_summary_groupby, cols= - "Permutation") #Reoerder table by the number of time it was significant
order_legend_pairwise <- paste0(Fisher_permute_summary_groupby$Taxa1,"_",Fisher_permute_summary_groupby$Taxa2) #Create vector to reorder scatter plot color legend
Fisher_permute_summary_groupby$Taxa1 <- gsub("_"," ",Fisher_permute_summary_groupby$Taxa1) #Remove underscore for aesthetic
Fisher_permute_summary_groupby$Taxa2 <- gsub("_"," ",Fisher_permute_summary_groupby$Taxa2) #Remove underscore for aesthetic
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test/bubble_plot_fisher_test.pdf", height =  5, width = 10)
plot <- ggplot(Fisher_permute_summary_groupby, aes(x=Taxa1, y=Taxa2, color=Mean_OddRatio, size=Permutation)) +
  geom_point() +
  scale_size(range = c(0, 15)) +
  theme_light() +
  scale_fill_gradient(low = "royalblue4", high = "red") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=15, hjust = 1, angle = 20, ))
print(plot)
dev.off()



#Scatter plot, Odds ratio and Pvalue
colours <- c("grey75","dodgerblue","coral1","chartreuse4","mAgenta","grey14","seagreen1","darkorange","purple","forestgreen","firebrick","aquamarine","lightpink")
Fisher_permute_summary$pairwise <- factor(Fisher_permute_summary$pairwise, level=c("NS",order_legend_pairwise))
bubble_plot <- ggplot(Fisher_permute_summary,aes(x=Odds_ratio, y=Pvalue_Adjust, color=pairwise, shape=pairwise)) +
  geom_point(size=5) +
  scale_shape_manual(values=c(1:20))+
  xlim(0,20) +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values=colours) +
  theme_light() +
  guides(col=guide_legend(ncol=2)) +
  geom_vline(xintercept=1, linetype='dashed') +
  geom_text(aes(x=1, label="\nOdds Ratio = 1", y=0.25), colour="red", angle=90, size=10) +
  geom_hline(yintercept=0.05, linetype='dashed') +
  geom_text(aes(y=0.07, label="\nAdj.Pvalue = 0.05", x=2.5), colour="red", angle=0, size=10) +
  theme(axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=20),
        legend.text = element_text(size=14),
        legend.position = 'bottom')


combined <- ggarrange(
  bubble_plot,                # First row with line plot
  # Second row with box and histogram
  ggarrange(plot,bar_plot,hist_plot, widths = c(0.6,0.2,0.3), ncol=3),
  nrow = 2, 
  labels = "A",       # Label of the line plot
  heights = c(1,0.5),
  widths = c(1,1)
  )
ggsave(combined,filename = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test/fisher_test_multipannel.pdf", width = 30, height = 20)

write.table(Fisher_permute_summary,"/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/bubble_plot_fisher_test.txt", sep="\t", row.names = FALSE)
write.table(Fisher_test_summary,"/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/fisher_test/fisher_test.txt", sep ="\t" )
Fisher_permute_summary <- fread("/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_btw_symbiont/bubble_plot_fisher_test.txt")


#Understand Infinite value
x <- subset(Fisher_permute_summary, Odds_ratio == "Inf")
x2 <- subset(x, Pvalue_Adjust<0.05 )
taxa1_taxa2_cont_table <- table(asv_tax_spl_rnk8[,"Parentodinium_sp._ASV132"],asv_tax_spl_rnk8[,"Cycloposthium_bipalmatum"])
fisher.test(taxa1_taxa2_cont_table)
dimnames(taxa1_taxa2_cont_table) <- list(Parentodinium_sp._ASV132 = c("Absence_Taxa1", "Presence_Taxa1"),
                                         Cycloposthium_bipalmatum = c("Absence_Taxa2", "Presence_Taxa2"))
table <- as.data.frame(taxa1_taxa2_cont_table)
colours <- c("orangered1","blue")
plot <- ggplot(table, aes(x=Cycloposthium_bipalmatum, y=Freq, fill=Parentodinium_sp._ASV132, color = colours )) +
  geom_bar(stat="identity", width=1, color="black", linetype = "solid") +  #height of the column equal the value, stacked
  scale_fill_manual(values=colours) +
  theme_bw()+
  labs(x ="Cycloposthium_bipalmatum" ) +
  #ggtitle(paste0("Total:",total_sample,", Sens:",sensitivity, ", Spec:",specificity,", Accu:",accuracy,", Pval(Fisher):",fisher_Pval))+
  theme(legend.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        plot.title=element_text(size=8,face="bold"))

#
x <- subset(Fisher_permute_summary, pairwise == "Blastocystis_Unknown_ASV45_Enteromonad_Unknown_ASV73")
x2 <- subset(x, Pvalue_Adjust>0.2 )

ggplot(x, aes(x=Pvalue_Adjust, y=Spl_taxa1)) +
  geom_point() +
  geom_vline(xintercept=0.05, linetype='dashed')

test <- cor.test(x$Pvalue_Adjust,x$Spl_taxa1)
summary(test)

taxa1_taxa2_cont_table <- table(asv_tax_spl_rnk8[,"Blastocystis_Unknown_ASV45"],asv_tax_spl_rnk8[,"Enteromonad_Unknown_ASV73"])
dimnames(taxa1_taxa2_cont_table) <- list(Taxa1 = c("Absence_Taxa1", "Presence_Taxa1"),
                                         Taxa2 = c("Absence_Taxa2", "Presence_Taxa2"))

colours <- c("orangered1","blue")
plot <- ggplot(taxa1_taxa2_cont_table, aes(x=Taxa1, y=Freq, fill=Taxa2, color = colours )) +
  geom_bar(stat="identity", width=1, color="black", linetype = "solid") +  #height of the column equal the value, stacked
  scale_fill_manual(values=colours) +
  theme_bw()+
  labs(y ="Taxa2" , x ="Taxa1" ) +
  #ggtitle(paste0("Total:",total_sample,", Sens:",sensitivity, ", Spec:",specificity,", Accu:",accuracy,", Pval(Fisher):",fisher_Pval))+
  theme(legend.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        plot.title=element_text(size=8,face="bold"))



########################################################

### Correlation between abundance
########################################################


########################################################

### Association between Presence and alpha diversity
########################################################

########################################################

### Association between Presence and beta diversity
########################################################

########################################################


################################################################################################################
###                  RELATIONSHIP BETWEEN INTESINAL EUKARYOTES and ENVIRONMENT AND HOST FACTORS              ###
################################################################################################################

########################################################
###                Beta diversity                    ###
########################################################

sort(sample_sums(project_data))
project_data_rarefied <- rarefy_even_depth(project_data, sample.size = 14418);project_data_rarefied

#Season
##project_data_rarefied <- project_data_rarefied %>%
#  subset_samples(Site != "Lokoue")

#Combine factor
sample_df <- data.frame(sample_data(project_data_rarefied))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Site*Age*Yaws*Social_Status*Season, data=sample_df, method="bray")
#Need to check if I used sequential or marhinal (Type 1 and type 2)
res.adonis.rarefied.bray


#Site
sample_df <- data.frame(sample_data(project_data_rarefied))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Site, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Site, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
## Create Colour Palettes ###
n <- length(unique(NMDS$gorilla_name))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-16]
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Site
sample_df <- data.frame(sample_data(project_data_rarefied))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Site, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Site, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Site_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Site)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Season
project_data_Season <- project_data_rarefied %>%
  subset_samples(Season != "NA")

sample_df <- data.frame(sample_data(project_data_Season))
project_bray.rarefied.bray <- phyloseq::distance(project_data_Season, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Site*Age*Yaws*Social_Status*Season, data=sample_df, method="bray")


#Season
sample_df <- data.frame(sample_data(project_data_Season))
project_bray.rarefied.bray <- phyloseq::distance(project_data_Season, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Season, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_Season, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_Season))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Season, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Season_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Season
sample_df <- data.frame(sample_data(project_data_Season))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_Season, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ Season, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
NMDS.jaccard <- ordinate(physeq = project_data_Season, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_Season))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Season_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = Season)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()


#control <- ggarrange(NMDS_Season_jaccard,NMDS_Season_bray, labels = c("A","B"))
#control2 <- ggarrange(control,alpha_Site,nrow = 2, labels = c("","C"))


#Remove data from Moba and Maya Sites (no information on dates)
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Moba") %>%
  subset_samples(Site !="Maya")

#Age
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Age, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Age, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Age_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Age)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","blue2","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Age
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Age, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Age, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Age_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Age)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  scale_color_manual(values=c("darkgreen","blue2","orangered3")) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()


#Social_Status
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Social_Status, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Social_Status, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Social_Status_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Social_Status)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Social_Status
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Social_Status, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Social_Status, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Social_Status_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Social_Status)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("orangered3","darkgreen")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()


#Yaws
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Yaws, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Yaws)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Yaws
sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
project_bray.rarefied.bray <- phyloseq::distance(project_data_rarefied_host_factor, method = "jaccard")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ Yaws, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
NMDS.bray <- ordinate(physeq = project_data_rarefied_host_factor, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_rarefied_host_factor))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$Yaws, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_Yaws_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Yaws)) + 
  geom_point(size=4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  ggtitle(as.character(title)) + 
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","orangered3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#Remove data from Moba and Maya Sites (no information on dates)
project_data_exp <- project_data %>% subset_samples(temporal == "yes") 

#gorilla name
sample_df <- data.frame(sample_data(project_data_exp))
project_bray.rarefied.bray <- phyloseq::distance(project_data_exp, method = "bray")
res.adonis.rarefied.bray <- adonis(project_bray.rarefied.bray ~ gorilla_name, data=sample_df, method="bray")
r2_treatment <- round(res.adonis.rarefied.bray$aov.tab$R2[1],3)
pval_treatment <- round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
NMDS.bray <- ordinate(physeq = project_data_exp, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data_exp))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$gorilla_name, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
## Create Colour Palettes ###
n <- length(unique(NMDS$gorilla_name))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-16]
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_gorilla_name_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = gorilla_name)) + 
  geom_point(size=4) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=col_vector) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()

#gorilla_name
sample_df <- data.frame(sample_data(project_data_exp))
project_jaccard.rarefied.bray <- phyloseq::distance(project_data_exp, method = "jaccard")
res.adonis.rarefied.jacard <- adonis(project_jaccard.rarefied.bray ~ gorilla_name, data=sample_df, method="jaccard")
r2_treatment <- round(res.adonis.rarefied.jacard$aov.tab$R2[1],3)
title <- paste0("R2=",r2_treatment," ","P=",pval_treatment)
pval_treatment <- round(res.adonis.rarefied.jacard$aov.tab$`Pr(>F)`[1],3)
NMDS.jaccard <- ordinate(physeq = project_data_exp, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(project_data_exp))
jaccard <- as.data.frame(NMDS.jaccard$points)
row.names(jaccard) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.jaccard1 <- jaccard$MDS1
NMDS$NMDS.jaccard2 <- jaccard$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$gorilla_name, NMDS$Yaws),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
NMDS_gorilla_name_jaccard <- ggplot(NMDS, aes(x=NMDS.jaccard1, y=NMDS.jaccard2, color = gorilla_name)) + 
  geom_point(size=4) +
  ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=col_vector) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_classic()



#combine plot
NMDS_individuals <- ggarrange(NMDS_Site_bray,NMDS_Site_jaccard,NMDS_Season_bray,NMDS_Season_jaccard,NMDS_Age_bray,NMDS_Age_jaccard,NMDS_Social_Status_bray,NMDS_Social_Status_jaccard,NMDS_Yaws_bray,NMDS_Yaws_jaccard,NMDS_gorilla_name_bray,NMDS_gorilla_name_jaccard,
                              labels = c("A","B","C","D","E","F","G","H","I","J","K","L"), nrow = 6, ncol = 2)
ggsave(NMDS_individuals,filename = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_with_gorilla_features/beta/interaction_with_features_NMDS_individual.pdf", width = 15, height = 25)

########################################################

########################################################
###                Alpha diversity                   ###
########################################################

#Site
my_comparisons <- list(c("Romani", "Lokoue"),c("Romani", "Maya"),c("Romani", "Moba"),
                       c("Lokoue", "Maya"),c("Lokoue", "Moba"),
                       c("Maya", "Moba"))

alpha_Site <- plot_richness(project_data_rarefied, x = "Site", color = "Site", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot()+
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  stat_compare_means(comparisons=my_comparisons, aes(label = paste0("p =", ..p.format..)))

#Season
project_data_rarefied_Season <- project_data_rarefied %>%
  subset_samples(Season != "NA")

my_comparisons <- list(c("Short-Dry", "Short-Rain"),c("Short-Dry", "Long-Dry"),c("Short-Dry", "Long-Rain"),
                       c("Short-Rain", "Long-Rain"),c("Short-Rain", "Long-Dry"),
                       c("Long-Rain", "Long-Dry"))

alpha_Season <- plot_richness(project_data_rarefied_Season, x = "Season", color = "Season", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot()+
  scale_color_manual(values=c("orangered3","darkgreen","blue2","pink")) +
  stat_compare_means(comparisons=my_comparisons, aes(label = paste0("p =", ..p.format..)))


#Age
project_data_rarefied_host_factor <- project_data_rarefied %>%
  subset_samples(Site !="Moba") %>%
  subset_samples(Site !="Maya")

my_comparisons <- list(c("BB", "SB"),c("BB", "sub"),c("SB", "sub"))

alpha_Age <- plot_richness(project_data_rarefied_host_factor, x = "Age", color = "Age", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot()+
  scale_color_manual(values=c("darkgreen","blue2","orangered3")) +
  stat_compare_means(comparisons=my_comparisons, aes(label = paste0("p =", ..p.format..)))

#Social status  
alpha_Social_Status <- plot_richness(project_data_rarefied_host_factor, x = "Social_Status", color = "Social_Status", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot()+
  scale_color_manual(values=c("orangered3","darkgreen")) +
  stat_compare_means(aes(label = paste0("p =", ..p.format..)))

#Yaws
alpha_Yaws <- plot_richness(project_data_rarefied_host_factor, x = "Yaws", color = "Yaws", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot()+
  scale_color_manual(values=c("darkgreen","orangered3")) +
  stat_compare_means(aes(label = paste0("p =", ..p.format..)))


#combine plot
alpha <- ggarrange(alpha_Site,alpha_Age,alpha_Social_Status,alpha_Yaws,alpha_Season, labels = c("A","B","C","D","E"), nrow = 3, ncol = 2)
ggsave(alpha,filename = "/Users/vincebilly/Desktop/vince/phd/statistic/gorilla/figure/interaction_with_gorilla_features/alpha/alpha_individual_allMetrics.pdf", width = 30, height = 20)

#Yaws per Site
alpha_Yaws_Site <- plot_richness(project_data, x = "Site", color = "Yaws", measures = c("Shannon")) +
  geom_boxplot(position=position_dodge2(padding=0.5, width = 0.5)) + 
  stat_compare_means(aes(label = paste0("p =", ..p.format..))) +
  scale_color_manual(values=c("darkgreen","orangered3"))
  


sample_df <- data.frame(sample_data(project_data_rarefied_host_factor))
otu <- otu_table(t(project_data_rarefied_host_factor))
shannon <- diversity(otu, index = "shannon", MARGIN = 1, base = exp(1))
otu[otu >0] <- 1
number_ASVs <- rowSums(otu)
meta_shannon <- cbind(sample_df,shannon,number_ASVs)
meta_shannon <- meta_shannon %>% subset(Age !="sub")

model <- lm(shannon ~ Yaws * Social_Status * Age * Site * Season, data = meta_shannon)
summary(model)

my_comparisons <- list(c("no", "yes"))


alpha_shannon_treatment_without <- ggplot(meta_shannon, aes(x = Site, y =shannon, color=Yaws)) + 
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8))+
  xlab("Site") +
  ylab("Shannon") +
  theme_bw() +
  stat_compare_means(aes(label = paste0("p =", ..p.format..))) +
  labs(title="Shannon") +
  scale_color_manual(values=c("darkgreen","orangered3")) +

lm <- lm(shannon ~ Yaws + Site + Social_Status + Age,  data = meta_shannon)
anova(lm)
TukeyHSD(aov(lm))

lm <- lm(shannon ~ Yaws + Site + Social_Status,  data = meta_shannon)
anova(lm)
TukeyHSD(aov(lm))

lm <- lm(shannon ~ Yaws * Site * Social_Status * Age,  data = meta_shannon)
anova(lm)
TukeyHSD(aov(lm))

########################################################

