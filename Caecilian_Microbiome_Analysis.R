#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(tidyverse)
library(microbiome)
library(ggplot2)
library(vegan)

setwd("C:/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/")

##Note: when running this code the metadata import kept cutting out the header and making the first row of data the header, so followed the advice at the end of this tutorial and added a #q2:types row
#https://forum.qiime2.org/t/headers-removed-when-importing-phyloseq-object-for-decontam/8075/8
CB <- qza_to_phyloseq(features = "Outputs/table-dada2-caecilian-paired-reads.qza",tree = "Outputs/rooted-tree-caecilian-paired.qza",taxonomy = "Outputs/taxonomy-caecilian-paired.qza",metadata = "Inputs/Caecilian_Micro_Metadata_2020.2.txt")
CB
#Check metadata
sample_data(CB)

################ Basic SUMMARY of data ################
## total number of sequences per sample, sorted in order
sort(sample_sums(CB))
#CB stands for caecilian bacteria

## lowest sample has 11,458 sequences and highest has 38,235
38235/11458
## Approx 3 fold difference, which isn't too bad
caecilians_final <- subset_samples(CB, Category == "caecilian")
sample_data(caecilians_final)
## total number of sequences per sample, sorted in order
sort(sample_sums(caecilians_final))

## Based on McMurdie & Holmes 2014 we will not rarefy because of the bias in introduces and since our library sizes differ by three fold then according to Weiss et al. 2017 rarefying would not change the false discovery rate. It only does so when library sizes are greater than ten fold differences

##Bactera, archaea and unassigned  
get_taxa_unique(caecilians_final, "Kingdom")

#How many archaeal phyla?
caecilian_Archae = subset_taxa(caecilians_final, Kingdom == "d__Archaea")
get_taxa_unique(caecilian_Archaea, "Phylum")
#"Nanoarchaeota"            "Methanobacteriota_A_1229" "Thermoproteota"           "Halobacteriota"           "Thermoplasmatota" 
#5

#How many bacterial phyla?
caecilian_Bacteria = subset_taxa(caecilians_final, Kingdom == "d__Bacteria")
get_taxa_unique(caecilian_Bacteria, "Phylum")
# [1] NA                          "Proteobacteria"            "Patescibacteria"           "Chloroflexota"             "Deinococcota"             
#[6] "Dependentiae"              "Armatimonadota"            "Fusobacteriota"            "Verrucomicrobiota"         "Spirochaetota"            
#[11] "Cyanobacteria"             "Bacteroidota"              "Firmicutes_D"              "Firmicutes_C"              "Firmicutes_A"             
#[16] "Firmicutes_B_370539"       "Actinobacteriota"          "Firmicutes_B_370543"       "Dormibacterota"            "Nitrospirota_A_437815"    
#[21] "Acidobacteriota"           "Bdellovibrionota_E"        "Desulfobacterota_B"        "Desulfobacterota_G_459544" "Desulfobacterota_I"       
#[26] "Myxococcota_A_473307"      "Desulfobacterota_G_459546" "Bdellovibrionota_C"        "Desulfobacterota_C"        "UBA10199"                 
#[31] "Planctomycetota"           "FCPU426"                   "Eisenbacteria"             "Firmicutes_G"              "Gemmatimonadota"          
#[36] "Firmicutes_B_370541"       "Elusimicrobiota"           "Sumerlaeota"               "Eremiobacterota"           "Methylomirabilota"        
#[41] "Chlamydiota"               "Myxococcota_A_437813"      "Campylobacterota"          "Fibrobacterota"
#44; 43 minus NA

library(dplyr)
get_taxa_unique(caecilians_final, "Phylum")
tax_table(caecilians_final)[,"Phylum"] %>% unique %>% na.exclude %>% length
#Another way to see 48 total phyla in caecilian samples

#### 352,739 total sequences
sum(taxa_sums(caecilians_final))

#### 3064 ASVs
caecilians_final

## total sequence count
sum(taxa_sums(caecilians_final))

##Top 15 phyla and genera
top15ph <- sort(tapply(taxa_sums(caecilians_final), tax_table(caecilians_final)[, "Phylum"],sum), TRUE)[1:15]
top15g <- sort(tapply(taxa_sums(caecilians_final), tax_table(caecilians_final)[, "Genus"], sum), TRUE)[1:15]
top15ph
top15g

##What percentage of total sequences are these top phyla?
##Proteobacteria
Proteo = subset_taxa(caecilians_final, Phylum == "Proteobacteria")
sum(taxa_sums(Proteo))
204433/352739
##Proteobacteria makes up 58% of all sequences

##Firmicutes
Firmicutes = subset_taxa(caecilians_final, Phylum %in% c("Firmicutes_D","Firmicutes_B_370539","Firmicutes_C","Firmicutes_A","Firmicutes_B_370543","Firmicutes_G","Firmicutes_B_370541"))
sum(taxa_sums(Firmicutes))
61269/352739
##Firmicutes makes up 17% of all sequences

##Bacteroidota
Bacteroidota = subset_taxa(caecilians_final, Phylum == "Bacteroidota")
sum(taxa_sums(Bacteroidota))
47713/352739
##Bacteroidota makes up 14% of all sequences

##Actinobacteriota
Actino = subset_taxa(caecilians_final, Phylum == "Actinobacteriota")
sum(taxa_sums(Actino))
20011/352739
##Actinobacteria makes up 6% of all sequences

##Acidobacteriota
Acid = subset_taxa(caecilians_final, Phylum == "Acidobacteriota")
sum(taxa_sums(Acid))
2129/352739
##Acidobacteriota makes up 0.6% of all sequences

####Soil samples only####
SB <- subset_samples(CB, Category == "soil")
sample_data(SB)
## total number of sequences per sample, sorted in order
sort(sample_sums(SB))

## lowest sample has 11,434 sequences and highest has 13,355
13355/11434
## Approx 1 fold difference

##Bactera, archaea and unassigned  
get_taxa_unique(SB, "Kingdom")

#How many archaeal phyla?
SB_Archae = subset_taxa(SB, Kingdom == "d__Archaea")
get_taxa_unique(SB_Archae, "Phylum")
#"Nanoarchaeota"            "Methanobacteriota_A_1229" "Thermoproteota"           "Halobacteriota"           "Thermoplasmatota" 
#5

#How many bacterial phyla?
SB_Bacteria = subset_taxa(SB, Kingdom == "d__Bacteria")
get_taxa_unique(SB_Bacteria, "Phylum")
# [1] NA                          "Proteobacteria"            "Patescibacteria"           "Chloroflexota"             "Deinococcota"             
#[6] "Dependentiae"              "Armatimonadota"            "Fusobacteriota"            "Verrucomicrobiota"         "Spirochaetota"            
#[11] "Cyanobacteria"             "Bacteroidota"              "Firmicutes_D"              "Firmicutes_C"              "Firmicutes_A"             
#[16] "Firmicutes_B_370539"       "Actinobacteriota"          "Firmicutes_B_370543"       "Dormibacterota"            "Nitrospirota_A_437815"    
#[21] "Acidobacteriota"           "Bdellovibrionota_E"        "Desulfobacterota_B"        "Desulfobacterota_G_459544" "Desulfobacterota_I"       
#[26] "Myxococcota_A_473307"      "Desulfobacterota_G_459546" "Bdellovibrionota_C"        "Desulfobacterota_C"        "UBA10199"                 
#[31] "Planctomycetota"           "FCPU426"                   "Eisenbacteria"             "Firmicutes_G"              "Gemmatimonadota"          
#[36] "Firmicutes_B_370541"       "Elusimicrobiota"           "Sumerlaeota"               "Eremiobacterota"           "Methylomirabilota"        
#[41] "Chlamydiota"               "Myxococcota_A_437813"      "Campylobacterota"          "Fibrobacterota"
#44; 43 minus NA

#### 37,644 total sequences
sum(taxa_sums(SB))
#### 3064 ASVs
SB
##Top 15 phyla and genera
top15phSB <- sort(tapply(taxa_sums(SB), tax_table(SB)[, "Phylum"],sum), TRUE)[1:15]
top15gSB <- sort(tapply(taxa_sums(SB), tax_table(SB)[, "Genus"], sum), TRUE)[1:15]
top15phSB
top15gSB

##What percentage of total sequences are these top phyla?
##Proteobacteria
ProteoSB = subset_taxa(SB, Phylum == "Proteobacteria")
sum(taxa_sums(ProteoSB))
18479/37644
##Proteobacteria makes up 49% of all sequences

##Acidobacteriota
AcidSB = subset_taxa(SB, Phylum == "Acidobacteriota")
sum(taxa_sums(AcidSB))
4654/37644
##Acidobacteriota makes up 12% of all sequences

##Bacteroidota
BacteroidotaSB = subset_taxa(SB, Phylum == "Bacteroidota")
sum(taxa_sums(BacteroidotaSB))
4314/37644
##Bacteroidota makes up 11% of all sequences

##Actinobacteriota
ActinoSB = subset_taxa(SB, Phylum == "Actinobacteriota")
sum(taxa_sums(ActinoSB))
3245/37644
##Actinobacteria makes up 9% of all sequences

##Firmicutes
FirmicutesSB = subset_taxa(SB, Phylum %in% c("Firmicutes_D","Firmicutes_C","Firmicutes_A"))
sum(taxa_sums(FirmicutesSB))
1487/37644
##Firmicutes makes up 4% of all sequences




#### Alpha Diversity Comparison
CB_alpha <- microbiome::alpha(CB, index = "all")
head(CB_alpha)
CB_meta <- sample_data(CB)
CB_meta$Observed <- CB_alpha$observed
CB_meta$Shannon <- CB_alpha$diversity_shannon

CB_filtered <- filter_taxa(CB, function(x) sum(x > 2) > 2, TRUE)
CB
#No change

library(phyloseq)
library(vegan)
library(RVAideMemoire)

#subset to only look at caecilian samples
CB_caecilian <- subset_samples(CB, Category == "caecilian")

##Are there significant differences by body location sampled?
#Jaccard
Jacc_CB <- phyloseq::distance(CB_caecilian,"jaccard")
df_CB <- as(sample_data(CB_caecilian), "data.frame")
df_CB <- na.omit(df_CB)
perm <- how(nperm = 999)
adonis2(Jacc_CB ~ Body_location, data = df_CB, permutations = perm)
#Not significantly different
#Bray-Curtis
bray_CB <- phyloseq::distance(CB_caecilian,"bray")
adonis2(bray_CB ~ Body_location, data = df_CB, permutations = perm)
#Not for bray-curtis either

##What about differences in the cutaneous microbiome between specimens?
#Jaccard
Jacc_CB <- phyloseq::distance(CB_caecilian,"jaccard")
adonis2(Jacc_CB ~ Caecilian_number, data = df_CB, permutations = perm)
#Significantly different
pairwise.perm.manova(Jacc_CB, df_CB$Caecilian_number, nperm = 999)
#Bray-Curtis
bray_CB <- phyloseq::distance(CB_caecilian,"bray")
adonis2(bray_CB ~ Caecilian_number, data = df_CB, permutations = perm)
#Significantly different
pairwise.perm.manova(bray_CB, df_CB$Caecilian_number, nperm = 999)

##Are caecilian samples significantly different from soil samples?
all_B <- subset_samples(CB, Category %in% c("caecilian","soil"))
df_all_B <- as(sample_data(all_B), "data.frame")
df_all_B <- na.omit(df_all_B)
perm <- how(nperm = 999)
#Jaccard
Jacc_all_B <- phyloseq::distance(all_B,"jaccard")
adonis2(Jacc_all_B ~ Caecilian_number, data = df_all_B, permutations = perm)
#Significantly different
pairwise.perm.manova(Jacc_CB, df_CB$Caecilian_number, nperm = 999)
#Bray-Curtis
bray_CB <- phyloseq::distance(CB_caecilian,"bray")
adonis2(bray_CB ~ Caecilian_number, data = df_CB, permutations = perm)
#Significantly different
pairwise.perm.manova(bray_CB, df_CB$Caecilian_number, nperm = 999)

bray_CB <- phyloseq::distance(CB_caecilian,"bray")
bray_ord <- ordinate(CB_caecilian,method = "PCoA","bray")
p <- plot_ordination(CB_caecilian, bray_ord, color = "Caecilian_number",shape = "Body_location") +
  geom_point(size = 5) +
  scale_shape_manual(values=c(21,22,23,24,25), name="Body location") +
  scale_color_manual(values = c("darkorchid4","goldenrod","darkgrey"), name = "Specimen number") +
  theme_q2r() +
  theme(text = element_text(size=20))

#Fixed scale_shape_manual issue using this tutorial https://github.com/joey711/phyloseq/issues/250
p
p$layers
p$layers <- p$layers[-1]
p
library(Cairo)
ggsave(filename = "/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/Figures/Bray_curtis_caecilians_no_env.eps",
       plot = print(p),
       device = "eps",
       dpi = 300,
       width = 15,
       height = 12,
       unit = "cm")
ggsave(filename = "/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/Figures/Bray_curtis_caecilians_no_env.tiff",
       plot = print(p),
       device = "tiff",
       dpi = 300,
       width = 15,
       height = 12,
       unit = "cm")

bray_CB <- phyloseq::distance(all_B,"bray")
bray_ord <- ordinate(all_B,method = "PCoA","bray")
p <- plot_ordination(all_B, bray_ord, color = "Caecilian_number",shape = "Body_location") +
  geom_point(size = 5) +
  scale_shape_manual(values=c(21,22,23,24,25), name="Body location") +
  scale_color_manual(values = c("darkorchid4","goldenrod","darkgrey"), name = "Specimen number") +
  theme_q2r() +
  theme(text = element_text(size=20))

#Fixed scale_shape_manual issue using this tutorial https://github.com/joey711/phyloseq/issues/250
p
p$layers
p$layers <- p$layers[-1]
p
library(Cairo)
ggsave(filename = "/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/Figures/Bray_curtis_caecilians_no_env.eps",
       plot = print(p),
       device = "eps",
       dpi = 300,
       width = 15,
       height = 12,
       unit = "cm")

#####################Stacked bar plot
#install.packages("devtools")
#devtools::install_github('bbc/bbplot')
library(bbplot)

#####***Make a stacked barplot with genus and phyla data by specimen and body location

#https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/index.html
###Stacked bar plot by Phylum
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
caecilians_phylum <- CB_caecilian %>% #Create new object
  tax_glom(taxrank = "Phylum") %>% #Condense by Phylum      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts
  arrange(Phylum) #Order Phyla


ggplot(caecilians_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

caecilians_genus <- CB_caecilian %>% #Create new object
  tax_glom(taxrank = "Genus") %>% #Condense by Genus      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts
  arrange(Genus) #Order Genera

g<- ggplot(caecilians_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

ggsave(filename = "/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/Figures/Stacked_barplot_genera.eps",
       plot = print(g),
       device = "eps",
       dpi = 300,
       width = 15,
       height = 12,
       unit = "cm")

#Run again with x axis sample names to be sure what order the samples are in
ggplot(caecilians_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))


##Stacked barplot with soil samples included
all_B <- subset_samples(CB, Category %in% c("caecilian","soil"))
sample_names(all_B)
all_genus <- all_B %>% #Create new object
  tax_glom(taxrank = "Genus") %>% #Condense by Genus      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts
  arrange(Genus) #Order Genera

ggplot(all_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category,space = "free_x", scales = "free_x")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

if(!"devtools" %in% installed.packages()){
  install.packages("devtools")
}
devtools::install_github("gmteunisse/fantaxtic")
devtools::install_github("gmteunisse/ggnested")
library(fantaxtic)
library(ggnested)
top_nested <- nested_top_taxa(all_B,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Genus",
                              n_top_taxa = 6, 
                              n_nested_taxa = 6)
nested_plot <- plot_nested_bar(ps_obj = top_nested$ps_obj,
                               top_level = "Phylum",
                               nested_level = "Genus")
nested_plot
ggsave(filename = "/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/Figures/Stacked_barplot_nested.eps",
       plot = print(nested_plot),
       device = "eps",
       dpi = 300,
       width = 50,
       height = 25,
       unit = "cm")


jacc_CB <- phyloseq::distance(CB, "jaccard")
jacc_ord <- ordinate(CB, method = "PCoA", "jaccard")
plot_ordination(CB, jacc_ord, color = "X2") +
  geom_point(size = 5) +
  theme_bw()


jacc_CB <- phyloseq::distance(all_B, "jaccard")
jacc_ord <- ordinate(all_B, method = "PCoA", "jaccard")
plot_ordination(all_B, jacc_ord, color = "X2") +
  geom_point(size = 5) +
  theme_bw()



#Set permutations
perm <- how(nperm = 999)
#Set blocks for stratification
setBlocks(perm) <- with(df_SB, Site)
#Run permutaion test
adonis2(bray_SB ~ Day, data = df_SB, permutations = perm)
