##Frog Microbiome Analysis
##Foam nest breeding frogs (Genus:Polypedates) from Brunei, Borneo
##Analysis performed April 2020

##Intall required packages.
##Needed up update BioConductor and install phyloseq specially.https://support.bioconductor.org/p/120724/
#if(!requireNamespace("BiocManager")){
#  install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")
#install.packages("ape")
#install.packages("vegan")
#install.packages("ggplot2")
#install.packages("tidyverse")
#install.packages("RVAideMemoire")
#install.packages("picante")
#install.packages("agricolae")

##Load libraries.
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)
library(dplyr)
library(RVAideMemoire)
library(picante)
library(plyr)
library(car)
library(agricolae)

##Set working directory.
setwd("/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Github_repo/")

########## INITIAL ANALYSES TO SHOW THAT RUN EFFECT WAS MINIMAL IN TERMS OF BETA-DIVERSITY############
### BIOM FILE ###
biom_sep_runs <- import_biom("biom_sep_runs/table-with-taxonomy.biom", parseFunction = parse_taxonomy_greengenes)
sample_names(biom_sep_runs)

### MAPPING FILE ####

map_sep_runs <- import_qiime_sample_data("map_sep_runs/Frog_Micro_Metadata_2020.2 - Sheet1 (2).tsv")
sample_names(map_sep_runs)

#####  TREE FILE  #####
tree = read.tree("tree/tree.nwk")

fnm <- merge_phyloseq(biom_sep_runs, tree, map_sep_runs)
fnm

### look at run effects

jaccard_run <- distance(fnm,"jaccard", binary = T)

jacc.ord_run <-ordinate(fnm, method = "PCoA", jaccard_run)

p_run <- plot_ordination(fnm, jacc.ord_run, color = "MiSeqRun")

p_run + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) 
#geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

df_fnm <- as(sample_data(fnm), "data.frame")

groups_run <- df_fnm[["MiSeqRun"]]
jacc_disprun <-betadisper(jaccard_run, groups_run, type=c("median"))
plot(jacc_disprun)
anova(jacc_disprun)
boxplot(jacc_disprun)

#Based on these results, then went back into QIIME2 and grouped samples by True_SampleID based on this code https://docs.qiime2.org/2020.2/plugins/available/feature-table/group/ (see also GitHub ReadMe Doc)

########################################Samples summed by True_SampleID##############################
#See QIIME2-2020.2 code for how these were grouped.
### BIOM FILE, RUNS MERGED #####
biom_grouped <- import_biom("biom_grouped/grouped-table-with-taxonomy.biom", parseFunction = parse_taxonomy_greengenes)
sample_names(biom_grouped)

#### MAPPING FILE ####
### We added alpha-diversity to Frog_Metadata_grouped.txt file, see at end of this script
map_grouped <- import_qiime_sample_data("map_grouped/Frog_Micro_Metadata_2020.2 - Sheet3.tsv")
sample_names(map_grouped)

#####  TREE FILE  #####
tree = read.tree("tree/tree.nwk")

fnm_g <- merge_phyloseq(biom_grouped,tree, map_grouped)

###Need to add alpha diversity metrics to mapping file###
## extracting your otu_table table from the phyloseq object
species_site <-as(otu_table(fnm_g), "matrix")

## transpose it to get in correct format
site_species <- t(species_site)

# prune the tree to only include tips from your otu_table (might be extra)
prunedTree <- prune.sample(site_species,tree)

#Changed include.root = F to include.root = T because I kept getting a warning code
# info on what this change might mean https://www.rdocumentation.org/packages/picante/versions/1.8.1/topics/pd
#Tutorial with include.root = T http://picante.r-forge.r-project.org/picante-intro.pdf 
PD <- pd(site_species, prunedTree, include.root = T)

#need to have both alpha and df having the same column info
PD$X.SampleID <- row.names(PD)

#Create dataframe
df_fnm_g <- as(sample_data(fnm_g), "data.frame")

#now merge
df_fnm_final <- merge(df_fnm_g, PD, by = "X.SampleID")

#Adding alpha-diversity to meta file 
# NOW, we will just read in this file from now on
write.table(df_fnm_final, "Frog_Metadata_Final_2020.txt", sep ="\t", row.names = T)

##Re-run with new mapping file to get final phyloseq object
### BIOM FILE, RUNS MERGED #####
biom_grouped <- import_biom("biom_grouped/grouped-table-with-taxonomy.biom", parseFunction = parse_taxonomy_greengenes)

#### MAPPING FILE ####
### We added alpha-diversity to Frog_Metadata_grouped.txt file
map_grouped_w_PD <- import_qiime_sample_data("Frog_Metadata_Final_2020.txt")

#####  TREE FILE  #####
tree = read.tree("tree/tree.nwk")

fnm_g_w_PD <- merge_phyloseq(biom_grouped,tree, map_grouped_w_PD)

############## FILTER OTUS ###################

## keep taxa that have more than two sequences in more than 1 sample 
# It has at least two sequences in at least 2 samples 

## from here: https://joey711.github.io/phyloseq/preprocess.html
## search filter_taxa
fnm_r_final = filter_taxa(fnm_g_w_PD, function(x) sum(x > 2) > 2, TRUE)


################ Basic SUMMARY of data ################

## total number of sequences per sample, sorted in order
sort(sample_sums(fnm_r_final))

## lowest sample has 5,213 sequences and highest has 28,676
## Approx 5 fold difference, which isn't too bad
28676/5213

## Based on McMurdie & Holmes 2014 we will not rarefy because of the bias in introduces and
## since our library sizes differ by five fold then according to Weiss et al. 2017 rarefying would 
## not change the false discovery rate. It only does so when library sizes are greater
## than ten fold differences

## 
get_taxa_unique(fnm_r_final, "Kingdom")

fnm_Archae = subset_taxa(fnm_r_final, Kingdom == "Archaea")

## these are your Archael phyla
get_taxa_unique(fnm_Archae, "Phylum")

## 28 bacterial and archael phyla (one NA, and 2 Archaeal phylum)
get_taxa_unique(fnm_r_final, "Phylum")
tax_table(fnm_r_final)[,"Phylum"] %>% unique %>% na.exclude %>% length

#### 1,377,542 total sequences
sum(taxa_sums(fnm_r_final))

#### 2850 sOTUs
fnm_r_final
#This says there are 0 taxa that have a sum of 0.
sum(taxa_sums(fnm_r_final) == 0)

## % Proteobacteria
Proteo = subset_taxa(fnm_r_final, Phylum == "Proteobacteria")

## total sequence count
sum(taxa_sums(fnm_r_final))

## Proteobacteria
sum(taxa_sums(Proteo))

## Proteobacteria makes up 58% of all sequences
803331/1377542

## Firmicutes
Firmicutes = subset_taxa(fnm_r_final, Phylum == "Firmicutes")
sum(taxa_sums(Firmicutes))
214289/1377542
##Firmicutes makes up 15% of all sequences

##Bacteroidetes
Bacteroidetes = subset_taxa(fnm_r_final, Phylum == "Bacteroidetes")
sum(taxa_sums(Bacteroidetes))
181077/1377542
##Bacteroidetes makes up 13% of all sequences

##% Cyanobacteria
Cyano = subset_taxa(fnm_r_final, Phylum == "Cyanobacteria")
sum(taxa_sums(Cyano))
31866/1377542
##Cyanobacteria makes up 2% of all sequences

##% Actinobacteria
Actino = subset_taxa(fnm_r_final, Phylum == "Actinobacteria")
sum(taxa_sums(Actino))
74877/1377542
##Actinobacteria makes up 5% of all sequences

##Verrucomicrobia
Verruc = subset_taxa(fnm_r_final, Phylum == "Verrucomicrobia")
sum(taxa_sums(Verruc))
26331/1377542
##Verrucomicrobia makes up 2% of all sequences

####  distribution of samples for reduced ####
#Create data frame from filtered reads for easier data manipulation.
df_fnm_2 <- as(sample_data(fnm_r_final), "data.frame")
readsumsdf = data.frame(nreads = sort(taxa_sums(fnm_r_final), TRUE), sorted = 1:ntaxa(fnm_r_final), type = "OTUs")

readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(fnm_r_final),  TRUE), sorted = 1:nsamples(fnm_r_final), type = "Samples"))

title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")


#######################ADULTS##########################################
### Beta Diversity ###
adults_mun <- subset_samples(fnm_r_final, Category == "adult" & Fert_status == "fertilized")
sample_names(adults_mun)

jaccard_adult_mun <- distance(adults_mun,"jaccard", binary = T)
jacc.ord_adult_mun <-ordinate(adults_mun, method = "PCoA", jaccard_adult_mun)

bray_adult_mun <- distance(adults_mun,"bray")
bray.ord_adult_mun <-ordinate(adults_mun, method = "PCoA", bray_adult_mun)

## PERMANOVA
## make dataframe for vegan to use for analysis
df_adult_mun <- as(sample_data(adults_mun), "data.frame")

############## Sex effects ############################
#Info on the strata addition to adonis function from https://rdrr.io/rforge/vegan/man/adonis.html. 
"The strata argument keeps groups intact for a particular hypothesis test where one does not want to permute the data among particular groups. For instance, strata = B causes permutations among levels of A but retains data within levels of B (no permutation among levels of B)."
## From Belden et al. 2016 Panamanian frog species host unique skin bacterial communities
## we accounted for the body location in the models using the 'strata' argument in the adonis
## function in vegan package (Oksanen et al. 2013) in R.
## we are using strata to nest body-location within adult category
adonis(jaccard_adult_mun ~ Sex, strata = df_adult_mun$Body_location, data = df_adult_mun)
adonis(bray_adult_mun ~ Sex, strata = df_adult_mun$Body_location, data = df_adult_mun)
##################  BODY LOCATION COMPARISONS  ###########
adonis(jaccard_adult_mun ~ Body_location, data = df_adult_mun)
adonis(bray_adult_mun ~ Body_location, data = df_adult_mun)
#################Adult by species comparisons ##########################
adonis(jaccard_adult_mun ~ Species, strata = df_adult_mun$Body_location, data = df_adult_mun)
adonis(bray_adult_mun ~ Species, strata = df_adult_mun$Body_location, data = df_adult_mun)

################## NEST comparison  ###################

###### BETA-diversity  ######
nest <- subset_samples(fnm_r_final, Category %in% c("nest-in","nest-out") & Fert_status == "fertilized")
sample_names(nest)
nest

jaccard_nest <- distance(nest,"jaccard", binary = T)
jacc.ord_nest <-ordinate(nest, method = "PCoA", jaccard_nest)

bray_nest <- distance(nest,"bray")
bray.ord_nest <-ordinate(nest, method = "PCoA", bray_nest)

## PERMANOVA
## make dataframe for vegan to use for analysis
df_nest <- as(sample_data(nest), "data.frame")
adonis(jaccard_nest ~ Nest_location,strata = df_nest$Nest_no, data = df_nest)
adonis(bray_nest ~ Nest_location, strata = df_nest$Nest_no,data = df_nest)
pairwise.perm.manova(jaccard_nest, df_nest$Nest_location, nperm = 999)

### ANOSIM to look at variation within and between
nest.ano <- anosim(jaccard_nest, df_nest$Nest_location)
summary(nest.ano)
plot(nest.ano, ylab = "Dissimilarity rank")
boxplot(nest.ano$dis.rank~nest.ano$class.vec)

#Think that nest insides are significantly different and nest outsides are not.
nest_in <- subset_samples(nest, Category == "nest-in")
sample_names(nest_in)
df_nest_in <- as(sample_data(nest_in), "data.frame")
jaccard_nest_in <- distance(nest_in,"jaccard", binary = T)
bray_nest_in <- distance(nest_in,"bray")
adonis(jaccard_nest_in ~ Species, data = df_nest_in, strata = df_nest_in$Nest_no)
adonis(bray_nest_in ~ Species, data = df_nest_in, strata = df_nest_in$Nest_no)
#Tried ANOSIM and was significant, not sure what this means.
nest.in.ano <- anosim(jaccard_nest_in, df_nest_in$Nest_no)
summary(nest.in.ano)

nest_out <- subset_samples(nest, Category == "nest-out")
sample_names(nest_out)
df_nest_out <- as(sample_data(nest_out), "data.frame")
jaccard_nest_out <- distance(nest_out,"jaccard", binary = T)
bray_nest_out <- distance(nest_out,"bray")
adonis(jaccard_nest_out ~ Species, data = df_nest_out, strata = df_nest_out$Nest_no)
adonis(bray_nest_out ~ Species, data = df_nest_out, strata = df_nest_out$Nest_no)

#Significant differences between nest insides based on PERMANOVA.
pairwise.perm.manova(jaccard_nest_in, df_nest_in$Species, nperm = 999)
pairwise.perm.manova(bray_nest_in, df_nest_in$Species, nperm = 999)

##### We now are going to move forward with cloacal microbiome  #####

## Removing dorsal, ventral swabs from adult, unfertilized adults and nest,
## sterile-swab (no replication), randomly select one swab from inside and
## one swab from outside

'%ni%' <- Negate('%in%')
fnm_cf <- subset_samples(fnm_r_final, Body_location %ni% c("dorsal", "ventral") & 
                           X.SampleID %ni% c ('AFCN2','AFDN2', 'AFVN2',  'AMC',
                                              'AMD', 'AMV','N2-1', 'N2-2', 'N2-3', 'SterileSwab', 'N1in1','N1in2',
                                              'N1out2','N1out3','N3in1','N3in2','N3out2','N3out3','N4in2',
                                              'N4in3','N4out2', 'N4out3', 'AMSecretion','NAR2'))

sample_names(fnm_cf)

##########################  MAIN ANALYSES  #####################
## Now, we just want all the adult cloaca, nest, tadpoles-nest, tadpoles-pond, water-enclosure,
## leaf, and water-terrarium
## frog nest microbiome - cloaca fertilized all
fnm_cfa <- subset_samples(fnm_cf, Category %in% c("adult", "nest-in","nest-out",'tadpole-nest', 
                                                  'tadpole-pond', 'water-enclosure', 'leaf', 'water-terrarium'))

fnm_cfa
#This says there are 26 taxa present which sum to 0.
sum(taxa_sums(fnm_cfa) == 0)
#Filter out any taxa that sum to 0.
filt_fnm_cfa = filter_taxa(fnm_cfa, function(x)sum(x) != 0, TRUE)
#Now we see that there are 2824 taxa by 58 samples in our OTU table.
filt_fnm_cfa

sample_names(fnm_cfa)
sort(sample_sums(fnm_cfa))

jaccard_cfa <- distance(fnm_cfa,"jaccard", binary = T)
jacc.ord_cfa <-ordinate(fnm_cfa, method = "PCoA", jaccard_cfa)

bray_cfa <- distance(fnm_cfa,"bray")
bray.ord_cfa <-ordinate(fnm_cfa, method = "PCoA", bray_cfa)

## PERMANOVA
## make dataframe for vegan to use for analysis
df_cfa <- as(sample_data(fnm_cfa), "data.frame")

#strata by next number
adonis(jaccard_cfa ~ Category, strata = df_cfa$Nest_no, data = df_cfa)
adonis(bray_cfa ~ Category, strata = df_cfa$Nest_no, data = df_cfa)
pairwise.perm.manova(jaccard_cfa, df_cfa$Category, nperm = 999)
pairwise.perm.manova(bray_cfa, df_cfa$Category, nperm = 999)

#######Looking to see if again the nest insides are significantly different from one another.
fnm_cfa_nest_in <- subset_samples(fnm_cfa, Category == "nest-in")
sample_names(fnm_cfa_nest_in)
df_cfa_nest_in <- as(sample_data(fnm_cfa_nest_in), "data.frame")
jaccard_cfa_nest_in <- distance(fnm_cfa_nest_in,"jaccard", binary = T)
bray_cfa_nest_in <- distance(fnm_cfa_nest_in,"bray")
adonis(jaccard_cfa_nest_in ~ Species, data = df_cfa_nest_in)
adonis(bray_cfa_nest_in ~ Species, data = df_cfa_nest_in)

#### Now let's look alpha-diversity for these 
########### summarize samples for SR and PD  #########
sum_data_SR <- ddply(df_cfa, c("Category"), summarise,
                     N    = length(Category),
                     mean = mean(SR),
                     median = median(SR),
                     sd   = sd(SR),
                     se   = sd / sqrt(N))
sum_data_SR

sum_data_PD <- ddply(df_cfa, c("Category"), summarise,
                     N    = length(Category),
                     mean = mean(PD),
                     median = median(PD),
                     sd   = sd(PD),
                     se   = sd / sqrt(N))
sum_data_PD



hist(df_cfa$SR)
hist(df_cfa$PD)

## >0.05 data is normally distributed
## close enough for species richness, fits within intervals
qqnorm(df_cfa$SR^2)
shapiro.test((df_cfa$SR)^2)
hist((df_cfa$SR)^2)

qqnorm(df_cfa$PD^2)
shapiro.test((df_cfa$PD)^2)
hist((df_cfa$PD)^2)

## > 0.05 variances are similar, meets assumption of ANOVA
# homogeniety of variance
leveneTest((SR^2) ~ Category, data = df_cfa)
leveneTest((PD^2) ~ Category, data = df_cfa)

aov_cfa <- aov(SR^2 ~ Category, data = df_cfa)
summary(aov_cfa)

TukeyHSD(aov_cfa)

aov_cfaPD <- aov(PD^2 ~ Category, data = df_cfa)
aov_cfaPD_2 <- aov(PD ~ Category, data = df_cfa)
summary(aov_cfaPD)

TukeyHSD(aov_cfaPD)

max(df_cfa$PD)
min(df_cfa$PD)
sd(df_cfa$PD)

#Tested with Kruskal-Walis test for just nest-in (data not normally distributed)
fnm_nest_cfa <- subset_samples(fnm_cfa, Category == "nest-in")
sample_names(fnm_nest_cfa)
df_fnm_nest_cfa <- as(sample_data(fnm_nest_cfa), "data.frame")
qqp((df_fnm_nest_cfa$SR^2), "norm")
shapiro.test((df_fnm_nest_cfa$SR)^2)
shapiro.test((df_fnm_nest_cfa$PD)^2)

kruskal.test(SR^2 ~ Nest_no, data = df_fnm_nest_cfa)
kruskal.test(PD^2 ~ Nest_no, data = df_fnm_nest_cfa)

#Nest-out
fnm_nest_out_cfa <- subset_samples(fnm_cfa, Category == "nest-out")
sample_names(fnm_nest_out_cfa)
df_fnm_nest_out_cfa <- as(sample_data(fnm_nest_out_cfa), "data.frame")
qqp((df_fnm_nest_out_cfa$SR^2), "norm")
shapiro.test((df_fnm_nest_out_cfa$SR)^2)
shapiro.test((df_fnm_nest_out_cfa$PD)^2)

## > 0.05 variances are similar, meets assumption of ANOVA
# homogeniety of variance
leveneTest((SR^2) ~ Nest_no, data = df_fnm_nest_out_cfa)
leveneTest((PD^2) ~ Category, data = df_cfa)

aov__fnm_nest_out_cfa <- aov(SR^2 ~ Nest_no, data = df_fnm_nest_out_cfa)
summary(aov__fnm_nest_out_cfa)

aov__fnm_nest_out_cfa <- aov(PD^2 ~ Nest_no, data = df_fnm_nest_out_cfa)
summary(aov__fnm_nest_out_cfa)

kruskal.test(SR^2 ~ Nest_no, data = df_fnm_nest_out_cfa)
kruskal.test(PD^2 ~ Nest_no, data = df_fnm_nest_out_cfa)


##Plot of differences in phylogenetic diversity between categories
tukey_result_1 <- HSD.test(aov_cfa, "Category", group = TRUE)
print(tukey_result_1)

tukey_result_2 <- HSD.test(aov_cfaPD_2, "Category", group = TRUE)
print(tukey_result_2)

# Plot result phylogenetic diversity
group_data <- tukey_result_2$groups[order(rownames(tukey_result_2$groups)),]

#Manuscript figure 3.
tiff("alpha_diversity_PD.tiff", units = "mm", width = 297, height = 170, res = 300)
my_plot <- ggplot(df_cfa, aes(x = Category, y = PD)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data),
                y = max(df_cfa$PD) + 1,
                label = group_data$groups),
            col = 'black',
            size = 8, vjust = 0) +
  geom_boxplot()+
  theme_bw() +
  theme_classic() +
  theme(text=element_text(size = 18),
        axis.title.x = element_blank()) +
  labs(y="Phylogenetic Diversity") +
  scale_x_discrete(labels=c("Adults", "Leaf", "Nest-Inside","Nest-Outside","Tadpole-Nest","Tadpole-Pond","Water-Pond","Water-Terrarium"))
my_plot
dev.off()

# Plot result species richness
group_data_SR <- tukey_result_1$groups[order(rownames(tukey_result_1$groups)),]
my_plot_SR <- ggplot(df_cfa, aes(x = Category, y = SR)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data),
                y = max(df_cfa$SR) + 1,
                label = group_data$groups),
            col = 'black',
            size = 10) +
  geom_boxplot()+
  theme_bw() +
  theme_classic() +
  theme(text=element_text(size = 25),
        axis.title.x = element_blank()) +
  labs(y="Species Richenss") +
  scale_x_discrete(labels=c("Adult", "Leaf", "Nest-In","Nest-Out","Tad-Nest","Tad-Pond","Water-Pond","Water-Terrarium"))
my_plot_SR

######################### Comparison of Tadpoles before/after Environmental Interaction###############
## frog nest microbiome - tadpole comparisons (biological replicates because each tadpole swabbed is a different individual)
fnm_tads <- subset_samples(fnm_cfa, Category %in% c('tadpole-nest', 'tadpole-pond'))

##only tadpoles from foam nest (FTad) and from pond (WTad)
sample_names(fnm_tads)

jaccard_tads <- distance(fnm_tads,"jaccard", binary = T)
jacc.ord_tads <-ordinate(fnm_tads, method = "PCoA", jaccard_tads)

bray_tads <- distance(fnm_tads,"bray")
bray.ord_tads <-ordinate(fnm_tads, method = "PCoA", bray_tads)

## PERMANOVA

## make dataframe for vegan to use for analysis
df_tads <- as(sample_data(fnm_tads), "data.frame")

#strata by nest number.
adonis(jaccard_tads ~ Category, strata = df_tads$Nest_no, data = df_tads)
adonis(jaccard_tads ~ Category, data = df_tads)
adonis(bray_tads ~ Category, strata = df_tads$Nest_no, data = df_tads)

pairwise.perm.manova(jaccard_tads, df_tads$Category, nperm = 999)
pairwise.perm.manova(bray_tads, df_tads$Category, nperm = 999)

## PERMDISP
groups_tads <- df_tads[["Category"]]
jacc_disp_tads <-betadisper(jaccard_tads, groups_tads, type=c("median"))
plot(jacc_disp_tads)
anova(jacc_disp_tads)
boxplot(jacc_disp_tads)
TukeyHSD(jacc_disp_tads)

bray_disp_tads <- betadisper(bray_tads, groups_tads, type = c("median"))
plot(bray_disp_tads)

#####################################Checking basic info for amphibs#######################
###########Just amphibs######################
## total number of sequences per sample, sorted in order
#Need all adult samples, not just the cloacal microbiome, so subsetting from more complete phyloseq object.
'%ni%' <- Negate('%in%')
fnm_taxa_plots <- subset_samples(fnm_r_final, X.SampleID %ni% c ('AFCN2','AFDN2', 'AFVN2',
                                                                 'AMC','AMD', 'AMV','N2-1', 
                                                                 'N2-2', 'N2-3', 'SterileSwab',
                                                                 'AMSecretion','NAR2',
                                                                 'H2Ocontrol1','H2Ocontrol2'
                                                                 ,'H2Ocontrol3','H20controlN4'
                                                                 ,'H2ON3','N1H2OFTads'))
fnm_amphibs <- subset_samples(fnm_taxa_plots, Category %in% c("adult","tadpole-nest","tadpole-pond"))
sample_data(fnm_amphibs)

#writing csv file for amphibian-specific microbiome comparison to global amphibian microbiome dataset (obtained by reaching out to corresponding author of Kueneman et al., 2019).
otus_amphibs <- otu_table(fnm_amphibs)
write.csv(otus_amphibs,file = 'otu_table_amphibs.csv')
tax_table_fnm_amphibs <- tax_table(fnm_amphibs)
write.csv(tax_table_fnm_amphibs,file = 'tax_table_amphibs.csv')
sample_names(sample_data(fnm_amphibs))

sum(tax_table(fnm_amphibs))
top15ph <- sort(tapply(taxa_sums(fnm_amphibs), tax_table(fnm_amphibs)[, "Phylum"],sum), TRUE)[1:15]
top15g <- sort(tapply(taxa_sums(fnm_amphibs), tax_table(fnm_amphibs)[, "Genus"], sum), TRUE)[1:15]
top15g
top15ph

get_taxa_unique(fnm_amphibs, "Kingdom")

fnm_Archae = subset_taxa(fnm_amphibs, Kingdom == "Archaea")

## these are your Archael phyla
get_taxa_unique(fnm_Archae, "Phylum")

## 28 bacterial and archael phyla (one NA, and 2 Archaeal phylum)
get_taxa_unique(fnm_amphibs, "Phylum")
tax_table(fnm_amphibs)[,"Phylum"] %>% unique %>% na.exclude %>% length

#### 512,532 total sequences
sum(taxa_sums(fnm_amphibs))

#This says there are 329 taxa present which sum to 0.
sum(taxa_sums(fnm_amphibs) == 0)
#Filter out any taxa that sum to 0.
filt_fnm_amphibs = filter_taxa(fnm_amphibs, function(x)sum(x) != 0, TRUE)
#Now we see that there are 2521 taxa by 36 samples in our OTU table.
filt_fnm_amphibs

## % Proteobacteria
Proteo = subset_taxa(fnm_amphibs, Phylum == "Proteobacteria")

## Proteobacteria
sum(taxa_sums(Proteo))

## Proteobacteria makes up 45% of all sequences
232924/512532

## Firmicutes
Firmicutes = subset_taxa(fnm_amphibs, Phylum == "Firmicutes")
sum(taxa_sums(Firmicutes))
94494/512532
##Firmicutes makes up 18% of all sequences

##Bacteroidetes
Bacteroidetes = subset_taxa(fnm_amphibs, Phylum == "Bacteroidetes")
sum(taxa_sums(Bacteroidetes))
83099/512532
##Bacteroidetes makes up 16% of all sequences

##% Cyanobacteria
Cyano = subset_taxa(fnm_amphibs, Phylum == "Cyanobacteria")
sum(taxa_sums(Cyano))
10062/512532
##Cyanobacteria makes up 2% of all sequences

##% Fusobacteria
Fuso = subset_taxa(fnm_amphibs, Phylum == "Fusobacteria")
sum(taxa_sums(Fuso))
1788/512532
##Fusobacteria makes up 0.3% of all sequences

##% Actinobacteria
Actino = subset_taxa(fnm_amphibs, Phylum == "Actinobacteria")
sum(taxa_sums(Actino))
33571/512532
##Actinobacteria makes up 7% of all sequences

##Verrucomicrobia
Verruc = subset_taxa(fnm_amphibs, Phylum == "Verrucomicrobia")
sum(taxa_sums(Verruc))
24011/512532
##Verrucomicrobia makes up 5% of all sequences

##Bordetella (Genus)
Bord = subset_taxa(fnm_amphibs, Genus == "Bordetella")
sum(taxa_sums(Bord))
15467/512532
##3% of sequences

#Acinetobacter
53147/512532
#10%

#Flavobacterium
22031/512532
#4%

#Pseudomonas
19143/512532
#4%
#Lactobacillus
11543/512532
#2%
## Now, we want to check for significant differences between species.

jaccard_amphibs <- distance(fnm_amphibs,"jaccard", binary = T)
bray_amphibs <- distance(fnm_amphibs,"bray")

## PERMANOVA

## make dataframe for vegan to use for analysis
df_amphibs <- as(sample_data(fnm_amphibs), "data.frame")

## From Belden et al. 2016? Panamanian frog species host unique skin bacterial communities
## we accounted for the nest_no in the models using the 'strata' argument in the adonis
## function in vegan package (Oksanen et al. 2013) in R.
adonis(jaccard_amphibs ~ Category, strata = df_amphibs$Species, data = df_amphibs)
adonis(jaccard_amphibs ~ Species, data = df_amphibs)

pairwise.perm.manova(jaccard_amphibs, df_amphibs$Species, nperm = 999)
pairwise.perm.manova(bray_amphibs, df_amphibs$Species, nperm = 999)

fnm_amphibs_adults <- subset_samples(fnm_taxa_plots, Category == "adult")
sample_names(fnm_amphibs_adults)
df_amphibs_adults <- as(sample_data(fnm_amphibs_adults), "data.frame")
jaccard_amphibs_adults <- distance(fnm_amphibs_adults,"jaccard", binary = T)
bray_amphibs_adults <- distance(fnm_amphibs_adults,"bray")
adonis(jaccard_amphibs_adults ~ Species, data = df_amphibs_adults)
adonis(bray_amphibs_adults ~ Species, data=df_amphibs_adults)
pairwise.perm.manova(jaccard_amphibs_adults, df_amphibs_adults$Species, nperm = 999)
pairwise.perm.manova(bray_amphibs_adults, df_amphibs_adults$Species, nperm = 999)
pairwise.perm.manova(jaccard_amphibs_adults, df_amphibs_adults$Sex, nperm = 999)
pairwise.perm.manova(bray_amphibs_adults, df_amphibs_adults$Sex, nperm = 999)
pairwise.perm.manova(jaccard_amphibs_adults, df_amphibs_adults$Body_location, nperm = 999)
pairwise.perm.manova(bray_amphibs_adults, df_amphibs_adults$Body_location, nperm = 999)


fnm_amphibs_tadnest <- subset_samples(fnm_taxa_plots, Category == "tadpole-nest")
sample_names(fnm_amphibs_tadnest)
df_amphibs_tadnest <- as(sample_data(fnm_amphibs_tadnest), "data.frame")
jaccard_amphibs_tadnest <- distance(fnm_amphibs_tadnest,"jaccard", binary = T)
bray_amphibs_tadnest <- distance(fnm_amphibs_tadnest,"bray")
adonis(jaccard_amphibs_tadnest ~ Species, data = df_amphibs_tadnest, strata = df_amphibs_tadnest$Nest_no)
adonis(bray_amphibs_tadnest ~ Species, data = df_amphibs_tadnest, strata = df_amphibs_tadnest$Nest_no)

fnm_amphibs_tadpond <- subset_samples(fnm_taxa_plots, Category == "tadpole-pond")
sample_names(fnm_amphibs_tadpond)
df_amphibs_tadpond <- as(sample_data(fnm_amphibs_tadpond), "data.frame")
jaccard_amphibs_tadpond <- distance(fnm_amphibs_tadpond,"jaccard", binary = T)
bray_amphibs_tadpond <- distance(fnm_amphibs_tadpond,"bray")
adonis(jaccard_amphibs_tadpond ~ Species, data = df_amphibs_tadpond, strata = df_amphibs_tadpond$Nest_no)
adonis(bray_amphibs_tadpond ~ Species, data = df_amphibs_tadpond, strata = df_amphibs_tadpond$Nest_no)

fnm_amphibs_nest_in <- subset_samples(fnm_taxa_plots, Category == "nest-in")
sample_names(fnm_amphibs_nest_in)
df_amphibs_nest_in <- as(sample_data(fnm_amphibs_nest_in), "data.frame")
jaccard_amphibs_nest_in <- distance(fnm_amphibs_nest_in,"jaccard", binary = T)
bray_amphibs_nest_in <- distance(fnm_amphibs_nest_in,"bray")
adonis(jaccard_amphibs_nest_in ~ Species, data = df_amphibs_nest_in, strata = df_amphibs_nest_in$Nest_no)
adonis(bray_amphibs_nest_in ~ Species, data = df_amphibs_nest_in, strata = df_amphibs_nest_in$Nest_no)
#Don't do pairwise if they aren't significant.
pairwise.perm.manova(jaccard_amphibs_nest_in, df_amphibs_nest_in$Species, nperm = 999)
pairwise.perm.manova(bray_amphibs_nest_in, df_amphibs_nest_in$Species, nperm = 999)


fnm_amphibs_nest_out <- subset_samples(fnm_taxa_plots, Category == "nest-out")
sample_names(fnm_amphibs_nest_out)
df_amphibs_nest_out <- as(sample_data(fnm_amphibs_nest_out), "data.frame")
jaccard_amphibs_nest_out <- distance(fnm_amphibs_nest_out,"jaccard", binary = T)
bray_amphibs_nest_out <- distance(fnm_amphibs_nest_out,"bray")
adonis(jaccard_amphibs_nest_out ~ Species, data = df_amphibs_nest_out)
adonis(bray_amphibs_nest_out ~ Species, data = df_amphibs_nest_out)
pairwise.perm.manova(jaccard_amphibs_nest_out, df_amphibs_nest_out$Species, nperm = 999)
pairwise.perm.manova(bray_amphibs_nest_out, df_amphibs_nest_out$Species, nperm = 999)


#Nests only
fnm_nest <- subset_samples(fnm_taxa_plots, Category %in% c("nest-in","nest-out"))

sample_names(fnm_nest)

get_taxa_unique(fnm_nest, "Kingdom")

fnm_Archae = subset_taxa(fnm_nest, Kingdom == "Archaea")

## these are your Archael phyla
get_taxa_unique(fnm_Archae, "Phylum")

## 28 bacterial and archael phyla (one NA, and 2 Archaeal phylum)
get_taxa_unique(fnm_nest, "Phylum")
tax_table(fnm_nest)[,"Phylum"] %>% unique %>% na.exclude %>% length

#### 244,169 total sequences
sum(taxa_sums(fnm_nest))

#### 2850 sOTUs
fnm_nest
#This says there are 1106 taxa present which sum to 0.
sum(taxa_sums(fnm_nest) == 0)
#Filter out any taxa that sum to 0.
filt_fnm_nest = filter_taxa(fnm_nest, function(x)sum(x) != 0, TRUE)
#Now we see that there are 1744 taxa by 18 samples in our OTU table.
filt_fnm_nest

top15ph_nest <- sort(tapply(taxa_sums(fnm_nest), tax_table(fnm_nest)[, "Phylum"],sum), TRUE)[1:15]
top15ph_nest

## % Proteobacteria
Proteo = subset_taxa(fnm_nest, Phylum == "Proteobacteria")

## Proteobacteria
sum(taxa_sums(Proteo))

## Proteobacteria makes up 67% of all sequences
164676/244169

## Firmicutes
Firmicutes = subset_taxa(fnm_nest, Phylum == "Firmicutes")
sum(taxa_sums(Firmicutes))
23789/244169
##Firmicutes makes up 10% of all sequences

##Bacteroidetes
Bacteroidetes = subset_taxa(fnm_nest, Phylum == "Bacteroidetes")
sum(taxa_sums(Bacteroidetes))
43597/244169
##Bacteroidetes makes up 18% of all sequences

##% Cyanobacteria
Cyano = subset_taxa(fnm_nest, Phylum == "Cyanobacteria")
sum(taxa_sums(Cyano))
854/244169
##Cyanobacteria makes up 0.3% of all sequences

##% Actinobacteria
Actino = subset_taxa(fnm_nest, Phylum == "Actinobacteria")
sum(taxa_sums(Actino))
9319/244169
##Actinobacteria makes up 4% of all sequences

#% Tenericutes
Ternero = subset_taxa(fnm_nest, Phylum == "Tenericutes")
sum(taxa_sums(Ternero))
908/244169
#Tenericutes makes up 0.4% of all sequences

##Verrucomicrobia
Verruc = subset_taxa(fnm_nest, Phylum == "Verrucomicrobia")
sum(taxa_sums(Verruc))
96/244169
##Verrucomicrobia makes up 0.03% of all sequences

#Environment only
fnm_environ <- subset_samples(fnm_taxa_plots, Category %in% c("leaf","water-terrarium","water-enclosure"))

sample_names(sample_data(fnm_environ))

get_taxa_unique(fnm_environ, "Kingdom")

fnm_Archae = subset_taxa(fnm_environ, Kingdom == "Archaea")

## these are your Archael phyla
get_taxa_unique(fnm_Archae, "Phylum")

## 28 bacterial and archael phyla (one NA, and 2 Archaeal phylum)
get_taxa_unique(fnm_environ, "Phylum")
tax_table(fnm_environ)[,"Phylum"] %>% unique %>% na.exclude %>% length

#### 381,901 total sequences
sum(taxa_sums(fnm_environ))

#### 2850 sOTUs
fnm_environ

## % Proteobacteria
Proteo = subset_taxa(fnm_environ, Phylum == "Proteobacteria")

## Proteobacteria
sum(taxa_sums(Proteo))

## Proteobacteria makes up 62% of all sequences
237210/381901

## Firmicutes
Firmicutes = subset_taxa(fnm_environ, Phylum == "Firmicutes")
sum(taxa_sums(Firmicutes))
58898/381901
##Firmicutes makes up 15% of all sequences

##Bacteroidetes
Bacteroidetes = subset_taxa(fnm_environ, Phylum == "Bacteroidetes")
sum(taxa_sums(Bacteroidetes))
38195/381901
##Bacteroidetes makes up 10% of all sequences

##% Cyanobacteria
Cyano = subset_taxa(fnm_environ, Phylum == "Cyanobacteria")
sum(taxa_sums(Cyano))
19770/381901
##Cyanobacteria makes up 5% of all sequences

##% Actinobacteria
Actino = subset_taxa(fnm_environ, Phylum == "Actinobacteria")
sum(taxa_sums(Actino))
19379/381901
##Actinobacteria makes up 5% of all sequences

##Verrucomicrobia
Verruc = subset_taxa(fnm_environ, Phylum == "Verrucomicrobia")
sum(taxa_sums(Verruc))
2137/381901
##Verrucomicrobia makes up 0.5% of all sequences

#############################Testing Other Replicates################################

############################### REP1 ##################################################
#Testing other replicates to see if we get the same results.
'%ni%' <- Negate('%in%')
fnm_cf_rep1 <- subset_samples(fnm_r_final, Body_location %ni% c("cloaca", "ventral") & 
                                X.SampleID %ni% c ('AFCN2','AFDN2', 'AFVN2',  'AMC',
                                                   'AMD', 'AMV','N2-1', 'N2-2', 'N2-3',
                                                   'SterileSwab', 'N1in2','N1in3',
                                                   'N1out1','N1out2','N3in2','N3in3','N3out1'
                                                   ,'N3out2','N4in1','N4in2','N4out1', 'N4out2',
                                                   'AMSecretion','NAR2'))
sample_names(fnm_cf_rep1)

## frog nest microbiome - cloaca fertilized all
fnm_cf_rep1 <- subset_samples(fnm_cf_rep1, Category %in% c("adult", "nest-in","nest-out",'tadpole-nest', 
                                                           'tadpole-pond', 'water-enclosure', 'leaf', 'water-terrarium'))



sample_names(fnm_cf_rep1)
sort(sample_sums(fnm_cf_rep1))


jaccard_cf_rep1 <- distance(fnm_cf_rep1,"jaccard", binary = T)
jacc.ord_cf_rep1 <-ordinate(fnm_cf_rep1, method = "PCoA", jaccard_cf_rep1)

p_cf_rep1 <- plot_ordination(fnm_cf_rep1, jacc.ord_cf_rep1, color = "Category", shape = "Nest_no")
p_cf_rep1 + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

uni_cf_rep1 <- distance(fnm_cf_rep1,"uunifrac")
## This is just telling us that there are more branches 5427 than the OTUs (2714)
## it's fine

uni.ord_cf_rep1 <-ordinate(fnm_cf_rep1, method = "PCoA", uni_cf_rep1)

p_cf_rep1 <- plot_ordination(fnm_cf_rep1, uni.ord_cf_rep1, color = "Category", shape = "Nest_no")
p_cf_rep1 + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 


bray_cf_rep1 <- distance(fnm_cf_rep1,"bray")

bray.ord_cf_rep1 <-ordinate(fnm_cf_rep1, method = "PCoA", bray_cf_rep1)

p_cf_rep1 <- plot_ordination(fnm_cf_rep1, bray.ord_cf_rep1, color = "Category", shape = "Nest_no")
p_cf_rep1 + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 


## PERMANOVA

## make dataframe for vegan to use for analysis
df_cf_rep1 <- as(sample_data(fnm_cf_rep1), "data.frame")

adonis(jaccard_cf_rep1 ~ Category, strata = df_cf_rep1$Nest_no, data = df_cf_rep1)
adonis(jaccard_cf_rep1 ~ Category, data = df_cf_rep1)
pairwise.perm.manova(jaccard_cf_rep1, df_cf_rep1$Category, nperm = 999)

adonis(uni_cf_rep1 ~ Category, strata = df_cf_rep1$Nest_no, data = df_cf_rep1)
pairwise.perm.manova(uni_cf_rep1, df_cf_rep1$Category, nperm = 999)

adonis(bray_cf_rep1 ~ Category, strata = df_cf_rep1$Nest_no, data = df_cf_rep1)
pairwise.perm.manova(bray_cf_rep1, df_cf_rep1$Category, nperm = 999)

## PERMDISP
groups_cf_rep1 <- df_cf_rep1[["Category"]]
jacc_disp_cf_rep1 <-betadisper(jaccard_cf_rep1, groups_cf_rep1, type=c("median"))
plot(jacc_disp_cf_rep1)
anova(jacc_disp_cf_rep1)
boxplot(jacc_disp_cf_rep1)
TukeyHSD(jacc_disp_cf_rep1)

fnm_cf_rep1_nest_in <- subset_samples(fnm_cf_rep1, Category == "nest-in")
sample_names(fnm_cf_rep1_nest_in)
df_cf_rep1_nest_in <- as(sample_data(fnm_cf_rep1_nest_in), "data.frame")
jaccard_cf_rep1_nest_in <- distance(fnm_cf_rep1_nest_in,"jaccard", binary = T)
bray_cf_rep1_nest_in <- distance(fnm_cf_rep1_nest_in,"bray")
adonis(jaccard_cf_rep1_nest_in ~ Species, data = df_cf_rep1_nest_in)
adonis(bray_cf_rep1_nest_in ~ Species, data = df_cf_rep1_nest_in)

groups_cf_rep1_nest_in <- df_cf_rep1_nest_in[["Nest_no"]]
jacc_disp_cf_rep1_nest_in <-betadisper(jaccard_cf_rep1_nest_in, groups_cf_rep1_nest_in, type=c("median"))
plot(jacc_disp_cf_rep1_nest_in)
anova(jacc_disp_cf_rep1_nest_in)
boxplot(jacc_disp_cf_rep1_nest_in)
TukeyHSD(jacc_disp_cf_rep1_nest_in)


#### Now let's look alpha-diversity for these 

########### summarize samples for SR and PD  #########
library(plyr)
sum_data_SR <- ddply(df_cf_rep1, c("Category"), summarise,
                     N    = length(Category),
                     mean = mean(SR),
                     median = median(SR),
                     sd   = sd(SR),
                     se   = sd / sqrt(N))
sum_data_SR

sum_data_PD <- ddply(df_cf_rep1, c("Category"), summarise,
                     N    = length(Category),
                     mean = mean(PD),
                     median = median(PD),
                     sd   = sd(PD),
                     se   = sd / sqrt(N))
sum_data_PD



hist(df_cf_rep1$SR)
hist(df_cf_rep1$PD)

## >0.05 data is normally distributed
## close enough for species richness, fits within intervals
qqp((df_cf_rep1$SR^2), "norm")
shapiro.test((df_cf_rep1$SR)^2)

hist((df_cf_rep1$SR)^2)

shapiro.test((df_cf_rep1$PD)^2)

## > 0.05 variances are similar, meets assumption of ANOVA
# homogeniety of variance
leveneTest((SR^2) ~ Category, data = df_cf_rep1)
leveneTest((PD^2) ~ Category, data = df_cf_rep1)

aov_cf_rep1 <- aov(SR^2 ~ Category, data = df_cf_rep1)
summary(aov_cf_rep1)

plot(SR^2 ~ Category, data = df_cf_rep1)
plot(SR ~ Category, data = df_cf_rep1)

TukeyHSD(aov_cf_rep1)

aov_cf_rep1PD <- aov(PD^2 ~ Category, data = df_cf_rep1)
aov_cf_rep1PD_2 <- aov(PD ~ Category, data = df_cf_rep1)
summary(aov_cf_rep1PD)

plot(PD ~ Category, data = df_cf_rep1)
plot(PD^2 ~ Category, data = df_cf_rep1)

TukeyHSD(aov_cf_rep1PD)

max(df_cf_rep1$PD)
min(df_cf_rep1$PD)
sd(df_cf_rep1$PD)

##Plot of differences in phylogenetic diversity between categories
tukey_result_1 <- HSD.test(aov_cf_rep1, "Category", group = TRUE)
print(tukey_result_1)

tukey_result_2 <- HSD.test(aov_cf_rep1PD_2, "Category", group = TRUE)
print(tukey_result_2)

# Plot result
group_data <- tukey_result_2$groups[order(rownames(tukey_result_2$groups)),]
my_plot <- ggplot(df_cf_rep1, aes(x = Category, y = PD)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data),
                y = max(df_cf_rep1$PD) + 1,
                label = group_data$groups),
            col = 'black',
            size = 10) +
  geom_boxplot()+
  theme(text=element_text(size = 20),
        axis.title.x = element_blank())
my_plot

#Tested with Kruskal-Walis test for just nest-in (data not normally distributed)
fnm_nest_in_rep1 <- subset_samples(fnm_cf_rep1, Category == "nest-in")
sample_names(fnm_nest_in_rep1)
df_fnm_nest_in_rep1 <- as(sample_data(fnm_nest_in_rep1), "data.frame")
qqp((df_fnm_nest_in_rep1$SR^2), "norm")
shapiro.test((df_fnm_nest_in_rep1$SR)^2)
shapiro.test((df_fnm_nest_in_rep1$PD)^2)


kruskal.test(SR^2 ~ Nest_no, data = df_fnm_nest_in_rep1)
kruskal.test(PD^2 ~ Nest_no, data = df_fnm_nest_in_rep1)

#Nest-out
fnm_nest_out_rep1 <- subset_samples(fnm_cf_rep1, Category == "nest-out")
sample_names(fnm_nest_out_rep1)
df_fnm_nest_out_rep1 <- as(sample_data(fnm_nest_out_rep1), "data.frame")
shapiro.test((df_fnm_nest_out_rep1$SR)^2)
shapiro.test((df_fnm_nest_out_rep1$PD)^2)


kruskal.test(SR^2 ~ Nest_no, data = df_fnm_nest_out_rep1)
kruskal.test(PD^2 ~ Nest_no, data = df_fnm_nest_out_rep1)
################################ REP2 #################################################
#Testing other replicates to see if we get the same results.
'%ni%' <- Negate('%in%')
fnm_cf_rep2 <- subset_samples(fnm_r_final, Body_location %ni% c("cloaca", "dorsal") & 
                                X.SampleID %ni% c ('AFCN2','AFDN2', 'AFVN2',  'AMC',
                                                   'AMD', 'AMV','N2-1', 'N2-2', 'N2-3',
                                                   'SterileSwab', 'N1in1','N1in3',
                                                   'N1out1','N1out3','N3in1','N3in3'
                                                   ,'N3out1','N3out3','N4in1','N4in3'
                                                   ,'N4out1', 'N4out3','AMSecretion','NAR2'
                                ))
sample_names(fnm_cf_rep2)
## frog nest microbiome - cloaca fertilized all
fnm_cf_rep2 <- subset_samples(fnm_cf_rep2, Category %in% c("adult", "nest-in","nest-out",'tadpole-nest', 
                                                           'tadpole-pond', 'water-enclosure', 'leaf', 'water-terrarium'))



sample_names(fnm_cf_rep2)
sort(sample_sums(fnm_cf_rep2))


jaccard_cf_rep2 <- distance(fnm_cf_rep2,"jaccard", binary = T)
jacc.ord_cf_rep2 <-ordinate(fnm_cf_rep2, method = "PCoA", jaccard_cf_rep2)

p_cf_rep2 <- plot_ordination(fnm_cf_rep2, jacc.ord_cf_rep2, color = "Category", shape = "Nest_no")
p_cf_rep2 + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

uni_cf_rep2 <- distance(fnm_cf_rep2,"uunifrac")
## This is just telling us that there are more branches 5427 than the OTUs (2714)
## it's fine

uni.ord_cf_rep2 <-ordinate(fnm_cf_rep2, method = "PCoA", uni_cf_rep2)

p_cf_rep2 <- plot_ordination(fnm_cf_rep2, uni.ord_cf_rep2, color = "Category", shape = "Nest_no")
p_cf_rep2 + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 


bray_cf_rep2 <- distance(fnm_cf_rep2,"bray")

bray.ord_cf_rep2 <-ordinate(fnm_cf_rep2, method = "PCoA", bray_cf_rep2)

p_cf_rep2 <- plot_ordination(fnm_cf_rep2, bray.ord_cf_rep2, color = "Category", shape = "Nest_no")
p_cf_rep2 + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 


## PERMANOVA

## make dataframe for vegan to use for analysis
df_cf_rep2 <- as(sample_data(fnm_cf_rep2), "data.frame")

adonis(jaccard_cf_rep2 ~ Category, strata = df_cf_rep2$Nest_no, data = df_cf_rep2)
adonis(jaccard_cf_rep2 ~ Category, data = df_cf_rep2)
pairwise.perm.manova(jaccard_cf_rep2, df_cf_rep2$Category, nperm = 999)

adonis(uni_cf_rep2 ~ Category, strata = df_cf_rep2$Nest_no, data = df_cf_rep2)
pairwise.perm.manova(uni_cf_rep2, df_cf_rep2$Category, nperm = 999)

adonis(bray_cf_rep2 ~ Category, strata = df_cf_rep2$Nest_no, data = df_cf_rep2)
pairwise.perm.manova(bray_cf_rep2, df_cf_rep2$Category, nperm = 999)

## PERMDISP
groups_cf_rep2 <- df_cf_rep2[["Category"]]
jacc_disp_cf_rep2 <-betadisper(jaccard_cf_rep2, groups_cf_rep2, type=c("median"))
plot(jacc_disp_cf_rep2)
anova(jacc_disp_cf_rep2)
boxplot(jacc_disp_cf_rep2)
TukeyHSD(jacc_disp_cf_rep2)


#### Now let's look alpha-diversity for these 

########### summarize samples for SR and PD  #########
sum_data_SR <- ddply(df_cf_rep2, c("Category"), summarise,
                     N    = length(Category),
                     mean = mean(SR),
                     median = median(SR),
                     sd   = sd(SR),
                     se   = sd / sqrt(N))
sum_data_SR

sum_data_PD <- ddply(df_cf_rep2, c("Category"), summarise,
                     N    = length(Category),
                     mean = mean(PD),
                     median = median(PD),
                     sd   = sd(PD),
                     se   = sd / sqrt(N))
sum_data_PD



hist(df_cf_rep2$SR)
hist(df_cf_rep2$PD)

## >0.05 data is normally distributed
## close enough for species richness, fits within intervals
qqp((df_cf_rep2$SR^2), "norm")
shapiro.test((df_cf_rep2$SR)^2)

hist((df_cf_rep2$SR)^2)

shapiro.test((df_cf_rep2$PD)^2)

## > 0.05 variances are similar, meets assumption of ANOVA
# homogeniety of variance
leveneTest((SR^2) ~ Category, data = df_cf_rep2)
leveneTest((PD^2) ~ Category, data = df_cf_rep2)

aov_cf_rep2 <- aov(SR^2 ~ Category, data = df_cf_rep2)
summary(aov_cf_rep2)

plot(SR^2 ~ Category, data = df_cf_rep2)
plot(SR ~ Category, data = df_cf_rep2)

TukeyHSD(aov_cf_rep2)

aov_cf_rep2PD <- aov(PD^2 ~ Category, data = df_cf_rep2)
aov_cf_rep2PD_2 <- aov(PD ~ Category, data = df_cf_rep2)
summary(aov_cf_rep2PD)

plot(PD ~ Category, data = df_cf_rep2)
plot(PD^2 ~ Category, data = df_cf_rep2)

TukeyHSD(aov_cf_rep2PD)

max(df_cf_rep2$PD)
min(df_cf_rep2$PD)
sd(df_cf_rep2$PD)

##Plot of differences in phylogenetic diversity between categories
tukey_result_2 <- HSD.test(aov_cf_rep2, "Category", group = TRUE)
print(tukey_result_2)

tukey_result_2 <- HSD.test(aov_cf_rep2PD_2, "Category", group = TRUE)
print(tukey_result_2)


# Plot result
group_data <- tukey_result_2$groups[order(rownames(tukey_result_2$groups)),]
my_plot <- ggplot(df_cf_rep2, aes(x = Category, y = PD)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data),
                y = max(df_cf_rep1$PD) + 1,
                label = group_data$groups),
            col = 'black',
            size = 10) +
  geom_boxplot()+
  theme(text=element_text(size = 20),
        axis.title.x = element_blank())
my_plot

#Tested with Kruskal-Walis test for just nest-in (data not normally distributed)
fnm_nest_in_rep2 <- subset_samples(fnm_cf_rep2, Category == "nest-in")
sample_names(fnm_nest_in_rep2)
df_fnm_nest_in_rep2 <- as(sample_data(fnm_nest_in_rep2), "data.frame")
qqp((df_fnm_nest_in_rep2$SR^2), "norm")
shapiro.test((df_fnm_nest_in_rep2$SR)^2)
shapiro.test((df_fnm_nest_in_rep2$PD)^2)


kruskal.test(SR^2 ~ Nest_no, data = df_fnm_nest_in_rep2)
kruskal.test(PD^2 ~ Nest_no, data = df_fnm_nest_in_rep2)

#Nest-out
fnm_nest_out_rep2 <- subset_samples(fnm_cf_rep2, Category == "nest-out")
sample_names(fnm_nest_out_rep2)
df_fnm_nest_out_rep2 <- as(sample_data(fnm_nest_out_rep2), "data.frame")
shapiro.test((df_fnm_nest_out_rep2$SR)^2)
shapiro.test((df_fnm_nest_out_rep2$PD)^2)


kruskal.test(SR^2 ~ Nest_no, data = df_fnm_nest_out_rep2)
kruskal.test(PD^2 ~ Nest_no, data = df_fnm_nest_out_rep2)

