#R Code for stacked bar plot figures
#Sarah McGrath-Blaser
#August 2020
#Run code in "Frog_microbiome_analysis_condensed.R" before running this

################### Creating Stacked Bar Plot Figure ###############################################
library(ggplot2)

#Taxanomic stacked bar plots 
#Get subset of samples with nest 2 removed but all other replicates
'%ni%' <- Negate('%in%')
fnm_taxa_plots <- subset_samples(fnm_r_final, X.SampleID %ni% c ('AFCN2','AFDN2', 'AFVN2',
                                                                 'AMC','AMD', 'AMV','N2-1', 
                                                                 'N2-2', 'N2-3', 'SterileSwab',
                                                                 'AMSecretion','NAR2',
                                                                 'H2Ocontrol1','H2Ocontrol2'
                                                                 ,'H2Ocontrol3','H20controlN4'
                                                                 ,'H2ON3','N1H2OFTads'))

sample_names(fnm_taxa_plots)
#Interested in arranging via Phylum
#Look at number of features per Phylum
table(tax_table(fnm_taxa_plots)[, "Phylum"], exclude = NULL)
#We see that there is a Phylum characterized as "NA"

top15ph <- sort(tapply(taxa_sums(fnm_taxa_plots), tax_table(fnm_taxa_plots)[, "Phylum"],sum), TRUE)[1:15]
top15ph

#This removes Phylum NA and any Phylum characterized by ambiguous annotation.
fnm_taxa_plots <- subset_taxa(fnm_taxa_plots, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
table(tax_table(fnm_taxa_plots)[, "Phylum"], exclude = NULL)
#Now we see that the NA is removed

top15ph <- sort(tapply(taxa_sums(fnm_taxa_plots), tax_table(fnm_taxa_plots)[, "Phylum"],sum), TRUE)[1:15]
top15ph


############
#install.packages("devtools")
#devtools::install_github('bbc/bbplot')
library(bbplot)
library(dplyr)

#https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/index.html
###Stacked bar plot by Phylum
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
###Yowza###
fnm_phylum <- fnm_taxa_plots %>% #Create new object
  tax_glom(taxrank = "Phylum", NArm=FALSE) %>% #Condense by Phylum      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts        
  arrange(Phylum) #Order Phyla

fnm_phylum <- with(fnm_phylum, fnm_phylum[order(sample_Species, Sex, Category),])
fnm_phylum
  
ggplot(fnm_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
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

ggplot(data = fnm_adult, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style()

#Check the labels to see what order the stacked bars are in.
ggplot(data = fnm_phylum[fnm_phylum$Category %in% c('adult','nest-in','nest-out','tadpole-nest','tadpole-pond'),], aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  facet_grid(~Category,scales="free") +
  theme(axis.text.x = element_text(angle = 90))
  
#High res image of the first row of the stacked barplots.
tiff("All_replicates_phylum_1.tiff", units = "mm", width = 170, height = 60, res = 300)

ggplot(data = fnm_phylum[fnm_phylum$Category %in% c('adult','nest-in','nest-out','tadpole-nest','tadpole-pond'),], aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(text = element_text(size=10),
        legend.position = "right", 
        legend.justification = "left",
        strip.text.x = element_blank(),
        legend.text=element_text(size=12)) +
  guides(color = guide_legend(override.aes = list(size = 0.1)))
dev.off()

#Rest of the samples in the second row.
tiff("All_replicates_phylum_2.tiff", units = "mm", width = 170, height = 60, res = 300)

ggplot(data = fnm_phylum[fnm_phylum$Category %in% c('leaf','water-terrarium','water-enclosure'),], aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(text = element_text(size=10),
        legend.position = "right", 
        legend.justification = "left",
        strip.text.x = element_blank(),
        legend.text=element_text(size=12)) +
  guides(color = guide_legend(override.aes = list(size = 0.1)))
dev.off()

#https://stackoverflow.com/questions/45611725/multiple-rows-in-facet-grid
#install.packages("gridExtra")
library(gridExtra)
############################# All Replicates Phylum Level #############################################
#Make each plot an object to use later in grid.arrange
tiff("All_replicates_phylum_1.tiff", units = "mm", width = 170, height = 200, res = 300)

phylum1 <- ggplot(data = fnm_phylum[fnm_phylum$Category %in% c('adult','leaf','nest-in','water-terrarium'),], aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(text = element_text(size=20),
        legend.position = "right", 
        legend.justification = "left",
        strip.text.x = element_blank(),
        legend.text=element_text(size=16)) +
  guides(color = guide_legend(override.aes = list(size = 0.1)))
phylum1
#function for getting the legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Make legend an object
legend_phylum <- get_legend(phylum1)

#Second plot
phylum2 <- ggplot(data = fnm_phylum[fnm_phylum$Category %in% 
                                      c('tadpole-nest','nest-out',
                                        'tadpole-pond','water-enclosure'),], 
                  aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category, scale="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "none",
        strip.text.x = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

#Remove legend from first plot
phylum1 <- phylum1 + theme(legend.position="none")

#Arrange in a nice fashion
#Note: order of objects will change plot and widths refers to width of items in columns.
grid.arrange(phylum1, legend_phylum, phylum2, ncol=2,nrow=2,layout_matrix = rbind(c(1,2), c(3,2)),
             widths = c(4,1))
dev.off()
############################# All replicates Genus #############################################################
#Genus taxonomic rank
#Interested in arranging via Genus
#Look at number of features per Genus
table(tax_table(fnm_taxa_plots)[, "Genus"], exclude = NULL)
#We see that there is a Phylum characterized as "NA"

#This removes Genus NA and any Genus characterized by ambiguous annotation.
fnm_taxa_plots <- subset_taxa(fnm_taxa_plots, !is.na(Genus) & !Genus %in% c("", "uncharacterized"))
table(tax_table(fnm_taxa_plots)[, "Genus"], exclude = NULL)
#Now we see that the NA is removed

fnm_genus <- fnm_taxa_plots %>% #Create new object
  tax_glom(taxrank = "Genus") %>% #Condense by Genus      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>% #Filter out low abundance counts    
  arrange(Genus) #Order Genera

fnm_genus$Genus <- as.character(fnm_genus$Genus) 
fnm_genus$Genus[fnm_genus$Abundance < 0.08] <- "Low Abundance"

#Look at just the adults to see how they are ordered (i.e., by body location sampled and by sex)
ggplot(data=fnm_genus[fnm_genus$Category %in% ('adult'),], aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "none") +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

ggplot(data=fnm_genus[fnm_genus$Category %in% ('nest-in'),], aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "none") +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

ggplot(data=fnm_genus[fnm_genus$Category %in% ('tadpole-nest'),], aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "none") +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

ggplot(fnm_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
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

#Make each plot an object to use later in grid.arrange
genus1 <- ggplot(data = fnm_genus[fnm_genus$Category %in% c('adult','leaf','nest-in','nest-out'),], aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category, scales="free") +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(text = element_text(size=20),
        legend.position = "right", 
        legend.justification = "left",
        strip.text.x = element_blank(),
        legend.text=element_text(size=16)) +
  guides(color = guide_legend(override.aes = list(size = 0.1)))

#function for getting the legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Make legend an object
legend_genus <- get_legend(genus1)

#Second plot
genus2 <- ggplot(data = fnm_genus[fnm_genus$Category %in% 
                                    c('tadpole-nest','water-terrarium',
                                      'tadpole-pond','water-enclosure'),], 
                 aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category, scale="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "none",
        strip.text.x = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

#Remove legend from first plot
genus1 <- genus1 + theme(legend.position="none")

#Arrange in a nice fashion
#Note: order of objects will change plot and widths refers to width of items in columns.
grid.arrange(genus1, legend_genus, genus2, ncol=2,nrow=2,layout_matrix = rbind(c(1,2), c(3,2)),
             widths = c(4,1))

#Other good tutorials
#https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#supervised_learning
#This one has really good info and visualizations!
#https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/05--plotting.html

################################Merged Plots###########################################
merged_fnm_taxa_plots = merge_samples(fnm_taxa_plots, "Category")
SD = merge_samples(sample_data(fnm_taxa_plots), "Category")

fnm_merged_phylum <- merged_fnm_taxa_plots %>% #Create new object
  tax_glom(taxrank = "Phylum") %>% #Condense by Phylum      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts        
  arrange(Phylum) #Order Phyla


ggplot(fnm_merged_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
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


#Genus taxonomic rank
fnm_merged_genus <- merged_fnm_taxa_plots %>% #Create new object
  tax_glom(taxrank = "Genus") %>% #Condense by Genus     
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts        
  arrange(Genus) #Order Genus

fnm_merged_genus$Genus <- as.character(fnm_merged_genus$Genus) 
fnm_merged_genus$Genus[fnm_merged_genus$Abundance < 0.05] <- "Low Abundance"

library(RColorBrewer)
display.brewer.all()

ggplot(fnm_merged_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=16)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

#Different colors with labels
ggplot(fnm_merged_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_fill_manual(values = colorRampPalette(brewer.pal(33,"Set3"))(33))

################################ Polypedates leucomystax ##################################################
fnm_leucomystax <- subset_samples(fnm_taxa_plots, Species == "leucomystax")
sample_names(fnm_leucomystax)

fnm_leuc_genus <- fnm_leucomystax %>% #Create new object
  tax_glom(taxrank = "Genus") %>% #Condense by Genus      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() #%>% #Reformat to long form instead of wide
  
fnm_leuc_genus$Genus <- as.character(fnm_leuc_genus$Genus) 
fnm_leuc_genus$Genus[fnm_leuc_genus$Abundance < 0.05] <- "Low Abundance"
fnm_leuc_genus

#library(scales)
#names(palette1) = fnm_leuc_genus$Genus
palette1 = c()
#print(palette1)
palette1['Low Abundance'] = 'grey'
palette1['Acinetobacter'] = 'yellowgreen'
palette1['Agrobacterium'] = 'yellow2'
palette1['Bacteroides'] = 'steelblue2'
palette1['Bdellovibrio'] = 'springgreen2'
palette1['Brenneria'] = 'violetred3'
palette1['Desulfovibrio'] = 'slategray2'
palette1['Devosia'] = 'yellow4'
palette1['Faecalibacterium'] = 'tan2'
palette1['Flavobacterium'] = 'thistle'
palette1['Janthinobacterium'] = 'sienna2'
palette1['Lactobacillus'] = 'turquoise3'
palette1['Methylobacillus'] = 'wheat3'
palette1['Methylotenera'] = 'tomato2'
palette1['Prevotella'] = 'slateblue3'
palette1['Pseudomonas'] = 'forestgreen'
palette1['PW3'] = 'deepskyblue4'
palette1['Sphingobacterium'] = 'orange2'
palette1['Staphylococcus'] = 'palevioletred3'
palette1['Streptococcus'] = 'goldenrod1'
palette1['Vogesella'] = 'darkseagreen2'
palette1['Achromobacter'] = 'darkkhaki'
palette1['Bordetella'] = 'red3'
palette1['Cetobacterium'] = 'darkslategray'
palette1['Chryseobacterium'] = 'lemonchiffon1'
palette1['Herbaspirillum'] = 'honeydew2'
palette1['Rhizobium'] = 'lightcoral'
palette1['Stenotrophomonas'] = 'khaki'
palette1['Aquitalea'] = 'ivory4'
palette1['Comamonas'] = 'mediumorchid3'
palette1['Elstera'] = 'lightsalmon3'
palette1['Enterobacter'] = 'lightseagreen'
palette1['Haloferula'] = 'lightsteelblue'
palette1['Novosphingobium'] = 'mistyrose3'
palette1['Parabacteroides'] = 'plum4'
#print(palette1)

tiff("P_leuc_genera.tiff", units = "mm", width = 170, height = 100, res = 300)

leuc <- ggplot(fnm_leuc_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_fill_manual(values=palette1) +
  facet_grid(~Category,scales="free") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=12),
        axis.text.x=element_blank())
leuc
dev.off()


################################### Polypedates macrotis ################################################
fnm_macrotis <- subset_samples(fnm_taxa_plots, Species == "macrotis")

fnm_mac_genus <- fnm_macrotis %>% #Create new object
  tax_glom(taxrank = "Genus") %>% #Condense by Genus      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() #%>% #Reformat to long form instead of wide

fnm_mac_genus$Genus <- as.character(fnm_mac_genus$Genus) 
fnm_mac_genus$Genus[fnm_mac_genus$Abundance < 0.05] <- "Low Abundance"

#set color palette to accommodate the number of genera
colourCount = length(unique(fnm_mac_genus$Genus))
getPalette = colorRampPalette(brewer.pal(9, "Set3"))

tiff("P_mac_genera.tiff", units = "mm", width = 170, height = 100, res = 300)
ggplot(fnm_mac_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_fill_manual(values=palette1) +
  facet_grid(~Category,scales="free") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=12),
        axis.text.x=element_blank())
dev.off()
################################# Polypedates otilophus ####################################################
fnm_otil <- subset_samples(fnm_taxa_plots, Species == "otilophus")

fnm_otil_genus <- fnm_otil %>% #Create new object
  tax_glom(taxrank = "Genus") %>% #Condense by Genus      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() #%>% #Reformat to long form instead of wide

fnm_otil_genus$Genus <- as.character(fnm_otil_genus$Genus) 
fnm_otil_genus$Genus[fnm_otil_genus$Abundance < 0.05] <- "Low Abundance"

#set color palette to accommodate the number of genera
#colourCount = length(unique(fnm_otil_genus$Genus))
#getPalette = colorRampPalette(brewer.pal(9, "Set3"))

tiff("P_otil_genera.tiff", units = "mm", width = 170, height = 100, res = 300)

ggplot(fnm_otil_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_fill_manual(values=palette1) +
  facet_grid(~Category,scales="free") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=12),
        axis.text.x=element_blank())

dev.off()


ggplot(data=fnm_otil_genus[fnm_otil_genus$Category %in% ('adult'),], aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "none") +
  guides(color = guide_legend(override.aes = list(size = 0.2)))







