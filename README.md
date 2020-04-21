# Amphib_microbiome_Brunei
Frog and caecilian microbiome data from Brunei, Borneo.  Collected and analyzed as part of my Master's work with James Madison University.

Miniconda and Ubuntu terminal already installed.

Downloaded latest version of QIIME2 into Ubuntu terminal.
mcgratse@DESKTOP-0PO1GR1:~$ wget https://data.qiime2.org/distro/core/qiime2-2020.2-py36-linux-conda.yml 

Created new environment to run latest QIIME2 version.
mcgratse@DESKTOP-0PO1GR1:~$ conda env create -n qiime2-2020.2 --file qiime2-2020.2-py36-linux-conda.yml 

Activated my new environment.
mcgratse@DESKTOP-0PO1GR1:~$ source activate qiime2-2020.2  

Made new directory for all working files within the pipeline.
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~$ mkdir amphib_micro_brunei   
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~$ cd amphib_micro_brunei/ 

Import sequencing files (from Illumina sequencing runs) into directory for visualization.

#Had to copy folder with sequences from Windows 10 to Ubuntu terminal
#Files can be accessed in Ubuntu command line from /mnt/c/Users/Sarah (this is where items are more perminently "mounted" to Ubuntu).
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ cp -r /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Frog_Microbiome_Run_1/Frog-45501458 ~/amphib_micro_brunei 

#Changed name of directory to make sense for when we bring in the second sequencing run later on.
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ mv Frog-45501458/ Frog_Microbiome_Run_1/                
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ ls                                                      
Frog_Microbiome_Run_1    

#Imported sequence data
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Frog_Microbiome_Run_1/ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end-frog1.qza       Imported Frog_Microbiome_Run_1/ as CasavaOneEightSingleLanePerSampleDirFmt to demux-paired-end-frog1.qza  

#Imported Illumina run 2 sequence data
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ cp -r /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Frog_Microbiome_Run_2/Frog_2/ ~/amphib_micro_brunei   

#Changed name
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ mv Frog_2/ Frog_Microbiome_Run_2/  

#Summarize demux data to view quality of sequence reads
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime demux summarize --i-data demux-paired-end-frog1.qza --o-visualization  demux-paird-end-frog1.qzv  Saved Visualization to: demux-paired-end-frog1.qzv              

(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime demux summarize --i-data demux-paired-end-frog2.qza --o-visualization  demux-paird-end-frog2.qzv  Saved Visualization to: demux-paired-end-frog2.qzv 

#Copied visualization files to folder on Windows system.
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ cp -r ~/amphib_micro_brunei/demux-paired-end-frog1.qzv /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Visualizations/                                                          (qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ cp -r ~/amphib_micro_brunei/demux-paired-end-frog2.qzv /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Visualizations/   

#Have to open view.qiime2.org then drag and drop files to view.
#Looked at interactive quality plots for quality scores.

#Lots more sequence counts with run 2 but better quality with run 1.

#Link to forum post I created asking about plot interpretation. Helps to understand how to understand these plots.
https://forum.qiime2.org/t/interactive-quality-plot-interpretation-and-colors/1843

#Imported the forward reads from Illumina runs Frog1 and Frog2 (did two runs because got really low coverage from the first run).
#Had to copy all R1 files from FrogRun1 and FrogRun2 into a new folder.  Renamed each file with R1 or R2 at the end of each SampleID name to identify which run each file came from.
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path ~/amphib_micro_brunei/Frog_Microbiome_Forward_Reads/ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-frog-forward-reads.qza    
Imported Frog_Microbiome_Forward_Reads/ as CasavaOneEightSingleLanePerSampleDirFmt to demux-frog-forward-reads.qza

(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime demux summarize --i-data demux-frog-forward-reads.qza --o-visualization demux-frog-forward-reads.qzv  Saved Visualization to: demux-frog-forward-reads.qzv  

(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ cp -r ~/amphib_micro_brunei/demux-frog-forward-reads.qzv /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Visualizations/  

##Based on this plot of my forward reads, I would like to trim to a length where >50% of my reads fall below q 20, therefore I will trim at sequence base 220.
##Forward primer (515F) is 19 nucleotides long so will trim 19 from the front (left).
##Next I filtered sequences based on quality score. I used the default command options provided via QIIME2 (minimum PHRED score =4, maximum number of consecutive low PHRED =3, maximum number of ambiguous (i.e., N) base calls =0).
##Link to code information https://docs.qiime2.org/2018.2/plugins/available/quality-filter/q-score/

##Link to documentation regarding these numbers https://qiita.ucsd.edu/static/doc/html/deblur_quality.html
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime quality-filter q-score --i-demux demux-frog-forward-reads.qza --o-filtered-sequences demux-filtered-frog-forward-reads.qza --o-filter-stats demux-filter-stats-frog-forward-reads.qza Saved SampleData[SequencesWithQuality] to: demux-filtered-frog-forward-reads.qza    Saved QualityFilterStats to: demux-filter-stats-frog-forward-reads.qza   

(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime deblur denoise-16S --i-demultiplexed-seqs demux-filtered-frog-forward-reads.qza --p-trim-length 220 --p-left-trim-len 19 --o-representative-sequences rep-seqs-deblur-frog-forward-reads.qza --p-sample-stats --o-stats deblur-frog-forward-read-stats.qza --o-table table-deblur-frog-forward-reads.qza                                  Saved FeatureTable[Frequency] to: table-deblur-frog-forward-reads.qza               
Saved FeatureData[Sequence] to: rep-seqs-deblur-frog-forward-reads.qza              
Saved DeblurStats to: deblur-frog-forward-read-stats.qza 

(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime metadata tabulate --m-input-file demux-filter-stats-frog-forward-reads.qza --o-visualization demux-filter-stats-frog-forward-reads.qzv                                                
Saved Visualization to: demux-filter-stats-frog-forward-reads.qzv                   
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime deblur visualize-stats --i-deblur-stats deblur-frog-forward-read-stats.qza --o-visualization deblur-frog-forward-read-stats.qzv                                                       
Saved Visualization to: deblur-frog-forward-read-stats.qzv  

(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime feature-table summarize --i-table table-deblur-frog-forward-reads.qza --o-visualization table-deblur-frog-forward-reads.qzv --m-sample-metadata-file '/mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Metadata/Frog_Micro_Metadata_2020.2 - Sheet1.tsv'       
Saved Visualization to: table-deblur-frog-forward-reads.qzv  
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ cp -r ~/amphib_micro_brunei/table-deblur-frog-forward-reads.qzv /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Visualizations/  

##Sequence count per sample ranges from 22958 to 2078.  Used the interactive sample detail to determine sampling depth (3,500).
##The rep-seqs visualization file allows you to interactively BLAST against the NCBI database.
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime feature-table tabulate-seqs --i-data rep-seqs-deblur-frog-forward-reads.qza --o-visualization rep-seqs-deblur-frog-forward-reads.qzv                                                  
Saved Visualization to: rep-seqs-deblur-frog-forward-reads.qzv                      
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ cp -r ~/amphib_micro_brunei/rep-seqs-deblur-frog-forward-reads.qzv /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Visualizations/ 


##Generating phylogenetic trees for diversity analyses.

##Multiple sequence alignment via mafft.
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime alignment mafft --i-sequences rep-seqs-deblur-frog-forward-reads.qza --o-alignment aligned-rep-seqs-frog-forward-reads.qza                                                            
Saved FeatureData[AlignedSequence] to: aligned-rep-seqs-frog-forward-reads.qza 

##Masked (filtered) alignment to remove highly variable sections.
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime alignment mask --i-alignment aligned-rep-seqs-frog-forward-reads.qza --o-masked-alignment masked-aligned-rep-seqs-frog-forward-reads.qza                                              
Saved FeatureData[AlignedSequence] to: masked-aligned-rep-seqs-frog-forward-reads.qza        

##Creating an unrooted tree via fasttree.
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs-frog-forward-reads.qza --o-tree unrooted-tree-frog-forward-reads.qza                                                         
Saved Phylogeny[Unrooted] to: unrooted-tree-frog-forward-reads.qza 

##Now rooting our tree to use in analyses (using a midpoint root).
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime phylogeny midpoint-root --i-tree unrooted-tree-frog-forward-reads.qza --o-rooted-tree rooted-tree-frog-forward-reads.qza                                                              
Saved Phylogeny[Rooted] to: rooted-tree-frog-forward-reads.qza  

##Before continuing with examination of alpha and beta diversity metrics, I want to visualize the alpha rarefaction plots to see if I captured the total diversity of my samples.
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime diversity alpha-rarefaction --i-table table-deblur-frog-forward-reads.qza --i-phylogeny rooted-tree-frog-forward-reads.qza --p-max-depth 4500 --m-metadata-file '/mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Metadata/Frog_Micro_Metadata_2020.2 - Sheet1.tsv' --o-visualization alpha-rarefaction-frog-forward-reads.qzv                     
Saved Visualization to: alpha-rarefaction-frog-forward-reads.qzv                    
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ cp -r ~/amphib_micro_brunei/alpha-rarefaction-frog-forward-reads.qzv /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Visualizations/    

##This visualization produces two plots.  Here is the information about what these plots mean.
“The visualization will have two plots. The top plot is an alpha rarefaction plot, and is primarily used to determine if the richness of the samples has been fully observed or sequenced. If the lines in the plot appear to “level out” (i.e., approach a slope of zero) at some sampling depth along the x-axis, that suggests that collecting additional sequences beyond that sampling depth would not be likely to result in the observation of additional features. If the lines in a plot don’t level out, this may be because the richness of the samples hasn’t been fully observed yet (because too few sequences were collected), or it could be an indicator that a lot of sequencing error remains in the data (which is being mistaken for novel diversity).

The bottom plot in this visualization is important when grouping samples by metadata. It illustrates the number of samples that remain in each group when the feature table is rarefied to each sampling depth. If a given sampling depth d is larger than the total frequency of a sample s (i.e., the number of sequences that were obtained for sample s), it is not possible to compute the diversity metric for sample s at sampling depth d. If many of the samples in a group have lower total frequencies than d, the average diversity presented for that group at d in the top plot will be unreliable because it will have been computed on relatively few samples. When grouping samples by metadata, it is therefore essential to look at the bottom plot to ensure that the data presented in the top plot is reliable.”

##Based on these plots I feel that the sequencing depth of 4500 captures the total diversity in the samples and increasing my sampling depth would not provide more accurate data for diversity calculations.

##Next I assigned taxonomy using a pre-trained classifier (Native-Bayes trained on Greengenes database 515F/806R region) downloaded from here…
https://docs.qiime2.org/2018.2/data-resources/#taxonomy-classifiers-for-use-with-q2-feature-classifier

(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ wget \              >   -O "gg-13-8-99-515-806-nb-classifier.qza" \                                     >   "https://data.qiime2.org/2020.2/common/gg-13-8-99-515-806-nb-classifier.qza"    Will not apply HSTS. The HSTS database must be a regular and non-world-writable file.                                                                                   ERROR: could not open HSTS store at '/home/mcgratse/.wget-hsts'. HSTS will be disabled.                                                                                 --2020-04-03 10:43:26--  https://data.qiime2.org/2020.2/common/gg-13-8-99-515-806-nb-classifier.qza                                                                     Resolving data.qiime2.org (data.qiime2.org)... 52.35.38.247                         Connecting to data.qiime2.org (data.qiime2.org)|52.35.38.247|:443... connected.     HTTP request sent, awaiting response... 302 FOUND                                   Location: https://s3-us-west-2.amazonaws.com/qiime2-data/2020.2/common/gg-13-8-99-515-806-nb-classifier.qza [following]                                                 --2020-04-03 10:43:27--  https://s3-us-west-2.amazonaws.com/qiime2-data/2020.2/common/gg-13-8-99-515-806-nb-classifier.qza                                              Resolving s3-us-west-2.amazonaws.com (s3-us-west-2.amazonaws.com)... 52.218.222.24  Connecting to s3-us-west-2.amazonaws.com (s3-us-west-2.amazonaws.com)|52.218.222.24|:443... connected.                                                                  HTTP request sent, awaiting response... 200 OK                                      Length: 28373581 (27M) [application/x-www-form-urlencoded]                          Saving to: ‘gg-13-8-99-515-806-nb-classifier.qza’                                                                                                                       gg-13-8-99-515-806-n 100%[======================>]  27.06M  10.6MB/s    in 2.6s                                                                                         2020-04-03 10:43:30 (10.6 MB/s) - ‘gg-13-8-99-515-806-nb-classifier.qza’ saved [28373581/28373581]    

(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime feature-classifier classify-sklearn --i-classifier gg-13-8-99-515-806-nb-classifier.qza --i-reads rep-seqs-deblur-frog-forward-reads.qza --o-classification taxonomy-frog-forward-reads.qza                Saved FeatureData[Taxonomy] to: taxonomy-frog-forward-reads.qza
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime metadata tabulate --m-input-file taxonomy-frog-forward-reads.qza --o-visualization taxonomy-frog-forward-reads.qzv        
Saved Visualization to: taxonomy-frog-forward-reads.qzv         
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ cp -r ~/amphib_micro_brunei/taxonomy-frog-forward-reads.qzv /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Visualizations/           


#Make grouped table.
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime feature-table group --i-table table-deblur-frog-forward-reads.qza --p-axis sample --m-metadata-file '/mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Metadata/Frog_Micro_Metadata_2020.2 - Sheet1 (1).tsv' --m-metadata-column True_SampleID --p-mode sum --o-grouped-table grouped-table-deblur-frog-forward-reads.qza                                                       
Saved FeatureTable[Frequency] to: grouped-table-deblur-frog-forward-reads.qza                                                   
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime feature-table summarize --i-table grouped-table-deblur-frog-forward-reads.qza --o-visualization grouped-table-deblur-frog-forward-reads.qzv                                               
Saved Visualization to: grouped-table-deblur-frog-forward-reads.qzv         

#Exported files to use in R
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime tools export --input-path table-deblur-frog-forward-reads.qza --output-path /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-table-deblur-frog-forward-reads 
Exported table-deblur-frog-forward-reads.qza as BIOMV210DirFmt to directory /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-table-deblur-frog-forward-reads                                                         
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime tools export --input-path rooted-tree-frog-forward-reads.qza --output-path /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-rooted-tree-frog-forward-reads   
Exported rooted-tree-frog-forward-reads.qza as NewickDirectoryFormat to directory /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-rooted-tree-frog-forward-reads                                                    
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ qiime tools export --input-path taxonomy-frog-forward-reads.qza --output-path /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-taxonomy-frog-forward-reads      Exported taxonomy-frog-forward-reads.qza as TSVTaxonomyDirectoryFormat to directory /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-taxonomy-frog-forward-reads    

(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ cp -r /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-taxonomy-frog-forward-reads/taxonomy.tsv /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-taxonomy-frog-forward-reads/biom-taxonomy.tsv   

#

#I do not think this is correct
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ biom add-metadata -i /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-table-deblur-frog-forward-reads/feature-table.biom -o /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-table-deblur-frog-forward-reads/table-with-taxonomy.biom --observation-metadata-fp '/mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/Metadata/Frog_Micro_Metadata_2020.2 - Sheet1 (1).tsv' --sc-separated taxonomy    

#Re-ran the code as this...
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ biom add-metadata -i /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-table-deblur-frog-forward-reads/feature-table.biom -o /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-table-deblur-frog-forward-reads/table-with-taxonomy.biom --observation-metadata-fp '/mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-taxonomy-frog-forward-reads/biom-taxonomy.tsv' --sc-separated taxonomy     

#Still was not working so found some more information on this forum...https://forum.qiime2.org/t/biom-add-metadata-problem/8302/15
#Before running this, I went to the biom-taxonomy.tsv file and changed the header to #OTUID (tab) taxonomy (tab) confidence.
#Then ran this code.
(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ biom add-metadata -i /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-table-deblur-frog-forward-reads/feature-table.biom -o /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-table-deblur-frog-forward-reads/table-with-taxonomy.biom --observation-metadata-fp '/mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-taxonomy-frog-forward-reads/biom-taxonomy.tsv' --sc-separated taxonomy                               

(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ biom convert -i /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-table-deblur-frog-forward-reads/table-with-taxonomy.biom  -o /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-table-deblur-frog-forward-reads/table-with-taxonomy.tsv --to-tsv --header-key taxonomy    

(qiime2-2020.2) mcgratse@DESKTOP-0PO1GR1:~/amphib_micro_brunei$ biom head -i /mnt/c/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/R_files/Exported_files/exported-table-deblur-frog-forward-reads/table-with-taxonomy.tsv                
# Constructed from biom file                                                                                            
#OTU ID AFCN1R1 AFCN1R2 AFCN2R1 AFCN2R2 AFCN3R1                                                                         237c27a2fe9ca9dabc0df88201ab7a87        1374.0  3486.0  0.0     0.0     0.0                                             e79a0e75fb6d179a719843ec80071e27        891.0   2222.0  0.0     0.0     0.0                                             82587105b8072b9fe258652b741704a8        207.0   564.0   5.0     25.0    4.0                                             64a6299a3845b8876341ea645bbed657        150.0   447.0   0.0     2.0     307.0                                           c3a1660be67fd87761a29eac85dbf344        136.0   367.0   0.0     0.0     2.0   

#Checked the table-with-taxonomy.tsv file in Notepad and it does have a taxonomy column added at the end of the column list.

#Brought it into R.  This fucking worked.  Holy fucking fuck.

