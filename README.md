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


