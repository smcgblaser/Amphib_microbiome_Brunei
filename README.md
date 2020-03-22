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
