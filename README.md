### Hello, this is a standalone (local) version of siSPOTr Appliaction by the Davidson Lab ###

### You can also access it online: https://davidsonlabchop.shinyapps.io/sispotrio/ ####

### Simply download this repository, and the libraries below and run it using RStudio (or your favorite R code editor). ###

### Prerequisites/Libraries Used:
R-4.4.0+
library(shiny)
library(DT)
library(shinythemes)
library(Biostrings)
library(data.table)
library(stringr)
library(shinyjs)
library(writexl)
library(shinyalert)

### All packages can be downloaded by running:
install.packages("PACKAGE_NAME") 
### Biostrings can be installed using bioconductor BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")

### Usage:
Simply open siSPOTR_basic_V3.R and (if running in RStudio) press RunApp, which will open the GUI 

### Input:
simple fasta/text file with fasta formating of your RNA guides/gene of interest (an example input file "test_seqs.fa" is in the repository)

For Example: 
>RNA guide 1
UGCAGGUA....AAGCACAAG
>RNA guide 2
UCACAAG.....CCCGACUUA


