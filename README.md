### Hello, this is a standalone version of siSPOTr Appliaction by the Davidson Lab ###

### No instalation necesary. Simply download this repository and run it using RStudio (or your favorite R code editor). ###

### Prerequisites/Libraries Used:
R
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

### Usage
### Simply open siSPOTR_basic_V3.R and (if running in RStudio) press RunApp 
