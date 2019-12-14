#'---
#' title: "Comparison of germline SNP calls between Line 6, Line7, and Arkansas lines"
#' author: Alec Steep
#' date: "`r format.Date( Sys.Date(), '%Y%m%d' )`"
#' output: 
#'     html_document:
#'         code_folding: hide
#'         toc: true
#'         highlight: zenburn
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(cache = FALSE)

#' ## Setup the Environment

#+ Setup Environment

################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Users/Alec/Documents/Bioinformatics/germplasm_collab/20191205_pilot-snp-compare_steep'
setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
BiocManager::install("ensemblVEP")
#install.packages("BSgenome.Ggallus.UCSC.galGal5")

# Load dependencies
pacs...man <- c("tidyverse","data.table","R.utils","ggpubr","plyr","readxl", "GenomicRanges","BSgenome","rtracklayer","ensemblVEP")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) 
})

############################################################
##### Functions ############################################
############################################################

# Make the 'not in' operator
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}

# Capture the Date
date <- format.Date( Sys.Date(), '%Y%m%d' )
auth <- "steep"

## explicit gc, then execute `expr` `n` times w/o explicit gc, return timings (~"Jenny Bryan, updating work of Winston Chang")
benchmark <- function(n = 1, expr, envir = parent.frame()) {
        expr <- substitute(expr)
        gc()
        map(seq_len(n), ~ system.time(eval(expr, envir), gcFirst = FALSE))
}

#' ## Load & Clean Data
#' ##### Data Files to Load:
#' * Arkansas SNP Loci
#' * Line 6, Line 7, and 6x7 F1 Genotype calls from 600K SNP Arrays
#' 
#'
#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# Line 6&7 and ARK SNPs
######################

# Load in the data
L67_ARK <- read.table(file = paste0('./data/',date,'_L67-ARK-df_',auth,'.txt'), 
                      header = TRUE, sep ='\t')
# As tibble
L67_ARK <- as_tibble(L67_ARK)

# Double check structure
str(L67_ARK)
summary(L67_ARK)


