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
#BiocManager::install("UpSetR")
#install.packages("doSNOW")

# Load dependencies
pacs...man <- c("tidyverse","data.table","R.utils","ggpubr","plyr","readxl", "GenomicRanges","BSgenome","rtracklayer","NAM","phylogram","fields","RColorBrewer","mapplots","LEA","adegenet","stringi","ggplot2","snow","doSNOW","parallel","UpSetR")
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

# Arkansas & Michigan SNPs
######################

# Load the annotated SNP file
sys_cmd <- paste0('ls -1 ./data/*_L67RP-annotations-quantitative_steep.txt | sort -r | head -n1')
L67RP_anne_file <- system(sys_cmd, intern = TRUE)
L67RP_anne_file <- str_sub(L67RP_anne_file,2,str_length(L67RP_anne_file)) # removes 'point' from absolute path
L67RP_anne_file <- paste0(WD,L67RP_anne_file)
# xls files
L67RP_anne <- read.table(L67RP_anne_file, 
                         sep = '\t', header = TRUE, check.names = FALSE) %>% as_tibble()

## Examine data
dim(L67RP_anne)
head(L67RP_anne)
str(L67RP_anne)
summary(L67RP_anne)

# Adjust factors
L67RP_anne$LINE_6 <- as.factor(L67RP_anne$LINE_6)
L67RP_anne$LINE_7 <- as.factor(L67RP_anne$LINE_7)
L67RP_anne$PRO <- as.factor(L67RP_anne$PRO)
L67RP_anne$REG <- as.factor(L67RP_anne$REG)

# EDA: REG-spec genes
#######################################
genes <- L67RP_anne %>% 
  filter(GT %in% c('0/1', '1/2')) %>% 
  filter(PRO != REG)%>% 
  dplyr::select(SYMBOL) %>% table() %>% sort() %>% names() %>% as.character()

write.table(genes, file = './data/REG-mut-genes_steep.txt', sep = '\n', quote = FALSE, row.names = FALSE)
########################################

# Genes of Interest
############################
# Select overlap for genes of interest
# Load the Arkansas SNP file
ark_file <- "./data/20191205_arkansas-snps-meta_unknown.xlsx"
ark_file <- paste0(WD,"/data/20191205_arkansas-snps-meta_unknown.xlsx")
# xls files
ark_df <- read_excel(ark_file, sheet = "Poultry Repro Genes")
# Select columns of interest
goi_df <- ark_df %>% dplyr::select("Gene")
# Get a vector of genes of interest
GOI <- goi_df$Gene %>% unique()

########################################################################
########## Visualization of SNPs with Upset ############################
########################################################################
# Create a matrix of just the genotype data
geno <- L67RP_anne[,11:ncol(L67RP_anne)] %>% transpose() %>% as.matrix()
# Grab a transposed df to double check naming orientation
L67RP_tp <- L67RP_anne[,11:ncol(L67RP_anne)] %>% transpose()

# Extract the sample names
SAMPLE <- names(L67RP_anne[,11:ncol(L67RP_anne)])

# Make the rownames to IDs for now
L67RP_anne$ID <- row.names(L67RP_anne)

# Select the appropriate columns
L67RP_upset <- L67RP_anne %>% dplyr::select(ID,LINE_6, LINE_7, REG, PRO)

# Remove the ID column
L67RP_anne <- L67RP_anne %>% dplyr::select(-ID)

# Replace 0.5 values
L67RP_upset[L67RP_upset == 0.5] <- 1

# Visualize intersections with UpSet
upset(data.frame(L67RP_upset), nsets = 4, number.angles = 30, point.size = 4, line.size = 2,
      mainbar.y.label = "Frequency", sets.x.label = "Variants Per LINE",
      text.scale = c(2.5, 2, 1.5, 1.4, 2, 2))

################################################################################
#####     K means cluster approach #############################################
################################################################################
# Code copied/adjusted from Kevin Falk: https://www.youtube.com/watch?v=cC6lMr8LgPQ&list=PLYX_up7DciW3r3VOmzRxKAyrZBB_eHLmy&index=9&t=0s

# Create a matrix of just the genotype data
#geno = as.matrix(L67RP_anne[,11:ncol(L67RP_anne)])
# transpose the genotypes to fit genind function
#geno_trans <- data.frame(transpose(as_tibble(geno)))
# Convert a data.frame of allele data to a genind object
#genind_geno <- df2genind(geno_trans, ploidy=2,sep='\t')
#genind_geno_b <- df2genind(data.frame(geno), ploidy=2,sep='\t')
# Perform k means with different values of k (here we'll use either 2,3, or 4)
#find.clusters(genind_geno, max.n =3, n.pca = 200, scale = FALSE)
# Increase k until no increase in fit (e.g. measure of fit interpreted by decrease in BIC)
#phenogrp <- find.clusters(x, max.n=2, n.pca=200, scale=TRUE)


################################################################################
#####     LD, Fst, and genetic distance params      ############################
################################################################################
# Code copied/adjusted from Kevin Falk: https://www.youtube.com/watch?v=VbHiYvseNdU&list=PLYX_up7DciW3r3VOmzRxKAyrZBB_eHLmy&index=7

# Create a matrix of just the genotype data
#geno <- L67RP_anne[,11:ncol(L67RP_anne)] %>% transpose() %>% as.matrix()
# Grab a transposed df to double check naming orientation
#L67RP_tp <- L67RP_anne[,11:ncol(L67RP_anne)] %>% transpose()

# Extract the sample names
#samples <- names(L67RP_anne[,11:ncol(L67RP_anne)])

# Examine SNP data
#head(geno)

################################################################################
#####     PCA of Genotypes      ################################################
################################################################################
# Code copied/adjusted from Kevin Falk: https://www.youtube.com/watch?v=dyy0nWOqKHI&list=PLYX_up7DciW3r3VOmzRxKAyrZBB_eHLmy&index=6

# Create a matrix of just the genotype data
geno <- L67RP_anne[,11:ncol(L67RP_anne)] %>% 
        transpose() %>% as.matrix()
# Grab a transposed df to double check naming orientation
L67RP_tp <- L67RP_anne[,11:ncol(L67RP_anne)] %>% transpose()

# Examine SNP data
#geno[1:4,1:20]

# Make sure there are no NA values
is.na(geno) %>% sum()

# Find columns with zero variance
length(which(apply(geno, 2, var)==0))

# Remove columns with zero variance
geno_pca <- geno[ , apply(geno, 2, var) != 0]

# Examine the data
#geno_pca[1:4,1:20]

# How many variants?
dim(geno_pca)

# Calculate principle components
PCAs <- prcomp(geno_pca, scale. = TRUE)
summary(PCAs)

# Loadings: rotated data -- the centered data X the rotation matrix
PCA_loadings = PCAs$x

# Extract the sample names
SAMPLE <- names(L67RP_anne[,11:ncol(L67RP_anne)])

# Merge the PCA loadings with annotation (meta) data
PCA_loadings <- as_tibble(cbind(SAMPLE, PCA_loadings))
names(PCA_loadings)[1] <- 'SAMPLE'

# Plot the PCA plot using ggplot2
ggplot(PCA_loadings, aes(x = PC1, y = PC2, color = SAMPLE, size = 2)) +
        geom_point(alpha=1) +
        labs(x = "PC1 (45.83%)", y = "PC2 (29.93%)")

################################################################################
#####     Sample Unique Variants      ##########################################
################################################################################

# Grab variants unique only to Regressor Line
LR_anne <- L67RP_anne %>% filter(LINE_6 == 0 & LINE_7 == 0 & PRO == 0)

# Vizulaize the data
LR_UQ <- LR_anne %>% dplyr::select(CHROM,POS,REF,ALT) %>% unique()
dim(LR_UQ)[1]

# Select appropriate columns
LR_anne %>% dplyr::select(-STRAND, -IMPACT) %>% data.frame()

# Summary
summary(LR_anne)


        