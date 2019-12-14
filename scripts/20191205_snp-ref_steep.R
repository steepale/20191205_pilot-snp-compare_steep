#'---
#' title: "Determine Gallus Gallus Reference Genome for Variant Calls"
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
# BiocManager::install("ade4")
#install.packages("ade4")

# Load dependencies
pacs...man <- c("tidyverse","data.table","R.utils","ggpubr","plyr","readxl", "GenomicRanges","BSgenome","BSgenome.Ggallus.UCSC.galGal4","BSgenome.Ggallus.UCSC.galGal5","rtracklayer")
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

# Generate a frequency of genotypes
################################################################################
count_genos <- function(x) {
        # Ensure x is a list
        if (!is.list(x))
                stop("input is not a list")
        # Collect the most common genotype per row
        common <- cat_mode(x) %>% unlist() %>% as.character()
        # Count the common genotype occurance
        row_vector <- x %>% unlist()
        sum(str_count(row_vector, common)) / (length(row_vector)-sum(str_count(row_vector,'---')))
}
################################################################################

# Calculate the most common categorical value from vector or list
################################################################################
cat_mode <- function(x) {
        uniqx <- unique(na.omit(x))
        uniqx[which.max(tabulate(match(x, uniqx)))]
}
################################################################################

# Function to load germline vcfs
################################################################################
load_normal_vcfs <- function(x) {
        # Ensure 'x' is a string
        if ( class(x) != "character" ) {
                stop("'x' must be a string", class.= FALSE)
        }
        # Assign the vcf file
        ND <- "/Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels/data/germline_snps/indv_samples/parents"
        vcf.gz <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_10K.g.vcf.gz')
        vcf <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_10K.g.vcf')
        # Gunzip the file
        gunzip(vcf.gz)
        # Load vcf data
        PL <- read.table(file = vcf, comment.char = "", check.names = FALSE, header = TRUE, sep = '\t')
        PL <- as_tibble(PL)
        # Adjust col name
        colnames(PL)[1] <- 'CHROM'
        SAMPLE <- colnames(PL)[10]
        # gzip the file
        gzip(vcf)
        
        # Extract the genotype information from the last column
        PL <- PL %>% separate(SAMPLE, c("GT", NA), sep = ':')
        # Select the columns of interest
        PL <- PL %>% dplyr::select(CHROM, POS, REF, ALT,GT)
        
        # Rename the allele columns
        colnames(PL)[5] <- paste0("GT")
        
        # Return the output
        as_tibble(PL)
}
################################################################################

# Function to load germline vcfs (total)
################################################################################
load_total_vcf <- function(x) {
        # Ensure 'x' is a string
        if ( class(x) != "character" ) {
                stop("'x' must be a string", class.= FALSE)
        }
        # Assign the vcf file
        ND <- "/Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels/data/germline_snps/indv_samples/parents"
        vcf.gz <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_10K.g.vcf.gz')
        vcf <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_10K.g.vcf')
        # Gunzip the file
        gunzip(vcf.gz)
        # Load vcf data
        PL <- read.table(file = vcf, comment.char = "", check.names = FALSE, header = TRUE, sep = '\t')
        PL <- as_tibble(PL)
        # Adjust col name
        colnames(PL)[1] <- 'CHROM'
        SAMPLE <- colnames(PL)[10]
        # gzip the file
        gzip(vcf)
        # Return the output
        as_tibble(PL)
}
################################################################################

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

# Chain file for liftover
######################
# Install a chain file for liftover
# Chain format: https://genome.ucsc.edu/goldenpath/help/chain.html
# Liftover file downloaded (20191205) from: http://hgdownload.soe.ucsc.edu/goldenPath/galGal4/liftOver/

# Arkansas SNPs
######################

# Load the Arkansas SNP file
ark_file <- "./data/20191205_arkansas-snps-meta_unknown.xlsx"
# xls files
ark_df <- read_excel(ark_file, sheet = "ARKANSAS with SNPs")
# Select columns of interest
ark_df <- ark_df %>% dplyr::select("Chromosome #","Strand Direction","MID","Ref Pos","Type","Ref Base","Called Base","Genotype")
# Adjust names
names(ark_df) <- c("CHROM","STRAND","LINE","POS","TYPE","REF","ALT","GENOTYPE")
# Select for SNPs & Remove Chromosomes with NA values
ark_df <- ark_df %>% filter(TYPE == 'SNP' & !is.na(CHROM))
# Adjust the strands to have proper annotation
ark_df <- ark_df %>% mutate(STRAND = case_when(STRAND == "forward strand" ~ "+", 
                                               STRAND == "reverse strand" ~ "-"))

# Adjust the genotype
table(ark_df$GENOTYPE)
ark_df <- ark_df %>% mutate(GT = case_when(GENOTYPE == "Reference" ~ "0/0", 
                                           GENOTYPE == "Hetero. Ref." ~ "0/1",
                                           GENOTYPE == "Homo. Variant" ~ "1/1",
                                           GENOTYPE == "Hetero. Not Ref." ~ "1/2"))

# Perform one-hot encoding of the LINE categorical variable: 
# https://stackoverflow.com/questions/52539750/r-how-to-one-hot-encoding-a-single-column-while-keep-other-columns-still
ark_df <- ark_df  %>% mutate(value = 1)  %>% spread(LINE, value,  fill = 0 )

# Remove the GENOTYPE and TYPE columns with GT
ark_df <- ark_df %>% dplyr::select(-GENOTYPE, -TYPE)

# Reorder the columns
ark_df <- ark_df %>% dplyr::select(CHROM, POS, REF, ALT, GT, REG, PRO, STRAND)

# Ensure unique
ark_df <- ark_df %>% unique()

# Examine data
dim(ark_df)
head(ark_df)
str(ark_df)
summary(ark_df)

# Load in Line 6 and Line 7 Data
######################
L6_df <- load_normal_vcfs('002683_Line-6')
L7_df <- load_normal_vcfs('002684_Line-7')

# Rename the all columns
L6_df$LINE <- "LINE_6"
L7_df$LINE <- "LINE_7"

# Perform a full-join
L6L7 <- full_join(L6_df, L7_df, by = c("CHROM","POS","REF","ALT","GT","LINE"), copy = TRUE, suffix = c(".x", ".y"))

# Remove heavy memory objects
rm(L6_df)
rm(L7_df)

# Perform one-hot encoding of the LINE categorical variable: 
L6L7 <- L6L7 %>% mutate(value = 1)  %>% spread(LINE, value,  fill = 0 )

# Ensure unique
L6L7 <- L6L7 %>% unique()

# Examine data
dim(L6L7)
head(L6L7)
str(L6L7)
summary(L6L7)


#' ## Determine the reference build
#'
#+ Determine reference

################################################################################
#####     Determine the Reference Build      ###################################
################################################################################

# ARKANSAS
##############

# Obtain the reference sequences (getSeq is vectorized)
ark_df$GG4_REF <- getSeq(BSgenome.Ggallus.UCSC.galGal4::Ggallus, names = paste0('chr',ark_df$CHROM), start=ark_df$POS, end=ark_df$POS,
                         strand=ark_df$STRAND, as.character=TRUE)

# Make sure all reference positions are true 
table(ark_df$REF == ark_df$GG4_REF)
# FALSE results from SNPs annotated as reverse strand (STRAND might refer to gene?)

# Realization, all of the Arkansas SNPs REF values are in relation to the positive strand

# Obtain the reference sequences (getSeq is vectorized)
ark_df$GG4_REF <- getSeq(BSgenome.Ggallus.UCSC.galGal4::Ggallus, names = paste0('chr',ark_df$CHROM), start=ark_df$POS, end=ark_df$POS,
                         strand="+", as.character=TRUE)
ark_df$GG5_REF <- getSeq(BSgenome.Ggallus.UCSC.galGal5::Ggallus, names = paste0('chr',ark_df$CHROM), start=ark_df$POS, end=ark_df$POS,
                         strand="+", as.character=TRUE)

# Make sure all reference positions are true 
table(ark_df$REF == ark_df$GG4_REF)
table(ark_df$REF == ark_df$GG5_REF)
# Chance the Rapper
7094/(21555+7094)

# Conclusion: Reference build is GG4 and all REFs refer to "+" strand

# Lines 6 and 7
######################

# Load the data into a GRanges object
L6L7_gr <- makeGRangesFromDataFrame(L6L7,
                                    keep.extra.columns = TRUE,
                                    ignore.strand=FALSE,
                                    seqnames.field='CHROM',
                                    start.field = 'POS',
                                    end.field = 'POS',
                                    strand.field = "+")

# Add 'chr' to CHROM names
seqlevelsStyle(L6L7_gr) = "UCSC"

# Convert to tibble
L6L7_rb <- as_tibble(L6L7_gr)
# Drop old REF columns
L6L7_rb <- L6L7_rb %>% dplyr::select(-end, -strand, -width)
# Adjust names
names(L6L7_rb)[1:2] <- c('CHROM','POS')

# Filter out any chromosome names beginning with "NT_", and custom adjust others
L6L7_rb <- L6L7_rb %>% filter(!grepl("^NT_",CHROM))
L6L7_rb$CHROM <- L6L7_rb$CHROM %>% str_replace_all("^Z", "chrZ")
L6L7_rb$CHROM <- L6L7_rb$CHROM %>% str_replace_all("^W", "chrW")
L6L7_rb$CHROM <- L6L7_rb$CHROM %>% str_replace_all("^LGE64", "chrLGE64")

# Collect Galgal5 reference sequences for each SNP
L6L7_rb$GG5_REF <- getSeq(BSgenome.Ggallus.UCSC.galGal5::Ggallus, names = L6L7_rb$CHROM, start=L6L7_rb$POS, end=L6L7_rb$POS, strand="+", as.character=TRUE)

# Make sure all reference positions are true 
table(L6L7_rb$REF == L6L7_rb$GG5_REF)

# Remove L6 and L7 reference calls
L6L7_rb <- L6L7_rb %>% dplyr::select(-GG5_REF)

# Conclusion: All line 6 and 7 calls refer to Galgal5, remove large object
rm(L6L7_rb)

################################################################################