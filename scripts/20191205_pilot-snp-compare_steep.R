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
        vcf.gz <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf.gz')
        vcf <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf')
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
        vcf.gz <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf.gz')
        vcf <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf')
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

# Perform one-hot encoding of the LINE categorical variable (adjust for homo and hetero vars): 
# https://stackoverflow.com/questions/52539750/r-how-to-one-hot-encoding-a-single-column-while-keep-other-columns-still
ark_df <- ark_df  %>% mutate(value = case_when(GT == "0/0" ~ 0, 
                                               GT == "0/1" ~ 0.5,
                                               GT == "1/1" ~ 1,
                                               GT == "1/2" ~ 1))  %>% 
        spread(LINE, value,  fill = 0 )

# To only include 0's and 1's
#ark_df <- ark_df  %>% mutate(value = 1)  %>% spread(LINE, value,  fill = 0 )


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
# L6L7 <- L6L7 %>% mutate(value = 1)  %>% spread(LINE, value,  fill = 0 )
L6L7 <- L6L7  %>% mutate(value = case_when(GT == "0/0" ~ 0, 
                                               GT == "0/1" ~ 0.5,
                                               GT == "1/1" ~ 1,
                                               GT == "1/2" ~ 1))  %>% 
        spread(LINE, value,  fill = 0 )

# Ensure unique
L6L7 <- L6L7 %>% unique()

# Examine data
dim(L6L7)
head(L6L7)
str(L6L7)
summary(L6L7)

#' ## Perform liftover from Galgal4 to Galgal5 on Arkansas SNPs
#'
#+ GG4 to GG5

################################################################################
#####     Perform Liftover      ################################################
################################################################################

# Save a copy of the GG4 coords
ark_df4 <- ark_df
#ark_df <- ark_df4

# Load the data into a GRanges object
ark_gr <- makeGRangesFromDataFrame(ark_df,
                                   keep.extra.columns = TRUE,
                                   ignore.strand=FALSE,
                                   seqnames.field='CHROM',
                                   start.field = 'POS',
                                   end.field = 'POS',
                                   strand.field = "+")
# Add 'chr' to CHROM names
seqlevelsStyle(ark_gr) = "UCSC"

# Decompress the chain file
gunzip('./data/galGal4ToGalGal5.over.chain.gz')
# Import the Chain file
GG4to5_chain <- import.chain('./data/galGal4ToGalGal5.over.chain')
# Compress the chain file
gzip('./data/galGal4ToGalGal5.over.chain')

# Perform liftover
ark_gr5 <- liftOver(ark_gr, GG4to5_chain)
ark_gr5 <- unlist(ark_gr5)
# Convert to tibble
ark_df <- as_tibble(ark_gr5)
# Drop old REF columns
ark_df <- ark_df %>% dplyr::select(-REF, -end, -strand, -width)

# Adjust names
names(ark_df)[1:2] <- c('CHROM','POS')

# Collect the GG5 References, ideally they will be the same as the GG4 references
ark_df$GG5_REF <- getSeq(BSgenome.Ggallus.UCSC.galGal5::Ggallus, names = ark_df$CHROM, start=ark_df$POS, end=ark_df$POS,strand="+", as.character=TRUE)

# Select columns and rename
ark_df <- ark_df %>% dplyr::select(CHROM, POS, GG5_REF, ALT, GT, REG, PRO)
names(ark_df)[3] <- 'REF'

# Remove duplicate values
ark_df <- ark_df %>% unique()

# Dimensions
dim(ark_df)

#' ## Perform an inner join to determine if any of the snps are shared across datasets
#'
#+ Inner join (Lines 6&7 and Arkansas)

################################################################################
#####     Perform an Inner Join      ###########################################
################################################################################

# Remove the 'chr' string to all chromosome names
ark_df$CHROM <- str_remove_all(ark_df$CHROM, 'chr')

# Perform an inner join
L67RP <- inner_join(L6L7, ark_df, by = c('CHROM','POS','REF'), copy = FALSE, suffix = c(".x", ".y"))

# Remove the large memory object
rm(L6L7)

# Adjust columns (.x, .y) *TODO: Generate dplyr-friendly function

# Generate NA values for rows that contain similar ALT and GT values

#ALT
# Iterate through the rows and combine and add columns depending on context
L67RP <- L67RP %>% mutate(ALT = ifelse( ALT.x == ALT.y, ALT.x, NA))

# Create duplicate dataframe with differing ALT values
df.x <- L67RP %>% filter(is.na(ALT))
df.y <- L67RP %>% filter(is.na(ALT))
# Adjust ALT values
df.x <- df.x %>% mutate(ALT = ALT.x)
df.y <- df.y %>% mutate(ALT = ALT.y)
# Combine the dataframes and adjust
L67RP <- L67RP %>% filter(!is.na(ALT)) %>% 
        bind_rows(df.x,df.y) %>% 
        dplyr::select(-ALT.x, -ALT.y)

#GT
# Iterate through the rows and combine and add columns depending on context
L67RP <- L67RP %>% mutate(GT = ifelse(GT.x == GT.y, GT.x, NA))
# Create duplicate dataframe with differing GT values
df.x <- L67RP %>% filter(is.na(GT))
df.y <- L67RP %>% filter(is.na(GT))
# Adjust GT values
df.x <- df.x %>% mutate(GT = GT.x)
df.y <- df.y %>% mutate(GT = GT.y)
# Combine the dataframes and adjust
L67RP <- L67RP %>% filter(!is.na(GT)) %>% 
        bind_rows(df.x,df.y) %>% 
        dplyr::select(-GT.x, -GT.y)

# Reorder the dataframe
L67RP <- L67RP %>% dplyr::select(CHROM,POS,REF,ALT,GT,LINE_6,LINE_7,REG,PRO)

# Remove any duplicate values
L67RP <- L67RP %>% unique()
dim(L67RP)

# Save the data to file
write.table(L67RP, file = paste0('./data/',date,'_L67RP_',auth,'.txt'),
          quote = FALSE, sep= "\t", row.names = FALSE)

# Extract SNPs from VCF format

# Load in Line 6 and Line 7 Data
######################
L6_df <- load_total_vcf('002683_Line-6')
L7_df <- load_total_vcf('002684_Line-7')

# Line 6
###########################
# Collect the intersection from variants of interest and line 6 in vcf format
L6_vcf <- semi_join(L6_df, L67RP, by =c('CHROM','POS','REF','ALT'))

# Adjust the column names
colnames(L6_vcf)[1] <- '#CHROM'

# Write the VCF file
write.table(L6_vcf, file = paste0('./data/',date,'_L6-ARK_',auth,'.vcf'),
            quote = FALSE, sep= "\t", row.names = FALSE)

# Line 7
############################

# Collect the intersection from variants of interest and line 6 in vcf format
L7_vcf <- semi_join(L7_df, L67RP, by =c('CHROM','POS','REF','ALT'))

# Adjust the column names
colnames(L7_vcf)[1] <- '#CHROM'

# Write the VCF file
write.table(L7_vcf, file = paste0('./data/',date,'_L7-ARK_',auth,'.vcf'),
            quote = FALSE, sep= "\t", row.names = FALSE)

#' ## Biology: What Genes are Mutated and What is Predicted Consequence?
#'
#+ VEP Annotation (Lines 6&7 and Arkansas intersections)

################################################################################
#####     VEP Annotation      ##################################################
################################################################################

# Load a docker image with the appropriate VEP release (92) and the galgal5 annotation information

# VEP Annotation performed in: "./scripts/20191206_vep-docker-example_steep.sh"


#' ## Load and Join VEP Annotations
#'
#+ Load and Join VEP Annotations

################################################################################
#####     Load and Join VEP Annotations      ###################################
################################################################################

# Line 6
#####################

# Load the VEP annotations
L6_file <- './data/20191208_L6-ARK-vep_steep.vcf'
L6 <- read.table(L6_file, sep = '\t', header = TRUE, skip = 2, comment.char = '', 
                 check.names = FALSE) %>% as_tibble()
# Rename columns
names(L6)[1] <- 'CHROM'
names(L6)[10] <- 'SAMPLE'
# Split the INFO column
L6 <- L6 %>% separate(INFO, into = c("INFO","CONSEQUENCE"), 
                      sep = "CSQ=", extra = "merge")
# Remove unneccassarry columns (learn to spell)
L6 <- L6 %>% dplyr::select(-QUAL, -FILTER, -FORMAT, -SAMPLE, -INFO)

# Splits the column by a delimiter and replicates rows (SPLIT_COL)
L6 <- L6 %>% mutate(CONSEQUENCE = strsplit(as.character(CONSEQUENCE), ",")) %>% 
        unnest(CONSEQUENCE)

# Split the CONSEQUENCE column
L6 <- L6 %>% separate(CONSEQUENCE, c("ALLELE", "CONSEQUENCE","IMPACT","SYMBOL","GENE_ID","FEATURE_TYPE","FEATURE","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_POS","CDS_POS","PRO_POS","AA","CODONS","VAR_EXIST","DIST","STRAND","FLAGS","SYMBOL_SOURCE","HGNC_ID","SIFT"), sep = "\\|")

# Select the columns of interest
L6 <- L6 %>% dplyr::select(CHROM,POS,REF,ALT,STRAND,CONSEQUENCE,IMPACT,SYMBOL,SIFT)

# Rename the all columns
#L6$LINE <- "LINE_6"

# Unique rows
L6 <- L6 %>% unique()


# Line 7
#####################

# Load the VEP annotations
L7_file <- './data/20191208_L7-ARK-vep_steep.vcf'
L7 <- read.table(L7_file, sep = '\t', header = TRUE, skip = 2, comment.char = '', 
                 check.names = FALSE) %>% as_tibble()
# Rename columns
names(L7)[1] <- 'CHROM'
names(L7)[10] <- 'SAMPLE'
# Split the INFO column
L7 <- L7 %>% separate(INFO, into = c("INFO","CONSEQUENCE"), 
                      sep = "CSQ=", extra = "merge")
# Remove unneccassarry columns (learn to spell)
L7 <- L7 %>% dplyr::select(-QUAL, -FILTER, -FORMAT, -SAMPLE, -INFO)

# Splits the column by a delimiter and replicates rows (SPLIT_COL)
L7 <- L7 %>% mutate(CONSEQUENCE = strsplit(as.character(CONSEQUENCE), ",")) %>% 
        unnest(CONSEQUENCE)

# Split the CONSEQUENCE column
L7 <- L7 %>% separate(CONSEQUENCE, c("ALLELE", "CONSEQUENCE","IMPACT","SYMBOL","GENE_ID","FEATURE_TYPE","FEATURE","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_POS","CDS_POS","PRO_POS","AA","CODONS","VAR_EXIST","DIST","STRAND","FLAGS","SYMBOL_SOURCE","HGNC_ID","SIFT"), sep = "\\|")

# Select the columns of interest
L7 <- L7 %>% dplyr::select(CHROM,POS,REF,ALT,STRAND,CONSEQUENCE,IMPACT,SYMBOL,SIFT)

# Rename the all columns
#L7$LINE <- "LINE_7"

# Unique rows
L7 <- L7 %>% unique()

# Perform a full-join
L6L7 <- full_join(L6, L7, copy = TRUE, suffix = c(".x", ".y"))

# Conver CHROM to factor
L6L7$CHROM <- as.factor(L6L7$CHROM)

# Incorporate sample to variant results with annotations
L67RP_anne <- left_join(L6L7, L67RP, by = c('CHROM','POS','REF','ALT'))

# Save the data to file and for subsequent analyses
write.table(L67RP_anne, file = paste0('./data/',date,'_L67RP-annotations-quantitative_',auth,'.txt'), quote = FALSE, sep= "\t", row.names = FALSE)













# Collect SNPs found only in L6 and R
L6R <- L67RP_anne %>% filter(
                (LINE_6 == 1 & REG == 1 & LINE_7 == 0 & PRO == 0)
        )
table(L6R$SYMBOL) %>% sort() %>% names()

L6R %>% filter(IMPACT %!in% c('MODIFIER','LOW'))


res <- L67RP_anne %>% filter( (SYMBOL == 'DNAH1') & 
                              ( IMPACT %!in% c('MODIFIER','LOW') ))
res %>% filter(!str_detect(SIFT,'tolerated'))
?str_detect
table(L67RP_anne$SYMBOL) %>% sort()



################################################################################
#####     Determine Genotype Similarity      ###################################
################################################################################

# Extract data frame of interest
#F1_df <- df %>% dplyr::select(PROBE_SET_ID,CHR,POS,AFFY_SNP_ID,F1_1, F1_2, F1_3, F1_4, F1_5, F1_6)

# Generate rows for number and percent similarity between columns
# Note: Row-oriented work is difficult in R--https://resources.rstudio.com/webinars/thinking-inside-the-box-you-can-do-that-inside-a-data-frame-april-jenny-bryan

F1_df <- df %>% dplyr::select(F1_1, F1_2, F1_3, F1_4, F1_5, F1_6)

# Convert rows into lists (20 seconds)
row_list <- F1_df %>% purrr::pmap(list)
# Count the genotypes (4 minutes)
F1_df$GENO_FREQ <- lapply(row_list,count_genos) %>% unlist()
# Calculate the percentage of similarity
F1_df %>% filter(GENO_FREQ < 1)



