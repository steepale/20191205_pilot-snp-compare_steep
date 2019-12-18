# An Exploration of Germline SNPs Across Chicken Lines

High confidence germline variants were compared between 4 chicken lines (2 from Michigan and 2 from Arkansas)--Line 6, Line 7, Progressor, and Regressor--in an attempt at a 'birds-eye view' of genotype distributions. The biological hypotheses guiding this pilot analysis were created by Beth Krehbiel and Harvey Blackburn, and the analysis was performed with the help of Alec Steep.  

## Information about Line 6 and Line 7
Line 6 and Line 7 were maintained at the Avian Disease and Oncology Lab in East Lansing, Michigan. They represent highly inbred lines that are resistant (Line 6) and susceptible (Line 7) to Marek's disease pathogensis, which is a lymphoproliferative disease driven by a herpesvirus. One might infer that these lines are resistant and susceptible to cancer in general. 

## Information about the Arkansas Progressor and Regressor Lines
The Progressor line demonstrates infertility while the regressor line demonstrates fertility.

# Order of analysis and explanation of scripts
* 20191205_pilot-snp-compare_steep.R  
    The initial data wrangling script. This script feeds in large portions of data from different sources and formats the data for proper sorting. As this script lays the foundation for all analyses, this script should be held to the highest scrutiny and errors should be searched for; that being said, ample attention was paid to the construction of this script.
    * 20191205_snp-ref_steep.R  
        Checks the reference genome builds associated with each dataset. This analysis concluded that the Arkansas lines were in reference to galgal4 and the Michigan lines were in reference to galgal5.
    * 20191206_vep-docker-example_steep.sh (Requires docker or singularity)  
        Allows for VEP annotation for variants of interest. VEP is a complex software and it notorious for its difficult installation. These scripts allow for use of vep in nearly any computing environment with a docker image.
* 20191209_pilot-snp-eda_steep.R  & 20191209_pilot-snp-eda_steep.html  
    The start of the global analysis of genotypes across samples. Identifies the intersection of genotypes and illustrates PCA analysis. The html is a friendly way of visualizing thess data.

# VEP Annotation Summaries
VEP generates global stats associated with annotated variants and puts them in a html file. Those files are:
* `./data/20191214_L6-ARK-vep_steep.vcf_summary.html`
* `./data/20191214_L7-ARK-vep_steep.vcf_summary.html`

# TODOs:
* Clean up mislabeled non-ref genos (1/2)
* Make sure that genotype assignments and quantitative representations of those genotypes are equivalent across all rows