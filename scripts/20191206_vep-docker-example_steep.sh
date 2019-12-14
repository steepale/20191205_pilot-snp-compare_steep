# Login to docker
docker login

# Run the docker image as a container
docker run -it \
-v /Users/Alec/Documents/Bioinformatics/germplasm_collab/20191205_pilot-snp-compare_steep:/Users/Alec/Documents/Bioinformatics/germplasm_collab/20191205_pilot-snp-compare_steep \
--rm --name vep_r steepale/20191205_vep-galgal5-r_steep:release_92

# Make sure VEP is in PATH

# Change to working directory
cd /Users/Alec/Documents/Bioinformatics/germplasm_collab/20191205_pilot-snp-compare_steep

# Run VEP
vep \
-i ./data/20191208_L6-ARK_steep.vcf \
-o ./data/20191208_L6-ARK-vep_steep.vcf \
--vcf \
--cache \
--offline \
--species gallus_gallus \
--force_overwrite \
--sift b

# Run VEP
vep \
-i ./data/20191208_L7-ARK_steep.vcf \
-o ./data/20191208_L7-ARK-vep_steep.vcf \
--vcf \
--cache \
--offline \
--species gallus_gallus \
--force_overwrite \
--sift b

# # Run VEP
vep \
-i ./data/002684_Line-7_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf \
-o ./data/002684_Line-7_germline_99.9_vep.g.vcf \
--vcf \
--cache \
--offline \
--species gallus_gallus \
--force_overwrite \
--sift b

vep \
-i ./data/002683_Line-6_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf.gz \
-o ./data/002683_Line-6_germline_99.9_vep.g.vcf \
--vcf \
--cache \
--offline \
--species gallus_gallus \
--force_overwrite \
--sift b

# Exit the docker image: files will be saved on your computer
exit
