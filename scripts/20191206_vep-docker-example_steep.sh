# Login to docker
docker login

WD=${PWD}

# Run the docker image as a container (make sure to be in the working directory -- not ./data/ but one up from that)
docker run -it \
-v ${PWD}:${PWD} \
--rm --name vep_r steepale/20191205_vep-galgal5-r_steep:release_92

# Make sure VEP is in PATH

# Change to working directory
WD='/Users/Alec/Documents/Bioinformatics/germplasm_collab/20191205_pilot-snp-compare_steep' # Chnage this to the one you were in own (AKA PWD from before)
cd ${WD}

# Note: This script does not work right now, it requires rerunning of prior scripts (cannot input vep annotated vcfs to these scripts)
# VEP annotation
for vcf in `ls -1 ./data/*_L*-ARK_steep.vcf`
do
# Adjust the output string
vcf_out=`echo ${vcf} | sed 's/-ARK/-ARK-vep/'`
#Perform VEP annotation, base annotation with SIFT added
vep \
-i ${vcf} \
-o ${vcf_out} \
--vcf \
--cache \
--offline \
--species gallus_gallus \
--force_overwrite \
--sift b
done

for vcf in `ls -1 ./data/*_L*-ARK_steep.vcf`
do
# Adjust the output string
vcf_out=`echo ${vcf} | sed 's/-ARK/-ARK-vep/'`
#Perform VEP annotation, base annotation with SIFT added
vep \
-i ${vcf} \
-o ${vcf_out} \
--vcf \
--cache \
--offline \
--species gallus_gallus \
--force_overwrite \
--sift b
done

# Exit the docker image: files will be saved on your computer
exit
