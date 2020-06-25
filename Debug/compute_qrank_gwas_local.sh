
quantiles="0.995,0.999"
pheno_file=~/Desktop/MendelianDiseaseProject/Analysis/Section_4/QRankGWAS/QRankDatasets/GBR/OMIM_ICD_125_QRankData.txt
phenotype_name='Bagged'
subject_id_col='IID'
bgen_file_path=/Volumes/WD_Backup/UKBiobank/genetic_data/imputed_genotypes/CHROM_3/ukb_imp_chr3_v3.bgen
#bgen_file_path=//ukbb/ukb_imp_chr14_v3.bgen
output_prefix=~/Desktop/MendelianDiseaseProject/QRankGWAS/Debug/tmp_Neurofibromatosis_testsnps
covariates='sex,age_normalized,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10'
subset=~/Desktop/MendelianDiseaseProject/Data/UKBiobank/filtered_subject_ids_caucasian_rel_excluded.txt
snps=test_snp_nf1.txt
maf=0.001
# variant_subset_file=~/Desktop/MendelianDiseaseProject/Analysis/Section_4/BuildVariantTables-4/RandomSubsetSNPs/CHROM_14/RandomSubset.txt
# variant_subset_file=test_snps.txt


python ../QRankGWAS.py \
    ${quantiles} \
    ${pheno_file} \
    ${phenotype_name} \
    ${subject_id_col} \
    ${bgen_file_path} \
    ${output_prefix} \
    --covariate_list ${covariates} \
    --subject_subset ${subset} \
    --maf ${maf} \
    --genetic_model 'Additive' \
    --print_freq=200 \
    --variant_subset=${snps}
    # --randomize
