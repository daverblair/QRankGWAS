
quantiles="0.25,0.5,0.75"
pheno_file=~/Desktop/MendelianDiseaseProject/Analysis/Section_4/QRankGWAS/QRankDatasets/GBR/OMIM_ICD_86_QRankData.txt
phenotype_name='Bagged'
subject_id_col='IID'
bgen_file_path=../../Sandbox/ukbb/ukb_imp_chr14_v3.bgen
output_prefix=~/Desktop/MendelianDiseaseProject/Analysis/Section_4/QRankGWAS/Debug/tmp_A1AT_testsnp_central
covariates='sex,age_normalized,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10'
subset=~/Desktop/MendelianDiseaseProject/Data/UKBiobank/filtered_subject_ids_caucasian_rel_excluded.txt
maf=0.001


python QRankGWAS.py \
    ${quantiles} \
    ${pheno_file} \
    ${phenotype_name} \
    ${subject_id_col} \
    ${bgen_file_path} \
    ${output_prefix} \
    --covariate_list ${covariates} \
    --subject_subset ${subset} \
    --maf ${maf} \
    --variant_subset test_snps.txt \
    --genetic_model 'Recessive'
