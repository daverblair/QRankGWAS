
quantiles="0.25,0.5,0.75"
pheno_file=/path/to/phenotypes
phenotype_name='***'
subject_id_col='***'
bgen_file_path=/path/to/bgen
output_prefix=test
covariates='age,sex'
maf=0.001

####optionial####
#subset=path/to/subset file
#snps=path/to/snps file


python -m ../run-QRankGWAS.py \
    ${quantiles} \
    ${pheno_file} \
    ${phenotype_name} \
    ${subject_id_col} \
    ${bgen_file_path} \
    ${output_prefix} \
    --covariate_list ${covariates} \
    --maf ${maf} \
    --print_freq=200
    # --subject_subset ${subset} \
    # --variant_subset=${snps}
