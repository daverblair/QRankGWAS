# QRankGWAS
 Software implementing the QRank method described in Song el al. Bioninformatics 2017. It was adapted for use within the UK Biobank using Python. It was designed to be used as a command-line tool. Note, the R QRank version and the python implementation will not yield identical results. The python version uses Iterative Weighted Least Squares to fit the null regression models, while the R version uses the simplex method. Therefore, the two implementations can produce slightly different p-values, but they are highly consistent.

Installation:

pip install QRankGWAS

Following installation, specific details regarding the software can be found by running the following command:

python -m QRankGWAS -h

usage: __main__.py [-h] [--covariate_file_path COVARIATE_FILE_PATH]
                   [--sample_file_path SAMPLE_FILE_PATH]
                   [--covariate_list COVARIATE_LIST]
                   [--subject_subset SUBJECT_SUBSET]
                   [--variant_subset VARIANT_SUBSET] [--maf MAF]
                   [--print_freq PRINT_FREQ]
                   [--model_param_tol MODEL_PARAM_TOL] [--null_model_only]
                   [--randomize]
                   quantiles phenotype_file_path phenotype_name subject_id_col
                   bgen_file_path output_prefix

Performs GWAS for quantitative phenotype using QRank method from Song et al.
Bioinformatics 2017. Designed for use on UKBiobank.

positional arguments:
  quantiles             Comma-sep list of quantiles for analysis. Recommended
                        max: 3 quantiles.
  phenotype_file_path   Specifies path to phentoype file. Expects tab-
                        delimitted data WITH header. One column must contain
                        subject ids.
  phenotype_name        string value that provides column name of phenotype
  subject_id_col        string value that provides column name of subject ids
  bgen_file_path        path to bgen file containing genotypes. Expects .bgi
                        index file with same prefix as well.
  output_prefix         prefix (including path) of output file

optional arguments:
  -h, --help            show this help message and exit
  --covariate_file_path COVARIATE_FILE_PATH
                        Optional covariate file path. If not provided, then
                        covariates (if given) will be read from phenotype
                        file.
  --sample_file_path SAMPLE_FILE_PATH
                        Path to .sample file for bgen dataset. If not
                        provided, will search in path of .bgen file for
                        .sample file with same prefix.
  --covariate_list COVARIATE_LIST
                        List of covariates to include into the model. Provided
                        as comma-sep list (no spaces). DO NOT include
                        intercept; this automatically included.
  --subject_subset SUBJECT_SUBSET
                        Text file containing subject ids to include into the
                        analysis. Header with subject id must be present.
                        Single column expected, but can contain other columns
                        as well (tab-delimitted).
  --variant_subset VARIANT_SUBSET
                        Text file containing rsids (not snpids, 1 per line, no
                        header) for a set of variants to analyze. Note: this
                        is effective only when analyzing subset of variants
                        approximately 1/10th of total. Otherwise, likely
                        faster to cycle through entire file.
  --maf MAF             Minor allele frequency to filter variants. Default is
                        0.0001.
  --print_freq PRINT_FREQ
                        Progress printing frequency. Default: Every 1000
                        variants.
  --model_param_tol MODEL_PARAM_TOL
                        Tolerance for fitting regression models. Default: 1e-6
  --null_model_only     Flag that indicates the program to compute only the
                        null models, and these output results (plus residuals)
                        will be stored in the indicated directory. GWAS
                        p-values will not be computed.
  --randomize           Flag that indicates that GWAS should be conducted over
                        randomized rank scores. This is useful for calibrating
                        null statistics for randomization test. Note,
                        randomization occurs once and is NOT unique per
                        variant.
