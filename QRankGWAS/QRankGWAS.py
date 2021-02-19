import os
import io
import argparse
import time
import copy
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2
from sklearn.utils import shuffle
from bgen.reader import BgenFile


__version__ = "0.0.9"


class QRank:

    def _computeNullRanks(self,tau, residuals):
        new_resid=np.clip(residuals.to_numpy(),a_min=0.0,a_max=None)
        new_resid[new_resid>0.0]=1.0
        return (tau-(1-new_resid)).reshape(-1,1)

    def __init__(self,phenotypes,covariate_matrix=None,quantiles=[0.25,0.5,0.75],intercept_included=False):

        assert isinstance(phenotypes, pd.DataFrame), "Expects pandas DataFrame object for phenotypes"
        self.phenotypes=phenotypes

        if covariate_matrix is not None:
            assert isinstance(covariate_matrix, pd.DataFrame), "Expects pandas DataFrame object for covariate matrix"
            self.covariate_matrix=covariate_matrix
            if intercept_included==False:
                self.covariate_matrix=self.covariate_matrix.assign(intercept=pd.Series([1]*covariate_matrix.shape[0]).values)

        else:
            self.covariate_matrix=pd.DataFrame({'intercept':np.ones(len(self.phenotypes))},index=self.phenotypes.index)

        self.quantiles=np.array(quantiles)

        self._base_model=None
        self.null_model_results={}
        self.null_ranks={}

        self.covariate_matrix_numpy_view=self.covariate_matrix.to_numpy()
        self.VN=None


    def FitAltModels(self,dosage_df,tol=1e-8,maxiter=1000,ci_level=0.05):
        if isinstance(dosage_df, pd.DataFrame) is False:
            new_df=pd.DataFrame(index=self.phenotypes.index)
            new_df['rsid_0']=dosage_df
            dosage_df=new_df

        alt_model=sm.QuantReg(self.phenotypes,pd.concat([self.covariate_matrix,dosage_df],axis=1),hasconst=True)
        betas=np.zeros(self.quantiles.shape[0])
        cis=np.zeros((self.quantiles.shape[0],2))
        for i,q in enumerate(self.quantiles):
            results=alt_model.fit(q=q,p_tol=tol,max_iter=maxiter)
            betas[i]=results.params[-1]
            cis[i]=results.conf_int(alpha=ci_level).loc[dosage_df.columns[0]]
        return betas,cis

    def FitNullModels(self,tol=1e-8,maxiter=1000,randomize=False):
        self._base_model=sm.QuantReg(self.phenotypes,self.covariate_matrix,hasconst=True)
        for i,q in enumerate(self.quantiles):
            results=self._base_model.fit(q=q,p_tol=tol,max_iter=maxiter)
            self.null_model_results[q]=copy.deepcopy(results)
            self.null_ranks[q]=self._computeNullRanks(q,results.resid)
        if randomize:
            rand_index=np.random.permutation(np.arange(self.phenotypes.shape[0]))
            for i,q in enumerate(self.quantiles):
                self.null_ranks[q]=self.null_ranks[q][rand_index]

        self.VN=np.zeros((self.quantiles.shape[0],self.quantiles.shape[0]))
        for i in range(self.quantiles.shape[0]):
            for j in range(self.quantiles.shape[0]):
                self.VN[i,j]=min(self.quantiles[i],self.quantiles[j])-self.quantiles[i]*self.quantiles[j]

    def ComputePValues(self,dosage):
        if len(dosage.shape)!=2:
            dosage=dosage.reshape(-1,1)

        lin_mod_x=np.linalg.lstsq(self.covariate_matrix_numpy_view,dosage,rcond=None)
        xstar=dosage-np.dot(self.covariate_matrix_numpy_view,lin_mod_x[0])
        SN=np.zeros(self.quantiles.shape)

        for i,q in enumerate(self.quantiles):
            SN[i]=np.sum(self.null_ranks[q]*xstar)

        VN2=self.VN*np.sum(xstar*xstar)
        pvals_each=chi2(1).sf((SN*SN)/np.diag(VN2))

        e=np.linalg.solve(np.linalg.cholesky(VN2).T,np.identity(VN2.shape[0]))

        SN2=np.dot(e.T,SN)
        pval_composite=chi2(self.quantiles.shape[0]).sf(np.sum(SN2*SN2))


        return pvals_each,pval_composite


class QRankGWAS:

    def __init__(self,bgen_file_path,phenotype_file_path,index_column_name,covariate_file_path=None,sample_file_path=None):

        """
        This software is meant to be called from the command line, so no documentation is included here. Note, the code here is a bit verbose, which was done in an attempt to minimize the number of function calls given the need to perform millions of calls. This could likely be optimized in a better way.


        """
        self.index_column_name=index_column_name

        assert os.path.isfile(bgen_file_path),"bgen file does not exist"

        if os.path.isfile(bgen_file_path+'.bgi') is False:
            print("Warning: No bgen index (.bgi) file provided in same directory as bgen file. Initial reading of the bgen is MUCH faster with index file. ")

        if sample_file_path is not None:
            assert os.path.isfile(sample_file_path),"sample file does not exist at provided location"
        else:
            sample_file_path=bgen_file_path.strip('bgen')+'sample'
            if os.path.isfile(sample_file_path) is False:
                raise FileNotFoundError("No sample file at {0:s}. A sample file must be provided.".format(sample_file_path))


        print('Reading bgen file from {0:s} using sample file {1:s}. If these seem like an error, kill program.'.format(bgen_file_path,sample_file_path))

        self.bgen_dataset=BgenFile(bgen_file_path,sample_path=sample_file_path)

        if os.path.isfile(phenotype_file_path):
            self.phenotype_dataset = pd.read_csv(phenotype_file_path,sep='\t',index_col=index_column_name)
        else:
            raise FileNotFoundError("No phenotype file at provided location")

        if covariate_file_path is not None:
            if os.path.isfile(covariate_file_path):
                self.covariate_dataset = pd.read_csv(covariate_file_path,sep='\t',index_col=index_column_name)
            else:
                raise FileNotFoundError("No covariate file at provided location")
        else:
            print("No covariate file provided. Will use phenotype file for covariates.\n",flush=True)
            self.covariate_dataset=self.phenotype_dataset


    def ConstructDataArrays(self,phenotype_name,covariate_cols=None,included_subjects=None):
        if included_subjects is None:
            self.included_subjects=self.phenotype_dataset.index.to_numpy()
        else:
            self.included_subjects=np.intersect1d(included_subjects,self.phenotype_dataset.index.to_numpy())

        self.Y=self.phenotype_dataset.loc[self.included_subjects][[phenotype_name]]
        if covariate_cols is not None:
            self.Z=self.covariate_dataset.loc[self.included_subjects][covariate_cols]
        else:
            self.Z=None

        sample_vals_np = np.array(self.bgen_dataset.samples,dtype=self.included_subjects.dtype)
        sample_vals_np_sorted=np.sort(sample_vals_np)
        sample_vals_np_idx_sorted=np.argsort(sample_vals_np)
        conv_dict=dict(zip(sample_vals_np_sorted,sample_vals_np_idx_sorted))
        self.included_subjects_bgen_idx=np.array([conv_dict[x] for x in self.included_subjects])



    def BuildQRank(self,quantiles,param_tol=1e-8, max_fitting_iter=5000,output_file_prefix=None,randomize=False):
        self.qrank=QRank(self.Y,covariate_matrix=self.Z,quantiles=quantiles)
        self.qrank.FitNullModels(tol=param_tol,maxiter=max_fitting_iter,randomize=randomize)
        if output_file_prefix is not None:
            residual_table=pd.DataFrame(index=self.included_subjects)

            for q in quantiles:
                residual_table['q.{0:g}.residuals'.format(q)]=self.qrank.null_model_results[q].resid
                with open(output_file_prefix+'.NullModelResults.{0:g}.txt'.format(q),'w') as model_file:
                    model_file.write(self.qrank.null_model_results[q].summary().as_text())
                    self.qrank.null_model_results[q].save(output_file_prefix+'.NullModel.{0:g}.pth'.format(q))
            residual_table.to_csv(output_file_prefix+'.NullModelResiduals.txt',sep='\t')


    def PerformGWASAdditive(self,output_file_prefix,maf_cutoff,print_freq=1000,variant_list=None):

        if variant_list is None:
            total_num_variants=len(self.bgen_dataset)
            variant_iterator=self.bgen_dataset
        elif len(variant_list) > 1000:
            print("Adjusting bgen index to drop excluded variants from the analysis. This may take several minutes up front.")
            all_rsids=self.bgen_dataset.rsids()
            rsid_table=pd.DataFrame({'rsid':all_rsids,'bgen_index':np.arange(len(all_rsids))})
            rsid_table.set_index('rsid',inplace=True,drop=False)
            rsid_table=rsid_table.drop(np.intersect1d(variant_list,rsid_table.index.to_numpy()))
            self.bgen_dataset.drop_variants(rsid_table['bgen_index'].to_list())

            total_num_variants=len(self.bgen_dataset)
            def variant_iterator_func(num_var):
                for x in range(num_var):
                    yield self.bgen_dataset[x]

            variant_iterator=variant_iterator_func(total_num_variants)
        else:
            # use a custom generator, load in real time
            #
            def variant_iterator_func(v_list):
                for x in v_list:
                    yield self.bgen_dataset.with_rsid(x)

            variant_iterator=variant_iterator_func(variant_list)

        with open(output_file_prefix+'.Additive.QRankGWAS.txt','w',buffering=io.DEFAULT_BUFFER_SIZE*10) as output_file:
            output_file.write('snpid\trsid\tchrom\tpos\tmaj\tmin\tmaf\t')
            output_file.write('\t'.join(['p.{0:g}'.format(x) for x in self.qrank.quantiles])+'\tp.comp\n')


            variant_counter=0
            avg_elapsed_time=0.0
            block_counter=0
            start=time.time()

            for variant in variant_iterator:
                if len(variant.alleles)==2:
                    dosage=variant.minor_allele_dosage[self.included_subjects_bgen_idx]
                    maf=dosage.sum()/(dosage.shape[0]*2.0)

                    if (maf>=maf_cutoff):
                        if (variant.alleles.index(variant.minor_allele)==1) and (maf<=0.5):
                            alleles=variant.alleles
                        else:
                            alleles=variant.alleles[::-1]



                        output_file.write('{0:s}'.format(variant.varid))
                        output_file.write('\t{0:s}'.format(variant.rsid))
                        output_file.write('\t{0:s}'.format(variant.chrom))
                        output_file.write('\t{0:d}'.format(variant.pos))
                        output_file.write('\t{0:s}'.format(alleles[0]))
                        output_file.write('\t{0:s}'.format(alleles[1]))
                        output_file.write('\t{0:.8g}'.format(maf))
                        pvals=self.qrank.ComputePValues(dosage)
                        for p in pvals[0]:
                            output_file.write('\t{0:.8g}'.format(p))
                        output_file.write('\t{0:.8g}'.format(pvals[1]))
                        output_file.write('\n')
                variant_counter+=1
                if (variant_counter) % print_freq==0:
                    end=time.time()
                    block_counter+=1
                    elapsed=end-start
                    print('Processed {0:d} of {1:d} variants ({2:.1f}% of total)'.format(variant_counter,total_num_variants,round((variant_counter/total_num_variants)*1000.0)/10.0),flush=True)
                    print('Elapsed time {0:.2f} sec'.format(elapsed))
                    avg_elapsed_time = ((avg_elapsed_time*(block_counter-1)+elapsed)/block_counter)
                    print('Estimated Total Time Required: {0:.2f} hours\n'.format(((total_num_variants/print_freq)*avg_elapsed_time)/3600))
                    start=time.time()



def main():
    parser = argparse.ArgumentParser(description='Performs GWAS for quantitative phenotype using QRank method from Song et al. Bioinformatics 2017. Designed for use on UKBiobank.')
    parser.add_argument("quantiles",help="Comma-sep list of quantiles for analysis. Recommended max: 3 quantiles.",type=str )
    parser.add_argument("phenotype_file_path",help="Specifies path to phentoype file. Expects tab-delimitted data WITH header. One column must contain subject ids.",type=str )
    parser.add_argument("phenotype_name",help="string value that provides column name of phenotype",type=str)
    parser.add_argument("subject_id_col",help="string value that provides column name of subject ids",type=str)
    parser.add_argument("bgen_file_path",help="path to bgen file containing genotypes. Expects .bgi index file with same prefix as well.",type=str)
    parser.add_argument("output_prefix",help="prefix (including path) of output file",type=str)
    parser.add_argument("--covariate_file_path",help="Optional covariate file path. If not provided, then covariates (if given) will be read from phenotype file.",type=str)
    parser.add_argument("--sample_file_path",help="Path to .sample file for bgen dataset. If not provided, will search in path of .bgen file for .sample file with same prefix.",type=str)
    parser.add_argument("--covariate_list",help="List of covariates to include into the model. Provided as comma-sep list (no spaces). DO NOT include intercept; this automatically included.",type=str)
    parser.add_argument("--subject_subset",help="Text file containing subject ids to include into the analysis. Header with subject id must be present. Single column expected, but can contain other columns as well (tab-delimitted).",type=str)
    parser.add_argument("--variant_subset",help="Text file containing rsids (not snpids, 1 per line, no header) for a set of variants to analyze. Note: this is effective only when analyzing subset of variants approximately 1/10th of total. Otherwise, likely faster to cycle through entire file.",type=str)
    parser.add_argument("--maf",help="Minor allele frequency to filter variants. Default is 0.0001.",type=float)
    parser.add_argument("--print_freq",help="Progress printing frequency. Default: Every 1000 variants.",type=int)
    parser.add_argument("--model_param_tol",help="Tolerance for fitting regression models. Default: 1e-6",type=float)
    parser.add_argument("--null_model_only",help="Flag that indicates the program  to compute only the null models, and these output results (plus residuals) will be stored in the indicated directory. GWAS p-values will not be computed.",action="store_true")
    parser.add_argument("--randomize",help="Flag that indicates that GWAS should be conducted over randomized rank scores. This is useful for calibrating null statistics for randomization test. Note, randomization occurs once and is NOT unique per variant.",action="store_true")
    args = parser.parse_args()



    quantiles=np.array(args.quantiles.split(','),dtype=np.float32)
    phenotype_file_path=args.phenotype_file_path
    phenotype_name=args.phenotype_name
    subject_id_col=args.subject_id_col
    bgen_file_path=args.bgen_file_path
    output_prefix=args.output_prefix

    #default to None, so can subsume value even if none
    covariate_file_path=args.covariate_file_path




    print('#'*20)
    print('Initiating QRankGWAS Analysis')
    print('Phenotype File: {0:s}'.format(phenotype_file_path))
    print('Quantiles: {0:s}'.format(args.quantiles))
    print('Phenotype Name: {0:s}'.format(phenotype_name))
    print('Subject ID Column: {0:s}'.format(subject_id_col))
    print('bgen dataset: {0:s}'.format(bgen_file_path))
    print('Output File: {0:s}.QRankGWAS.txt'.format(output_prefix))
    print('#'*20+'\n')

    sample_file_path=args.sample_file_path
    if sample_file_path is None:
        print('Sample file not provided. Will attempt to read one from same directory as bgen file.\n')

    if args.covariate_list is not None:
        covariate_list=args.covariate_list.split(',')
        print("Covariates: "+ args.covariate_list+'\n')
    else:
        print("No Covariates included into the model. Quantile regression will be performed with intercept only.\n")
        covariate_list=None

    if args.maf is not None:
        print('MAF Filter: {0:g}.\n'.format(args.maf))
        maf_cutoff=args.maf
    else:
        print('MAF Filter: Using default of 0.0001.\n')
        maf_cutoff=0.0001

    if args.subject_subset is not None:
        included_subjects=pd.read_csv(args.subject_subset,sep='\t',index_col=subject_id_col).index.to_numpy()
        print("Subset of subject ids read from {0:s}.\n".format(args.subject_subset))
    else:
        included_subjects=None

    if args.variant_subset is not None:
        included_variants=pd.read_csv(args.variant_subset,sep='\t',header=None,names=['rsid'])['rsid'].values
        print("Subset of variant ids read from {0:s}.\n".format(args.variant_subset))
    else:
        included_variants=None


    if args.print_freq is not None:
        print_freq=args.print_freq
    else:
        print_freq=1000

    if args.model_param_tol is not None:
        model_param_tol=args.model_param_tol
    else:
        model_param_tol=1e-6


    if args.randomize is True:
        print('Randomization invoked. This will generate a null distribution of p-values.\n')

    print("Step 1: Reading bgen, phenotype, and covariate files.\n",flush=True)
    gwas=QRankGWAS(bgen_file_path,phenotype_file_path,subject_id_col,covariate_file_path=covariate_file_path)
    #
    print("Step 2: Constructing phenotype and covariate data arrays.\n",flush=True)
    gwas.ConstructDataArrays(phenotype_name,covariate_cols=covariate_list,included_subjects=included_subjects)


    if args.null_model_only is False:
        print("Step 3: Inferring Null Quantile Regression models.\n",flush=True)

        gwas.BuildQRank(quantiles,param_tol=model_param_tol, max_fitting_iter=5000,randomize=args.randomize)

        print("Step 4: Performing GWAS using additive genetic model. Will print update every {0:d} variants\n".format(print_freq),flush=True)



        gwas.PerformGWASAdditive(output_prefix,maf_cutoff,print_freq=print_freq,variant_list=included_variants)
    else:
        assert args.randomize is False,"Randomization has no effect on null model inference. Please remove this option if fitting null model only. "
        print("Step 3: Inferring Null Quantile Regression models only. No GWAS will be performed.\n",flush=True)
        gwas.BuildQRank(quantiles,param_tol=model_param_tol, max_fitting_iter=5000,output_file_prefix=output_prefix)

    print("Analysis successfully completed!")
