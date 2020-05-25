import numpy as np
import os
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2
from bgen.reader import BgenFile
import io
import argparse
import time
import copy

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
                self.covariate_matrix['intercept']=np.ones(len(self.covariate_matrix))
        else:
            self.covariate_matrix=pd.DataFrame({'intercept':np.ones(len(self.phenotypes))},index=self.phenotypes.index)

        self.quantiles=np.array(quantiles)

        self._base_model=None
        self.null_model_results={}
        self.null_ranks={}

        self.covariate_matrix_numpy_view=self.covariate_matrix.to_numpy()
        self.VN=None


    def FitNullModels(self,tol=1e-8,maxiter=1000):
        self._base_model=sm.QuantReg(self.phenotypes,self.covariate_matrix,hasconst=True)
        for i,q in enumerate(self.quantiles):
            results=self._base_model.fit(q=q,p_tol=tol,max_iter=maxiter)

            self.null_ranks[q]=self._computeNullRanks(q,results.resid)
            self.null_model_results[q]=copy.deepcopy(results)


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
            self.included_subjects=included_subjects



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



    def BuildQRank(self,quantiles,output_file_prefix,param_tol=1e-8, max_fitting_iter=5000):
        self.qrank=QRank(self.Y,covariate_matrix=self.Z,quantiles=quantiles)
        self.qrank.FitNullModels(tol=param_tol,maxiter=max_fitting_iter)
        for q in quantiles:
            with open(output_file_prefix+'.NullModel.{0:g}.txt'.format(q),'w') as model_file:
                model_file.write(self.qrank.null_model_results[q].summary().as_text())

    def PerformGWASAdditive(self,output_file_prefix,maf_cutoff,print_freq=1000,variant_list=None):

        with open(output_file_prefix+'.QRankGWAS.txt','w',buffering=io.DEFAULT_BUFFER_SIZE*10) as output_file:
            output_file.write('snpid\trsid\tchrom\tpos\tmaj\tmin\tmaf\t')
            output_file.write('\t'.join(['p.{0:g}'.format(x) for x in self.qrank.quantiles])+'\tp.comp\n')

            variant_counter=0
            avg_elapsed_time=0.0
            block_counter=0
            start=time.time()

            if variant_list is None:
                total_num_variants=len(self.bgen_dataset)
                for variant in self.bgen_dataset:
                    if len(variant.alleles)==2:
                        dosage=variant.minor_allele_dosage[self.included_subjects_bgen_idx]
                        maf=dosage.sum()/(dosage.shape[0]*2.0)

                        if (maf>=maf_cutoff):
                            if variant.alleles.index(variant.minor_allele)!=1:
                                alleles=variant.alleles[::-1]
                            else:
                                alleles=variant.alleles


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

            else:
                total_num_variants=len(variant_list)
                failed_list=[]
                for rsid in variant_list:
                    try:
                        variant=self.bgen_dataset.with_rsid(rsid)
                        if len(variant.alleles)==2:
                            dosage=variant.minor_allele_dosage[self.included_subjects_bgen_idx]
                            maf=dosage.sum()/(dosage.shape[0]*2.0)

                            if (maf>=maf_cutoff):
                                if variant.alleles.index(variant.minor_allele)!=1:
                                    alleles=variant.alleles[::-1]
                                else:
                                    alleles=variant.alleles


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
                        else:
                            raise ValueError("Variant {0:s} is not biallelic".format(rsid))
                    except ValueError:
                        failed_list+=[rsid]

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

                print('Failed to process {0:d} variants due to missing ids/not biallelic'.format(len(failed_list)))
                print('Failed rsids written to {0:s}.FailedSNPs.txt\n'.format(output_file_prefix))
                with open(output_file_prefix+'.FailedSNPs.txt','w') as failed_file:
                    failed_file.write('\n'.join(failed_list)+'\n')

    def PerformGWASDominant(self,output_file_prefix,maf_cutoff,print_freq=1000,variant_list=None):

        with open(output_file_prefix+'.QRankGWAS.txt','w',buffering=io.DEFAULT_BUFFER_SIZE*10) as output_file:
            output_file.write('snpid\trsid\tchrom\tpos\tmaj\tmin\tmaf\t')
            output_file.write('\t'.join(['p.{0:g}'.format(x) for x in self.qrank.quantiles])+'\tp.comp\n')

            variant_counter=0
            avg_elapsed_time=0.0
            block_counter=0
            start=time.time()

            if variant_list is None:
                total_num_variants=len(self.bgen_dataset)
                for variant in self.bgen_dataset:
                    if len(variant.alleles)==2:
                        probs=variant.probabilities[self.included_subjects_bgen_idx]
                        dosage=np.sum(probs[:,1:],axis=1)
                        maf=variant.minor_allele_dosage[self.included_subjects_bgen_idx].sum()/(dosage.shape[0]*2)

                        if (maf>=maf_cutoff):
                            if variant.alleles.index(variant.minor_allele)!=1:
                                alleles=variant.alleles[::-1]
                            else:
                                alleles=variant.alleles


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

            else:
                total_num_variants=len(variant_list)
                failed_list=[]
                for rsid in variant_list:
                    try:
                        variant=self.bgen_dataset.with_rsid(rsid)
                        if len(variant.alleles)==2:
                            probs=variant.probabilities[self.included_subjects_bgen_idx]
                            dosage=np.sum(probs[:,1:],axis=1)
                            maf=variant.minor_allele_dosage[self.included_subjects_bgen_idx].sum()/(dosage.shape[0]*2)

                            if (maf>=maf_cutoff):
                                if variant.alleles.index(variant.minor_allele)!=1:
                                    alleles=variant.alleles[::-1]
                                else:
                                    alleles=variant.alleles


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
                        else:
                            raise ValueError("Variant {0:s} is not biallelic".format(rsid))
                    except ValueError:
                        failed_list+=[rsid]

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

                print('Failed to process {0:d} variants due to missing ids/not biallelic'.format(len(failed_list)))
                print('Failed rsids written to {0:s}.FailedSNPs.txt\n'.format(output_file_prefix))
                with open(output_file_prefix+'.FailedSNPs.txt','w') as failed_file:
                    failed_file.write('\n'.join(failed_list)+'\n')

    def PerformGWASRecssive(self,output_file_prefix,maf_cutoff,print_freq=1000,variant_list=None):

        with open(output_file_prefix+'.QRankGWAS.txt','w',buffering=io.DEFAULT_BUFFER_SIZE*10) as output_file:
            output_file.write('snpid\trsid\tchrom\tpos\tmaj\tmin\tmaf\t')
            output_file.write('\t'.join(['p.{0:g}'.format(x) for x in self.qrank.quantiles])+'\tp.comp\n')

            variant_counter=0
            avg_elapsed_time=0.0
            block_counter=0
            start=time.time()

            if variant_list is None:
                total_num_variants=len(self.bgen_dataset)
                for variant in self.bgen_dataset:
                    if len(variant.alleles)==2:
                        probs=variant.probabilities[self.included_subjects_bgen_idx]
                        dosage=probs[:,2]
                        maf=variant.minor_allele_dosage[self.included_subjects_bgen_idx].sum()/(dosage.shape[0]*2)

                        if (maf>=maf_cutoff):
                            if variant.alleles.index(variant.minor_allele)!=1:
                                alleles=variant.alleles[::-1]
                            else:
                                alleles=variant.alleles


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

            else:
                total_num_variants=len(variant_list)
                failed_list=[]
                for rsid in variant_list:
                    try:
                        variant=self.bgen_dataset.with_rsid(rsid)
                        if len(variant.alleles)==2:
                            probs=variant.probabilities[self.included_subjects_bgen_idx]
                            dosage=probs[:,2]
                            maf=variant.minor_allele_dosage[self.included_subjects_bgen_idx].sum()/(dosage.shape[0]*2)

                            if (maf>=maf_cutoff):
                                if variant.alleles.index(variant.minor_allele)!=1:
                                    alleles=variant.alleles[::-1]
                                else:
                                    alleles=variant.alleles


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
                        else:
                            raise ValueError("Variant {0:s} is not biallelic".format(rsid))
                    except ValueError:
                        failed_list+=[rsid]

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

                print('Failed to process {0:d} variants due to missing ids/not biallelic'.format(len(failed_list)))
                print('Failed rsids written to {0:s}.FailedSNPs.txt\n'.format(output_file_prefix))
                with open(output_file_prefix+'.FailedSNPs.txt','w') as failed_file:
                    failed_file.write('\n'.join(failed_list)+'\n')


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Performs GWAS for quantitative phenotype using QRank method from Song et al. Bioinformatic 2017. Designed for use on UKBiobank')
    parser.add_argument("quantiles",help="Comma-sep list of quantiles for analysis. Recommended max: 3 quantiles.",type=str )
    parser.add_argument("phenotype_file_path",help="Specifies path to phentoype file. Expects tab-delimitted data WITH header. One column must contain subject ids.",type=str )
    parser.add_argument("phenotype_name",help="string value that provides column name of phenotype",type=str)
    parser.add_argument("subject_id_col",help="string value that provides column name of subject ids",type=str)
    parser.add_argument("bgen_file_path",help="path to bgen file containing genotypes",type=str)
    parser.add_argument("output_prefix",help="prefix (including path) of output file",type=str)
    parser.add_argument("--covariate_file_path",help="Optional covariate file path. If not provided, then covariates (if given) will be read from phenotype file.",type=str)
    parser.add_argument("--sample_file_path",help="Path to .sample file for bgen dataset. If not provided, will search in path of .bgen file for .sample file with same prefix.",type=str)
    parser.add_argument("--covariate_list",help="List of covariates to include into the model. Provided as comma-sep list (no spaces)",type=str)
    parser.add_argument("--subject_subset",help="Text file containing subject ids to include into the analysis. Header with subject id must be present. Single column expected, but can contain other columns as well (tab-delimitted).",type=str)
    parser.add_argument("--variant_subset",help="Text file containing rsids (not snpids, 1 per line, no header) for a set of variants to analyze. Note: this is effective only when analyzing subset of variants approximately 1/10th of total. Otherwise, likely faster to cycle through entire file.",type=str)
    parser.add_argument("--maf",help="Minor allele frequency to filter variants. Default is 0.0001.",type=float)
    parser.add_argument("--print_freq",help="Progress printing frequency. Default: Every 1000 variants.",type=int)
    parser.add_argument("--null_model_tol",help="Tolerance for fitting null models. Default: 1e-6",type=float)
    parser.add_argument("--genetic_model",help="Genetic model for GWAS. Must be in ['Additive','Recessive','Dominant']. Default: 'Additive'",type=str)
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
        included_variants=pd.read_csv(args.variant_subset,sep='\t',header=None,names=['rsid'])
        print("Subset of variant ids read from {0:s}.\n".format(args.variant_subset))
    else:
        included_variants=None


    if args.print_freq is not None:
        print_freq=args.print_freq
    else:
        print_freq=1000

    if args.null_model_tol is not None:
        null_model_tol=args.print_freq
    else:
        null_model_tol=1e-6

    if args.genetic_model is not None:
        assert args.genetic_model in ['Additive','Dominant','Recessive'],"Genetic model must be in ['Additive','Recessive','Dominant']"
        genetic_model=args.genetic_model
    else:
        genetic_model='Additive'


    print("Step 1: Reading bgen, phenotype, and covariate files.\n",flush=True)
    gwas=QRankGWAS(bgen_file_path,phenotype_file_path,subject_id_col,covariate_file_path=covariate_file_path)
    #
    print("Step 2: Constructing phenotype and covariate data arrays.\n",flush=True)
    gwas.ConstructDataArrays(phenotype_name,covariate_cols=covariate_list,included_subjects=included_subjects)

    print("Step 3: Inferring Null Quantile Regression models.\n",flush=True)
    gwas.BuildQRank(quantiles,output_prefix,param_tol=null_model_tol, max_fitting_iter=5000)
    # #
    print("Step 4: Performing GWAS using {0:s} genetic model. Will print update every {1:d} variants\n".format(genetic_model,print_freq),flush=True)

    if genetic_model=='Dominant':
        if included_variants is not None:
            gwas.PerformGWASDominant(output_prefix,maf_cutoff=maf_cutoff,print_freq=print_freq,variant_list=included_variants['rsid'].values)
        else:
            gwas.PerformGWASDominant(output_prefix,maf_cutoff=maf_cutoff,print_freq=print_freq)
    elif genetic_model=='Recessive':
        if included_variants is not None:
            gwas.PerformGWASRecssive(output_prefix,maf_cutoff=maf_cutoff,print_freq=print_freq,variant_list=included_variants['rsid'].values)
        else:
            gwas.PerformGWASRecssive(output_prefix,maf_cutoff=maf_cutoff,print_freq=print_freq)
    else:
        if included_variants is not None:
            gwas.PerformGWASAdditive(output_prefix,maf_cutoff=maf_cutoff,print_freq=print_freq,variant_list=included_variants['rsid'].values)
        else:
            gwas.PerformGWASAdditive(output_prefix,maf_cutoff=maf_cutoff,print_freq=print_freq)
    print("Successfully completed GWAS.\n",flush=True)
