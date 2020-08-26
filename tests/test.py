"""
This code compare R and Python implementations directly. It was used for debugging.

"""



if __name__=='__main__':
    import sys
    from QRankGWAS import QRank
    import numpy as np
    import pandas as pd


    N=10000

    np.random.seed(123)
    intercept=-0.5
    p=0.05
    q=1.0-p

    beta=0.1
    alpha=np.array([0.7,0.1])

    X = np.random.multinomial(1,pvals=[p*p,2.0*p*q,q*q],size=N)

    Z = np.hstack([np.random.binomial(1,0.5,size=(N,1)),np.random.normal(0.0,1.0,size=(N,1))])

    effects=beta*np.sum(X*np.arange(3),axis=1,keepdims=True)+np.sum(alpha*Z,axis=1,keepdims=True)
    Y=intercept+effects+(1.0+effects)*np.random.normal(0.0,1.0,size=(N,1))

    pheno=pd.DataFrame({'phenotype':Y.ravel()},index=np.arange(Y.shape[0]))
    covariates=pd.DataFrame({'sex':Z[:,0],'age':Z[:,1]},index=np.arange(Y.shape[0]))

    quants=np.array([0.9,0.95,0.99])

    dosage=np.sum(X*np.arange(3),axis=1)

    dosage_df=pd.DataFrame(index=pheno.index)
    dosage_df['Minor Allele Dosage']=dosage

    qrank=QRank(pheno,covariate_matrix=covariates,quantiles=quants)
    qrank.FitNullModels(tol=1e-8)
    p=qrank.ComputePValues(dosage)
    betas,ci=qrank.FitAltModels(dosage_df)

    print("Rank P's: {0:g},{1:g},{2:g}".format(*p[0]))
    print("Composite P: {0:g}".format(p[1]))
