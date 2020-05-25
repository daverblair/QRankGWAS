import numpy as np
import statsmodels.api as sm
from scipy.stats import chi2
import pandas as pd



if __name__=='__main__':
    from QRankGWAS import QRank
    N=100000

    intercept=-0.1
    p=0.05
    q=1.0-p

    beta=0.05
    alpha=np.array([0.7,0.1])

    X = np.random.multinomial(1,pvals=[p*p,2.0*p*q,q*q],size=N)

    Z = np.hstack([np.random.binomial(1,0.5,size=(N,1)),np.random.normal(0.0,1.0,size=(N,1))])

    tmp=beta*np.sum(X*np.arange(3),axis=1,keepdims=True)+np.sum(alpha*Z,axis=1,keepdims=True)
    Y=intercept+tmp+(1.0+tmp)*np.random.normal(0.0,1.0,size=(N,1))

    pheno=pd.DataFrame({'phenotype':Y.ravel()},index=np.arange(Y.shape[0]))
    covariates=pd.DataFrame({'sex':Z[:,0],'age':Z[:,1]},index=np.arange(Y.shape[0]))

    quants=np.array([0.9,0.95,0.99])

    dosage=np.sum(X*np.arange(3),axis=1)

    qrank=QRank(pheno,covariate_matrix=covariates,quantiles=quants)
    qrank.FitNullModels()
    p=qrank.ComputePValues(dosage)
    # #subset=np.random.choice(np.arange(N),9500,replace=False)
    # #
    #
    from rpy2.robjects.packages import importr
    from rpy2.robjects import numpy2ri
    numpy2ri.activate()
    qrank_r=importr("QRank")
    r_output=qrank_r.QRank(Y,dosage,Z,quants)

    print("Python Rank P's: {0:g},{1:g},{2:g}".format(*p[0]))
    print("Python Composite P: {0:g}".format(p[1]))

    print("R Rank P's: {0:g},{1:g},{2:g}".format(*np.array(r_output[1])))
    print("R Composite P: {0:g}".format(np.array(r_output[0][0])))
