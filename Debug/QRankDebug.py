import numpy as np
import statsmodels.api as sm
from scipy.stats import chi2
import pandas as pd



if __name__=='__main__':
    import sys
    sys.path.append("..")
    from QRankGWAS import QRank

    import matplotlib.pyplot as plt
    from matplotlib import cm
    import seaborn as sns

    sns.set(context='talk',color_codes=True,style='ticks',font='Arial',font_scale=1,rc={'axes.linewidth':5,"font.weight":"bold",'axes.labelweight':"bold",'xtick.major.width':4,'xtick.minor.width': 2})
    cmap = cm.get_cmap('viridis', 12)
    color_list=[cmap(x) for x in [0.0,0.1,0.25,0.5,0.75,0.9,1.0]]
    grey_color=(0.5490196078431373, 0.5490196078431373, 0.5490196078431373)
    red_color = '#b42f2f'

    N=100000

    intercept=-0.1
    p=0.005
    q=1.0-p

    beta=-0.5
    alpha=np.array([0.7,0.1])

    X = np.random.multinomial(1,pvals=[p*p,2.0*p*q,q*q],size=N)

    Z = np.hstack([np.random.binomial(1,0.5,size=(N,1)),np.random.normal(0.0,1.0,size=(N,1))])

    tmp=beta*np.sum(X*np.arange(3),axis=1,keepdims=True)+np.sum(alpha*Z,axis=1,keepdims=True)
    Y=intercept+tmp+(1.0+tmp)*np.random.normal(0.0,1.0,size=(N,1))

    pheno=pd.DataFrame({'phenotype':Y.ravel()},index=np.arange(Y.shape[0]))
    covariates=pd.DataFrame({'sex':Z[:,0],'age':Z[:,1]},index=np.arange(Y.shape[0]))

    quants=np.array([0.9,0.95,0.99])

    dosage=np.sum(X*np.arange(3),axis=1)

    dosage_df=pd.DataFrame(index=pheno.index)
    dosage_df['Minor Allele Dosage']=dosage

    qrank=QRank(pheno,covariate_matrix=covariates,quantiles=quants)
    qrank.FitNullModels()
    p=qrank.ComputePValues(dosage)
    betas,ci=qrank.FitAltModels(dosage_df)
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

    print("Quantile Effect Params: {0:g},{1:g},{2:g}".format(*betas))
    print("Quantile Effects 95CI: ({0:s}),({1:s}),({2:s})".format(*[','.join(list(np.array(x,dtype=np.str))) for x in ci]))
    fig_data=pd.concat([pheno,dosage_df],axis=1)
    sns.boxenplot(x='Minor Allele Dosage',y='phenotype',data=fig_data,k_depth='proportion')
    plt.show()
