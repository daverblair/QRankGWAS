# QRankGWAS
 Software implementing the QRank method described in Song el al. Bioninformatics 2017. It was adapted for use within the UK Biobank using Python. It was designed to be used as a command-line tool. Note, the R QRank version and the python implementation will not yield identical results. The python version uses Iterative Weighted Least Squares to fit the null regression models, while the R version uses the simplex method. Therefore, the two implementations can produce slightly different p-values, but they are generally highly correlated.

Installation:

pip install QRankGWAS

Following installation, specific details regarding the software can be found by running the following command:

python -m QRankGWAS -h
