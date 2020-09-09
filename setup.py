#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 09:11:51 2019

@author: davidblair
"""

###

import setuptools
import re

version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('QRankGWAS/QRankGWAS.py').read(),
    re.M).group(1)


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="QRankGWAS",
    version=version,
    author="David Blair",
    author_email="david.blair@ucsf.edu",
    description="Python implementation of the QRank method described in Song et al Bioninformatics 2017.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/daverblair/QRankGWAS",
    packages=["QRankGWAS"],
    entry_points = {
        "console_scripts": ['QRankGWAS = QRankGWAS.QRankGWAS:main']
        },
    install_requires=[
        'argparse>=1.1',
        'numpy>=1.19.0',
        'pandas>=1.0.5',
        'statsmodels>=0.11.1',
        'scipy>=1.5.2',
        'scikit-learn>=0.22.1',
        'bgen==1.2.7'
        ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
