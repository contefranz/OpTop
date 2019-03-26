MacOS/Linux: [![Build Status](https://travis-ci.org/contefranz/OpTop.svg?branch=master)](https://travis-ci.org/contefranz/OpTop)
Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/contefranz/OpTop/branch/master?svg=true)](https://ci.appveyor.com/api/projects/status/contefranz/OpTop/branch/master)

# OpTop: detect the optimal number of topics from a pool of LDA models

## Overview

**OpTop** is the direct consequence of the paper by Lewis and Grossetti (2019). 
The paper is in an advanced draft and is about to be submitted to a journal.

Latent Dirichelet Allocation (LDA) was developed by Blei, Ng, and Jordan in 2003
[Blei et al., (2003)] and is based on 
the idea that a corpus can be represented by a set of topics. LDA has been used 
extensively in computational linguistics, is replicable, and is automated so it 
cannot be influenced by researcher prejudice. LDA uses a likelihood approach to 
discover clusters of text, namely topics that frequently appear in a corpus.

One of the open challenges in topic modeling is to rigorously determine the 
optimal number of topics for a corpus. Since there are no well-defined heuristic 
approaches, researchers rely on iterative trial-and-error approaches.
A standard approach is to determine which specification is the least perplexed 
by the test sets. Perplexity is based on the intuition that a high degree of 
similarity, identified as a low level of perplexity, can be used to determine 
the appropriate number of topics [Blei et al., (2003); Hornik and Gr&uuml;n, (2011)].

**OpTop** introduces a set of parametric tests to find the optimal number of topics 
in a collection of LDA models. **OpTop** also includes 
several tests to explore topic stability and redundancy. 

## Installation

The package is not on CRAN yet. You can install the development version as follows:
``` r
# Install the development version from Github:
devtools::install_github( "contefranz/OpTop" )
```

## Functions

All the procedures described in the paper will be implemented in this package.
The package is in early alpha stage and contains two functions:

* `get_topic_models()`: handy function to immediately get the list of topic models
the user wants to process from a specified environment;

* `word_proportions()`: computes word proportions from a `corpus` object created 
by __quanteda__ [Benoit et al. (2018)];

* `optimal_topic()`: implements _Test 1_ from the methodological paper 
[Lewis and Grossetti (2019)].

* `topic_stability()`: implements _Test 2_ from the methodological paper 
[Lewis and Grossetti (2019)].

More functions which implement the other tests are to come in future releases.

## Bug Reporting

Bugs and issues can be reported at
[https://github.com/contefranz/OpTop/issues](https://github.com/contefranz/OpTop/issues).

## Bibliography

1. Lewis, C. and Grossetti, F. (2019 - forthcoming): _A Statistical Approach
for Optimal Topic Model Identification_.
2. Blei, D. M., Ng, A. Y., and Jordan, M. I. (2003). _Latent Dirichlet Allocation_.
Journal of Machine Learning Research, 3(Jan):993–1022.
3. Benoit K., Watanabe K., Wang H., Nulty P., Obeng A., M&uuml;ller S., Matsuo A.
(2018): _`quanteda`: An R package for the
quantitative analysis of textual data_. Journal of Open Source Software, 3(30), 774. doi: 10.21105/joss.00774
(URL: http://doi.org/10.21105/joss.00774), URL: https://quanteda.io)

***
  
