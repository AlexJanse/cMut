# cMut
cMut (clustered mutations) is an R package that allows you to find clustered mutations, find mutational  
patterns you're looking for and validate the results.

## Purpose of cMut
# What can the cMut package do? 
A subset of the mutations that happen in the human genome occur very close to each other. 
Such local enrichment of mutations are often called mutation clusters. It has been hypothesized that 
clustered mutations occur by mechanisms that are different from the mechanisms of non-clustered mutations, 
making this group particularly interesting to study.
There are mutational influences that tend to cause mutations of very specific types. One example is the 
human protein APOBEC3A, which causes C > T mutations at nucleotide positions that are preceded by CT or TT 
nucleotides and followed by an A nucleotide at high specificity. APOBEC3A is suspected to cause a subset 
of clustered mutations.
The cMut package analyses clustered SNV mutations for the presence of custom mutation patterns and can 
calculate whether these patterns are enriched. 

## Get Started
The following instruction will tell you what you need to use cMut and how to install it.

### Prerequisites
The software that you need:

[R](https://www.r-project.org/) (at least version 3.5)

Highly recommended extra software:

[RStudio](https://www.rstudio.com/)

### Install
Start R (or RStudio) and use the following command:
```
install.packages("cMut")
```
Please also checkout the Vignette for intructions how to use this package:
```
vignette("analysis_of_clusterpattterns",package = "cMut")
```