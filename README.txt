# cMut
cMut (clustered mutations) is an R package that allows you to find clustered mutations, find mutational  
patterns you're looking for and validate the results.

## Purpose of cMut
When learning about DNA, we've learned that all the cells in our body that contain DNA have the same 
DNA throughout the body. Although that is mostly true, when we get older the more variations arise in 
different cells due exposure to mutators and different speed of decay of the telomeres. 
The variations are the best displayed when looking at the mutations in DNA sequence (sequence of 
nucleotides A, C, G and T in different orders and frequencies). These mutations can exist in many 
shapes and sizes; from single nucleotide variant (SNV) till large copy number variants (CNV).  
In this package we're interested in finding and annotating clustered SNV (cSNV a.k.a. multi-nucleotide 
mutation; MNM). We also try to link them with known mutator patterns and check if the enrichment is 
valid.  

Therefore the following questions arise:
* **Bioinformatics:**
  Is there an enrichment of clustered mutations that match with certain mutation patterns?
* **Biology:**
  Which mutators are responsible for the clustered mutations?

With this package we hope to answer the bioinformatics question and therefore give a clue about the 
answer for the biology question.

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
