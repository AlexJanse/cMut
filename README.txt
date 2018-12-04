# cMUT
cMut (clustered mutations) is an R package that allows you to find clustered mutations, find mutational 
patterns you're looking for and validate the results.

## History
_De Novo_ mutations (DNM) are DNA mutations that can only be found in the offspring and not in the parents.

Previous research has found that 2 á 3 percent of these DNMs are located within clusters (2 or more 
mutations within 20 Kb range) ([Goldmann et al. 2016](https://www.nature.com/articles/ng.3597); [GoldMann et al 2018](https://www.nature.com/articles/s41588-018-0071-6)).
This is much more than the 0.11% you would expected if these clusters arrive by chance alone ([Besenbacher et al. 2016](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006315)).
Also the increased number of transitions ( G <> C, A <> T ) in clustered DNMs (cDNMs) compared with 
non-clustered DNMs (ncDNMs) is another evidence that cDNMs most likely arise from different processes than ncDNMs 
([Besenbacher et al. 2016](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006315)).

Therefore the following questions arise:
* **Bioinformatics:**
  Is there an enrichment of clustered de novo mutations that match with certain mutation patterns?
* **Biology:**
  Which mutators are responsible for the clustered de novo mutations?

With this package we hope to answer the bioinformatics question and therefore give a clue about the answer 
for the biology question.

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
