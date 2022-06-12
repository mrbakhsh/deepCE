# Inference of Interactome and Underying Protein Complexes from Co-fractionation MS-Data via Deep Learning
Co-fractionation mass spectrometry (CF-MS) has emerged as a powerful approach 
for cell-wide identification of protein complexes that performs nearly all
biological processes and cellular activities. While computational tools based 
on traditional machine learning is widely applied to analyze protein complexes, 
it remains unclear if it is the most suitable choice for analyzing CF/MS 
dataset. Here, we introduce deepCE, a R/Bioconductor package that employs deep 
learning approaches for automated scoring of co-fractionation mass spectrometry 
(CF-MS) and sophisticated clustering procedure for network inference of 
underlying complexes. 


## Installation
To install the stable version from Bioconductor:

```r
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") 
}
BiocManager::install("deepCE")
```
To view documentation for the version of this package installed in your system, 
start R and enter:

```r
browseVignettes("deepCE")
```

To install the development version in `R`, run:
  
```r
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}
devtools::install_github("mrbakhsh/deepCE")
```

## Required Software
Although, this is a R package, users will need to 
install [Python](https://www.python.org/downloads/) as well. In addition, 
tensorflow and keras should be installed in Python and 
they should be connected to R, before installing this package.
