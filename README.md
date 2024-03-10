# Tutorial_16S_metabarcoding_analyses

## Overall workflow


## Package installation
First, some packages need to be installed.

For some of packages, the installation needs to go through BiocManager
```
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
```

For the whole workflow, two main packages are necessary: 

- [phyloseq](https://joey711.github.io/phyloseq/)
  
- [microbiome](https://microbiome.github.io/tutorials/)

```
BiocManager::install("phyloseq")
BiocManager::install("microbiome")
```



BiocManager::install("decontam") #for decontamination of the dataset

#for discriminant analyses
BiocManager::install("microbiomeMarker")


install.packages("randomcoloR")
install.packages("ggplot2")
install.packages("stringi")
install.packages("rlang")
install.packages(“vegan”)
install.packages("agricolae")
install.packages(“RVAideMemoire”)
install.packages(“ape”)
install.packages(“dplyr”)
install.packages("cowplot")
install.packages("stringi")
c
