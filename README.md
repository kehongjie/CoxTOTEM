# CoxTOTEM
A Cox model based two-stage variable selection method for the detection of survival associated biomarkers with multiple genomic studies. In the first stage, it performs sure screening with multiple studies in high-dimensional Cox model, and in the second stage, it penalizes the partial likelihood with a group lasso penalty to select the final set of variables in all studies simultaneously. 

# Installation
To install the `CoxTOTEM` package, you will first need to install `devtools` package and then execute the following code: 
```
devtools::install_github('kehongjie/CoxTOTEM')
```

# Data
The `PanGyn.RData` in the `data` folder has RNA-seq data and clinical data containing the survival outcomes of the five Pan-Gyn cancer types, including high-grade serous ovarian cystadenocarcinoma (OV), uterine corpus endometrial carcinoma (UCEC), cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC), uterine carcinosarcoma (UCS), and invasive breast carcinoma (BRCA). The data comes from The Cancer Genome Atlas (TCGA) projects. To access to the data, simply run the following code: 
````
library(CoxTOTEM)
data(PanGyn)
``` 
