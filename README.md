## KBAC Statistic Implementation
This is repository for the R implementation of [KBAC statistic](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1001156) for association testing.

### Methodology
This program implements the KBAC statistic in Liu and Leal 2010((Dajiang J. Liu and Suzanne M. Leal (2010). **A Novel Adaptive Method for the Analysis of Next-Generation Sequencing Data to Detect Complex Trait Associations with Rare Variants Due to Gene Main Effects and Interactions**. *PLoS Genetics*)). It carries out case-control association testing for rare variants for whole exome association studies. Briefly, consider a gene of length \(N\) which harbors \(M\) rare variants. Genotype on the \(M\) variant sites & the disease status (case/control) are known for each individual. The program takes as input the \(M\)-site genotype and disease status (case/control) data files, and computes a \(p\) value indicating the significance of association. In order to speed up permutation testing we use an "adaptive" approach to obtain \(p\) values.

### Latest Version

You can install the latest version of ''KBAC'' via the R script below:

```r
if (!require("devtools", character.only=TRUE, quietly=TRUE)) {
   install.packages("devtools")
}
devtools::install_github("gaow/kbac")
```

#### Dataset

*  Example dataset [[download](http://tigerwang.org/downloads/kbac-dataset.tar.gz)]

### ChangeLog

*  **2016.04.09** Move source code to [github](http://github.com/gaow/kbac). No updates made to the package.
*  **2013.02.22** Output the weights for the original genotype patterns, in response to user requests
*  **2011.05.27** GSL portability: I keep updating the package trying to port the GSL libraries into the program for a speedy hypergeometric routine. Complexities thus arise (various dependency problems) that I keep resolving. Please let me know if you fail to compile the package on your computer
*  **2011.05.18** Improved R/CPP interface
*  **2011.04.12** A more proper implementation of two-sided test
*  **2010.12.02** Initial release

## Documentation

### Note
In this R implementation

*  Only hyper-geometric kernel implemented
*  Only single candidate region analysis is supported (no covariates allowed)
	*  The length of candidate region reported in ''KBAC'' log is the number of variant sites analyzed (excluding sites with zero MAF or MAF above given threshold)
*  Allows for testing of both protective and deleterious mutations (set *alternative = 2*); Genotypes codings are 0 or 1 or 2. Invalid codings will be re-coded as 0 (wild-type). Phenotype codings are 0 for ctrls, 1 for cases
*  Genotype should consist only of the rare variants of interest -- synonymous, non-polymophic sites and common variants (MAF > 0.01) must be excluded;
*  Variant sites must be SNPs with numeric coding, no missing data allowed.

The purpose of this R packages to demonstration the KBAC methodology. It is not designed for analysis of real-world next-generation sequencing data, for which I suggest [VAT](http://varianttools.sf.net/VAT) as a resort.

### Input

Genotype for each locus are coded as *0* for wild type, *1* for heterozygotes and *2* for homozygotes, e.g.

| STATUS | M1 | M2 | M3 | M4 | M5 | M6 | M7 | M8 | M9 | M10 | M11 | M12 | 
| ------ | -- | -- | -- | -- | -- | -- | -- | -- | -- | --- | --- | --- | 
| 0      | 0  | 0  | 0  | 0  | 0  | 0  | 0  | 0  | 0  | 0   | 0   | 0   | 
| 0      | 0  | 0  | 0  | 0  | 1  | 0  | 0  | 0  | 0  | 0   | 0   | 0   | 
| 1      | 0  | 0  | 0  | 0  | 0  | 0  | 0  | 0  | 0  | 0   | 0   | 0   | 

### Example

#### Getting started
Install the package, start R, and type the following to load the package and read the usage

```r
library("KBAC")
?KbacTest
```

####  A quick demonstration
```r
casectrl.dat <- read.table("phengen.dat", skip = 1)
# Set parameters and use the KbacTest() function to obtain p-value
alpha <- 0.05
num.permutation <- 3000
quiet <- 1
alternative <- 1
maf.upper.bound <- 0.05
kbac.pvalue <- KbacTest(casectrl.dat, alpha, num.permutation, quiet, maf.upper.bound, alternative)
print(kbac.pvalue)
# To evaluate test at small alpha we need huge number of permutations. Adaptive approach is thus necessary.
kbac.pvalue <- KbacTest(casectrl.dat, 0.00001, 1000000, quiet, maf.upper.bound, alternative)
print(kbac.pvalue)
# Not using adaptive p-value calculation; will take longer time
kbac.pvalue <- KbacTest(casectrl.dat, 9, 1000000, quiet, maf.upper.bound, alternative)
print(kbac.pvalue)
```

##### Result

If you set ''quiet = F'' then you will see the screen output like this:

```
    Number of each unique individual genotype patterns (totaling 16 patterns excluding wildtype):
    3, 4, 4, 1, 29, 1, 16, 6, 40, 19, 1, 1, 51, 1, 1, 1,
    
    Unique genotype patterns weights (model 1):
    0.124812 0.687688 0.312312 1 0.98844 0.5 0.962177 0.656485 0.925138 0.821517 0.5 0.5 0.922296 0.5 0.5 1
```

Take the first genotype pattern for example. There are 3 individuals having this genotype pattern in the input data, **all of whom are controls**.

```r
> casectrl.dat[c(507,525,549),]
      V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21
    507  0  0  1  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0
    525  0  0  1  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0
    549  0  0  1  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0
      V22 V23 V24 V25 V26 V27 V28
    507   0   0   0   0   0   0   0
    525   0   0   0   0   0   0   0
    549   0   0   0   0   0   0   0
```

To see how the weight *0.124812* is tabulated, consider the ''phyper(x,m,n,k)'' notation in R where

*  x: number of cases (for model 1) having this genotype pattern
*  m: total number of cases
*  n: total number of ctrls
*  k: total number of samples having this genotype pattern

Then the weight is computed by ''phyper(0,1000,1000,3)'' which is *0.124812*

#### Compare with CMC method

An R script for the CMC method Li and Leal 2008((Bingshan Li and Suzanne M. Leal (2008). **Methods for Detecting Associations with Rare Variants for Common Diseases: Application to Analysis of Sequence Data**. *The American Journal of Human Genetics*)) via Fisher's exact test is provided for comparison purpose:

```r

arg <- commandArgs()

Cmc <- function(fn) {
    pgdata <- as.matrix(read.table(fn, as.is=T, skip = 1))
    y <- pgdata[,1];
    x `<- ((rowSums(as.matrix(pgdata[,-1])))>`0);
    m <- matrix(nrow=2,ncol=2);
    m[1,1] <- sum(x==1 & y==1);
    m[1,2] <- sum(x==0 & y==1);
    m[2,1] <- sum(x==1 & y==0);
    m[2,2] <- sum(x==0 & y==0);
    print(m);
    stat <- fisher.test(m)$p.value;
    return(stat);
}

Cmc(arg[3])
```

To use this script:

```bash
R --no-save phengen.dat < cmc.R
```
