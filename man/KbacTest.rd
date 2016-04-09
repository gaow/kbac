 \name{KbacTest}
 \alias{KbacTest}
 \title{KBAC statistic implementation}
 \description{
 This program implements the KBAC statistic in [Liu and Leal 2010]. It carries out case-control association testing for rare variants for whole exome association studies. Briefly, consider a gene of length n which harbors m rare variants. Genotype on the m variant sites & the disease status (case/control) are known for each individual. The program takes as input the m-site genotype and disease status (case/control) data files, and computes a p-value indicating the significance of association. In order to speed up permutation testing we use an "adaptive" approach to obtain p-values.
}

\usage{
KbacTest(casectrl.dat, alpha = NULL, num.permutation = NULL, quiet = T, maf.upper.bound = 1.0, alternative = 1)
}

\arguments{
      \item{casectrl.dat}{a numeric matrix with first column having disease status '0' or '1' and the rest columns codes the locus genotype as '0', '1', and '2'. DO NOT allow for missing data. }
      \item{alpha}{size of test, or the significant level. This feature will be useful in adaptive p-value calculation. If you do not want to use adaptive p-value, set alpha = 999 (or any number greater than 1.0).}
      \item{num.permutation}{number of permutations for p-value calculation.}
      \item{quiet}{when quiet = 0 the screen output would contain a summary of the KBAC test; otherwise only the p-value will be printed on screen.}
      \item{maf.upper.bound}{The upper bound of the MAF to be included in analysis. MAF is calculated based on observed sample. This can be arbitary although it is usually defined as 0.01 for analysis of rare variants.}
      \item{alternative}{Set alternative = 1 for test of deleterious variants, = 2 for test of both deleterious and protective variants. Please note that this is different from the "one/two-sided" definition in the KBAC paper.} 
}

\value{
	\item{pvalue}{the p-value of test.}
}

\details{
    ....
}

\author{Gao Wang | wangow@gmail.com}

\references{
Liu DJ,  Leal SM, 2010 A Novel Adaptive Method for the Analysis of Next-Generation Sequencing Data to Detect Complex Trait Associations with Rare Variants Due to Gene Main Effects and Interactions. PLoS Genet 6(10): e1001156. doi:10.1371/journal.pgen.1001156  
}

\examples{
# Load the package
library("KBAC")
?KbacTest

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
}
