
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GxE

The goal of GxE is to enable easy estimation of gene-environment (GxE)
interaction from individual-level polygenic scores (GRS) and outcome
phenotypes alone.

## Installation

You can install the latest version of GxE from
[GitHub](https://github.com/) with:

``` r
install.packages( "devtools" )
library( devtools )
install_github( "zkutalik/GRSxE_software", subdir = "Rcode" )
```

## Example

The main function is `estimate_gxe`, which requires the outcome
phenotype of interest `y` (corrected for any relevant covariates) and
the polygenic score, `GRS`. This returns a list of results.

`xopt` and `xopt0` are vectors, each containing 4 parameters, namely
`alpha1` (linear genetic effect), `alpha2` (quadratic genetic effect),
`beta` (environmental effect), and `gamma` (GxE interaction effect).
While the `xopt` vector contains these parameters estimated in the real
data, `xopt0` contains parameters estimated using a counterfeit GRS.

`SExopt`, `SExopt0`, `Pxopt`, and `Pxopt0` are similarly structured but
contain the respective standard errors and p-values obtained from
bootstrapping.

`tdiff` contains the t-statistic for the difference between the results
from the real and the counterfeit GRS.

``` r
library(GxE)

gxe  =  estimate_gxe( y, GRS )
print( gxe )
```

In addition, the `ukb_estimate_gxe` is a helper function to allow easy
analysis of UK Biobank data (requires access to a local copy of
individual-level UK Biobank genetic and phenotypic data). It allows
covariates to be specified and corrected for before estimating GxE (age,
age^2, sex, and the top 10 genetic PCs are used as covariates by
default).

``` r
gxe  =  ukb_gxe_interaction( phenotype_name = '21001-0.0',
                             ukb_filename   = 'uk_biobank/pheno/ukb21067.csv',
                             bgens_path     = 'uk_biobank/imp',
                             snps           = snps,
                             sqc_filename   = 'uk_biobank/geno/ukb_sqc_v2.txt',
                             fam_filename   = 'uk_biobank/plink/ukb1638_cal_chr1_v2_s488366.fam' )
print( gxe )
```

`get_betas_from_neale` is another helper function included to
automatically extract and prune SNPs and their effect sizes from files
as those provided by the [Neale
Lab](http://www.nealelab.is/uk-biobank/).

``` r
snps = get_betas_from_neale( neale_filename = '21001_irnt.gwas.imputed_v3.both_sexes.tsv.gz',
                             variants_filename = 'variants.tsv.gz' )
```
