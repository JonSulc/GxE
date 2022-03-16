
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
the polygenic score, `grs`. This returns a list containing a separate
list for the real data (`real_data`) and the simulated GRS (`fake_grs`).
Each of these lists contain a vector, `coefficients`, with the average
estimates obtained for each of the 4 parameters, namely `alpha1` (linear
genetic effect), `alpha2` (quadratic genetic effect), `beta`
(environmental effect), and `gamma` (GxE interaction effect), as well as
vectors for the corresponding standard errors (`se`) and p-values (`p`)
and a matrix with the individual estimates of each bootstrap permutation
(`individual_estimates`).

If `simulate_phenotype` was used, the returned list `fake_phenotype`
contains similar information for the simulated phenotype, without
`alpha2` as quadratic genetic effects are accounted for in the phenotype
simulation. Additionally, this list contains the `skewness` and
`kurtosis` values which most closely match the data, the root mean
square difference between the real and fake phenotypes (`rms_diff`), the
main linear genetic effect (`alp`), and the estimates of the genetic
effects from linear to seventh power (`fY_coefficients`).

In addition, the `t_real_fgrs` contains the t-statistic for the
difference between the results from the real and fake GRS.

``` r
library(GxE)

gxe  =  estimate_gxe( y, grs )
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
