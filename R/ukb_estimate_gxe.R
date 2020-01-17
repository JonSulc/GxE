#' Estimate contribution of GRSxE to variance in outcome phenotype in UK Biobank
#'
#' @name ukb_estimate_gxe
#'
#' @param phenotype_name A string containing the column name of the phenotype to
#'   use as outcome, e.g. '21001-0.0'
#' @param ukb_filename A string or character vector containing the path to the
#'   file(s) containing the UK Biobank data
#' @param bgens_path path to the folder containing the \code{*.bgen} files with
#'   the genetic data
#' @param snps a character vector of SNP RSIDs or a data.frame-type structure
#'   containing a column named 'rsid' with SNP RSIDs, as well as an optional
#'   column named 'beta' containing the effect estimates of the alternate allele
#'   on the outcome
#' @param covariate_names a string or character vector containing the names of
#'   any covariates to adjust for before estimating GxE
#' @param covariate_factor_names a string or character vector containing the
#'   names of covariates to be considered as factors when adjusting for them
#' @param correct_age2_sex logical indicating whether to adjust for age, age^2,
#'   and sex (in addition to any other covariates)
#' @param sample_ids a data.frame containing the individual ids under the column
#'   name 'eid'. Any other included columns (e.g. genetic PCs) will be adjusted
#'   for. Required if \code{sqc_filename} or \code{fam_filename} are missing
#' @param sqc_filename path to the UK Biobank sqc file (e.g. ukb_sqc_v2.txt).
#'   Required if \code{sample_ids} is missing
#' @param fam_filename path to any \code{*.fam} file containing the sample ids
#'   for the values in \code{sqc_filename}. Required if \code{sample_ids} is
#'   missing
#' @param imp_sample_filename path to a \code{*.sample} file containing the
#'   sample ids for the genetic data. The default is the first \code{*.sample}
#'   in the \code{bgens_path} folder
#' @param ids_to_remove vector containing the ids of any samples to remove
#'   before analysis
#' @param npcs number of genetic principal components to correct for before
#'   analysis
#' @param betas vector of effects of SNPs on outcome. Ignored if \code{snps}
#'   contains a column named 'beta'. Must be in the same order as \code{snps}
#' @param sim_num Number of permutations for bootstrap and fake GRSs
#'
#' @importFrom data.table fread
#' @importFrom dplyr bind_cols mutate_at tibble
#' @importFrom parallel mclapply
#' @importFrom rbgen bgen.load
#' @export
#'
#' @return \code{ukb_estimate_gxe} returns a list containing parameter estimates
#'   for alpha1, alpha2, beta, and gamma (\code{xopt}), their standard error
#'   (\code{SExopt}), and the associated p-values (\code{Pxopt}). The
#'   corresponding estimates for the fake GRS are stored in the \code{xopt0},
#'   \code{SExopt0}, and \code{Pxopt0}, respectively.
#'
#'   In addition, the \code{Xopt} and \code{Xopt0} matrices contain the
#'   estimates for each bootstrap and fake GRS.
#'
#'   The \code{tdiff} contains the t-statistic for the difference between the
#'   real data estimates and those from the fake GRS.
#' @examples
#' \dontrun{
#'   library( GxE )
#'   # Load GWAS results from the Neale lab
#'   snps = get_betas_from_neale( neale_filename = '21001_irnt.gwas.imputed_v3.both_sexes.tsv.gz',
#'                                variants_filename = 'variants.tsv.gz' )
#'
#'   # UK Biobank files not provided
#'   gxe  =  ukb_gxe_interaction( phenotype_name = '21001-0.0',
#'                                ukb_filename   = 'uk_biobank/pheno/ukb21067.csv',
#'                                bgens_path     = 'uk_biobank/imp',
#'                                snps           = snps,
#'                                sqc_filename   = 'uk_biobank/geno/ukb_sqc_v2.txt',
#'                                fam_filename   = 'uk_biobank/plink/ukb1638_cal_chr1_v2_s488366.fam' )
#'
#'   print( gxe )
#' }

library( rbgen )
library( parallel )
library( data.table )
library( dplyr )

.get_ukb_imputed  =  function( path,
                               rsids,
                               samples_filter ){
    imputed_data  =  mclapply( list.files( path, pattern = '[.]bgen$', full.names = TRUE ),
               function( filename ) {
                   single_imputed  =  bgen.load( filename, rsids = rsids )$data[ , samples_filter, ]
                   single_imputed[ , , 2 ] + 2 * single_imputed[ , , 3 ]
               } )
    do.call( rbind, imputed_data )
}

.get_ukb_sqc  =  function( sqc_filename,
                           fam_filename,
                           ids_to_remove = c(),
                           npcs = 10 ){
    sqc_columns  =  colnames( fread( sqc_filename, nrows = 1 ) )
    if( length( sqc_columns ) == 68 ) {
        sqc_data  =  fread( sqc_filename,
                            select = c( 4, # Batch
                                        24,
                                        25:(25 + npcs) ))
    } else {
        sqc_data  =  fread( sqc_filename,
                            select = c( 2, # Batch
                                        22,
                                        23:(23 + npcs) ))
    }

    sqc_data  =  bind_cols( fread( fam_filename,
                                   select = c( 1 ) ),
                            sqc_data )
    colnames( sqc_data )  =  c( 'eid',
                                'batch',
                                'white',
                                'unrelated',
                                paste0( 'pc', 1:npcs ) )

    sqc_data  =  sqc_data[ !(sqc_data$eid %in% ids_to_remove) &
                               sqc_data$white == 1 &
                               sqc_data$unrelated == 1, ]
    sqc_data  =  sqc_data[ , -c(3:4) ]

    sqc_data$batch = as.factor(sqc_data$batch)
    sqc_data
}

ukb_estimate_gxe  =  function( phenotype_name,
                               ukb_filename,
                               bgens_path,
                               snps,
                               covariate_names        = NULL,
                               covariate_factor_names = NULL,
                               correct_age2_sex       = TRUE,
                               sample_ids             = NULL,
                               sqc_filename           = NULL,
                               fam_filename           = NULL,
                               imp_sample_filename    = list.files( bgens_path,
                                                                    pattern = '[.]sample$',
                                                                    full.names = TRUE )[1],
                               ids_to_remove          = NULL,
                               npcs                   = 10,
                               betas                  = NULL,
                               sim_num                = 100,
                               ... ){
    covariate_names  =  unique( c( covariate_names, covariate_factor_names ) )

    if (is.null( sample_ids )) {
        sample_ids  =  .get_ukb_sqc( sqc_filename, fam_filename, ids_to_remove, npcs )
    }

    if (correct_age2_sex) {
        covariate_names  =  c( covariate_names, '21022-0.0', '31-0.0' )
    }

    ukb_data  =  do.call( cbind,
                          lapply( ukb_filename,
                                  function(filename) fread( filename,
                                                            select = c( 'eid',
                                                                        phenotype_name,
                                                                        covariate_names ) ) ) )
    ukb_data  =  ukb_data[ , c( 'eid', phenotype_name, covariate_names ), with = FALSE ]

    colnames( ukb_data )[2]  =  'phenotype'
    ukb_data  =  ukb_data[ , lapply( .SD, as.numeric ) ]
    if (correct_age2_sex) {
        ukb_data$age2  =  ukb_data$'21022-0.0'^2
    }

    ukb_data  =  left_join( sample_ids, ukb_data )
    ukb_data  =  ukb_data[ !is.na(ukb_data$phenotype), ]
    if (!is.null( covariate_factor_names )) {
        ukb_data  =  mutate_at( ukb_data, covariate_factor_names, as.factor )
    }

    phenotype_residuals  =  lm( phenotype ~ .,
                                data = ukb_data[ , -1 ],
                                na.action = na.exclude )$residuals
    phenotype_residuals  =  tibble( eid = ukb_data$eid,
                                    phenotype = c( scale( phenotype_residuals ) ) )

    bgens_path  =  sub('/$', '', bgens_path)
    imp_samples  =  fread( imp_sample_filename )
    imp_samples  =  imp_samples[ -1, ]
    samples_to_keep  =  imp_samples$ID_1 %in% phenotype_residuals$eid

    if (is.null( dim( snps ))){
        snps  =  tibble( rsid = snps )
    }

    imputed_data  =  .get_ukb_imputed( bgens_path, snps$rsid, samples_to_keep )

    phenotype_residuals  =   left_join( phenotype_residuals,
                                        imp_samples[ , 'ID_1', drop = FALSE ],
                                        by = c( 'eid' = 'ID_1' ))

    if (any( duplicated( rownames( imputed_data ) ) )) {
        snps = snps[ !(snps$rsid %in% rownames( imputed_data )[ duplicated( rownames( imputed_data ) ) ]), ]
    }
    snps  =  snps[ snps$rsid %in% rownames( imputed_data ), ]
    imputed_data  =  imputed_data[ snps$rsid, ]

    if (!( 'beta' %in% colnames( snps ) )){
        snps$beta  =  betas
    }

    grs  =  scale ( t( imputed_data ) %*% as.matrix( snps$beta ) )

    estimate_gxe( phenotype_residuals$phenotype, grs, sim_num, ... )
}







