library( rbgen )
library( parallel )
library( data.table )
library( dplyr )

source( 'estimate_gxe.R' )

get_ukb_imputed  =  function( path,
                              rsids,
                              samples_filter ){
    imputed_data  =  mclapply( list.files( path, pattern = '[.]bgen$', full.names = TRUE ),
               function( filename ) {
                   single_imputed  =  bgen.load( filename, rsids = rsids )$data[ , samples_filter, ]
                   single_imputed[ , , 2 ] + 2 * single_imputed[ , , 3 ]
               } )
    do.call( rbind, imputed_data )
}

get_ukb_sqc  =  function( sqc_filename,
                          geno_filename,
                          ids_to_remove = c(),
                          npcs = 10 ){
    sqc_columns  =  colnames( fread( sqc_filename, nrows = 0 ) )
    if( length( sqc_columns ) == 68 ) {
        sqc_data  =  fread( sqc_filename,
                            select = c( 4, 24:(25 + npcs) ))
    } else {
        sqc_data  =  fread( sqc_filename,
                            select = c( 2, 22:(23 + npcs) ))
    }

    sqc_data  =  bind_cols( fread( geno_filename,
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

get_betas_from_neale  =  function( neale_filename,
                                   variants_filename,
                                   threshold      = 5e-8,
                                   prune_distance = 5e5 ) {
    neale_stats  =  fread( paste( 'zcat', neale_filename ),
                           select = c( 'variant', 'low_confidence_variant', 'beta', 'pval' ) )
    neale_stats  =  filter( neale_stats,
                            low_confidence_variant == 'false',
                            pval < threshold )
    neale_stats  =  select( neale_stats, -low_confidence_variant )
    neale_stats  =  left_join( neale_stats,
                               fread( paste( 'zcat', variants_filename ),
                                      select = c( 'variant', 'rsid', 'chr', 'pos' ) ) )
    neale_stats  =  do.call( rbind,
             mclapply( 1:22,
                       function( chromosome ){
                           data  =  neale_stats[ neale_stats$chr == chromosome, ]
                           data  =  arrange( data, pval )
                           final_data  =  data[ 0, ]

                           while (nrow( data ) != 0) {
                               final_data  =  bind_rows( final_data, data[ 1, ] )
                               new_position  =  data$pos[ 1 ]
                               data  =  filter( data, abs( pos - new_position ) > prune_distance )
                           }
                           final_data
                       } ) )
    neale_stats[ , c( 'rsid', 'beta' ) ]
}

ukb_gxe_interaction  =  function( phenotype_name,
                                  ukb_filename,
                                  bgens_path,
                                  snps,
                                  covariate_names        = NULL,
                                  covariate_factor_names = NULL,
                                  correct_age2_sex       = TRUE,
                                  sample_ids             = NULL,
                                  sqc_filename           = NULL,
                                  geno_filename          = NULL,
                                  imp_sample_filename    = list.files( bgens_path,
                                                                       pattern = '[.]sample$',
                                                                       full.names = TRUE )[1],
                                  ids_to_remove          = NULL,
                                  npcs                   = 10,
                                  betas                  = NULL ){
    covariate_names  =  unique( c( covariate_names, covariate_factor_names ) )

    if (is.null( sample_ids )) {
        sample_ids  =  get_ukb_sqc( sqc_filename, geno_filename, ids_to_remove, npcs )
    }

    if (correct_age2_sex) {
        covariate_names  =  c( covariate_names, '21022-0.0', '31-0.0' )
    }

    ukb_data  =  do.call( cbind,
                          lapply( ukb_filename,
                                  function(x) fread( ukb_filename,
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

    imputed_data  =  get_ukb_imputed( bgens_path, snps$rsid, samples_to_keep )

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

    estimate_gxe( phenotype_residuals$phenotype, grs )
}







