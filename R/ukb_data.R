require( tidyverse )
require( rbgen )

ukb_filename  =  '../pheno/ukb21067.csv'
neale_filename  =  '../neale_files/21001.assoc.tsv.gz'
imp_sample_filename  =  '../imp/ukb1638_imp_chr1_v2_s487398.sample'
bgens_path  =  '../imp'
npcs  =  10

get_unrelated  =  function( input,
                            output ){
    unrelated_individuals  =  read_delim( input,
                                          delim = ' ',
                                          col_names = c( 'id', paste0( 'pc', 1:40 ) ) )

    save( unrelated_individuals, file = paste0( output, '.RData' ) )

    write_tsv( unrelated_individuals, paste0( output, '.tsv' ) )
}
# get_unrelated( '/data/sgg3/eleonora/projects/CNV_UKBB/QUALITY_SCORE/PHENO_FILES/europeanID_pca.txt.ID2keep',
#                'data/ukb_unrelated_pcs1-40' )

load( '../data/ukb_unrelated_pcs1-40.RData' )

phenotype_name  =  '21001-0.0'
names( phenotype_name )  =  'phenotype'

covariate_names  =  c( 'eid', '21003-0.0', '31-0.0', '22000-0.0' )

covariate_factor_names  =  c( '31-0.0', '22000-0.0' )

columns_to_keep  =  spec_csv( ukb_filename )

columns_to_keep  =  columns_to_keep$cols[ names( columns_to_keep$cols )
                                          %in% c( phenotype_name, covariate_names ) ]

ukb_data  =  read_csv( ukb_filename, col_types = do.call( cols_only, columns_to_keep ) ) %>%
    left_join( unrelated_individuals[ , 1:( npcs+1 ) ], ., by = c( 'id' = 'eid' ) ) %>%
    mutate_at( covariate_factor_names, as.factor ) %>%
    rename( !!!phenotype_name ) %>%
    filter( !is.na( phenotype_name ) )

inv_rank_normalise  =  function(x){
    res  =  rank( x, na.last = 'keep' )
    qnorm( res / ( length( res ) + 0.5 ) )
}

phenotype_residuals  =  lm( phenotype_name ~ .,
                            data = ukb_data[ , -1 ],
                            na.action = na.exclude ) %>%
    resid %>%
    tibble( id = ukb_data$id, phenotype = . ) # %>%
    # filter( !is.na( phenotype ) )

neale_assoc  =  read_tsv( neale_filename ) %>%
    filter( pval < 5e-8 ) %>%
    mutate( chromosome = as.numeric( str_match( variant, '^[0-9]+' ) ),
            position = as.numeric( str_match( variant, '^[0-9]+:([0-9]+)' )[ , 2 ] ) )

filter_within_chromosome  =  function( data, chr_number = NULL ){
    if (!is.null( chr_number )) {
        data  =  data %>%
            filter( chromosome == chr_number )
    }
    data  =  data %>%
        arrange( pval )
    final_data  =  data[ 0, ]

    while (nrow( data ) != 0) {
        final_data  =  bind_rows( final_data, data[ 1, ] )
        new_position  =  data$position[ 1 ]
        data  =  data[ -1, ] %>%
            filter( abs( position - new_position ) > 5e5 )
    }
    final_data
}

final_assoc  =  neale_assoc[ 0, ]

for (chr_number in 1:22) {
    final_assoc  =  final_assoc %>%
        bind_rows( filter_within_chromosome( neale_assoc, chr_number ) )
}

imp_samples  =  read_delim( imp_sample_filename, delim = ' ' )[ -1, ]

get_imputed_data  =  function( filename, rsids, samples_filter ){
    imputed_data  =  bgen.load( filename, rsids = rsids )$data[ , samples_filter, ]
    # imputed_data$samples  =  imputed_data$samples[ samples_filter ]
    imputed_data[ , , 2 ] + 2 * imputed_data[ , , 3 ]
}

get_all_imputed_data  =  function( path, rsids, samples_filter ){
    imputed_data  =  NULL
    for (filename in list.files( path, pattern = '[.]bgen$', full.names = TRUE )) {
        print( filename )
        if (is.null( imputed_data )) {
            imputed_data  =  get_imputed_data( filename, rsids, samples_filter )
        }
        else {
            imputed_data  =  get_imputed_data( filename, rsids, samples_filter ) %>%
                # as.tibble %>%
                rbind( imputed_data, . )
        }
    }
    imputed_data
}

samples_to_keep  =  imp_samples$ID_1 %in% phenotype_residuals$id

imputed_data  =  get_all_imputed_data( bgens_path, final_assoc$rsid, samples_to_keep )

phenotype_residuals  =  imp_samples[ , 'ID_1' ] %>%
    inner_join( phenotype_residuals, ., by = c( 'id' = 'ID_1' ) )

colnames( imputed_data )  =  phenotype_residuals$id



phenotype_residuals$phenotype  =  scale( phenotype_residuals$phenotype )

save( imputed_data, file = 'imputed_data.RData' )
save( phenotype_residuals, file = 'phenotype_residuals.RData' )

for (snp in 1:nrow( final_assoc )) {
    final_assoc[ snp, 'ukbeta' ]  =  lm( phenotype_residuals$phenotype ~ imputed_data[ snp, ] ) %>%
        coef %>%
        '['( 2 ) %>%
        unname
}

save( final_assoc, file = 'final_assoc.RData' )

grs  =  t( imputed_data ) %*% as.matrix( final_assoc[ , 'ukbeta', drop = FALSE ] ) %>%
    scale

# nlm( function( x ) IA_fit_<xyz>( x, phenotype_residuals$phenotype, grs ), rep( 0, <3-5> ) )






