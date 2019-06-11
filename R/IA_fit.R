require( tidyverse )
require( rbgen )

IA_fit  =  function( abc, y, g ){
    if (is.null( dim( abc ) )) {
        abc  =  array( abc, dim = c( 3, length( abc ) / 3 ) )
    }
    if (is.null( dim( y ) )) {
        y  =  array( y, dim = c( length( y ), 1 ) )
    }
    if (is.null( dim( g ) )) {
        g  =  array( g, dim = c( length( g ), 1 ) )
    }

    a  =  abc[ 1, , drop = FALSE ]
    b  =  abc[ 2, , drop = FALSE ]
    c  =  abc[ 3, , drop = FALSE ]

    s2  =  1 - a^2 - b^2 - c^2
    d1  =  array( 1, dim = c( length( y ), 1 ) )
    d2  =  array( 1, dim = c( 1, length( a ) ) )

    v0  =  0.5 * colSums( ( y %*% d2 - g %*% a )^2 / ( ( d1 %*% b + g %*% c )^2 + d1 %*% s2 )
                      + suppressWarnings( log( ( d1 %*% b + g %*% c )^2 + d1 %*% s2 ) ) )
    v  =  0.5 * colSums( ( y %*% d2 )^2 )

    sel  =  s2 > 0
    v[ sel ]  =  v0[ sel ]
    v
}

IA_fit_GEcorr  =  function( abc, y, g ){
    if (is.null( dim( abc ) )) {
        abc  =  array( abc, dim = c( 4, length( abc ) / 4 ) )
    }
    if (is.null( dim( y ) )) {
        y  =  array( y, dim = c( length( y ), 1 ) )
    }
    if (is.null( dim( g ) )) {
        g  =  array( g, dim = c( length( g ), 1 ) )
    }

    a  =  abc[ 1, ]
    b  =  abc[ 2, ]
    c  =  abc[ 3, ]
    d  =  abc[ 4, ]

    s2  =  1 - a^2 - b^2 - c^2 * ( 1 + d^2 ) - 2 * a * b * d
    d1  =  array( 1, dim = c( length( y ), 1 ) )
    d2  =  array( 1, dim = c( 1, length( a ) ) )
    g2  =  g^2

    av   =  a + b * d
    a2v  =  c * d
    bv   =  b * suppressWarnings( sqrt( 1 - d^2 ) )
    cv   =  c * suppressWarnings( sqrt( 1 - d^2 ) )

    v0  =  0.5 * colSums( ( y %*% d2 - g %*% av - ( g2 - 1 ) %*% a2v )^2
                          / ( ( d1 %*% bv + g %*% cv )^2 + d1 %*% s2 )
                          + suppressWarnings( log( ( d1 %*% bv + g %*% cv )^2 + d1 %*% s2 ) ) )
    v   =  0.5 * colSums( ( y %*% d2 )^2 )

    sel  =  s2 > 0 & s2 < 1 & a^2 < 1 & b^2 < 1 & c^2 < 1 & d^2 < 1
    v[ sel ]  =  v0[ sel ]
    v
}

IA_fit_pow  =  function( abc, y, g ){
    if (is.null( dim( abc ) )) {
        abc  =  array( abc, dim = c( 4, length( abc ) / 4 ) )
    }
    if (is.null( dim( y ) )) {
        y  =  array( y, dim = c( length( y ), 1 ) )
    }
    if (is.null( dim( g ) )) {
        g  =  array( g, dim = c( length( g ), 1 ) )
    }

    a   =  abc[ 1, , drop = FALSE ]
    b   =  abc[ 2, , drop = FALSE ]
    c   =  abc[ 3, , drop = FALSE ]
    d2  =  abc[ 4, , drop = FALSE ]
    d1  =  suppressWarnings( sqrt( 1 - 2 * d2^2 ) )

    s2  =  1 -
        (d1 * a)^2 - 2 * ( d2 * a^2 )^2 -
        ( d1 * b )^2 -
        ( d1 * c + 2 * d2 * a * b )^2 -
        3 * ( 2 * d2 * a * c )^2 -
        2 * ( d1 * c + 2 * d2 * a * b ) * ( 2 * d2 * a * c)
    s2  =  s2[ 1, 1 ]
    h1  =  array( 1, dim = c( length( y ), 1 ) )
    h2  =  array( 1, dim = c( 1, length( a ) ) )
    g2  =  g^2

    c_g    =  d1 * a
    c_g2   =  d2 * a^2
    c_e    =  d1 * b
    c_ge   =  d1 * c + 2 * d2 * a * b
    c_g2e  =  2 * d2 * a * c

    yy  =  y %*% h2 - ( g %*% c_g + ( g2 - 1 ) %*% c_g2 )
    sig2  =  ( h1 %*% c_e + g %*% c_ge + g2 %*% c_g2e )^2 + s2

    v0  =  0.5 * colSums( yy^2 / sig2
                          + suppressWarnings( log( sig2 ) ) )
    v  =  0.5 * colSums( ( y %*% h2 )^2 )

    sel  =  s2 > 0 & a^2 + b^2 + c^2 <= 1 & d2^2 <= 0.5
    v[ sel ]  =  v0[ sel ]
    v
}

IA_fit_GEcorr_pow  =  function( abc, y, g ){
    if (is.null( dim( abc ) )) {
        abc  =  array( abc, dim = c( 5, length( abc ) / 5 ) )
    }
    if (is.null( dim( y ) )) {
        y  =  array( y, dim = c( length( y ), 1 ) )
    }
    if (is.null( dim( g ) )) {
        g  =  array( g, dim = c( length( g ), 1 ) )
    }

    a   =  abc[ 1, , drop = FALSE ]
    b   =  abc[ 2, , drop = FALSE ]
    c   =  abc[ 3, , drop = FALSE ]
    d   =  abc[ 4, , drop = FALSE ]
    d2  =  abc[ 5, , drop = FALSE ]
    d1  =  suppressWarnings( sqrt( 1 - 2 * d2^2 ) )

    h1  =  array( 1, dim = c( length( y ), 1 ) )
    h2  =  array( 1, dim = c( 1, length( a ) ) )
    g2  =  g^2
    g3  =  g^3
    g4  =  g^4

    c_g    =  d1 * (a + b * d)
    c_g2   =  d1 * c * d + d2 * (a + b * d)^2
    c_g3   =  2 * d2 * (a + b * d) * c * d
    c_g4   =  d2 * (c * d)^2
    c_e    =  d1 * b * suppressWarnings( sqrt( 1 - d^2 ) )
    c_ge   =  (d1 * c + 2 * d2 * (a + b * d) * b) * suppressWarnings( sqrt( 1 - d^2 ) )
    c_g2e  =  2 * d2 * (a + 2 * b * d) * c * suppressWarnings( sqrt( 1 - d^2 ) )
    c_g3e  =  2 * d2 * c^2 * d * suppressWarnings( sqrt( 1 - d^2 ) )
    c_0    =  -( c_g2 + 3 * c_g4 )

    s2  =  1 -
        ( c_g^2
          + 3   * c_g2^2
          + 15  * c_g3^2
          + 105 * c_g4^2
          +       c_e^2
          +       c_ge^2
          + 3   * c_g2e^2
          + 15  * c_g3e^2
          + 2   * c_e     * c_g2e
          + 6   * c_ge    * c_g3e
          + 6   * c_g     * c_g3
          + 30  * c_g2    * c_g4
          - (c_g2 + 3 * c_g4)^2 )

    yy  =  (  y %*% h2
            - h1 %*% c_0
            - g  %*% c_g
            - g2 %*% c_g2
            - g3 %*% c_g3
            - g4 %*% c_g4 )

    vy  =  (  h1 %*% c_e
            + g  %*% c_ge
            + g2 %*% c_g2e
            + g3 %*% c_g3e )^2 + h1 %*% s2

    v0  =  0.5 * colSums( yy^2 / vy
                          + suppressWarnings( log( vy ) ) )
    v  =  0.5 * colSums( ( y %*% h2 )^2 )

    sel  =  s2 > 0 &
            s2 < 1 &
            d2^2 < 0.5 &
            a^2 < 1 &
            b^2 < 1 &
            c^2 < 1 &
            d^2 < 1 &
            (a + b * d)^2 + b^2 * (1 - d^2) + 3 * (c * d)^2 + c^2 * (1 - d^2) < 1
    v[ sel ]  =  v0[ sel ]
    v
}

IA_fit_gam0_pow  =  function( abc, y, g ){
    if (is.null( dim( abc ) )) {
        abc  =  array( abc, dim = c( 3, length( abc ) /3 ) )
    }

    a   =  abc[ 1, , drop = FALSE ]
    b   =  abc[ 2, , drop = FALSE ]
    # c   =  0
    # d   =  0
    d2  =  abc[ 3, , drop = FALSE ]
    d1  =  suppressWarnings( sqrt( 1 - 2 * d2^2 ) )

    h1  =  array( 1, dim = c( length( y ), 1 ) )
    h2  =  array( 1, dim = c( 1, length( a ) ) )
    g2  =  g^2
    g3  =  g^3
    g4  =  g^4

    c_g    =  d1 * a
    c_g2   =  d2 * a^2
    # c_g3   =  0
    # c_g4   =  0
    c_e    =  d1 * b
    c_ge   =  2 * d2 * a * b
    # c_g2e  =  0
    # c_g3e  =  0
    c_0    =  -c_g2

    s2  =  1 -
        ( c_g^2
          + 3   * c_g2^2
          +       c_e^2
          +       c_ge^2
          - c_g2^2 )

    yy  =  (  y %*% h2
              - h1 %*% c_0
              - g  %*% c_g
              - g2 %*% c_g2 )

    vy  =  (  h1 %*% c_e
              + g  %*% c_ge )^2 + h1 %*% s2

    v0  =  0.5 * colSums( yy^2 / vy
                          + suppressWarnings( log( vy ) ) )
    v  =  0.5 * colSums( ( y %*% h2 )^2 )

    sel  =  s2 > 0 &
        s2 < 1 &
        d2^2 < 0.5 &
        a^2 < 1 &
        b^2 < 1 &
        a^2 + b^2 < 1
    v[ sel ]  =  v0[ sel ]
    v
}

bic  =  function( n_parameters,
                  sample_size,
                  loglikelihood ){
    2 * n_parameters * log( sample_size ) - 2 * loglikelihood
}

aic  =  function( n_parameters,
                  loglikelihood ){
    2 * n_parameters - 2 * loglikelihood
}

se_bootstrap  =  function( ia_function,
                           n_parameters,
                           phenotype,
                           genotype,
                           n_simulations = 100 ){
    abc  =  array( NA, dim = c( n_simulations, n_parameters ) )
    for (j in 1:n_simulations) {
        ix  =  sample( length( phenotype ), length( phenotype ), replace = TRUE )
        abc[ j, ]  =  nlm( function( x ) do.call( ia_function, list( x, phenotype[ ix, ], genotype[ ix, ] ) ),
                           rep( 0, n_parameters ) ) %>%
            '$'( 'estimate' )
    }
    apply( abc, 2, sd )
}

# Wrapper

single_gxe  =  function( ia_function,
                         n_parameters,
                         phenotype,
                         genotype,
                         n_simulations = 100 ){
    results  =  nlm( function( x ) do.call( ia_function, list( x, phenotype, genotype ) ),
                     rep( 0, n_parameters ) )
    # return(results)
    results$se  =  se_bootstrap( ia_function, n_parameters, phenotype, genotype, n_simulations )
    results$pvalues  =  2 * pnorm( -abs( results$estimate ) / results$se )
    results$aic  =  aic( length( results$estimate ), -results$minimum )
    results
}

gxe_interaction  =  function( phenotype,
                              genotype,
                              n_simulations = 100 ){
    if (is.null( dim( phenotype ) )) {
        phenotype  =  array( phenotype, dim = c( length( phenotype ), 1 ) )
    }
    if (is.null( dim( genotype ) )) {
        genotype  =  array( genotype, dim = c( length( genotype ), 1 ) )
    }

    models  =  list( gxe          =  c( 'IA_fit',            3 ),
                     pow          =  c( 'IA_fit_gam0_pow',   3 ),
                     gxe_cor      =  c( 'IA_fit_GEcorr',     4 ),
                     gxe_pow      =  c( 'IA_fit_pow',        4 ),
                     gxe_cor_pow  =  c( 'IA_fit_GEcorr_pow', 5 ) )

    results  =  list()

    for (model_name in names( models )) {
        results[[ model_name ]]  =  single_gxe( models[[ model_name ]][ 1 ],
                                                models[[ model_name ]][ 2 ],
                                                phenotype,
                                                genotype,
                                                n_simulations )
    }

    results$simple      =  lm( phenotype ~ genotype )
    results$simple$aic  =  AIC( results$simple )
    results
}

get_imputed_data  =  function( filename, rsids, samples_filter ){
    imputed_data  =  bgen.load( filename, rsids = rsids )$data[ , samples_filter, ]
    imputed_data[ , , 2 ] + 2 * imputed_data[ , , 3 ]
}

get_all_imputed_data  =  function( path, rsids, samples_filter ){
    imputed_data  =  NULL
    for (filename in list.files( path, pattern = '[.]bgen$', full.names = TRUE )) {
        print( filename )
        if (is.null( imputed_data )) {
            imputed_data  =  get_imputed_data( filename, rsids, samples_filter )
        } else {
            imputed_data  =  get_imputed_data( filename, rsids, samples_filter ) %>%
                rbind( imputed_data, . )
        }
    }
    imputed_data
}

get_unrelated_pcs  =  function( sqc_filename,
                                geno_filename,
                                ids_to_remove = c(),
                                npcs = 10 ){
    sqc_columns  =  spec_delim( sqc_filename, delim = ' ', col_names = FALSE, guess_max = 1 )
    if( length( sqc_columns$cols ) < 68 ){
        sqc_data  =  read_delim( sqc_filename, delim = ' ',
                                 col_types = paste( c( '_',
                                                       'c',
                                                       rep( '_', 19 ),   # Ignore
                                                       'll',             # White british, unrelated
                                                       rep( 'd', npcs ),   # Genetic principal components
                                                       rep( '_', 43 - npcs  ) ), # Ignore
                                                    collapse = '' ) )
    } else {
        sqc_data  =  read_delim( sqc_filename,
                                 delim = ' ',
                                 col_names = FALSE,
                                 col_types = paste( c( rep( '_', 3 ),
                                                       'c',
                                                       rep( '_', 19 ),   # Ignore
                                                       'll',             # White british, unrelated
                                                       rep( 'd', npcs ),   # Genetic principal components
                                                       rep( '_', 43 - npcs  ) ), # Ignore
                                                    collapse = '' ) )
    }

    # sqc_data  =  sqc_data %>%
        # filter_at( 2:3, all_vars( . ) ) %>%
        # select( c( 1:npcs + 3 ) )

    sqc_data  =  read_delim( geno_filename,
                             delim = ' ',
                             col_names = FALSE,
                             col_types = paste( c( 'i',              # Subject ID
                                                   rep( '_', 5  ) ), # Ignore
                                                collapse = '' ) ) %>%
        bind_cols( sqc_data ) %>%
        filter( !(.[[ 1 ]] %in% ids_to_remove) ) %>%
        filter_at( 3:4, all_vars( . ) ) %>%
        select( -c( 3:4 ) )
    sqc_data[[ 2 ]]  =  as.factor( sqc_data[[ 2 ]] )
    colnames( sqc_data )  =  c( 'eid', 'batch', paste0( 'pc', 1:npcs ) )
    sqc_data
}

ukb_gxe_interaction  =  function( phenotype_name,
                                  ukb_filename,
                                  bgens_path,
                                  snps,
                                  covariate_names,
                                  covariate_factor_names,
                                  sample_ids = NULL,
                                  sqc_filename = NULL,
                                  geno_filename = NULL,
                                  ids_to_remove = c(),
                                  npcs = 10,
                                  covariate_filename = NULL,
                                  betas = NULL ){
    names( phenotype_name )  =  'phenotype'

    if (!is.null( sqc_filename )) {
        if (is.null( geno_filename ))
            stop( "Using the SQC file requires providing a *.fam file providing the sample IDs" )
        if (!is.null( sample_ids )) {
            sample_ids  =  get_unrelated_pcs( sqc_filename, geno_filename, ids_to_remove, npcs ) %>%
                left_join( tibble( eid = sample_ids ) )
        } else {
            sample_ids  =  get_unrelated_pcs( sqc_filename, geno_filename, ids_to_remove, npcs )
        }
    }

    if (is.null( covariate_filename )) {
        columns_to_keep  =  spec_csv( ukb_filename )
        columns_to_keep  =  columns_to_keep$cols[ names( columns_to_keep$cols )
                                                  %in% c( 'eid', phenotype_name, covariate_names ) ]
    } else {
        columns_to_keep  =  spec_csv( ukb_filename )
        columns_to_keep  =  columns_to_keep$cols[ names( columns_to_keep$cols )
                                                  %in% c( 'eid', phenotype_name ) ]
        covariate_columns  =  spec_csv( covariate_filename )
        covariate_columns  =  covariate_columns$cols[ names( covariate_columns$cols )
                                                      %in% c( 'eid', covariate_names ) ]
    }

    ukb_data  =  read_csv( ukb_filename, col_types = do.call( cols_only, columns_to_keep ) ) %>%
        rename( !!!phenotype_name )

    if (!is.null( sample_ids )) {
        ukb_data  =  ukb_data %>%
            left_join( as.tibble( sample_ids ), . )
    }

    if (!is.null( covariate_filename )) {
        ukb_data  =  read_csv( covariate_filename, col_types = do.call( cols_only, covariate_columns ) ) %>%
            left_join( ukb_data, . )
    }

    ukb_data  =  ukb_data %>%
        filter( !is.na( phenotype ) ) %>%
        mutate_at( covariate_factor_names, as.factor )

    phenotype_residuals  =  lm( phenotype ~ .,
                                data = ukb_data[ , -1 ],
                                na.action = na.exclude ) %>%
        resid %>%
        scale %>%
        c %>%
        tibble( eid = ukb_data$eid, phenotype = . )

    bgens_path  =  sub('/$', '', bgens_path)
    imp_samples  =  list.files( bgens_path, pattern = '[.]sample$', full.names = TRUE )[ 1 ] %>%
        read_delim( ., delim = ' ' )
    imp_samples  =  imp_samples[ -1, ]
    samples_to_keep  =  imp_samples$ID_1 %in% phenotype_residuals$eid

    if (is.null( dim( snps )))
        snps   =  tibble( rsid = snps )

    imputed_data  =  get_all_imputed_data( bgens_path, snps$rsid, samples_to_keep )

    phenotype_residuals  =  imp_samples[ , 'ID_1' ] %>%
        inner_join( phenotype_residuals, ., by = c( 'eid' = 'ID_1' ))

    # colnames( imputed_data )  =  phenotype_residuals$eid

    # phenotype_residuals$phenotype  =  scale( phenotype_residuals$phenotype )

    if (is.null( betas )){
        for (snp in 1:nrow( snps )) {
            snps[ snp, 'beta' ]  =  lm( phenotype_residuals$phenotype ~ imputed_data[ snp, ] )%>%
                coef %>%
                '['( 2 ) %>%
                unname
        }
        betas = snps[ , 'beta' ]
    }

    grs  =  t( imputed_data ) %*% as.matrix( betas ) %>%
        scale

    results  =  gxe_interaction( phenotype_residuals$phenotype, grs )

    # aics  =  sapply( results, function( x ) x$aic )

    # smallest  =  which( aics == min( aics ) )
    #
    # texts  =  c(
    #     gxe          =  c( "The data suggests that the genetic and environmental factors are uncorrelated and the outcome is untransformed, however this may be due to insufficient power.\n" ),
    #     pow          =  c( "The data suggests that the genetic and environmental factors are uncorrelated and the outcome is untransformed, however this may be due to insufficient power.\n" ),
    #     gxe_cor      =  c( "The data suggests that that the genetic and environmental factors are correlated, but the outcome is untransformed (or power is insufficient to detect it).\n" ),
    #     gxe_pow      =  c( "The data suggests that the outcome is observed on a transformed scale, but the genetic and environmental factors are not correlated (possibly due to insufficient power).\n" ),
    #     gxe_cor_pow  =  c( "The data supports the hypothesis that the genetic and environmental factors are correlated, and that the outcome is observed on a transformed scale.\n" )
    # )
    #
    # cat( texts[ smallest ] )

    results$betas  =  betas
    results$phenotype  =  phenotype_residuals$phenotype
    results$grs  =  grs

    results
}


ukb_filename  =  '../pheno/ukb21067.csv'
phenotype_name  =  '21001-0.0'
imp_sample_filename  =  '../imp/ukb1638_imp_chr1_v2_s487398.sample'
bgens_path  =  '../imp'
npcs  =  10
covariate_names  =  c( 'eid', '21003-0.0', '31-0.0', '22000-0.0' )
covariate_factor_names  =  c( '31-0.0', '22000-0.0' )
sqc_filename   =  '../geno/ukb_sqc_v2.txt'
geno_filename  =  '../plink/ukb1638_cal_chr1_v2_s488366.fam'

neale_filename  =  '../neale_files/21001.assoc.tsv.gz'
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
snps = final_assoc$rsid

test = ukb_gxe_interaction( phenotype_name,
                            ukb_filename,
                            bgens_path,
                            snps,
                            covariate_names,
                            covariate_factor_names,
                            # sample_ids = NULL,
                            sqc_filename  = sqc_filename,
                            geno_filename = geno_filename,
                            # ids_to_remove = c(),
                            npcs = npcs
                            # covariate_filename = NULL,
                            # betas = NULL
                            )














