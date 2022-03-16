source( 'R/IA_fit.R'  )

.create_zs  =  function( ys,
                         y_sorted,
                         ys_order = order( ys ) ){
  Ys  =  outer( ys,  1:7, function( y, power ) y^power )
  betas  =  lm( y_sorted ~ Ys[ ys_order, ] )$coefficients

  zs  =  sweep( Ys, 2, betas[-1], '*' )
  scale( rowSums( zs ) + betas[1] )
}

simulate_fY  =  function( phenotypes,
                          grs,
                          skewness,
                          kurtosis,
                          sim_num = 100 ){
  if (!requireNamespace( "RNOmni", quietly = TRUE )
      | !requireNamespace( "PearsonDS", quietly = TRUE )) {
    stop("Packages 'RNOmni' and 'PearsonDS' are required for the fake phenotype simulation.",
         call. = FALSE)
  }
  if (skewness^2 > kurtosis-1) {
    stop( "Invalid combination of skewness and kurtosis" )
  }
  phenotypes  =  scale( phenotypes )
  grs         =  scale( grs )

  if (is.null( dim( phenotypes ) )) {
    phenotypes  =  array( phenotypes, dim = c( length( phenotypes ), 1 ) )
  }
  if (is.null( dim( grs ) )) {
    grs  =  array( grs, dim = c( length( grs ), 1 ) )
  }
  z    =  matrix( RNOmni::rankNorm( phenotypes[ , 1 ] ), ncol = 1 )

  a1  =  cor( phenotypes, grs )[1]

  thY =  optim( c( a1, 0.1, 0 ), IA_fit, gr = NULL, y = phenotypes, grs = grs )$par

  y_sorted  =  sort( phenotypes )

  a1s  =  matrix( a1 + seq( -0.25, 0.25, by = 0.01 ), nrow = 1 )
  noi  =  matrix( PearsonDS::rpearson( length( phenotypes ),
                            moments = c( 0, 1, skewness, kurtosis ) ),
                  ncol = 1 )

  yS   =  grs %*% a1s + noi %*% sqrt( 1-a1s^2 )

  zS   =  apply( yS, 2, .create_zs, y_sorted )

  ts  =  apply( zS, 2, function( z, grs ) lm( grs ~ z-1 )$coefficients, grs )

  a1_best  =  a1s[ which.min( abs( ts-a1 ) ) ]

  noi  =  matrix( PearsonDS::rpearson( length( phenotypes ) * sim_num,
                            moments = c( 0, 1, skewness, kurtosis ) ),
                  ncol = sim_num )

  yS  =  a1_best * matrix( grs, nrow = length(grs), ncol = sim_num ) +
         noi * sqrt( 1-a1_best^2 )
  mean_ys  =  apply( apply( yS, 2, sort ), 1, mean )
  mean_Ys  =  outer( mean_ys,  1:7, function( phenotypes, power ) phenotypes^power )
  betas  =  lm( y_sorted ~ mean_Ys )$coefficients
  betas  =  setNames( betas, c( 'Intercept', paste0( 'y', 1:7 ) ) )

  return( list( fphenotypes = apply( yS, 2, .create_zs, y_sorted ),
                alp = a1_best,
                f_betas = betas) )
}








