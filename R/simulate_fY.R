# library( R.matlab )
library( RNOmni ) # IRNT
library( PearsonDS ) # Random number generator for Pearson distribution
# library( numDeriv ) # Hessian

source( 'R/IA_fit.R'  )

# data_dir  =  '/data/sgg2/zoltan/project/Project_tailPDF'
# phenotype  =  '21001'
#
# tmp  =  readMat( sprintf( '%s/data-mat/ukbb/data_ukbb%s_v2.mat',
#                           data_dir,
#                           phenotype ) )
# y   =  tmp$y0
# grs  =  tmp$GRS
# rm( tmp )
#
# sel  =  abs(y0) < 7
# y0   =  scale( y0[ sel ] )
# grs  =  scale( grs[ sel ] )
# z0   =  matrix( rankNorm( y0[ , 1 ] ), ncol = 1 )

.create_zs  =  function( ys,
                         y_sorted,
                         skewness,
                         kurtosis ){
  y_order  =  order( ys )
  Ys  =  outer( ys,  1:7, function( y, power ) y^power )
  betas  =  lm( y_sorted ~ Ys[ y_order, ] )$coefficients

  zs  =  sweep( Ys, 2, betas[-1], '*' )
  scale( rowSums( zs ) + betas[1] )
}

simulate_fY  =  function( y,
                          grs,
                          skewness,
                          kurtosis,
                          sim_num = 100 ){
  y    =  scale( y )
  grs  =  scale( grs )

  if (is.null( dim( y ) )) {
    y  =  array( y, dim = c( length( y ), 1 ) )
  }
  if (is.null( dim( grs ) )) {
    grs  =  array( grs, dim = c( length( grs ), 1 ) )
  }
  z    =  matrix( rankNorm( y[ , 1 ] ), ncol = 1 )

  ####### Check for a1 = corrz( y0, GRS )
  # a1   =  summary( lm( y0~grs ) )$coef[ 'grs', 'Estimate' ]
  a1  =  cor( y, grs )[1]

  thY =  optim( c( a1, 0.1, 0 ), IA_fit, gr = NULL, y = y, grs = grs )$par
  # thZ =  optim( thY,             IA_fit, gr = NULL, y = z, grs = grs )$par

  y_sorted  =  sort( y )

  a1s  =  matrix( a1 + seq( -0.25, 0.25, by = 0.01 ), nrow = 1 )
  try({
    noi  =  matrix( rpearson( length( y ),
                              moments = c( 0, 1, skewness, kurtosis ) ),
                    ncol = 1 )
  }, silent = TRUE )
  if (!exists( 'noi' )) {
    return()
  }

  yS   =  grs %*% a1s + noi %*% sqrt( 1-a1s^2 )

  zS   =  apply( yS, 2, .create_zs, y_sorted, skewness, kurtosis )

  ts  =  apply( zS, 2, function( z, grs ) lm( grs ~ z-1 )$coefficients, grs )

  a1_best  =  a1s[ which.min( abs( ts-a1 ) ) ]

  noi  =  matrix( rpearson( length( y ) * sim_num,
                            moments = c( 0, 1, skewness, kurtosis ) ),
                  ncol = sim_num )

  yS  =  a1_best * grs %*% matrix( 1, ncol = sim_num ) + noi * sqrt( 1-a1_best^2 )
  # zS  =  apply( yS, 2, .create_zs, y_sorted, skewness, kurtosis )

  return( apply( yS, 2, .create_zs, y_sorted, skewness, kurtosis ) )

  # qual  =  sqrt(mean( (apply( zS, 2, function( z, y_sorted ) sort(z) - y_sorted, y_sorted ))^2 ))
  #
  # thYs = thYs_SE = thZs = thZs_SE  =  matrix( 0, nrow = 3, ncol = sim_num )
  #
  # for (simulation_n in 1:sim_num) {
  #   minimum  =  optim( c( a1, 0.1, 0 ),
  #                      IA_fit,
  #                      gr = NULL,
  #                      # list( y = zS[ , simulation_n ], g = grs ),
  #                      y = zS[ , simulation_n ],
  #                      grs = grs,
  #                      hessian = TRUE )
  #   thYs[ , simulation_n ]  =  minimum$par
  #   thYs_SE[ , simulation_n ]  =  sqrt( diag( solve( minimum$hessian ) ) )
  #
  #   # tt  =  rankNorm( zS[ , simulation_n ] )
  #   # minimum  =  optim( c( a1, 0.1, 0 ),
  #   #                    IA_fit,
  #   #                    gr = NULL,
  #   #                    y = tt,
  #   #                    grs = grs,
  #   #                    hessian = TRUE )
  #   # thZs[ , simulation_n ]  =  minimum$par
  #   # thZs_SE[ , simulation_n ]  =  sqrt( diag( solve( minimum$hessian ) ) )
  # }
  #
  # list( thYs    = thYs,
  #       thYs_SE = thYs_SE,
  #       # thZs    = thZs,
  #       # thZs_SE = thZs_SE,
  #       qual    = qual )
}








