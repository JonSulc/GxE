library( parallel )

source( 'regress_NOp.R' )
source( 'IA_fit.R' )
source( 'IA_fit_G2.R' )

.single_gxe  =  function( y,
                          GRS,
                          params ){
    params  =  nlm( function( x ) IA_fit( x, y, GRS ),
                    params )$estimate
    nlm( function( x ) IA_fit_G2( x, y, GRS ),
         c( params[ 1 ], 0, params[ 2:3 ] ) )$estimate
}

estimate_gxe  =  function( y,
                           GRS,
                           sim_num = 100 ){
    if (is.null( dim( y ) )) {
        y = as.matrix( y )
    }

    ao_het  =  regress_NOp( y,
                            cbind( matrix( 1, nrow = length( y ) ),
                                   GRS,
                                   GRS^2 ) )$x
    Xopt  =  mclapply( 1:sim_num,
                       function(x) {
                           ix  =  sample( length(y), length(y), replace = TRUE )
                           .single_gxe( y[   ix, , drop = FALSE ],
                                        GRS[ ix, , drop = FALSE ],
                                        c( ao_het[ 2 ], 0.1, 0 ) )
                       } )

    fGRS  =  simulate_fGRS( y, GRS, sim_num )
    Xopt0  =  mclapply( as.data.frame( fGRS ),
                        function(x) {
                            .single_gxe( y,
                                         as.matrix( x ),
                                         c( ao_het[ 2 ], 0.1, 0 ) )
                        } )
    Xopt    =  do.call( cbind, Xopt  )
    Xopt0   =  do.call( cbind, Xopt0 )
    xopt0   =  apply( Xopt0, 1, mean )
    xopt    =  apply( Xopt,  1, mean )
    SExopt0 =  apply( Xopt0, 1, sd   )
    SExopt  =  apply( Xopt,  1, sd   )
    Pxopt0  =  2 * pnorm( -abs( xopt0 / SExopt0 ) )
    Pxopt   =  2 * pnorm( -abs( xopt  / SExopt  ) )
    tdiff   =  ( xopt - xopt0 ) / sqrt( SExopt^2 + SExopt0^2 )

    list( Xopt    =  Xopt,
          Xopt0   =  Xopt0,
          xopt0   =  xopt0,
          xopt    =  xopt,
          SExopt0 =  SExopt0,
          SExopt  =  SExopt,
          Pxopt0  =  Pxopt0,
          Pxopt   =  Pxopt,
          tdiff   =  tdiff )
}
