simulate_fGRS = function( y0, GRS, sim_num ){
    b1   =  mean( y0 * GRS )
    b2   =  mean( y0 * (GRS^2) )
    mu2  =  mean( y0^2 )
    mu3  =  mean( y0^3 )
    mu4  =  mean( y0^4 )
    mu5  =  mean( y0^5 )
    noi  =  t( scale( t( matrix( rnorm( length(y0) * sim_num ), ncol = sim_num ) ) ) )

    A  =  mu3^3 - 2 * mu4 * mu3 + mu5
    B  =  b1 - mu4 * b1 + mu3^2 * b1
    D  =  b1^2 - 2 * mu4 * b1^2 +
          2 * mu3^2 * b1^2 + mu4^2 * b1^2 +
          mu3^3 * b2 -
          2 * mu4 * mu3 * b2 -
          mu5 * mu3 * b1^2 + mu5 * b2

    if (D > 0) {
        a2  =  ( B + sqrt( D ) ) / A
        a0  =  -a2 * mu2
        a1  =  ( b1 - a2 * mu3 )
        varG0  =  a0^2 + a1^2 + a2^2 * mu4 + 2 * a0 * a2 + 2 * a1 * a2 * mu3

        G0  =  a0 + a1 * y0 + a2 * y0^2
        fGRS  =  G0 %*% matrix( 1, ncol = sim_num ) +
                 sqrt( 1-varG0 ) * noi
    } else {
        fGRS  =  b1 * y0 * matrix( 1, ncol = sim_num ) + sqrt( 1-b1^2 ) * noi
    }

}




