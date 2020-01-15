IA_fit_G2  =  function( abc, y, grs ){
    if (is.null( dim( abc ) )) {
        abc  =  array( abc, dim = c( 4, length( abc ) / 4 ) )
    }
    if (is.null( dim( y ) )) {
        y  =  array( y, dim = c( length( y ), 1 ) )
    }
    if (is.null( dim( grs ) )) {
        grs  =  array( grs, dim = c( length( grs ), 1 ) )
    }

    a  =  abc[ 1, ]
    b  =  abc[ 2, ]
    c  =  abc[ 3, ]
    d  =  abc[ 4, ]

    s2  =  1 - a^2 - 2*b^2 - c^2 - d^2
    d1  =  array( 1, dim = c( length( y ), 1 ) )
    d2  =  array( 1, dim = c( 1, length( a ) ) )
    g2  =  grs^2

    v0  =  0.5 * colSums( ( y %*% d2 - grs %*% a - (g2-1) %*% b )^2
                          / ( ( d1 %*% c + grs %*% d )^2 + d1 %*% s2 )
                          + suppressWarnings( log( ( d1 %*% c + grs %*% d )^2 + d1 %*% s2 ) ) )
    v  =  0.5 * colSums( ( y %*% d2 )^2 )

    sel  =  s2 > 0
    v[ sel ]  =  v0[ sel ]
    v
}
