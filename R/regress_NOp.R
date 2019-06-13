regress_NOp  =  function( y, X, d_f = NULL ){
    # Model: y = X*beta
    if (is.null( dim(y) )) {
        y = matrix( y )
    }

    nonsel  =  which( apply( cbind( y, X ), 1, function(x) all(!is.na(x)) ) )
    y  =  y[ nonsel, , drop = FALSE ]
    if (diff(range(y)) == 0) {
        return()
    }
    X  =  X[ nonsel, ]

    NA_field    =  which( apply( apply( X, 2, range ), 2, diff ) == 0 )
    noNA_field  =  which( apply( apply( X, 2, range ), 2, diff ) != 0 )

    if (length(NA_field) != 0) {
        noNA_field  =  c( NA_field[ 1 ], noNA_field )
    } else {
        y  =  y - mean( y )
    }

    x  =  matrix( NA, nrow = ncol(X), ncol = 1 )
    se =  x
    X  =  X[ , noNA_field ]

    Z  =  solve( t(X) %*% X )
    tmp_x  =  Z %*% ( t(X) %*% y )
    x[ noNA_field ]  =  tmp_x

    if (is.null( d_f )) {
        d_f  =  length( tmp_x )
    }

    err  =  y - X %*% tmp_x
    sig  =  sqrt( sum(err^2) / (length( y ) - d_f) )
    tmp_se  =  sig * sqrt( diag( Z ) )
    se[ noNA_field ]  =  tmp_se
    Cov  =  sig^2 * Z

    sel  =  which( Im( se ) != 0 )
    x[  sel ]  =  NA
    se[ sel ]  =  NA

    list( x = x, se = se, Cov = Cov )
}
