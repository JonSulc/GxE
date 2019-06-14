#' Estimate contribution of GRSxE to variance in outcome phenotype
#'
#' @name estimate_gxe
#'
#' @param y Outcome phenotype
#' @param GRS Polygenic score
#' @param sim_num Number of permutations for bootstrap and fake GRSs
#'
#' @return \code{estimate_gxe} returns a list containing parameter estimates for
#'   alpha1, alpha2, beta, and gamma (\code{xopt}), their standard error
#'   (\code{SExopt}), and the associated p-values (\code{Pxopt}). The
#'   corresponding estimates for the fake GRS are stored in the \code{xopt0},
#'   \code{SExopt0}, and \code{Pxopt0}, respectively.
#'
#'   In addition, the \code{Xopt} and \code{Xopt0} matrices contain the
#'   estimates for each bootstrap and fake GRS.
#'
#'   The \code{tdiff} contains the t-statistic for the difference between the
#'   real data estimates and those from the fake GRS.
#' @importFrom dplyr bind_cols
#' @importFrom parallel mclapply
#' @importFrom stats lm na.exclude nlm pnorm sd setNames
#' @export
#'
#' @examples
#' library( GxE )
#' library( PearsonDS )
#'
#' m    = 100       # number of genetic markers
#' n    = 1e4       # sample size
#' a0   = sqrt(.05) # linear effect of GRS on y
#' b1   = sqrt(.3)  # linear effect of E on y
#' c1   = sqrt(.05) # interaction effect
#' d1   = 0         # correlation between E and GRS
#' skwE = 0         # skewness of E
#' krtE = 3         # kurtosis of E
#' skwN = 0         # skewness of the noise
#' krtN = 3         # kurtosis of the noise
#' pow  = 1         # transformation power
#'
#' if (krtE <= skwE^2 + 1 | krtN <= skwN^2 + 1) {
#'     stop( 'Skew and kurtosis values not compatible (kurtosis > skew^2 + 1 not satisfied)' )
#' }
#'
#' sim_num = 1e2 # number of bootstrap / fake GRS
#'
#' # Simulation
#'
#' maf  =  matrix( runif( m ), nrow = 1 )
#' G    =  apply( maf, 2, function(x) rbinom( n, 2, x ) )
#' eff  =  matrix( runif( m ), ncol = 1 )
#' eff  =  eff / sqrt( sum( eff^2 ) )
#' GRS  =  scale( G %*% eff )
#'
#' E    =  d1 * GRS + sqrt( 1 - d1^2 ) * matrix( rpearson( n, moments = c( 0, 1, skwE, krtE ) ) )
#' noi  =  matrix( rpearson( n, moments = c( 0, 1, skwN, krtN ) ) )
#' sig  =  sqrt( 1 - a0^2 - b1^2 - c1^2 * (1 + d1^2) - 2 * a0 * b1 * d1 )
#'
#' z    =  scale( a0 * GRS + b1 * E + c1 * GRS * E + sig * noi )
#'
#' Ftrans  =  function( s, p1, p2 ) ( (s-p1)^p2 - 1 ) / p2
#'
#' if (pow != 0) {
#'     y  =  scale( Ftrans( z, min(z)-1e-5, pow ) )
#' } else {
#'     y  =  scale( log( z - min(z)+1e-5 ) )
#' }
#'
#' sel  =  which( abs( y ) > 10 )
#'
#' while (length( sel ) > 0) {
#'     y[ sel ] = NA
#'     y  =  scale( y )
#'     sel  =  which( abs( y ) > 10 )
#' }
#'
#' # Estimate interaction effect for GRS
#' estimate_gxe( y, GRS, sim_num )

library( parallel )

source( 'R/regress_NOp.R' )
source( 'R/IA_fit.R' )
source( 'R/IA_fit_G2.R' )

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
    parameter_names = c( 'alpha1', 'alpha2', 'beta', 'gamma' )

    Xopt    =  do.call( cbind, Xopt  )
    rownames( Xopt  )  =  parameter_names
    Xopt0   =  do.call( cbind, Xopt0 )
    rownames( Xopt0 )  =  parameter_names
    xopt0   =  apply( Xopt0, 1, mean )
    xopt    =  apply( Xopt,  1, mean )
    SExopt0 =  apply( Xopt0, 1, sd   )
    SExopt  =  apply( Xopt,  1, sd   )
    Pxopt0  =  2 * pnorm( -abs( xopt0 / SExopt0 ) )
    Pxopt   =  2 * pnorm( -abs( xopt  / SExopt  ) )
    tdiff   =  ( xopt - xopt0 ) / sqrt( SExopt^2 + SExopt0^2 )

    list( Xopt    =  Xopt,
          Xopt0   =  Xopt0,
          xopt0   =  setNames( xopt0,   parameter_names ),
          xopt    =  setNames( xopt,    parameter_names ),
          SExopt0 =  setNames( SExopt0, parameter_names ),
          SExopt  =  setNames( SExopt,  parameter_names ),
          Pxopt0  =  setNames( Pxopt0,  parameter_names ),
          Pxopt   =  setNames( Pxopt,   parameter_names ),
          tdiff   =  setNames( tdiff,   parameter_names ) )
}
