#' Estimate contribution of GRSxE to variance in outcome phenotype
#'
#' @name estimate_gxe
#'
#' @param y Outcome phenotype
#' @param grs Polygenic score
#' @param sim_num Number of permutations for bootstrap, fake GRSs, and fake phenotypes
#' @param simulate_phenotype Generate pseudo-phenotype fY
#' @param skewness_range Possible values of skewness to test for in the simulated phenotype
#' @param kurtosis_range possible values of kurtosis to test for in the simulated phenotype
#' @param max_sd Threshold to remove outliers
#' @param use_rslurm Use the rslurm package for parallelization
#' @param rslurm_jobname Name of the rslurm job
#' @param rslurm_wait Wait for the rslurm job to complete
#' @param slurm_parameters List of parameters to pass to Slurm
#'
#' @return \code{estimate_gxe} returns a list containing parameter estimates for
#'   alpha1, alpha2, beta, and gamma (\code{xopt}), their standard error
#'   (\code{SExopt}), and the associated p-values (\code{Pxopt}). The
#'   corresponding estimates for the fake GRS are stored in the \code{xopt_fGRS},
#'   \code{SExopt_fGRS}, and \code{Pxopt_fGRS}, respectively.
#'
#'   In addition, the \code{Xopt} and \code{Xopt_fGRS} matrices contain the
#'   estimates for each bootstrap and fake GRS, respectively.
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
#' grs  =  scale( G %*% eff )
#'
#' E    =  d1 * grs + sqrt( 1 - d1^2 ) * matrix( rpearson( n, moments = c( 0, 1, skwE, krtE ) ) )
#' noi  =  matrix( rpearson( n, moments = c( 0, 1, skwN, krtN ) ) )
#' sig  =  sqrt( 1 - a0^2 - b1^2 - c1^2 * (1 + d1^2) - 2 * a0 * b1 * d1 )
#'
#' z    =  scale( a0 * grs + b1 * E + c1 * grs * E + sig * noi )
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
#' estimate_gxe( y, grs, sim_num )

library( parallel )

# source( 'R/regress.R' )
source( 'R/IA_fit.R' )
source( 'R/IA_fit_G2.R' )
source( 'R/simulate_fY.R' )

.single_gxe  =  function( y,
                          grs,
                          params ){
    params  =  optim( params,
                      IA_fit,
                      y = y, grs = grs )$par
    optim( c( params[ 1 ], 0, params[ 2:3 ] ),
           IA_fit_G2,
           y = y, grs = grs )$par
}

estimate_gxe  =  function( y,
                           grs,
                           sim_num = 100,
                           simulate_phenotype = FALSE,
                           skewness_range = seq( -1.6, 1.6, by = 0.2 ),
                           kurtosis_range = c( 2:4 ),
                           max_sd = 7,
                           use_rslurm = TRUE,
                           rslurm_jobname = 'estimate_gxe',
                           rslurm_wait = TRUE,
                           slurm_parameters = list() ){
    if (is.null( dim( y ) )) {
        y = as.matrix( y )
    }

    ao_het  =  lm( y ~ grs + I( grs^2 ) )$coefficients[ 2 ]
    Xopt  =  mclapply( 1:sim_num,
                       function(x) {
                           ix  =  sample( length(y), length(y), replace = TRUE )
                           .single_gxe( y[   ix, , drop = FALSE ],
                                        grs[ ix, , drop = FALSE ],
                                        c( ao_het, 0.1, 0 ) )
                       } )

    fGRS        =  simulate_fGRS( y, grs, sim_num )
    Xopt_fGRS  =  mclapply( as.data.frame( fGRS ),
                            function(fgrs) {
                                .single_gxe( y,
                                             as.matrix( fgrs ),
                                             c( ao_het, 0.1, 0 ) )
                            } )

    parameter_names = c( 'alpha1', 'alpha2', 'beta', 'gamma' )

    Xopt        =  do.call( cbind, Xopt  )
    Xopt_fGRS   =  do.call( cbind, Xopt_fGRS )
    rownames( Xopt      )  =  parameter_names
    rownames( Xopt_fGRS )  =  parameter_names
    xopt        =  apply( Xopt,  1, mean )
    xopt_fGRS   =  apply( Xopt_fGRS, 1, mean )
    SExopt      =  apply( Xopt,  1, sd   )
    SExopt_fGRS =  apply( Xopt_fGRS, 1, sd   )
    Pxopt       =  2 * pnorm( -abs( xopt  / SExopt  ) )
    Pxopt_fGRS  =  2 * pnorm( -abs( xopt_fGRS / SExopt_fGRS ) )
    tdiff       =  ( xopt - xopt_fGRS ) / sqrt( SExopt^2 + SExopt_fGRS^2 )
    results  =  list()

    if (simulate_phenotype) {
        cor_y_grs  =  cor( y, grs )[1]

        find_optimal_fY  =  function( skewness,
                                      kurtosis ){
            fY  =  simulate_fY( y,
                                grs,
                                skewness = skewness,
                                kurtosis = kurtosis,
                                sim_num = sim_num )
            if (is.null( fY )) {
                return( list( thYe     = rep( NA, 3 ),
                              thYe_SE  = rep( NA, 3 ),
                              skewness = skewness,
                              kurtosis = kurtosis,
                              qual     = NA ) )
            }
            qual  =  sqrt(mean(
                (apply( fY, 2,function( fy, y_sorted ) sort(fy) - y_sorted, sort(y) ))^2
            ) )

            thYs = thYs_SE = thZs = thZs_SE  =  matrix( 0, nrow = 3, ncol = sim_num )

            for (simulation_n in 1:sim_num) {
                minimum  =  optim( c( cor_y_grs, 0.1, 0 ),
                                   IA_fit,
                                   gr = NULL,
                                   y = fY[ , simulation_n ],
                                   grs = grs ) # ,
                                   # hessian = TRUE )
                thYs[    , simulation_n ]  =  minimum$par
                # thYs_SE[ , simulation_n ]  =  sqrt( diag( solve( minimum$hessian ) ) )
            }

            list( # thYs     = thYs,
                  # thYs_SE  = thYs_SE,
                  thYe     = apply( thYs, 1, mean, na.rm = TRUE ),
                  thYe_SE  = apply( thYs, 1, sd,   na.rm = TRUE ) / sqrt( sim_num ),
                  skewness = skewness,
                  kurtosis = kurtosis,
                  qual     = qual )
        }
        parameters  =  expand.grid( skewness = skewness_range,
                                    kurtosis = kurtosis_range )
        parameters  =  parameters[ parameters$skewness^2 < parameters$kurtosis - 1, ]

        y0  =  scale( y )
        keep  =  abs(y0) < max_sd
        y0  =  y0[ keep, , drop = FALSE ]
        grs0  =  scale( grs[ keep, , drop = FALSE ] )

        minimum  =  optim( c( cor_y_grs, 0.1, 0 ),
                           IA_fit,
                           gr = NULL,
                           y = y0, grs = grs0,
                           hessian = TRUE )
        thY     =  minimum$par
        thY_SE  =  sqrt( diag( solve( minimum$hessian ) ) )

        if (use_rslurm) {
            if(!require(rslurm))
                stop( "Requires rslurm package." )
            simulate_fY_job  =  slurm_apply( find_optimal_fY,
                                             params = parameters,
                                             add_objects = c( 'y',
                                                              'grs',
                                                              'simulate_fY',
                                                              'cor_y_grs',
                                                              'sim_num',
                                                              'max_sd',
                                                              'IA_fit',
                                                              '.create_zs' ),
                                             jobname = rslurm_jobname )
            # if (rslurm_wait) {
                fY_results  =  get_slurm_out( simulate_fY_job, wait = TRUE )
#             } else {
#                 print( sprintf( "Slurm script submitted. Once completed, run:
# get_fY_results( '%s', nodes = %s ) to see the results.",
#                             rslurm_jobname,
#                             simulate_fY_job$nodes ) )
#             }

        } else {
            parameters  =  as.data.frame( t( parameters ) )
            fY_results  =  mclapply( parameters, function( param ){
                find_optimal_fY( param[1], param[2] )
            } )
            # get_fY_results( fY_results = fY_results )
        }
        score  =  sapply( fY_results, function( fY_result ){
            sum(abs(fY_result$thYe - thY) / thY_SE)
        } )
        results  =  fY_results[[ which.min( score ) ]]
        results$krA  =  results$skewness^2 + results$kurtosis
    }

    c( list( Xopt        =  Xopt,
             Xopt_fGRS   =  Xopt_fGRS,
             xopt        =  setNames( xopt,        parameter_names ),
             xopt_fGRS   =  setNames( xopt_fGRS,   parameter_names ),
             SExopt      =  setNames( SExopt,      parameter_names ),
             SExopt_fGRS =  setNames( SExopt_fGRS, parameter_names ),
             Pxopt       =  setNames( Pxopt,       parameter_names ),
             Pxopt_fGRS  =  setNames( Pxopt_fGRS,  parameter_names ),
             tdiff       =  setNames( tdiff,       parameter_names ) ),
       results )
}
