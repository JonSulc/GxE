#' Estimate contribution of GRSxE to variance in outcome phenotype
#'
#' @name estimate_gxe
#'
#' @param phenotypes Outcome phenotype
#' @param grs Polygenic score
#' @param sim_num Number of permutations for bootstrap, fake GRSs, and fake phenotypes
#' @param simulate_phenotype Generate pseudo-phenotype fY
#' @param skewness_range Possible values of skewness to test for in the simulated phenotype
#' @param kurtosis_range possible values of kurtosis to test for in the simulated phenotype
#' @param max_sd Threshold to remove outliers
#' @param use_rslurm Use the rslurm package for parallelization
#' @param rslurm_jobname Name of the rslurm job
#' @param rslurm_suffix Append suffix to rslurm_jobname
#' @param rslurm_overwrite Remove any rslurm files with the same job name in the directory
#' @param slurm_options List of parameters to pass to Slurm
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

estimate_gxe  =  function( phenotypes,
                           grs,
                           sim_num = 100,
                           simulate_phenotype = FALSE,
                           # The rest is only used if simulate_phenotype is TRUE
                           skewness_range = seq( 3, 3, by = 0.2 ),
                           k_range = c( 2:4 ), # kurtosis = skewness^2 + k
                           max_sd = 7,
                           use_rslurm = simulate_phenotype
                             & 'rslurm' %in% rownames(installed.packages()),
                           rslurm_jobname = 'estimate_gxe',
                           slurm_options = list() ) {
    if (is.null( dim( phenotypes ) )) {
        phenotypes = as.matrix( phenotypes )
    }

    ao_het  =  lm( phenotypes ~ grs + I( grs^2 ) )$coefficients[ 2 ]
    individual_coefficients  =  mclapply( 1:sim_num,
                       function(x) {
                           ix  =  sample( length(phenotypes),
                                          length(phenotypes),
                                          replace = TRUE )
                           .single_gxe( phenotypes[   ix, , drop = FALSE ],
                                        grs[ ix, , drop = FALSE ],
                                        c( ao_het, 0.1, 0 ) )
                       } )

    fGRS        =  simulate_fGRS( phenotypes, grs, sim_num )
    individual_coefficients_fgrs  =  mclapply( as.data.frame( fGRS ),
                            function(fgrs) {
                                .single_gxe( phenotypes,
                                             as.matrix( fgrs ),
                                             c( ao_het, 0.1, 0 ) )
                            } )

    coef_names = c( 'alpha1', 'alpha2', 'beta', 'gamma' )

    individual_coefficients      =  do.call( cbind, individual_coefficients  )
    individual_coefficients_fgrs =  do.call( cbind, individual_coefficients_fgrs )
    rownames( individual_coefficients      )  =  coef_names
    rownames( individual_coefficients_fgrs )  =  coef_names
    coefficients         =  apply( individual_coefficients,      1, mean )
    coefficients_fgrs    =  apply( individual_coefficients_fgrs, 1, mean )
    se_coefficients      =  apply( individual_coefficients,      1, sd   )
    se_coefficients_fgrs =  apply( individual_coefficients_fgrs, 1, sd   )
    p_coefficients       =  2 * pnorm( -abs( coefficients      / se_coefficients      ) )
    p_coefficients_fgrs  =  2 * pnorm( -abs( coefficients_fgrs / se_coefficients_fgrs ) )
    t_real_fgrs          =  ( coefficients - coefficients_fgrs ) / sqrt( se_coefficients^2 + se_coefficients_fgrs^2 )
    results  =  list()

    if (simulate_phenotype) {
        cor_y_grs  =  cor( phenotypes, grs )[1]

        find_optimal_fY  =  function( skewness,
                                      kurtosis ){
            fY_full  =  simulate_fY( phenotypes,
                                     grs,
                                     skewness = skewness,
                                     kurtosis = kurtosis,
                                     sim_num = sim_num )
            fphenotype  =  fY_full$fphenotype
            qual  =  sqrt(mean(
                (apply( fphenotype, 2,function( fy, y_sorted ) sort(fy) - y_sorted, sort(phenotypes) ))^2
            ) )

            # thYs = matrix( 0, nrow = 3, ncol = sim_num )

            # for (simulation_n in 1:sim_num) {
            #     minimum  =  optim( c( cor_y_grs, 0.1, 0 ),
            #                        IA_fit,
            #                        gr = NULL,
            #                        y = fphenotype[ , simulation_n ],
            #                        grs = grs )
            #     thYs[    , simulation_n ]  =  minimum$par
            # }
            thYs  =  apply( fphenotype, 2, function( y ){
                optim( c( cor_y_grs, 0.1, 0 ),
                       IA_fit,
                       gr = NULL,
                       y = y,
                       grs = grs )$par
            } )

            list( coefficients = rowMeans( thYs, na.rm = TRUE ),
                  se           = apply( thYs, 1, sd,   na.rm = TRUE ) / sqrt( sim_num ),
                  skewness     = skewness,
                  kurtosis     = kurtosis,
                  rms_diff     = qual,
                  alp          = fY_full$alp,
                  fY_coefficients = fY_full$f_betas )
        }
        parameters  =  expand.grid( skewness = skewness_range,
                                    kurtosis = k_range )
        parameters$kurtosis = parameters$skewness^2 + parameters$kurtosis

        y0   =  scale( phenotypes )
        keep =  abs(y0) < max_sd
        y0   =  y0[ keep, , drop = FALSE ]
        grs0 =  scale( grs[ keep, , drop = FALSE ] )

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
                                             nodes = nrow( parameters ),
                                             cpus_per_node = 1,
                                             add_objects = c( 'y',
                                                              'grs',
                                                              'simulate_fY',
                                                              'cor_y_grs',
                                                              'sim_num',
                                                              'max_sd',
                                                              'IA_fit',
                                                              '.create_zs' ),
                                             libPaths = '/home/josulc/miniconda3/lib/R/library',
                                             slurm_options = slurm_options,
                                             jobname = rslurm_jobname )
            fY_results  =  get_slurm_out( simulate_fY_job, wait = TRUE )

        } else {
            parameters  =  as.data.frame( t( parameters ) )
            fY_results  =  mclapply( parameters, function( param ){
                find_optimal_fY( param[1], param[2] )
            } )
        }
        score  =  sapply( fY_results, function( fY_result ){
            sum(abs(fY_result$coefficients - thY) / thY_SE)
        } )
        results  =  fY_results[[ which.min( score ) ]]
    }

    c( list(
        real_phenotype = list(
            coefficients = setNames( coefficients, coef_names ),
            se = setNames( se_coefficients, coef_names ),
            p  = setNames( p_coefficients,  coef_names ),
            individual_estimates = individual_coefficients
        ),
        fake_grs = list(
            coefficients = setNames( coefficients_fgrs, coef_names ),
            se = setNames( se_coefficients_fgrs, coef_names ),
            p  = setNames( p_coefficients_fgrs,  coef_names ),
            individual_estimates = individual_coefficients_fgrs
        ),
        fake_phenotype = results,
        t_real_fgrs = setNames( t_real_fgrs,       coef_names )
    ) )
}
