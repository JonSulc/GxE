#' Load and prune SNP effect sizes from Neale lab files
#'
#' @name get_betas_from_neale
#'
#' @param neale_filename path to the GWAS summary statistics file
#' @param variants_filename path to the variants information file
#' @param threshold p-value threshold to filter SNPs. Default is 5*10^-8
#' @param prune_distance distance (in base pairs) to determine independence of
#'   SNPs. Default is 5*10^5 (500 kb)
#'
#' @return Returns a \code{data.frame} containing the RSIDs and beta effect
#' estimates of the pruned SNPs, with column names 'rsid' and 'beta'
#' @importFrom data.table fread
#' @importFrom dplyr filter select left_join arrange
#' @importFrom parallel mclapply
#' @export
#'
#' @examples
#'   library( GxE )
#'   download.file( paste0( 'https://www.dropbox.com/s/gulqmkaighh8w3b/',
#'                          '21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz?raw=1' ),
#'                  '21001_irnt.gwas.imputed_v3.both_sexes.tsv.gz',
#'                  method = 'libcurl' )
#'   download.file( paste0( 'https://www.dropbox.com/s/puxks683vb0omeg/',
#'                          'variants.tsv.bgz?raw=1' ),
#'                  'variants.tsv.gz',
#'                  method = 'libcurl' )
#'   snps = get_betas_from_neale( neale_filename = '21001_irnt.gwas.imputed_v3.both_sexes.tsv.gz',
#'                                variants_filename = 'variants.tsv.gz' )
#'


library( data.table )
library( dplyr )
library( parallel )

get_betas_from_neale  =  function( neale_filename,
                                   variants_filename,
                                   threshold      = 5e-8,
                                   prune_distance = 5e5 ) {

    neale_stats  =  fread( cmd = paste( 'zcat', neale_filename ),
                           select = c( 'variant', 'low_confidence_variant', 'beta', 'pval' ) )
    neale_stats  =  filter( neale_stats,
                            # low_confidence_variant == 'false',
                            (!low_confidence_variant | low_confidence_variant == 'false'),
                            pval < threshold )
    neale_stats  =  select( neale_stats, -low_confidence_variant )
    neale_stats  =  left_join( neale_stats,
                               fread( cmd = paste( 'zcat', variants_filename ),
                                      select = c( 'variant', 'rsid', 'chr', 'pos' ) ) )
    neale_stats  =  do.call( rbind,
                             mclapply( 1:22,
                                       function( chromosome ){
                                           data  =  neale_stats[ neale_stats$chr == chromosome, ]
                                           data  =  arrange( data, pval )
                                           final_data  =  data[ 0, ]

                                           while (nrow( data ) != 0) {
                                               final_data  =  rbind( final_data, data[ 1, ] )
                                               new_position  =  data$pos[ 1 ]
                                               data  =  filter( data, abs( pos - new_position ) > prune_distance )
                                           }
                                           final_data
                                       } ) )
    neale_stats[ , c( 'rsid', 'beta' ) ]
}