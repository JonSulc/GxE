library( GxE )
library( PearsonDS )

m    = 100       # number of genetic markers
n    = 1e4       # sample size
a0   = sqrt(.05) # linear effect of GRS on y
b1   = sqrt(.3)  # linear effect of E on y
c1   = sqrt(.05) # interaction effect
d1   = 0         # correlation between E and GRS
skwE = 0         # skewness of E
krtE = 3         # kurtosis of E
skwN = 0         # skewness of the noise
krtN = 3         # kurtosis of the noise
pow  = 1         # transformation power

if (krtE <= skwE^2 + 1 | krtN <= skwN^2 + 1) {
    stop( 'Skew and kurtosis values not compatible (kurtosis > skew^2 + 1 not satisfied)' )
}

sim_num = 1e2 # number of bootstrap / fake GRS

# Simulation

maf  =  matrix( runif( m ), nrow = 1 )
G    =  apply( maf, 2, function(x) rbinom( n, 2, x ) )
eff  =  matrix( runif( m ), ncol = 1 )
eff  =  eff / sqrt( sum( eff^2 ) )
GRS  =  scale( G %*% eff )

E    =  d1 * GRS + sqrt( 1 - d1^2 ) * matrix( rpearson( n, moments = c( 0, 1, skwE, krtE ) ) )
noi  =  matrix( rpearson( n, moments = c( 0, 1, skwN, krtN ) ) )
sig  =  sqrt( 1 - a0^2 - b1^2 - c1^2 * (1 + d1^2) - 2 * a0 * b1 * d1 )

z    =  scale( a0 * GRS + b1 * E + c1 * GRS * E + sig * noi )

Ftrans  =  function( s, p1, p2 ) ( (s-p1)^p2 - 1 ) / p2

if (pow != 0) {
    y  =  scale( Ftrans( z, min(z)-1e-5, pow ) )
} else {
    y  =  scale( log( z - min(z)+1e-5 ) )
}

sel  =  which( abs( y ) > 10 )

while (length( sel ) > 0) {
    y[ sel ] = NA
    y  =  scale( y )
    sel  =  which( abs( y ) > 10 )
}

# Estimate interaction effect for GRS
print( estimate_gxe( y, GRS, sim_num ) )







