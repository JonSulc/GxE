source( 'ukb_gxe.R' )

# Column name for outcome phenotype, here BMI on first visit
phenotype_name  =  '21001-0.0'

# Filename (or vector of filenames) containing any UK Biobank phenotypes
# to use as outcome or covariates to correct for, including age and sex
ukb_filename  =  '../uk_biobank/pheno/ukb21067.csv'

# Path to the imputed files
bgens_path  =  '../uk_biobank/imp'

# GWAS results from the Neale lab
neale_filename = '../uk_biobank/pipeline/neale_files/both_sexes/body/21001_irnt.gwas.imputed_v3.both_sexes.tsv.gz'
variants_filename = '../uk_biobank/pipeline/neale_files/variants.tsv.gz'

snps  =  get_betas_from_neale( neale_filename,
                               variants_filename )

# UK Biobank SQC file and any corresponding *.fam file for IDs
sqc_filename   =  '../uk_biobank/geno/ukb_sqc_v2.txt'
geno_filename  =  '../uk_biobank/plink/ukb1638_cal_chr1_v2_s488366.fam'

gxe  =  ukb_gxe_interaction( phenotype_name,
                             ukb_filename  = ukb_filename,
                             bgens_path    = bgens_path,
                             snps          = snps,
                             sqc_filename  = sqc_filename,
                             geno_filename = geno_filename )

print( gxe )
