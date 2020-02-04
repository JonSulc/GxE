library( GxE )

download.file( 'https://www.dropbox.com/s/gulqmkaighh8w3b/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz?raw=1',
               '21001_irnt.gwas.imputed_v3.both_sexes.tsv.gz',
               method = 'libcurl' )
download.file( 'https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?raw=1',
               'variants.tsv.gz',
               method = 'libcurl' )

snps = get_betas_from_neale( neale_filename = '23100_irnt.gwas.imputed_v3.both_sexes.tsv.gz',
                             variants_filename = 'variants.tsv.gz' )

gxe  =  ukb_estimate_gxe( phenotype_name = '23100-0.0',
                          ukb_filename   = 'uk_biobank/pheno/ukb21067.csv',
                          bgens_path     = 'uk_biobank/imp',
                          snps           = snps,
                          sqc_filename   = 'uk_biobank/geno/ukb_sqc_v2.txt',
                          fam_filename   = 'uk_biobank/plink/ukb1638_cal_chr1_v2_s488366.fam' )

print( gxe )
