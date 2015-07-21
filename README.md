# snpMiner2
The snpMiner2 is a gene annotation script written in Perl that will assist you with annotation of variable call format (vcf) files in bacterial genomes.  
The current version of program is Alpha. This is a worked version without help information. The current version annotate SNPs and does not work with INDELs and intergenic substitutions.

The first you should create database file by gb2db.pl script:

./gb2db.pl -o <db_name> -i <genbank_fileformat_file> -snp <snp_collection_file>

snp_collection file consists SNPs which you want to check in process of annotation. The structure of snp_collection_file is tab-delimeted fields:

SNP_set1
	rv0557	455g>c	test
	rv0557	457c>g	S
	rv0557	532c>g	Cameroon
	rv0557	457c>g	S
	rv0557	221c>t	M.africanum 1a/1b
SNP_set2
	-	489935_G>C	test
SNP_set3
	Rv3833	V105I	test
SNP_set4
	-	761155_Asp749Glu	test
SNP_set5
	Rv0667	codon450	Rv0667_rpoB_531_RIF
intergenic_set1
	-	igpos1673425	inhA_minus15_INH


