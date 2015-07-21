# snpMiner2
The snpMiner2 is a gene annotation script written in Perl that will assist you with annotation of variable call format (vcf) files in bacterial genomes.  
The current version of program is Alpha. This is a worked version without help information. The current version annotate SNPs and does not work with INDELs and intergenic substitutions.

The snpMiner2.pl depends on Bioperl and  following modules:

```
Bio::SeqIO;
Bio::Perl;
Bio::DB::GenBank;
Pod::Usage;
Getopt::Long;
File::Temp qw/ tempfile tempdir /;
Storable qw(freeze thaw);
Data::Dumper;
Carp;
IO::Compress::RawDeflate qw(rawdeflate $RawDeflateError);
IO::Uncompress::RawInflate qw(rawinflate $RawInflateError);
List::Compare;
Cwd;
List::Util qw(first max maxstr min minstr reduce shuffle sum);
```


The first you should to create database file by gb2db.pl script:

```
./gb2db.pl -o <db_name> -i <genbank_fileformat_file> -snp <snp_collection_file>

```

If you does not have genbank file, you can use -id <organism_id> key to download it automatically.

```
./gb2db.pl -o <db_name> -id <organism_id> -snp <snp_collection_file>

```
snp_collection file consists SNPs which you want to check in process of annotation. The structure of snp_collection_file is tab-delimeted fields:

```
Name<end_of_line>
<tab><locus><tab><SNP_notation><tab><type><tab><description>
```

The SNP_notation field supports following types:

```
1. 455g>c (SNP position in gene. Needs a locus name.)
2. 489935_G>C (SNP position in genome)
3. V105I (AA position in gene. Needs a locus name.)
4. 761155_Asp749Glu (AA position in genome and gene)
5. codon450 (codon position in gene. Needs a locus name.)
6. igpos1673425 (position in genome for intergenic SNPs)
```
Example of snp_collection_file's fields:
```
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
```


Usage snpMiner2:

The numbers in action name are different levels of output information.

annotation of vcf file

```
./snpMiner2.pl -db <db_name> -action annotation -vcf <vcf_file> > output_file
./snpMiner2.pl -db <db_name> -action annotation.1 -vcf <vcf_file> > output_file
./snpMiner2.pl -db <db_name> -action annotation.11 -vcf <vcf_file> > output_file
./snpMiner2.pl -db <db_name> -action annotation.2 -vcf <vcf_file> > output_file
./snpMiner2.pl -db <db_name> -action annotation.3 -vcf <vcf_file> > output_file
```
find uniq SNP for all vcf files in currend directory:

```
./snpMiner2.pl -db <db_name> -action uniq > output_file
./snpMiner2.pl -db <db_name> -action uniq.1 > output_file
./snpMiner2.pl -db <db_name> -action uniq.2 > output_file
./snpMiner2.pl -db <db_name> -action uniq.3 > output_file
./snpMiner2.pl -db <db_name> -action uniq.4 > output_file
./snpMiner2.pl -db <db_name> -action uniq.5 > output_file
./snpMiner2.pl -db <db_name> -action uniq.6 > output_file
./snpMiner2.pl -db <db_name> -action uniq.7 > output_file

```

check SNP 

```
./snpMiner2.pl -db <db_name> -action check_snp -snp_list <snp_collection_name>

```
