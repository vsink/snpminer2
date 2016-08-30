# snpMiner2

The snpMiner2 is a gene annotation script written in Perl that will assist you with annotation of variable call format (vcf) files in bacterial genomes.  
The current version of program is Beta.


PREREQUISITES

In order to install and use this package you will need Perl version
5.8.1 or better.  Some modules within this package depend on other
packages that are distributed separately from Perl.  We recommend that
you have the following packages installed before you install
snpMiner2:


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
  
HOW TO WORK:

The first you should to create database file by gb2db.pl script:

```
./gb2db.pl -o <db_name> -i <genbank_fileformat_file> -snp <snp_collection_file>

```

If you does not have genbank file, you can use -id <organism_id> key to download it automatically.

```
./gb2db.pl -o <db_name> -id <organism_id> -snp <snp_collection_file>

```
snp_collection file is optional and consists SNPs which you want to check in process of annotation. The structure of snp_collection_file is tab-delimeted fields:

```
Name<end_of_line>
<tab><locus><tab><SNP_notation><tab><<description>
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
./snpMiner2.pl -db <db_name> -action annotation -indel -vcf <vcf_file> > output_file
./snpMiner2.pl -db <db_name> -action annotation -vcf all > output_file
./snpMiner2.pl -db <db_name> -action annotation -lvl 2 -vcf all > output_file
./snpMiner2.pl -db <db_name> -action annotation -lvl 3 -vcf all > output_file
./snpMiner2.pl -db <db_name> -action annotation -lvl 4 -vcf all > output_file
./snpMiner2.pl -db <db_name> -action annotation -lvl 5 -vcf all > output_file
./snpMiner2.pl -db <db_name> -action annotation -lvl 6 -vcf all > output_file
./snpMiner2.pl -db <db_name> -action annotation -debug -vcf all > output_file

```
unique SNPs for all vcf files in current directory:

```
./snpMiner2.pl -db <db_name> -action uniq > output_file
./snpMiner2.pl -db <db_name> -action uniq -lvl 1 > output_file
./snpMiner2.pl -db <db_name> -action uniq -lvl 2 > output_file
./snpMiner2.pl -db <db_name> -action uniq -lvl 3 > output_file
./snpMiner2.pl -db <db_name> -action uniq -lvl 4 > output_file
./snpMiner2.pl -db <db_name> -action uniq -lvl 5 > output_file
./snpMiner2.pl -db <db_name> -action uniq -lvl 6 > output_file
```

check SNPs 

```
./snpMiner2.pl -db <db_name> -action check_snp -snp_list <snp_collection_name>

```

information about gene:

```
./snpMiner2.pl -db tb -debug  -action info -locus Rv3919c -snp "276A>C"

```

it will create a report like that:
```
--------------------------------------------------
Locus: Rv3919c (4407528 - 4408202) 
Product: rRNA small subunit methyltransferase G
Note: glucose-inhibited division protein B
Gene length (bp): 675
GC (%): 66.1
..................................................
DNA :
  ATGTCTCCGATCGAGCCCGCGGCGTCTGCGATCTTCGGACCGCGGCTTGGCCTTGCTCGGCGGTACGCCGAAGCGTTGGCGGGACCCGGTGTGGAGCGGGGGCTGGTGGGACCCCGCGAAGTCGGTAGGCTATGGGACCGGCATCTACTGAACTGCGCCGTGATCGGTGAGCTCCTCGAACGCGGTGACCGGGTCGTGGATATCGGTAGCGGAGCCGGGTTGCCGGGCGTGCCATTGGCGATAGCGCGGCCGGACCTCCAGGTAGTTCTCCTAGAACCGCTACTGCGCCGCACCGAGTTTCTTCGAGAGATGGTGACAGATCTGGGCGTGGCCGTTGAGATCGTGCGGGGGCGCGCCGAGGAGTCCTGGGTGCAGGACCAATTGGGCGGCAGCGACGCTGCGGTGTCACGGGCGGTGGCCGCGTTGGACAAGTTGACGAAATGGAGCATGCCGTTGATACGGCCGAACGGGCGAATGCTCGCCATCAAAGGCGAGCGGGCTCACGACGAAGTACGGGAGCACCGGCGTGTGATGATCGCATCGGGCGCGGTTGATGTCAGGGTGGTGACATGTGGCGCGAACTATTTGCGTCCGCCCGCGACCGTGGTGTTCGCACGACGTGGAAAGCAGATCGCCCGAGGGTCGGCACGGATGGCGAGTGGAGGGACGGCGTGA
..................................................
Protein:
  MSPIEPAASAIFGPRLGLARRYAEALAGPGVERGLVGPREVGRLWDRHLLNCAVIGELLERGDRVVDIGSGAGLPGVPLAIARPDLQVVLLEPLLRRTEFLREMVTDLGVAVEIVRGRAEESWVQDQLGGSDAAVSRAVAALDKLTKWSMPLIRPNGRMLAIKGERAHDEVREHRRVMIASGAVDVRVVTCGANYLRPPATVVFARRGKQIARGSARMASGGTAX
Molecular mass (Da) approx. : 24013.8914
Protein length (aa): 224
..................................................
SNP:
  ATGTCTCCGATCGAGCCCGCGGCGTCTGCGATCTTCGGACCGCGGCTTGGCCTTGCTCGGCGGTACGCCGAAGCGTTGGCGGGACCCGGTGTGGAGCGGGGGCTGGTGGGACCCCGCGAAGTCGGTAGGCTATGGGACCGGCATCTACTGAACTGCGCCGTGATCGGTGAGCTCCTCGAACGCGGTGACCGGGTCGTGGATATCGGTAGCGGAGCCGGGTTGCCGGGCGTGCCATTGGCGATAGCGCGGCCGGACCTCCAGGTAGTTCTCCTAGA[A>C]CCGCTACTGCGCCGCACCGAGTTTCTTCGAGAGATGGTGACAGATCTGGGCGTGGCCGTTGAGATCGTGCGGGGGCGCGCCGAGGAGTCCTGGGTGCAGGACCAATTGGGCGGCAGCGACGCTGCGGTGTCACGGGCGGTGGCCGCGTTGGACAAGTTGACGAAATGGAGCATGCCGTTGATACGGCCGAACGGGCGAATGCTCGCCATCAAAGGCGAGCGGGCTCACGACGAAGTACGGGAGCACCGGCGTGTGATGATCGCATCGGGCGCGGTTGATGTCAGGGTGGTGACATGTGGCGCGAACTATTTGCGTCCGCCCGCGACCGTGGTGTTCGCACGACGTGGAAAGCAGATCGCCCGAGGGTCGGCACGGATGGCGAGTGGAGGGACGGCGTGA
```

If you will use snpMiner2 in your work, please, cite us:

```
Viacheslav V. Sinkov. (2016). snpMiner2: vcf files annotation tool. Zenodo. 10.5281/zenodo.51052
```


COPYRIGHT

  © 2014-2016 Viacheslav V. Sinkov. The laboratory of Social and Epidemic Infections, 
"Scientific Centre for Family Health and Human Reproduction Problems" 
16, Timiriazeva st., Irkutsk, Russia.
  
This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

Enjoy!
