#!/usr/bin/env perl
eval {
    require Pod::Usage;
    require Getopt::Long;
    require File::Temp;
    require File::Basename;
    require Storable;
    require Data::Dumper;
    require Carp;
    require IO::Compress::RawDeflate;
    require IO::Uncompress::RawInflate;
    require List::Compare;
    require Cwd;
    require List::Util;

};
if ($@) {
    die("You have to install all required modules!\n");
}
use Pod::Usage;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;
use File::Basename;
use Storable qw(freeze thaw);
use Data::Dumper;
use Carp;
use IO::Compress::RawDeflate qw(rawdeflate $RawDeflateError);
use IO::Uncompress::RawInflate qw(rawinflate $RawInflateError);
use List::Compare;
use Cwd;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;
use warnings;
use Time::HiRes qw(usleep ualarm gettimeofday tv_interval);
no warnings qw/uninitialized/;

use constant {
    TANG_TH      => 0.4, #treshold for Tang index
    GD_TH        => 100,
    BLOSSUM45_TH => 0,
    PAM30_TH     => 0
};

# use Benchmark qw(:all) ;
# my $t0 = Benchmark->new;

my $GenomeSequenceBioSeq      = "";        # BioSeq объект с геномом
# my $LocusInfoReader_res_codon_seq   = "";
my $tree_codon_sequence = "";
my $db_version          = "";
my $SNPListExists     = 0;
my $RefGenomeName     = "";
my $OrganismName       = "";
my $version             = "0.3.1";
my $about_msg
    = "Program: snpMiner2 (Tool for analysis of variable call format files)\n";

#-------------------------------
my %database;
my %genetic_code;
my %tang_dic;
my %deltaH_dic;
my %GD_dic;
my %B45_dic;
my %PAM30_dic;
my %aa_hash;
my %snp_list_h;
my %snp_igpos_list_h;
my %aa_mw;
my %intergenic;
my @to_circos_a;
my %options;    # options hash
my %snp_per_gene_h;

GetOptions(
    \%options,    'db=s',    'action=s', 'snp_list=s',
    'vcf=s',      "about!",  'help',     'snp_th=s',
    'list=s',     "DP=s",    "target=s", "tang_th=s",
    "locus=s",    "prosite", "o=s",      "snp=s",
    "seq_type=s", "indel",   "lvl=s",    "debug",
    "intergenic", "toR",     "codon"
);

&show_help( "verbose" => 99, -utf8 ) if $options{help};

if ( $options{action} eq "rf" ) {    #extract db as dump file
    &load_db( $options{db} );
    &write_dump;
}

elsif ( $options{action} eq "check_snp" )
{                                    #check existance of snp from database

    if ( $options{snp_list} ne "" && $options{db} ne "" ) {
        &load_db( $options{db} );

        if ( exists $snp_list_h{ $options{snp_list} } ) {

            &use_snp_classification_coding;
        }
        if ($options{snp_list}  eq "?") {

          while (my ($key, $value) = each %snp_list_h) {
            print $key, "\n";
        }
      }

        elsif ( exists $snp_igpos_list_h{ $options{snp_list} } ) {

            # print "Yes, this is intergenic positions!!!\n";
            &use_snp_classification_intergenic( $options{snp_list} );
          }
        else {
            print "The SNP classification name is not found!\n";
        }
    }
    else {
        print "Use -snp_list key to set SNP classification name!\n";
    }

}

elsif ( $options{action} eq "check_list" ) {

    &load_db( $options{db} );

    &compare_list_by_snp( $options{vcf} );

}
elsif ( $options{action} eq "uniq" ) {
    $options{lvl} = 1 if $options{lvl} eq "";
    &load_db( $options{db} );

    &open_multiple_vcfs;

}
elsif ( $options{action} eq "t4" ) {

    &load_db( $options{db} );

    &test4( $options{vcf} );

}
<<<<<<< HEAD
elsif ( $options{action} eq "t5" || $options{action} eq "t5.1" || $options{action} eq "make_seq_aa" ) {

    &load_db( $options{db} );

    &make_seq_aa( $options{vcf} );

}

elsif ( $options{action} eq "t6" || $options{action} eq "make_seq_nc") {

    &load_db( $options{db} );

    &make_seq_nc( $options{vcf} );
=======
elsif ($options{action} eq "t5"
    || $options{action} eq "t5.1"
    || $options{action} eq "make_seq_aa" )
{

    &load_db( $options{db} );

    &make_seq_aa( $options{vcf} );
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e

}

elsif ( $options{action} eq "t6" || $options{action} eq "make_seq_nc" ) {

    &load_db( $options{db} );

    &make_seq_nc( $options{vcf} );

}

elsif ( $options{action} eq "annotation"
    and $options{vcf} ne ""
    and $options{vcf} ne "all" )
{
    $options{lvl} = 1 if $options{lvl} eq "";
    if (   $options{lvl} == 1
        || $options{lvl} == 2
        || $options{lvl} == 3
        || $options{lvl} == 4
        || $options{lvl} == 5
        || $options{lvl} == 6 )
    {
        &load_db( $options{db} );
        &open_single_vcf( $options{vcf}, 1 );
    }

    #
}
elsif ( $options{action} eq "annotation" and $options{vcf} eq "all" ) {
    $options{lvl} = 1 if $options{lvl} eq "";
    if (   $options{lvl} == 1
        || $options{lvl} == 2
        || $options{lvl} == 3
        || $options{lvl} == 4
        || $options{lvl} == 5
        || $options{lvl} == 6 )
    {
        &load_db( $options{db} );
        my @files = glob('*.vcf');
        foreach my $file (@files) {
            print "-" x 60 . "\n$file\n" . "-" x 60 . "\n" if $options{debug};
            open_single_vcf( $file, 1 );
        }
    }

    #
}

elsif ( $options{action} eq "make_seq" and $options{seq_type} eq "axt" ) {

    &load_db( $options{db} );

    &make_axt_sequence( $options{vcf} );

}

elsif ( $options{action} eq "info" ) {
    if ( $options{locus} ne "" ) {
        &load_db( $options{db} );
        &get_information( $options{locus} );

    }
    else {
        print "Use -locus <LOCUS NAME> option. \n";
    }

}

if ( $options{about} ) {
    print "-" x 60
        . "\nThe snpMiner2 ($version) was written by Viacheslav V.Sinkov, 2015\n\t\thttps://github.com/vsink/snpminer2\n"
        . "-" x 60 . "\n";
}

# if ( $options{examples} ) {
#     print "-" x 60
#         . "\n Annotation of the vcf file:\n./snpMiner2.pl -db db_name -action annotation -vcf vcf_file_name.vcf\n"
#         . "-" x 60 . "\n";
# }

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub write_dump {

    # write loaded db as hash-dump ({db_name}.dump)
    my $buffer = "";
    rawinflate $options{db} => \$buffer
        or die "rawinflate failed: $RawInflateError\n";
    my $hash_ref = thaw($buffer);
    %database = %{$hash_ref};

    open my $FH, '>', $options{db} . '.dump';
    print $FH Dumper( \%database );
    close $FH;

    # print

}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub codon2aa {

    # translate codon to amino-acid (three letter)
    my $codon = shift;

    $codon = uc $codon;

    if ( exists $genetic_code{$codon} && $genetic_code{$codon} ne "STOP" ) {
        return $genetic_code{$codon};
    }
    if ( exists $genetic_code{$codon} && $genetic_code{$codon} eq "STOP" ) {
        return "X";
    }
    else {

        return "BAD_CODON";

    }
}

sub aa_decode {

    #convert three letter AA to short nomenclature (one letter)
    my $aa_input = shift;

    if ( exists $aa_hash{$aa_input} ) {
        return $aa_hash{$aa_input};
    }
    elsif ( $aa_input eq "X" ) {
        return "X";
    }

    # else
    # {
    #     return "X";
    # }

}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub aa_mw_calc {

    #convert three letter AA to short nomenclature (one letter)
    my $aa_input = shift;

    if ( exists $aa_mw{$aa_input} ) {
        return $aa_mw{$aa_input};
    }

}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub load_db {

    #load genome database created by gb2db.pl
    my $buffer      = "";
    my $loc_in_file = shift;
    rawinflate $loc_in_file => \$buffer
        or die "rawinflate failed: $RawInflateError\n";
    my $hash_ref = thaw($buffer);
    %database     = %{ $hash_ref->{annotations} };    # CDS and gene sequences
    %intergenic   = %{ $hash_ref->{intergenic} };     # CDS and gene sequences
    %tang_dic     = %{ $hash_ref->{dict}{tang} };     # Tang index values
    %deltaH_dic   = %{ $hash_ref->{dict}{deltaH} };   # Delta H values
    %GD_dic       = %{ $hash_ref->{dict}{GD} };       # Grantham values
    %B45_dic      = %{ $hash_ref->{dict}{B45} };      # BLOSUM45 values
    %PAM30_dic    = %{ $hash_ref->{dict}{PAM30} };    # PAM30 values
    %genetic_code = %{ $hash_ref->{dict}{codon} };    #codon table
    %aa_hash      = %{ $hash_ref->{dict}{aa} };       #amino acid table
    %aa_mw        = %{ $hash_ref->{dict}{aa_mw} };    #amino acid table

    #------------------------------------
    $db_version      = $database{options}{version};
    $RefGenomeName = $database{options}{reference};
    $OrganismName   = $database{options}{organism_name};
    $SNPListExists
        = scalar $database{options}{snp_list_exists}; #snp_list existance flag
    my $igpos_list_exists = scalar $database{options}{igpos_list_exists}
        ;                                             #snp_list existance flag
    if ( $SNPListExists == 1 ) {
        %snp_list_h = %{ $hash_ref->{dict}{snp_list}{coding} };    # snp_list
        if ( $igpos_list_exists == 1 ) {
            %snp_igpos_list_h = %{ $hash_ref->{dict}{snp_list}{intergenic} }
                ;    # intergenic positions list
        }
    }

}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

<<<<<<< HEAD
# sub annotate_vcf_new {
#     my $t0 = [gettimeofday];
#     my %vcf_data = &open_single_vcf( $options{vcf}, 1 );
#     my %snp_per_gene_h;

#     foreach my $key ( sort { $a <=> $b } keys %vcf_data ) {
#         next if $vcf_data{$key}{locus} eq "";
#         print
#             "$vcf_data{$key}{locus}\t$key\t$vcf_data{$key}{snp}\t$vcf_data{$key}{ref_codon}/$vcf_data{$key}{alt_codon}\t$vcf_data{$key}{ref_aa_short}$vcf_data{$key}{aa_pos}$vcf_data{$key}{alt_aa_short}\t$vcf_data{$key}{snp_type}\t$vcf_data{$key}{product}\n"
#             if $options{lvl} == 1 and $vcf_data{$key}{_type} eq "snp";
#         print
#             "$vcf_data{$key}{locus}\t$key\t$vcf_data{$key}{snp}\t$vcf_data{$key}{ref_codon}/$vcf_data{$key}{alt_codon}\t$vcf_data{$key}{ref_aa_short}$vcf_data{$key}{aa_pos}$vcf_data{$key}{alt_aa_short}\t$vcf_data{$key}{snp_type}\t$vcf_data{$key}{tang}\t$vcf_data{$key}{gd}\t$vcf_data{$key}{blossum45}\t$vcf_data{$key}{pam30}\t$vcf_data{$key}{AACI}\t$vcf_data{$key}{product}\n"
#             if $options{lvl} == 2 and $vcf_data{$key}{_type} eq "snp";
#         print
#             "$vcf_data{$key}{locus}\t$key\t$vcf_data{$key}{snp}\t$vcf_data{$key}{ref_codon}/$vcf_data{$key}{alt_codon}\t$vcf_data{$key}{ref_aa_short}$vcf_data{$key}{aa_pos}$vcf_data{$key}{alt_aa_short}\t$vcf_data{$key}{snp_type}\t$vcf_data{$key}{tang}\t$vcf_data{$key}{gd}\t$vcf_data{$key}{blossum45}\t$vcf_data{$key}{pam30}\t$vcf_data{$key}{AACI}\t$vcf_data{$key}{product}\n"
#             if $options{lvl} == 3
#             and $vcf_data{$key}{AACI} ne "----"
#             and $vcf_data{$key}{_type} eq "snp";
#         print
#             "$vcf_data{$key}{locus}\t$key\t$vcf_data{$key}{snp}\t$vcf_data{$key}{ref_codon}/$vcf_data{$key}{alt_codon}\t$vcf_data{$key}{ref_aa_short}$vcf_data{$key}{aa_pos}$vcf_data{$key}{alt_aa_short}\t$vcf_data{$key}{snp_type}\t$vcf_data{$key}{tang}\t$vcf_data{$key}{gd}\t$vcf_data{$key}{blossum45}\t$vcf_data{$key}{pam30}\t$vcf_data{$key}{AACI}\t$vcf_data{$key}{product}\n"
#             if $options{lvl} == 4
#             and $vcf_data{$key}{AACI} eq "++++"
#             and $vcf_data{$key}{_type} eq "snp";

#         #---------------------------------------------------------
#         if ( $options{indel} ) {
#             print
#                 "$vcf_data{$key}{locus}\t$vcf_data{$key}{'snp_genome_pos'}\t$vcf_data{$key}{'snp_gene_pos'}$vcf_data{$key}{'short_indel_note'}$vcf_data{$key}{'indel'}\t$vcf_data{$key}{'full_formated_indel'}\t$vcf_data{$key}{'frame_shift'}\t$vcf_data{$key}{'snp_type'}\t$vcf_data{$key}{'product'}\n"
#                 if $options{lvl} == 1 and $vcf_data{$key}{_type} eq "indel";
#             print
#                 "$vcf_data{$key}{locus}\t$vcf_data{$key}{'snp_genome_pos'}\t$vcf_data{$key}{'snp_gene_pos'}$vcf_data{$key}{'short_indel_note'}$vcf_data{$key}{'indel'}\t$vcf_data{$key}{'full_formated_indel'}\t$vcf_data{$key}{'frame_shift'}\t$vcf_data{$key}{'snp_type'}\t$vcf_data{$key}{tang}\t$vcf_data{$key}{gd}\t$vcf_data{$key}{blossum45}\t$vcf_data{$key}{pam30}\t$vcf_data{$key}{AACI}\t$vcf_data{$key}{product}\n"
#                 if $options{lvl} == 2 and $vcf_data{$key}{_type} eq "indel";
#             print
#                 "$vcf_data{$key}{locus}\t$vcf_data{$key}{'snp_genome_pos'}\t$vcf_data{$key}{'snp_gene_pos'}$vcf_data{$key}{'short_indel_note'}$vcf_data{$key}{'indel'}\t$vcf_data{$key}{'full_formated_indel'}\t$vcf_data{$key}{'frame_shift'}\t$vcf_data{$key}{'snp_type'}\t$vcf_data{$key}{tang}\t$vcf_data{$key}{gd}\t$vcf_data{$key}{blossum45}\t$vcf_data{$key}{pam30}\t$vcf_data{$key}{AACI}\t$vcf_data{$key}{product}\n"
#                 if $options{lvl} == 3
#                 and $vcf_data{$key}{AACI} ne "----"
#                 and $vcf_data{$key}{_type} eq "indel";
#         }

#         # }

#         #---------------------------------------------------------
#         $snp_per_gene_h{ $vcf_data{$key}{locus} }{$key}
#             = "$vcf_data{$key}{snp}\t$vcf_data{$key}{'ref_codon'}/$vcf_data{$key}{'alt_codon'}\t$vcf_data{$key}{'ref_aa_short'}$vcf_data{$key}{'aa_pos'}$vcf_data{$key}{'alt_aa_short'}\t$vcf_data{$key}{'snp_type'}"
#             if $options{lvl} == 5 and $vcf_data{$key}{_type} eq "snp";

#         $snp_per_gene_h{ $vcf_data{$key}{locus} }{$key}
#             = "$vcf_data{$key}{snp}\t$vcf_data{$key}{'ref_codon'}/$vcf_data{$key}{'alt_codon'}\t$vcf_data{$key}{'ref_aa_short'}$vcf_data{$key}{'aa_pos'}$vcf_data{$key}{'alt_aa_short'}\t$vcf_data{$key}{'snp_type'}\t$vcf_data{$key}{tang}\t$vcf_data{$key}{gd}\t$vcf_data{$key}{blossum45}\t$vcf_data{$key}{pam30}\t$vcf_data{$key}{AACI}"
#             if $options{lvl} == 6 and $vcf_data{$key}{_type} eq "snp";
#         if ( $options{indel} ) {
#             $snp_per_gene_h{ $vcf_data{$key}{locus} }{$key}
#                 = "$vcf_data{$key}{'snp_gene_pos'}$vcf_data{$key}{'short_indel_note'}$vcf_data{$key}{'indel'}\t$vcf_data{$key}{'full_formated_indel'}\t$vcf_data{$key}{'frame_shift'}\t$vcf_data{$key}{'snp_type'}"
#                 if $options{lvl} == 5 and $vcf_data{$key}{_type} eq "indel";
#             $snp_per_gene_h{ $vcf_data{$key}{locus} }{$key}
#                 = "$vcf_data{$key}{'snp_gene_pos'}$vcf_data{$key}{'short_indel_note'}$vcf_data{$key}{'indel'}\t$vcf_data{$key}{'full_formated_indel'}\t$vcf_data{$key}{'frame_shift'}\t$vcf_data{$key}{'snp_type'}\t$vcf_data{$key}{tang}\t$vcf_data{$key}{gd}\t$vcf_data{$key}{blossum45}\t$vcf_data{$key}{pam30}\t$vcf_data{$key}{AACI}"
#                 if $options{lvl} == 6 and $vcf_data{$key}{_type} eq "indel";

#         }

#         #------------------------------------------------------------
#     }
#     if ( $options{lvl} == 5 or $options{lvl} == 6 ) {
#         foreach my $key ( sort keys %snp_per_gene_h ) {

#             my $key2_count = 0;
#             foreach my $key2 ( sort keys %{ $snp_per_gene_h{$key} } ) {
#                 $key2_count++;
#                 if ( $key2_count == 1 ) {
#                     print "$key" . "\t"
#                         . $database{$key}{'product'} . "\t"
#                         . "\n\t$key2\t$snp_per_gene_h{$key}{$key2}\n";

#                 }
#                 else {
#                     print "\t$key2\t$snp_per_gene_h{$key}{$key2}\n";
#                 }
#             }
#         }
#     }

#     my $elapsed = tv_interval( $t0, [gettimeofday] );
#     print "elapsed:\t$elapsed\n";

# }

# sub annotate_vcf {
#     my $locInputFile   = shift;
#     my $AACI          = "";
#     my $loc_AACI_flag = "----";
#     my $loc_AACI      = 0;
#     my $locRefGenomeName  = "";
#     my $DP            = "";
#     my $count         = 0;

# # print
# # "locus\tsnp_pos\tgene_pos\tcodons\taa_pos\ttype\tu_index\tGD\tscore\tDP\tproduct\n";
#     open( my $fh, "<", $locInputFile )
#         or croak "cannot open < $locInputFile: $!";

#     while (<$fh>) {
#         my $str = "$_";
#         next if $str =~ /^#/;
#         $str =~ /.*($RefGenomeName)/x;
#         $locRefGenomeName = $1;
#         $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;
#         if ( length($2) == 1 && length($3) == 1 ) {
#             my %LocusInfoReader = &get_locus_info( $1, $3 );
#             if ( $locRefGenomeName ne $RefGenomeName ) {
#                 print
#                     "REFERENCE GENOME: $RefGenomeName NOT FOUND! ($locRefGenomeName)\n";
#                 exit;
#             }

#  # print "$LocusInfoReader{'alt_aa_long'}\t!!!!!\n" if $LocusInfoReader{'alt_aa_long'} eq "BAD_CODON";
#             if ($LocusInfoReader{'locus'} ne ""

#                 # && $LocusInfoReader{'alt_aa_long'} ne 'STOP'
#                 # && $LocusInfoReader{'ref_aa_long'} ne 'STOP'
#                 && $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
#                 && $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON'
#                 )
#             {
#                 $DP = $4;
#                 if ( $LocusInfoReader{'snp_type'} eq "missense" ) {
#                     my $loc_tang_index
#                         = calcutalte_tang_index( $LocusInfoReader{'ref_aa_long'},
#                         $LocusInfoReader{'alt_aa_long'} );

#                     # my $loc_deltaH
#                     #     = calcutalte_deltaH( $LocusInfoReader{'ref_aa_long'},
#                     #     $LocusInfoReader{'alt_aa_long'} );
#                     my $loc_GD
#                         = calculate_grantham_matrix( $LocusInfoReader{'ref_aa_long'},
#                         $LocusInfoReader{'alt_aa_long'} );
#                     my $blosum
#                         = calculate_blosum45_matrix( $LocusInfoReader{'ref_aa_long'},
#                         $LocusInfoReader{'alt_aa_long'} );
#                     my $pam
#                         = &calculate_pam30_matrix( $LocusInfoReader{'ref_aa_long'},
#                         $LocusInfoReader{'alt_aa_long'} );

#                     # my $titv=&calcutalte_tv_ti($LocusInfoReader{'ref'},
#                     #     $LocusInfoReader{'alt'} );
#                     #  $titv="-" if $titv eq "";
#                     if ( $loc_tang_index < 0.5 ) {
#                         $loc_AACI++;
#                         substr $loc_AACI_flag, 0, 1, "+";
#                     }

#                     # if ( $loc_deltaH > 1.22 ) {
#                     #     $loc_AACI++;
#                     #     substr $loc_AACI_flag, 1, 1, "+";
#                     # }
#                     if ( $loc_GD > 100 ) {
#                         $loc_AACI++;
#                         substr $loc_AACI_flag, 1, 1, "+";
#                     }
#                     if ( $blosum < 0 ) {
#                         $loc_AACI++;
#                         substr $loc_AACI_flag, 2, 1, "+";
#                     }
#                     if ( $pam < 0 ) {
#                         $loc_AACI++;
#                         substr $loc_AACI_flag, 3, 1, "+";
#                     }

#                     $AACI
#                         = "$loc_tang_index\t$loc_GD\t$blosum\t$pam\t$loc_AACI_flag\t";
#                 }

#                 else {
#                     $AACI = "-\t-\t-\t-\t$loc_AACI_flag\t";
#                 }
#                 if ( $options{action} eq "annotation.1" ) {

#                     # if ( $loc_AACI_flag eq "+++" ) {
#                     print
#                         "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$AACI$LocusInfoReader{'product'}\n";

#                     # }
#                 }

#        # elsif ( $options{action} eq "annotation.2" ) {
#        #     print
#        #         "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t\t$LocusInfoReader{'product'}\n";
#        # }
#                 elsif ( $options{action} eq "annotation" ) {

#                     # print  color 'bold green';
#                     # if ( $loc_AACI_flag ne "---" ) {

#                     print
#                         "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$LocusInfoReader{'product'}\n";

#                     # }
#                 }
#                 elsif ( $options{action} eq "annotation.11" ) {
#                     print
#                         "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$AACI$LocusInfoReader{'product'}\n"
#                         if $loc_AACI_flag ne "----";
#                 }
#                 elsif ( $options{action} eq "annotation.12" ) {
#                     print
#                         "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$AACI$LocusInfoReader{'product'}\n"
#                         if $loc_AACI_flag eq "++++";
#                 }

#                 $loc_AACI      = 0;
#                 $loc_AACI_flag = "----";

#             }

#         }
#         elsif ( length($3) > length($2) and $options{indel} ) {
#             my %LocusInfoReader = &get_indel_info( $1, $2, $3 );
#             next if $LocusInfoReader{locus} eq "";
#             print $LocusInfoReader{locus} . "\t"
#                 . $LocusInfoReader{snp_genome_pos} . "\t"
#                 . $LocusInfoReader{snp_gene_pos} . "ins"
#                 . $LocusInfoReader{indel} . "\t"
#                 . $LocusInfoReader{ref} . "/"
#                 . $LocusInfoReader{formated_indel} . "\t"
#                 . $LocusInfoReader{frame_shift} . "\t"
#                 . $LocusInfoReader{snp_type}
#                 . "\t-\t-\t-\t-\t"
#                 . $LocusInfoReader{product} . "\n";

#   # print $LocusInfoReader{full_formated_indel} .  "\t-\t-\t-\t-\t" . $LocusInfoReader{product} ."\n";

# # my $alt=$3;
# # my $snp_pos=$1;
# # my $ref=$2;
# # my %LocusInfoReader = &get_locus_info($snp_pos);
# # next if $LocusInfoReader{locus} eq "";
# # $alt=~/^($ref)(.*)/;
# # my $ref_pos=$LocusInfoReader{snp_gene_pos}+length($ref);
# # my $alt_pos=$ref_pos+length($2);
# # print $LocusInfoReader{locus}."\t$snp_pos\t$ref_pos\_$alt_pos" . "ins$2\tref:$ref\talt:$alt\tb:$1\ta:$2\n";
# # # exit;
# # # my %LocusInfoReader = &get_locus_info( $1,$3);
# # # print
# # #             "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t\t\t\t\t\tinsertion\t\t\n";
#         }
#         elsif ( length($3) < length($2) and $options{indel} ) {
#             my %LocusInfoReader = &get_indel_info( $1, $2, $3 );
#             next if $LocusInfoReader{locus} eq "";
#             print $LocusInfoReader{locus} . "\t"
#                 . $LocusInfoReader{snp_genome_pos} . "\t"
#                 . $LocusInfoReader{snp_gene_pos} . "del"
#                 . $LocusInfoReader{indel} . "\t"
#                 . $LocusInfoReader{formated_indel} . "/"
#                 . $LocusInfoReader{ref} . "\t"
#                 . $LocusInfoReader{frame_shift} . "\t"
#                 . $LocusInfoReader{snp_type}
#                 . "\t-\t-\t-\t-\t"
#                 . $LocusInfoReader{product} . "\n";

#         }

#     }
#     close $fh;

# }

=======
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
sub get_intergen_info {

    # my $locSnpPosInGenome    = shift;
    # my $loc_alt        = shift;
    my ( $locSnpPosInGenome, $loc_alt, $enable_calcs )
        = @_;    #gets values of abs. snp pos,
    my $locGeneLocusName = "";
    my %LocusInfoProvider;
    my $snpPosInGene;
    my $typeOfSNP    = "";
    my $gene_length = 0;
    my ( $tang, $gd, $blossum45, $pam30 );
    my $loc_AACI_flag = "----";

    foreach my $key ( sort keys %intergenic )
    {            #looking for positions of intergenic snps

        if (  $locSnpPosInGenome > scalar $intergenic{$key}{'start'}
            && $locSnpPosInGenome < scalar $intergenic{$key}{'end'} )
        {
            $locGeneLocusName = $key;
            last;

        }

    }

    if ( $locGeneLocusName ne "" ) {
        my $locStartOfGene = $intergenic{$locGeneLocusName}{"start"};
        my $locEndOfGene   = $intergenic{$locGeneLocusName}{"end"};

      #     my $locStrandDirection       = $database{$locGeneLocusName}{"strand"};
      #     my $locTranslationTable = $database{$locGeneLocusName}{"transl_table"};
      #     my $loc_product      = $database{$locGeneLocusName}{"product"};
      #     my $loc_note         = $database{$locGeneLocusName}{"note"};
        my $loc_nuc_seq = $intergenic{$locGeneLocusName}{"sequence"};
        my $igen_genes  = $intergenic{$locGeneLocusName}{"igen_genes"};

        #     if ( $locStrandDirection == 1 ) {Rv1132
        $snpPosInGene = ( ( $locSnpPosInGenome - $locStartOfGene ) + 1 );

        #     }
        #     elsif ( $locStrandDirection == -1 ) {
        #         $snpPosInGene = ( ( $locEndOfGene - $locSnpPosInGenome ) + 1 );
        #         $loc_alt =~ tr/ACGT/TGCA/;
        #     }
        #     if ( $locEndOfGene > $locStartOfGene ) {
        $gene_length = ( $locEndOfGene - $locStartOfGene ) + 1;

        #     }
        #     elsif ( $locStartOfGene > $locEndOfGene ) {
        #         $gene_length = ( $locStartOfGene - $locEndOfGene ) + 1;
        #     }
        my $loc_ref
            = substr( $loc_nuc_seq, $snpPosInGene - 1, 1 );    #ref letter
        my $codon_nbr = int( ( $snpPosInGene - 1 ) / 3 ) + 1;  #codon position
        my @codons = unpack( "(A3)*", $loc_nuc_seq );          # nbr of codons
        my $res_codon = lc( $codons[ $codon_nbr - 1 ] );
        my $loc_codon_pos = ( ( $codon_nbr * 3 ) - $snpPosInGene );
        if ( $loc_codon_pos == 2 ) {
            $loc_codon_pos = 0;
        }
        elsif ( $loc_codon_pos == 0 ) {
            $loc_codon_pos = 2;
        }

        substr( $res_codon, $loc_codon_pos, 1 )
            = uc( substr( $res_codon, $loc_codon_pos, 1 ) );

        my $alt_res_codon = $res_codon;
        substr $alt_res_codon, $loc_codon_pos, 1, $loc_alt;

<<<<<<< HEAD
        # my $ref_aa = codon2aa($res_codon);
        # my $alt_aa = codon2aa($alt_res_codon);
        # if ( $ref_aa eq $alt_aa and $alt_aa ne "X" ) {
        #     $typeOfSNP = "synonymous";
        # }
        # elsif ( $ref_aa ne $alt_aa and $alt_aa ne "X" ) {
        #     $typeOfSNP = "missense";
        # }
        # elsif ( $ref_aa ne $alt_aa and $alt_aa eq "X" ) {
        #     $typeOfSNP = "nonsense";

        # }

        #     if ( $enable_calcs eq 1 and $typeOfSNP eq "missense" ) {
        #         $tang = calcutalte_tang_index( $ref_aa, $alt_aa );
        #         $gd = calculate_grantham_matrix( $ref_aa, $alt_aa );
        #         $blossum45 = calculate_blosum45_matrix( $ref_aa, $alt_aa );
        #         $pam30 = calculate_pam30_matrix( $ref_aa, $alt_aa );
        #         if ( $tang < TANG_TH ) {

        #             substr $loc_AACI_flag, 0, 1, "+";
        #         }

        #         if ( $gd > GD_TH ) {

        #             substr $loc_AACI_flag, 1, 1, "+";
        #         }
        #         if ( $blossum45 < BLOSSUM45_TH ) {

        #             substr $loc_AACI_flag, 2, 1, "+";
        #         }
        #         if ( $pam30 < PAM30_TH ) {

        #             substr $loc_AACI_flag, 3, 1, "+";
        #         }

        #     }
        #     elsif ( $enable_calcs eq 1 and $typeOfSNP eq "synonymous" ) {
        #         $tang          = "-";
        #         $gd            = "-";
        #         $blossum45     = "-";
        #         $pam30         = "-";
        #         $loc_AACI_flag = "----";
        #     }
        #     elsif ( $enable_calcs eq 1 and $typeOfSNP eq "nonsense" ) {
        #         $tang          = "-";
        #         $gd            = "-";
        #         $blossum45     = "-";
        #         $pam30         = "-";
        #         $loc_AACI_flag = "----";
        #     }

        %LocusInfoProvider = (
            "locus" => $locGeneLocusName,
=======
        %locus_info_h = (
            "locus" => $loc_locus_name,
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e

            "snp_genome_pos" => $locSnpPosInGenome,
            "start"          => $locStartOfGene,
            "end"            => $locEndOfGene,

            # "strand"         => $locStrandDirection,
            # "product"        => $loc_product,
            "sequence" => $loc_nuc_seq,

            # "snp_gene_pos"   => $snpPosInGene,
            "alt" => $loc_alt,
            "snp" => "$locSnpPosInGenome$loc_ref>$loc_alt",
            "ref" => $loc_ref,

            # "aa_pos"    => $codon_nbr,
            "ref_codon" => $res_codon,

            # "snp_type"       => $typeOfSNP,
            # "ref_aa_long"    => $ref_aa,
            # "alt_aa_long"    => $alt_aa,
            # "ref_aa_short"   => &aa_decode($ref_aa),
            # "alt_aa_short"   => &aa_decode($alt_aa),
            "alt_codon"  => $alt_res_codon,
            "igen_genes" => $igen_genes,
            "product"    => "intergenic region"

                # "note"           => $loc_note,
                # "gene_length"    => $gene_length
                # "protein_length" => scalar(@codons) - 1,
                # "tang"           => $tang,
                # "gd"             => $gd,
                # "blossum45"      => $blossum45,
                # "pam30"          => $pam30,
                # "AACI"           => $loc_AACI_flag,
        );
        return %LocusInfoProvider;

        # $loc_AACI_flag = "----";

    }
    else {
        return;
    }

    # }
}

sub open_single_vcf {
    my ( $locInputFile, $enable_calcs ) = @_;
    my $locRefGenomeName = "";
    my %results;
    # my $loc_type = "";
    my %locus_info;
    my $loc_AACI_flag = "----";
    my %file_hash;
    my $vcfParserQuery = qr/^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/;
    my $t0    = [gettimeofday];
    open( my $fh, "<", $locInputFile )
        or croak "cannot open < $locInputFile: $!";
        #**** HEADERS DISABLED
    # if ($options{lvl} == 1 && $options{debug}) {
    # print "locus\tpos.\tsnp\tcodon\tAA\ttype\t\tdesc.\n";
    # }
    # elsif ($options{lvl} == 2 || $options{lvl} == 3 || $options{lvl} == 4) {
    #   print "locus\tpos.\tsnp\tcodon\tAA\ttype\ttang\tGD\tblossumm45\tpam30\tscore\tdesc.\n" if $options{debug};
    # }
    while (<$fh>) {
        my $str = "$_";
        next if $str =~ /^#/;
        $str =~ /.*($RefGenomeName)/x;
        $locRefGenomeName = $1;
        if ( $locRefGenomeName ne $RefGenomeName ) {
            print "REFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
            exit;
        }

        $str =~$vcfParserQuery;
        if ( length($2) == 1 and length($3) == 1 ) {

            &print_locus_info( $1, $3 );
        }
        elsif ( length($2) != 1 and length($3) != 1 and $options{indel} ) {
            &print_indel_info( $1, $2, $3 );
        }
    }
    close $fh;

    &print_grouped_info;

    # print Dumper @to_circos_results;
    # exit;

    my $elapsed = tv_interval( $t0, [gettimeofday] );
    print "-" x 50 . "\n\tTime elapsed:\t$elapsed\n"
        if $options{debug} and $options{vcf} ne "all";

}

sub print_locus_info {

    # my $locSnpPosInGenome    = shift;
    # my $loc_alt        = shift;
    my ( $locSnpPosInGenome, $loc_alt ) = @_;
    my $locGeneLocusName = "";
    my %LocusInfoProvider;
    my $snpPosInGene;
    my $typeOfSNP    = "";
    my $gene_length = 0;
    my ( $tang, $gd, $blossum45, $pam30 );
    my $loc_AACI_flag = "----";
    my $snp           = "";
    my $ref_aa_short;
    my $alt_aa_short;

    # my %snp_per_gene_h;
    my $enable_calcs = "";

    if ( $options{lvl} eq 1 || $options{lvl} eq 5 ) {
        $enable_calcs = 0;
    }
    else {
        $enable_calcs = 1;
    }

    foreach my $key ( keys %database ) {
        if (   $locSnpPosInGenome >= scalar $database{$key}{'start'}
            && $locSnpPosInGenome <= scalar $database{$key}{'end'} )
        {
            $locGeneLocusName = $key;
            last;
        }

    }

    if ( $locGeneLocusName ne "" ) {
        my $locStartOfGene        = $database{$locGeneLocusName}{"start"};
        my $locEndOfGene          = $database{$locGeneLocusName}{"end"};
        my $locStrandDirection       = $database{$locGeneLocusName}{"strand"};
        my $locTranslationTable = $database{$locGeneLocusName}{"transl_table"};
        my $loc_product      = $database{$locGeneLocusName}{"product"};
        my $loc_note         = $database{$locGeneLocusName}{"note"};
        my $loc_nuc_seq      = $database{$locGeneLocusName}{"sequence"};

        if ( $locStrandDirection == 1 ) {
            $snpPosInGene = ( ( $locSnpPosInGenome - $locStartOfGene ) + 1 );
        }
        elsif ( $locStrandDirection == -1 ) {
            $snpPosInGene = ( ( $locEndOfGene - $locSnpPosInGenome ) + 1 );
            $loc_alt =~ tr/ACGT/TGCA/;
        }
        if ( $locEndOfGene > $locStartOfGene ) {
            $gene_length = ( $locEndOfGene - $locStartOfGene ) + 1;
        }
        elsif ( $locStartOfGene > $locEndOfGene ) {
            $gene_length = ( $locStartOfGene - $locEndOfGene ) + 1;
        }
        my $loc_ref = substr( $loc_nuc_seq, $snpPosInGene - 1, 1 )
            ;    #reference target letter
        my $codon_nbr = int( ( $snpPosInGene - 1 ) / 3 )
            + 1;    #number of target codon from start
        my @codons = unpack( "(A3)*", $loc_nuc_seq );     # codons of gene
        my $res_codon = lc( $codons[ $codon_nbr - 1 ] )
            ;    # letters from target codon in low case
        my $loc_codon_pos = ( ( $codon_nbr * 3 ) - $snpPosInGene )
            ;    # position letter in codon (0,1,2)
         # print "$loc_ref\n$codon_nbr\n$res_codon\n$loc_codon_pos\n" . scalar(@codons) . "\n";
         # exit;

        if ( $loc_codon_pos == 2 ) {
            $loc_codon_pos = 0;
        }
        elsif ( $loc_codon_pos == 0 ) {
            $loc_codon_pos = 2;
        }

        substr( $res_codon, $loc_codon_pos, 1 )
            = uc( substr( $res_codon, $loc_codon_pos, 1 ) );
        my $alt_res_codon = $res_codon;
        substr $alt_res_codon, $loc_codon_pos, 1, $loc_alt;
        my $ref_aa = codon2aa($res_codon);
        my $alt_aa = codon2aa($alt_res_codon);
        if (    $ref_aa eq $alt_aa
            and $alt_aa ne "X"
            and $alt_aa ne "BAD_CODON"
            and $ref_aa ne "BAD_CODON" )
        {
            $typeOfSNP = "synonymous";
        }
        elsif ( $ref_aa ne $alt_aa
            and $alt_aa ne "X"
            and $alt_aa ne "BAD_CODON"
            and $ref_aa ne "BAD_CODON" )
        {
            $typeOfSNP = "missense";
        }
        elsif ( $ref_aa ne $alt_aa
            and $alt_aa eq "X"
            and $alt_aa ne "BAD_CODON"
            and $ref_aa ne "BAD_CODON" )
        {
            $typeOfSNP = "nonsense";

        }

        if ( $enable_calcs eq 1 and $typeOfSNP eq "missense" ) {
            $tang = calcutalte_tang_index( $ref_aa, $alt_aa );
            $gd = calculate_grantham_matrix( $ref_aa, $alt_aa );
            $blossum45 = calculate_blosum45_matrix( $ref_aa, $alt_aa );
            $pam30 = calculate_pam30_matrix( $ref_aa, $alt_aa );
            if ( $tang < TANG_TH ) {

                substr $loc_AACI_flag, 0, 1, "+";
            }

            if ( $gd > GD_TH ) {

                substr $loc_AACI_flag, 1, 1, "+";
            }
            if ( $blossum45 < BLOSSUM45_TH ) {

                substr $loc_AACI_flag, 2, 1, "+";
            }
            if ( $pam30 < PAM30_TH ) {

                substr $loc_AACI_flag, 3, 1, "+";
            }

        }
        elsif ( $enable_calcs eq 1 and $typeOfSNP eq "synonymous" ) {
            $tang          = "-";
            $gd            = "-";
            $blossum45     = "-";
            $pam30         = "-";
            $loc_AACI_flag = "----";
        }
        elsif ( $enable_calcs eq 1 and $typeOfSNP eq "nonsense" ) {
            $tang          = "-";
            $gd            = "-";
            $blossum45     = "-";
            $pam30         = "-";
            $loc_AACI_flag = "----";
        }

        $snp          = "$snpPosInGene$loc_ref>$loc_alt";
        $ref_aa_short = &aa_decode($ref_aa);
        $alt_aa_short = &aa_decode($alt_aa);

        # if ($options{lvl} ==7){
        # my $sub_seq = $loc_nuc_seq;
        # substr( $sub_seq, $snpPosInGene - 1, 1 ) = "[$loc_ref>$loc_alt]";
        # my $sub_res=substr($sub_seq, $snpPosInGene - 6,16);
        # }

        print
            "$locGeneLocusName\t$locSnpPosInGenome\t$snp\t$res_codon/$alt_res_codon\t$ref_aa_short$codon_nbr$alt_aa_short\t$typeOfSNP\t$loc_product\n"
            if $options{lvl} == 1;

        print
            "$locGeneLocusName\t$locSnpPosInGenome\t$snp\t$res_codon/$alt_res_codon\t$ref_aa_short$codon_nbr$alt_aa_short\t$typeOfSNP\t$tang\t$gd\t$blossum45\t$pam30\t$loc_AACI_flag\t$loc_product\n"
            if $options{lvl} == 2;
        print
            "$locGeneLocusName\t$locSnpPosInGenome\t$snp\t$res_codon/$alt_res_codon\t$ref_aa_short$codon_nbr$alt_aa_short\t$typeOfSNP\t$tang\t$gd\t$blossum45\t$pam30\t$loc_AACI_flag\t$loc_product\n"
            if $options{lvl} == 3
            and $typeOfSNP eq "missense";

        # and $loc_AACI_flag ne "----";

        print
            "$locGeneLocusName\t$locSnpPosInGenome\t$snp\t$res_codon/$alt_res_codon\t$ref_aa_short$codon_nbr$alt_aa_short\t$typeOfSNP\t$tang\t$gd\t$blossum45\t$pam30\t$loc_AACI_flag\t$loc_product\n"
            if $options{lvl} == 4
            and $loc_AACI_flag eq "++++";

        $snp_per_gene_h{$locGeneLocusName}{$locSnpPosInGenome}
            = "\t$snp\t$res_codon/$alt_res_codon\t$ref_aa_short$codon_nbr$alt_aa_short\t$typeOfSNP"
            if $options{lvl} == 5;

        $snp_per_gene_h{$locGeneLocusName}{$locSnpPosInGenome}
            = "\t$snp\t$res_codon/$alt_res_codon\t$ref_aa_short$codon_nbr$alt_aa_short\t$typeOfSNP\t$tang\t$gd\t$blossum45\t$pam30\t$loc_AACI_flag"
            if $options{lvl} == 6;

<<<<<<< HEAD
        # if ( $options{lvl} == 5 or $options{lvl} == 6 ) {

        #     foreach my $key ( sort keys %snp_per_gene_h ) {

        #         my $key2_count = 0;

        #         foreach my $key2 ( sort keys %{ $snp_per_gene_h{$key} } ) {
        #             $key2_count++;
        #             if ( $key2_count == 1 ) {
        #                 print "$key" . "\t"
        #                     . $database{$key}{'product'} . "\t"
        #                     . "\n\t$key2$snp_per_gene_h{$key}{$key2}\n";

        #             }
        #             else {
        #                 print "\t$key2$snp_per_gene_h{$key}{$key2}\n";
        #             }
        #         }

        #     }
        # }

        # %LocusInfoProvider = (
        #     "locus"          => $locGeneLocusName,
        #     "snp_genome_pos" => $locSnpPosInGenome,
        #     "start"          => $locStartOfGene,
        #     "end"            => $locEndOfGene,
        #     "strand"         => $locStrandDirection,
        #     "product"        => $loc_product,
        #     "sequence"       => $loc_nuc_seq,
        #     "snp_gene_pos"   => $snpPosInGene,
        #     "alt"            => $loc_alt,
        #     "snp"            => "$snpPosInGene$loc_ref>$loc_alt",
        #     "ref"            => $loc_ref,
        #     "aa_pos"         => $codon_nbr,
        #     "ref_codon"      => $res_codon,
        #     "snp_type"       => $typeOfSNP,
        #     "ref_aa_long"    => $ref_aa,
        #     "alt_aa_long"    => $alt_aa,
        #     "ref_aa_short"   => &aa_decode($ref_aa),
        #     "alt_aa_short"   => &aa_decode($alt_aa),
        #     "alt_codon"      => $alt_res_codon,
        #     "product"        => $loc_product,
        #     "note"           => $loc_note,
        #     "gene_length"    => $gene_length,
        #     "protein_length" => scalar(@codons) - 1,
        #     "tang"    => $tang,
        #     "gd"    => $gd,
        #     "blossum45"    => $blossum45,
        #     "pam30"    => $pam30,
        #     "AACI"    => $loc_AACI_flag,

        # );
        # return %LocusInfoProvider;
=======
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
        $loc_AACI_flag = "----";

    }
    elsif ( $locGeneLocusName eq "" and $options{intergenic} ) {
        my %LocusInfoReader = &get_intergen_info( $locSnpPosInGenome, $loc_alt );
        if ( $LocusInfoReader{locus} ne "" ) {
            print
                "$LocusInfoReader{igen_genes}\t$locSnpPosInGenome\t$LocusInfoReader{snp}\t$LocusInfoReader{ref_codon}\/$LocusInfoReader{alt_codon}\t-\t-\t$LocusInfoReader{product} ($LocusInfoReader{igen_genes})\n";

        }
    }

    # else {
    #     return;
    # }

}

sub print_grouped_info {

    if ( $options{lvl} == 5 or $options{lvl} == 6 ) {

        foreach my $key ( sort keys %snp_per_gene_h ) {

            my $key2_count = 0;

            foreach my $key2 ( sort keys %{ $snp_per_gene_h{$key} } ) {
                $key2_count++;
                if ( $key2_count == 1 ) {
                    print "$key" . "\t"
                        . $database{$key}{'product'} . "\t"
                        . "\n\t$key2$snp_per_gene_h{$key}{$key2}\n";

                }
                else {
                    print "\t$key2$snp_per_gene_h{$key}{$key2}\n";
                }
            }

        }
    }

    undef %snp_per_gene_h;
}

sub print_indel_info {

    # my $locSnpPosInGenome    = shift;
    # my $loc_alt        = shift;
    my ( $locSnpPosInGenome, $ref, $alt ) = @_;
    my $locGeneLocusName = "";

    # my %snp_per_gene_h;
    my $snpPosInGene;
    my $typeOfSNP            = "";
    my $gene_length         = 0;
    my $indel               = "";
    my $indel_bp_count      = 0;
    my $frame_shift         = "";
    my $formated_indel      = "";
    my $short_indel_note    = "";
    my $full_formated_indel = "";
    my $loc_AACI_flag       = "----";
    my ( $tang, $gd, $blossum45, $pam30 );

    # my $common="";
    foreach my $key ( keys %database ) {
        if (   $locSnpPosInGenome >= scalar $database{$key}{'start'}
            && $locSnpPosInGenome <= scalar $database{$key}{'end'} )
        {
            $locGeneLocusName = $key;
            last;
        }

    }

    if ( $locGeneLocusName ne "" ) {
        my $locStartOfGene        = $database{$locGeneLocusName}{"start"};
        my $locEndOfGene          = $database{$locGeneLocusName}{"end"};
        my $locStrandDirection       = $database{$locGeneLocusName}{"strand"};
        my $locTranslationTable = $database{$locGeneLocusName}{"transl_table"};
        my $loc_product      = $database{$locGeneLocusName}{"product"};
        my $loc_note         = $database{$locGeneLocusName}{"note"};
        my $loc_nuc_seq      = $database{$locGeneLocusName}{"sequence"};

        if ( $locStrandDirection == 1 ) {
            $snpPosInGene = ( ( $locSnpPosInGenome - $locStartOfGene ) + 1 );
        }
        elsif ( $locStrandDirection == -1 ) {
            $snpPosInGene = ( ( $locEndOfGene - $locSnpPosInGenome ) + 1 );

            $alt =~ tr/ACGT/TGCA/;
            $ref =~ tr/ACGT/TGCA/;

        }
        if ( $locEndOfGene > $locStartOfGene ) {
            $gene_length = ( $locEndOfGene - $locStartOfGene ) + 1;
        }
        elsif ( $locStartOfGene > $locEndOfGene ) {
            $gene_length = ( $locStartOfGene - $locEndOfGene ) + 1;
        }

        if ( length($alt) > length($ref) ) {
            $typeOfSNP = "insertion";

<<<<<<< HEAD
            # $alt =~ /^(?<left>$ref)\W+(?<right>.*)/;
            # $alt =~ /^\S+\W+\d+\W+(?<left>\w+)\W+(?<right>\w+)/;
            # ^\S+\W+\d+\W+(\w+)\W+(\w+)

            # $indel=uc($2);

# if ($locSnpPosInGenome eq "2704884")
# {
# print $+{left} . "\t" . $+{right} . "\n";
#     exit;
# }
# next  if $+{left} eq "" || $+{right} eq "";
# print "$options{vcf}\t$locGeneLocusName\t$alt\tleft:$+{left}\tright:$+{right}\n" if  $+{left} eq "" || $+{right} eq "";

            # $indel = uc( $+{right} );
=======
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
            $indel = uc($alt);
            if ( $locStrandDirection == -1 ) {
                $indel = reverse($indel);
                $ref   = reverse($ref);
            }

            # $indel_bp_count=length($indel)/3;
            $indel_bp_count = length($indel) % 3;

            # print length($indel). "\t$indel_bp_count\n";

            if ( $indel_bp_count == 0 ) {
                $frame_shift = "no_frame_shift";
            }
            elsif ( $indel_bp_count != 0 ) {
                $frame_shift = "frame_shift";
            }

            $formated_indel = lc($ref) . uc($indel);

        }
        elsif ( length($alt) < length($ref) ) {
            $typeOfSNP = "deletion";

            # $ref =~ /^(?<left>$alt)(?<right>.*)/;

            # $locGeneLocusName = "" if $+{left} eq "" || $+{right} eq "";
            # $indel = uc( $+{left} );
            $indel = $ref;
            if ( $locStrandDirection == -1 ) {

                $ref =~ tr/ACGTacgt/TGCAtgca/;
                $indel =~ tr/ACGTacgt/TGCAtgca/;
                $indel =~ reverse($indel);
                $ref =~ reverse($ref);
            }

            $indel_bp_count = length($indel) % 3;

            if ( $indel_bp_count == 0 ) {
                $frame_shift = "no_frame_shift";
            }
            elsif ( $indel_bp_count != 0 ) {
                $frame_shift = "frame_shift";
            }
            $formated_indel = lc($ref) . uc($indel);

        }

<<<<<<< HEAD
        # %LocusInfoProvider = (
        #     "locus"            => $locGeneLocusName,
        #     "snp_genome_pos"   => $locSnpPosInGenome,
        #     "start"            => $locStartOfGene,
        #     "end"              => $locEndOfGene,
        #     "strand"           => $locStrandDirection,
        #     "product"          => $loc_product,
        #     "sequence"         => $loc_nuc_seq,
        #     "snp_gene_pos"     => $snpPosInGene,
        #     "snp_type"         => $typeOfSNP,
        #     "note"             => $loc_note,
        #     "gene_length"      => $gene_length,
        #     "indel"            => $indel,
        #     "ref"              => lc($ref),
        #     "formated_indel"   => $formated_indel,
        #     "frame_shift"      => $frame_shift,
        #     "short_indel_note" => substr( $typeOfSNP, 0, 3 )

        #         # "protein_length" => scalar(@codons) - 1,

        # );

        if ( $typeOfSNP eq "insertion" ) {
=======
        if ( $snp_type eq "insertion" ) {
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
            $full_formated_indel
                = lc($ref) . "/" . lc($ref) . uc($indel);

        }
        elsif ( $typeOfSNP eq "deletion" ) {
            $full_formated_indel
                = lc($ref) . uc($indel) . "/" . lc($ref);

        }

        $tang          = "-";
        $gd            = "-";
        $blossum45     = "-";
        $pam30         = "-";
        $loc_AACI_flag = "----";
        print "$locGeneLocusName\t$locSnpPosInGenome\t$snpPosInGene"
            . substr( $typeOfSNP, 0, 3 )
            . "$indel\t$full_formated_indel\t$frame_shift\t$typeOfSNP\t$loc_product\n"
            if $options{lvl} == 1;
        print "$locGeneLocusName\t$locSnpPosInGenome\t$snpPosInGene"
            . substr( $typeOfSNP, 0, 3 )
            . "$indel\t$full_formated_indel\t$frame_shift\t$typeOfSNP\t$tang\t$gd\t$blossum45\t$pam30\t$loc_AACI_flag\t$loc_product\n"
            if $options{lvl} == 2;

        $snp_per_gene_h{$locGeneLocusName}{$locSnpPosInGenome}
            = "\t$snpPosInGene"
            . substr( $typeOfSNP, 0, 3 )
            . "$indel\t$frame_shift\t$typeOfSNP"
            if $options{lvl} == 5;

        $snp_per_gene_h{$locGeneLocusName}{$locSnpPosInGenome}
            = "\t$snpPosInGene"
            . substr( $typeOfSNP, 0, 3 )
            . "$indel\t$full_formated_indel\t$frame_shift\t$typeOfSNP\t$tang\t$gd\t$blossum45\t$pam30\t$loc_AACI_flag"
            if $options{lvl} == 6;

<<<<<<< HEAD
        # if ( $options{lvl} == 5 or $options{lvl} == 6 ) {
        #     foreach my $key ( sort keys %snp_per_gene_h ) {

        #         my $key2_count = 0;
        #         foreach my $key2 ( sort keys %{ $snp_per_gene_h{$key} } ) {
        #             $key2_count++;
        #             if ( $key2_count == 1 ) {
        #                 print "$key" . "\t"
        #                     . $database{$key}{'product'} . "\t"
        #                     . "\n\t$key2\t$snp_per_gene_h{$key}{$key2}\n";

        #             }
        #             else {
        #                 print "\t$key2\t$snp_per_gene_h{$key}{$key2}\n";
        #             }
        #         }
        #     }
        # }

        #     "locus"            => $locGeneLocusName,
        #     "snp_genome_pos"   => $locSnpPosInGenome,
        #     "start"            => $locStartOfGene,
        #     "end"              => $locEndOfGene,
        #     "strand"           => $locStrandDirection,
        #     "product"          => $loc_product,
        #     "sequence"         => $loc_nuc_seq,
        #     "snp_gene_pos"     => $snpPosInGene,
        #     "snp_type"         => $typeOfSNP,
        #     "note"             => $loc_note,
        #     "gene_length"      => $gene_length,
        #     "indel"            => $indel,
        #     "ref"              => lc($ref),
        #     "formated_indel"   => $formated_indel,
        #     "frame_shift"      => $frame_shift,
        #     "short_indel_note" => substr( $typeOfSNP, 0, 3 )
        # # return %LocusInfoProvider;
=======
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
        $loc_AACI_flag = "----";

    }

    # else {
    #     return;
    # }
}

sub get_indel_info {

    # my $locSnpPosInGenome    = shift;
    # my $loc_alt        = shift;
    my ( $locSnpPosInGenome, $ref, $alt ) = @_;
    my $locGeneLocusName = "";
    my %snp_per_gene_h;
    my $snpPosInGene;
    my $typeOfSNP            = "";
    my $gene_length         = 0;
    my $indel               = "";
    my $indel_bp_count      = 0;
    my $frame_shift         = "";
    my $formated_indel      = "";
    my $short_indel_note    = "";
    my $full_formated_indel = "";
    my $loc_AACI_flag       = "----";
    my ( $tang, $gd, $blossum45, $pam30 );
    my %LocusInfoProvider;

    # my $common="";
    foreach my $key ( sort keys %database ) {
        if (   $locSnpPosInGenome >= scalar $database{$key}{'start'}
            && $locSnpPosInGenome <= scalar $database{$key}{'end'} )
        {
            $locGeneLocusName = $key;
            last;
        }

    }

    if ( $locGeneLocusName ne "" ) {
        my $locStartOfGene        = $database{$locGeneLocusName}{"start"};
        my $locEndOfGene          = $database{$locGeneLocusName}{"end"};
        my $locStrandDirection       = $database{$locGeneLocusName}{"strand"};
        my $locTranslationTable = $database{$locGeneLocusName}{"transl_table"};
        my $loc_product      = $database{$locGeneLocusName}{"product"};
        my $loc_note         = $database{$locGeneLocusName}{"note"};
        my $loc_nuc_seq      = $database{$locGeneLocusName}{"sequence"};

        if ( $locStrandDirection == 1 ) {
            $snpPosInGene = ( ( $locSnpPosInGenome - $locStartOfGene ) + 1 );
        }
        elsif ( $locStrandDirection == -1 ) {
            $snpPosInGene = ( ( $locEndOfGene - $locSnpPosInGenome ) + 1 );

            $alt =~ tr/ACGT/TGCA/;
            $ref =~ tr/ACGT/TGCA/;

        }
        if ( $locEndOfGene > $locStartOfGene ) {
            $gene_length = ( $locEndOfGene - $locStartOfGene ) + 1;
        }
        elsif ( $locStartOfGene > $locEndOfGene ) {
            $gene_length = ( $locStartOfGene - $locEndOfGene ) + 1;
        }

        if ( length($alt) > length($ref) ) {
            $typeOfSNP = "insertion";

<<<<<<< HEAD
            # $alt =~ /^(?<left>$ref)(?<right>.*)/;

            # $indel=uc($2);

            # if ($locSnpPosInGenome eq "2704884")
            # {
            # print $+{left} . "\t" . $+{right} . "\n";
            #     exit;
            # }

            # $locGeneLocusName = "" if $+{left} eq "" || $+{right} eq "";

=======
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
            $indel = $alt;
            if ( $locStrandDirection == -1 ) {
                $indel = reverse($indel);
                $ref   = reverse($ref);
            }

            # $indel_bp_count=length($indel)/3;
            $indel_bp_count = length($indel) % 3;

            # print length($indel). "\t$indel_bp_count\n";

            if ( $indel_bp_count == 0 ) {
                $frame_shift = "no_frame_shift";
            }
            elsif ( $indel_bp_count != 0 ) {
                $frame_shift = "frame_shift";
            }

            $formated_indel = lc($ref) . uc($indel);

        }
        elsif ( length($alt) < length($ref) ) {
            $typeOfSNP = "deletion";

            # $ref =~ /^(?<left>$alt)(?<right>.*)/;

            # $locGeneLocusName = "" if $+{left} eq "" || $+{right} eq "";
            $indel = $ref;
            if ( $locStrandDirection == -1 ) {

                $ref =~ tr/ACGTacgt/TGCAtgca/;
                $indel =~ tr/ACGTacgt/TGCAtgca/;
                $indel =~ reverse($indel);
                $ref =~ reverse($ref);
            }

            $indel_bp_count = length($indel) % 3;

            if ( $indel_bp_count == 0 ) {
                $frame_shift = "no_frame_shift";
            }
            elsif ( $indel_bp_count != 0 ) {
                $frame_shift = "frame_shift";
            }
            $formated_indel = lc($ref) . uc($indel);

        }

        if ( $typeOfSNP eq "insertion" ) {
            $full_formated_indel
                = lc($ref) . "/" . lc($ref) . uc($indel);

        }
        elsif ( $typeOfSNP eq "deletion" ) {
            $full_formated_indel
                = lc($ref) . uc($indel) . "/" . lc($ref);

        }

        $tang          = "-";
        $gd            = "-";
        $blossum45     = "-";
        $pam30         = "-";
        $loc_AACI_flag = "----";

        my $standard_report
            = "$locGeneLocusName\t$locSnpPosInGenome\t$snpPosInGene"
            . substr( $typeOfSNP, 0, 3 )
            . "$indel\t$full_formated_indel\t$frame_shift\t$typeOfSNP\t$loc_product\n";
        my $full_report
            = "$locGeneLocusName\t$locSnpPosInGenome\t$snpPosInGene"
            . substr( $typeOfSNP, 0, 3 )
            . "$indel\t$full_formated_indel\t$frame_shift\t$typeOfSNP\t$tang\t$gd\t$blossum45\t$pam30\t$loc_AACI_flag\t$loc_product\n";

        %LocusInfoProvider = (
            "locus"            => $locGeneLocusName,
            "snp_genome_pos"   => $locSnpPosInGenome,
            "start"            => $locStartOfGene,
            "end"              => $locEndOfGene,
            "strand"           => $locStrandDirection,
            "product"          => $loc_product,
            "sequence"         => $loc_nuc_seq,
            "snp_gene_pos"     => $snpPosInGene,
            "snp_type"         => $typeOfSNP,
            "note"             => $loc_note,
            "gene_length"      => $gene_length,
            "indel"            => $indel,
            "ref"              => lc($ref),
            "formated_indel"   => $formated_indel,
            "frame_shift"      => $frame_shift,
            "short_indel_note" => substr( $typeOfSNP, 0, 3 ),
            "tang"             => $tang,
            "gd"               => $gd,
            "blossum45"        => $blossum45,
            "pam30"            => $pam30,
            "AACI"             => $loc_AACI_flag,
            "standard_report"  => $standard_report,
            "full_report"      => $full_report

                # "protein_length" => scalar(@codons) - 1,

        );
        return %LocusInfoProvider;

        #     "locus"            => $locGeneLocusName,
        #     "snp_genome_pos"   => $locSnpPosInGenome,
        #     "start"            => $locStartOfGene,
        #     "end"              => $locEndOfGene,
        #     "strand"           => $locStrandDirection,
        #     "product"          => $loc_product,
        #     "sequence"         => $loc_nuc_seq,
        #     "snp_gene_pos"     => $snpPosInGene,
        #     "snp_type"         => $typeOfSNP,
        #     "note"             => $loc_note,
        #     "gene_length"      => $gene_length,
        #     "indel"            => $indel,
        #     "ref"              => lc($ref),
        #     "formated_indel"   => $formated_indel,
        #     "frame_shift"      => $frame_shift,
        #     "short_indel_note" => substr( $typeOfSNP, 0, 3 )
        # # return %LocusInfoProvider;
        $loc_AACI_flag = "----";
    }

    else {
        return;
    }
}

<<<<<<< HEAD
# sub get_indel_info_old {

#     # my $locSnpPosInGenome    = shift;
#     # my $loc_alt        = shift;
#     my ( $locSnpPosInGenome, $ref, $alt ) = @_;
#     my $locGeneLocusName = "";
#     my %LocusInfoProvider;
#     my $snpPosInGene;
#     my $typeOfSNP         = "";
#     my $gene_length      = 0;
#     my $indel            = "";
#     my $indel_bp_count   = 0;
#     my $frame_shift      = "";
#     my $formated_indel   = "";
#     my $short_indel_note = "";

#     # my $common="";
#     foreach my $key ( keys %database ) {
#         if (   $locSnpPosInGenome > scalar $database{$key}{'start'}
#             && $locSnpPosInGenome < scalar $database{$key}{'end'} )
#         {
#             $locGeneLocusName = $key;
#             last;
#         }

#     }

#     if ( $locGeneLocusName ne "" ) {
#         my $locStartOfGene        = $database{$locGeneLocusName}{"start"};
#         my $locEndOfGene          = $database{$locGeneLocusName}{"end"};
#         my $locStrandDirection       = $database{$locGeneLocusName}{"strand"};
#         my $locTranslationTable = $database{$locGeneLocusName}{"transl_table"};
#         my $loc_product      = $database{$locGeneLocusName}{"product"};
#         my $loc_note         = $database{$locGeneLocusName}{"note"};
#         my $loc_nuc_seq      = $database{$locGeneLocusName}{"sequence"};

#         if ( $locStrandDirection == 1 ) {
#             $snpPosInGene = ( ( $locSnpPosInGenome - $locStartOfGene ) + 1 );
#         }
#         elsif ( $locStrandDirection == -1 ) {
#             $snpPosInGene = ( ( $locEndOfGene - $locSnpPosInGenome ) + 1 );

#             $alt =~ tr/ACGT/TGCA/;
#             $ref =~ tr/ACGT/TGCA/;

#         }
#         if ( $locEndOfGene > $locStartOfGene ) {
#             $gene_length = ( $locEndOfGene - $locStartOfGene ) + 1;
#         }
#         elsif ( $locStartOfGene > $locEndOfGene ) {
#             $gene_length = ( $locStartOfGene - $locEndOfGene ) + 1;
#         }

#         if ( length($alt) > length($ref) ) {
#             $typeOfSNP = "insertion";
#             $alt =~ /^(?<left>$ref)(?<right>.*)/;

#             # $indel=uc($2);

#             # if ($locSnpPosInGenome eq "2704884")
#             # {
#             # print $+{left} . "\t" . $+{right} . "\n";
#             #     exit;
#             # }
#             $locGeneLocusName = "" if $+{left} eq "" || $+{right} eq "";
#             $indel = uc( $+{right} );
#             if ( $locStrandDirection == -1 ) {
#                 $indel = reverse($indel);
#                 $ref   = reverse($ref);
#             }

#             # $indel_bp_count=length($indel)/3;
#             $indel_bp_count = length($indel) % 3;

#             # print length($indel). "\t$indel_bp_count\n";

#             if ( $indel_bp_count == 0 ) {
#                 $frame_shift = "no_frame_shift";
#             }
#             elsif ( $indel_bp_count != 0 ) {
#                 $frame_shift = "frame_shift";
#             }

#             $formated_indel = lc($ref) . uc($indel);

#         }
#         elsif ( length($alt) < length($ref) ) {
#             $typeOfSNP = "deletion";

#             $ref =~ /^(?<left>$alt)(?<right>.*)/;

#             $locGeneLocusName = "" if $+{left} eq "" || $+{right} eq "";
#             $indel = uc( $+{left} );
#             if ( $locStrandDirection == -1 ) {
#                 $indel = reverse($indel);
#                 $ref   = reverse($ref);
#             }

#             $indel_bp_count = length($indel) % 3;

#             if ( $indel_bp_count == 0 ) {
#                 $frame_shift = "no_frame_shift";
#             }
#             elsif ( $indel_bp_count != 0 ) {
#                 $frame_shift = "frame_shift";
#             }
#             $formated_indel = lc($ref) . uc($indel);

#         }

#         %LocusInfoProvider = (
#             "locus"            => $locGeneLocusName,
#             "snp_genome_pos"   => $locSnpPosInGenome,
#             "start"            => $locStartOfGene,
#             "end"              => $locEndOfGene,
#             "strand"           => $locStrandDirection,
#             "product"          => $loc_product,
#             "sequence"         => $loc_nuc_seq,
#             "snp_gene_pos"     => $snpPosInGene,
#             "snp_type"         => $typeOfSNP,
#             "note"             => $loc_note,
#             "gene_length"      => $gene_length,
#             "indel"            => $indel,
#             "ref"              => lc($ref),
#             "formated_indel"   => $formated_indel,
#             "frame_shift"      => $frame_shift,
#             "short_indel_note" => substr( $typeOfSNP, 0, 3 )

#                 # "protein_length" => scalar(@codons) - 1,

#         );
#         if ( $typeOfSNP eq "insertion" ) {
#             $LocusInfoProvider{full_formated_indel}
#                 = lc($ref) . "/" . lc($ref) . uc($indel);

#         }
#         elsif ( $typeOfSNP eq "deletion" ) {
#             $LocusInfoProvider{full_formated_indel}
#                 = lc($ref) . uc($indel) . "/" . lc($ref);

#         }

#         return %LocusInfoProvider;
#     }
#     else {
#         return;
#     }
# }

=======
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub use_snp_classification_intergenic {
    my $snp_class_name = shift;
    my @files          = glob('*.vcf');
    my $class_found    = 0;
    my $class_query    = "";
    my $query_locus    = "";
    my $query_snp      = "";
    my @class_a;
    my $snp_count        = 0;
    my $tab              = "";
    my $res              = "";
    my $locRefGenomeName     = "";
    my $loc_snp_notation = "";
    my $position_type    = "";

    foreach my $file (@files) {
        print "$file";
        my @res_a;
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*(?<ref_genome>$RefGenomeName)/x;
            if (   $+{ref_genome} eq "$RefGenomeName"
                && $+{ref_genome} ne "" )
            {
                $str =~ /^\S+\W(\d+)/x;
                $query_locus = "IGPOS$1";

                # print "$query_locus\n";
                # print $snp_class_name . "\n";
                if ( exists $snp_igpos_list_h{$snp_class_name}{$query_locus} )
                {

                    $loc_snp_notation
                        = $snp_igpos_list_h{$snp_class_name}{$query_locus}
                        {'note'};
                    push( @res_a, $loc_snp_notation );
                    $res = join( ",", sort @res_a );
                    $class_found = 1;
                }
            }
            elsif ($+{ref_genome} ne "$RefGenomeName"
                && $+{ref_genome} ne "" )
            {
                print "\t\tREFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
                last;
            }

<<<<<<< HEAD
     #         $locRefGenomeName = $+{ref_genome};
     #         $str =~ /^\S+\W+(?<genome_pos>\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

            #         # $position_type= &check_position($+{genome_pos});
            #         if ( length($2) == 1 && length($3) == 1 ) {

#             # print "!\n" if ($1 =="1673425");
#             # print &check_position($1). "\n";
#             if ( $locRefGenomeName ne $RefGenomeName ) {
#                 print
#                     "\t\tREFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
#                 last;
#             }

            #             print $+{genome_pos} . "\n";

            #             # $class_query
            #             #     = uc( $LocusInfoReader{'locus'} . "_" . $LocusInfoReader{'snp'} );

            #         }
=======
            #         $loc_ref_name = $+{ref_genome};
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
        }
        if ( $class_found == 0 ) {
            print "\t---";
        }
        $class_query = "";
        $file        = "";
        $class_found = 0;
        $snp_count   = 0;
        print "\t$res\n";
        $res = "";
        close $fh;

    }

}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub use_snp_classification_coding {
    my @files       = glob('*.vcf');
    my $class_found = 0;
    my $class_query = "";
    my $query_locus = "";
    my $query_snp   = "";
    my @class_a;
    my $snp_count        = 0;
    my $tab              = "";
    my $res              = "";
    my $locRefGenomeName     = "";
    my $loc_snp_notation = "";
    my $position_type    = "";

    foreach my $file (@files) {
        print "$file";
        my @res_a;
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*(?<ref_genome>$RefGenomeName)/x;
            $locRefGenomeName = $+{ref_genome};
            $str =~ /^\S+\W+(?<genome_pos>\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

            # $position_type= &check_position($+{genome_pos});
            if ( length($2) == 1 && length($3) == 1 ) {

                # print "!\n" if ($1 =="1673425");
                # print &check_position($1). "\n";
                my %LocusInfoReader = &get_locus_info( $+{genome_pos}, $3 );
                if ( $locRefGenomeName ne $RefGenomeName ) {
                    print
                        "\t\tREFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
                    last;
                }
                if ($LocusInfoReader{'locus'} ne ""

                    # && $LocusInfoReader{'alt_aa_long'} ne 'STOP'
                    # && $LocusInfoReader{'ref_aa_long'} ne 'STOP'
                    && $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
                    && $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON'
                    )
                {
                    # $class_query
                    #     = uc( $LocusInfoReader{'locus'} . "_" . $LocusInfoReader{'snp'} );
                    $query_locus = uc( $LocusInfoReader{'locus'} );
                    $query_snp   = uc( $LocusInfoReader{'snp'} );

                    # print "$class_query\n";
                    $loc_snp_notation = &load_coding_snp(
                        $options{snp_list},
                        $position_type,
                        $query_locus,
                        $query_snp,
                        $LocusInfoReader{'snp_genome_pos'} . "_"
                            . $LocusInfoReader{'ref'} . ">"
                            . $LocusInfoReader{'alt'},
                        $LocusInfoReader{'ref_aa_short'}
                            . $LocusInfoReader{'aa_pos'}
                            . $LocusInfoReader{'alt_aa_short'},
                        uc(       $LocusInfoReader{'ref_aa_long'}
                                . $LocusInfoReader{'aa_pos'}
                                . $LocusInfoReader{'alt_aa_long'}
                        ),
                        $LocusInfoReader{'snp'}

                    );

      # if (exists $snp_list_h{$options{snp_list}}{$query_locus}{$query_snp} )
                    if ( $loc_snp_notation ne "" )

                    {

                        # if ( $snp_count == 0 ) {

                 # print "\t$class_query\t"
                 #     . $snp_list_h{$options{snp_list}}{$query_locus}
                 #     {$query_snp}{'note'} . "\n";
                 # push(@res_a,  $snp_list_h{$options{snp_list}}{$query_locus}
                 # {$query_snp}{'note'});
                        push( @res_a, $loc_snp_notation );

                 # }
                 # else {
                 # $tab = "\t\t" x $snp_count;
                 # print "\t$tab$class_query\t"
                 #     . $snp_list_h{$options{snp_list}}{$query_locus}
                 #     {$query_snp}{'note'} . "\n";
                 # push(@res_a,  $snp_list_h{$options{snp_list}}{$query_locus}
                 # {$query_snp}{'note'});
                 # push( @res_a, $loc_snp_notation );

                        # }

                        # }

                        $class_found = 1;
                        $snp_count++;
                        $res = join( ",", sort @res_a );

                        # last;
                    }

                }

            }

        }
        close $fh;

        if ( $class_found == 0 ) {
            print "\t---";
        }
        $class_query = "";
        $file        = "";
        $class_found = 0;
        $snp_count   = 0;
        print "\t$res\n";
        $res = "";

    }

}

# sub save_circos_files {
#     my (@to_circos_data) = @_;
#     my $to_circos_str    = "";
#     my $to_circos_output = "";

#     if ( scalar(@to_circos_a) != 0 ) {
#         foreach my $i ( 0 .. scalar(@to_circos_a) - 1 ) {

#             # print "$to_circos_a[$i]\n";
#             $to_circos_str .= $to_circos_a[$i] . ";"

#         }
#         my %to_circos_h
#             = map { split /\:/, $_ } ( split /;/, $to_circos_str );
#         if (    $to_circos_h{chr} ne ""
#             and $to_circos_h{name} ne ""
#             and $to_circos_h{bcolor} ne ""
#             and $to_circos_h{color} ne ""
#             and $to_circos_h{out} ne "" )
#         {

#             # print "x" x 60 . "\n";
#             open( my $fh, '>', "$to_circos_h{out}.karyotype.txt" )
#                 or die "Cannot open $to_circos_h{out}.karyotype.txt $!";

#             print $fh
#                 "chr - $to_circos_h{'chr'} $to_circos_h{'chr'} 1 $database{options}{genome_length} $to_circos_h{'bcolor'}\n";

#             foreach my $i ( 0 .. scalar(@to_circos_data) - 1 ) {
#                 $to_circos_data[$i]
#                     =~ /^(\w+)\s+(\d+)\s+(\S+)/;

# my ($start,$end)=($2-1,$2+1);
#         print $fh "band $to_circos_h{'chr'} $1 $1_$3 $start $end $to_circos_h{'color'}\n";

#             }

#             # print "x" x 60 . "\n";

#             close $fh;
#         }
#     }
# }

<<<<<<< HEAD
# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

# sub annotate_vcf_formatted {
#     my $locInputFile = shift;
#     my %snp_per_gene_h;
#     my %gene_mutations_h;

#     # my @snp_per_gene_a;
#     my $locRefGenomeName  = "";
#     my $loc_AACI_flag = "---";
#     open( my $fh, "<", $locInputFile )
#         or croak "cannot open < $locInputFile: $!";
#     while (<$fh>) {
#         my $str = "$_";
#         $str =~ /.*($RefGenomeName)/x;
#         $locRefGenomeName = $1;
#         $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;
#         if ( length($2) == 1 && length($3) == 1 ) {
#             my %LocusInfoReader = &get_locus_info( $1, $3 );
#             if ( $locRefGenomeName ne $RefGenomeName ) {
#                 print "REFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
#                 exit;
#             }

#  # print "$LocusInfoReader{'alt_aa_long'}\t!!!!!\n" if $LocusInfoReader{'alt_aa_long'} eq "BAD_CODON";
#             if ($LocusInfoReader{'locus'} ne ""

#                 # && $LocusInfoReader{'alt_aa_long'} ne 'STOP'
#                 # && $LocusInfoReader{'ref_aa_long'} ne 'STOP'
#                 && $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
#                 && $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON'
#                 )
#             {

#                 # $snp_per_gene_h{$LocusInfoReader{'locus'}}=$LocusInfoReader{'snp'};
#                 my $locGeneLocusName = $LocusInfoReader{'locus'};
#                 my $loc_snp        = $LocusInfoReader{'snp'};
#                 my $locSnpPosInGenome    = $LocusInfoReader{'snp_genome_pos'};
#                 if ( $options{action} eq "annotation.2" ) {
#                     $snp_per_gene_h{$locGeneLocusName}{$locSnpPosInGenome}
#                         = "$loc_snp\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}";

#                 }
#                 elsif ( $options{action} eq "annotation.3" ) {

#                     $loc_AACI_flag
#                         = &calc_aa_change_info( $LocusInfoReader{'ref_aa_long'},
#                         $LocusInfoReader{'alt_aa_long'}, $LocusInfoReader{'snp_type'} );

#                     $snp_per_gene_h{$locGeneLocusName}{$locSnpPosInGenome}
#                         = "$loc_snp\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$loc_AACI_flag";
#                 }

#             }

#         }

#     }
#     close $fh;

#     # print Dumper (%snp_per_gene_h);

#     foreach my $key ( sort keys %snp_per_gene_h ) {

#         my $key2_count = 0;
#         foreach my $key2 ( sort keys %{ $snp_per_gene_h{$key} } ) {

#             # foreach my $key3 ( keys %{ $hash{$key}{$key2} } ) {
#             # $value = $hash{$key}{$key2}->{$key3};
#             $key2_count++;
#             if ( $key2_count == 1 ) {
#                 print "$key" . "\t"
#                     . $database{$key}{'product'} . "\t"
#                     . "\n\t$key2\t$snp_per_gene_h{$key}{$key2}\n";

#             }
#             else {
#                 print "\t$key2\t$snp_per_gene_h{$key}{$key2}\n";

#             }

#             # .
#             # .
#             # Do something with $value
#             # .
#             # .
#             # .
#             # }
#         }
#     }

# }

=======
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
sub make_axt_sequence {

    # my $locInputFile = shift;

    my @files = glob('*.vcf');
    my %snp_per_gene_h;
    my %gene_mutations_h;

    # my @snp_per_gene_a;
    my $locRefGenomeName  = "";
    my $loc_AACI_flag = "---";
    my $nuc_seq       = $database{ $options{locus} }{sequence};

# print "$RefGenomeName\n" . $database{$options{locus}}{sequence} . "_" . $database{$options{locus}}{sequence}  . "\n\n";
    foreach my $locInputFile (@files) {

        open( my $fh, "<", $locInputFile )
            or croak "cannot open < $locInputFile: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*($RefGenomeName)/x;
            $locRefGenomeName = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;
            if ( length($2) == 1 && length($3) == 1 ) {
                my %LocusInfoReader = &get_locus_info( $1, $3 );
                if ( $locRefGenomeName ne $RefGenomeName ) {
                    print "REFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
                    exit;
                }

 # print "$LocusInfoReader{'alt_aa_long'}\t!!!!!\n" if $LocusInfoReader{'alt_aa_long'} eq "BAD_CODON";
 # print $options{locus};
 # print $LocusInfoReader{locus} . "\t" . $options{locus} . "\n";
                if ($LocusInfoReader{'locus'} ne ""

                    # && $LocusInfoReader{'alt_aa_long'} ne 'STOP'
                    # && $LocusInfoReader{'ref_aa_long'} ne 'STOP'
                    && $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
                    && $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON'
                    && uc( $LocusInfoReader{locus} ) eq uc( $options{locus} )

                    # && $LocusInfoReader{locus} eq $options{locus}
                    )
                {

                    # $snp_per_gene_h{$LocusInfoReader{'locus'}}=$LocusInfoReader{'snp'};

                    my $locGeneLocusName = $LocusInfoReader{'locus'};
                    my $loc_snp        = $LocusInfoReader{'snp'};
                    my $locSnpPosInGenome    = $LocusInfoReader{'snp_genome_pos'};
                    $snp_per_gene_h{$locGeneLocusName}{$locSnpPosInGenome}
                        = "$loc_snp\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}";

                    $gene_mutations_h{$locGeneLocusName}{$locSnpPosInGenome}{pos}
                        = $LocusInfoReader{snp_gene_pos};
                    $gene_mutations_h{$locGeneLocusName}{$locSnpPosInGenome}{ref}
                        = $LocusInfoReader{ref};

                }
                else {
                    next;
                }

            }

        }
        close $fh;

        # print Dumper (%snp_per_gene_h);

        foreach my $key ( sort keys %snp_per_gene_h ) {

            next if uc($key) ne uc( $options{locus} );
            my $key2_count = 0;
            foreach my $key2 ( sort keys %{ $snp_per_gene_h{$key} } ) {

                # foreach my $key3 ( keys %{ $hash{$key}{$key2} } ) {
                # $value = $hash{$key}{$key2}->{$key3};
                $key2_count++;

                # if ( $key2_count == 1 ) {
                #     # print "$key" . "\t"
                #     #     . $database{$key}{'product'} . "\t"
                #     #     . "\n\t$key2\t$snp_per_gene_h{$key}{$key2}\n";
                # my $pos=$gene_mutations_h{$key}{$key2}{pos};
                # my $ref=$gene_mutations_h{$key}{$key2}{ref};
                # substr( $nuc_seq, $pos-1, 1)="$ref";
                # # print "\n$nuc_seq\n";

                # }
                # else {
                # print "\t$key2\t$snp_per_gene_h{$key}{$key2}\n";
                my $pos = $gene_mutations_h{$key}{$key2}{pos};
                my $ref = $gene_mutations_h{$key}{$key2}{ref};
                substr( $nuc_seq, $pos - 1, 1 ) = lc($ref);

                # print "\n$nuc_seq\n";

                # }

                # .
                # .
                # Do something with $value
                # .
                # .
                # .
                # }
            }

        }

    }
    print "$RefGenomeName\_"
        . length($nuc_seq)
        . "\n$nuc_seq\n"
        . $database{ $options{locus} }{sequence} . "\n\n";

}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub print_unique_snp_info {
    my $input_hash         = shift;
    my $count_f            = shift;
    my %all_snp_h          = %{$input_hash};
    my $t0                 = [gettimeofday];
    my $files_per_snp_size = 0;
    my @results_a;
    my ( $miss_stat, $syn_stat, $frame_shift_stat, $no_frame_shift_stat,
        $insertion_stat, $deletion_stat, $intergenic_stat )
        = ( 1, 1, 1, 1, 1, 1 );
    my %snp_per_gene_h;
    my %indel_h;

    #----------------------------------

    foreach my $key ( sort keys %all_snp_h ) {
        my $files_per_snp_size;
        my @filenames_a = split /\t/, $all_snp_h{$key}{'file'};
        my @alts        = split /\t/, $all_snp_h{$key}{'alt_list'};
        my @DPs         = split /\t/, $all_snp_h{$key}{'DP'};
        my $alt           = $all_snp_h{$key}{'alt'};
        my $ref           = $all_snp_h{$key}{'ref'};
        my $loc_AACI_flag = "";
        $files_per_snp_size = @filenames_a;
        my $DP_min         = min @DPs;
        my $DP_max         = max @DPs;
        my $DP_average_raw = ( sum @DPs ) / $files_per_snp_size;
        my $DP_average     = sprintf( "%.1i", $DP_average_raw );

        if (    length($alt) == 1
            and $files_per_snp_size == $count_f
            and length($ref) == 1 )
        {

            my %LocusInfoReader = &get_locus_info( $key, $alt, 1 );
            $miss_stat++ if $LocusInfoReader{snp_type} eq "missense";
            $syn_stat++  if $LocusInfoReader{snp_type} eq "synonymous";

            if (    $options{lvl} == 1
                and $LocusInfoReader{locus} ne ""
                and $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
                and $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON' )
            {

                push( @results_a,
                    "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$LocusInfoReader{'product'}\n"
                );

            }
            elsif ( $options{lvl} == 2
                and $LocusInfoReader{locus} ne ""
                and $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
                and $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON' )
            {

                push( @results_a,
                    "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$LocusInfoReader{'product'}\t$files_per_snp_size/$count_f\n"
                );
            }
            elsif ( $options{lvl} == 3
                and $LocusInfoReader{locus} ne ""
                and $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
                and $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON' )
            {

                push( @results_a,
                    "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$LocusInfoReader{tang}\t$LocusInfoReader{gd}\t$LocusInfoReader{blossum45}\t$LocusInfoReader{pam30}\t$LocusInfoReader{AACI}\t$LocusInfoReader{'product'}\n"
                );
            }
            elsif ( $options{lvl} == 4
                and $LocusInfoReader{locus} ne ""
                and $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
                and $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON' )
            {

                push( @results_a,
                    "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$LocusInfoReader{'product'}\t"
                        . "DP:"
                        . $DP_min . "-"
                        . $DP_max
                        . " (Avg: $DP_average)"
                        . "\n" );
            }
            elsif ( $options{lvl} == 5
                and $LocusInfoReader{locus} ne ""
                and $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
                and $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON' )
            {

                $snp_per_gene_h{ $LocusInfoReader{'locus'} }{ $LocusInfoReader{'snp_genome_pos'} }
                    = "\t$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$LocusInfoReader{'product'}"

            }
            elsif ( $options{lvl} == 6
                and $LocusInfoReader{locus} ne ""
                and $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
                and $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON' )
            {

                $snp_per_gene_h{ $LocusInfoReader{'locus'} }{ $LocusInfoReader{'snp_genome_pos'} }
                    = "\t$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$LocusInfoReader{tang}\t$LocusInfoReader{gd}\t$LocusInfoReader{blossum45}\t$LocusInfoReader{pam30}\t$LocusInfoReader{AACI}\t$LocusInfoReader{'product'}"

            }

            elsif ( $options{lvl} == 1
                and $LocusInfoReader{locus} eq ""
                and $options{intergenic} )
            {

                my %ig = &get_intergen_info( $key, $alt );
                my $next_end = scalar( $ig{end} ) - 1;
                my $snp_dist_next_gene
                    = $next_end - scalar( $ig{snp_genome_pos} );
                if ( $ig{locus} ne "" ) {
                    $intergenic_stat++;
                    push( @results_a,
                        "$ig{igen_genes}\t$ig{snp_genome_pos}\t$ig{snp}\t$ig{ref_codon}\/$ig{alt_codon}\t-\t-\t$ig{product} ($ig{igen_genes}) $snp_dist_next_gene\n"
                    );
                }

            }

        }

        elsif (length($alt) > 1
            && $options{indel}
            && $files_per_snp_size == $count_f
            || length($ref) > 1
            && $options{indel}
            && $files_per_snp_size == $count_f )
        {

            %indel_h = &get_indel_info( $key, $ref, $alt );
            $insertion_stat++   if $indel_h{snp_type} eq "insertion";
            $deletion_stat++    if $indel_h{snp_type} eq "deletion";
            $frame_shift_stat++ if $indel_h{frame_shift} eq "frame_shift";
            $no_frame_shift_stat++
                if $indel_h{frame_shift} eq "no_frame_shift";
            if ( $options{lvl} == 1 ) {
                push( @results_a, $indel_h{standard_report} );
            }

        }

    }

    if (   $options{lvl} == 1
        || $options{lvl} == 2
        || $options{lvl} == 3
        || $options{lvl} == 4 )
    {
        print sort { lc($a) cmp lc($b) } @results_a;

        @results_a = undef;

    }
    elsif ( $options{lvl} == 5 || $options{lvl} == 6 ) {
        foreach my $key ( sort keys %snp_per_gene_h ) {
            my $key2_count = 0;
            foreach my $key2 ( sort keys %{ $snp_per_gene_h{$key} } ) {

                # foreach my $key3 ( keys %{ $hash{$key}{$key2} } ) {
                # $value = $hash{$key}{$key2}->{$key3};
                $key2_count++;
                if ( $key2_count == 1 ) {
                    print "$key" . "\t"
                        . $database{$key}{'product'} . "\t"
                        . "\n\t$key2\t$snp_per_gene_h{$key}{$key2}\n";
                }
                else {
                    print "\t$key2\t$snp_per_gene_h{$key}{$key2}\n";
                }

            }
        }

    }

    #----------------------------------
    if ( $options{debug} ) {
        my $elapsed = tv_interval( $t0, [gettimeofday] );
        print "-" x 50
            . "\n\tSynonymous:\t"
            . $syn_stat
            . "\n\tMissense:\t"
            . $miss_stat;
        if ( $options{indel} ) {
            print "\n\tInsertions:\t"
                . $insertion_stat
                . "\n\tDeletions:\t"
                . $deletion_stat
                . "\n\tFrame shifts:\t"
                . $frame_shift_stat
                . "\n\tNo frame shifts:\t"
                . $no_frame_shift_stat;
        }
        if ( $options{intergenic} ) {
            print "\n\tIntergenic:\t" . $intergenic_stat;

        }
        print "\n\tFiles:\t$count_f" . "\n\tTime elapsed:\t$elapsed\n";

    }
}

sub open_multiple_vcfs {
    my @files = glob('*.vcf');
    my %all_snp_h;
    my $count_f = @files;

    # my %toR_h;
    my $out_header = "";
    my @toR_a;

    open( my $out_fh, '>', 'toR.csv' ) if $options{toR};

    foreach my $file (@files) {

        $out_header = "$file";

        # print "$file";

        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*($RefGenomeName)/x;
            next if $str =~ /^#/;
            my $locRefGenomeName = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

            # if ( length($2) == 1 && length($3) == 1 ) {

            # my %LocusInfoReader = &get_locus_info( $1, $3 );
            if ( $locRefGenomeName ne $RefGenomeName ) {
                print
                    "$file\tREFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
                last;
            }
            $all_snp_h{$1}{'file'}     .= "$file\t";
            $all_snp_h{$1}{'alt_list'} .= "$3\t";
            $all_snp_h{$1}{'alt'} = "$3";
            $all_snp_h{$1}{'ref'} = "$2";
            $all_snp_h{$1}{'DP'} .= "$4\t";

            if ( $options{toR} and length($3) == 1 and length($2) == 1 ) {

                $out_header .= " $1";

                # print $out_fh "$1\t$file\t$1_$2>$3\n";
            }

            # }
            # $toR_h{$file}=$out_header if $options{toR};

        }
        close $fh;

        # push( @toR_a, "$out_header\n" )
        # @toR_swap=split (/\t/,$out_header);
        push( @toR_a, "$out_header\n" );

        # if ( $options{toR}) {

        #     # $out_header .= "$file\t";
        #     print $out_fh Dumper %all_snp_h;
        # }

    }

    &print_unique_snp_info( \%all_snp_h, $count_f );
    if ( $options{toR} ) {

        foreach my $val (@toR_a) {
            print $out_fh $val;
        }
    }
}

<<<<<<< HEAD
# sub find_uniq_genes {
#     my @files       = glob('*.vcf');
#     my $count_f     = @files;       # считаем их количество
#     my $class_found = 0;
#     my $class_query = "";
#     my $query_locus = "";
#     my $query_snp   = "";
#     my @class_a;
#     my $snp_count    = 0;
#     my $tab          = "";
#     my $res          = "";
#     my $locRefGenomeName = "";
#     my %all_snp_h;
#     my @results_a;

#     my $start = time;

#     foreach my $file (@files) {

#         # print "$file";
#         open( my $fh, "<", $file )
#             or croak "cannot open < $file: $!";
#         while (<$fh>) {
#             my $str = "$_";
#             $str =~ /.*($RefGenomeName)/x;
#             next if $str =~ /^#/;
#             $locRefGenomeName = $1;
#             $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

#             # if ( length($2) == 1 && length($3) == 1 ) {

#             # my %LocusInfoReader = &get_locus_info( $1, $3 );
#             if ( $locRefGenomeName ne $RefGenomeName ) {
#                 print
#                     "$file\tREFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
#                 last;
#             }
#             $all_snp_h{$1}{'file'}     .= "$file\t";
#             $all_snp_h{$1}{'alt_list'} .= "$3\t";
#             $all_snp_h{$1}{'alt'} = "$3";
#             $all_snp_h{$1}{'DP'} .= "$4\t";

#             # }

#         }
#         close $fh;

#     }
#     foreach my $key ( sort keys %all_snp_h ) {
#         my $files_per_snp_size;
#         my @filenames_a = split /\t/, $all_snp_h{$key}{'file'};
#         my @alts        = split /\t/, $all_snp_h{$key}{'alt_list'};
#         my @DPs         = split /\t/, $all_snp_h{$key}{'DP'};
#         my $alt           = $all_snp_h{$key}{'alt'};
#         my $loc_AACI_flag = "";
#         $files_per_snp_size = @filenames_a;
#         my $DP_min         = min @DPs;
#         my $DP_max         = max @DPs;
#         my $DP_average_raw = ( sum @DPs ) / $files_per_snp_size;
#         my $DP_average     = sprintf( "%.1i", $DP_average_raw );

#         # print "$key\n\t"
#         #     . $all_snp_h{$key}{alt_list} . "\t"
#         #     . $all_snp_h{$key}{file} . "\n";

#         if ( $files_per_snp_size == $count_f ) {
#             my %LocusInfoReader = &get_locus_info( $key, $alt );

#             if ($LocusInfoReader{'locus'} ne ""

#                 # && $LocusInfoReader{'alt_aa_long'} ne 'STOP'
#                 # && $LocusInfoReader{'ref_aa_long'} ne 'STOP'
#                 && $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
#                 && $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON'
#                 )
#             {
#                 $loc_AACI_flag
#                     = &calc_aa_change_info( $LocusInfoReader{'ref_aa_long'},
#                     $LocusInfoReader{'alt_aa_long'}, $LocusInfoReader{'snp_type'} );
#                 if ( $options{action} eq "uniq" ) {
#                     push( @results_a,
#                         "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$LocusInfoReader{'product'}\t$files_per_snp_size/$count_f\n"
#                     );
#                 }
#                 elsif ( $options{action} eq "uniq.1" ) {
#                     push( @results_a,
#                         "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$LocusInfoReader{'product'}\n"
#                     );

#                 }
#                 elsif ( $options{action} eq "uniq.6" ) {
#                     push( @results_a,
#                         "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$loc_AACI_flag$LocusInfoReader{'product'}\n"
#                     );

#                 }
#                 elsif ( $options{action} eq "uniq.7" ) {
#                     push( @results_a,
#                         "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$loc_AACI_flag$LocusInfoReader{'product'}\t"
#                             . "DP:"
#                             . $DP_min . "-"
#                             . $DP_max
#                             . " (Avg: $DP_average)"
#                             . "\n" );

#                 }
#                 elsif ($options{action} eq "uniq.8"
#                     || $options{action} eq "uniq.8.1" )
#                 {
#                     ( my $AACI_res )
#                         = $loc_AACI_flag
#                         =~ /^\S+\W+\S+\W+\S+\W+\S+\s+(\W{4})/;

#                     if (    $AACI_res eq "++++"
#                         and $options{action} eq "uniq.8" )
#                     {
#                         push( @results_a,
#                             "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$loc_AACI_flag$LocusInfoReader{'product'}\n"
#                         );
#                     }
#                     elsif ( $AACI_res eq "++++"
#                         and $options{action} eq "uniq.8.1" )
#                     {
#                         push( @results_a,
#                             "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$LocusInfoReader{'product'}\n"
#                         );
#                     }

#                 }
#                 elsif ( $options{action} eq "uniq.2" ) {
#                     push( @results_a,
#                         "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'product'}\n"
#                     );

#                 }

#             }

#         }

#     }

#     print sort @results_a;
#     undef @results_a;

#     # my $duration = time - $start;
#     # print "Execution time: $duration s\n";
# }

# # ------------------------------------------------------------
# #
# #
# # ------------------------------------------------------------

# sub find_uniq_genes_formated {
#     my @files       = glob('*.vcf');
#     my $count_f     = @files;       # считаем их количество
#     my $class_found = 0;
#     my $class_query = "";
#     my $query_locus = "";
#     my $query_snp   = "";
#     my @class_a;
#     my $snp_count    = 0;
#     my $tab          = "";
#     my $res          = "";
#     my $locRefGenomeName = "";
#     my %all_snp_h;
#     my @results_a;
#     my %snp_per_gene_h;

#     my $start = time;

#     foreach my $file (@files) {

#         # print "$file";
#         open( my $fh, "<", $file )
#             or croak "cannot open < $file: $!";
#         while (<$fh>) {
#             my $str = "$_";
#             $str =~ /.*($RefGenomeName)/x;
#             $locRefGenomeName = $1;
#             $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

#             if ( length($2) == 1 && length($3) == 1 ) {

#                 # my %LocusInfoReader = &get_locus_info( $1, $3 );
#                 if ( $locRefGenomeName ne $RefGenomeName ) {
#                     print
#                         "$file\tREFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
#                     last;
#                 }
#                 $all_snp_h{$1}{'file'}     .= "$file\t";
#                 $all_snp_h{$1}{'alt_list'} .= "$3\t";
#                 $all_snp_h{$1}{'alt'} = "$3";

#             }

#         }
#         close $fh;

#     }
#     foreach my $key ( sort keys %all_snp_h ) {
#         my $files_per_snp_size;
#         my @filenames_a = split /\t/, $all_snp_h{$key}{'file'};
#         my @alts        = split /\t/, $all_snp_h{$key}{'alt_list'};
#         my $alt = $all_snp_h{$key}{'alt'};
#         $files_per_snp_size = @filenames_a;
#         if ( $files_per_snp_size == $count_f ) {
#             my %LocusInfoReader = &get_locus_info( $key, $alt );
#             if ($LocusInfoReader{'locus'} ne ""

#                 # && $LocusInfoReader{'alt_aa_long'} ne 'STOP'
#                 # && $LocusInfoReader{'ref_aa_long'} ne 'STOP'
#                 && $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
#                 && $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON'
#                 )
#             {
# # push( @results_a,
# #     "$LocusInfoReader{'locus'}\t$LocusInfoReader{'snp_genome_pos'}\t$LocusInfoReader{'snp'}\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$LocusInfoReader{'product'}\t$files_per_snp_size/$count_f\n"
# # );
#                 my $locGeneLocusName = $LocusInfoReader{'locus'};
#                 my $loc_snp        = $LocusInfoReader{'snp'};
#                 my $locSnpPosInGenome    = $LocusInfoReader{'snp_genome_pos'};
#                 my $loc_AACI_flag
#                     = &calc_aa_change_info( $LocusInfoReader{'ref_aa_long'},
#                     $LocusInfoReader{'alt_aa_long'}, $LocusInfoReader{'snp_type'} );
#                 if ( $options{action} eq "uniq.3" ) {
#                     $snp_per_gene_h{$locGeneLocusName}{$locSnpPosInGenome}
#                         = "\t$locSnpPosInGenome\t$loc_snp\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$loc_AACI_flag\t$files_per_snp_size/$count_f";
#                 }
#                 elsif ( $options{action} eq "uniq.4" ) {
#                     $snp_per_gene_h{$locGeneLocusName}{$locSnpPosInGenome}
#                         = "\t$locSnpPosInGenome\t$loc_snp\t$LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'}\t$LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'}\t$LocusInfoReader{'snp_type'}\t$loc_AACI_flag";
#                 }

#             }

#         }

#     }

#     # print sort @results_a;
#     # print Dumper(%snp_per_gene_h);
#     foreach my $key ( sort keys %snp_per_gene_h ) {
#         my $key2_count = 0;
#         foreach my $key2 ( sort keys %{ $snp_per_gene_h{$key} } ) {

#             # foreach my $key3 ( keys %{ $hash{$key}{$key2} } ) {
#             # $value = $hash{$key}{$key2}->{$key3};
#             $key2_count++;
#             if ( $key2_count == 1 ) {
#                 print "$key" . "\t"
#                     . $database{$key}{'product'} . "\t"
#                     . "\n\t$key2\t$snp_per_gene_h{$key}{$key2}\n";
#             }
#             else {
#                 print "\t$key2\t$snp_per_gene_h{$key}{$key2}\n";
#             }

#             # .
#             # .
#             # Do something with $value
#             # .
#             # .
#             # .
#             # }
#         }
#     }

#     my $duration = time - $start;

#     # print "Execution time: $duration s\n";

# }

=======
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
# # ------------------------------------------------------------
# #
# #
# ------------------------------------------------------------

sub compare_files_by_snp {
    my @files       = glob('*.vcf');
    my $count_f     = @files;       # считаем их количество
    my $class_found = 0;
    my $class_query = "";
    my $query_locus = "";
    my $query_snp   = "";
    my @class_a;
    my $snp_count    = 0;
    my $tab          = "";
    my $res          = "";
    my $locRefGenomeName = "";
    my %all_snp_h;
    my @results_a;
    my %results_h;

    # my %results_files_h;
    # my %snp2seq_h;
    my $filename = "";
    my $start    = time;

    # my $th=10;
    # $th=(@files/2);
    $options{snp_th} = 20 if ( $options{snp_th} eq "" );
    foreach my $file (@files) {
        $filename = $file;

        # print "$file";
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*($RefGenomeName)/x;
            $locRefGenomeName = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

            if ( length($2) == 1 && length($3) == 1 ) {

                # my %LocusInfoReader = &get_locus_info( $1, $3 );
                if ( $locRefGenomeName ne $RefGenomeName ) {
                    print
                        "$file\tREFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
                    last;
                }
                $all_snp_h{$1}{'pos'} = "$1";
                $all_snp_h{$1}{'file'} .= "$file,";

                # $all_snp_h{$1}{'full'}     .= "$1\t$3\t$file\n";
                $all_snp_h{$1}{'alt_list'} .= "$3,";
                $all_snp_h{$1}{'alt'} = "$3";

            }

        }

        close $fh;
    }
    foreach my $key ( sort keys %all_snp_h ) {
        my $files_per_snp_size;
        my @filenames_a = split /\,/, $all_snp_h{$key}{'file'};

        my $lc           = List::Compare->new( \@files, \@filenames_a );
        my @intersection = $lc->get_symmetric_difference();
        my @alts         = split /\,/, $all_snp_h{$key}{'alt_list'};
        my $alt          = $all_snp_h{$key}{'alt'};
        my $target_a_size;
        my @target_a;
        $files_per_snp_size = @filenames_a;
        $results_h{$key}{'nf'} = join( ',', @intersection );
        $results_h{$key}{'nf_count'} = @intersection;

        if ( $files_per_snp_size >= $options{'snp_th'} ) {
            $results_h{$key}{'count'}     = $files_per_snp_size;
            $results_h{$key}{'file_list'} = $all_snp_h{$key}{'file'};
            $results_h{$key}{'alt'}       = $alt;

    # $results_h{$key}{'full'}=$all_snp_h{$key}{'full'};
    # $results_files_h{$files_per_snp_size}{'files'}=$all_snp_h{$key}{'file'};
    # $results_files_h{$files_per_snp_size}{'pos'}=$key;
    # # $results_files_h{$files_per_snp_size}{$key}=$all_snp_h{$key}{'file'} ;
        }

    }

    foreach my $key ( sort keys %results_h ) {
        my %LocusInfoReader = &get_locus_info( $key, $results_h{$key}{'alt'} );
        if ($LocusInfoReader{'locus'} ne ""

            # && $LocusInfoReader{'alt_aa_long'} ne 'STOP'
            # && $LocusInfoReader{'ref_aa_long'} ne 'STOP'
            && $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
            && $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON'
            )
        {
            # if ( $options{color} == 1 ) {

  # print colored ['bold yellow'], "FOUND: " . $results_h{$key}{'count'} . "/"
  #     . scalar(@files);
            my @filenames_a = split /\,/, $all_snp_h{$key}{'file'};

            if ( $options{target} ne "" ) {
                if ( my ($matched) = grep $_ eq $options{target},
                    @filenames_a )
                {
                    print "\n" . "-" x 50;

                    # . "\nSNP_POS:\n\t $key";
                    print "\nLOCUS:\n\t";
                    print
                        "$LocusInfoReader{'locus'} $key $LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'} $LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'} $LocusInfoReader{'snp_type'} $LocusInfoReader{'product'} \n";

                    # print "***: $matched\n";
                    print "FOUND_IN_W_TARGET ("
                        . $results_h{$key}{'count'} . "/"
                        . scalar(@files)
                        . "):\n\t";
                    print $results_h{$key}{'file_list'} . "\n";
                    if ( $results_h{$key}{'nf_count'} > 0 ) {
                        print "NOT_FOUND_IN ("
                            . $results_h{$key}{'nf_count'} . "/"
                            . scalar(@files)
                            . "):\n\t";
                        print $results_h{$key}{'nf'};
                    }
                }
            }
            else {
                print "\n" . "-" x 50;

                # . "\nSNP_POS:\n\t $key";
                print "\nLOCUS:\n\t";
                print
                    "$LocusInfoReader{'locus'} $key $LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'} $LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'} $LocusInfoReader{'snp_type'} $LocusInfoReader{'product'} \n";
                print "FOUND_IN ("
                    . $results_h{$key}{'count'} . "/"
                    . scalar(@files)
                    . "):\n\t";
                print $results_h{$key}{'file_list'} . "\n";
                if ( $results_h{$key}{'nf_count'} > 0 ) {
                    print "NOT_FOUND_IN ("
                        . $results_h{$key}{'nf_count'} . "/"
                        . scalar(@files)
                        . "):\n\t";
                    print $results_h{$key}{'nf'};
                }
            }

# }
# else {
#     print "\n"
#         . "-" x 50
#         . "\nLOCUS:\n\t $LocusInfoReader{'locus'} $key $LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'} $LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'} $LocusInfoReader{'snp_type'} $LocusInfoReader{'product'} \n"
#         . "FOUND_IN ("
#         . $results_h{$key}{'count'} . "/"
#         . scalar(@files)
#         . "):\n\t"
#         . $results_h{$key}{'file_list'} . "\n";
#     print "NOT_FOUND_IN ("
#         . $results_h{$key}{'nf_count'} . "/"
#         . scalar(@files)
#         . "):\n\t"
#         . $results_h{$key}{'nf'}
#         if $results_h{$key}{'nf_count'} > 0;
# }

           # my @loc_files_list = split( " ", $results_h{$key}{'file_list'} );
           # foreach my $file_key ( sort @files ) {
           #     foreach my $file_key1 ( sort @loc_files_list ) {
           #         print "$file_key\t$file_key1";
           #     }
           # }

            # print $results_h{$key}{'full'};

        }
    }

# foreach my $key ( sort keys %results_files_h ) {
#     print "$key\n\t$results_files_h{$key}{'pos'}\n\t\t$results_files_h{$key}{'files'}\n";
# }
    my $duration = time - $start;

    # print "Execution time: $duration s\n";

}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub compare_list_by_snp {
    my @files       = glob('*.vcf');
    my $count_f     = @files;       # считаем их количество
    my $class_found = 0;
    my $class_query = "";
    my $query_locus = "";
    my $query_snp   = "";
    my @class_a;
    my $snp_count    = 0;
    my $tab          = "";
    my $res          = "";
    my $locRefGenomeName = "";
    my %all_snp_h;
    my @results_a;
    my %results_h;
    my @filelist_a;
    my $filelist_joined = "";
    my $files_joied     = "";
    my $dir             = getcwd;

    # my %results_files_h;
    # my %snp2seq_h;
    my $filename = "";
    my $start    = time;

    # my $th=10;
    # $th=(@files/2);
    $options{snp_th} = 20 if ( $options{snp_th} eq "" );
    if ( $options{list} eq "" ) {
        print
            "Check list is not found! Use option --list <FILE> to load filelist!\n";
        print $options{list};
        exit;
    }

    open( my $fh, "<", $options{list} )
        or croak "cannot open file!";

    while (<$fh>) {
        my ($loc_fname) = $_ =~ /(^\S+)/;

        # chomp;
        if ( -e "$dir/$loc_fname" ) {
            print uc($loc_fname) . "\t[ ";
            print "OK";
            print " ]\n";
            push( @filelist_a, $loc_fname );
        }
        else {
            print uc($loc_fname) . "\t[ ";
            print "NOT_FOUND";
            print " ]\n";
        }
    }
    close $fh;
    $filelist_joined = join( ',', @filelist_a );
    $files_joied     = join( ',', @files );

    # exit;

    foreach my $file (@files) {
        $filename = $file;

        # print "$file";
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*($RefGenomeName)/x;
            $locRefGenomeName = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

            if ( length($2) == 1 && length($3) == 1 ) {

                # my %LocusInfoReader = &get_locus_info( $1, $3 );
                if ( $locRefGenomeName ne $RefGenomeName ) {
                    print
                        "$file\tREFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
                    last;
                }
                $all_snp_h{$1}{'pos'} = "$1";
                $all_snp_h{$1}{'file'} .= "$file,";

                # $all_snp_h{$1}{'full'}     .= "$1\t$3\t$file\n";
                $all_snp_h{$1}{'alt_list'} .= "$3,";
                $all_snp_h{$1}{'alt'} = "$3";

            }

        }

        close $fh;
    }
    foreach my $key ( sort keys %all_snp_h ) {
        my $files_per_snp_size;
        my @filenames_a = split /\,/, $all_snp_h{$key}{'file'};
        my @alts        = split /\,/, $all_snp_h{$key}{'alt_list'};
        my $alt          = $all_snp_h{$key}{'alt'};
        my $lc           = List::Compare->new( \@filelist_a, \@filenames_a );
        my $lc_all       = List::Compare->new( \@filenames_a, \@files );
        my @intersection = $lc->get_intersection();
        my @complement_files = $lc_all->get_complement();
        $files_per_snp_size = @filenames_a;
        $results_h{$key}{'intersection'} = join( ',', @intersection );
        $results_h{$key}{'intersection_count'} = @intersection;
        $results_h{$key}{'not_equal'}       = join( ',', @complement_files );
        $results_h{$key}{'not_equal_count'} = @complement_files;

        if ( $files_per_snp_size >= $options{'snp_th'} ) {
            $results_h{$key}{'count'}     = $files_per_snp_size;
            $results_h{$key}{'file_list'} = $all_snp_h{$key}{'file'};
            $results_h{$key}{'alt'}       = $alt;

    # $results_h{$key}{'full'}=$all_snp_h{$key}{'full'};
    # $results_files_h{$files_per_snp_size}{'files'}=$all_snp_h{$key}{'file'};
    # $results_files_h{$files_per_snp_size}{'pos'}=$key;
    # # $results_files_h{$files_per_snp_size}{$key}=$all_snp_h{$key}{'file'} ;
        }

    }

    foreach my $key ( sort keys %results_h ) {
        my %LocusInfoReader = &get_locus_info( $key, $results_h{$key}{'alt'} );
        if ($LocusInfoReader{'locus'} ne ""

            # && $LocusInfoReader{'alt_aa_long'} ne 'STOP'
            # && $LocusInfoReader{'ref_aa_long'} ne 'STOP'
            && $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
            && $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON'
            )
        {
            # if ( $options{color} == 1 ) {
            print "\n" . "-" x 50;

            # . "\nSNP_POS:\n\t $key";
            print "\nLOCUS:\n\t";
            print
                "$LocusInfoReader{'locus'} $key $LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'} $LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'} $LocusInfoReader{'snp_type'} $LocusInfoReader{'product'} \n";

  # print colored ['bold yellow'], "FOUND: " . $results_h{$key}{'count'} . "/"
  #     . scalar(@files);
            print "FILES (" . scalar(@files) . "):\n\t";
            print "$files_joied\n";
            print "FILE_LIST (" . scalar(@filelist_a) . "):\n\t";
            print "$filelist_joined\n";
            print "FOUND ("
                . $results_h{$key}{'intersection_count'} . "/"
                . scalar(@filelist_a) . "/"
                . scalar(@files)
                . "):\n\t";
            print $results_h{$key}{'intersection'} . "\n";
            print "NOT_FOUND_IN_FILES ("
                . $results_h{$key}{'not_equal_count'} . "/"
                . scalar(@files)
                . "):\n\t";
            print $results_h{$key}{'not_equal'};

            # if ($results_h{$key}{'nf_count'} > 0 ) {

# }
# }
# else {
#     print "\n"
#         . "-" x 50
#         . "\nLOCUS:\n\t $LocusInfoReader{'locus'} $key $LocusInfoReader{'ref_codon'}/$LocusInfoReader{'alt_codon'} $LocusInfoReader{'ref_aa_short'}$LocusInfoReader{'aa_pos'}$LocusInfoReader{'alt_aa_short'} $LocusInfoReader{'snp_type'} $LocusInfoReader{'product'} \n"
#         . "FOUND_IN ("
#         . $results_h{$key}{'count'} . "/"
#         . scalar(@files)
#         . "):\n\t"
#         . $results_h{$key}{'file_list'} . "\n";
#     print "NOT_FOUND_IN ("
#         . $results_h{$key}{'nf_count'} . "/"
#         . scalar(@files)
#         . "):\n\t"
#         . $results_h{$key}{'nf'}
#         if $results_h{$key}{'nf_count'} > 0;
# }

           # my @loc_files_list = split( " ", $results_h{$key}{'file_list'} );
           # foreach my $file_key ( sort @files ) {
           #     foreach my $file_key1 ( sort @loc_files_list ) {
           #         print "$file_key\t$file_key1";
           #     }
           # }

            # print $results_h{$key}{'full'};

        }
    }

# foreach my $key ( sort keys %results_files_h ) {
#     print "$key\n\t$results_files_h{$key}{'pos'}\n\t\t$results_files_h{$key}{'files'}\n";
# }
    my $duration = time - $start;

    # print "Execution time: $duration s\n";

}

sub test4 {
    my @files       = glob('*.vcf');
    my $class_found = 0;
    my $class_query = "";
    my $query_locus = "";
    my $query_snp   = "";
    my @class_a;
    my $snp_count        = 0;
    my $tab              = "";
    my $res              = "";
    my $locRefGenomeName     = "";
    my $loc_snp_notation = "";
    my %seq_per_file_h;
    my $start = time;

    foreach my $file (@files) {

        # print "$file";
        my @res_a;
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*($RefGenomeName)/x;
            $locRefGenomeName = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;
            if ( length($2) == 1 && length($3) == 1 ) {

                # my %LocusInfoReader = &get_locus_info( $1, $3 );
                if ( $locRefGenomeName ne $RefGenomeName ) {
                    print
                        "\t\tREFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
                    last;
                }

                $seq_per_file_h{$file}{'seq'} .= $3;

                # }

            }

        }
        close $fh;

    }

    foreach my $key ( sort keys %seq_per_file_h ) {

        print ">$key\n" . $seq_per_file_h{$key}{'seq'} . "\n";
    }
    my $duration = time - $start;

    # print "Execution time: $duration s\n";
}

sub make_seq_aa {
    my @files       = glob('*.vcf');
    my $class_found = 0;
    my $class_query = "";
    my $query_locus = "";
    my $query_snp   = "";
    my @class_a;
    my $snp_count        = 0;
    my $tab              = "";
    my $res              = "";
    my $locRefGenomeName     = "";
    my $loc_snp_notation = "";
    my %seq_per_file_h;
    my %seq_ref_h;
    my %all_pos_h;
    my $start          = time;
    my $char_count     = 0;
    my $check_ref_used = 0;
    my $pb             = ".";

<<<<<<< HEAD
    print "Reading files:\n"  if $options{o} ne "";
=======
    print "Reading files:\n" if $options{o} ne "";
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
    open( my $fh_out, ">", $options{o} ) if $options{o} ne "";

    foreach my $file (@files) {

<<<<<<< HEAD
        print "$file\n"  if $options{o} ne "";
        my @res_a;
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*($RefGenomeName)/x;
            $locRefGenomeName = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;
            if ( length($2) == 1 && length($3) == 1 ) {

                my %LocusInfoReader = &get_locus_info( $1, $3 );
                if ( $locRefGenomeName ne $RefGenomeName ) {
                    print
                        "\t\tREFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
                    last;
                }
                if ($LocusInfoReader{'locus'} ne ""

                    # && $LocusInfoReader{'alt_aa_long'} ne 'STOP'
                    # && $LocusInfoReader{'ref_aa_long'} ne 'STOP'
                    && $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
                    && $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON'
                    && $4 > 10
                    )
                {
                    $char_count++;
                    $seq_per_file_h{$file}{'seq'} .= $LocusInfoReader{"alt_aa_short"};

       # $seq_per_file_h{$file}{$LocusInfoReader{"snp_genome_pos"}}= $LocusInfoReader{"alt_aa_short"};
                    $seq_per_file_h{$file}{'aa_pos'}
                        .= $char_count . ":"
                        . $LocusInfoReader{"locus"} . ":"
                        . $LocusInfoReader{"snp_genome_pos"} . ":"
                        . $LocusInfoReader{"ref_aa_short"}
                        . $LocusInfoReader{"aa_pos"}
                        . $LocusInfoReader{"alt_aa_short"} . ":"
                        . $4 . "|";
                    $seq_per_file_h{$file}{ $LocusInfoReader{"snp_genome_pos"} }
                        = $LocusInfoReader{"alt_aa_short"};

                    # if ($check_ref_used ==0){
                    $seq_ref_h{ $LocusInfoReader{"snp_genome_pos"} }
                        = $LocusInfoReader{"ref_aa_short"};
                    $all_pos_h{ $LocusInfoReader{"snp_genome_pos"} }
                        = $LocusInfoReader{"snp_genome_pos"};

                    # }

                }

            }

        }
        close $fh;
        $char_count = 0;
        local $| = 1;

        # my $pb =".";
        # print "\b$pb"; if $options{debug};
        # $pb .= ".";

    }
    print "done!\n"  if $options{o} ne "";

    # system("clear");
    # print "-" x 50 . "\n";

    foreach my $file (@files) {
        my $res_seq = "";
        foreach my $key ( sort { $a <=> $b } keys %all_pos_h ) {
            if ( exists $seq_per_file_h{$file}{$key} ) {
                $res_seq .= $seq_per_file_h{$file}{$key};
            }
            else {
                $res_seq .= $seq_ref_h{$key};
            }

        }

        print ">" . basename( $file, ".vcf" ) . "\n$res_seq\n"
            if $options{o} eq "";
        print $fh_out ">" . basename( $file, ".vcf" ) . "\n$res_seq\n"
            if $options{o} ne "";

    }

    # ---------------------------------------------------------
    # my $ref_seq = "";
    # foreach my $key ( sort { $a <=> $b } keys %all_pos_h ) {
    #     $ref_seq .= $seq_ref_h{$key};
    # }
    # print ">$RefGenomeName\t" . length($ref_seq) . "\n$ref_seq\n";
    # ----------------------------------------------------------

    # print Dumper(sort {$a <=> $b} %seq_ref_h);
    # if ($check_ref_used == 0){
    # print ">$RefGenomeName\n";
    # foreach my $key ( sort { $a <=> $b} keys  %seq_ref_h ) {
    #     print $seq_ref_h{$key};

    # }
    # print "\n";
    # $check_ref_used=1;
    # }
    #-----------------------------------------------
    # foreach my $key ( sort keys %seq_per_file_h ) {

    #     print ">$key\n" . $seq_per_file_h{$key}{'seq'} . "\n";
    # }

    if ( $options{action} eq "t5.1" ) {
        print "-" x 50 . "\n";
        foreach my $key ( sort keys %seq_per_file_h ) {

            print ">$key\n" . $seq_per_file_h{$key}{'aa_pos'} . "\n";

        }
        my $duration = time - $start;
    }

    print "Work was finished. The $options{o} file was saved.\n"
        if $options{o} ne "";

    # print "Execution time: $duration s\n";

}


sub make_seq_nc {
    my @files       = glob('*.vcf');
    my $class_found = 0;
    my $class_query = "";
    my $query_locus = "";
    my $query_snp   = "";
    my @class_a;
    my $snp_count        = 0;
    my $tab              = "";
    my $res              = "";
    my $locRefGenomeName     = "";
    my $loc_snp_notation = "";
    my %seq_per_file_h;
    my %seq_ref_h;
    my %all_pos_h;
    my $start          = time;
    my $char_count     = 0;
    my $check_ref_used = 0;
    my $pb             = ".";

    print "Reading files:\n"  if $options{o} ne "";
    open( my $fh_out, ">", $options{o} ) if $options{o} ne "";

    foreach my $file (@files) {

        print "$file\n"  if $options{o} ne "";
=======
        print "$file\n" if $options{o} ne "";
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
        my @res_a;
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*($RefGenomeName)/x;
            $locRefGenomeName = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;
            if ( length($2) == 1 && length($3) == 1 ) {

                my %LocusInfoReader = &get_locus_info( $1, $3 );
                if ( $locRefGenomeName ne $RefGenomeName ) {
                    print
                        "\t\tREFERENCE GENOME: $RefGenomeName NOT FOUND!!!\n";
                    last;
                }
                if ($LocusInfoReader{'locus'} ne ""

                    # && $LocusInfoReader{'alt_aa_long'} ne 'STOP'
                    # && $LocusInfoReader{'ref_aa_long'} ne 'STOP'
                    && $LocusInfoReader{'alt_aa_long'} ne 'BAD_CODON'
                    && $LocusInfoReader{'ref_aa_long'} ne 'BAD_CODON'
                    && $4 > 10
                    )
                {
                    $char_count++;
                    $seq_per_file_h{$file}{'seq'} .= $LocusInfoReader{"alt"};

       # $seq_per_file_h{$file}{$LocusInfoReader{"snp_genome_pos"}}= $LocusInfoReader{"alt_aa_short"};
                    $seq_per_file_h{$file}{'aa_pos'}
                        .= $char_count . ":"
                        . $LocusInfoReader{"locus"} . ":"
                        . $LocusInfoReader{"snp_genome_pos"} . ":"
                        . $LocusInfoReader{"ref_aa_short"}
                        . $LocusInfoReader{"aa_pos"}
                        . $LocusInfoReader{"alt_aa_short"} . ":"
                        . $4 . "|";
                    $seq_per_file_h{$file}{ $LocusInfoReader{"snp_genome_pos"} }
                        = $LocusInfoReader{"alt"};

                    # if ($check_ref_used ==0){
                    $seq_ref_h{ $LocusInfoReader{"snp_genome_pos"} }
                        = $LocusInfoReader{"ref"};
                    $all_pos_h{ $LocusInfoReader{"snp_genome_pos"} }
                        = $LocusInfoReader{"snp_genome_pos"};

                    # }

                }

            }

        }
        close $fh;
        $char_count = 0;
        local $| = 1;

        # my $pb =".";
        # print "\b$pb"; if $options{debug};
        # $pb .= ".";

    }
<<<<<<< HEAD
    print "done!\n"  if $options{o} ne "";
=======
    print "done!\n" if $options{o} ne "";
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e

    # system("clear");
    # print "-" x 50 . "\n";

    foreach my $file (@files) {
        my $res_seq = "";
        foreach my $key ( sort { $a <=> $b } keys %all_pos_h ) {
            if ( exists $seq_per_file_h{$file}{$key} ) {
                $res_seq .= $seq_per_file_h{$file}{$key};
            }
            else {
                $res_seq .= $seq_ref_h{$key};
            }

        }

        print ">" . basename( $file, ".vcf" ) . "\n$res_seq\n"
            if $options{o} eq "";
        print $fh_out ">" . basename( $file, ".vcf" ) . "\n$res_seq\n"
            if $options{o} ne "";

    }

    # ---------------------------------------------------------
    # my $ref_seq = "";
    # foreach my $key ( sort { $a <=> $b } keys %all_pos_h ) {
    #     $ref_seq .= $seq_ref_h{$key};
    # }
    # print ">$RefGenomeName\t" . length($ref_seq) . "\n$ref_seq\n";
    # ----------------------------------------------------------

    # print Dumper(sort {$a <=> $b} %seq_ref_h);
    # if ($check_ref_used == 0){
    # print ">$RefGenomeName\n";
    # foreach my $key ( sort { $a <=> $b} keys  %seq_ref_h ) {
    #     print $seq_ref_h{$key};

    # }
    # print "\n";
    # $check_ref_used=1;
    # }
    #-----------------------------------------------
    # foreach my $key ( sort keys %seq_per_file_h ) {

    #     print ">$key\n" . $seq_per_file_h{$key}{'seq'} . "\n";
    # }

    if ( $options{action} eq "t5.1" ) {
        print "-" x 50 . "\n";
        foreach my $key ( sort keys %seq_per_file_h ) {

            print ">$key\n" . $seq_per_file_h{$key}{'aa_pos'} . "\n";

        }
        my $duration = time - $start;
    }

    print "Work was finished. The $options{o} file was saved.\n"
<<<<<<< HEAD
=======
        if $options{o} ne "";

    # print "Execution time: $duration s\n";

}

sub make_seq_nc {
    my @files       = glob('*.vcf');
    my $class_found = 0;
    my $class_query = "";
    my $query_locus = "";
    my $query_snp   = "";
    my @class_a;
    my $snp_count        = 0;
    my $tab              = "";
    my $res              = "";
    my $loc_ref_name     = "";
    my $loc_snp_notation = "";
    my %seq_per_file_h;
    my %seq_ref_h;
    my %all_pos_h;
    my $start          = time;
    my $char_count     = 0;
    my $check_ref_used = 0;
    my $pb             = ".";
    my $file_count     = 0;
    print "Total files: " . scalar(@files) . "\n" if $options{o} ne "";
    open( my $fh_out, ">", $options{o} ) if $options{o} ne "";

    foreach my $file (@files) {
        $file_count++;
        print "$file_count) $file\n" if $options{o} ne "";
        my @res_a;
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*($ref_genome_name)/x;
            $loc_ref_name = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;
            if ( length($2) == 1 && length($3) == 1 ) {

                my %tmp = &get_locus_info( $1, $3 );
                if ( $loc_ref_name ne $ref_genome_name ) {
                    print
                        "\t\tREFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
                    last;
                }
                if ($tmp{'locus'} ne ""

                    # && $tmp{'alt_aa_long'} ne 'STOP'
                    # && $tmp{'ref_aa_long'} ne 'STOP'
                    && $tmp{'alt_aa_long'} ne 'BAD_CODON'
                    && $tmp{'ref_aa_long'} ne 'BAD_CODON'
                    && $4 > 10
                    )
                {
                    $char_count++;
                    if ( $options{codon} ) {
                        $seq_per_file_h{$file}{'seq'} .= $tmp{"alt_codon"};
                        $seq_per_file_h{$file}{'aa_pos'}
                            .= $char_count . ":"
                            . $tmp{"locus"} . ":"
                            . $tmp{"snp_genome_pos"} . ":"
                            . $tmp{"ref_aa_short"}
                            . $tmp{"aa_pos"}
                            . $tmp{"alt_aa_short"} . ":"
                            . $4 . "|";
                        $seq_per_file_h{$file}{ $tmp{"snp_genome_pos"} }
                            = $tmp{"alt_codon"};

                        # if ($check_ref_used ==0){
                        $seq_ref_h{ $tmp{"snp_genome_pos"} }
                            = $tmp{"ref_codon"};
                        $all_pos_h{ $tmp{"snp_genome_pos"} }
                            = $tmp{"snp_genome_pos"};

                    }
                    else {
                        $seq_per_file_h{$file}{'seq'} .= $tmp{"alt"};
                        $seq_per_file_h{$file}{'aa_pos'}
                            .= $char_count . ":"
                            . $tmp{"locus"} . ":"
                            . $tmp{"snp_genome_pos"} . ":"
                            . $tmp{"ref_aa_short"}
                            . $tmp{"aa_pos"}
                            . $tmp{"alt_aa_short"} . ":"
                            . $4 . "|";
                        $seq_per_file_h{$file}{ $tmp{"snp_genome_pos"} }
                            = $tmp{"alt"};

                        # if ($check_ref_used ==0){
                        $seq_ref_h{ $tmp{"snp_genome_pos"} }
                            = $tmp{"ref"};
                        $all_pos_h{ $tmp{"snp_genome_pos"} }
                            = $tmp{"snp_genome_pos"};

                    }

       # $seq_per_file_h{$file}{$tmp{"snp_genome_pos"}}= $tmp{"alt_aa_short"};

                    # }

                }

            }

        }
        close $fh;
        $char_count = 0;
        local $| = 1;

        # my $pb =".";
        # print "\b$pb"; if $options{debug};
        # $pb .= ".";

    }
    print "done!\n" if $options{o} ne "";

    # system("clear");
    # print "-" x 50 . "\n";

    foreach my $file (@files) {
        my $res_seq = "";
        foreach my $key ( sort { $a <=> $b } keys %all_pos_h ) {
            if ( exists $seq_per_file_h{$file}{$key} ) {
                $res_seq .= $seq_per_file_h{$file}{$key};
            }
            else {
                $res_seq .= $seq_ref_h{$key};
            }

        }

        print ">" . basename( $file, ".vcf" ) . "\n$res_seq\n"
            if $options{o} eq "";
        print $fh_out ">" . basename( $file, ".vcf" ) . "\n$res_seq\n"
            if $options{o} ne "";

    }

    # ---------------------------------------------------------
    # my $ref_seq = "";
    # foreach my $key ( sort { $a <=> $b } keys %all_pos_h ) {
    #     $ref_seq .= $seq_ref_h{$key};
    # }
    # print ">$ref_genome_name\t" . length($ref_seq) . "\n$ref_seq\n";
    # ----------------------------------------------------------

    # print Dumper(sort {$a <=> $b} %seq_ref_h);
    # if ($check_ref_used == 0){
    # print ">$ref_genome_name\n";
    # foreach my $key ( sort { $a <=> $b} keys  %seq_ref_h ) {
    #     print $seq_ref_h{$key};

    # }
    # print "\n";
    # $check_ref_used=1;
    # }
    #-----------------------------------------------
    # foreach my $key ( sort keys %seq_per_file_h ) {

    #     print ">$key\n" . $seq_per_file_h{$key}{'seq'} . "\n";
    # }

    if ( $options{action} eq "t5.1" ) {
        print "-" x 50 . "\n";
        foreach my $key ( sort keys %seq_per_file_h ) {

            print ">$key\n" . $seq_per_file_h{$key}{'aa_pos'} . "\n";

        }
        my $duration = time - $start;
    }

    print "Work was finished. The $options{o} file was saved.\n"
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e
        if $options{o} ne "";

    # print "Execution time: $duration s\n";

}



# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub get_locus_info {

    # my $locSnpPosInGenome    = shift;
    # my $loc_alt        = shift;
    my ( $locSnpPosInGenome, $loc_alt, $enable_calcs ) = @_;
    my $locGeneLocusName = "";
    my %LocusInfoProvider;
    my $snpPosInGene;
    my $typeOfSNP    = "";
    my $gene_length = 0;
    my ( $tang, $gd, $blossum45, $pam30 );
    my $loc_AACI_flag = "----";

    foreach my $key ( sort keys %database ) {

        if (   $locSnpPosInGenome >= scalar $database{$key}{'start'}
            && $locSnpPosInGenome <= scalar $database{$key}{'end'} )
        {
            $locGeneLocusName = $key;
            last;

        }

    }

    if ( $locGeneLocusName ne "" ) {
        my $locStartOfGene        = $database{$locGeneLocusName}{"start"};
        my $locEndOfGene          = $database{$locGeneLocusName}{"end"};
        my $locStrandDirection       = $database{$locGeneLocusName}{"strand"};
        my $locTranslationTable = $database{$locGeneLocusName}{"transl_table"};
        my $loc_product      = $database{$locGeneLocusName}{"product"};
        my $loc_note         = $database{$locGeneLocusName}{"note"};
        my $loc_nuc_seq      = $database{$locGeneLocusName}{"sequence"};

        if ( $locStrandDirection == 1 ) {
            $snpPosInGene = ( ( $locSnpPosInGenome - $locStartOfGene ) + 1 );
        }
        elsif ( $locStrandDirection == -1 ) {
            $snpPosInGene = ( ( $locEndOfGene - $locSnpPosInGenome ) + 1 );
            $loc_alt =~ tr/ACGT/TGCA/;
        }
        if ( $locEndOfGene > $locStartOfGene ) {
            $gene_length = ( $locEndOfGene - $locStartOfGene ) + 1;
        }
        elsif ( $locStartOfGene > $locEndOfGene ) {
            $gene_length = ( $locStartOfGene - $locEndOfGene ) + 1;
        }
        my $loc_ref = substr( $loc_nuc_seq, $snpPosInGene - 1, 1 );
        my $codon_nbr = int( ( $snpPosInGene - 1 ) / 3 ) + 1;
        my @codons        = unpack( "(A3)*", $loc_nuc_seq );
        my $res_codon     = lc( $codons[ $codon_nbr - 1 ] );
        my $loc_codon_pos = ( ( $codon_nbr * 3 ) - $snpPosInGene );
        if ( $loc_codon_pos == 2 ) {
            $loc_codon_pos = 0;
        }
        elsif ( $loc_codon_pos == 0 ) {
            $loc_codon_pos = 2;
        }

        substr( $res_codon, $loc_codon_pos, 1 )
            = uc( substr( $res_codon, $loc_codon_pos, 1 ) );
        my $alt_res_codon = $res_codon;
        substr $alt_res_codon, $loc_codon_pos, 1, $loc_alt;
        my $ref_aa = codon2aa($res_codon);
        my $alt_aa = codon2aa($alt_res_codon);
        if (    $ref_aa eq $alt_aa
            and $alt_aa ne "X"
            and $alt_aa ne "BAD_CODON"
            and $ref_aa ne "BAD_CODON" )
        {
            $typeOfSNP = "synonymous";
        }
        elsif ( $ref_aa ne $alt_aa
            and $alt_aa ne "X"
            and $alt_aa ne "BAD_CODON"
            and $ref_aa ne "BAD_CODON" )
        {
            $typeOfSNP = "missense";
        }
        elsif ( $ref_aa ne $alt_aa
            and $alt_aa eq "X"
            and $alt_aa ne "BAD_CODON"
            and $ref_aa ne "BAD_CODON" )
        {
            $typeOfSNP = "nonsense";

        }

        if ( $enable_calcs eq 1 and $typeOfSNP eq "missense" ) {
            $tang = calcutalte_tang_index( $ref_aa, $alt_aa );
            $gd = calculate_grantham_matrix( $ref_aa, $alt_aa );
            $blossum45 = calculate_blosum45_matrix( $ref_aa, $alt_aa );
            $pam30 = calculate_pam30_matrix( $ref_aa, $alt_aa );
            if ( $tang < TANG_TH ) {

                substr $loc_AACI_flag, 0, 1, "+";
            }

            if ( $gd > GD_TH ) {

                substr $loc_AACI_flag, 1, 1, "+";
            }
            if ( $blossum45 < BLOSSUM45_TH ) {

                substr $loc_AACI_flag, 2, 1, "+";
            }
            if ( $pam30 < PAM30_TH ) {

                substr $loc_AACI_flag, 3, 1, "+";
            }

        }
        elsif ( $enable_calcs eq 1 and $typeOfSNP eq "synonymous" ) {
            $tang          = "-";
            $gd            = "-";
            $blossum45     = "-";
            $pam30         = "-";
            $loc_AACI_flag = "----";
        }
        elsif ( $enable_calcs eq 1 and $typeOfSNP eq "nonsense" ) {
            $tang          = "-";
            $gd            = "-";
            $blossum45     = "-";
            $pam30         = "-";
            $loc_AACI_flag = "----";
        }

        %LocusInfoProvider = (
            "locus"          => $locGeneLocusName,
            "snp_genome_pos" => $locSnpPosInGenome,
            "start"          => $locStartOfGene,
            "end"            => $locEndOfGene,
            "strand"         => $locStrandDirection,
            "product"        => $loc_product,
            "sequence"       => $loc_nuc_seq,
            "snp_gene_pos"   => $snpPosInGene,
            "alt"            => $loc_alt,
            "snp"            => "$snpPosInGene$loc_ref>$loc_alt",
            "ref"            => $loc_ref,
            "aa_pos"         => $codon_nbr,
            "ref_codon"      => $res_codon,
            "snp_type"       => $typeOfSNP,
            "ref_aa_long"    => $ref_aa,
            "alt_aa_long"    => $alt_aa,
            "ref_aa_short"   => &aa_decode($ref_aa),
            "alt_aa_short"   => &aa_decode($alt_aa),
            "alt_codon"      => $alt_res_codon,
            "product"        => $loc_product,
            "note"           => $loc_note,
            "gene_length"    => $gene_length,
            "protein_length" => scalar(@codons) - 1,
            "tang"           => $tang,
            "gd"             => $gd,
            "blossum45"      => $blossum45,
            "pam30"          => $pam30,
            "AACI"           => $loc_AACI_flag,

        );
        return %LocusInfoProvider;
        $loc_AACI_flag = "----";
    }
    else {
        return;
    }
}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub get_locus_name {
    my $locSnpPosInGenome    = shift;
    my $locGeneLocusName = "";

    foreach my $key ( keys %database ) {
        if (   $locSnpPosInGenome >= scalar $database{$key}{'start'}
            && $locSnpPosInGenome <= scalar $database{$key}{'end'} )
        {
            $locGeneLocusName = $key;
            last;
        }

    }

    if ( $locGeneLocusName ne "" ) {
        return $locGeneLocusName;
    }
}

sub check_position {
    my $locSnpPosInGenome = shift;
    my $type_pos    = "";

    foreach my $key ( keys %database ) {
        if (   $locSnpPosInGenome >= scalar $database{$key}{'start'}
            && $locSnpPosInGenome <= scalar $database{$key}{'end'} )
        {
            $type_pos = "coding";
            last;
        }
        else {
            $type_pos = "intergenic";
        }

    }

    if ( $type_pos ne "" ) {
        return $type_pos;
    }
}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub calcutalte_tv_ti {
    my $nc_ref     = uc(shift);
    my $nc_alt     = uc(shift);
    my (%tvti_dic) = (
        "AG" => "ti",
        "TC" => "ti",
        "AC" => "tv",
        "AT" => "tv",
        "GC" => "tv",
        "GT" => "tv",
    );
    my $nc_query = $nc_ref . $nc_alt;

    if ( $tvti_dic{$nc_query} ne "" ) {
        return $tvti_dic{$nc_query};

    }
    else {
        return scalar $tvti_dic{ reverse($nc_query) };

    }

}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub calcutalte_tang_index {

    my $aa_ref = &aa_decode(shift);
    my $aa_alt = &aa_decode(shift);

    # my $aa_query = $aa_ref . $aa_alt;

    if ( $tang_dic{ $aa_ref . $aa_alt } ne "" ) {
        return $tang_dic{ $aa_ref . $aa_alt };

    }
    else {
        return scalar $tang_dic{ $aa_alt . $aa_ref };

    }

}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub calcutalte_deltaH {

    my $aa_ref = &aa_decode(shift);
    my $aa_alt = &aa_decode(shift);

    if ( $deltaH_dic{ $aa_ref . $aa_alt } ne "" ) {
        return $deltaH_dic{ $aa_ref . $aa_alt };

    }
    else {
        return scalar $deltaH_dic{ $aa_alt . $aa_ref };

    }

}

# sub calcutalte_sneath {

#     # my ($aa_input) = shift =~ /(?<ref>\w+)\t(?<alt>\w+)/;
#     my $aa_ref = &aa_decode(shift);
#     my $aa_alt = &aa_decode(shift);
#     my $query  = "$aa_ref$aa_alt";

#     if ( $sneath_dic{$query} ne "" ) {
#         return $sneath_dic{$query};

#     }
#     else {
#         return $sneath_dic{ reverse($query) };

#     }

# }

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub calculate_grantham_matrix {
    my $aa_ref = &aa_decode(shift);
    my $aa_alt = &aa_decode(shift);

    # print @AAs;

    # print values $LocusInfoReader{matrix}{R};
    if ( $aa_ref ne "" && $aa_alt ne "" ) {

        # return $GD_dic->{$aa_ref}->{$aa_alt};
        return $GD_dic{$aa_ref}{$aa_alt};

        # exit;

        # print $LocusInfoReader->{GD}->{$aa_ref}->{$aa_alt} . "!!!!!!\n";
        # exit;
    }
}

sub calculate_blosum45_matrix {
    my $aa_ref = &aa_decode(shift);
    my $aa_alt = &aa_decode(shift);

    # print @AAs;

    # print values $LocusInfoReader{matrix}{R};
    if ( $aa_ref ne "" && $aa_alt ne "" ) {

        # return $GD_dic->{$aa_ref}->{$aa_alt};
        return $B45_dic{$aa_ref}{$aa_alt};

        # exit;

        # print $LocusInfoReader->{GD}->{$aa_ref}->{$aa_alt} . "!!!!!!\n";
        # exit;
    }
}

sub calculate_pam30_matrix {
    my $aa_ref = &aa_decode(shift);
    my $aa_alt = &aa_decode(shift);

    if ( $aa_ref ne "" && $aa_alt ne "" ) {

        return $PAM30_dic{$aa_ref}{$aa_alt};

    }
}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub calc_aa_change_info {

    # my $ref_aa_long   = shift;
    # my $alt_aa_long   = shift;
    # my $type          = shift;
    my ( $ref_aa_long, $alt_aa_long, $type ) = @_;
    my $loc_AACI_flag = "----";
    my $loc_AACI;
    my $AACI;
    my $loc_tang_th;
    if ( $options{tang_th} ne "" ) {
        $loc_tang_th = $options{tang_th};
    }
    else {
        $loc_tang_th = 0.5;
    }
    if ( $type eq "missense" ) {
        my $loc_tang_index
            = calcutalte_tang_index( $ref_aa_long, $alt_aa_long );

        # my $loc_deltaH = calcutalte_deltaH( $ref_aa_long, $alt_aa_long );
        my $loc_GD = calculate_grantham_matrix( $ref_aa_long, $alt_aa_long );
        my $blosum = calculate_blosum45_matrix( $ref_aa_long, $alt_aa_long );
        my $pam = calculate_pam30_matrix( $ref_aa_long, $alt_aa_long );
        if ( $loc_tang_index < $loc_tang_th ) {
            $loc_AACI++;
            substr $loc_AACI_flag, 0, 1, "+";
        }

        # if ( $loc_deltaH > 1.22 ) {
        #     $loc_AACI++;
        #     substr $loc_AACI_flag, 1, 1, "+";
        # }
        if ( $loc_GD > 100 ) {
            $loc_AACI++;
            substr $loc_AACI_flag, 1, 1, "+";
        }
        if ( $blosum < 0 ) {
            $loc_AACI++;
            substr $loc_AACI_flag, 2, 1, "+";
        }
        if ( $pam < 0 ) {
            $loc_AACI++;
            substr $loc_AACI_flag, 3, 1, "+";
        }

        # return $loc_AACI_flag;

        $AACI = "$loc_tang_index\t$loc_GD\t$blosum\t$pam\t$loc_AACI_flag\t";
        return $AACI;
    }
    else {
        # return "---";

        $AACI = "-\t-\t-\t-\t$loc_AACI_flag\t";
        return $AACI;
    }

}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub get_snp_notation {
    my $mut = shift;
    my %result;
    $mut =~ /^(\d+)(\D)>(\D)/x;
    if ( $1 ne "" && $2 ne "" && $3 ne "" ) {

        # return "nc_gene";
        $result{type} = 'nc_gene';
        $result{pos}  = $1;
        $result{ref}  = $2;
        $result{alt}  = $3;
        return %result;
    }
    else {
        $mut =~ /^(\D)(\d+)(\D)/x;
        if ( $1 ne "" && $2 ne "" && $3 ne "" ) {

            # return "aa";
            $result{type} = 'aa';
            $result{pos}  = $2;
            $result{ref}  = $1;
            $result{alt}  = $3;
            return %result;
        }
        else {

            $mut =~ /^(\d+)_(\w)>(\w)/x;
            if ( $1 ne "" && $2 ne "" && $3 ne "" ) {

                # return "nc_genome";
                $result{type} = 'nc_genome';
                $result{pos}  = $1;
                $result{ref}  = $2;
                $result{alt}  = $3;
                return %result;
            }
            else {
                $mut =~ /^^(\w{3})(\d+)(\w{3})/x;
                if ( $1 ne "" && $2 ne "" && $3 ne "" ) {

                    # return "aa_long";
                    $result{type} = 'aa_long';
                    $result{pos}  = $2;
                    $result{ref}  = $1;
                    $result{alt}  = $3;
                    return %result;
                }
                else {

                    $result{type} = 'bad';
                    return %result;
                }

            }

        }
    }
}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub load_coding_snp {

    # my $snp_class_name = shift;
    # my $position_type  = shift;
    # my $query_locus    = shift;
    # my $query_snp      = shift;
    # my $position       = shift;
    # my $aa             = shift;
    # my $aa_long        = uc(shift);
    # my $snp            = shift;

    my ( $snp_class_name, $position_type, $query_locus, $query_snp,
        $position, $aa, $aa_long, $snp )
        = @_;

    my ($position_extracted) = $position =~ /^(\d+)/;
    my ($aa_reformat)        = $aa_long =~ /^(\D{3})(\d+)(\D{3})/;
    my ($aa_codon_pos)       = $aa_long =~ /^\D{3}(\d+)\D{3}/;
    my $aa_long_genome_pos   = "$position_extracted\_$1>$3";
    my $pos_check            = "POS$position_extracted";
    my $codon_check          = "CODON$aa_codon_pos";

    # print "$aa_long_genome_pos\n";

    # if ($position eq "489935_G>C"){print "!!!!!\n"};
    # my $snp_class_type = $snp_list_type_h{$options{snp_list}};
    if ( exists $snp_list_h{$snp_class_name}{$query_locus}{$query_snp} ) {
        return $snp_list_h{$snp_class_name}{$query_locus}{$query_snp}{'note'}
            . "[$query_locus:$query_snp]";
    }
    elsif ( exists $snp_list_h{$snp_class_name}{$position} ) {
        return $snp_list_h{$snp_class_name}{$position}{'note'};
    }
    elsif ( exists $snp_list_h{$snp_class_name}{$pos_check} ) {
        return $snp_list_h{$snp_class_name}{$pos_check}{'note'}
            . "\ ($query_locus\_$position_extracted\_$snp\_$aa_long)";
    }
    elsif ( exists $snp_list_h{$snp_class_name}{$aa_long_genome_pos} ) {
        return $snp_list_h{$snp_class_name}{$aa_long_genome_pos}{'note'};
    }
    if ( exists $snp_list_h{$snp_class_name}{$query_locus}{$aa} ) {
        return $snp_list_h{$snp_class_name}{$query_locus}{$aa}{'note'};
    }
    if ( exists $snp_list_h{$snp_class_name}{$query_locus}{$codon_check} ) {
        return $snp_list_h{$snp_class_name}{$query_locus}{$codon_check}
            {'note'} . "\ ($query_locus\_$position\_$snp\_$aa_long)";
    }
    if ( exists $snp_list_h{$snp_class_name}{$query_locus}{$aa_long} ) {
        return $snp_list_h{$snp_class_name}{$query_locus}{$aa_long}{'note'};
    }

  # my $snp_class_type = $snp_list_type_h{$options{snp_list}};
  # if ( $snp_class_type eq "nc_gene"
  #     && exists $snp_list_h{$snp_class_name}{$query_locus}{$query_snp} )
  # {
  #     return $snp_list_h{$snp_class_name}{$query_locus}{$query_snp}{'note'};
  # }
  # elsif ( $snp_class_type eq "nc_genome"
  #     && exists $snp_list_h{$snp_class_name}{$position} )
  # {
  #     return $snp_list_h{$snp_class_name}{$position}{'note'};
  # }
  # if ( $snp_class_type eq "aa"
  #     && exists $snp_list_h{$snp_class_name}{$query_locus}{$aa} )
  # {
  #     return $snp_list_h{$snp_class_name}{$query_locus}{$aa}{'note'};
  # }
  # if ( $snp_class_type eq "aa_long"
  #     && exists $snp_list_h{$snp_class_name}{$query_locus}{$aa_long} )
  # {
  #     return $snp_list_h{$snp_class_name}{$query_locus}{$aa_long}{'note'};
  # }
}

sub show_help {
    print "Version: $version by Viacheslav Sinkov (vsinkov_at_gmail.com)\n\n";
    pod2usage( "verbose" => 99, -utf8 );
}

sub get_information {
    my ($locus_name) = @_;
    my $gene_length = 0;
    my $snp = $options{snp} if $options{snp} ne "";
    if ( $locus_name ne "" && $database{$locus_name}{"end"} ne "" ) {
        my $nuc_seq = $database{$locus_name}{"sequence"};

        my @codons = unpack( "(A3)*", $nuc_seq );
        my @aminoAcids    = map { codon2aa $_} @codons;
        my @aa_one_letter = map { aa_decode $_} @aminoAcids;
        my $aa_seq = join '', @aa_one_letter;
        my @aaMw  = map { aa_mw_calc $_} @aminoAcids;
        my $start = $database{$locus_name}{"start"};
        my $end   = $database{$locus_name}{"end"};

        my $rna = $nuc_seq;
        $rna =~ tr/T/U/;
        if ( $end > $start ) {
            $gene_length = ( $end - $start ) + 1;
        }
        elsif ( $start > $end ) {
            $gene_length = ( $start - $end ) + 1;
        }
        my $count = 0;
        for ( my $i = 0; $i < $gene_length; $i++ ) {
            my $sub = substr( $nuc_seq, $i, 1 );
            if ( $sub =~ /G|C/i ) {
                $count++;
            }
        }

        my $GC = sprintf( "%.1f", $count * 100 / $gene_length );

        # my $left= substr( $nuc_seq, $res{pos}-1, $res{pos}-1 );

 # my @aminoAcids = map { exists $aacode{$_} ? $aacode{$_} : "?$_?" } @codons;

        print "-" x 50
            . "\nLocus: $locus_name ($start - $end) "
            . "\nProduct: "
            . $database{$locus_name}{"product"}
            . "\nNote: "
            . $database{$locus_name}{"note"}
            . "\nGene length (bp): "
            . $gene_length
            . "\nGC (%): " . "$GC\n"
            . "." x 50
            . "\nDNA :\n  "
            . $nuc_seq

            # . "\nRNA :\n  "
            # . "$rna\n"
            . "\n"
            . "." x 50
            . "\nProtein:\n  "
            . $aa_seq
            . "\nMolecular mass (Da) approx. : "
            . sum(@aaMw)
            . "\nProtein length (aa): "
            . scalar( length($aa_seq) - 1 ) . "\n";
        if ( $options{snp} ne "" ) {
            my %res = &get_snp_notation($snp);
            my $sub = $nuc_seq;
            substr( $sub, $res{pos} - 1, 1 ) = "[$res{ref}>$res{alt}]";
            print "." x 50 . "\nSNP:\n  $sub\n";
        }
        if ( $options{prosite} ) {
            &procite_patterns_parsing($aa_seq);
        }

    }

    else {
        print "Locus $locus_name not found!\n";
    }
}

sub procite_patterns_parsing {
#
# This function based on original script: parse_prosite created by James Tisdall, 2001
#
    my ($protein)    = @_;
    my $prosite_file = "prosite.dat";
    my $dir          = getcwd;
    my $record       = '';
    if ( -e "$dir/$prosite_file" ) {
        print "." x 50 . "\nProsite domains:\n";

        open( my $fh, $prosite_file )
            or die "Cannot open PROSITE file $prosite_file";
        #
        # set input separator to termination line //
        $/ = "//\n";

        while ( $record = <$fh> ) {

            #
            # Parse the PROSITE record into its "line types"
            #
            my %line_types = prosite_to_hash($record);

            #
            # Skip records without an ID (the first record)
            #
            defined $line_types{'ID'} or next;

            #
            # Skip records that are not PATTERN
            # (such as MATRIX or RULE)
            #
            $line_types{'ID'} =~ /PATTERN/ or next;

            #
            # Get the ID of this record
            #
            my $id = $line_types{'ID'};
            my $de = $line_types{'DE'};
            my $ac = $line_types{'AC'};
            $de =~ s/^DE   //;
            $ac =~ s/^AC   //;
            $ac =~ s/;//;
            $id =~ s/^ID   //;
            $id =~ s/; .*//;

            #
            # Get the PROSITE pattern from the PA line(s)
            #
            my $pattern = $line_types{'PA'};

            # Remove the PA line type tag(s)
            $pattern =~ s/PA   //g;

            #
            # Calculate the Perl regular expression
            # from the PROSITE pattern
            #
            my $regexp = prosite_get_pattern($pattern);

            #
            # Find the PROSITE regular expression patterns
            # in the protein sequence, and report
            #
            while ( $protein =~ /$regexp/g ) {
                my $res_protein    = $protein;
                my $pattern_length = length($&);
                my $position       = ( pos $protein ) - $pattern_length + 1;
                my $match          = $&;
                my $end_pos        = $position + ( length($match) - 1 );
                $res_protein =~ s/$match/\_\_$match\_\_/g;
                print "$de Found $id at positions $position:$end_pos\n";
                print "   $res_protein\n";
                print "   ID:   $ac\n";
                print "   match:   $match\n";
                print "   pattern: $pattern\n";
                print "   regexp:  $regexp\n\n";

            }

        }
    }
    else {
        print
            "prosite.dat file not found! Please, download the latest version of prosite.dat from ftp://ftp.expasy.org/databases/prosite/ and copy it into working directory!\n";
        exit;
    }
}

sub prosite_get_pattern {

# This function was copied from script parse_prosite created by James Tisdall, 2001

    my ($pattern) = @_;

    # print "@_\t";
    $pattern =~ s/{/[^/g;

    # print "$pattern\n";
    $pattern =~ tr/cx}()<>\-\./C.]{}^$/d;

    $pattern =~ s/\[G\$\]/(G|\$)/;

    # Return PROSITE pattern translated to Perl regular expression
    return $pattern;
}

#
# Parse a PROSITE record into "line types" hash
#
sub prosite_to_hash {

# This function was copied from script parse_prosite created by James Tisdall, 2001
#
# Collect the PROSITE record
#
    my ($record) = @_;

    #
    # Initialize the hash
    #   key   = line type
    #   value = lines
    #
    my %line_types_hash = ();

    #
    # Split the PROSITE record to an array of lines
    #
    my @records = split( /\n/, $record );

    #
    # Loop through the lines of the PROSITE record
    #
    foreach my $line (@records) {

        #
        # Extract the 2-character name
        # of the line type
        #
        my $line_type = substr( $line, 0, 2 );

        #
        # Append the line to the hash
        # indexed by this line type
        #
        ( defined $line_types_hash{$line_type} )
            ? ( $line_types_hash{$line_type} .= $line )
            : ( $line_types_hash{$line_type} = $line );
    }

    #
    # Return the hash
    #
    return %line_types_hash;
}

__END__


=encoding utf8


Usage:   snpMiner2 -db <database name> -action <type of> <command> <options>

=head1 -db

The select of database filename created by gb2db.pl script.

=over 20

=back

=head1 actions

=over 20

=item B<annotation >

Description....


=item B<uniq>

Description....

=item B<check_snp >

Description....

=item B<aa_cons>

Description....

=item B<make_seq_aa>  

Description....

=item B<make_seq_nc>  


=item --codon  

Description....


Description....


=back

=head1 commands

=over 20

<<<<<<< HEAD
=item B<-vcf> <vcf file>
=======
=item B<-vcf> <vcf file>   

Description....

=item B<-lvl> <1..6>   

Description....

=item B<-indel>    
>>>>>>> a1827220d481b5c1534176827f8ddb76a4b3f37e

Description....

=item B<-lvl> <1..6>

Description....

=item B<-indel>

Description....

=item B<--snp_list> <snp_list name>

Description....

=item B<--snp_list> <?>

Get a list of all snp list names.


=back
=cut
