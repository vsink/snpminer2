#!/usr/bin/env perl
use Pod::Usage;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;
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
no warnings qw/uninitialized/;

# use Benchmark qw(:all) ;
# my $t0 = Benchmark->new;

my $genome_seq_obj      = "";        # BioSeq объект с геномом
my $tmp_res_codon_seq   = "";
my $tree_codon_sequence = "";
my $db_version          = "";
my $snp_list_exists     = 0;
my $ref_genome_name     = "";
my $organism_name       = "";
my $version             = "0.2.2";
my $about_msg
    = "Program: snpMiner2 (Tool for analysis of variable call format files)\n";

#-------------------------------
my %database;
my %genetic_code;
my %tang_dic;
my %deltaH_dic;
my %GD_dic;
my %aa_hash;
my %snp_list_h;
my %snp_igpos_list_h;
my %aa_mw;

my %options;    # options hash

GetOptions(
    \%options, 'db=s',    'action=s', 'snp_list=s',
    'vcf=s',   "about!",  'help',     'snp_th=s',
    'list=s',  "DP=s",    "target=s", "tang_th=s",
    "locus=s", "prosite", "o=s"
);

&show_help( "verbose" => 99, -utf8 ) if $options{help};

if ( $options{action} eq "rf" ) {
    &load_db( $options{db} );
    &write_dump;
}
elsif ($options{action} eq "annotation" && $options{vcf} ne ""
    || $options{action} eq "annotation.1"  && $options{vcf} ne ""
    || $options{action} eq "annotation.11" && $options{vcf} ne "" )
{

    &load_db( $options{db} );

    &annotate_vcf( $options{vcf} );

}
elsif ( $options{action} eq "check_snp" ) {

    if ( $options{snp_list} ne "" && $options{db} ne "" ) {
        &load_db( $options{db} );

        #         foreach my $key (keys %snp_list_h){
        #     print "$key\n";
        # }
        # exit;
        if ( exists $snp_list_h{ $options{snp_list} } ) {

            &use_snp_classification_coding;
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
        print "Use -cn key to set SNP classification name!\n";
    }

}
elsif ($options{action} eq "annotation.2" && $options{vcf} ne ""
    || $options{action} eq "annotation.3" && $options{vcf} ne "" )
{

    &load_db( $options{db} );

    &annotate_vcf_formatted( $options{vcf} );

}
elsif ($options{action} eq "uniq"
    || $options{action} eq "uniq.1"
    || $options{action} eq "uniq.2"
    || $options{action} eq "uniq.6"
    || $options{action} eq "uniq.7" )
{

    &load_db( $options{db} );

    &find_uniq_genes( $options{vcf} );

}
elsif ($options{action} eq "uniq.3"
    || $options{action} eq "uniq.4" )
{

    &load_db( $options{db} );

    &find_uniq_genes_formated( $options{vcf} );

}
elsif ( $options{action} eq "uniq.5" ) {

    &load_db( $options{db} );

    &compare_files_by_snp( $options{vcf} );

}
elsif ( $options{action} eq "check_list" ) {

    &load_db( $options{db} );

    &compare_list_by_snp( $options{vcf} );

}
elsif ( $options{action} eq "t4" ) {

    &load_db( $options{db} );

    &test4( $options{vcf} );

}
elsif ( $options{action} eq "t5" || $options{action} eq "t5.1" ) {

    &load_db( $options{db} );

    &test5( $options{vcf} );

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
        . "\nThe snpMiner2 was written by Viacheslav V.Sinkov, 2015\n\t\twebsite: http://snpminer.ru\n"
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

    if ( exists $genetic_code{$codon} ) {
        return $genetic_code{$codon};
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
    %tang_dic     = %{ $hash_ref->{dict}{tang} };     # Tang index values
    %deltaH_dic   = %{ $hash_ref->{dict}{deltaH} };   # Delta H values
    %GD_dic       = %{ $hash_ref->{dict}{GD} };       # Grantham values
    %genetic_code = %{ $hash_ref->{dict}{codon} };    #codon table
    %aa_hash      = %{ $hash_ref->{dict}{aa} };       #amino acid table
    %aa_mw        = %{ $hash_ref->{dict}{aa_mw} };    #amino acid table

    #------------------------------------
    $db_version      = $database{options}{version};
    $ref_genome_name = $database{options}{reference};
    $organism_name   = $database{options}{organism_name};
    $snp_list_exists
        = scalar $database{options}{snp_list_exists}; #snp_list existance flag
    my $igpos_list_exists = scalar $database{options}{igpos_list_exists}
        ;                                             #snp_list existance flag
    if ( $snp_list_exists == 1 ) {
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
sub annotate_vcf {
    my $loc_file_in   = shift;
    my $AACI          = "";
    my $loc_AACI_flag = "--";
    my $loc_AACI      = 0;
    my $loc_ref_name  = "";
    my $DP            = "";
    print
        "locus\tsnp_pos\tgene_pos\tcodons\taa_pos\ttype\tu_index\tGD\tscore\tDP\tproduct\n";
    open( my $fh, "<", $loc_file_in )
        or croak "cannot open < $loc_file_in: $!";

    while (<$fh>) {
        my $str = "$_";
        $str =~ /.*($ref_genome_name)/x;
        $loc_ref_name = $1;
        $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;
        if ( length($2) == 1 && length($3) == 1 ) {
            my %tmp = &get_locus_info( $1, $3 );
            if ( $loc_ref_name ne $ref_genome_name ) {
                print "REFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
                exit;
            }

 # print "$tmp{'alt_aa_long'}\t!!!!!\n" if $tmp{'alt_aa_long'} eq "BAD_CODON";
            if (   $tmp{'locus'} ne ""
                # && $tmp{'alt_aa_long'} ne 'STOP'
                # && $tmp{'ref_aa_long'} ne 'STOP'
                && $tmp{'alt_aa_long'} ne 'BAD_CODON'
                && $tmp{'ref_aa_long'} ne 'BAD_CODON' )
            {
                $DP = $4;                
                if ( $tmp{'snp_type'} eq "missense" ) {
                    my $loc_tang_index
                        = calcutalte_tang_index( $tmp{'ref_aa_long'},
                        $tmp{'alt_aa_long'} );

                    # my $loc_deltaH
                    #     = calcutalte_deltaH( $tmp{'ref_aa_long'},
                    #     $tmp{'alt_aa_long'} );
                    my $loc_GD
                        = calculate_grantham_matrix( $tmp{'ref_aa_long'},
                        $tmp{'alt_aa_long'} );

                    # my $titv=&calcutalte_tv_ti($tmp{'ref'},
                    #     $tmp{'alt'} );
                    #  $titv="-" if $titv eq "";
                    if ( $loc_tang_index < 0.5 ) {
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
                    $AACI = "$loc_tang_index\t$loc_GD\t$loc_AACI_flag\t$DP\t";
                }

                else {
                    $AACI = "-\t-\t$loc_AACI_flag\t$DP\t";
                }
                if ( $options{action} eq "annotation.1" ) {

                    # if ( $loc_AACI_flag eq "+++" ) {
                    print
                        "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t$tmp{'snp'}\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$AACI$tmp{'product'}\n";

                    # }
                }

       # elsif ( $options{action} eq "annotation.2" ) {
       #     print
       #         "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t\t$tmp{'product'}\n";
       # }
                elsif ( $options{action} eq "annotation" ) {

                    # print  color 'bold green';
                    # if ( $loc_AACI_flag ne "---" ) {

    print
        "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t$tmp{'snp'}\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$tmp{'product'}\n";
# }
                }
                elsif ( $options{action} eq "annotation.11" ) {
                    print
                        "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t$tmp{'snp'}\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$AACI$tmp{'product'}\n"
                        if $loc_AACI_flag ne "--";
                }

                $loc_AACI      = 0;
                $loc_AACI_flag = "--";

            }

        }

    }
    close $fh;

}

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
    my $loc_ref_name     = "";
    my $loc_snp_notation = "";
    my $position_type    = "";

    foreach my $file (@files) {
        print "$file";
        my @res_a;
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*(?<ref_genome>$ref_genome_name)/x;
            if (   $+{ref_genome} eq "$ref_genome_name"
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
            elsif ($+{ref_genome} ne "$ref_genome_name"
                && $+{ref_genome} ne "" )
            {
                print "\t\tREFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
                last;
            }

     #         $loc_ref_name = $+{ref_genome};
     #         $str =~ /^\S+\W+(?<genome_pos>\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

            #         # $position_type= &check_position($+{genome_pos});
            #         if ( length($2) == 1 && length($3) == 1 ) {

#             # print "!\n" if ($1 =="1673425");
#             # print &check_position($1). "\n";
#             if ( $loc_ref_name ne $ref_genome_name ) {
#                 print
#                     "\t\tREFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
#                 last;
#             }

            #             print $+{genome_pos} . "\n";

            #             # $class_query
            #             #     = uc( $tmp{'locus'} . "_" . $tmp{'snp'} );

            #         }
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
    my $loc_ref_name     = "";
    my $loc_snp_notation = "";
    my $position_type    = "";

    foreach my $file (@files) {
        print "$file";
        my @res_a;
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*(?<ref_genome>$ref_genome_name)/x;
            $loc_ref_name = $+{ref_genome};
            $str =~ /^\S+\W+(?<genome_pos>\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

            # $position_type= &check_position($+{genome_pos});
            if ( length($2) == 1 && length($3) == 1 ) {

                # print "!\n" if ($1 =="1673425");
                # print &check_position($1). "\n";
                my %tmp = &get_locus_info( $+{genome_pos}, $3 );
                if ( $loc_ref_name ne $ref_genome_name ) {
                    print
                        "\t\tREFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
                    last;
                }
                if (   $tmp{'locus'} ne ""
                    # && $tmp{'alt_aa_long'} ne 'STOP'
                    # && $tmp{'ref_aa_long'} ne 'STOP'
                    && $tmp{'alt_aa_long'} ne 'BAD_CODON'
                    && $tmp{'ref_aa_long'} ne 'BAD_CODON' )
                {
                    # $class_query
                    #     = uc( $tmp{'locus'} . "_" . $tmp{'snp'} );
                    $query_locus = uc( $tmp{'locus'} );
                    $query_snp   = uc( $tmp{'snp'} );

                    # print "$class_query\n";
                    $loc_snp_notation = &load_coding_snp(
                        $options{snp_list},
                        $position_type,
                        $query_locus,
                        $query_snp,
                        $tmp{'snp_genome_pos'} . "_"
                            . $tmp{'ref'} . ">"
                            . $tmp{'alt'},
                        $tmp{'ref_aa_short'}
                            . $tmp{'aa_pos'}
                            . $tmp{'alt_aa_short'},
                        uc(       $tmp{'ref_aa_long'}
                                . $tmp{'aa_pos'}
                                . $tmp{'alt_aa_long'}
                        ),
                        $tmp{'snp'}

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

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub annotate_vcf_formatted {
    my $loc_file_in = shift;
    my %snp_per_gene_h;

    # my @snp_per_gene_a;
    my $loc_ref_name  = "";
    my $loc_AACI_flag = "---";
    open( my $fh, "<", $loc_file_in )
        or croak "cannot open < $loc_file_in: $!";
    while (<$fh>) {
        my $str = "$_";
        $str =~ /.*($ref_genome_name)/x;
        $loc_ref_name = $1;
        $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;
        if ( length($2) == 1 && length($3) == 1 ) {
            my %tmp = &get_locus_info( $1, $3 );
            if ( $loc_ref_name ne $ref_genome_name ) {
                print "REFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
                exit;
            }

 # print "$tmp{'alt_aa_long'}\t!!!!!\n" if $tmp{'alt_aa_long'} eq "BAD_CODON";
            if (   $tmp{'locus'} ne ""
                # && $tmp{'alt_aa_long'} ne 'STOP'
                # && $tmp{'ref_aa_long'} ne 'STOP'
                && $tmp{'alt_aa_long'} ne 'BAD_CODON'
                && $tmp{'ref_aa_long'} ne 'BAD_CODON' )
            {

                # $snp_per_gene_h{$tmp{'locus'}}=$tmp{'snp'};
                my $loc_locus_name = $tmp{'locus'};
                my $loc_snp        = $tmp{'snp'};
                my $loc_snp_pos    = $tmp{'snp_genome_pos'};
                if ( $options{action} eq "annotation.2" ) {
                    $snp_per_gene_h{$loc_locus_name}{$loc_snp_pos}
                        = "$loc_snp\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}";
                }
                elsif ( $options{action} eq "annotation.3" ) {

                    $loc_AACI_flag
                        = &calc_aa_change_info( $tmp{'ref_aa_long'},
                        $tmp{'alt_aa_long'}, $tmp{'snp_type'} );

                    $snp_per_gene_h{$loc_locus_name}{$loc_snp_pos}
                        = "$loc_snp\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$loc_AACI_flag";
                }
            }

        }

    }
    close $fh;

    # print Dumper (%snp_per_gene_h);

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

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub find_uniq_genes {
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
    my $loc_ref_name = "";
    my %all_snp_h;
    my @results_a;

    my $start = time;

    foreach my $file (@files) {

        # print "$file";
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*($ref_genome_name)/x;
            $loc_ref_name = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

            if ( length($2) == 1 && length($3) == 1 ) {

                # my %tmp = &get_locus_info( $1, $3 );
                if ( $loc_ref_name ne $ref_genome_name ) {
                    print
                        "$file\tREFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
                    last;
                }
                $all_snp_h{$1}{'file'}     .= "$file\t";
                $all_snp_h{$1}{'alt_list'} .= "$3\t";
                $all_snp_h{$1}{'alt'} = "$3";
                $all_snp_h{$1}{'DP'} .= "$4\t";

            }

        }
        close $fh;

    }
    foreach my $key ( sort keys %all_snp_h ) {
        my $files_per_snp_size;
        my @filenames_a = split /\t/, $all_snp_h{$key}{'file'};
        my @alts        = split /\t/, $all_snp_h{$key}{'alt_list'};
        my @DPs         = split /\t/, $all_snp_h{$key}{'DP'};
        my $alt           = $all_snp_h{$key}{'alt'};
        my $loc_AACI_flag = "";
        $files_per_snp_size = @filenames_a;
        my $DP_min         = min @DPs;
        my $DP_max         = max @DPs;
        my $DP_average_raw = ( sum @DPs ) / $files_per_snp_size;
        my $DP_average     = sprintf( "%.1i", $DP_average_raw );

        if ( $files_per_snp_size == $count_f ) {
            my %tmp = &get_locus_info( $key, $alt );

            if (   $tmp{'locus'} ne ""
                # && $tmp{'alt_aa_long'} ne 'STOP'
                # && $tmp{'ref_aa_long'} ne 'STOP'
                && $tmp{'alt_aa_long'} ne 'BAD_CODON'
                && $tmp{'ref_aa_long'} ne 'BAD_CODON' )
            {
                $loc_AACI_flag
                    = &calc_aa_change_info( $tmp{'ref_aa_long'},
                    $tmp{'alt_aa_long'}, $tmp{'snp_type'} );
                if ( $options{action} eq "uniq" ) {
                    push( @results_a,
                        "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t$tmp{'snp'}\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$tmp{'product'}\t$files_per_snp_size/$count_f\n"
                    );
                }
                elsif ( $options{action} eq "uniq.1" ) {
                    push( @results_a,
                        "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t$tmp{'snp'}\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$tmp{'product'}\n"
                    );

                }
                elsif ( $options{action} eq "uniq.6" ) {
                    push( @results_a,
                        "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t$tmp{'snp'}\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$loc_AACI_flag$tmp{'product'}\t\n"
                    );

                }
                elsif ( $options{action} eq "uniq.7" ) {
                    push( @results_a,
                        "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t$tmp{'snp'}\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$loc_AACI_flag$tmp{'product'}\t"
                            . "DP:"
                            . $DP_min . "-"
                            . $DP_max
                            . " (Avg: $DP_average)"
                            . "\n" );

                }
                elsif ( $options{action} eq "uniq.2" ) {
                    push( @results_a,
                        "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t$tmp{'product'}\n"
                    );

                }

            }

        }

    }

    print sort @results_a;

    # my $duration = time - $start;
    # print "Execution time: $duration s\n";
}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub find_uniq_genes_formated {
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
    my $loc_ref_name = "";
    my %all_snp_h;
    my @results_a;
    my %snp_per_gene_h;

    my $start = time;

    foreach my $file (@files) {

        # print "$file";
        open( my $fh, "<", $file )
            or croak "cannot open < $file: $!";
        while (<$fh>) {
            my $str = "$_";
            $str =~ /.*($ref_genome_name)/x;
            $loc_ref_name = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

            if ( length($2) == 1 && length($3) == 1 ) {

                # my %tmp = &get_locus_info( $1, $3 );
                if ( $loc_ref_name ne $ref_genome_name ) {
                    print
                        "$file\tREFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
                    last;
                }
                $all_snp_h{$1}{'file'}     .= "$file\t";
                $all_snp_h{$1}{'alt_list'} .= "$3\t";
                $all_snp_h{$1}{'alt'} = "$3";

            }

        }
        close $fh;

    }
    foreach my $key ( sort keys %all_snp_h ) {
        my $files_per_snp_size;
        my @filenames_a = split /\t/, $all_snp_h{$key}{'file'};
        my @alts        = split /\t/, $all_snp_h{$key}{'alt_list'};
        my $alt = $all_snp_h{$key}{'alt'};
        $files_per_snp_size = @filenames_a;
        if ( $files_per_snp_size == $count_f ) {
            my %tmp = &get_locus_info( $key, $alt );
            if (   $tmp{'locus'} ne ""
                # && $tmp{'alt_aa_long'} ne 'STOP'
                # && $tmp{'ref_aa_long'} ne 'STOP'
                && $tmp{'alt_aa_long'} ne 'BAD_CODON'
                && $tmp{'ref_aa_long'} ne 'BAD_CODON' )
            {
# push( @results_a,
#     "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t$tmp{'snp'}\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$tmp{'product'}\t$files_per_snp_size/$count_f\n"
# );
                my $loc_locus_name = $tmp{'locus'};
                my $loc_snp        = $tmp{'snp'};
                my $loc_snp_pos    = $tmp{'snp_genome_pos'};
                my $loc_AACI_flag
                    = &calc_aa_change_info( $tmp{'ref_aa_long'},
                    $tmp{'alt_aa_long'}, $tmp{'snp_type'} );
                if ( $options{action} eq "uniq.3" ) {
                    $snp_per_gene_h{$loc_locus_name}{$loc_snp_pos}
                        = "\t$loc_snp_pos\t$loc_snp\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$loc_AACI_flag\t$files_per_snp_size/$count_f";
                }
                elsif ( $options{action} eq "uniq.4" ) {
                    $snp_per_gene_h{$loc_locus_name}{$loc_snp_pos}
                        = "\t$loc_snp_pos\t$loc_snp\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$loc_AACI_flag";
                }

            }

        }

    }

    # print sort @results_a;
    # print Dumper(%snp_per_gene_h);
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

            # .
            # .
            # Do something with $value
            # .
            # .
            # .
            # }
        }
    }

    my $duration = time - $start;

    # print "Execution time: $duration s\n";

}

# ------------------------------------------------------------
#
#
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
    my $loc_ref_name = "";
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
            $str =~ /.*($ref_genome_name)/x;
            $loc_ref_name = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

            if ( length($2) == 1 && length($3) == 1 ) {

                # my %tmp = &get_locus_info( $1, $3 );
                if ( $loc_ref_name ne $ref_genome_name ) {
                    print
                        "$file\tREFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
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
        my %tmp = &get_locus_info( $key, $results_h{$key}{'alt'} );
        if (   $tmp{'locus'} ne ""
            # && $tmp{'alt_aa_long'} ne 'STOP'
            # && $tmp{'ref_aa_long'} ne 'STOP'
            && $tmp{'alt_aa_long'} ne 'BAD_CODON'
            && $tmp{'ref_aa_long'} ne 'BAD_CODON' )
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
                        "$tmp{'locus'} $key $tmp{'ref_codon'}/$tmp{'alt_codon'} $tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'} $tmp{'snp_type'} $tmp{'product'} \n";

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
                    "$tmp{'locus'} $key $tmp{'ref_codon'}/$tmp{'alt_codon'} $tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'} $tmp{'snp_type'} $tmp{'product'} \n";
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
#         . "\nLOCUS:\n\t $tmp{'locus'} $key $tmp{'ref_codon'}/$tmp{'alt_codon'} $tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'} $tmp{'snp_type'} $tmp{'product'} \n"
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
    my $loc_ref_name = "";
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
            $str =~ /.*($ref_genome_name)/x;
            $loc_ref_name = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;

            if ( length($2) == 1 && length($3) == 1 ) {

                # my %tmp = &get_locus_info( $1, $3 );
                if ( $loc_ref_name ne $ref_genome_name ) {
                    print
                        "$file\tREFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
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
        my %tmp = &get_locus_info( $key, $results_h{$key}{'alt'} );
        if (   $tmp{'locus'} ne ""
            # && $tmp{'alt_aa_long'} ne 'STOP'
            # && $tmp{'ref_aa_long'} ne 'STOP'
            && $tmp{'alt_aa_long'} ne 'BAD_CODON'
            && $tmp{'ref_aa_long'} ne 'BAD_CODON' )
        {
            # if ( $options{color} == 1 ) {
            print "\n" . "-" x 50;

            # . "\nSNP_POS:\n\t $key";
            print "\nLOCUS:\n\t";
            print
                "$tmp{'locus'} $key $tmp{'ref_codon'}/$tmp{'alt_codon'} $tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'} $tmp{'snp_type'} $tmp{'product'} \n";

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
#         . "\nLOCUS:\n\t $tmp{'locus'} $key $tmp{'ref_codon'}/$tmp{'alt_codon'} $tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'} $tmp{'snp_type'} $tmp{'product'} \n"
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
    my $loc_ref_name     = "";
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
            $str =~ /.*($ref_genome_name)/x;
            $loc_ref_name = $1;
            $str =~ /^\S+\W+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)/x;
            if ( length($2) == 1 && length($3) == 1 ) {

                # my %tmp = &get_locus_info( $1, $3 );
                if ( $loc_ref_name ne $ref_genome_name ) {
                    print
                        "\t\tREFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
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

sub test5 {
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

    foreach my $file (@files) {

        # print "$file";
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
                if (   $tmp{'locus'} ne ""
                    # && $tmp{'alt_aa_long'} ne 'STOP'
                    # && $tmp{'ref_aa_long'} ne 'STOP'
                    && $tmp{'alt_aa_long'} ne 'BAD_CODON'
                    && $tmp{'ref_aa_long'} ne 'BAD_CODON'
                    && $4 > 10 )
                {
                    $char_count++;
                    $seq_per_file_h{$file}{'seq'} .= $tmp{"alt_aa_short"};

       # $seq_per_file_h{$file}{$tmp{"snp_genome_pos"}}= $tmp{"alt_aa_short"};
                    $seq_per_file_h{$file}{'aa_pos'}
                        .= $char_count . ":"
                        . $tmp{"locus"} . ":"
                        . $tmp{"snp_genome_pos"} . ":"
                        . $tmp{"ref_aa_short"}
                        . $tmp{"aa_pos"}
                        . $tmp{"alt_aa_short"} . ":"
                        . $4 . "|";
                    $seq_per_file_h{$file}{ $tmp{"snp_genome_pos"} }
                        = $tmp{"alt_aa_short"};

                    # if ($check_ref_used ==0){
                    $seq_ref_h{ $tmp{"snp_genome_pos"} }
                        = $tmp{"ref_aa_short"};
                    $all_pos_h{ $tmp{"snp_genome_pos"} }
                        = $tmp{"snp_genome_pos"};

                    # }

                }

            }

        }
        close $fh;
        $char_count = 0;
        print "$file\t[ ";
        print "OK";
        print " ]\n";

    }
    print "-" x 50 . "\n";

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
        print ">$file\t" . length($res_seq) . "\n$res_seq\n";
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

    # print "Execution time: $duration s\n";
}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub get_locus_info {

    # my $loc_snp_pos    = shift;
    # my $loc_alt        = shift;
    my ( $loc_snp_pos, $loc_alt ) = @_;
    my $loc_locus_name = "";
    my %locus_info_h;
    my $snp_gene_pos;
    my $snp_type = "";

    foreach my $key ( keys %database ) {
        if (   $loc_snp_pos > scalar $database{$key}{'start'}
            && $loc_snp_pos < scalar $database{$key}{'end'} )
        {
            $loc_locus_name = $key;
            last;
        }

    }

    if ( $loc_locus_name ne "" ) {
        my $loc_start        = $database{$loc_locus_name}{"start"};
        my $loc_end          = $database{$loc_locus_name}{"end"};
        my $loc_strand       = $database{$loc_locus_name}{"strand"};
        my $loc_transl_table = $database{$loc_locus_name}{"transl_table"};
        my $loc_product      = $database{$loc_locus_name}{"product"};
        my $loc_note         = $database{$loc_locus_name}{"note"};
        my $loc_nuc_seq      = $database{$loc_locus_name}{"sequence"};

        if ( $loc_strand == 1 ) {
            $snp_gene_pos = ( ( $loc_snp_pos - $loc_start ) + 1 );
        }
        elsif ( $loc_strand == -1 ) {
            $snp_gene_pos = ( ( $loc_end - $loc_snp_pos ) + 1 );
            $loc_alt =~ tr/ACGT/TGCA/;
        }
        my $loc_ref = substr( $loc_nuc_seq, $snp_gene_pos - 1, 1 );
        my $codon_nbr = int( ( $snp_gene_pos - 1 ) / 3 ) + 1;
        my @codons        = unpack( "(A3)*", $loc_nuc_seq );
        my $res_codon     = lc( $codons[ $codon_nbr - 1 ] );
        my $loc_codon_pos = ( ( $codon_nbr * 3 ) - $snp_gene_pos );
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
        if ( $ref_aa eq $alt_aa and $alt_aa ne "STOP") {
            $snp_type = "synonymous";
        }
        elsif ( $ref_aa ne $alt_aa and $alt_aa ne "STOP" ) {
            $snp_type = "missense";
        }
        elsif ( $ref_aa ne $alt_aa  and $alt_aa eq "STOP" ) {            
            $snp_type = "nonsense";
            
        }

        %locus_info_h = (
            "locus"          => $loc_locus_name,
            "snp_genome_pos" => $loc_snp_pos,
            "start"          => $loc_start,
            "end"            => $loc_end,
            "strand"         => $loc_strand,
            "product"        => $loc_product,
            "sequence"       => $loc_nuc_seq,
            "snp_gene_pos"   => $snp_gene_pos,
            "alt"            => $loc_alt,
            "snp"            => "$snp_gene_pos$loc_ref>$loc_alt",
            "ref"            => $loc_ref,
            "aa_pos"         => $codon_nbr,
            "ref_codon"      => $res_codon,
            "snp_type"       => $snp_type,
            "ref_aa_long"    => $ref_aa,
            "alt_aa_long"    => $alt_aa,
            "ref_aa_short"   => &aa_decode($ref_aa),
            "alt_aa_short"   => &aa_decode($alt_aa),
            "alt_codon"      => $alt_res_codon,
            "product"        => $loc_product,
            "note"           => $loc_note

        );
        return %locus_info_h;
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
    my $loc_snp_pos    = shift;
    my $loc_locus_name = "";

    foreach my $key ( keys %database ) {
        if (   $loc_snp_pos > scalar $database{$key}{'start'}
            && $loc_snp_pos < scalar $database{$key}{'end'} )
        {
            $loc_locus_name = $key;
            last;
        }

    }

    if ( $loc_locus_name ne "" ) {
        return $loc_locus_name;
    }
}

sub check_position {
    my $loc_snp_pos = shift;
    my $type_pos    = "";

    foreach my $key ( keys %database ) {
        if (   $loc_snp_pos > scalar $database{$key}{'start'}
            && $loc_snp_pos < scalar $database{$key}{'end'} )
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

    # print values $tmp{matrix}{R};
    if ( $aa_ref ne "" && $aa_alt ne "" ) {

        # return $GD_dic->{$aa_ref}->{$aa_alt};
        return $GD_dic{$aa_ref}{$aa_alt};

        # exit;

        # print $tmp->{GD}->{$aa_ref}->{$aa_alt} . "!!!!!!\n";
        # exit;
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
    my $loc_AACI_flag = "--";
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

        # return $loc_AACI_flag;

        $AACI = "$loc_tang_index\t$loc_GD\t$loc_AACI_flag\t";
        return $AACI;
    }
    else {
        # return "---";

        $AACI = "-\t-\t$loc_AACI_flag\t";
        return $AACI;
    }

}

# ------------------------------------------------------------
#
#
# ------------------------------------------------------------

sub check_snp_notation {
    my $mut = shift;
    $mut =~ /^(\d+)(\D)>(\D)/x;
    if ( $1 ne "" && $2 ne "" && $3 ne "" ) {
        return "nc_gene";
    }
    else {
        $mut =~ /^(\D)(\d+)(\D)/x;
        if ( $1 ne "" && $2 ne "" && $3 ne "" ) {
            return "aa";
        }
        else {

            $mut =~ /^(\d+)_(\w)>(\w)/x;
            if ( $1 ne "" && $2 ne "" && $3 ne "" ) {
                return "nc_genome";
            }
            else {
                $mut =~ /^^(\w{3})(\d+)(\w{3})/x;
                if ( $1 ne "" && $2 ne "" && $3 ne "" ) {
                    return "aa_long";
                }
                else {

                    return "bad";
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
        return $snp_list_h{$snp_class_name}{$query_locus}{$query_snp}{'note'};
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
            . "\nDNA :\n"
            . $nuc_seq
            . "\nRNA :\n"
            . "$rna\n"
            . "." x 50
            . "\nProtein:\n"
            . $aa_seq
            . "\nMolecular mass (Da) approx. : "
            . sum(@aaMw)
            . "\nProtein length (aa): "
            . length($aa_seq) . "\n";
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
        print "." x 50 . "\n";

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
                $res_protein =~ s/$match/\_\_$match\_\_/g;
                print "$de Found $id at position $position\n";
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

=back

=head1 commands

=over 20

=item B<--vcf> <vcf file>   

Description....

=item B<--snp_list> <snp_list name>   

Description....

=back
=cut
