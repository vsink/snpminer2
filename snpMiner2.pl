#!/usr/bin/env perl
use Pod::Usage;
use strict;
use warnings;
no warnings qw/uninitialized/;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;
use Storable qw(freeze thaw);
use Data::Dumper;
use Carp;
use IO::Compress::RawDeflate qw(rawdeflate $RawDeflateError);
use IO::Uncompress::RawInflate qw(rawinflate $RawInflateError);

# use Benchmark qw(:all) ;
# my $t0 = Benchmark->new;

my $genome_seq_obj      = "";    # BioSeq объект с геномом
my $tmp_res_codon_seq   = "";
my $tree_codon_sequence = "";
my $db_version          = "";
my $snp_list_exists     = 0;
my $ref_genome_name     = "";
my $organism_name       = "";
my $version="0.2.1";
my $about_msg="Program: snpMiner2 (Tool for analysis of variable call format files)\n";

#-------------------------------
my %database;
my %genetic_code;
my %tang_dic;
my %deltaH_dic;
my %GD_dic;
my %aa_hash;
my %snp_list_h;
my %snp_igpos_list_h;


# my %snps_h;
my %options;

GetOptions( \%options, 'db=s', 'action=s', 'snp_list=s', 'vcf=s', "about!", 'help' );

if ( $options{action} eq "rf" ) {
    &load_db( $options{db} );
    &write_dump;
}
elsif ($options{action} eq "annotation" && $options{vcf} ne ""
    || $options{action} eq "annotation.1" && $options{vcf} ne "" )
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

            &use_snp_classification;
        }
        elsif ( exists $snp_igpos_list_h{$options{snp_list} } ) {

            print "Yes, this is intergenic positions!!!\n";
        }
        else {
            print "The SNP classification name is not found!\n";
        }
    }
    else {
        print "Use -cn key to set SNP classification name!\n";
    }

}
elsif ( $options{action} eq "annotation.2" && $options{vcf} ne "" || $options{action} eq "annotation.3" && $options{vcf} ne "" ) {

    &load_db( $options{db} );

    &annotate_vcf_formatted( $options{vcf} );

}
elsif ($options{action} eq "uniq"
    || $options{action} eq "uniq.1"
    || $options{action} eq "uniq.2" )
{

    &load_db( $options{db} );

    &find_uniq_genes( $options{vcf} );

}
elsif ( $options{action} eq "uniq.3" ) {

    &load_db( $options{db} );

    &find_uniq_genes_formated( $options{vcf} );

}
elsif ( $options{action} eq "t3" ) {

    &load_db( $options{db} );

    &test3( $options{vcf} );

}
elsif ( $options{action} eq "t4" ) {

    &load_db( $options{db} );

    &test4( $options{vcf} );

}
elsif ( $options{action} eq "t5" || $options{action} eq "t5.1" ) {

    &load_db( $options{db} );

    &test5( $options{vcf} );

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
    my ($aa_input) = shift;

    if ( exists $aa_hash{$aa_input} ) {
        return $aa_hash{$aa_input};
    }

}

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

    #------------------------------------
    $db_version      = $database{options}{version};
    $ref_genome_name = $database{options}{reference};
    $organism_name   = $database{options}{organism_name};
    $snp_list_exists = scalar $database{options}{snp_list_exists};
    if ( $snp_list_exists == 1 ) {
        %snp_list_h = %{ $hash_ref->{dict}{snp_list}{coding} };
        %snp_igpos_list_h = %{ $hash_ref->{dict}{snp_list}{intergenic} };
        # %snp_list_type_h = %{ $hash_ref->{dict}{snp_list_type} };
    }

}

sub annotate_vcf {
    my $loc_file_in   = shift;
    my $AACI          = "";
    my $loc_AACI_flag = "---";
    my $loc_AACI      = 0;
    my $loc_ref_name  = "";

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
                && $tmp{'alt_aa_long'} ne 'STOP'
                && $tmp{'ref_aa_long'} ne 'STOP'
                && $tmp{'alt_aa_long'} ne 'BAD_CODON'
                && $tmp{'ref_aa_long'} ne 'BAD_CODON' )
            {

                if ( $tmp{'snp_type'} eq "missense" ) {
                    my $loc_tang_index
                        = calcutalte_tang_index( $tmp{'ref_aa_long'},
                        $tmp{'alt_aa_long'} );
                    my $loc_deltaH
                        = calcutalte_deltaH( $tmp{'ref_aa_long'},
                        $tmp{'alt_aa_long'} );
                    my $loc_GD
                        = calculate_grantham_matrix( $tmp{'ref_aa_long'},
                        $tmp{'alt_aa_long'} );
                    if ( $loc_tang_index < 0.5 ) {
                        $loc_AACI++;
                        substr $loc_AACI_flag, 0, 1, "+";
                    }
                    if ( $loc_deltaH > 1.22 ) {
                        $loc_AACI++;
                        substr $loc_AACI_flag, 1, 1, "+";
                    }
                    if ( $loc_GD > 100 ) {
                        $loc_AACI++;
                        substr $loc_AACI_flag, 2, 1, "+";
                    }
                    $AACI
                        = "$loc_tang_index\t$loc_deltaH\t$loc_GD\t$loc_AACI_flag\t";
                }
                else {
                    $AACI = "-\t-\t-\t$loc_AACI_flag\t";
                }
                if ( $options{action} eq "annotation.1" ) {
                    print
                        "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t$tmp{'snp'}\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$AACI$tmp{'product'}\n";
                }

       # elsif ( $options{action} eq "annotation.2" ) {
       #     print
       #         "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t\t$tmp{'product'}\n";
       # }
                elsif ( $options{action} eq "annotation" ) {
                    print
                        "$tmp{'locus'}\t$tmp{'snp_genome_pos'}\t$tmp{'snp'}\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$tmp{'product'}\n";
                }

                $loc_AACI      = 0;
                $loc_AACI_flag = "---";

            }

        }

    }
    close $fh;

}

sub use_snp_classification {
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
    my $position_type="";

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
            if ( length($2) == 1 && length($3) == 1) {
                # print "!\n" if ($1 =="1673425");
                # print &check_position($1). "\n";
                my %tmp = &get_locus_info( $+{genome_pos}, $3 );
                if ( $loc_ref_name ne $ref_genome_name ) {
                    print
                        "\t\tREFERENCE GENOME: $ref_genome_name NOT FOUND!!!\n";
                    last;
                }
                if (   $tmp{'locus'} ne ""
                    && $tmp{'alt_aa_long'} ne 'STOP'
                    && $tmp{'ref_aa_long'} ne 'STOP'
                    && $tmp{'alt_aa_long'} ne 'BAD_CODON'
                    && $tmp{'ref_aa_long'} ne 'BAD_CODON' )
                {
                    # $class_query
                    #     = uc( $tmp{'locus'} . "_" . $tmp{'snp'} );
                    $query_locus = uc( $tmp{'locus'} );
                    $query_snp   = uc( $tmp{'snp'} );

                    # print "$class_query\n";
                    $loc_snp_notation = &get_snp_notation(
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
                        $tmp{'ref_aa_long'}
                            . $tmp{'aa_pos'}
                            . $tmp{'alt_aa_long'},
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
                && $tmp{'alt_aa_long'} ne 'STOP'
                && $tmp{'ref_aa_long'} ne 'STOP'
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

                    $loc_AACI_flag  = &calc_aa_change_info( $tmp{'ref_aa_long'},
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
                && $tmp{'alt_aa_long'} ne 'STOP'
                && $tmp{'ref_aa_long'} ne 'STOP'
                && $tmp{'alt_aa_long'} ne 'BAD_CODON'
                && $tmp{'ref_aa_long'} ne 'BAD_CODON' )
            {
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
                && $tmp{'alt_aa_long'} ne 'STOP'
                && $tmp{'ref_aa_long'} ne 'STOP'
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
                $snp_per_gene_h{$loc_locus_name}{$loc_snp_pos}
                    = "\t$loc_snp_pos\t$loc_snp\t$tmp{'ref_codon'}/$tmp{'alt_codon'}\t$tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'}\t$tmp{'snp_type'}\t$loc_AACI_flag\t$files_per_snp_size/$count_f";

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
    print "Execution time: $duration s\n";

}

sub test3 {
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
                $all_snp_h{$1}{'file'} .= "$file\t";

                # $all_snp_h{$1}{'full'}     .= "$1\t$3\t$file\n";
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

        if ( $files_per_snp_size > 20 ) {
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
            && $tmp{'alt_aa_long'} ne 'STOP'
            && $tmp{'ref_aa_long'} ne 'STOP'
            && $tmp{'alt_aa_long'} ne 'BAD_CODON'
            && $tmp{'ref_aa_long'} ne 'BAD_CODON' )
        {
            print
                "$key\n\t$tmp{'locus'} $tmp{'ref_codon'}/$tmp{'alt_codon'} $tmp{'ref_aa_short'}$tmp{'aa_pos'}$tmp{'alt_aa_short'} $tmp{'snp_type'} $tmp{'product'} \n\t\t"
                . $results_h{$key}{'count'} . "/"
                . scalar(@files)
                . "\n\t\t\t"
                . $results_h{$key}{'file_list'} . "\n";

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
    print "Execution time: $duration s\n";

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
                    && $tmp{'alt_aa_long'} ne 'STOP'
                    && $tmp{'ref_aa_long'} ne 'STOP'
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
        print "$file...OK\n";

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

sub get_locus_info {
    my $loc_snp_pos    = shift;
    my $loc_alt        = shift;
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
        if ( $ref_aa eq $alt_aa ) {
            $snp_type = "synonymous";
        }
        elsif ( $ref_aa ne $alt_aa ) {
            $snp_type = "missense";
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
    my $loc_snp_pos    = shift;    
    my $type_pos="";

    foreach my $key ( keys %database ) {
        if ($loc_snp_pos > scalar $database{$key}{'start'}
            && $loc_snp_pos < scalar $database{$key}{'end'} )
        {
            $type_pos="coding";
            last;
        }
        else
        {
            $type_pos="intergenic";
        }

    }

    if ( $type_pos ne "" ) {
        return $type_pos;
    }
}

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

sub calc_aa_change_info {
    my $ref_aa_long   = shift;
    my $alt_aa_long   = shift;
    my $type          = shift;
    my $loc_AACI_flag = "---";
    my $loc_AACI;
    my $AACI;
    if ( $type eq "missense" ) {
        my $loc_tang_index
            = calcutalte_tang_index( $ref_aa_long, $alt_aa_long );
        my $loc_deltaH = calcutalte_deltaH( $ref_aa_long, $alt_aa_long );
        my $loc_GD = calculate_grantham_matrix( $ref_aa_long, $alt_aa_long );
        if ( $loc_tang_index < 0.5 ) {
            $loc_AACI++;
            substr $loc_AACI_flag, 0, 1, "+";
        }
        if ( $loc_deltaH > 1.22 ) {
            $loc_AACI++;
            substr $loc_AACI_flag, 1, 1, "+";
        }
        if ( $loc_GD > 100 ) {
            $loc_AACI++;
            substr $loc_AACI_flag, 2, 1, "+";
        }

        # return $loc_AACI_flag;

        $AACI = "$loc_tang_index\t$loc_deltaH\t$loc_GD\t$loc_AACI_flag\t";
        return $AACI;
    }
    else {
        # return "---";

        $AACI = "-\t-\t-\t$loc_AACI_flag\t";
        return $AACI;
    }

}

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

sub get_snp_notation {
    my $snp_class_name = shift;
    my $position_type=shift;
    my $query_locus    = shift;
    my $query_snp      = shift;
    my $position       = shift;
    my $aa             = shift;
    my $aa_long        = uc(shift);
    my $snp            = shift;

    my ($position_extracted) = $position =~ /^(\d+)/;
    my ($aa_reformat)        = $aa_long =~ /^(\D{3})(\d+)(\D{3})/;
    my ($aa_codon_pos)        = $aa_long =~ /^\D{3}(\d+)\D{3}/;
    my $aa_long_genome_pos   = "$position_extracted\_$1>$3";
    my $pos_check            = "POS$position_extracted";
    my $codon_check="CODON$aa_codon_pos";
    

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
        return $snp_list_h{$snp_class_name}{$query_locus}{$codon_check}{'note'}
        . "\ ($query_locus\_$snp\_$aa_long)";
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


