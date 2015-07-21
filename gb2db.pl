#!/usr/bin/env perl
use Bio::SeqIO;
use Bio::Perl;
use Bio::DB::GenBank;
use Pod::Usage;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;
use Storable qw(freeze thaw);
use Data::Dumper;
use IO::Compress::RawDeflate qw(rawdeflate $RawDeflateError);
use IO::Uncompress::RawInflate qw(rawinflate $RawInflateError);
use strict;

my $in_file       = "";
my $out_file      = "";
my $snp_list      = "";
my $version       = "0.2.170415";
my $line          = "-" x 50;
my $strain_name   = "";
my $organism_name = "";
my $organism_id   = "";
my $remote        = 0;

GetOptions(
    'i=s'   => \$in_file,
    'o=s'   => \$out_file,
    'id=s'  => \$organism_id,
    'snp=s' => \$snp_list

        # 'c=s' => \$codon_file,
);

if ( $in_file ne "" and $out_file ne "" ) {
    &create_db_file;
}
elsif ( $in_file eq "" and $out_file ne "" and $organism_id ne "" ) {

    #
    my $gb_file = "$organism_id.gb";
    if ( -e $gb_file ) {
        print
            "The file $gb_file already exists. If you want to use it press Y. If you want to download it again press N. Press X to exit.\n";
        my $input = uc(<STDIN>);
        chomp $input;
        if ( $input eq "N" || $input eq "NO" ) {
            $remote = 1;
            &create_db_file;
        }
        elsif ( $input eq "Y" || $input eq "YES" ) {
            $remote  = 0;
            $in_file = $gb_file;
            &create_db_file;
        }

    }
    else {
        $remote = 1;
        &create_db_file;
    }

}

else {
    print
        "$line\nUsage: gb2db.pl -i GENBANK_FILE -o OUT_FILE\n$line\nver.$version\n";
}

sub create_db_file {
    my (%annotations, $gene_feature, $start,     $end,
        $tag,         $locus_tag,    $gene,      $gene_synonym,
        $protein_id,  $product,      $note,      $translation_table,
        $seqobj,      $feat,         $tag_value, $nucleotide,
        $count_genes, $loc_strand,   %dict,      $seqio_object
    );
    if ( $remote == 0 ) {
        print "Parsing of $in_file file...\n";
        $seqio_object = Bio::SeqIO->new( -file => $in_file );
    }
    elsif ( $remote == 1 ) {
        print "Downloading $organism_id file from GenBank server...\n";

        my $genBank
            = new Bio::DB::GenBank; # This object knows how to talk to GenBank
        my $seq = $genBank->get_Seq_by_acc($organism_id)
            ;                       # get a record by accession
        my $out_file   = "$organism_id.gb";
        my $seq_object = Bio::SeqIO->new(
            -format  => 'genbank',
            '-file'  => ">$out_file",
            -verbose => -1
        );
        $seq_object->write_seq($seq);
        print "The file saved as $out_file\n";
        $seqio_object = Bio::SeqIO->new( -file => $out_file );

        # my $seqio_object = Bio::SeqIO->new(-format => 'genbank');

        # while ( my $seq = $in->next_seq() ) {
        #     $out->write_seq($seq);
        #     }
        # $seqio_object->write_seq($seq);

    }

    my $seq_object = $seqio_object->next_seq;

    # my $genome_seq     = $seq_object->seq();
    my $genome_seq_obj = Bio::Seq->new( -seq => $seq_object->seq() );
    my $display_id     = $seq_object->display_id();
    my $genome_length  = $seq_object->length();
    my $nuc_seq        = "";
    my $rev_seq        = "";

    print "-" x 50 . "\n";
    for my $feat_object ( $seq_object->get_SeqFeatures ) {

        if ( $feat_object->strand != 0 ) {
            $loc_strand = $feat_object->strand, "\n";    # -1,1
        }
        if ( $feat_object->primary_tag eq "source" ) {
            for my $val ( $feat_object->get_tag_values('organism') ) {

                $organism_name = $val;
                print "Organism: $organism_name\n";
            }
            for my $val ( $feat_object->get_tag_values('strain') ) {
                $strain_name = $val;
                print "Strain: $strain_name\n";
            }

        }

        if (   $feat_object->primary_tag eq "CDS"
            || $feat_object->primary_tag eq "rRNA"
            || $feat_object->primary_tag eq "tRNA" )
        {
            $start = $feat_object->location->start;
            $end   = $feat_object->location->end;
            $count_genes++;

            # print "start: $start\nstop: $end\n";
            # my $gene_length = ( $end - $start ) + 1;
            # my $gene_seq    = substr( $genome_seq, $start - 1,
            #         $gene_length );
            if ( $feat_object->has_tag('gene') ) {
                for my $val ( $feat_object->get_tag_values('gene') ) {

                    # print "gene: ", $val, "\n";
                    $gene = $val;

                    # $count_genes++;

                    # e.g. 'NDP', from a line like '/gene="NDP"'
                }
            }
            if ( $feat_object->has_tag('locus_tag') ) {
                for my $val ( $feat_object->get_tag_values('locus_tag') ) {

                    # print "locus_tag: ", $val, "\n";
                    $locus_tag = $val;

                    # e.g. 'NDP', from a line like '/gene="NDP"'
                }
            }
            if ( $feat_object->has_tag('note') ) {
                for my $val ( $feat_object->get_tag_values('note') ) {

                    # print "note: ", $val, "\n";
                    $note = $val;

                    # e.g. 'NDP', from a line like '/gene="NDP"'
                }

            }

            if ( $feat_object->has_tag('product') ) {
                for my $val ( $feat_object->get_tag_values('product') ) {

                    # print "product: ", $val, "\n";
                    $product = $val;

                    # e.g. 'NDP', from a line like '/gene="NDP"'
                }
            }

            if ( $feat_object->has_tag('transl_table') ) {
                for my $val ( $feat_object->get_tag_values('transl_table') ) {

                    # print "translation:\n", $val, "\n";
                    $translation_table = $val;

                    # e.g. 'NDP', from a line like '/gene="NDP"'
                }
            }

            if ( $loc_strand == 1 ) {
                $nuc_seq = $genome_seq_obj->subseq( $start, $end );
            }
            elsif ( $loc_strand == -1 ) {
                $nuc_seq
                    = revcom( $genome_seq_obj->subseq( $start, $end ) )
                    ->seq();
            }

            # $reverse_complement = revcom( $sequence );

            $annotations{$locus_tag}{gene} = $gene;

            # $annotations{$locus_tag}{gene_synonym} = $gene_synonym;
            $annotations{$locus_tag}{product} = $product;

            # $annotations{$locus_tag}{protein_id}   = $protein_id;
            $annotations{$locus_tag}{start}        = $start;
            $annotations{$locus_tag}{end}          = $end;
            $annotations{$locus_tag}{note}         = $note;
            $annotations{$locus_tag}{transl_table} = $translation_table;
            $annotations{$locus_tag}{strand}       = $loc_strand;
            $annotations{$locus_tag}{sequence}     = $nuc_seq;

            # $annotations{$locus_tag}{nc}       = $gene_seq;

            #TANG
            #-------------------------------------

# %hash = map { get_a_key_for($_) => $_ } @array;
# %new_hash = map { new_key($_) => new_value($hash{$_}) } keys %hash;
# my @tmp_tang= qw /ST QP EA LH RL VI SG SC KM GC SA AS QH RS RP PL NS VL RT EG RC DE RH IM VF NY IL AP QL EK SW NT KN LW DG SF YF RQ PH IF DV EQ SP TI SI CF LM AV LF GV NI TA DN SL RG CW RK TM KI EV CY KQ TP HY SY RW NH KT DA RI GW GA VM DH RM DY HF EH YQ IQ AR W E G Q Y/;
# my %tmp_tang_h= map {tang($_) => $_} @tmp_tang;

            ${dict}{tang}{ST} = 2.490;
            ${dict}{tang}{QP} = 1.377;
            ${dict}{tang}{EA} = 0.906;
            ${dict}{tang}{LH} = 0.560;
            ${dict}{tang}{RL} = 0.414;
            ${dict}{tang}{VI} = 2.415;
            ${dict}{tang}{SG} = 1.360;
            ${dict}{tang}{SC} = 0.852;
            ${dict}{tang}{KM} = 0.559;
            ${dict}{tang}{GC} = 0.414;
            ${dict}{tang}{SA} = 2.380;
            ${dict}{tang}{AS} = 2.380;
            ${dict}{tang}{QH} = 1.351;
            ${dict}{tang}{RS} = 0.850;
            ${dict}{tang}{RP} = 0.559;
            ${dict}{tang}{PL} = 0.388;
            ${dict}{tang}{NS} = 2.053;
            ${dict}{tang}{VL} = 1.329;
            ${dict}{tang}{RT} = 0.827;
            ${dict}{tang}{EG} = 0.553;
            ${dict}{tang}{RC} = 0.382;
            ${dict}{tang}{DE} = 2.033;
            ${dict}{tang}{RH} = 1.317;
            ${dict}{tang}{IM} = 0.827;
            ${dict}{tang}{VF} = 0.548;
            ${dict}{tang}{NY} = 0.378;
            ${dict}{tang}{IL} = 1.726;
            ${dict}{tang}{AP} = 1.288;
            ${dict}{tang}{QL} = 0.805;
            ${dict}{tang}{EK} = 0.548;
            ${dict}{tang}{SW} = 0.375;
            ${dict}{tang}{NT} = 1.695;
            ${dict}{tang}{KN} = 1.075;
            ${dict}{tang}{LW} = 0.793;
            ${dict}{tang}{DG} = 0.548;
            ${dict}{tang}{SF} = 0.365;
            ${dict}{tang}{YF} = 1.649;
            ${dict}{tang}{RQ} = 1.045;
            ${dict}{tang}{PH} = 0.784;
            ${dict}{tang}{IF} = 0.545;
            ${dict}{tang}{DV} = 0.361;
            ${dict}{tang}{EQ} = 1.634;
            ${dict}{tang}{SP} = 1.039;
            ${dict}{tang}{TI} = 0.750;
            ${dict}{tang}{SI} = 0.540;
            ${dict}{tang}{CF} = 0.321;
            ${dict}{tang}{LM} = 1.601;
            ${dict}{tang}{AV} = 1.017;
            ${dict}{tang}{LF} = 0.732;
            ${dict}{tang}{GV} = 0.539;
            ${dict}{tang}{NI} = 0.321;
            ${dict}{tang}{TA} = 1.587;
            ${dict}{tang}{DN} = 1.015;
            ${dict}{tang}{SL} = 0.725;
            ${dict}{tang}{RG} = 0.534;
            ${dict}{tang}{CW} = 0.271;
            ${dict}{tang}{RK} = 1.583;
            ${dict}{tang}{TM} = 1.007;
            ${dict}{tang}{KI} = 0.688;
            ${dict}{tang}{EV} = 0.506;
            ${dict}{tang}{CY} = 0.268;
            ${dict}{tang}{KQ} = 1.466;
            ${dict}{tang}{TP} = 1.001;
            ${dict}{tang}{HY} = 0.665;
            ${dict}{tang}{SY} = 0.503;
            ${dict}{tang}{RW} = 0.263;
            ${dict}{tang}{NH} = 1.382;
            ${dict}{tang}{KT} = 0.989;
            ${dict}{tang}{DA} = 0.657;
            ${dict}{tang}{RI} = 0.490;
            ${dict}{tang}{GW} = 0.242;
            ${dict}{tang}{GA} = 1.379;
            ${dict}{tang}{VM} = 0.986;
            ${dict}{tang}{DH} = 0.560;
            ${dict}{tang}{RM} = 0.470;
            ${dict}{tang}{DY} = 0.241;
            ${dict}{tang}{HF} = '';
            ${dict}{tang}{EH} = '';
            ${dict}{tang}{YQ} = '';
            ${dict}{tang}{IQ} = '';
            ${dict}{tang}{AR} = '';
            ${dict}{tang}{W}  = '';
            ${dict}{tang}{E}  = '';
            ${dict}{tang}{G}  = '';
            ${dict}{tang}{Q}  = '';
            ${dict}{tang}{Y}  = '';

            # Sneath dic
            #----------------------------------------
            ${dict}{sneath}{LI} = '0.889';
            ${dict}{sneath}{LV} = '0.785';
            ${dict}{sneath}{LG} = '0.380';
            ${dict}{sneath}{LA} = '0.643';
            ${dict}{sneath}{LP} = '0.432';
            ${dict}{sneath}{LQ} = '0.501';
            ${dict}{sneath}{LN} = '0.506';
            ${dict}{sneath}{LM} = '0.515';
            ${dict}{sneath}{LT} = '0.432';
            ${dict}{sneath}{LS} = '0.411';
            ${dict}{sneath}{LC} = '0.398';
            ${dict}{sneath}{LE} = '0.333';
            ${dict}{sneath}{LD} = '0.387';
            ${dict}{sneath}{LK} = '0.492';
            ${dict}{sneath}{LR} = '0.360';
            ${dict}{sneath}{LY} = '0.347';
            ${dict}{sneath}{LF} = '0.570';
            ${dict}{sneath}{LW} = '0.368';
            ${dict}{sneath}{LH} = '0.450';
            ${dict}{sneath}{IV} = '0.843';
            ${dict}{sneath}{IG} = '0.371';
            ${dict}{sneath}{IA} = '0.588';
            ${dict}{sneath}{IP} = '0.419';
            ${dict}{sneath}{IQ} = '0.453';
            ${dict}{sneath}{IN} = '0.456';
            ${dict}{sneath}{IM} = '0.494';
            ${dict}{sneath}{IT} = '0.493';
            ${dict}{sneath}{IS} = '0.360';
            ${dict}{sneath}{IC} = '0.348';
            ${dict}{sneath}{IE} = '0.366';
            ${dict}{sneath}{ID} = '0.338';
            ${dict}{sneath}{IK} = '0.477';
            ${dict}{sneath}{IR} = '0.342';
            ${dict}{sneath}{IY} = '0.266';
            ${dict}{sneath}{IF} = '0.487';
            ${dict}{sneath}{IW} = '0.287';
            ${dict}{sneath}{IH} = '0.368';
            ${dict}{sneath}{VG} = '0.437';
            ${dict}{sneath}{VA} = '0.675';
            ${dict}{sneath}{VP} = '0.473';
            ${dict}{sneath}{VQ} = '0.416';
            ${dict}{sneath}{VN} = '0.395';
            ${dict}{sneath}{VM} = '0.465';
            ${dict}{sneath}{VT} = '0.551';
            ${dict}{sneath}{VS} = '0.439';
            ${dict}{sneath}{VC} = '0.430';
            ${dict}{sneath}{VE} = '0.239';
            ${dict}{sneath}{VD} = '0.279';
            ${dict}{sneath}{VK} = '0.419';
            ${dict}{sneath}{VR} = '0.307';
            ${dict}{sneath}{VY} = '0.199';
            ${dict}{sneath}{VF} = '0.380';
            ${dict}{sneath}{VW} = '0.195';
            ${dict}{sneath}{VH} = '0.3';
            ${dict}{sneath}{GA} = '0.659';
            ${dict}{sneath}{GP} = '0.499';
            ${dict}{sneath}{GQ} = '0.163';
            ${dict}{sneath}{GN} = '0.19';
            ${dict}{sneath}{GM} = '0.149';
            ${dict}{sneath}{GT} = '0.396';
            ${dict}{sneath}{GS} = '0.323';
            ${dict}{sneath}{GC} = '0.29';
            ${dict}{sneath}{GE} = '0.049';
            ${dict}{sneath}{GD} = '0.015';
            ${dict}{sneath}{GK} = '0.309';
            ${dict}{sneath}{GR} = '0.149';
            ${dict}{sneath}{GY} = '0.163';
            ${dict}{sneath}{GF} = '0.259';
            ${dict}{sneath}{GW} = '0.138';
            ${dict}{sneath}{GH} = '0.183';
            ${dict}{sneath}{AP} = '0.533';
            ${dict}{sneath}{AQ} = '0.356';
            ${dict}{sneath}{AN} = '0.28';
            ${dict}{sneath}{AM} = '0.421';
            ${dict}{sneath}{AT} = '0.417';
            ${dict}{sneath}{AS} = '0.477';
            ${dict}{sneath}{AC} = '0.578';
            ${dict}{sneath}{AE} = '0.159';
            ${dict}{sneath}{AD} = '0.156';
            ${dict}{sneath}{AK} = '0.426';
            ${dict}{sneath}{AR} = '0.315';
            ${dict}{sneath}{AY} = '0.214';
            ${dict}{sneath}{AF} = '0.356';
            ${dict}{sneath}{AW} = '0.224';
            ${dict}{sneath}{AH} = '0.32';
            ${dict}{sneath}{PQ} = '0.168';
            ${dict}{sneath}{PN} = '0.172';
            ${dict}{sneath}{PM} = '0.265';
            ${dict}{sneath}{PT} = '0.330';
            ${dict}{sneath}{PS} = '0.321';
            ${dict}{sneath}{PC} = '0.318';
            ${dict}{sneath}{PE} = '0.003';
            ${dict}{sneath}{PD} = '0.015';
            ${dict}{sneath}{PK} = '0.295';
            ${dict}{sneath}{PR} = '0.155';
            ${dict}{sneath}{PY} = '0.179';
            ${dict}{sneath}{PF} = '0.282';
            ${dict}{sneath}{PW} = '0.211';
            ${dict}{sneath}{PH} = '0.172';
            ${dict}{sneath}{QN} = '0.589';
            ${dict}{sneath}{QM} = '0.699';
            ${dict}{sneath}{QT} = '0.340';
            ${dict}{sneath}{QS} = '0.501';
            ${dict}{sneath}{QC} = '0.482';
            ${dict}{sneath}{QE} = '0.685';
            ${dict}{sneath}{QD} = '0.492';
            ${dict}{sneath}{QK} = '0.545';
            ${dict}{sneath}{QR} = '0.561';
            ${dict}{sneath}{QY} = '0.368';
            ${dict}{sneath}{QF} = '0.459';
            ${dict}{sneath}{QW} = '0.353';
            ${dict}{sneath}{QH} = '0.406';
            ${dict}{sneath}{NM} = '0.518';
            ${dict}{sneath}{NT} = '0.488';
            ${dict}{sneath}{NS} = '0.581';
            ${dict}{sneath}{NC} = '0.485';
            ${dict}{sneath}{NE} = '0.578';
            ${dict}{sneath}{ND} = '0.637';
            ${dict}{sneath}{NK} = '0.401';
            ${dict}{sneath}{NR} = '0.427';
            ${dict}{sneath}{NY} = '0.391';
            ${dict}{sneath}{NF} = '0.34';
            ${dict}{sneath}{NW} = '0.316';
            ${dict}{sneath}{NH} = '0.459';
            ${dict}{sneath}{MT} = '0.409';
            ${dict}{sneath}{MS} = '0.48';
            ${dict}{sneath}{MC} = '0.612';
            ${dict}{sneath}{ME} = '0.402';
            ${dict}{sneath}{MD} = '0.292';
            ${dict}{sneath}{MK} = '0.482';
            ${dict}{sneath}{MR} = '0.522';
            ${dict}{sneath}{MY} = '0.307';
            ${dict}{sneath}{MF} = '0.465';
            ${dict}{sneath}{MW} = '0.355';
            ${dict}{sneath}{MH} = '0.345';
            ${dict}{sneath}{TS} = '0.668';
            ${dict}{sneath}{TC} = '0.485';
            ${dict}{sneath}{TE} = '0.218';
            ${dict}{sneath}{TD} = '0.254';
            ${dict}{sneath}{TK} = '0.224';
            ${dict}{sneath}{TR} = '0.258';
            ${dict}{sneath}{TY} = '0.285';
            ${dict}{sneath}{TF} = '0.254';
            ${dict}{sneath}{TW} = '0.176';
            ${dict}{sneath}{TH} = '0.208';
            ${dict}{sneath}{SC} = '0.613';
            ${dict}{sneath}{SE} = '0.312';
            ${dict}{sneath}{SD} = '0.330';
            ${dict}{sneath}{SK} = '0.285';
            ${dict}{sneath}{SR} = '0.317';
            ${dict}{sneath}{SY} = '0.354';
            ${dict}{sneath}{SF} = '0.380';
            ${dict}{sneath}{SW} = '0.243';
            ${dict}{sneath}{SH} = '0.342';
            ${dict}{sneath}{CE} = '0.221';
            ${dict}{sneath}{CD} = '0.243';
            ${dict}{sneath}{CK} = '0.269';
            ${dict}{sneath}{CR} = '0.324';
            ${dict}{sneath}{CY} = '0.223';
            ${dict}{sneath}{CF} = '0.289';
            ${dict}{sneath}{CW} = '0.188';
            ${dict}{sneath}{CH} = '0.288';
            ${dict}{sneath}{ED} = '0.84';
            ${dict}{sneath}{EK} = '0.435';
            ${dict}{sneath}{ER} = '0.382';
            ${dict}{sneath}{EY} = '0.261';
            ${dict}{sneath}{EF} = '0.219';
            ${dict}{sneath}{EW} = '0.086';
            ${dict}{sneath}{EH} = '0.201';
            ${dict}{sneath}{DK} = '0.248';
            ${dict}{sneath}{DR} = '0.236';
            ${dict}{sneath}{DY} = '0.287';
            ${dict}{sneath}{DF} = '0.172';
            ${dict}{sneath}{DW} = '0.028';
            ${dict}{sneath}{DH} = '0.2';
            ${dict}{sneath}{KR} = '0.733';
            ${dict}{sneath}{KY} = '0.285';
            ${dict}{sneath}{KF} = '0.381';
            ${dict}{sneath}{KW} = '0.297';
            ${dict}{sneath}{KH} = '0.421';
            ${dict}{sneath}{RY} = '0.407';
            ${dict}{sneath}{RF} = '0.339';
            ${dict}{sneath}{RW} = '0.288';
            ${dict}{sneath}{RH} = '0.396';
            ${dict}{sneath}{YF} = '0.729';
            ${dict}{sneath}{YW} = '0.565';
            ${dict}{sneath}{YH} = '0.504';
            ${dict}{sneath}{FW} = '0.741';
            ${dict}{sneath}{FH} = '0.605';
            ${dict}{sneath}{WH} = '0.484';

            #-------------------------------------

            ${dict}{deltaH}{PS} = '2.56';
            ${dict}{deltaH}{SY} = '2.83';
            ${dict}{deltaH}{TI} = '2.53';
            ${dict}{deltaH}{ED} = '0.01';
            ${dict}{deltaH}{IN} = '2.96';
            ${dict}{deltaH}{PT} = '2.16';
            ${dict}{deltaH}{ST} = '0.40';
            ${dict}{deltaH}{TA} = '0.19';
            ${dict}{deltaH}{HR} = '0.67';
            ${dict}{deltaH}{VF} = '0.96';
            ${dict}{deltaH}{PA} = '1.97';
            ${dict}{deltaH}{SI} = '2.93';
            ${dict}{deltaH}{ZQ} = '2.32';
            ${dict}{deltaH}{HD} = '0.86';
            ${dict}{deltaH}{VA} = '1.06';
            ${dict}{deltaH}{PL} = '0.18';
            ${dict}{deltaH}{SF} = '2.61';
            ${dict}{deltaH}{LH} = '1.02';
            ${dict}{deltaH}{HN} = '1.39';
            ${dict}{deltaH}{VL} = '0.73';
            ${dict}{deltaH}{PQ} = '2.50';
            ${dict}{deltaH}{SA} = '0.59';
            ${dict}{deltaH}{LM} = '1.12';
            ${dict}{deltaH}{MR} = '0.57';
            ${dict}{deltaH}{TK} = '1.06';
            ${dict}{deltaH}{PH} = '1.20';
            ${dict}{deltaH}{SL} = '2.38';
            ${dict}{deltaH}{LR} = '1.69';
            ${dict}{deltaH}{RC} = '0.08';
            ${dict}{deltaH}{TM} = '0.86';
            ${dict}{deltaH}{PR} = '1.87';
            ${dict}{deltaH}{SR} = '0.69';
            ${dict}{deltaH}{LW} = '0.58';
            ${dict}{deltaH}{RW} = '2.27';
            ${dict}{deltaH}{TR} = '0.29';
            ${dict}{deltaH}{GS} = '0.04';
            ${dict}{deltaH}{SN} = '0.03';
            ${dict}{deltaH}{KQ} = '1.40';
            ${dict}{deltaH}{DN} = '0.53';
            ${dict}{deltaH}{TN} = '0.44';
            ${dict}{deltaH}{GV} = '1.69';
            ${dict}{deltaH}{SC} = '0.61';
            ${dict}{deltaH}{KE} = '0.95';
            ${dict}{deltaH}{CW} = '2.35';
            ${dict}{deltaH}{AE} = '0.08';
            ${dict}{deltaH}{GA} = '0.63';
            ${dict}{deltaH}{SW} = '2.96';
            ${dict}{deltaH}{KM} = '0.20';
            ${dict}{deltaH}{IV} = '1.28';
            ${dict}{deltaH}{AD} = '0.09';
            ${dict}{deltaH}{GE} = '0.55';
            ${dict}{deltaH}{YF} = '0.22';
            ${dict}{deltaH}{KR} = '0.77';
            ${dict}{deltaH}{IF} = '0.32';
            ${dict}{deltaH}{VE} = '1.14';
            ${dict}{deltaH}{GR} = '0.73';
            ${dict}{deltaH}{YH} = '1.47';
            ${dict}{deltaH}{KN} = '1.49';
            ${dict}{deltaH}{IL} = '0.55';
            ${dict}{deltaH}{VM} = '0.39';
            ${dict}{deltaH}{GD} = '0.54';
            ${dict}{deltaH}{YD} = '2.33';
            ${dict}{deltaH}{QE} = '0.45';
            ${dict}{deltaH}{IK} = '1.47';
            ${dict}{deltaH}{VD} = '1.15';
            ${dict}{deltaH}{GC} = '0.65';
            ${dict}{deltaH}{YN} = '2.86';
            ${dict}{deltaH}{QH} = '1.30';
            ${dict}{deltaH}{IM} = '1.67';
            ${dict}{deltaH}{FL} = '0.23';
            ${dict}{deltaH}{GW} = '3.00';
            ${dict}{deltaH}{YC} = '2.22';
            ${dict}{deltaH}{QR} = '0.63';
            ${dict}{deltaH}{YR} = '2.24';
            ${dict}{deltaH}{FC} = '2.00';

            #------------------------------------
            ${dict}{codon}{TCA} = 'Ser';  # Serine
            ${dict}{codon}{TCC} = 'Ser';  # Serine
            ${dict}{codon}{TCG} = 'Ser';  # Serine
            ${dict}{codon}{TCT} = 'Ser';  # Serine
            ${dict}{codon}{TTC} = 'Phe';  # Phenylalanine
            ${dict}{codon}{TTT} = 'Phe';  # Phenylalanine
            ${dict}{codon}{TTA} = 'Leu';  # Leucineзамена кодонов
            ${dict}{codon}{TTG} = 'Leu';  # Leucine
            ${dict}{codon}{TAC} = 'Tyr';  # Tyrosine
            ${dict}{codon}{TAT} = 'Tyr';  # Tyrosine
            ${dict}{codon}{TAA} = 'STOP'; # Ochre
            ${dict}{codon}{TAG} = 'STOP'; # Amber
            ${dict}{codon}{TGC} = 'Cys';  # Cysteine
            ${dict}{codon}{TGT} = 'Cys';  # Cysteine
            ${dict}{codon}{TGA} = 'STOP'; # Opal
            ${dict}{codon}{TGG} = 'Trp';  # Tryptophan
            ${dict}{codon}{CTA} = 'Leu';  # Leucine
            ${dict}{codon}{CTC} = 'Leu';  # Leucine
            ${dict}{codon}{CTG} = 'Leu';  # Leucine
            ${dict}{codon}{CTT} = 'Leu';  # Leucine
            ${dict}{codon}{CCA} = 'Pro';  # Proline
            ${dict}{codon}{CCC} = 'Pro';  # Proline
            ${dict}{codon}{CCG} = 'Pro';  # Proline
            ${dict}{codon}{CCT} = 'Pro';  # Proline
            ${dict}{codon}{CAC} = 'His';  # Histidine
            ${dict}{codon}{CAT} = 'His';  # Histidine
            ${dict}{codon}{CAA} = 'Gln';  # Glutamine
            ${dict}{codon}{CAG} = 'Gln';  # Glutamine
            ${dict}{codon}{CGA} = 'Arg';  # Arginine
            ${dict}{codon}{CGC} = 'Arg';  # Arginine
            ${dict}{codon}{CGG} = 'Arg';  # Arginine
            ${dict}{codon}{CGT} = 'Arg';  # Arginine
            ${dict}{codon}{ATA} = 'Ile';  # Isoleucine
            ${dict}{codon}{ATC} = 'Ile';  # Isoleucine
            ${dict}{codon}{ATT} = 'Ile';  # Isoleucine
            ${dict}{codon}{ATG} = 'Met';  # Methionine
            ${dict}{codon}{ACA} = 'Thr';  # Threonine
            ${dict}{codon}{ACC} = 'Thr';  # Threonine
            ${dict}{codon}{ACG} = 'Thr';  # Threonine
            ${dict}{codon}{ACT} = 'Thr';  # Threonine
            ${dict}{codon}{AAC} = 'Asn';  # Asparagine
            ${dict}{codon}{AAT} = 'Asn';  # Asparagine
            ${dict}{codon}{AAA} = 'Lys';  # Lysine
            ${dict}{codon}{AAG} = 'Lys';  # Lysine
            ${dict}{codon}{AGC} = 'Ser';  # Serine
            ${dict}{codon}{AGT} = 'Ser';  # Serine
            ${dict}{codon}{AGA} = 'Arg';  # Arginine
            ${dict}{codon}{AGG} = 'Arg';  # Arginine
            ${dict}{codon}{GTA} = 'Val';  # Valine
            ${dict}{codon}{GTC} = 'Val';  # Valine
            ${dict}{codon}{GTG} = 'Val';  # Valine
            ${dict}{codon}{GTT} = 'Val';  # Valine
            ${dict}{codon}{GCA} = 'Ala';  # Alanine
            ${dict}{codon}{GCC} = 'Ala';  # Alanine
            ${dict}{codon}{GCG} = 'Ala';  # Alanine
            ${dict}{codon}{GCT} = 'Ala';  # Alanine
            ${dict}{codon}{GAC} = 'Asp';  # Aspartic Acid
            ${dict}{codon}{GAT} = 'Asp';  # Aspartic Acid
            ${dict}{codon}{GAA} = 'Glu';  # Glutamic Acid
            ${dict}{codon}{GAG} = 'Glu';  # Glutamic Acid
            ${dict}{codon}{GGA} = 'Gly';  # Glycine
            ${dict}{codon}{GGC} = 'Gly';  # Glycine
            ${dict}{codon}{GGG} = 'Gly';  # Glycine
            ${dict}{codon}{GGT} = 'Gly';  # Glycine
                 #------------------------------------
            ${dict}{aa}{Gly} = 'G';
            ${dict}{aa}{Leu} = 'L';
            ${dict}{aa}{Tyr} = 'Y';
            ${dict}{aa}{Ser} = 'S';
            ${dict}{aa}{Glu} = 'E';
            ${dict}{aa}{Gln} = 'Q';
            ${dict}{aa}{Asp} = 'D';
            ${dict}{aa}{Asn} = 'N';
            ${dict}{aa}{Phe} = 'F';
            ${dict}{aa}{Ala} = 'A';
            ${dict}{aa}{Lys} = 'K';
            ${dict}{aa}{Arg} = 'R';
            ${dict}{aa}{His} = 'H';
            ${dict}{aa}{Cys} = 'C';
            ${dict}{aa}{Val} = 'V';
            ${dict}{aa}{Pro} = 'P';
            ${dict}{aa}{Trp} = 'W';
            ${dict}{aa}{Ile} = 'I';
            ${dict}{aa}{Met} = 'M';
            ${dict}{aa}{Thr} = 'T';

            #------------------------------------

            ${dict}{aa_mw}{Ala} = '71.0788';
            ${dict}{aa_mw}{Arg} = '156.1875';
            ${dict}{aa_mw}{Asn} = '114.1038';
            ${dict}{aa_mw}{Asp} = '115.0886';
            ${dict}{aa_mw}{Cys} = '103.1388';
            ${dict}{aa_mw}{Glu} = '129.1155';
            ${dict}{aa_mw}{Gln} = '128.1307';
            ${dict}{aa_mw}{Gly} = '57.0519';
            ${dict}{aa_mw}{His} = '137.1411';
            ${dict}{aa_mw}{Ile} = '113.1594';
            ${dict}{aa_mw}{Leu} = '113.1594';
            ${dict}{aa_mw}{Lys} = '128.1741';
            ${dict}{aa_mw}{Met} = '131.1926';
            ${dict}{aa_mw}{Phe} = '147.1766';
            ${dict}{aa_mw}{Pro} = '97.1167';
            ${dict}{aa_mw}{Ser} = '87.0782';
            ${dict}{aa_mw}{Thr} = '101.1051';
            ${dict}{aa_mw}{Trp} = '186.2132';
            ${dict}{aa_mw}{Tyr} = '163.1760';
            ${dict}{aa_mw}{Val} = '99.1326';


            #------------------------------------
            
            my @grantham_matrix = qw/
                0 112 111 126 195  91 107  60  86  94  96 106  84 113  27  99  58 148 112  64
                112   0  86  96 180  43  54 125  29  97 102  26  91  97 103 110  71 101  77  96
                111  86   0  23 139  46  42  80  68 149 153  94 142 158  91  46  65 174 143 133
                126  96  23   0 154  61  45  94  81 168 172 101 160 177 108  65  85 181 160 152
                195 180 139 154   0 154 170 159 174 198 198 202 196 205 169 112 149 215 194 192
                91  43  46  61 154   0  29  87  24 109 113  53 101 116  76  68  42 130  99  96
                107  54  42  45 170  29   0  98  40 134 138  56 126 140  93  80  65 152 122 121
                60 125  80  94 159  87  98   0  98 135 138 127 127 153  42  56  59 184 147 109
                86  29  68  81 174  24  40  98   0  94  99  32  87 100  77  89  47 115  83  84
                94  97 149 168 198 109 134 135  94   0   5 102  10  21  95 142  89  61  33  29
                96 102 153 172 198 113 138 138  99   5   0 107  15  22  98 145  92  61  36  32
                106  26  94 101 202  53  56 127  32 102 107   0  95 102 103 121  78 110  85  97
                84  91 142 160 196 101 126 127  87  10  15  95   0  28  87 135  81  67  36  21
                113  97 158 177 205 116 140 153 100  21  22 102  28   0 114 155 103  40  22  50
                27 103  91 108 169  76  93  42  77  95  98 103  87 114   0  74  38 147 110  68
                99 110  46  65 112  68  80  56  89 142 145 121 135 155  74   0  58 177 144 124
                58  71  65  85 149  42  65  59  47  89  92  78  81 103  38  58   0 128  92  69
                148 101 174 181 215 130 152 184 115  61  61 110  67  40 147 177 128   0  37  88
                112  77 143 160 194  99 122 147  83  33  36  85  36  22 110 144  92  37   0  55
                64  96 133 152 192  96 121 109  84  29  32  97  21  50  68 124  69  88  55   0
                /;

            my @GD_AAs = qw/A R N D C Q E G H I L K M F P S T W Y V/;

            # my %GD_dic;

            my $GD_num = @GD_AAs;

            for ( my $i = 0; $i < $GD_num; $i++ ) {
                for ( my $j = 0; $j < $GD_num; $j++ ) {

                    # $GD_dic->{GD}->{ $GD_AAs[$i] }->{ $GD_AAs[$j] }
                    # = $grantham_matrix[ ( $i * $GD_num ) + $j ];
                    ${dict}{GD}->{ $GD_AAs[$i] }->{ $GD_AAs[$j] }
                        = $grantham_matrix[ ( $i * $GD_num ) + $j ];
                }
            }

            # $annotations{$locus_tag}{nucleotide}   = $nucleotide;
            $gene         = undef;
            $locus_tag    = undef;
            $gene_synonym = undef;
            $protein_id   = undef;
        }
    }

    # $annotations{options}{genome_seq} = $genome_seq;
    $annotations{options}{version}       = $version;
    # $annotations{options}{logo}          = $logo;
    $annotations{options}{reference}     = $display_id;
    $annotations{options}{strain_name}   = $strain_name;
    $annotations{options}{organism_name} = $organism_name;
    print "Genome length: $genome_length bp\n";
    print "Reference id: $display_id\n";

    #-----------------------------------------------
    if ( $snp_list ne "" ) {
        my $loc_snp_count = 0;
        my ( $col1, $col2, $col3, $col4 );
        my $snp_list_name = "";
        my $snp_notation  = "";

        $annotations{options}{snp_list_exists} = 1;
        open( my $fh, "<", $snp_list )
            or die "cannot open < $snp_list: $!";
        while (<$fh>) {
            next if substr( $_, 0, 1 ) eq "#";
            my $str = "$_";

            chomp($str);
            ( $col1, $col2, $col3, $col4 ) = split( "\t", $str );

            # print "$snp_list_name\n$col2\t$col3\n";
            # push(@snp_list_a,$)
            $snp_list_name = $col1 if ( $col1 ne "" );
            next if ( $col2 eq "" );
            if (   $snp_list_name ne ""
                && $col3 ne ""
                && $col4 ne "" )
            {
                $snp_notation = check_snp_notation($col3);

                # print "$snp_notation\t$col3\t$col4\n";
                if ( $snp_notation eq "nc_gene" && $col2 ne "" ) {
                    ${dict}{snp_list}{coding}{$snp_list_name}{ uc($col2) }
                        { uc($col3) }{'note'} = $col4;

                }
                elsif ( $snp_notation eq "nc_genome" ) {
                    ${dict}{snp_list}{coding}{$snp_list_name}{ uc($col3) }{'note'}
                        = $col4;

                }
                elsif ( $snp_notation eq "pos" ) {

                    # my ($pos_val)=$col3=~/(\d+)/;
                    ${dict}{snp_list}{coding}{$snp_list_name}{ uc($col3) }{'note'}
                        = $col4;

                    # print "$snp_list_name\t$col3\t$col4\t$pos_val\n";

                }
                elsif ( $snp_notation eq "codon" ) {

                    # my ($pos_val)=$col3=~/(\d+)/;
                    ${dict}{snp_list}{coding}{$snp_list_name}{ uc($col2) }
                        { uc($col3) }{'note'} = $col4;

                }
                elsif ( $snp_notation eq "igpos" ) {

                    # my ($pos_val)=$col3=~/(\d+)/;                    
                    ${dict}{snp_list}{intergenic}{$snp_list_name}
                        { uc($col3) }{'note'} = $col4;
                         $annotations{options}{igpos_list_exists} = 1;

                }
                elsif ( $snp_notation eq "aa_genome" ) {
                    ${dict}{snp_list}{coding}{$snp_list_name}{ uc($col3) }{'note'}
                        = $col4;

                }
                elsif ( $snp_notation eq "aa_genome_long" ) {
                    ${dict}{snp_list}{coding}{$snp_list_name}{ uc($col3) }{'note'}
                        = $col4;
                }
                elsif ( $snp_notation eq "aa" && $col2 ne "" ) {
                    ${dict}{snp_list}{coding}{$snp_list_name}{ uc($col2) }
                        { uc($col3) }{'note'} = $col4;
                }
                elsif ( $snp_notation eq "aa_long" ) {
                    ${dict}{snp_list}{coding}{$snp_list_name}{ uc($col2) }
                        { uc($col3) }{'note'} = $col4;
                }
                $loc_snp_count++;

            }
            else {
                print "Bad format of $snp_list file!\n";
                $annotations{options}{snp_list_exists} = 0;
                $annotations{options}{igpos_list_exists} = 0;
            }

       # $str =~ /^(\w+)\t(\S+)\t(.*)/;
       # if ( $col1 eq "" ) {
       #     print "Bad format of $snp_list file!\n";
       #     $annotations{options}{snp_list_exists} = 0;
       # }
       # else {
       #     ${dict}{snp_list}{ uc($col2) . "_" . uc($col3) }{'note'} = $col4;
       #     $loc_snp_count++;

            #     # last;
            # }

            # my $loc_pos=$1;
            # my $loc_ref=$2;
            # my $loc_alt=$3;

        }

        close $fh;
        print "SNP's processed: $loc_snp_count\n";

    }
    else {
        $annotations{options}{snp_list_exists} = 0;
        print "SNP's processed: file not set\n";
    }

    #-----------------------------------------------

    # $annotations{options}{codon_usage} = "no";

    my %hoh = ( annotations => \%annotations, dict => \%dict );

    # store \%hoh, 'HoH.sto';

    print "CDS processed: $count_genes\n";
    print "Filename: $out_file\n";
    print "Version: $version\n";

    print "-" x 50 . "\n";

    my $tmp_annotations = freeze( \%hoh );

    # my $tmp_tang=

    my $buffer = "";
    my $fh     = tempfile();
    open $fh, '>', \$buffer;
    print $fh $tmp_annotations;

    # print $fh $tmp_annotations;
    rawdeflate \$buffer => $out_file
        or die "rawdeflate failed: $RawDeflateError\n";

}

# sub codon_usage_db {
#     my $loc_read_file = "";
#     open( FILE, $codon_usage )
#         or die "Couldn't open file $codon_usage, $!"

#         #             ;    #открываем каждый файл
#         while (<FILE>) {    # ищем строки со снипами
#         $loc_read_file =~ /locus_tag="(\w+)".*\n(.*)/;
#         print "$loc_read_file\n";

#     }

#     # body...
# }

# @read_codon_info="";
# if ($codon_file) {
#     print "$codon_file!\n";
#     open( my $fh, "<", $codon_file )
#         or die "cannot open < input.txt: $!";

#     while ( <$fh> ) {
#             # $_ =~ /locus_tag=\W(?<locus>\w+).*\W(?<codon_freq>\d+.*)/;
#             $_ =~ /locus_tag=\W(?<locus>\w+)/;
#             $_ =~ /(?<codon_freq>.*^\d+.*)/m;
#             print "$+{locus}  $+{codon_freq}\n";
#              # push (@read_codon_info, $1 . " " . $2);
#             }

# }

# foreach $val (@read_codon_info){
#     print "$val\n";
# }
sub check_snp_notation {
    my $mut = shift;
    $mut =~ /^(\d+)(\D)>(\D)/x;    # 100A>G
    if ( $1 ne "" && $2 ne "" && $3 ne "" ) {
        return "nc_gene";
    }
    else {
        $mut =~ /^(\D)(\d+)(\D)/x;    # L234R
        if ( $1 ne "" && $2 ne "" && $3 ne "" ) {
            return "aa";
        }
        else {

            $mut =~ /^(\d+)_(\w)>(\w)/x;    # 10000_A>G
            if ( $1 ne "" && $2 ne "" && $3 ne "" ) {
                return "nc_genome";
            }
            else {

                $mut =~ /^(\d+)_(\w{3})>(\w{3})/x;    # 10000_Ser>Leu
                if ( $1 ne "" && $2 ne "" && $3 ne "" ) {
                    return "aa_genome_long";
                }
                else {

                    $mut =~ /^(pos\d+)/x;             # pos10000
                    if ( $1 ne "" ) {
                        return "pos";
                    }
                    else {

                        $mut =~ /^(codon\d+)/x;       # codon100
                        if ( $1 ne "" ) {
                            return "codon";
                        }
                        else {

                            $mut =~ /^igpos(\d+)/x;    # igpos10000 intergenic
                            if ( $1 ne "" ) {
                                return "igpos";
                            }
                            else {
                                $mut =~ /^(\D{3})(\d+)(\D{3})/x;    #Asp749Glu
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
            }

        }
    }
}
