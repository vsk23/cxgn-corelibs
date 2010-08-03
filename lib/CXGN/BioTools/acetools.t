#! /usr/bin/perl

use Test::More tests=>12;



use strict;
use warnings;


use CXGN::BioTools::AceTools;




use FileHandle;


## Directories ## - have them preset to something?
my $data_dir ="/home/vsk23/cxgn/data/ace_files/Nspecies_snp.ace";  # ACE FILE
my $index_dir="/home/vsk23/cxgn/data/index_files/Nspecies_snp_index.db";  #DO i need this?
my $contig_id = "step3_c99";   # test contig id
my $contig_id2 = "step3_c10"; # test2 contig id
my $index_dir2= "/home/vsk23/cxgn/data/ace_files/hashtables/Nspecies_snpgff_indx.db";
my $new_gff_dir = "/home/vsk23/cxgn/data/gff_files/Nspecies_snpgff";
my $tab_dir = "/home/vsk23/cxgn/data/tab_files/Nspecies.tab";


## TESTING build_ace_hash_file ##
my $tester_hash =CXGN::BioTools::AceTools-> build_index_file($data_dir, $index_dir);

## TESTING convert_ace_to_gff ##
CXGN::BioTools::AceTools->convert_ace_to_gff($data_dir,$new_gff_dir, $index_dir, $tab_dir,$tester_hash,"/home/vsk23/cxgn/data/tabindex.db");

## TESTING get_location_by_id ###
my $byte_loc =CXGN::BioTools::AceTools-> get_location_by_id($contig_id,$tester_hash);
print "$byte_loc \n";
my $byte_loc2 = CXGN::BioTools::AceTools-> get_location_by_id($contig_id2, $tester_hash);


## TESTING get_seq_by_location ##

my $contig_info = CXGN::BioTools::AceTools->get_contig_information($byte_loc, $data_dir, $contig_id);
my $contig_info2 = CXGN::BioTools::AceTools->get_contig_information($byte_loc2, $data_dir, $contig_id2);


my $test_fh = FileHandle->new;

open($test_fh, ">> /home/vsk23/cxgn/data/test_info.txt");

while (my ($key,$value) = each (%$contig_info)){

    print $test_fh "Key: $key => Value: $value \n";
}

close($test_fh);

is ($contig_info->{"Ntab_s99612_sc"}, "[2280-2492]", "Ntor_s99612 sequence coordinates");
like ($contig_info->{"Ntab_s99612"},qr/GACAAACATG/, "Ntor_c7271 sequence");
is ($contig_info->{"Ntor_c12034_sc"},"[3219-3899]","Ntor_c12034 sequence coordinates");
like($contig_info->{"Ntor_c12034"},qr/CTGATATTGAT/, "Ntor_c20342 sequence");
is ($contig_info->{"Ntor_c6490_sc"}, "[2109-2627]","Nsyl_c501 sequence coordinates");
like ($contig_info->{"Ntor_c6490"},qr/GTAACAAAGTG/, "Nsyl_c501 sequence");


is($contig_info2->{"Ntab_c1965_rc"},"[1-3293]","Ntab_c1965 read coordinates");
like($contig_info2->{"Ntab_c1965"},qr/TATAG/,"Ntab_c1965 sequence");
is($contig_info2->{"Nsyl_c1114_rc"},"[1-3039]","Nsyl_c1114 read coordinates");
like($contig_info2->{"Nsyl_c1114"},qr/GCAGT/,"Nsyl_c1114 sequence");
is($contig_info2->{"Ntor_c20377_rc"},"[1-1276]","Ntor_c203777 read coordinates");
like($contig_info2->{"Ntor_c20377"},qr/TGAA*/,"Ntor_c20377 sequence");


