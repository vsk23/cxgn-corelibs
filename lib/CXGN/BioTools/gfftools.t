#! /usr/bin/perl

use Test::More tests=>12;



use strict;
use warnings;
use FileHandle;

use CXGN::BioTools::AceTools;
use CXGN::BioTools::GffTools;



my $tab_dir = "/home/vsk23/cxgn/data/tab_files/Nspecies.tab";
my $gff_data_dir = "/home/vsk23/cxgn/data/gff_files/Nspecies_snpgff_final.gff";#
my $gff_index_dir = "/home/vsk23/cxgn/data/index_files/Nspecies_snpgff_final_index.db";

#TESTING build_gff_index_file#
my $index_hash =CXGN::BioTools::GffTools-> build_index_file($gff_data_dir, $gff_index_dir);
print "size of returned has is". keys(%$index_hash)."\n";

#TESTING get_location_by_id#
my $location = CXGN::BioTools::GffTools-> get_location_by_id("step3_c10005",$index_hash);
print "location of step3_c10005 is : ".$location."\n";

#TESTING get_seq_by_location#
my $contig_info = CXGN::BioTools::GffTools->get_contig_information($index_hash,$location,"step3_c10005",$gff_data_dir);

while ( my ($key, $value) = each(%$contig_info) ) {
              
    if($key=~m/_sc/){
                   print "$key => $value\n";
    }
}

#TESTING INFORMATION IN THE HASH#
is ($contig_info->{"Nsyl_c4415_sc"}, "[1-554]", "Nsyl_c4415 seq coords");
like($contig_info->{"Nsyl_c4415"},qr/TATTTT/, "Nsyl_c4415 sequence");
is ($contig_info->{"Nsyl_c4415_READCOORDS"},"[1-554]","Nsyl_c4415 read coords");
is ($contig_info->{"Nsyl_c4415_SAMPLE"},"tomamtoes", "Nsyl_c4415 sample");
is ($contig_info->{"Nsyl_c4415_PARENT"}, "step3_c10005","Nsyl_c4415 parent");
is ($contig_info->{"Ntab_s96186_sc"},"[1-424]", "Ntab_s96186 sequence coordinates");
is ($contig_info->{"Ntab_s96186_READCOORDS"},"[1-424]", "Ntab_s96186 read coordinates");
is ($contig_info->{"Ntab_s96186_PARENT"},"step3_c10005", "Ntab_s96186 parent");
is ($contig_info->{"Ntab_s96186_SAMPLE"},"tomamtoes","Ntab_c96186 sample");
like($contig_info->{"Ntab_s96186"},qr/ACCGC/,"Ntab_s96186 sequence");
like($contig_info->{"step3_c10005"}, qr/CGTTAA/, "step3_c10005 sequence");
is ($contig_info->{"step3_c10005_sc"},"[1-669]","step3_c10005 coordinates");

