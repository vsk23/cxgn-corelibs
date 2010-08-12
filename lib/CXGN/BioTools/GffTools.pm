#! /usr/bin/perl

=head1 NAME

CXGN::BioTools::GffTools - Perl module to manipulate Gff FILES.

=head1 SYNOPSIS

    # Module loading
    use CXGN::BioTools::GffTools;

    # Module using methods

        CXGN::BioTools::GffTools -> method_name(args_here);

=head1 DESCRIPTION

CXGN::BioTools::GffTools was made to help parse and process GFF files
in a (hopefully) faster and more efficient way then currently available.
This module seeks to provide misc. tools to the user in order to make
the processing of GFF  FILES more convenient. ( Parsing, etc..)

=cut
=head1 AUTHOR - Vera Kutsenko

vsk23@cornell.edu

=cut


package CXGN::BioTools::GffTools;


require Exporter;
         @ISA = qw(Exporter);
         @EXPORT_OK = qw(munge frobnicate);  # symbols to export on request


use strict;
use warnings;

use FileHandle;
use Fcntl;   # For O_RDWR, O_CREAT, etc.
use DB_File;


=head2 build_index_file

    Title   : build_index_file
    Usage   : $tied_hash_ref = CXGN::BioTools::GffTools->build_index_file($data_dir,$index_dir);
    Function: ties a hash with key(contig_id) -> value (byte location in data gff file) to a file.Opens up the two 
              handles with one referring to the gff file ($data_handle) and the other to the indexfile with extension
              .db ($index_handle). The $data_handle goes through the gff file line by line and stores all the contigs/reads
              and their byte locations in the file.
    Returns : hash reference
    Args    : -data_file : file path of the GFF file
              -index_file: file path of index file or file path that you want the index to be created with
=cut
sub build_index_file{
my ($self,$data_file,$index_file) = @_;
my ($data_handle, $index_handle) = FileHandle->new;
my  $contig_id="";
my  $offset= 0;
 
if (-f $index_file){
     tie(my %index_hash, 'DB_File',"$index_file",O_RDONLY, 0666)
         or die "Couldn't tie DB file $index_file: $!; Check permissions or existance! Aborting";
       return \%index_hash;
} else {
   $data_handle-> open("< $data_file") or die "Can't open data handle: $!\n";
   tie(my  %index_hash, 'DB_File',"$index_file", O_CREAT|O_WRONLY, 0666)
        or die "Couldn't tie DB file $index_file: $!; Check permissions or existance! Aborting";
   while ( <$data_handle>) {
        ## STORE CONTIG ANNOTATION LOCATION ##
         if ($_ =~ m/^(\S+)\t?(\S+)\t?contig?\t?/){
             $index_hash{$1}=$offset;
           } 
          ## STORE CONTIG/READ SEQUENCE LOCATIONS ##
	 if($_ =~ m/^>(\S+)/){
             $index_hash{$1."_seq"} = $offset;
           }
       $offset = tell($data_handle);
   }
return \%index_hash;
}
}

=head2 get_location_by_id

    Title   : get_location_by_id
    Usage   : $location = CXGN::BioTools::GffTools->get_location_by_id($id, $index_hash_ref);
    Function: Obtains the byte location of the contig_id given by $id in the gff file.
    Returns : a value
    Args    : -id : contig_id whose byte location needs to be returned
              -index_hash_ref: tied hash_ref that contains the contig_id keys with byte loc. values.
=cut
sub get_location_by_id{
my ($self,$id,$index_hash_ref) = @_;
if(exists $index_hash_ref->{$id}){
   return $index_hash_ref->{$id};	
} else {
print "key like that does not exist!"
 }
}

=head2 get_contig_information

    Title   : get_contig_information
    Usage   : $data_hash = CXGN::BioTools::GffTools->get_contig_information($index_hash,$byte_location, $contig_id, $data_dir);
    Function: Obtains all the information ( Read IDs & Sequences .. Contig ID && Sequence.. Coordinates) 
              for the given $contig_id, retrieves and stores in a hash_ref. - ONLY ONE CONTIG - NOT ALL CONTIGS IN A GIVEN ACE FILE -
             Format of everything is stored in hash:
             - contigid_sc   : sequence coords of contig.
             - contigid      : sequence of contig
             - readid        : sequence of readid
             - readid_sc     : sequence coordinates of read_id if avail.
            GROUP TAGS: tag=value...
             - readid_tag     : value.
                          
                      e.g.if tag=value was PARENT=parenthere then key=value
                             readid_PARENT = parenthere;
                          if tag=value was SAMPLE=sample then key=value
                             readid_SAMPLE = sample; 
    Returns : a $hash_reference
    Args    : -byte_location : byte location of the $contig_id in $data_dir
              -index_hash_ref    : hash that contains byte location of contigs/reads
              -data_file         : filename of the GFF FILE
              -contig_id         : contig id to be located
=cut

sub get_contig_information{
my ($self,$index_hash_ref,$location,$contig_id,$data_file) = @_;
my $data_handle= FileHandle->new;
my $print_handle = FileHandle->new;
my %output_hash;

$data_handle-> open("< $data_file") or die "Can't open data handle: $!\n";
seek($data_handle, $location,0);

## obtain all the reads whose parent is this contig id
while(defined(fileno($data_handle))){
     my $line =(<$data_handle>); 
   #  print"$line \n";	
   #  print "$contig_id \n";
     if(!(eof($data_handle))) { # if not the end of the file/undef. line
        if($line =~ m/$contig_id\t?(\S+)\t?contig\t?(\d?)\t?(\d+)\t?/) { #if contig line
            $location = tell($data_handle);
	    $output_hash{$contig_id."_sc"} = "[".$2."-".$3."]";
            seek($data_handle, $index_hash_ref->{$contig_id."_seq"},0);
            CXGN::BioTools::GffTools->print_sequence(\%output_hash,$data_handle, $contig_id, $location);
          }        
	elsif($line=~m/$contig_id?/){
	#    print " contig id matched \n";
            my @line = split('\t',$line);
          # if($line=~m/^(\S+)\t?(\S+)\t?(\S+)\t?(\d?)\t?(\d+)\t?\.\t?\+\t?\.\t?(\S*)/){ # if read
	 #  print " found read annotation\n";
	     $location = tell($data_handle);
	     my $read_id = $line[0];
             $output_hash{$read_id."_sc"} = "[".$line[3]."-".$line[4]."]";
               #GROUP: contains information in tag=value format#
             my @group = split(';',$line[8]);
             print "GROUP INFO: $line[8] \n";
	     foreach my $info(@group){
		 my @split = split('=',$info);
                 $output_hash{$read_id."_".$split[0]}= $split[1];
	     } 
            print $index_hash_ref->{$read_id."_seq"}."\n";           
             seek($data_handle,$index_hash_ref->{$read_id."_seq"},0);
	   CXGN::BioTools::GffTools-> print_sequence(\%output_hash,$data_handle, $read_id, $location);#}
 
      }else{ # reached end of information for this contig.
	close($data_handle);
      }
     }
}
my $seq_location = CXGN::BioTools::GffTools->get_location_by_id($contig_id."_seq", $index_hash_ref);
return \%output_hash;
}


=head2 print_sequence

    Title   : print_sequence
    Usage   : $data_hash = CXGN::BioTools::GffTools->get_seq_by_location($data_handle, $contig_id, $print_handle, $location);
    Function: Obtains all the information ( Read IDs & Sequences .. Contig ID && Sequence.. Coordinates) 
              for the given $contig_id, retrieves and stores in a hash_ref.
    Returns : a $hash_reference
    Args    : -data_handle   : byte location of the $contig_id in $data_dir
              -print_handle  : directory & basename of the ACE FILE
              -contig_id     : contig id to be located
              -location      :
=cut
sub print_sequence{ 
    my ($self,$output_hash,$data_handle, $contig_id, $location)=@_;
    my $end="false";
    my $line;
    my $countwo=0;
    my $seq="";

   while($end eq "false")  {
         $line = <$data_handle>;
         if($line){
            if(($line =~ m/^>/)){
                  $countwo++;
                 }
            else {
		$seq .= $line;
                 } 
            if($countwo >1) { 
               $end = "true";
               $countwo=0;
               seek($data_handle, $location,0);
	    }
	 }

   }
    $output_hash->{$contig_id} = $seq;
}
