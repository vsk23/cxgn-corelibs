#! /usr/bin/perl

=head1 NAME

CXGN::BioTools:AceTools - Perl module to manipulate ACE FILES.

=head1 SYNOPSIS

    # Module loading
    use CXGN::BioTools::AceTools;

    # Module using methods

        CXGN::BioTools::AceTools -> method_name(args_here);

=head1 DESCRIPTION

CXGN::BioTools::AceTools was made to help parse and process ACE files
in a (hopefully) faster and more efficient way then currently available.
This module seeks to provide misc. tools to the user in order to make
the processing of ACE FILES more convenient. ( Converting, Parsing, etc..)

=cut
=head1 AUTHOR - Vera Kutsenko

vsk23@cornell.edu

=cut

package CXGN::BioTools::AceTools;


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
    Usage   : $tied_hash_ref =CXGN::BioTools::AceTools-> build_index_file($data_dir,$index_dir);
    Function: ties a hash with key(contig_id) -> value (byte location in data file) to a file.Opens up the two 
              handles with one referring to the ace file ($data_handle) and the other to the indexfile with extension
              .db ($index_handle). The $data_handle goes through the ace file line by line; when it reaches a line that starts with "CO"/"RD" it
              obtains its current byte location and stores this as the value for the key(which is the contig_id) in %index_hash which will be 
              returned as a hashref.
    Returns : hash reference
    Args    : -data_dir : directory with filename of the ACE file
              -index_dir: directory of index file / or directory with filename that you want the index to be created with
=cut

sub build_index_file{
    
my ($self,$data_file,$index_file) = @_;
my $contig_id="";
my $offset = 0;
my ($data_handle,$index_handle) = FileHandle->new;
if (-f $index_file){
     tie(my %index_hash, 'DB_File',"$index_file", O_RDONLY, 0777)
         or die "Couldn't tie DB file $index_file: $!; aborting";
     return \%index_hash;
} else {
     $data_handle-> open("< $data_file") or die "Can't open data handle: $!\n";
      tie(my  %index_hash, 'DB_File',"$index_file", O_CREAT|O_WRONLY, 0777)
          or die "Couldn't tie DB file $index_file: $!; aborting";
      while ( <$data_handle>) {
             if ($_ =~ m/^CO\s+(\S+)\s+/){
                 $contig_id = $1;
                 $index_hash{$contig_id}=$offset;
               } 
              $offset = tell($data_handle);
          }
	 return \%index_hash;
}

}


=head2 get_location_by_id ($id, $index_hash_ref)

    Title   : CXGN::BioTools::AceTools->get_location_by_id
    Usage   : $location = get_location_by_id($id, $index_hash_ref);
    Function: Obtains the byte location of the contig_id given by $id in the ace file.
    Returns : a value
    Args    : -id : contig_id whose byte location needs to be returned
              -index_hash_ref: tied hash_ref that contains the contig_id keys with byte loc. values.
=cut
sub get_location_by_id{     
     my ($self, $id,$index_hash_ref) = @_;   
     if(exists $index_hash_ref->{$id}){
          return $index_hash_ref->{$id};	
        } else {
          print "key like that does not exist!"
        }
       
}

=head2 get_contig_information ($byte_location, $data_file, $contig_id)

    Title   : CXGN::BioTools::AceTools->get_contig_information($byte_location, $data_file, $contig_id)
    Usage   : $data_hash = get_contig_information($byte_location, $data_file, $contig_id);
    Function: Obtains all the information ( Read IDs & Sequences .. Contig ID && Sequence.. Coordinates) 
              for the given $contig_id, retrieves and stores in a hash_ref.
    Returns : a $hash_reference ($id_and_seq_data in the method)
              Keys:(Abstract)                 Values:
                   contig_id_identifier            Program that created the assembly data
                   contig_id_length                Length of contig
                   contig_id                       Sequence
                   read_id                         Sequence of read
                   read_id_sc                      Sequence coordinates of this read
                   read_id_rc                      Read coordinates of this read
    Args    : -byte_location   : byte location of the $contig_id in $data_dir
              -dataf_file      : directory & basename of the ACE FILE
              -contig_id       : contig id to be located
=cut
sub get_contig_information{

my ($self, $byte_location, $data_file, $contig_id) = @_;
my ($seq_coord, $read_coord, $testline, $read_id);
my ($sequence, $identifier) = "";
my %id_and_seq_data;
my $stop = "false";
my $contig_length=0;
my $end= "false";   
my $read_length=0;
my $data_handle= FileHandle->new;

$data_handle-> open("< $data_file") or die "Can't open data handle: $!\n";
seek($data_handle, $byte_location,0);
## Get Contig and Read IDs and reads
while($stop eq "false" && defined(fileno($data_handle))){        
       my $line = <$data_handle>;
       if(!eof($data_handle)){
           #Set all identifiers#
	   if($line=~m/$contig_id\s?\S+\s?(\S+)\s?/){
	       $id_and_seq_data{$contig_id."_identifier"} = $1;
             
	   }
           if($line=~m/^CO\s+(\S+)\s+(\S+)\s+/ && tell($data_handle) != -1){
              if (!($1 eq $contig_id)){close ($data_handle);}
	      $id_and_seq_data{$1."_length"} = $2;
              $end="false";             
              $identifier = "CO";
	     }
	   if(defined(fileno($data_handle))){
              if($line=~m/^RD\s+(\S+)\s+(\d+)/){
	          $identifier = "RD";
	          $read_id = $1; 
		  $read_length=$2;
	          $end= "false";
	        }
             # check the appropriate tags 
             # 1. If identifier == CO, then store the lines as value for appropriate hash key
             # 2. If identifier == RD, then store the sequence as a value for read_id as key.
             # 3. If identifier == QA, then store the coordinates as values for appropriate keys.
             if( !($line =~ m/^(AS|BQ|BS|CT{|WA|DS|RT{)\s+/)){
                  if($identifier eq "CO" && !($line =~ m/^QA|CO/) && ($end eq "false")){
                      $id_and_seq_data{$contig_id} .= "".$line; 
                    }
                  if($identifier eq "RD" && !($line =~ m/^(QA|RD)/) && ($end eq "false")){
	              $id_and_seq_data{$read_id} .= "".$line;
                    }            
                  if($line=~m/^QA\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/){
	             $id_and_seq_data{$read_id."_rc"} = "[".$3."-".$4."]";
                     $id_and_seq_data{$read_id."_end"} = $id_and_seq_data{$read_id."_start"}+$read_length-1;
                     $id_and_seq_data{$read_id."_sc"} = "[".$id_and_seq_data{$read_id."_start"}."-".$id_and_seq_data{$read_id."_end"}."]";
                    }
		  if($line=~m/^AF\s?(\S+)\s?\S?\s?(\d*)/){                      
                      my $read_id = $1;
                      $id_and_seq_data{$read_id."_start"} = $2;		  
		  }
             } else {
                $end="true"; #reached end of sequence for specified read/contig
               }
	  }
	   
       }else{
         $stop = "true"; #reached end of parsing all NEEDED information
       }
}        
return \%id_and_seq_data;

}

=head2 convert_ace_to_gff

    Title   : convert_ace_to_gff
    Usage   : CXGN::BioTools::AceTools->convert_ace_to_gff($file_in, $file_out, $ace_index_file, $tab_file, $index_hash_ref);
    Function: Efficiently converts(as opposed to other converter found online) ACE file format ($file_in) 
              into a GFF format ($gff_file_out). Output includes the
              following tags: Contig_ID/Read_ID, program with which this was created,begin & end ends, PARENT ID for reads,
              READ COORDS, SAMPLE....
    Returns : NOTHING . Writes to a .gff file.
    Args    : -file_in           : ACE FILE to be converted
              -gff_file_out      : GFF FILE to be written to -> does not have to exist, will be created otherwise
              -index_file        : indexing file. (If non-existant, create using build_index_hash)
              -tab_file          : tab file with 1st column: read ID, and 2nd column: sample
              -index_hash_ref    : hash ref with the byte locations (If non-existant, create using get_contig_info)
=cut

sub convert_ace_to_gff{
    
my($self, $file_in, $gff_name, $index_file, $tab_file,$index_hash,$tab_index_file) = @_;
my ($gff_annotation_handle,$gff_fasta_handle,$tab_handle) = FileHandle->new;
my %tab_hash;
my $fhline;
if(!(-e $gff_name."_final.gff")){
    print "GFF FILE DOES NOT EXIST..CONVERTING...\n";
## TIE TAB TO A HASH ##
if(-f $tab_index_file){
    open($tab_handle, ">$tab_file") or die $!;
       tie (%tab_hash, 'DB_File', $tab_index_file, O_RDONLY,0666) or die "Couldn't tie existing hash";
}else{
     open($tab_handle,"<$tab_file") or die $!;
     tie(%tab_hash, 'DB_File',$tab_index_file, O_CREAT|O_WRONLY, 0666)
            or die "Couldn't tie DB file $index_file: $!; aborting";
       while(<$tab_handle>){
             if ($_ =~ m/(\S+)\t?\s?(\S+)/){
                 $tab_hash{$1} = $2;
               }
       }
}

## CONVERT ACE TO GFF ##
while ( my ($contig_id, $byte) = each(%$index_hash) ) {
        my $data_hash =CXGN::BioTools::AceTools-> get_contig_information($byte,$file_in,$contig_id);
        
         ## PRINT CONTIG INFO TO GFF FILE ##
	open($gff_fasta_handle, ">>$gff_name"."_fas.gff") or die $!; 
        open ($gff_annotation_handle, ">>$gff_name"."_ann.gff") or die $!; # create if not there

	 my $gff_line = "$contig_id\t".$data_hash->{$contig_id."_identifier"}."\tcontig\t1\t".$data_hash->{$contig_id."_length"}."\t.\t+\t.\tID=$contig_id\n";

	 print $gff_annotation_handle "$gff_line";
         print $gff_fasta_handle ">$contig_id \n".$data_hash->{$contig_id}."\n"; # print to second part of gff file in FASTA format
         delete $data_hash->{$contig_id."_length"};
         delete $data_hash->{$contig_id};
         
         ## PRINT READS OF PREVIOUS CONTIG TO GFF FILE ##
while ( my ($key, $value) = each(%$data_hash) ){
         if($key =~ m/(\S+)_sc/){
	     my $read_id = $1;
             print $gff_fasta_handle ">$read_id\n".$data_hash->{$read_id}."\n"; # append read seq to 2nd part of gff FASTA format
	     my $read_sc = $data_hash->{$read_id."_sc"};
             if($read_sc =~ m/\S?(\d+)\-(\d+)\S?/){
	     $fhline= "$read_id\t".$data_hash->{$contig_id."_identifier"}."\tread\t".$1."\t".$2."\t.\t+\t.\tID=$read_id;PARENT=$contig_id;READCOORDS=".$data_hash->{$read_id."_rc"}.";";
	     }
             if($tab_hash{$read_id}){
                 print $gff_annotation_handle $fhline."SAMPLE=".$tab_hash{$read_id}.";\n";
	        } else { 
                 print $gff_annotation_handle $fhline."\n";
                }
	 }
 
 }
     close ($gff_annotation_handle);
     close ($gff_fasta_handle);
     close ($tab_handle);
 }
## CONCATENATE THE TWO FILES TOGETHER INTO THE FINAL GFF FILE ##
my $concatenation_command = "cat $gff_name"."_ann.gff"." $gff_name"."_fas.gff"." > "."$gff_name"."_final.gff";
system($concatenation_command);
unlink("$gff_name"."_ann.gff");
unlink("$gff_name"."_fas.gff");
}else{
    print "GFF FILE EXISTS. NO CONVERSION OCCURED";
}
}

1;
