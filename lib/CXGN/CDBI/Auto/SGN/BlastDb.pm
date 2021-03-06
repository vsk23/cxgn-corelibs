package CXGN::CDBI::Auto::SGN::BlastDb;
# This class is autogenerated by cdbigen.pl.  Any modification
# by you will be fruitless.

=head1 DESCRIPTION

CXGN::CDBI::Auto::SGN::BlastDb - object abstraction for rows in the sgn.blast_db table.

Autogenerated by cdbigen.pl.

=head1 DATA FIELDS

  Primary Keys:
      blast_db_id

  Columns:
      blast_db_id
      file_base
      title
      type
      source_url
      lookup_url
      update_freq
      info_url
      index_seqs
      blast_db_group_id
      web_interface_visible
      description

  Sequence:
      none

=cut

use base 'CXGN::CDBI::Class::DBI';
__PACKAGE__->table( 'sgn.blast_db' );

our @primary_key_names =
    qw/
      blast_db_id
      /;

our @column_names =
    qw/
      blast_db_id
      file_base
      title
      type
      source_url
      lookup_url
      update_freq
      info_url
      index_seqs
      blast_db_group_id
      web_interface_visible
      description
      /;

__PACKAGE__->columns( Primary => @primary_key_names, );
__PACKAGE__->columns( All     => @column_names,      );


=head1 AUTHOR

cdbigen.pl

=cut

###
1;#do not remove
###
