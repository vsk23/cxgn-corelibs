package CXGN::CDBI::Auto::SGN::TempMapCorrespondence;
# This class is autogenerated by cdbigen.pl.  Any modification
# by you will be fruitless.

=head1 DESCRIPTION

CXGN::CDBI::Auto::SGN::TempMapCorrespondence - object abstraction for rows in the sgn.temp_map_correspondence table.

Autogenerated by cdbigen.pl.

=head1 DATA FIELDS

  Primary Keys:
      tmc_id

  Columns:
      tmc_id
      old_map_id
      map_version_id

  Sequence:
      none

=cut

use base 'CXGN::CDBI::Class::DBI';
__PACKAGE__->table( 'sgn.temp_map_correspondence' );

our @primary_key_names =
    qw/
      tmc_id
      /;

our @column_names =
    qw/
      tmc_id
      old_map_id
      map_version_id
      /;

__PACKAGE__->columns( Primary => @primary_key_names, );
__PACKAGE__->columns( All     => @column_names,      );


=head1 AUTHOR

cdbigen.pl

=cut

###
1;#do not remove
###
