package CXGN::CDBI::Auto::SGN::MarkerDerivedFrom;
# This class is autogenerated by cdbigen.pl.  Any modification
# by you will be fruitless.

=head1 DESCRIPTION

CXGN::CDBI::Auto::SGN::MarkerDerivedFrom - object abstraction for rows in the sgn.marker_derived_from table.

Autogenerated by cdbigen.pl.

=head1 DATA FIELDS

  Primary Keys:
      marker_derived_dummy_id

  Columns:
      marker_derived_dummy_id
      marker_id
      derived_from_source_id
      id_in_source

  Sequence:
      none

=cut

use base 'CXGN::CDBI::Class::DBI';
__PACKAGE__->table( 'sgn.marker_derived_from' );

our @primary_key_names =
    qw/
      marker_derived_dummy_id
      /;

our @column_names =
    qw/
      marker_derived_dummy_id
      marker_id
      derived_from_source_id
      id_in_source
      /;

__PACKAGE__->columns( Primary => @primary_key_names, );
__PACKAGE__->columns( All     => @column_names,      );


=head1 AUTHOR

cdbigen.pl

=cut

###
1;#do not remove
###
