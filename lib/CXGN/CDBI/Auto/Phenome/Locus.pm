package CXGN::CDBI::Auto::Phenome::Locus;
# This class is autogenerated by cdbigen.pl.  Any modification
# by you will be fruitless.

=head1 DESCRIPTION

CXGN::CDBI::Auto::Phenome::Locus - object abstraction for rows in the phenome.locus table.

Autogenerated by cdbigen.pl.

=head1 DATA FIELDS

  Primary Keys:
      locus_id

  Columns:
      locus_id
      locus_name
      locus_symbol
      original_symbol
      gene_activity
      locus_notes
      obsolete
      sp_person_id
      create_date
      modified_date
      description
      linkage_group
      lg_arm
      common_name_id
      updated_by

  Sequence:
      none

=cut

use base 'CXGN::CDBI::Class::DBI';
__PACKAGE__->table( 'phenome.locus' );

our @primary_key_names =
    qw/
      locus_id
      /;

our @column_names =
    qw/
      locus_id
      locus_name
      locus_symbol
      original_symbol
      gene_activity
      locus_notes
      obsolete
      sp_person_id
      create_date
      modified_date
      description
      linkage_group
      lg_arm
      common_name_id
      updated_by
      /;

__PACKAGE__->columns( Primary => @primary_key_names, );
__PACKAGE__->columns( All     => @column_names,      );


=head1 AUTHOR

cdbigen.pl

=cut

###
1;#do not remove
###
