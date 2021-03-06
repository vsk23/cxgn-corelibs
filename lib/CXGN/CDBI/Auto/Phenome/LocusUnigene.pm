package CXGN::CDBI::Auto::Phenome::LocusUnigene;
# This class is autogenerated by cdbigen.pl.  Any modification
# by you will be fruitless.

=head1 DESCRIPTION

CXGN::CDBI::Auto::Phenome::LocusUnigene - object abstraction for rows in the phenome.locus_unigene table.

Autogenerated by cdbigen.pl.

=head1 DATA FIELDS

  Primary Keys:
      locus_unigene_id

  Columns:
      locus_unigene_id
      locus_id
      unigene_id
      obsolete
      sp_person_id
      create_date
      modified_date

  Sequence:
      none

=cut

use base 'CXGN::CDBI::Class::DBI';
__PACKAGE__->table( 'phenome.locus_unigene' );

our @primary_key_names =
    qw/
      locus_unigene_id
      /;

our @column_names =
    qw/
      locus_unigene_id
      locus_id
      unigene_id
      obsolete
      sp_person_id
      create_date
      modified_date
      /;

__PACKAGE__->columns( Primary => @primary_key_names, );
__PACKAGE__->columns( All     => @column_names,      );


=head1 AUTHOR

cdbigen.pl

=cut

###
1;#do not remove
###
