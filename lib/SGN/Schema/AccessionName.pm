package SGN::Schema::AccessionName;

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("Core");
__PACKAGE__->table("accession_names");
__PACKAGE__->add_columns(
  "accession_name_id",
  {
    data_type => "bigint",
    default_value => "nextval('accession_names_accession_name_id_seq'::regclass)",
    is_auto_increment => 1,
    is_nullable => 0,
    size => 8,
  },
  "accession_name",
  {
    data_type => "character varying",
    default_value => undef,
    is_nullable => 1,
    size => 255,
  },
  "accession_id",
  {
    data_type => "bigint",
    default_value => undef,
    is_foreign_key => 1,
    is_nullable => 1,
    size => 8,
  },
);
__PACKAGE__->set_primary_key("accession_name_id");
__PACKAGE__->might_have(
  "accession",
  "SGN::Schema::Accession",
  { "foreign.accession_name_id" => "self.accession_name_id" },
);
__PACKAGE__->belongs_to(
  "accession",
  "SGN::Schema::Accession",
  { accession_id => "accession_id" },
  { join_type => "LEFT" },
);


# Created by DBIx::Class::Schema::Loader v0.04999_07 @ 2009-09-04 13:21:55
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:G9e5N3E4Irk71oOf99ETXg


# You can replace this text with custom content, and it will be preserved on regeneration
1;
