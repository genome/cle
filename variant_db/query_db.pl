#! /gsc/bin/perl

use strict;
use warnings;

use DBI;

unless (@ARGV == 1) {
    die "Provide mutation type(somatic or germline)";
}

my $type = $ARGV[0];
my $tablename = $type.'_variant';
my $dbfile = '/gscmnt/gc3042/cle_validation/CLE_variant_database/sqlite_variant_DB/cle_variants.sqlite';

my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', { AutoCommit => 0, RaiseError => 1 })
    or main->fatal_msg("Can't connect to db ($dbfile): " . $DBI::errstr);

my @columns = qw(
    assay_version
    pipeline_version
    date
    process_id
    normal_sample
    chromosome
    start
    stop
    reference
    variant
    variant_type
    transcript_name
    trv_type
    trv_type_category
    c_position
    amino_acid_change
    default_gene_name
    removed_manual_review
);
    
if ($type eq 'somatic') {
    splice (@columns, 5, 0, 'tumor_sample');
}

my $column_str = join ',', @columns;

my $sth = $dbh->prepare(
    qq(
        select $column_str
        from $tablename
    )
);
$sth->execute;
my $outs = $sth->fetchall_arrayref;

if ($outs) {
    for my $out (@$outs) {
        print join "\t", @$out, "\n";
    }
}
else {
    die "Failed to query sqlite DB $dbfile\n";
}
