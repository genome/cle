#! /gsc/bin/perl

use strict;
use warnings;

use DBI;
use File::Spec;

#For cle varaint DB, 1 db file, three tables

my @valid_names = qw(somatic_variant idt_somatic_variant myeloseq_variant myeloseq_variant_count);
my $valid_names = join ",", @valid_names;

die "Provide valid table name" unless @ARGV and @ARGV == 1;

my $tablename = $ARGV[0];
die "table_name must be one of $valid_names" unless grep{$_ eq $tablename}@valid_names;

my $db_file = '/gscmnt/gc3042/cle_validation/CLE_variant_database/sqlite_variant_DB/cle_variants.sqlite';

unless (-s $db_file) {
    die "db_file $db_file is not valid.";
}

my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", '', '', { AutoCommit => 0, RaiseError => 1 })
    or main->fatal_msg("Can't connect to db ($db_file): " . $DBI::errstr);

my $sth = $dbh->prepare(
    qq(
        select *
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
    die "Failed to query sqlite DB $db_file for $tablename\n";
}
