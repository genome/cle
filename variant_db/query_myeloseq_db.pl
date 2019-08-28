#! /gsc/bin/perl

use strict;
use warnings;

use DBI;

my $tablename = 'myeloseq_variant';
my $dbfile = '/gscmnt/gc3042/cle_validation/CLE_variant_database/sqlite_variant_DB/cle_myeloseq_variants.sqlite';

my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', { AutoCommit => 0, RaiseError => 1 })
    or main->fatal_msg("Can't connect to db ($dbfile): " . $DBI::errstr);

my @columns = qw(
    assay_version
    pipeline_version
    date
    sample
    chromosome
    position
    reference
    variant
    variant_callers
    TAMP
    SAMP
    CVAF
    NR
    NV
    consequence
    gene_symbol
    exon
    intron
    feature_type
    feature
    HGVSc
    HGVsp
    HGNC_ID
    MAX_AF
    MYELOSEQ_TCGA_AC
    MYELOSEQ_MDS_AC
);
    
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
