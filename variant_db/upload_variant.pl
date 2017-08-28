#! /usr/bin/perl

use strict;
use warnings;

use above 'Genome';
use DBI;

unless (@ARGV == 3) {
    die "Provide three arguments in order: mutation type(somatic or germline), process id, upload file path";
}

my ($type, $process_id, $file) = @ARGV;

die "germline mode is not supported now" if $type eq 'germline';

my $dbfile = '/gscmnt/gc3042/cle_validation/sqlite_variant_DB/cle_variants.sqlite';

my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', { AutoCommit => 0, RaiseError => 1 })
    or main->fatal_msg("Can't connect to db ($dbfile): " . $DBI::errstr);

#create the database
unless (-s $dbfile) {
    create_table($dbh, 'somatic');
    create_table($dbh, 'germline');
}

my $fh = Genome::Sys->open_file_for_reading($file);
my $tablename = $type.'_variant';

while (my $line = $fh->getline) {
    chomp $line;
    next if $line =~ /^Man/;
    my @columns = split /\t/, $line;
    my $sth;

    if ($type eq 'somatic') {
        die "Expect 62 columns from $type upload file" unless @columns == 62;
        $sth = $dbh->prepare("insert into $tablename (assay_version, pipeline_version, date, process_id, normal_sample, tumor_sample, chromosome, start, stop, reference, variant, variant_type, transcript_name, trv_type, trv_type_category, c_position, amino_acid_change, default_gene_name, ensembl_gene_id, inSegDup, AML_RMG, rsid, dbSNP_caf, dbSNP_max_alt_af, onTarget, MeetsMinDepthCutoff, min_coverage_observed, max_normal_vaf_observed, max_tumor_vaf_observed, variant_callers, variant_caller_count, removed_manual_review) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)") or die $DBI::errstr;
    }
    elsif ($type eq 'germline') {
        $sth = $dbh->prepare("insert into $tablename (assay_version, pipeline_version, date, process_id, normal_sample, chromosome, start, stop, reference, variant, variant_type, transcript_name, trv_type, trv_type_category, c_position, amino_acid_change, default_gene_name, ensembl_gene_id, inSegDup, AML_RMG, rsid, dbSNP_caf, dbSNP_max_alt_af, onTarget, MeetsMinDepthCutoff, min_coverage_observed, NHLBI_All_MAF, NHLBI_AA_MAF, NHLBI_EU_MAF, removed_manual_review) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)") or die $DBI::errstr;
    }
    my @info = get_info(@columns);
    $sth->execute_array({}, @info) or main->fatal_msg($DBI::errstr);
    $dbh->commit;
}

$fh->close;



sub create_table {
    my ($dbh, $type) = @_;
    my $table_name = $type.'_variant';
    my $sth;

    if ($type eq 'somatic') {
        $sth = $dbh->prepare("create table $table_name (id integer primary key, assay_version text, pipeline_version text, date text, process_id text, normal_sample text, tumor_sample text, chromosome text, start integer, stop integer, reference text, variant text, variant_type text, transcript_name text, trv_type text, trv_type_category text, c_position text, amino_acid_change text, default_gene_name text, ensembl_gene_id text, inSegDup boolean, AML_RMG boolean, rsid text, dbSNP_caf text, dbSNP_max_alt_af text, onTarget boolean, MeetsMinDepthCutoff boolean, min_coverage_observed integer, max_normal_vaf_observed float, max_tumor_vaf_observed float, variant_callers text, variant_caller_count integer, removed_manual_review boolean)");
    }
    elsif ($type eq 'germline') {
        $sth = $dbh->prepare("create table $table_name (id integer primary key, assay_version text, pipeline_version text, date text, process_id text, normal_sample text, chromosome text, start integer, stop integer, reference text, variant text, variant_type text, transcript_name text, trv_type text, trv_type_category text, c_position text, amino_acid_change text, default_gene_name text, ensembl_gene_id text, inSegDup boolean, AML_RMG boolean, rsid text, dbSNP_caf text, dbSNP_max_alt_af text, onTarget boolean, MeetsMinDepthCutoff boolean, min_coverage_observed integer, NHLBI_All_MAF numeric, NHLBI_AA_MAF numeric, NHLBI_EU_MAF numeric, removed_manual_review boolean)");
    }
    else {
        die "Type must be either somatic or germline";
    }
    $sth->execute or main->fatal_msg($DBI::errstr);
    $dbh->commit;

    #$sth = $dbh->prepare("create unique index variant on $table_name(id)");
    #$sth->execute or main->fatal_msg($DBI::errstr);
    #$dbh->commit;
}

sub get_info {
    my @cols = @_;

    my $process = Genome::Process->get($process_id);
    die "Failed to get process $process_id" unless $process;

    my $assay_version = 'NGv3 + RMG-1';
    my $version = get_pipeline_version($process->software_revision);
    my $date = $process->ended_at;
    my $normal_sample = $process->normal_sample->name;
    my $tumor_sample  = $process->tumor_sample->name;

    my @info;
    #somatic 32 columns, germline 30 columns
    if ($type eq 'somatic') {
        @info = ($assay_version, $version, $date, $process_id, $normal_sample, $tumor_sample);
        push @info, map{$cols[$_]}(1..20, 30..34);
        my $review = $cols[0] eq 'x' ? 1 : 0;
        push @info, $review;
    }
    return @info;
}


sub get_pipeline_version {
    my $revision = shift;
    my ($snap_shot_version) = $revision =~ /snapshots\/genome-(\d+)\//;

    if ($snap_shot_version >= 3714 and $snap_shot_version < 3722) {
        return '1.3';
    }
    elsif ($snap_shot_version >= 3722 and $snap_shot_version < 3734) {
        return '1.4';
    }
    elsif ($snap_shot_version >= 3734) {
        return '1.5';
    }
    else {
        die 'Failed to get pipeline version for snapshot_version genome-'.$snap_shot_version;
    }
}
