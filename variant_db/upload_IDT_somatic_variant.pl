#! /usr/bin/perl

use strict;
use warnings;

use DBI;
use POSIX;
use Genome;
use File::Basename;

my @cwl_dirs = @ARGV;

my $dir = '/gscmnt/gc13015/cle/IDT_somatic_exome_assay';
my $dbfile = '/gscmnt/gc3042/cle_validation/CLE_variant_database/sqlite_variant_DB/cle_variants.sqlite';

my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', { AutoCommit => 0, RaiseError => 1 })
    or main->fatal_msg("Can't connect to db ($dbfile): " . $DBI::errstr);

my $assay_version = '2.0';
my $pipeline_version = '1.5';
my $tablename = 'idt_somatic_variant';
#create the database
unless (check_table()) {
    create_table();
}

unless (@cwl_dirs) {
    print "Go through the whole IDT somatic cwl dirs\n";
    @cwl_dirs = glob($dir."/CI-*");
}

for my $cwl_dir (@cwl_dirs) {
    my @case_dirs = glob($cwl_dir."/H_*");

    for my $case_dir (@case_dirs) {
        next unless -d $case_dir;
        my $case_name = basename $case_dir;

        my $info = File::Spec->join($cwl_dir, $case_name.'.info');
        unless (-s $info) {
            die "There is no valid $info for $case_name";
        }

        my $n_str = `grep normal $info`;
        my $t_str = `grep tumor $info`;

        my ($normal_sample) = $n_str =~ /^\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/;
        my ($tumor_sample)  = $t_str =~ /^\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/;

        unless ($normal_sample and $tumor_sample) {
            die "Either normal_sample $normal_sample or tumor_sample $tumor_sample is not valid";
        }

        my @tsv = glob($case_dir."/variants.annotated.tsv*");
        unless (@tsv and @tsv == 1) {
            warn "There are none or more than 1 variants.annotated.tsv for $case_name";
            next;
        }
        my $date = POSIX::strftime("%Y-%m-%d %H:%M:%S", localtime((stat $tsv[0])[9]));
        my $fh = Genome::Sys->open_file_for_reading($tsv[0]);
        
        while (my $line = $fh->getline) {
            chomp $line;
            next if $line =~ /^(CHROM|This|#This)/;
            my @columns = split /\t/, $line;
            
            my $col_ct = scalar @columns;
            die "Expect 28 columns from myeloseq upload file" unless $col_ct =~ /^(28)$/;

            $columns[25] =~ s/HGNC://;

            my @info = ($assay_version, $pipeline_version, $date, $normal_sample, $tumor_sample, @columns);
            my $sth = $dbh->prepare("insert into $tablename (assay_version, pipeline_version, date, normal_sample, tumor_sample, chromosome, position, rs_id, reference, variant, variant_callers, normal_gt, normal_ad, normal_af, normal_dp, tumor_gt, tumor_ad, tumor_af, tumor_dp, consequence, symbol, feature_type, feature, HGVSc, HGVSp, cDNA_position, cds_position, protein_position, amino_acids, codons, HGNC_ID, existing_variation, gnomADe_af) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)") or die $DBI::errstr;
    
            $sth->execute_array({}, @info) or main->fatal_msg($DBI::errstr);
            $dbh->commit;
        }
        $fh->close;
    }
    print "$cwl_dir is done\n";
}


sub check_table {
    my $sth = $dbh->prepare("select name from sqlite_master where type='table' and name='$tablename'");
    $sth->execute or main->fatal_msg($DBI::errstr);
    my @tables = $sth->fetchrow_array();
    if (@tables and @tables == 1) {
        return 1;
    }
    else {
        return 0;
    }
}

sub create_table {
    my $sth = $dbh->prepare("create table $tablename (id integer primary key, assay_version text, pipeline_version text, date text, normal_sample text, tumor_sample text, chromosome text, position integer, rs_id text, reference text, variant text, variant_callers text, normal_gt text, normal_ad text, normal_af float, normal_dp integer, tumor_gt text, tumor_ad text, tumor_af float, tumor_dp integer, consequence text, symbol text, feature_type text, feature text, HGVSc text, HGVSp text, cDNA_position integer, cds_position integer, protein_position integer, amino_acids text, codons text, HGNC_ID integer, existing_variation text, gnomADe_af text)");
    $sth->execute or main->fatal_msg($DBI::errstr);
    $dbh->commit;
}
