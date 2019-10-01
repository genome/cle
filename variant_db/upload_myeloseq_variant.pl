#! /usr/bin/perl

use strict;
use warnings;

use DBI;
use POSIX;
use Genome;
use File::Basename;

my @batch_dirs = @ARGV;

my $dir = '/gscmnt/gc13016/cle/54f8f7b915cb472aa183c721307369ab_scratch_space/myeloseq/RUN/wdl_out';
my $dbfile = '/gscmnt/gc3042/cle_validation/CLE_variant_database/sqlite_variant_DB/cle_variants.sqlite';

my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', { AutoCommit => 0, RaiseError => 1 })
    or main->fatal_msg("Can't connect to db ($dbfile): " . $DBI::errstr);

my $tablename = 'myeloseq_variant';
#create the database
unless (check_table()) {
    create_table();
}

unless (@batch_dirs) {
    print "Go through the whole myeloseq RUN Batch dirs\n";
    @batch_dirs = glob($dir."/Batch*");
}

my @excludes = map{'Batch'.$_}qw(29 2_version1_1 30 30_2ndRun 39 47_Manual 5.2 9_new_test2);

BATCH: for my $batch_dir (@batch_dirs) {
    my $base = basename $batch_dir;
    if (grep{$base eq $_}@excludes) {
        print "$base is not valid\n";
        next;
    }
    
    opendir(my $dir_h, $batch_dir);
    CASE: for my $case (readdir $dir_h) {
        next if $case =~ /^\./;
        next if $case =~ /^cromwell\-/;
        my $case_dir = File::Spec->join($batch_dir, $case);
        next unless -d $case_dir;
        next if -l $case_dir;

        my ($case_name) = $case =~ /^(\S+lib\d+)_[ATCG]{8}/;
        
        my $tsv = File::Spec->join($case_dir, $case_name.'.variants_annotated.tsv');
        my $report = File::Spec->join($case_dir, $case_name.'.variant_report.txt');
   
        my $fh = Genome::Sys->open_file_for_reading($tsv);
        my $date = POSIX::strftime("%Y-%m-%d %H:%M:%S", localtime((stat $tsv)[9]));
    
        my $key='"'.'Myeloseq Assay Version'.'"';
        my $str = `grep $key $report`;
        my ($version) = $str =~ /Version\s+(\S+)/;
        unless ($version) {
            die "$base : $case does not have valid version";
        }

        TSV: while (my $line = $fh->getline) {
            chomp $line;
            next if $line =~ /^CHROM/;
            my @columns = split /\t/, $line;
            next unless $columns[3] eq "PASS";
            
            my $col_ct = scalar @columns;
            die "Expect 25 columns from myeloseq upload file" unless $col_ct =~ /^(25)$/;

            my @info = ($version, $version, $date, $case_name);
            push @info, map{$columns[$_]}(0..1, 4..8, 10..24);

            my $sth = $dbh->prepare("insert into $tablename (assay_version, pipeline_version, date, sample, chromosome, position, reference, variant, variant_callers, TAMP, SAMP, CVAF, NR, NV, consequence, gene_symbol, exon, intron, feature_type, feature, HGVSc, HGVSp, HGNC_ID, MAX_AF, MYELOSEQ_TCGA_AC, MYELOSEQ_MDS_AC) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)") or die $DBI::errstr;
    
            $sth->execute_array({}, @info) or main->fatal_msg($DBI::errstr);
            $dbh->commit;
        }
        $fh->close;
    }
    closedir $dir_h;
    print "$base is done\n";
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
    #For now not store ID, FILTER, VAFTYPE columns
    my $sth = $dbh->prepare("create table $tablename (id integer primary key, assay_version text, pipeline_version text, date text, sample text, chromosome text, position integer, reference text, variant text, variant_callers text, TAMP integer, SAMP integer, CVAF float, NR integer, NV integer, consequence text, gene_symbol text, exon text, intron text, feature_type text, feature text, HGVSc text, HGVSp text, HGNC_ID integer, MAX_AF float, MYELOSEQ_TCGA_AC integer, MYELOSEQ_MDS_AC integer)");
    $sth->execute or main->fatal_msg($DBI::errstr);
    $dbh->commit;

    #$sth = $dbh->prepare("create unique index variant on $table_name(id)");
    #$sth->execute or main->fatal_msg($DBI::errstr);
    #$dbh->commit;
}
