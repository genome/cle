#! /usr/bin/perl

#Copyright (C) 2015 Feiyu Du <fdu@genome.wustl.edu>
#              and Washington University The Genome Institute

#This script is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY or the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#GNU General Public License for more details.


use strict;
use warnings;

umask 002;

use IO::File;
use File::Spec;

die "Provide sample_list, gene_list and data_dir" unless @ARGV == 3;
my ($sample_list, $gene_list, $dir) = @ARGV;

my $gene_fh = IO::File->new($gene_list);
my $gene_names = get_gene_names($gene_fh);

my $sample_fh = IO::File->new($sample_list);

while (my $name = $sample_fh->getline) {
    chomp $name;

    for my $docm ('', 1) {
        my $mid_path = "discovery/SeqCap\ EZ\ Human\ Exome\ v3.0\ +\ RMG-1\ pooled\ probes/merged/full/";
        my $out_name = $name;
        if ($docm) {
            $mid_path = "docm/SeqCap\ EZ\ Human\ Exome\ v3.0\ +\ RMG-1\ pooled\ probes/merged/docm";
            $out_name .= '_docm';
        }
        my $report   = File::Spec->join($dir, $name, $mid_path, 'report.txt');
        my $out_file = File::Spec->join($dir, $out_name.'_gene.out');
        parse_file($report, $gene_names, $out_file, $docm);
        print "$out_name done\n";
    }
}
$sample_fh->close;


sub parse_file {
    my ($in_file, $names, $out_file, $docm) = @_;
    my $in_fh  = IO::File->new($in_file);
    my $out_fh = IO::File->new(">$out_file");

    my $index = $docm ? 9 : 11;

    while (my $line = $in_fh->getline) {
        if ($line =~ /^chromosome_name/) {
            $out_fh->print($line);
            next;
        }
        my @columns = split /\t/, $line;
        if ($names->{$columns[$index]}) {
            if ($docm) {
                if ($columns[15] > 5 and $columns[16] > 1) {
                    $out_fh->print($line);
                }
            }
            else {
                $out_fh->print($line);
            }
        }
    }
    $in_fh->close;
    $out_fh->close;
    return;
}

sub get_gene_names {
    my ($fh) = @_;
    my $names;
    while (my $name = $fh->getline) {
        chomp $name;
	$name =~ s/\s+$//;
        $names->{$name} = 1;
    }
    $fh->close;
    return $names;
}

