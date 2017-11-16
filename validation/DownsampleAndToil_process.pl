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
use Genome;

my $dir = '/gscmnt/gc13016/cle/54f8f7b915cb472aa183c721307369ab_scratch_space/IDT_validation/downsample';
my $cram_dir = '/gscmnt/gc13016/cle/54f8f7b915cb472aa183c721307369ab_scratch_space/IDT_validation/wdl_cram';

my $fh = Genome::Sys->open_file_for_reading($ARGV[0]);
my $template_json = File::Spec->join($dir, 'template.json');

while (my $pair = $fh->getline) {
    chomp $pair;
    my ($t_sample, $n_sample) = split /\s+/, $pair;
    my $id;

    for my $ratio qw(0.85 0.7 0.55 0.4) {
        $id = join '_', $t_sample, 'x', $n_sample, $ratio;
        my $id_dir = $dir.'/'.$id;
        Genome::Sys->create_directory($id_dir);

        my $tumor_cram  = get_cram_for_sample($cram_dir, $t_sample);
        my $normal_cram = get_cram_for_sample($cram_dir, $n_sample);

        my $json = File::Spec->join($id_dir, $id.'.json');
        my $json_fh = Genome::Sys->open_file_for_writing($json);

        my $in_json_fh = Genome::Sys->open_file_for_reading($template_json);

        while (my $l = $in_json_fh->getline) {
            if ($l =~ /TumorCramPath/) {
                $l =~ s/TumorCramPath/$tumor_cram/;
            }
            elsif ($l =~ /NormalCramPath/) {
                $l =~ s/NormalCramPath/$normal_cram/;
            }
            elsif ($l =~ /OutputDirPath/) {
                $l =~ s/OutputDirPath/$id_dir/;
            }
            elsif ($l =~ /DownsampleRatioValue/) {
                $l =~ s/DownsampleRatioValue/$ratio/;
            }
            elsif ($l =~ /TumorLabelName/) {
                $l =~ s/TumorLabelName/$t_sample/;
            }
            elsif ($l =~ /NormalLabelName/) {
                $l =~ s/NormalLabelName/$n_sample/;
            }

            $json_fh->print($l);
        }
        $in_json_fh->close;
        $json_fh->close;
        
        my $out = File::Spec->join($dir, 'log', $id.'.out');
        my $err = File::Spec->join($dir, 'log', $id.'.err');

        `bsub -oo $out -eo $err -q research-hpc -a "docker(sleongmgi/cromwell:develop-with-mysql)" /usr/bin/java -Dconfig.file=application.conf -jar /cromwell/cromwell.jar run DownsampleAndToil.wdl $json`;
        print "$id done\n";
    }
}
    
exit;


sub get_cram_for_sample {
    my $dir = shift;
    my $sample_name = shift;
    
    my $sample_cram_dir = $dir . '/'. $sample_name;

    unless (-d $sample_cram_dir) { die('Missing sample cram dir:'. $sample_cram_dir); }

    my $sample_cram_path = $sample_cram_dir .'/'. $sample_name .'.cram';
    unless (-s $sample_cram_path) { die('Missing sample cram file: '. $sample_cram_path); }

    unless (-s $sample_cram_path .'.crai') { die('Missing sample cram index: '. $sample_cram_path .'.crai'); }
    my $sample_crai_path = $sample_cram_dir .'/'. $sample_name .'.crai';
    unless (-l $sample_crai_path or -s $sample_crai_path) {
        print "create symlink for $sample_crai_path\n";
        Genome::Sys->create_symlink($sample_cram_path.'.crai', $sample_crai_path);
    }
    return $sample_cram_path;
}

