workflow ProcessHaloplexHS {

  File SampleSheet
  # sample sheet has this structure:
  # index  sample    other stuff

  Array[Array[String]] inputData = read_tsv(SampleSheet)

  String IlluminaDir

  String TargetBed
  String AmpliconBed  
  String CoverageBed

  String JobGroup
  String OutputDir

  Array[String] Adapters = ["NNNGATCGGAAGAGCACACGTCTGAACTCCAGTCAC","NNNAGATCGGAAGAGCGTCGTGTAGGGAAA"]

  String Reference = "/gscuser/dspencer/refdata/GRCh37/all_sequences.fa"
  String ReferenceIndex = "/gscuser/dspencer/refdata/GRCh37/all_sequences.fa.fai"
  String Dictionary = "/gscuser/dspencer/refdata/GRCh37/all_sequences.dict"

  call barcode_demux {
       input: Dir=IlluminaDir, #Fastqs=get_fastq_files.fastq_files,
              SampleSheet=SampleSheet,
	      SampleIndexMM=0,
              jobGroup=JobGroup
  }

  call prepare_samples {
       input: SampleSheet=SampleSheet,
              Fastq1=barcode_demux.read1_fastqs,
              Fastq2=barcode_demux.read2_fastqs,
              jobGroup=JobGroup
  }

  scatter (samples in inputData){
  	 call trim_reads {
               input: Index=samples[0],
	       	      SampleSheetFile=prepare_samples.sample_sheet,
                      Adapters=Adapters,
                      jobGroup=JobGroup
	 }
   	 call align_barcode_and_sort_reads {
	       input: Fastq=trim_reads.fastq_file,
		      refFasta=Reference,
		      jobGroup=JobGroup
         }
         call consensus_bam {
               input: Bam=align_barcode_and_sort_reads.bam_file,
		      BamIndex=align_barcode_and_sort_reads.bam_index,
                      TargetBed=TargetBed,
                      refFasta=Reference,
                      jobGroup=JobGroup
         }

	 call haloplex_qc {
	      input: refFasta=Reference,
	      	     AlignedBam=align_barcode_and_sort_reads.bam_file,
	      	     ConsensusBam=consensus_bam.bam_file,
		     TargetBed=TargetBed,
		     AmpliconBed=AmpliconBed,
		     CoverageBed=CoverageBed,
		     Name=samples[1],
		     DemuxFile=prepare_samples.sample_sheet,
		     jobGroup=JobGroup
	 }

	 call run_varscan {
              input: Bam=consensus_bam.bam_file,
                     BamIndex=consensus_bam.bam_index,
                     CoverageBed=CoverageBed,
                     refFasta=Reference,
		     Name=samples[1],
                     jobGroup=JobGroup
         }


	 call run_platypus {
             input: Bam=consensus_bam.bam_file,
                    BamIndex=consensus_bam.bam_index,
                    CoverageBed=CoverageBed,
		    Name=samples[1],
                    refFasta=Reference,
                    jobGroup=JobGroup
        }

	 call run_pindel_region as run_pindel_flt3itd {
              input: Bam=consensus_bam.bam_file,
	             BamIndex=consensus_bam.bam_index,
	             Reg='13:28608124-28608453',
	             refFasta=Reference,
	             Name=samples[1],
	             jobGroup=JobGroup
	}

	call combine_variants {
             input: VarscanSNV=run_varscan.varscan_snv_file,
                    VarscanIndel=run_varscan.varscan_indel_file,
                    PindelITD=run_pindel_flt3itd.pindel_vcf_file,
                    Platypus=run_platypus.platypus_vcf_file,
                    Bam=consensus_bam.bam_file,
                    BamIndex=consensus_bam.bam_index,
                    refFasta=Reference,
                    Name=samples[1],
                    jobGroup=JobGroup
	}
	
	call gather_files {
	      input: OutputFiles=[trim_reads.trimstats_file,
			   align_barcode_and_sort_reads.bam_file,
	      	           align_barcode_and_sort_reads.bam_index,
			   consensus_bam.bam_file,
			   consensus_bam.bam_index,
			   haloplex_qc.coverage_qc_file,
			   haloplex_qc.coverage_file,
			   haloplex_qc.amplicon_counts,
			   haloplex_qc.gc_length_plot,
			   haloplex_qc.coverage_plot,
			   run_varscan.varscan_snv_file,
			   run_varscan.varscan_indel_file,
			   run_pindel_flt3itd.pindel_vcf_file,
			   run_platypus.platypus_vcf_file,
			   combine_variants.combined_vcf_file],
		     OutputDir=OutputDir,
		     SubDir=samples[1] + "_" + samples[0],
		     jobGroup=JobGroup
	 }
  }
}

task barcode_demux {
     String Dir
     String SampleSheet 
     Int SampleIndexMM
     String jobGroup

     command <<<
     	     I1=$(ls ${Dir}/Fastq/*I1*.fastq.gz) && \
     	     I2=$(ls ${Dir}/Fastq/*I2*.fastq.gz) && \
     	     R1=$(ls ${Dir}/Fastq/*R1*.fastq.gz) && \
     	     R2=$(ls ${Dir}/Fastq/*R2*.fastq.gz) && \
	     cut -f 1 ${SampleSheet} > index_list.txt && \
     	     python /usr/bin/demultiplexer_hja_bd.py \
	     --mismatches ${SampleIndexMM} --index-read-file $I1 --tag-file index_list.txt \
	     --read1-file $R1 --read2-file $R2 --read3-file $I2
     >>>
     runtime {
             docker_image: "dhspence/docker-haloplex_demux:latest"
             cpu: "1"
             memory_gb: "12"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=12000]"
             job_group: jobGroup 
     }
     output {
	    Array[File] read1_fastqs = glob("*.1.fq.gz") 
	    Array[File] read2_fastqs = glob("*.2.fq.gz")
	    File unknown_read1 = "unknown.1.fq.gz" 
	    File unknown_read2 = "unknown.2.fq.gz" 
     }
}

task prepare_samples {
     File SampleSheet
     Array[File] Fastq1
     Array[File] Fastq2
     String jobGroup

     command <<<
             /bin/cat ${write_tsv(Fastq1)} > 1.tmp.txt
             /bin/cat ${write_tsv(Fastq2)} > 2.tmp.txt
             /usr/bin/perl -e 'open(R1,"1.tmp.txt"); @r1 = <R1>; \
                  chomp @r1; close R1;\
                  open(R2,"2.tmp.txt"); @r2 = <R2>; \
                  chomp @r2; close R2; \
                  open(SS,"${SampleSheet}");
                  while(<SS>){
                      chomp;
                      my @l = split("\t",$_);
                      my $r1 = (grep /$l[0].1.fq.gz/, @r1)[0];
                      my $r2 = (grep /$l[0].2.fq.gz/, @r2)[0];
		      my $persamplereads1 = `gunzip -c $r1 | wc -l`;
		      chomp $persamplereads1;
		      my $persamplereads2 = `gunzip -c $r2 | wc -l`;
		      chomp $persamplereads2;
                      print join("\t",$l[0],$l[1],$r1,$r2,($persamplereads1 / 4) + ($persamplereads2 / 4)),"\n";
                 }
                 close SS;
		 my $r1 = (grep /unknown.1.fq.gz/, @r1)[0];
                 my $r2 = (grep /unknown.2.fq.gz/, @r2)[0];
		 my $persamplereads1 = `gunzip -c $r1 | wc -l`;
		 chomp $persamplereads1;		      
		 my $persamplereads2 = `gunzip -c $r2 | wc -l`;
   		 chomp $persamplereads2;
		 print join("\t","unknown","unknown",$r1,$r2,($persamplereads1 / 4) + ($persamplereads2 / 4)),"\n"' > sample_sheet.txt
     >>>
     runtime {
             docker_image: "ubuntu:xenial"
             cpu: "1"
             memory_gb: "4"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=4000]"
             job_group: jobGroup
     }
     output {
            File sample_sheet = "sample_sheet.txt"
	    Array[Array[String]] sample_data = read_tsv("sample_sheet.txt")
     }
}

task trim_reads {
     String Index
     File SampleSheetFile
     Array[String] Adapters
     String jobGroup

     command {
     	export PYTHONPATH=/opt/cutadapt/lib/python2.7/site-packages/ && \
	/opt/cutadapt/bin/cutadapt --interleaved -u 4 -m 40 -a ${Adapters[0]} -A ${Adapters[1]} -o trimmed.fq \
	$(/bin/grep ${Index} ${SampleSheetFile} | cut -f 3) $(/bin/grep ${Index} ${SampleSheetFile} | cut -f 4) > trimstats.txt && \
	gzip -c trimmed.fq > trimmed.fq.gz
     }
     runtime {
             docker_image: "dhspence/docker-cutadapt:latest"
             cpu: "1"
             memory_gb: "8"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=8000]"
             job_group: jobGroup 
     }
     output {
	    File fastq_file = "trimmed.fq.gz"
	    File trimstats_file = "trimstats.txt"
     }
}

task align_barcode_and_sort_reads {
     String Fastq
     String refFasta
     String jobGroup
     
     command <<<     	     
	     (set -eo pipefail && /usr/local/bin/bwa mem -M -t 8 -p ${refFasta} ${Fastq} \
	     -R "@RG\tPU:XX\tID:1\tLB:YY\tSM:Sample\tPL:ILLUMINA\tCN:WUGSC" | \
	     /usr/bin/perl -ane '{ if (/^@/) { print join("\t", @F),"\n"; 
	      } elsif ($F[0] =~ /:([ACGT]{8,10})$/){ 
	       push @F, "X0:Z:$1"; 
	       $F[0] =~ s/:[ACGT]{8,10}$//; 
	       print join("\t", @F), "\n"; 
	      }
	     }' | /usr/local/bin/samtools view -b -S /dev/stdin | \
	     /usr/local/bin/samtools sort -m 1G -O bam - > "aligned_sorted.bam") && \
	     /usr/local/bin/samtools index aligned_sorted.bam
     >>>
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/tagged-alignment:2"
             cpu: "8"
             memory_gb: "32"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=32000]"
             job_group: jobGroup
     }
     output {
     	    File bam_file = "aligned_sorted.bam"
     	    File bam_index = "aligned_sorted.bam.bai"
     }
}

task haloplex_qc {
     String refFasta
     String DemuxFile
     String AlignedBam
     String ConsensusBam 
     String TargetBed
     String AmpliconBed
     String CoverageBed
     String Name
     Int? Mincov1
     Int? Mincov2
     String jobGroup

     command {
     	     /usr/bin/perl /usr/local/bin/CalculateCoverageQC.072617.pl ${refFasta} ${AmpliconBed} ${TargetBed} ${CoverageBed} \
	     ${default=50 Mincov1} ${default=20 Mincov2} ${DemuxFile} ${AlignedBam} ${ConsensusBam} ${Name}
     }
     runtime {
     	     docker_image: "dhspence/docker-haloplexqc:latest"
             cpu: "1"
             memory_gb: "16"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=16000]"
             job_group: jobGroup
     } 
     output {
     	    File coverage_qc_file = "${Name}.qc.txt"
     	    File coverage_file = "${Name}.coverage.txt"
     	    File amplicon_counts = "${Name}.amplicon_counts.txt"
     	    File coverage_plot = "${Name}.coverage_qc.pdf"
     	    File gc_length_plot = "${Name}.gc_length_qc.pdf"
     }
}

task run_varscan {
     File Bam
     File BamIndex
     Int? MinCov
     Float? MinFreq
     Int? MinReads
     String CoverageBed
     String refFasta
     String Name
     String jobGroup

     command <<<
             /usr/bin/samtools mpileup -f ${refFasta} -l ${CoverageBed} ${Bam} > /tmp/mpileup.out && \
	     java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2snp /tmp/mpileup.out --min-coverage ${default=8 MinCov} --min-reads2 ${default=5 MinReads} \
	     	  --min-var-freq ${default="0.01" MinFreq} --output-vcf > ${Name}.snv.vcf && \
	     java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2indel /tmp/mpileup.out --min-coverage ${default=8 MinCov} --min-reads2 ${default=5 MinReads} \
	     	  --min-var-freq ${default="0.01" MinFreq} --output-vcf > ${Name}.indel.vcf
     >>>

     runtime {
             docker_image: "mgibio/varscan-cwl:v2.4.2-samtools1.3.1"
             cpu: "1"
             memory_gb: "16"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=16000]"
             job_group: jobGroup
     }
     output {
            File varscan_snv_file = "${Name}.snv.vcf"
            File varscan_indel_file = "${Name}.indel.vcf"
     }
}

task run_pindel_region {
    File Bam
    File BamIndex
    String Reg
    Int? Isize
    Int? MinReads
    String refFasta
    String Name
    String jobGroup

    command <<<
        (set -eo pipefail && /usr/local/bin/samtools view ${Bam} ${Reg} | /opt/pindel-0.2.5b8/sam2pindel - /tmp/in.pindel ${default=250 Isize} tumor 0 Illumina-PairEnd) && \
	/usr/bin/pindel -f ${refFasta} -p /tmp/in.pindel -c ${Reg} -o /tmp/out.pindel && \
	/usr/bin/pindel2vcf -P /tmp/out.pindel -G -r ${refFasta} -e ${default=3 MinReads} -R GRCh37 -d GRCh37 -v ${Name}.pindel.vcf
    >>>

    runtime {
        docker_image: "mgibio/pindel2vcf-cwl:0.6.3"
        cpu: "1"
        memory_gb: "16"
        queue: "research-hpc"
        resource: "rusage[gtmp=10, mem=16000]"
        job_group: jobGroup
    }
    output {
        File pindel_vcf_file = "${Name}.pindel.vcf"
    }
}

task consensus_bam {
     File Bam
     File BamIndex
     String TargetBed
     String refFasta
     String jobGroup

     command {
             /usr/bin/java -Xmx6g \
	     -jar /opt/gatk/public/external-example/target/external-example-1.0-SNAPSHOT.jar \
	     -T WalkerTRConsensus -I ${Bam} -mmq 20 \
	     -L ${TargetBed} -dcov 500000 \
	     -cbam consensus.bam \
	     -R ${refFasta} && \
	     /usr/bin/samtools index consensus.bam
     }

     runtime {
             docker_image: "dhspence/docker-walker"
             cpu: "1"
             memory_gb: "8"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=8000]"
             job_group: jobGroup
     }
     output {
            File bam_file = "consensus.bam"
            File bam_index = "consensus.bam.bai"
     }
}

task run_platypus {
     File Bam
     File BamIndex
     String CoverageBed
     String? DocmVcf
     String Name
     String refFasta
     String jobGroup

     command <<<
     	     /usr/bin/awk '{ print $1":"$2+1"-"$3; }' ${CoverageBed} > "regions.txt" && \
     	     python /opt/platypus/Platypus_0.8.1/Platypus.py callVariants --bamFiles=${Bam} \
	     	 --regions regions.txt --refFile=${refFasta} \
	         --nCPU 1 ${"--source=" + DocmVcf} --output="${Name}.platypus.vcf" \
	         --minReads 5 --filterDuplicates=0 --minFlank=10 --assemble 1 --filterReadPairsWithSmallInserts=0 \
	         --trimOverlapping=0 --trimSoftClipped=0 --filterReadsWithDistantMates=0 --filterReadsWithUnmappedMates=0
     >>>
     runtime {
             docker_image: "dhspence/docker-platypus:cle"
             cpu: "4"
             memory_gb: "16"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=16000]"
             job_group: jobGroup
     }
     output {
     	    File platypus_vcf_file = "${Name}.platypus.vcf"
     }
}

task combine_variants {
     String VarscanSNV
     String VarscanIndel
     String PindelITD
     String Platypus
     String Bam
     String BamIndex
     String refFasta
     String Name
     String jobGroup

     command {
     	     /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T CombineVariants -R ${refFasta} --variant:varscanIndel ${VarscanIndel} \
	     --variant:varscanSNV ${VarscanSNV} --variant:Platypus ${Platypus} --variant:PindelITD ${PindelITD} -o /tmp/out.vcf --genotypemergeoption UNIQUIFY && \
	     /usr/bin/java -Xmx16g -jar /opt/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R ${refFasta} --variant /tmp/out.vcf -o combined.vcf && \
	     python /usr/local/bin/addAmpliconInfoAndCountReads.py -r ${refFasta} combined.vcf ${Bam} ${Name} > ${Name}.combined_and_tagged.vcf
     }
     runtime {
     	     docker_image: "dhspence/docker-amplicon-readcount:latest"
             cpu: "1"
             memory_gb: "10"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=10000]"
             job_group: jobGroup
     }
     output {
            File combined_vcf_file = "${Name}.combined_and_tagged.vcf"
     }

}

task gather_files {
     Array[String] OutputFiles
     String OutputDir
     String? SubDir
     String jobGroup

     command {
	     if [[ ${SubDir} != "" ]] && [[ ! -e ${OutputDir}/${SubDir} ]]; then
	     	mkdir ${OutputDir}/${SubDir}
	     fi
             /bin/mv -f -t ${OutputDir}/${SubDir} ${sep=" " OutputFiles}
     }
     runtime {
             docker_image: "ubuntu:xenial"
             queue: "research-hpc"
             job_group: jobGroup
     }
     output {
            String out = stdout()
     }
}

task return_object {
     Array[Object] obj
     command {
     	     cat ${write_objects(obj)} > "obj.tsv"
     }

     output {
     	    File results = "obj.tsv"
     }
}