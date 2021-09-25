workflow FastqToTargetedAlignedBam {

  String BaitBedFile = "/storage1/fs1/gtac-mgi/Active/CLE/reference/production_reference_GRCh38DH/accessory_bed_files/TCP_bait_hg38.bed"
  String TargetBedFile = "/storage1/fs1/gtac-mgi/Active/CLE/reference/production_reference_GRCh38DH/accessory_bed_files/TCP_target_hg38.bed"
  String ContBedFile = "/storage1/fs1/gtac-mgi/Active/CLE/reference/production_reference_GRCh38DH/accessory_bed_files/haplotyping_snps.010918.hg38.bed"
  String Reference = "/storage1/fs1/gtac-mgi/Active/CLE/reference/production_reference_GRCh38DH/reference/all_sequences.fa"
  String Dictionary = "/storage1/fs1/gtac-mgi/Active/CLE/reference/production_reference_GRCh38DH/reference/all_sequences.dict"

  Array[String] Read1s 
  Array[String] Read2s
  String JobGroup
  String TMPDIR
  String OutputDir
  String FinalLabel
  String ReadGroup 
  String Queue

  call bed_to_interval_list as BaitIntervals {
       input: bed=BaitBedFile,
       	      dict=Dictionary,
	      jobGroup=JobGroup,
              queue=Queue
  }

  call bed_to_interval_list as TargetIntervals {
       input: bed=TargetBedFile,
              dict=Dictionary,
	      jobGroup=JobGroup,
              queue=Queue
  }

  call cat_fastqs as cat_fastqs_read1 {
       input: readNum="1",
              fastqs=Read1s,
              finalLabel=FinalLabel,
	      jobGroup=JobGroup,
              queue=Queue
  }
 
  call cat_fastqs as cat_fastqs_read2 {
       input: readNum="2",
              fastqs=Read2s,
              finalLabel=FinalLabel,
	      jobGroup=JobGroup,
              queue=Queue
  }

  call align_and_tag {
       input: fastqs=[cat_fastqs_read1.Fastq,cat_fastqs_read2.Fastq],
              refFasta=Reference,
              jobGroup=JobGroup,
              readGroup=ReadGroup,
              queue=Queue
  } 

  call name_sort {
       input: tmp=TMPDIR,
              alignedBam=align_and_tag.TaggedBam,
              jobGroup=JobGroup,
              queue=Queue
  }

  call mark {
       input: tmp=TMPDIR,
              mergedBam=name_sort.SortedBam,
              label=FinalLabel,
              jobGroup=JobGroup,
              queue=Queue
  }

  call run_haplotect {
        input: Bam=mark.SortedBam,
            Bed=ContBedFile,
            Name=FinalLabel,
            MinReads=50,
            refFasta=Reference,
            queue=Queue,
            jobGroup=JobGroup
  }

  call collect_alignment_metrics {
       input: in=mark.SortedBam,
              ref=Reference,
              jobGroup=JobGroup,
              queue=Queue
  }

  call collect_gc_metrics {
       input: in=mark.SortedBam,
              ref=Reference,
              jobGroup=JobGroup,
              queue=Queue
  }

  call collect_insert_metrics {
       input: in=mark.SortedBam,
              ref=Reference,
              jobGroup=JobGroup,
              queue=Queue
  }

  call collect_hs_metrics {
       input: in=mark.SortedBam,
              ref=Reference,
	      baits=BaitIntervals.intervals,
	      targets=TargetIntervals.intervals,
              jobGroup=JobGroup,
              queue=Queue
  }

  call flagstat {
       input: in=mark.SortedBam,
              jobGroup=JobGroup,
              queue=Queue
  }

  call bamutil {
       input: in=mark.SortedBam,
              jobGroup=JobGroup,
              queue=Queue
  }

  call remove_file {
       input: file=align_and_tag.TaggedBam,
	      order_by=name_sort.SortedBam,
              jobGroup=JobGroup,
              queue=Queue
  }

  call gather_result as gather_the_rest {
       input: files=[cat_fastqs_read1.Fastq,
                     cat_fastqs_read2.Fastq,
                     mark.SortedBam,
	             mark.SortedBamIndex,
                     collect_alignment_metrics.alignMetrics,
                     collect_gc_metrics.gcOut,
                     collect_gc_metrics.gcSum,
                     collect_gc_metrics.gcPDF,
                     collect_insert_metrics.isOut,
                     collect_insert_metrics.isPDF,
		     collect_hs_metrics.hsMetrics,
		     collect_hs_metrics.perTargetCoverage,
                     flagstat.fsOut,
                     bamutil.bamutilOut,
                     mark.MarkMetrics,
                     run_haplotect.out_file,
                     run_haplotect.sites_file,
                     ],
              dir=OutputDir,
              jobGroup=JobGroup,
              queue=Queue

  }
}

task bed_to_interval_list {
     String bed
     String dict
     String jobGroup
     String queue

     command {
     	      /usr/bin/java -Xmx4g -jar /usr/picard/picard.jar BedToIntervalList I=${bed} O="interval_list" SD=${dict}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/sort-mark-duplicates:2"
             cpu: "1"
             memory: "4 G"
             queue: queue
             job_group: jobGroup
     }
     output {
     	    File intervals = "interval_list"
     }
}

task cat_fastqs {
     String readNum
     Array[String] fastqs
     String finalLabel
     String jobGroup
     String queue

     command {
            /bin/cat ${sep=" " fastqs} > "${finalLabel}.${readNum}.fastq.gz"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
             queue: queue
             job_group: jobGroup
     }
     output {
            File Fastq = "${finalLabel}.${readNum}.fastq.gz"
     }
}

task align_and_tag {
     Array[String] fastqs
     String refFasta
     String readGroup
     String jobGroup
     String queue

     command {
            (set -eo pipefail && /usr/local/bin/bwa mem -K 100000000 -t 8 -Y -R "${readGroup}" ${refFasta} ${sep=" " fastqs} | /usr/local/bin/samblaster -a --addMateTags | /usr/local/bin/samtools view -b -S /dev/stdin > "refAlign.bam")
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/tagged-alignment:2"
             cpu: "8"
             memory: "20 G"
             queue: queue
             job_group: jobGroup
     }
     output {
            File TaggedBam = "refAlign.bam"
     }
}

task name_sort {
     String alignedBam
     String tmp
     String jobGroup
     String queue

     command {
             /usr/local/bin/sambamba sort -t 8 -m 18G -n --tmpdir=${tmp} -o "NameSorted.bam" ${alignedBam}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/sambamba-0.6.4:1"
             cpu: "8"
             memory: "20 G"
             queue: queue
             job_group: jobGroup
     }
     output {
             File SortedBam = "NameSorted.bam"
     }
}

task mark {
     String mergedBam
     String jobGroup
     String tmp
     String label
     String queue

     command {
             (set -eo pipefail && /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar MarkDuplicates I=${mergedBam} O=/dev/stdout ASSUME_SORT_ORDER=queryname METRICS_FILE=mark_dups_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | /usr/local/bin/sambamba sort -t 8 -m 18G --tmpdir=${tmp} -o "${label}.bam" /dev/stdin)
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/sort-mark-duplicates:2"
             cpu: "8"
             memory: "50 G"
             queue: queue
             job_group: jobGroup
     }
     output {
             File SortedBam = "${label}.bam"
             File SortedBamIndex = "${label}.bam.bai"
             File MarkMetrics = "mark_dups_metrics.txt"
     }
     
}

task run_haplotect {
     String Bam
     String Bed
     String Name   
     String refFasta
     String queue
     String jobGroup

     Int? MinReads

     command <<<
         awk -v OFS="\t" '{print $1,$2-1,$3;}' ${Bed} > /tmp/haplotect.bed && \
         /usr/bin/java -Xmx6g \
         -jar /opt/gatk/public/external-example/target/external-example-1.0-SNAPSHOT.jar \
         -T Haplotect -R ${refFasta} \
         -mmq 20 -mbq 20 -dcov 20000 -unif 0 -gstol 0.001 \
         -mr ${default=100 MinReads} \
         -I ${Bam} \
         -htp ${Bed} \
         -L /tmp/haplotect.bed \
         -outPrefix ${Name} && \
         sort -u "${Name}.multihaploci.txt" > "haplotectloci.txt" && \
         /usr/bin/perl /opt/CalculateContamination.pl "${Name}.txt" "haplotectloci.txt" > "haplotect.txt"
     >>>

     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/docker-cle-contamination:1"
             cpu: "1"
             memory: "8 G"
             queue: queue
             job_group: jobGroup
     }
     output {
            File out_file = "haplotect.txt"
            File sites_file = "haplotectloci.txt"
     }
}

task collect_alignment_metrics {
     String in
     String ref
     String jobGroup
     String queue

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=${ref} INPUT=${in} OUTPUT="alignment_summary.txt" ASSUME_SORTED=true
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory: "18 G"
             queue: queue
             job_group: jobGroup
     }
     output {
            File alignMetrics = "alignment_summary.txt"
     }
}

task collect_gc_metrics {
     String in
     String ref
     String jobGroup
     String queue

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectGcBiasMetrics REFERENCE_SEQUENCE=${ref} INPUT=${in} OUTPUT="GC_bias.txt" SUMMARY_OUTPUT="GC_bias_summary.txt" CHART_OUTPUT="GC_bias_chart.pdf" ASSUME_SORTED=true
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory: "18 G"
             queue: queue
             job_group: jobGroup
     }
     output {
            File gcOut = "GC_bias.txt"
            File gcSum = "GC_bias_summary.txt"
            File gcPDF = "GC_bias_chart.pdf"
     }     
}

task collect_insert_metrics {
     String in
     String ref
     String jobGroup
     String queue

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectInsertSizeMetrics INPUT=${in} OUTPUT="insert_size_summary.txt" HISTOGRAM_FILE="insert_size.pdf" ASSUME_SORTED=true
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory: "18 G"
             queue: queue
             job_group: jobGroup
     }
     output {
            File isOut = "insert_size_summary.txt"
            File isPDF = "insert_size.pdf"
     }
     
}
task collect_hs_metrics {
     String in
     String ref
     String baits
     String targets
     String jobGroup
     String queue

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectHsMetrics REFERENCE_SEQUENCE=${ref} INPUT=${in} \
	     		   BAIT_INTERVALS=${baits} TARGET_INTERVALS=${targets} \
			   MINIMUM_MAPPING_QUALITY=1 MINIMUM_BASE_QUALITY=1 \
			   OUTPUT="hs_metric_summary.txt" PER_TARGET_COVERAGE="per_target_coverage.txt"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory: "18 G"
             queue: queue
             job_group: jobGroup
     }
     output {
            File hsMetrics = "hs_metric_summary.txt"
	    File perTargetCoverage = "per_target_coverage.txt"
     }
     
}

task flagstat {
     String in
     String jobGroup
     String queue

     command {
             /usr/local/bin/samtools flagstat ${in} > "flagstat.out"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/tagged-alignment:2"
             cpu: "1"
             memory: "10 G"
             queue: queue
             job_group: jobGroup
     }
     output {
            File fsOut = "flagstat.out"
    }
}

task bamutil {
     String in
     String jobGroup
     String queue

     command {
             /usr/local/bin/bam stats --noPhoneHome --in ${in} --phred --excludeFlags 3844 2> bamutil_stats.txt
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/bamutil:2"
             cpu: "1"
             memory: "10 G"
             queue: queue
             job_group: jobGroup
     }
     output {
            File bamutilOut = "bamutil_stats.txt"
    }
}

task gather_result {
     String dir
     Array[String] files
     String jobGroup
     String queue

     command {
             /bin/mv -t ${dir} ${sep=" " files}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
             queue: queue
             job_group: jobGroup
     }
     output {
            String out = stdout()
     }
}

task remove_file {
     String file
     String? order_by
     String jobGroup
     String queue

     command {
             /bin/rm ${file}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
             queue: queue
             job_group: jobGroup
     }
     output {
            String out = stdout()
     }
}
