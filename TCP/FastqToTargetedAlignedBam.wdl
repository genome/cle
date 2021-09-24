workflow FastqToTargetedAlignedBam {

  String BaitBedFile = "/gscmnt/gc2709/info/production_reference_GRCh38DH/accessory_bed_files/TCP_bait_hg38.bed"
  String TargetBedFile = "/gscmnt/gc2709/info/production_reference_GRCh38DH/accessory_bed_files/TCP_target_hg38.bed"
  String ContBedFile = "/gscmnt/gc2709/info/production_reference_GRCh38DH/accessory_bed_files/haplotyping_snps.010918.hg38.bed"
  String Reference = "/gscmnt/gc2709/info/production_reference_GRCh38DH/reference/all_sequences.fa"
  String Dictionary = "/gscmnt/gc2709/info/production_reference_GRCh38DH/reference/all_sequences.dict"

  Array[String] Read1s 
  Array[String] Read2s
  String JobGroup
  String TMPDIR
  String OutputDir
  String FinalLabel
  String ReadGroup 
  Int Priority
  String Project
  String Queue

  Int bed_to_interval_list_priority_1 = Priority + 1
  String bed_to_interval_list_priority_value_1 = " " + bed_to_interval_list_priority_1
  call bed_to_interval_list as BaitIntervals {
       input: bed=BaitBedFile,
       	      dict=Dictionary,
	      jobGroup=JobGroup,
              priority=bed_to_interval_list_priority_value_1,
              queue=Queue,
              project=Project
  }

  Int bed_to_interval_list_priority_2 = Priority + 2
  String bed_to_interval_list_priority_value_2 = " " + bed_to_interval_list_priority_2
  call bed_to_interval_list as TargetIntervals {
       input: bed=TargetBedFile,
              dict=Dictionary,
	      jobGroup=JobGroup,
              priority=bed_to_interval_list_priority_value_2,
              queue=Queue,
              project=Project
  }

  Int cat_fastqs_priority_3 = Priority + 3
  String cat_fastqs_priority_value_3 = " " + cat_fastqs_priority_3
  call cat_fastqs as cat_fastqs_read1 {
       input: readNum="1",
              fastqs=Read1s,
              finalLabel=FinalLabel,
	      jobGroup=JobGroup,
              priority=cat_fastqs_priority_value_3,
              queue=Queue,
              project=Project
  }
 
  Int cat_fastqs_priority_4 = Priority + 4
  String cat_fastqs_priority_value_4 = " " + cat_fastqs_priority_4
  call cat_fastqs as cat_fastqs_read2 {
       input: readNum="2",
              fastqs=Read2s,
              finalLabel=FinalLabel,
	      jobGroup=JobGroup,
              priority=cat_fastqs_priority_value_4,
              queue=Queue,
              project=Project
  }

  Int align_and_tag_priority_5 = Priority + 5
  String align_and_tag_priority_value_5 = " " + align_and_tag_priority_5
  call align_and_tag {
       input: fastqs=[cat_fastqs_read1.Fastq,cat_fastqs_read2.Fastq],
              refFasta=Reference,
              jobGroup=JobGroup,
              readGroup=ReadGroup,
              priority=align_and_tag_priority_value_5,
              queue=Queue,
              project=Project
  } 

  Int name_sort_priority_6 = Priority + 6
  String name_sort_priority_value_6 = " " + name_sort_priority_6
  call name_sort {
       input: tmp=TMPDIR,
              alignedBam=align_and_tag.TaggedBam,
              jobGroup=JobGroup,
              priority=name_sort_priority_value_6,
              queue=Queue,
              project=Project
  }

  Int mark_priority_7 = Priority + 7
  String mark_priority_value_7 = " " + mark_priority_7
  call mark {
       input: tmp=TMPDIR,
              mergedBam=name_sort.SortedBam,
              label=FinalLabel,
              jobGroup=JobGroup,
              priority=mark_priority_value_7,
              queue=Queue,
              project=Project
  }

  Int contamination_priority = Priority + 8 
  String contamination_priority_value = " " + contamination_priority
  call run_haplotect {
        input: Bam=mark.SortedBam,
            Bed=ContBedFile,
            Name=FinalLabel,
            MinReads=50,
            refFasta=Reference,
            queue=Queue,
            jobGroup=JobGroup
  }

  Int collect_alignment_metrics_priority_8 = Priority + 8
  String collect_alignment_metrics_priority_value_8 = " " + collect_alignment_metrics_priority_8
  call collect_alignment_metrics {
       input: in=mark.SortedBam,
              ref=Reference,
              jobGroup=JobGroup,
              priority=collect_alignment_metrics_priority_value_8,
              queue=Queue,
              project=Project
  }

  Int collect_gc_metrics_priority_9 = Priority + 9
  String collect_gc_metrics_priority_value_9 = " " + collect_gc_metrics_priority_9
  call collect_gc_metrics {
       input: in=mark.SortedBam,
              ref=Reference,
              jobGroup=JobGroup + 8,
              priority=collect_gc_metrics_priority_value_9,
              queue=Queue,
              project=Project
  }

  Int collect_insert_metrics_priority_10 = Priority + 10
  String collect_insert_metrics_priority_value_10 = " " + collect_insert_metrics_priority_10
  call collect_insert_metrics {
       input: in=mark.SortedBam,
              ref=Reference,
              jobGroup=JobGroup,
              priority=collect_insert_metrics_priority_value_10,
              queue=Queue,
              project=Project
  }

  Int collect_hs_metrics_priority_11 = Priority + 11
  String collect_hs_metrics_priority_value_11 = " " + collect_hs_metrics_priority_11
  call collect_hs_metrics {
       input: in=mark.SortedBam,
              ref=Reference,
	      baits=BaitIntervals.intervals,
	      targets=TargetIntervals.intervals,
              jobGroup=JobGroup,
              priority=collect_hs_metrics_priority_value_11,
              queue=Queue,
              project=Project
  }

  Int flagstat_priority_12 = Priority + 12
  String flagstat_priority_value_12 = " " + flagstat_priority_12
  call flagstat {
       input: in=mark.SortedBam,
              jobGroup=JobGroup,
              priority=flagstat_priority_value_12,
              queue=Queue,
              project=Project
  }

  Int bamutil_priority_13 = Priority + 13
  String bamutil_priority_value_13 = " " + bamutil_priority_13
  call bamutil {
       input: in=mark.SortedBam,
              jobGroup=JobGroup,
              priority=bamutil_priority_value_13,
              queue=Queue,
              project=Project
  }

  Int remove_file_priority_14 = Priority + 14
  String remove_file_priority_value_14 = " " + remove_file_priority_14
  call remove_file {
       input: file=align_and_tag.TaggedBam,
	      order_by=name_sort.SortedBam,
              jobGroup=JobGroup,
              priority=remove_file_priority_value_14,
              queue=Queue,
              project=Project
  }

  Int gather_result_priority_15 = Priority + 15
  String gather_result_priority_value_15 = " " + gather_result_priority_15
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
              priority=gather_result_priority_value_15,
              queue=Queue,
              project=Project

  }
}

task bed_to_interval_list {
     String bed
     String dict
     String jobGroup
     String priority
     String project
     String queue

     command {
     	      /usr/bin/java -Xmx4g -jar /usr/picard/picard.jar BedToIntervalList I=${bed} O="interval_list" SD=${dict}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/sort-mark-duplicates:2"
             cpu: "1"
             memory_gb: "4"
             queue: queue
             resource: "rusage[gtmp=10, mem=4000]"
             job_group: jobGroup
             priority: priority
             project: project
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
     String priority
     String project
     String queue

     command {
            /bin/cat ${sep=" " fastqs} > "${finalLabel}.${readNum}.fastq.gz"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
             queue: queue
             job_group: jobGroup
             priority: priority
             project: project
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
     String priority
     String project
     String queue

     command {
            (set -eo pipefail && /usr/local/bin/bwa mem -K 100000000 -t 8 -Y -R "${readGroup}" ${refFasta} ${sep=" " fastqs} | /usr/local/bin/samblaster -a --addMateTags | /usr/local/bin/samtools view -b -S /dev/stdin > "refAlign.bam")
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/tagged-alignment:2"
             cpu: "8"
             memory_gb: "20"
             queue: queue
             resource: "rusage[gtmp=10, mem=20000]"
             job_group: jobGroup
             priority: priority
             project: project
     }
     output {
            File TaggedBam = "refAlign.bam"
     }
}

task name_sort {
     String alignedBam
     String tmp
     String jobGroup
     String priority
     String project
     String queue

     command {
             /usr/local/bin/sambamba sort -t 8 -m 18G -n --tmpdir=${tmp} -o "NameSorted.bam" ${alignedBam}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/sambamba-0.6.4:1"
             cpu: "8"
             memory_gb: "20"
             queue: queue
             resource: "rusage[gtmp=10, mem=20000]"
             job_group: jobGroup
             priority: priority
             project: project
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
     String priority
     String project
     String queue

     command {
             (set -eo pipefail && /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar MarkDuplicates I=${mergedBam} O=/dev/stdout ASSUME_SORT_ORDER=queryname METRICS_FILE=mark_dups_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | /usr/local/bin/sambamba sort -t 8 -m 18G --tmpdir=${tmp} -o "${label}.bam" /dev/stdin)
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/sort-mark-duplicates:2"
             cpu: "8"
             memory_gb: "50"
             queue: queue
             resource: "rusage[gtmp=10, mem=50000]"
             job_group: jobGroup
             priority: priority
             project: project
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
             memory_gb: "8"
             queue: queue
             resource: "rusage[gtmp=10, mem=8000]"
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
     String priority
     String project
     String queue

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=${ref} INPUT=${in} OUTPUT="alignment_summary.txt" ASSUME_SORTED=true
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory_gb: "18"
             queue: queue
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup
             priority: priority
             project: project
     }
     output {
            File alignMetrics = "alignment_summary.txt"
     }
}

task collect_gc_metrics {
     String in
     String ref
     String jobGroup
     String priority
     String project
     String queue

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectGcBiasMetrics REFERENCE_SEQUENCE=${ref} INPUT=${in} OUTPUT="GC_bias.txt" SUMMARY_OUTPUT="GC_bias_summary.txt" CHART_OUTPUT="GC_bias_chart.pdf" ASSUME_SORTED=true
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory_gb: "18"
             queue: queue
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup
             priority: priority
             project: project
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
     String priority
     String project
     String queue

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectInsertSizeMetrics INPUT=${in} OUTPUT="insert_size_summary.txt" HISTOGRAM_FILE="insert_size.pdf" ASSUME_SORTED=true
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory_gb: "18"
             queue: queue
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup
             priority: priority
             project: project
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
     String priority
     String project
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
             memory_gb: "18"
             queue: queue
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup
             priority: priority
             project: project
     }
     output {
            File hsMetrics = "hs_metric_summary.txt"
	    File perTargetCoverage = "per_target_coverage.txt"
     }
     
}

task flagstat {
     String in
     String jobGroup
     String priority
     String project
     String queue

     command {
             /usr/local/bin/samtools flagstat ${in} > "flagstat.out"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/tagged-alignment:2"
             cpu: "1"
             memory_gb: "10"
             queue: queue
             resource: "rusage[gtmp=10, mem=10000]"
             job_group: jobGroup
             priority: priority
             project: project
     }
     output {
            File fsOut = "flagstat.out"
    }
}

task bamutil {
     String in
     String jobGroup
     String priority
     String project
     String queue

     command {
             /usr/local/bin/bam stats --noPhoneHome --in ${in} --phred --excludeFlags 3844 2> bamutil_stats.txt
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/bamutil:2"
             cpu: "1"
             memory_gb: "10"
             queue: queue
             resource: "rusage[gtmp=10, mem=10000]"
             job_group: jobGroup
             priority: priority
             project: project
     }
     output {
            File bamutilOut = "bamutil_stats.txt"
    }
}

task gather_result {
     String dir
     Array[String] files
     String jobGroup
     String priority
     String project
     String queue

     command {
             /bin/mv -t ${dir} ${sep=" " files}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
             queue: queue
             job_group: jobGroup
             priority: priority
             project: project
     }
     output {
            String out = stdout()
     }
}

task remove_file {
     String file
     String? order_by
     String jobGroup
     String priority
     String project
     String queue

     command {
             /bin/rm ${file}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
             queue: queue
             job_group: jobGroup
             priority: priority
             project: project
     }
     output {
            String out = stdout()
     }
}
