workflow FastqToTargetedAlignedBam {

  String BaitBedFile = "/gscmnt/gc2709/info/production_reference_GRCh38DH/CLE/IDTExome/xgen-exome-research-panel-probes.bed"
  String TargetBedFile = "/gscmnt/gc2709/info/production_reference_GRCh38DH/CLE/IDTExome/xgen-exome-research-panel-targets.bed"
  
  String Reference = "/gscmnt/gc2709/info/production_reference_GRCh38DH/reference/all_sequences.fa"
  String Dictionary = "/gscmnt/gc2709/info/production_reference_GRCh38DH/reference/all_sequences.dict"

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

  call fastqc as fastqc_read1 {
       input: in=cat_fastqs_read1.Fastq,
              readNum="1",
              finalLabel=FinalLabel,
              jobGroup=JobGroup,
              queue=Queue
  }

  call fastqc as fastqc_read2 {
       input: in=cat_fastqs_read2.Fastq,
              readNum="2",
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

  call collect_alignment_metrics {
       input: in=mark.SortedBam,
              finalLabel=FinalLabel,
              ref=Reference,
              jobGroup=JobGroup,
              queue=Queue
  }

  call collect_insert_metrics {
       input: in=mark.SortedBam,
              finalLabel=FinalLabel,
              ref=Reference,
              jobGroup=JobGroup,
              queue=Queue
  }

  call collect_hs_metrics {
       input: in=mark.SortedBam,
              finalLabel=FinalLabel,
              ref=Reference,
	      baits=BaitIntervals.intervals,
	      targets=TargetIntervals.intervals,
              jobGroup=JobGroup,
              queue=Queue
  }

  call flagstat {
       input: in=mark.SortedBam,
              finalLabel=FinalLabel,
              jobGroup=JobGroup,
              queue=Queue
  }

  call mapQ  {
       input: in=mark.SortedBam,
              label=FinalLabel,
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
                     collect_insert_metrics.isOut,
                     collect_insert_metrics.isPDF,
		     collect_hs_metrics.hsMetrics,
		     collect_hs_metrics.perTargetCoverage,
                     flagstat.fsOut,
                     fastqc_read1.fastqcHTML,
                     fastqc_read1.fastqcDATA,
                     fastqc_read2.fastqcHTML,
                     fastqc_read2.fastqcDATA,
                     mapQ.mapqOUT,
                     mark.MarkMetrics,
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
             memory_gb: "4"
             queue: queue
             resource: "rusage[gtmp=10, mem=4000]"
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
            /bin/cat ${sep=" " fastqs} > "${finalLabel}.read${readNum}.fastq.gz"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
             queue: queue
             job_group: jobGroup
     }
     output {
            File Fastq = "${finalLabel}.read${readNum}.fastq.gz"
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
             memory_gb: "20"
             queue: queue
             resource: "rusage[gtmp=10, mem=20000]"
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
             memory_gb: "20"
             queue: queue
             resource: "rusage[gtmp=10, mem=20000]"
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
             (set -eo pipefail && /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar MarkDuplicates I=${mergedBam} O=/dev/stdout ASSUME_SORT_ORDER=queryname METRICS_FILE=mark_dups_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | /usr/local/bin/sambamba sort -t 8 -m 18G --tmpdir=${tmp} -o "${label}.bam" /dev/stdin) && \
             /bin/sed -i "s/NameSorted.bam/${label}.bam/" mark_dups_metrics.txt
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/sort-mark-duplicates:2"
             cpu: "8"
             memory_gb: "50"
             queue: queue
             resource: "rusage[gtmp=10, mem=50000]"
             job_group: jobGroup
     }
     output {
             File SortedBam = "${label}.bam"
             File SortedBamIndex = "${label}.bam.bai"
             File MarkMetrics = "mark_dups_metrics.txt"
     }
     
}

task collect_alignment_metrics {
     String in
     String finalLabel
     String ref
     String jobGroup
     String queue

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=${ref} INPUT=${in} OUTPUT="${finalLabel}.alignment_summary.txt" ASSUME_SORTED=true
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory_gb: "18"
             queue: queue
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup
     }
     output {
            File alignMetrics = "alignment_summary.txt"
     }
}

task collect_insert_metrics {
     String in
     String finalLabel
     String ref
     String jobGroup
     String queue

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectInsertSizeMetrics INPUT=${in} OUTPUT="${finalLabel}.insert_size_summary.txt" HISTOGRAM_FILE="${finalLabel}.insert_size.pdf" ASSUME_SORTED=true
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory_gb: "18"
             queue: queue
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup
     }
     output {
            File isOut = "insert_size_summary.txt"
            File isPDF = "insert_size.pdf"
     }
     
}

task collect_hs_metrics {
     String in
     String finalLabel
     String ref
     String baits
     String targets
     String jobGroup
     String queue

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectHsMetrics REFERENCE_SEQUENCE=${ref} INPUT=${in} \
	     		   BAIT_INTERVALS=${baits} TARGET_INTERVALS=${targets} \
			   MINIMUM_MAPPING_QUALITY=1 MINIMUM_BASE_QUALITY=1 \
			   OUTPUT="${finalLabel}.hs_metric_summary.txt" PER_TARGET_COVERAGE="${finalLabel}.per_target_coverage.txt"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory_gb: "18"
             queue: queue
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup
     }
     output {
            File hsMetrics = "hs_metric_summary.txt"
	    File perTargetCoverage = "per_target_coverage.txt"
     }
     
}

task flagstat {
     String in
     String finalLabel
     String jobGroup
     String queue

     command {
             /usr/local/bin/samtools flagstat ${in} > "${finalLabel}.flagstat.out"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/tagged-alignment:2"
             cpu: "1"
             memory_gb: "10"
             queue: queue
             resource: "rusage[gtmp=10, mem=10000]"
             job_group: jobGroup
     }
     output {
            File fsOut = "${finalLabel}.flagstat.out"
    }
}

task fastqc {
     String in
     String readNum
     String finalLabel
     String jobGroup
     String queue

     command {
             /usr/local/bin/fastqc --extract --outdir ./ ${in} && \
             /bin/mv "${finalLabel}.read${readNum}_fastqc/fastqc_data.txt" ${finalLabel}.read${readNum}_fastqc_data.txt
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/fdu/fastqc:1"
             cpu: "1"
             memory_gb: "16"
             queue: queue
             resource: "rusage[gtmp=10, mem=16000]"
             job_group: jobGroup
     }
     output {
            File fastqcHTML = "${finalLabel}.read${readNum}_fastqc.html"
            File fastqcDATA = "${finalLabel}.read${readNum}_fastqc_data.txt"
    }
}

task mapQ {
     String in
     String label
     String jobGroup
     String queue
     
     command <<<
            /usr/local/bin/samtools view ${in}|/usr/bin/awk -v q="'" '{s+=$5; c++;} END {print "# title: ",q"Average Mapping Quality"q; print "# description: ",q"Average Mapping Quality"q;print "# section: ",q"Custom Data File"q; print "# format: ",q"tsv"q; print "# plot_type: ",q"bargraph"q; print "# pconfig:"; print "#    id:", q"custom_bargraph_w_header"q; print "#    ylab: ",q"Average MapQ"q; print "MapQ_Average",s/c;}' > "${label}_mapQ_mqc.txt"
             
     >>>
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/tagged-alignment:2"
             cpu: "1"
             memory_gb: "16"
             queue: queue
             resource: "rusage[gtmp=10, mem=16000]"
             job_group: jobGroup
     }
     output {
            File mapqOUT = "${label}_mapQ_mqc.txt"
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
