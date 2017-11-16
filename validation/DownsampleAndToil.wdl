workflow DownsampleAndToil {

  File YmlTemplate="/gscmnt/gc2709/info/production_reference_GRCh38DH/CLE/IDTExome/detect_variants_template_IDT.yml"
  String Cwlpath="/opt/cancer-genomics-workflow/detect_variants/detect_variants.cwl" #within its docker container
  String ToilContainer = "registry.gsc.wustl.edu/genome/docker-cle-cwl:5"

  String BaitBedFile = "/gscmnt/gc2709/info/production_reference_GRCh38DH/accessory_bed_files/xgen-exome-research-panel-probes.bed"
  String TargetBedFile = "/gscmnt/gc2709/info/production_reference_GRCh38DH/accessory_bed_files/xgen-exome-research-panel-targets.bed"
  String Reference = "/gscmnt/gc2709/info/production_reference_GRCh38DH/reference/all_sequences.fa"
  String Dictionary = "/gscmnt/gc2709/info/production_reference_GRCh38DH/reference/all_sequences.dict"

  String TumorCram
  String NormalCram
  String DownsampleRatio
  String JobGroup
  String TumorLabel
  String NormalLabel
  String OutputDir

  call convert_to_bam as convert_tumor_bam {
       input: refFasta=Reference,
              cram=TumorCram,
              label=TumorLabel,
              jobGroup=JobGroup
  }

  call convert_to_bam as convert_normal_bam {
       input: refFasta=Reference,
              cram=NormalCram,
              label=NormalLabel,
              jobGroup=JobGroup
  }

  call downsample_bam as downsample_tumor {
       input: inputBam=convert_tumor_bam.bam,
              ratio=DownsampleRatio,
              label=TumorLabel,
              jobGroup=JobGroup
  }

  call downsample_bam as downsample_normal {
        input: inputBam=convert_normal_bam.bam,
               ratio=DownsampleRatio,
               label=NormalLabel,
               jobGroup=JobGroup
  }

  call remove_file as rm_tumor_bam {
       input: file=convert_tumor_bam.bam,
              order_by=downsample_tumor.outputBam,
              jobGroup=JobGroup
  }

  call remove_file as rm_normal_bam {
       input: file=convert_normal_bam.bam,
              order_by=downsample_normal.outputBam,
              jobGroup=JobGroup
  }

  call bed_to_interval_list as BaitIntervals {
       input: bed=BaitBedFile,
              dict=Dictionary,
              jobGroup=JobGroup
  }

  call bed_to_interval_list as TargetIntervals {
       input: bed=TargetBedFile,
              dict=Dictionary,
              jobGroup=JobGroup
  }

  call collect_alignment_metrics as normal_alignment_metrics {
       input: in=downsample_normal.outputBam,
              ref=Reference,
              label=NormalLabel,
              jobGroup=JobGroup
  }

  call collect_gc_metrics as normal_gc_metrics {
       input: in=downsample_normal.outputBam,
              ref=Reference,
              label=NormalLabel,
              jobGroup=JobGroup
  }

  call collect_insert_metrics as normal_insert_metrics {
       input: in=downsample_normal.outputBam,
              ref=Reference,
              label=NormalLabel,
              jobGroup=JobGroup
  }

  call collect_hs_metrics as normal_hs_metrics {
       input: in=downsample_normal.outputBam,
              ref=Reference,
              baits=BaitIntervals.intervals,
              targets=TargetIntervals.intervals,
              label=NormalLabel,
              jobGroup=JobGroup
  }

  call flagstat as normal_flagstat {
       input: in=downsample_normal.outputBam,
              label=NormalLabel,
              jobGroup=JobGroup
  }


  call collect_alignment_metrics as tumor_alignment_metrics {
       input: in=downsample_tumor.outputBam,
              ref=Reference,
              label=TumorLabel,
              jobGroup=JobGroup
  }

  call collect_gc_metrics as tumor_gc_metrics {
       input: in=downsample_tumor.outputBam,
              ref=Reference,
              label=TumorLabel,
              jobGroup=JobGroup
  }

  call collect_insert_metrics as tumor_insert_metrics {
       input: in=downsample_tumor.outputBam,
              ref=Reference,
              label=TumorLabel,
              jobGroup=JobGroup
  }

  call collect_hs_metrics as tumor_hs_metrics {
       input: in=downsample_tumor.outputBam,
              ref=Reference,
              baits=BaitIntervals.intervals,
              targets=TargetIntervals.intervals,
              label=TumorLabel,
              jobGroup=JobGroup
  }

  call flagstat as tumor_flagstat {
       input: in=downsample_tumor.outputBam,
              label=TumorLabel,
              jobGroup=JobGroup
  }


  call convert_to_cram as convert_tumor_cram {
       input: refFasta=Reference,
              bam=downsample_tumor.outputBam,
              label=TumorLabel,
              jobGroup=JobGroup
  }

  call convert_to_cram as convert_normal_cram {
       input: refFasta=Reference,
              bam=downsample_normal.outputBam,
              label=NormalLabel,
              jobGroup=JobGroup
  }

  call prep_yaml {
       input: tumor=convert_tumor_cram.cram,
              normal=convert_normal_cram.cram,
              template=YmlTemplate,
              jobGroup=JobGroup
  }

  call run_toil {
       input: 
             Yml=prep_yaml.InputYml,
             Cwlpath=Cwlpath,
             Rundir=OutputDir,
             jobGroup=JobGroup,
             container=ToilContainer
  }

  call gather_result as gather_the_rest {
       input: files=[run_toil.final_report,
                     convert_normal_cram.cram,
                     convert_normal_cram.cram_index,
                     convert_tumor_cram.cram,
                     convert_tumor_cram.cram_index,
                     normal_alignment_metrics.alignMetrics,
                     normal_gc_metrics.gcOut,
                     normal_gc_metrics.gcSum,
                     normal_gc_metrics.gcPDF,
                     normal_insert_metrics.isOut,
                     normal_insert_metrics.isPDF,
                     normal_hs_metrics.hsMetrics,
                     normal_hs_metrics.perTargetCoverage,
                     normal_flagstat.fsOut,
                     tumor_alignment_metrics.alignMetrics,
                     tumor_gc_metrics.gcOut,
                     tumor_gc_metrics.gcSum,
                     tumor_gc_metrics.gcPDF,
                     tumor_insert_metrics.isOut,
                     tumor_insert_metrics.isPDF,
                     tumor_hs_metrics.hsMetrics,
                     tumor_hs_metrics.perTargetCoverage,
                     tumor_flagstat.fsOut
                     ],
              dir=OutputDir,
              jobGroup=JobGroup
  }

  call remove_file as rm_down_tumor_bam {
       input: file=downsample_tumor.outputBam,
       order_by=gather_the_rest.out,
       jobGroup=JobGroup
  }
 
  call remove_file as rm_down_normal_bam {
       input: file=downsample_normal.outputBam,
       order_by=gather_the_rest.out,
       jobGroup=JobGroup
  }
}

task bed_to_interval_list {
     String bed
     String dict
     String jobGroup

     command {
           /usr/bin/java -Xmx4g -jar /usr/picard/picard.jar BedToIntervalList I=${bed} O="interval_list" SD=${dict}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/sort-mark-duplicates:2"
             cpu: "1"
             memory_gb: "4"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=4000]"
             job_group: jobGroup
     }
     output {
         File intervals = "interval_list"
     }
}

task convert_to_bam {
     String refFasta
     String cram
     String label
     String jobGroup

     command {
             /usr/local/bin/samtools view -b -T ${refFasta} ${cram} > "${label}.bam"; /usr/local/bin/samtools index "${label}.bam"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/samtools-1.3.1-2:2"
             cpu: "1"
             memory_gb: "20"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=20000]"
             job_group: jobGroup
     }
     output {
            File bam = "${label}.bam"
            File bam_index = "${label}.bam.bai"
    }
}


task convert_to_cram {
     String refFasta
     String bam
     String label
     String jobGroup

     command {
             /usr/local/bin/samtools view -C -T ${refFasta} ${bam} > "${label}.cram" && \
             /usr/local/bin/samtools index "${label}.cram" && \
             /usr/local/bin/samtools index "${label}.cram" "${label}.crai"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/samtools-1.3.1-2:2"
             cpu: "1"
             memory_gb: "20"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=20000]"
             job_group: jobGroup
     }
     output {
            File cram = "${label}.cram"
            File cram_index = "${label}.cram.crai"
    }

}

task downsample_bam {
     String inputBam
     String ratio
     String label
     String jobGroup
     
     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar DownsampleSam INPUT=${inputBam} OUTPUT="${label}.downsample.bam" P=${ratio}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory_gb: "18"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup 
     }
     output {
            File outputBam = "${label}.downsample.bam"
     }
}

task collect_alignment_metrics {
     String in
     String ref
     String label
     String jobGroup

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=${ref} INPUT=${in} OUTPUT="${label}.alignment_summary.txt" ASSUME_SORTED=true
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory_gb: "18"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup 
     }
     output {
            File alignMetrics = "${label}.alignment_summary.txt"
     }
}

task collect_gc_metrics {
     String in
     String ref
     String label
     String jobGroup

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectGcBiasMetrics REFERENCE_SEQUENCE=${ref} INPUT=${in} OUTPUT="${label}.GC_bias.txt" SUMMARY_OUTPUT="${label}.GC_bias_summary.txt" CHART_OUTPUT="${label}.GC_bias_chart.pdf" ASSUME_SORTED=true
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory_gb: "18"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup 
     }
     output {
            File gcOut = "${label}.GC_bias.txt"
            File gcSum = "${label}.GC_bias_summary.txt"
            File gcPDF = "${label}.GC_bias_chart.pdf"
     }     
}

task collect_insert_metrics {
     String in
     String ref
     String label
     String jobGroup

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectInsertSizeMetrics INPUT=${in} OUTPUT="${label}.insert_size_summary.txt" HISTOGRAM_FILE="${label}.insert_size.pdf" ASSUME_SORTED=true
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory_gb: "18"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup 
     }
     output {
            File isOut = "${label}.insert_size_summary.txt"
            File isPDF = "${label}.insert_size.pdf"
     }
     
}
task collect_hs_metrics {
     String in
     String ref
     String baits
     String targets
     String label
     String jobGroup

     command {
             /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectHsMetrics REFERENCE_SEQUENCE=${ref} INPUT=${in} \
             BAIT_INTERVALS=${baits} TARGET_INTERVALS=${targets} \
             MINIMUM_MAPPING_QUALITY=1 MINIMUM_BASE_QUALITY=1 \
             OUTPUT="${label}.hs_metric_summary.txt" PER_TARGET_COVERAGE="${label}.per_target_coverage.txt"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
             cpu: "2"
             memory_gb: "18"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=18000]"
             job_group: jobGroup 
     }
     output {
            File hsMetrics = "${label}.hs_metric_summary.txt"
            File perTargetCoverage = "${label}.per_target_coverage.txt"
     }
     
}

task flagstat {
     String in
     String label
     String jobGroup

     command {
             /usr/local/bin/samtools flagstat ${in} > "${label}.flagstat.out"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/tagged-alignment:2"
             cpu: "1"
             memory_gb: "10"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=10000]"
             job_group: jobGroup 
     }
     output {
            File fsOut = "${label}.flagstat.out"
    }
}

task prep_yaml {
     String tumor
     String normal
     String template
     String jobGroup

     command {
             /usr/bin/printf 'tumor_cram:%s\n  class: File%s\n  path: ${tumor}%s\nnormal_cram:%s\n  class: File%s\n  path: ${normal}%s\n' >> "detect_variant.yml";
             /bin/cat ${template} >> "detect_variant.yml"
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
             queue: "research-hpc"
             job_group: jobGroup
     }
     output {
            File InputYml = "detect_variant.yml"
     }
}

task run_toil {
    File Yml
    String Cwlpath
    String Rundir
    String jobGroup
    String container

    command {
        /bin/mkdir ${Rundir}/toil_outdir; /bin/mkdir ${Rundir}/workDir; \
        export LSB_SUB_ADDITIONAL="docker(${container})" && \
        export LSB_DEFAULTQUEUE=research-hpc && \
        cwltoil --disableCaching --logLevel=DEBUG --logFile=${Rundir}/toil.log --outdir=${Rundir}/toil_outdir --workDir ${Rundir}/workDir --jobStore ${Rundir}/jobStore --batchSystem lsf ${Cwlpath} ${Yml}
    }
    runtime {
             docker_image: container
             cpu: "1"
             memory_gb: "16"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=16000]"
             job_group: jobGroup
     }
     output {
            File final_report = "${Rundir}/toil_outdir/variants.annotated.tsv"
     }
}

task gather_result {
     String dir
     Array[String] files
     String jobGroup

     command {
             /bin/mv -t ${dir} ${sep=" " files}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
             queue: "research-hpc"
             job_group: jobGroup 
     }
     output {
            String out = stdout()
     }
}

task remove_file {
     String file
     String order_by
     String jobGroup

     command {
             /bin/rm ${file}
     }
     runtime {
             docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
             queue: "research-hpc"
             job_group: jobGroup 
     }
     output {
            String out = stdout()
     }
}


