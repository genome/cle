workflow RunVep {
    File InputVcf
    String Name="test"
    String JobGroup="/fdu/haloplextest"
    String OutputDir="/gscuser/dspencer/projects/wdltest/haloplex/vep_test/output"
    String Reference="/gscuser/dspencer/refdata/GRCh37/all_sequences.fa"

    call run_vep {
        input: 
            CombineVcf=InputVcf,
            Reference=Reference,
            Name=Name,
            jobGroup=JobGroup
    }

    call bgzip {
        input:
            file=run_vep.vep_output,
            jobGroup=JobGroup
    }

    call index {
        input:
            file=bgzip.gzip_vcf,
            jobGroup=JobGroup
    }

    call variants_to_table {
        input:
            reference=Reference,
            input_vcf=index.index_vcf,
            Name=Name,
            jobGroup=JobGroup
    }

    call add_vep_annotation_to_table {
        input:
            input_tsv=variants_to_table.out_tsv,
            input_vcf=index.index_vcf,
            Name=Name,
            jobGroup=JobGroup
    }

    call gather_files {
        input: 
            OutputFiles=[run_vep.vep_output,add_vep_annotation_to_table.annotated_tsv],
            OutputDir=OutputDir,
            jobGroup=JobGroup
    }
}

task run_vep {
    File CombineVcf
    String Reference
    String Name
    String jobGroup

    command {
        /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
        --format vcf --vcf --plugin Downstream --plugin Wildtype --symbol --term SO --flag_pick -o ${Name}.vep_annotation.out \
        -i ${CombineVcf} --hgvs --fasta ${Reference} --offline --cache --af_exac --dir /gscmnt/gc2764/cad/fdu/fdu_scratch/ensembl_vep87_GRCh37
    }
    runtime {
             docker_image: "mgibio/cle"
             cpu: "1"
             memory_gb: "10"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=10000]"
             job_group: jobGroup
     }
     output {
            File vep_output = "${Name}.vep_annotation.out"
     }
}

task bgzip {
    String file
    String jobGroup

    command {
        /opt/htslib/bin/bgzip -c ${file} > ${file}.vcf.gz
    }
    runtime {
             docker_image: "mgibio/cle"
             queue: "research-hpc"
             job_group: jobGroup
     }
     output {
            File gzip_vcf = "${file}.vcf.gz"
     }
}

task index {
    String file
    String jobGroup

    command {
        /usr/bin/tabix -p vcf ${file}
    }
    runtime {
             docker_image: "mgibio/cle"
             queue: "research-hpc"
             job_group: jobGroup
    }
    output {
           File index_vcf = "${file}"
    }
}

task variants_to_table {
    String reference
    String input_vcf
    String Name
    String jobGroup

    command {
        /usr/bin/java -Xmx4g -jar /opt/GenomeAnalysisTK.jar -T VariantsToTable \
        -R ${reference} --variant ${input_vcf} -o ${Name}.variants.tsv \
        -F CHROM -F POS -F ID -F REF -F ALT -F set -F AC -F AF \
        -GF GT -GF AD -GF NR -GF NV -GF TAMP -GF SAMP -GF VAFTYPE -GF AMPS -GF CVAF
    }
    runtime {
             docker_image: "mgibio/cle"
             cpu: "1"
             memory_gb: "10"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=10000]"
             job_group: jobGroup
     }
     output {
            File out_tsv = "${Name}.variants.tsv"
     }
}

task add_vep_annotation_to_table {
    String input_tsv
    String input_vcf
    String Name
    String jobGroup

    command {
        /usr/bin/python /usr/bin/add_annotations_to_table_helper.py ${input_tsv} ${input_vcf} Consequence,SYMBOL,Feature_type,Feature,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,HGNC_ID,ExAC_AF,ExAC_Adj_AF,ExAC_AFR_AF,ExAC_AMR_AF,ExAC_EAS_AF,ExAC_FIN_AF,ExAC_NFE_AF,ExAC_OTH_AF,ExAC_SAS_AF,CLIN_SIG,SOMATIC,PHENO ./
    }
    runtime {
             docker_image: "mgibio/cle"
             cpu: "1"
             memory_gb: "10"
             queue: "research-hpc"
             resource: "rusage[gtmp=10, mem=10000]"
             job_group: jobGroup
     }
     output {
            File annotated_tsv = "variants.annotated.tsv"
     }
}

task gather_files {
     Array[String] OutputFiles
     String OutputDir
     String? SubDir
     String jobGroup

     command {
             /bin/mv -f -t ${OutputDir}/ ${sep=" " OutputFiles}
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


