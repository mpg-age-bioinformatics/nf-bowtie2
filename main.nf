#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.run_type}" == "r2d2" ]] || [[ "${params.run_type}" == "raven" ]] ; 
      then
        cd ${params.image_folder}
        if [[ ! -f flexbar-3.5.sif ]] ;
          then
            singularity pull bowtie2-2.4.2.sif docker://index.docker.io/mpgagebioinformatics/bowtie2:2.4
        fi

    fi


    if [[ "${params.run_type}" == "local" ]] ; 
      then
        docker pull mpgagebioinformatics/bowtie2:2.4.2
        
    fi

    """

}

// process genome_collector {
//   stageInMode 'symlink'
//   stageOutMode 'move'

//   output:
//     val "finished", emit: get_genome_status

//   when:
//     ( ! file("${params.genomes}/${params.organism}/${params.release}/${params.organism}.${params.release}.genome").exists() )
  
//   script:
//     """
//     target_folder=/genomes/${params.organism}/${params.release}/

//     if [[ ! -e \$target_folder ]] ; then mkdir -p \$target_folder ; fi

//     cd \$target_folder

//     if [[ ! -e ${params.organism}.${params.release}.gtf ]] ; 
//       then
//         curl -#O ${params.url_gtf} && gtf=`basename ${params.url_gtf}` || gtf=`curl -#l ${params.url_gtf} | grep "gtf" | grep -v abinitio | grep -v chr` && curl -#O ${params.url_gtf}/\$gtf
//         if [[ "\$gtf" == *".gz" ]] ; then unpigz -p ${task.cpus} \$gtf ; gtf=\${gtf%.gz} ; fi
//         mv \$gtf ${params.organism}.${params.release}.gtf

//         grep -v -i 'biotype "rRNA' ${params.organism}.${params.release}.gtf | grep -v -i "Mt_rRNA" | grep -v -i srrna > ${params.organism}.${params.release}.no.rRNA.gtf
//     fi

//     if [[ ! -e ${params.organism}.${params.release}.fa ]] ; 
//       then
//         curl -#O ${params.url_dna} && dna=\$(basename ${params.url_dna} ) || dna=""
//         if [[ ! -f \$dna ]] ;
//           then 
//             dna=\$(curl -#l ${params.url_dna} | grep .dna.toplevel.fa.gz)
//             curl -#O ${params.url_dna}/\$dna
//         fi
//         if [[ "\$dna" == *".gz" ]] ; then unpigz \$dna ; dna=\${dna%.gz} ; fi
//         mv \$dna ${params.organism}.${params.release}.fa
        
//     fi

//     if [[ ! -e ${params.organism}.${params.release}.genome ]] ;
//       then
//         samtools faidx ${params.organism}.${params.release}.fa
//         awk '{print \$1"\t"\$2}' ${params.organism}.${params.release}.fa.fai > ${params.organism}.${params.release}.genome
//     fi

//     """
// }



process mapping {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    tuple val(pair_id), path(fastq)

  output:
    val pair_id

  when:
    (  ! file("${params.project_folder}/bowtie2_output/${pair_id}.sam").exists() )
    
  script:

  def single = fastq instanceof Path

  if ( single ) {
    """
        mkdir -p /workdir/bowtie2_output
        bowtie2 -p ${task.cpus} -x ${params.GENOME} -U /trimmed_raw/${pair_id}.fastq.gz -S /workdir/bowtie2_output/${pair_id}.sam
    """
  } 
  else { 
    """
      mkdir -p /workdir/bowtie2_output
      bowtie2 -p ${task.cpus} -x ${params.GENOME} -1 /trimmed_raw/${pair_id}_1.fastq.gz -2 /trimmed_raw/${pair_id}_2.fastq.gz -S /workdir/bowtie2_output/${pair_id}.sam
    """
  }
}


workflow images {
  main:
    get_images()
}

// workflow get_genome {
//   main:
//     genome_collector()
// }

workflow {
  main:
    // Channel
    //   .fromFilePairs( "${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz", size: -1 )
    //   .ifEmpty { error "Cannot find any reads matching: ${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz" }
    //   .set { read_files }
    if ( params.seq == "single" ) {
      suffix="${params.suffix_single}"
    } else {
      suffix="${params.suffix_paired}"
    }
    read_files=Channel.fromFilePairs( "${params.trimmed_raw}/*${suffix}", size: -1 )
    read_files.view()
    mapping(read_files)
}