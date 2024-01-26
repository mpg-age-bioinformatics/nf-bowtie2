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
        if [[ ! -f bowtie2-2.4.2.sif ]] ;
          then
            singularity pull bowtie2-2.4.2.sif docker://index.docker.io/mpgagebioinformatics/bowtie2:2.4.2
        fi
        if [[ ! -f samtools-1.16.1.sif ]] ;
          then
            singularity pull samtools-1.16.1.sif docker://index.docker.io/mpgagebioinformatics/samtools:1.16.1
        fi
        if [[ ! -f picard-3.0.0.sif ]] ;
          then
            singularity pull picard-3.0.0.sif docker://index.docker.io/mpgagebioinformatics/picard:3.0.0
        fi
    fi


    if [[ "${params.run_type}" == "local" ]] ; 
      then
        docker pull mpgagebioinformatics/bowtie2:2.4.2
        docker pull mpgagebioinformatics/samtools:1.16.1  
        docker pull mpgagebioinformatics/picard:3.0.0
    fi

    """

}
process bowtie2_indexer {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.genomes}/${params.organism}/${params.release}/toplevel_bowtie2/index.fa").exists() ) 
  
  script:
    """

    target_folder=${params.genomes}/${params.organism}/${params.release}/toplevel_bowtie2

    if [[ ! -e \$target_folder ]] ; then mkdir -p \$target_folder ; fi

    cd \$target_folder

    ln -s ../${params.organism}.${params.release}.fa index.fa
    
    bowtie2-build index.fa index.fa
    """
}


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
        bowtie2 -p ${task.cpus} -x ${params.genomes}/${params.organism}/${params.release}/toplevel_bowtie2/index.fa -U /trimmed_raw/${pair_id}.fastq.gz -S /workdir/bowtie2_output/${pair_id}.sam
    """
  } 
  else { 
    """
      mkdir -p /workdir/bowtie2_output
      bowtie2 -p ${task.cpus} -x ${params.genomes}/${params.organism}/${params.release}/toplevel_bowtie2/index.fa -1 /trimmed_raw/${pair_id}_1.fastq.gz -2 /trimmed_raw/${pair_id}_2.fastq.gz -S /workdir/bowtie2_output/${pair_id}.sam
    """
  }
}


process remove_mito {
  stageInMode 'symlink'
  stageOutMode 'move'
  
  input:
    tuple val(pair_id), path(fastq)

  output:
    val pair_id

  when:
      (  ! file("${params.project_folder}/bowtie2_output/${pair_id}.mito.txt").exists() )

  script:
    def single = fastq instanceof Path

    """
    echo ${single}

    cd /workdir/bowtie2_output/

    if [ ${single} = true ]; then

      if [ ${params.remove_mito} = "yes" ] ; then
        awk '{if (\$3 != "MT" ) print }' ${pair_id}.sam | samtools view -@ 10 -q 10 -bS - > ${pair_id}.bam
      else
        samtools view -@ 10 -q 10 -f 2 -bS ${pair_id}.sam > ${pair_id}.bam
      fi
    
    else

      if [ ${params.remove_mito} = "yes" ] ; then
        awk '{if (\$3 != "MT" ) print }' ${pair_id}.sam | samtools view -@ 10 -q 10 -f 2 -bS - > ${pair_id}.bam
      else
        samtools view -@ 10 -q 10 -bS ${pair_id}.sam > ${pair_id}.bam
      fi
      
    fi

samtools sort -@ 10 -o ${pair_id}.ss.bam ${pair_id}.bam
samtools index ${pair_id}.ss.bam
awk '{if (\$3 != "MT" ) print }' ${pair_id}.sam | wc -l > ${pair_id}.mito.txt
awk '{if (\$3 == "MT" ) print }' ${pair_id}.sam | wc -l >> ${pair_id}.mito.txt

    """
}  

    
process remove_duplicates {
  stageInMode 'symlink'
  stageOutMode 'move'
  
  input:
    val input_file

  when:
      (  ! file("${input_file}".replaceAll('.ss.bam', '.md_metrics.txt')).exists() )

  script:
  """
  cd /workdir/bowtie2_output/
  mkdir -p /workdir/tmp

  input_file_name=\$(basename ${input_file} .ss.bam)

  java -Xmx16g -jar /PICARD/picard.jar MarkDuplicates I=${input_file}\
	O=\${input_file_name}.md.bam CREATE_INDEX=true \
	TMP_DIR=/workdir/tmp M=\${input_file_name}.md_metrics.txt REMOVE_DUPLICATES=true
  """

}

process samtools_flagstat {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val input_file

  when:
    ( ! file("${input_file}".replaceAll('.bam', '.flagstat.txt')).exists() ) 

  script:
    """
    cd /workdir/bowtie2_output/
    input_file_name=\$(basename ${input_file} .bam)
    samtools flagstat ${input_file} > \${input_file_name}.flagstat.txt
    """
} 


process qc_count {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val input_file

  when:
    ( ! file("${input_file}".replaceAll('.sam', '.count.txt')).exists() ) 

  script:
    """
    cd /workdir/bowtie2_output/
    input_file_name=\$(basename ${input_file} .sam)
    samtools view -@ 10 -c ${input_file} > \${input_file_name}.count.txt
    samtools view -@ 10 -q 10 -c ${input_file} >> \${input_file_name}.count.txt
    samtools view -@ 10 -q 10 -f 2 -c ${input_file} >> \${input_file_name}.count.txt
    """
} 



workflow images {
  main:
    get_images()
}

workflow index {
  main:
    bowtie2_indexer()
}


workflow align{
  main:
    // Channel
    //   .fromFilePairs( "${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz", size: -1 )
    //   .ifEmpty { error "Cannot find any reads matching: ${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz" }
    //   .set { read_files }
    if ( params.seq == "single" ) {

      // if ( 'suffix_single' in params.keySet() ) {
        suffix=params.suffix_single
      // } else {
      //   suffix=".fastq.gz"
      // }

    } else {
      // if ( 'suffix_paired' in params.keySet() ) {
        suffix=params.suffix_paired
      // } else {
      //   suffix="_{1,2}.fastq.gz"
      // }
    }
    
    read_files=Channel.fromFilePairs( "${params.bowtie_raw_data}/*${suffix}", size: -1 )
    // read_files.view()
    mapping(read_files)
    // flagstat( mapping.out.collect(), read_files )
}

workflow mito {
  main:
      if ( params.seq == "single" ) {
        suffix=params.suffix_single
    } else {
        suffix=params.suffix_paired
    }
    read_files=Channel.fromFilePairs( "${params.bowtie_raw_data}/*${suffix}", size: -1 )
    remove_mito(read_files)
}

workflow picard {
  main:
    data=Channel.fromPath( "${params.project_folder}/bowtie2_output/*.ss.bam" )
    remove_duplicates(data)
}

workflow flagstat {
  main:
    data=Channel.fromPath( ["${params.project_folder}/bowtie2_output/*.bam"] )
    samtools_flagstat(data)
    
}

workflow qccount {
    data=Channel.fromPath( ["${params.project_folder}/bowtie2_output/*.sam"] )
    qc_count(data)
}
