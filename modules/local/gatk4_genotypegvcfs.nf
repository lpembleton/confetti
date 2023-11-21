process GATK4_GENOTYPEGVCFS {
    tag "genotyping $id"
    label 'process_high_single'

    input:
    path gvcf
    path tbi
    path fasta
    path fai
    path dict
    tuple val(id), path(genomicsdb)

    output:
    tuple val(id), path("${params.outPrefix}.${id}.vcf.gz")                 , emit: vcf
    tuple val(id), path("${params.outPrefix}.${id}.vcf.gz.tbi")             , emit: vcfidx
    path  "versions.yml"                            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 22
    """

    #chmod -R 777 my_database/

    gatk --java-options "-Xmx${avail_mem}g" GenotypeGVCFs \\
        -R ${fasta} \\
        -V gendb://${genomicsdb} \\
        --allow-old-rms-mapping-quality-annotation-data \\
        -O ${params.outPrefix}.${id}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS

    """
}