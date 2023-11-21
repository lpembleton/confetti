process GATK4_HAPLOTYPECALLER {
    tag "$meta.name"
    label 'process_multi'

    publishDir "$params.outdir/variant_calling/gvcfs", pattern: '*gvcf.gz*', mode: 'copy'

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path('*gvcf.gz')       , emit: gvcf
    tuple val(meta), path('*gvcf.gz.tbi')   , emit: tbi
    path "versions.yml"                     , emit: versions

    script:
    def avail_mem = 20
    def prefix = task.ext.prefix ?: "${meta.name}"
    """
    echo "${bam.join('\n')}" > bam.list
    gatk --java-options "-Xmx${avail_mem}g" HaplotypeCaller \\
        -I bam.list \\
        --output ${prefix}.gvcf.gz \\
        -ERC GVCF \\
        --reference ${fasta}
    
    gatk IndexFeatureFile \\
        --input ${prefix}.gvcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}