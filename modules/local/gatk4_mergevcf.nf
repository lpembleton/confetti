process GATK4_MERGEVCFS {
    label 'process_medium'

    publishDir "$params.outdir/variant_calling", pattern: '*.vcf.gz', mode: 'copy'
    publishDir "$params.outdir/variant_calling", pattern: '*.tbi', mode: 'copy'


    input:
    path(vcf)
    path(dict)

    output:
    path('*.vcf.gz'), emit: vcf
    path('*.tbi'), emit: tbi
    path  "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def input_list = vcf.collect{ "--INPUT $it"}.join(' ')
    def reference_command = dict ? "--SEQUENCE_DICTIONARY $dict" : ""
    def avail_mem = 20

    """
    gatk --java-options "-Xmx${avail_mem}g" MergeVcfs \\
        $input_list \\
        --OUTPUT ${params.outPrefix}.vcf.gz \\
        $reference_command \\
        --TMP_DIR . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}