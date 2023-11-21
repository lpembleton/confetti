process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*bam"), path("*.bai") , optional:true, emit: bam
    path  "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ 1 \\
        $args \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}