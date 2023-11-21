process BCFTOOLS_SORT {
    tag "$id"
    label 'process_medium'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("*sorted.vcf.gz")  , emit: vcf

    script:


    """
    bcftools \\
        sort \\
        --output ${params.outPrefix}.${id}.split-multi.sorted.vcf.gz \\
        --temp-dir . \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

}