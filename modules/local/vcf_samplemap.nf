process VCF_SAMPLEMAP {
    label 'process_single'

    input:
    path(bam)

    output:
    path('*.sample_map'),   emit: sample_map
    path "versions.yml",    emit: versions

    script:
    """

    touch cohort.sample_map
    for f in *.gvcf.gz; do
        echo -e "\$(vcf-query -l \$f)\t\$f" >> cohort.sample_map
    done
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/(VCFtools)//' | sed 's/[()]//g')
    END_VERSIONS
    """

}