process GATK4_FILTERVARIANTS {
    label 'process_single'

    publishDir "$params.outdir/variant_calling", pattern: '*.vcf.gz', mode: 'copy'
    publishDir "$params.outdir/variant_calling", pattern: '*.tbi', mode: 'copy'

    input:
    path vcf
    path tbi
    path fasta
    path fai
    path dict

    output:
    tuple path('*split-multi.snps.vcf.gz'), path('*split-multi.snps.vcf.gz.tbi')          , emit: snp_vcf
    tuple path('*split-multi.snps.filtered.vcf.gz'), path('*split-multi.snps.filtered.vcf.gz.tbi'), emit: filtered_snp_vcf
    path "versions.yml"                                                                     , emit: versions

    script:
    """
    
    gatk LeftAlignAndTrimVariants \\
        -R ${fasta} \\
        -V ${vcf} \\
        -O ${params.outPrefix}.split-multi.vcf.gz \\
        --split-multi-allelics

    gatk SelectVariants \\
        -R ${fasta} \\
        -V ${params.outPrefix}.split-multi.vcf.gz \\
        --select-type-to-include SNP \\
        -O ${params.outPrefix}.split-multi.snps.vcf.gz

    gatk VariantFiltration \\
        -R ${fasta} \\
        -V ${params.outPrefix}.split-multi.snps.vcf.gz \\
        -filter "QD < 2.0" --filter-name "QD2" \\
        -filter "QUAL < 30.0" --filter-name "QUAL30" \\
        -filter "MQ < 40.0" --filter-name "MQ40" \\
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        -O ${params.outPrefix}.split-multi.snps.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}