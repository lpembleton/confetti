process GATK4_DICT {
    tag "$fasta"
    label 'process_single'

    input:
    path fasta

    output:
    path "*dict"        , emit: dict
    path "versions.yml" , emit: versions

    script: 
    def avail_mem = 3
    """    
    gatk --java-options "-Xmx${avail_mem}g" CreateSequenceDictionary \\
        --REFERENCE $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}