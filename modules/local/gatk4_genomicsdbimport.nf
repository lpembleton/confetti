process GATK4_GENOMICSDBIMPORT {

    input:
    path(vcf)
    path(tbi)
    path (sample_map)
    path region
    path fasta
    path fai
    path dict

    output:
    //tuple val(meta), path(my_database)        , optional:true, emit: genomicsdb
    //path(my_database), emit: genomicsdb
    tuple env(FOO), path(my_database), emit: genomicsdb
    path "versions.yml"                                    , emit: versions
    //env FOO, emit: ids


    script:
    def args = task.ext.args ?: ''
    def avail_mem = 80
    """

    FOO=\$(echo $region | sed 's/.bed//g' | sed 's/subregion_//g')
    
    gatk BedToIntervalList -I $region -O regions.interval_list -SD *.dict
    
    gatk --java-options "-Xmx${avail_mem}g" GenomicsDBImport \\
        --sample-name-map $sample_map \\
        --batch-size 50 \\
        --genomicsdb-vcf-buffer-size 16384000 \\
        --genomicsdb-workspace-path my_database \\
        -L regions.interval_list

    chmod -R 777 my_database/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS

    """
}