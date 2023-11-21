/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check mandatory parameters
if (params.input) { csv_file = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

 
if (params.reference == null) error "Please specify a reference genome with --reference"
// demo reference  can be found in docs/example-data/Lolium_2.6.1_Chromosomes_Scaffolds_chr2_0-10M.fasta



// print parameters


log.info """\
    ======================================================================
    C O N F E T T I
    A   R E S T R I C T I O N   E N Z Y M E   G B S   P I P E L I N E
    ======================================================================
    samplesheet: ${params.input}
    reference: ${params.reference}
    output directory: ${params.outdir}
    genomic regions: ${params.regions}
    output prefix: ${params.outPrefix}
    ======================================================================
    """
    .stripIndent()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP } from '../modules/local/fastp'
include { FASTQC } from '../modules/local/fastqc'
include { SAMTOOLS_FAIDX } from '../modules/local/samtools_faidx'
include { SAMTOOLS_INDEX } from '../modules/local/samtools_index'
include { BWA_INDEX } from '../modules/local/bwa_index'
include { BWA_MEM } from '../modules/local/bwa_mem'
include { SAMTOOLS_STATS } from '../modules/local/samtools_stats'
include { MOSDEPTH } from '../modules/local/mosdepth'
include { MULTIQC } from '../modules/local/multiqc'
include { GATK4_DICT } from '../modules/local/gatk4_dict'
include { GATK4_HAPLOTYPECALLER } from '../modules/local/gatk4_haplotypecaller'
include { GATK4_GENOTYPEGVCFS } from '../modules/local/gatk4_genotypegvcfs'
include { GATK4_FILTERVARIANTS } from '../modules/local/gatk4_filtervariants'
include { GATK4_GENOMICSDBIMPORT } from '../modules/local/gatk4_genomicsdbimport'
include { VCF_SAMPLEMAP } from '../modules/local/vcf_samplemap'
include { AWK_SPLITBED } from '../modules/local/awk_splitbed'
include { GATK4_MERGEVCFS } from '../modules/local/gatk4_mergevcf'
include { BCFTOOLS_SORT } from '../modules/local/bcftools_sort'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/local/custom/custom_dumpsoftwareversions'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow CONFETTI {
    // To gather all QC reports for MultiQC
    versions = Channel.empty()
    // To gather used softwares versions for MultiQC
    reports = Channel.empty()

    input_sample = extract_csv(file(csv_file))

    // PREPROCESSING

    // Read quality and adapter trimming
    //FASTP(input_sample)
    // Gather used softwares versions and reports
    //versions = versions.mix(FASTP.out.versions)
    //reports = reports.mix(FASTP.out.json.collect{ meta, json -> json })
    //reports = reports.mix(FASTP.out.html.collect{ meta, html -> html })
 

    // QC on trimmed reads
    //FASTQC(FASTP.out.reads)
    FASTQC(input_sample)
    // Gather used softwares versions and reports
    versions = versions.mix(FASTQC.out.versions)
    reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })

    // MAPPING

    // Prepare reference fasta index
    SAMTOOLS_FAIDX(params.reference)
    BWA_INDEX(params.reference)
    // Gather used softwares versions and reports
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    versions = versions.mix(BWA_INDEX.out.versions)


    // Map reads to reference fasta
    sort = true // sort after mapping
    //BWA_MEM(FASTP.out.reads, BWA_INDEX.out.index, sort)
    BWA_MEM(input_sample, BWA_INDEX.out.index, sort)
    // Gather used softwares versions and reports
    versions = versions.mix(BWA_MEM.out.versions)


    // get bam stats
    SAMTOOLS_STATS(BWA_MEM.out.bam, BWA_MEM.out.bai)
    // Gather used softwares versions and reports
    versions = versions.mix(SAMTOOLS_STATS.out.versions)
    reports = reports.mix(SAMTOOLS_STATS.out.stats)


    // WGS coverage stats
    MOSDEPTH(BWA_MEM.out.bam, BWA_MEM.out.bai)
    // Gather used softwares versions and reports
    versions = versions.mix(MOSDEPTH.out.versions)    
    reports = reports.mix(MOSDEPTH.out.global_txt)
    reports = reports.mix(MOSDEPTH.out.summary_txt)

    // VARIANT CALLING

    // Create GATK reference sequence dictionary
    GATK4_DICT(params.reference)

    // GATK4 variant calling
    GATK4_HAPLOTYPECALLER(BWA_MEM.out.bam.map{ meta, bam -> [ meta.subMap('name'), bam ] }.groupTuple(), BWA_MEM.out.bai.map{ meta, bai -> [ meta.subMap('name'), bai ] }.groupTuple(), params.reference, SAMTOOLS_FAIDX.out.fai, GATK4_DICT.out.dict)

    // GATK4 genotype caller
    // lines_from_file = Channel.fromPath(params.regions).splitText()

    AWK_SPLITBED(params.regions)


    VCF_SAMPLEMAP(GATK4_HAPLOTYPECALLER.out.gvcf.map { tuple -> tuple[1] }.collect())

    GATK4_GENOMICSDBIMPORT(GATK4_HAPLOTYPECALLER.out.gvcf.map { tuple -> tuple[1] }.collect(), GATK4_HAPLOTYPECALLER.out.tbi.map { tuple -> tuple[1] }.collect(), VCF_SAMPLEMAP.out.sample_map, AWK_SPLITBED.out.subregion.flatten(), params.reference, SAMTOOLS_FAIDX.out.fai, GATK4_DICT.out.dict)
 
    GATK4_GENOTYPEGVCFS(GATK4_HAPLOTYPECALLER.out.gvcf.map { tuple -> tuple[1] }.collect(), GATK4_HAPLOTYPECALLER.out.tbi.map { tuple -> tuple[1] }.collect(), params.reference, SAMTOOLS_FAIDX.out.fai, GATK4_DICT.out.dict, GATK4_GENOMICSDBIMPORT.out.genomicsdb)

    BCFTOOLS_SORT(GATK4_GENOTYPEGVCFS.out.vcf)

    GATK4_MERGEVCFS(BCFTOOLS_SORT.out.vcf.map { tuple -> tuple[1] }.collect(), GATK4_DICT.out.dict)

    // GATK4 basic variant filter (GATK standards)
    GATK4_FILTERVARIANTS(GATK4_MERGEVCFS.out.vcf, GATK4_MERGEVCFS.out.tbi, params.reference, SAMTOOLS_FAIDX.out.fai, GATK4_DICT.out.dict)


    // REPORTING
    version_yaml = Channel.empty()
    CUSTOM_DUMPSOFTWAREVERSIONS(versions.unique().collectFile(name: 'collated_versions.yml'))
    version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()

    // Multiqc
    //multiqc_files = Channel.empty()
    //multiqc_files = multiqc_files.mix(reports.collect())
    //MULTIQC(multiqc_files.collect())

    multiqc_files = Channel.empty()
    multiqc_files = multiqc_files.mix(version_yaml)
    multiqc_files = multiqc_files.mix(reports.collect().ifEmpty([]))

    MULTIQC(multiqc_files.collect())

    multiqc_report = MULTIQC.out.report.toList()
    versions = versions.mix(MULTIQC.out.versions)
    

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Function to extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file) {

    Channel.of(csv_file).splitCsv(header: true)
        // Retrieves number of lanes by grouping together by name and id and counting how many entries there are for this combination
        .map{ row ->
            //sample_count_all++
            if (!(row.name && row.seqid)) {
                error("Missing field in csv file header. The csv file must have fields named 'name' and 'seqid'.")
            }
            else if (row.name.contains(" ") || row.seqid.contains(" ")) {
                error("Invalid value in csv file. Values for 'name' and 'seqid' can not contain space.")
            }
            [ [ row.name.toString(), row.seqId.toString() ], row ]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [ rows, size ]
        }.transpose()
        .map{ row, num_lanes -> // from here do the usual thing for csv parsing

        def meta = [:]

        if (row.name) meta.name = row.name.toString()
        if (row.seqId)  meta.seqid  = row.seqid.toString()

        if (row.seq_type) meta.seq_type = row.seq_type.toString()
        else meta.seq_type = 'paired'
       

        meta.id         = "${row.name}_X_${row.seqid}".toString()
        def fastq_1     = file(row.fastq_1, checkIfExists: true)
        if(row.fastq_2) fastq_2     = file(row.fastq_2, checkIfExists: true)

        def flowcell    = flowcellLaneFromFastq(fastq_1)
        // Don't use a random element for ID, it breaks resuming
        // the below read_group needs to be revised once all meta params are accessible
        def read_group  = "\"@RG\\tID:${flowcell}.${meta.seqid}\\tSM:${meta.name}\\tDS:${params.reference}\""

        //meta.num_lanes  = num_lanes.toInteger()
        meta.read_group = read_group.toString()
        meta.data_type  = 'fastq'

        if (row.fastq_2) return [ meta, [fastq_1, fastq_2] ]
        else {
            return [ meta, [fastq_1] ]
        }
    }
}


// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
