[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

# confetti
<img align="right" src="docs/images/confetti.svg" height="100">



## Introduction

**confetti** is a nextflow pipeline for processes demultiplexed restriction enzyme genotyping-by-sequencing data. It aligns reads to a reference, detects single nucleotide variants and calls genotypes.

Although the pipeline has been written for the standard single end GBS read data as described by [Elshire et al. 2011](https://doi.org/10.1371/journal.pone.0019379) and derivatives of, the mapping modules have been developed with the option for pair end reads. It just hasn't been tested yet... 


## Pipeline Summary:

1. Data Preparation:
     - Read the CSV sample sheet to retrieve the required information for each sample, including optional sequencing metadata.
     - Download the FastQ.gz files associated with each sample.

2. Sequence Read Trimming:
     - Utilize ([`fastp`](https://github.com/OpenGene/fastp)) to perform quality-based trimming and adapter removal on the sequencing reads.
     - Generate trimmed FastQ files as output.

3. Quality Control:
     - Use ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) to perform quality control analysis on the trimmed FastQ files.
     - Generate quality control reports and metrics.

4. Alignment to Reference Genome:
     - Utilize ([`BWA-MEM`](https://bio-bwa.sourceforge.net/)) to align the trimmed sequencing reads to a reference genome.
     - Generate aligned BAM files as output.
     - Calculate depth and coverage stats with ([`Samtools`](https://www.htslib.org/)) and ([`mosdepth`](https://github.com/brentp/mosdepth))

5. Variant Calling:
     - Employ GATK ([`GATK`](https://gatk.broadinstitute.org/hc/en-us)) to identify single nucleotide variants (SNVs), and perform genotype calling.
     - Generate variant call format (VCF) files containing the detected variants.

6. Post-processing and Analysis:
     - Perform additional post-processing steps, such as filtering variants based on quality, coverage, or other criteria.
     - Output the filtered variants in as a VCF file.

7. Results and Reports:
     - Generate multi-sample quality control reports using ([`MultiQC`](http://multiqc.info/)) to aggregate FastQC data and other metrics from each step.
     - Output summary reports or visualizations to provide an overview of the processed samples individuals and as an aggregate.

## Input Requirements:
The pipeline expects a CSV samplesheet as input, which should contain the sample name, sequence id, read type (single/paired) and the path to the read1 fastq file and path to the read2 fastq file (optional). It should like something similar to:

```csv
name,seqid,seq_type,fastq1
sample01,B01-W2-sample01,single,HQW5_B01-W2-sample01_S1_L001_R1_001.fastq.gz
sample02,B02-W2-sample02,single,HQW5_B02-W2-sample02_S1_L001_R1_001.fastq.gz
```
*Note the column names are important*

Multiple entries can have the same sample `name`, however `seqid` must be unique

## Usage


```bash
nextflow run main.nf \
   -profile <local/awsbatch/awsbatchspot> \
   --input samplesheet.csv \
   --reference <path-to-fasta-reference> \
   --regions <path-to-regions-bed> \
   --outPrefix <output-prefix>
```


## Pipeline Output
The pipeline outputs genetic variants and associated genotype calls in vcf and gvcf format along with a MultiQC report.

## Credits

confetti (the nf pipeline, not the party store mess) was originally written by [`LWPembleton`](https://github.com:lpembleton).

A lot of inspiration and structure was taken from the Nextflow documentation, the fantastic nf-core community and modules.

> **Nextflow enables reproducible computational workflows.**
> 
> Paolo Di Tommaso, Maria Chatzou, Evan Floden, Pablo Prieto Barja, Emilio Palumbo & Cedric Notredame.
> 
> P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:10.1038/nbt.3820

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).


## Citations


