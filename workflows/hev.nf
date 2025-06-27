/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BLAST_BLASTN } from '../modules/nf-core/blast/blastn/main'
include { BLAST_MAKEBLASTDB } from '../modules/nf-core/blast/makeblastdb/main'
include { BLASTPARSE } from '../modules/local/blastparse/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FLYE } from '../modules/nf-core/flye/main'
include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/kraken2/kraken2/main'
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PRINSEQPLUSPLUS        } from '../modules/nf-core/prinseqplusplus/main' 
include { SAMTOOLS_BAM2FQ        } from '../modules/nf-core/samtools/bam2fq/main' 
include { SAMTOOLS_BAM2FQ as SAMTOOLS_BAM2FQ_MAPPED        } from '../modules/nf-core/samtools/bam2fq/main' 
include { SPADES } from '../modules/nf-core/spades/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_hev_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HEV {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Make BLAST database from HEV reference genomes
    //
    BLAST_MAKEBLASTDB (
        [ [id:"blast_db"], file(params.references) ] // Add empty meta.id map before the references file
    )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions.first())

    //
    // MODULE: Run Samtools bam2fq
    //
    SAMTOOLS_BAM2FQ (
        ch_samplesheet,
        false // Do not split reads into pairs
    )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        SAMTOOLS_BAM2FQ.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    
    //
    // MODULE: Run Prinseq++ to filter fastq reads on length
    //
    PRINSEQPLUSPLUS (
        SAMTOOLS_BAM2FQ.out.reads
    )
    ch_versions = ch_versions.mix(PRINSEQPLUSPLUS.out.versions.first())

    //
    // MODULE: Run Kraken2 to classify reads
    //
    KRAKEN2_KRAKEN2(
        PRINSEQPLUSPLUS.out.reads,
        file(params.kraken_all_db),
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

    //
    // MODULE: Map reads to a set of HEV reference genomes to keep mapped reads
    //
    MINIMAP2_ALIGN(
        PRINSEQPLUSPLUS.out.good_reads,
        [[], file(params.references)],
        true, // Save as bam file
        'bai', // Index extension for bam files
        false,
        false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    //
    // MODULE: Convert mapped reads to fastq
    //
    SAMTOOLS_BAM2FQ_MAPPED (
        MINIMAP2_ALIGN.out.bam,
        false // Do not split reads into pairs
    )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ_MAPPED.out.versions.first())

    //
    // MODULE: Run Spades to assemble reads
    //
    FLYE (
        SAMTOOLS_BAM2FQ_MAPPED.out.reads,
        "--nano-raw" // Flye mode for nanopore reads
    )
    ch_versions = ch_versions.mix(FLYE.out.versions.first())

    //
    // MODULE: Blast assembled output against HEV reference genomes
    //
    BLAST_BLASTN (
        FLYE.out.fasta,
        BLAST_MAKEBLASTDB.out.db
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    //
    // MODULE: Parse the blast output
    //
    BLASTPARSE (
        BLAST_BLASTN.out.txt.join(FLYE.out.fasta), // tuple val(meta), path(blast_out), path(scaffolds)
        file(params.references), // path to the same references used in the Blastn
    )
    ch_versions = ch_versions.mix(BLASTPARSE.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'hev_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
