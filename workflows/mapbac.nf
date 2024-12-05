/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mapbac_pipeline'
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.reference, params.multiqc_config,
                           params.shortread_qc_adapterlist, params.multiqc_logo, 
                           params.multiqc_methods_description ]
                            
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if ( params.input ) {
    ch_input = file(params.input, checkIfExists: true)
} else {
    error("Input samplesheet not specified")
}

if (params.reference) { ch_reference =  Channel.fromPath(params.reference) } else { exit 1, 'Reference sequence FASTA not specified!' }

// Modify reference channel to include meta data
ch_reference_meta = ch_reference.map{ it -> [[id:it[0].baseName], it] }.collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { FASTQSCANPARSE as FASTQSCANPARSE_TRIM      } from '../modules/local/fastqscanparse/main'
include { FASTQSCANPARSE as FASTQSCANPARSE_SUBSAMPLE } from '../modules/local/fastqscanparse/main'
include { SEQTK_COMP                                 } from '../modules/local/seqtk_comp/main'
include { SEQTK_PARSE                                } from '../modules/local/seqtk_parse/main'
include { ALIGNPSEUDOGENOMES                         } from '../modules/local/alignpseudogenomes/main'

include { SHORTREAD_PREPROCESSING                    } from '../subworkflows/local/shortread_preprocessing'
include { LONGREAD_PREPROCESSING                     } from '../subworkflows/local/longread_preprocessing'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// NF-CORE MODULES/PLUGINS
//
include { BWAMEM2_INDEX                          } from '../modules/nf-core/bwamem2/index/main'
include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { CAT_FASTQ as MERGE_RUNS                } from '../modules/nf-core/cat/fastq/main'
include { FASTQSCAN as FASTQSCAN_TRIM            } from '../modules/nf-core/modules/fastqscan/main'
include { RASUSA                                 } from '../modules/nf-core/rasusa/main'
include { FASTQSCAN as FASTQSCAN_SUBSAMPLE       } from '../modules/nf-core/modules/fastqscan/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SHORT } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_LONG  } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX                         } from '../modules/nf-core/samtools/index/main'
include { PICARD_COLLECTMULTIPLEMETRICS          } from '../modules/nf-core/picard/collectmultiplemetrics/main'
include { SAMTOOLS_FAIDX                         } from '../modules/nf-core/samtools/faidx/main'
include { BCFTOOLS_FILTER                        } from '../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_CONSENSUS                     } from '../modules/nf-core/bcftools/consensus/main'
include { SNPSITES                               } from '../modules/nf-core/snpsites/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'

//
// NF-CORE SUBWORKFLOWS
//
include { BAM_MARKDUPLICATES_PICARD                   } from '../subworkflows/nf-core/bam_markduplicates_picard/main'
include { BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS } from '../subworkflows/nf-core/bam_variant_calling_sort_freebayes_bcftools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAPBAC {

    take:
    samplesheet  // channel: samplesheet read in from --input
    ch_reference // channel: path(reference.fasta)

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    // Validate input files and create separate channels for FASTQ, FASTA, and Nanopore data
    ch_input = samplesheet
        .map { meta, run_accession, instrument_platform, fastq_1, fastq_2, fasta ->
            meta.run_accession = run_accession
            meta.instrument_platform = instrument_platform

            // Define single_end based on the conditions
            meta.single_end = ( fastq_1 && !fastq_2 && instrument_platform != 'OXFORD_NANOPORE' )

            // Define is_fasta based on the presence of fasta
            meta.is_fasta = fasta ? true : false

            if ( !meta.is_fasta && !fastq_1 ) {
                error("ERROR: Please check input samplesheet: entry `fastq_1` doesn't exist!")
            }
            if ( meta.instrument_platform == 'OXFORD_NANOPORE' && fastq_2 ) {
                error("Error: Please check input samplesheet: for Oxford Nanopore reads entry `fastq_2` should be empty!")
            }
            if ( meta.single_end && fastq_2 ) {
                error("Error: Please check input samplesheet: for single-end reads entry `fastq_2` should be empty")
            }
            return [ meta, run_accession, instrument_platform, fastq_1, fastq_2, fasta ]
        }
        .branch { meta, run_accession, instrument_platform, fastq_1, fastq_2, fasta ->
            fastq: meta.single_end || fastq_2
                return [ meta + [ type: "short" ], fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]
            nanopore: instrument_platform == 'OXFORD_NANOPORE' && !meta.is_fasta
                meta.single_end = true
                return [ meta + [ type: "long" ], [ fastq_1 ] ]
            fasta_short: meta.is_fasta && instrument_platform == 'ILLUMINA'
                meta.single_end = true
                return [ meta + [ type: "short" ], [ fasta ] ]
            fasta_long: meta.is_fasta && instrument_platform == 'OXFORD_NANOPORE'
                meta.single_end = true
                return [ meta + [ type: "long" ], [ fasta ] ]
        }

    // Merge ch_input.fastq and ch_input.nanopore into a single channel
    ch_input_for_fastqc = ch_input.fastq.mix( ch_input.nanopore )
    
    /*
        MODULE: Run Minimap2 index on reference
    */

    BWAMEM2_INDEX (
        ch_reference_meta
    )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.first())
     
    /*
        MODULE: Run FastQC
    */
    if ( !params.skip_preprocessing_qc ) {
    
        FASTQC (
            ch_input_for_fastqc
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    /*
        SUBWORKFLOW: PERFORM PREPROCESSING
    */

    if ( params.perform_shortread_qc ) {
        ch_shortreads_preprocessed = SHORTREAD_PREPROCESSING ( ch_input.fastq, adapterlist ).reads
        ch_versions = ch_versions.mix( SHORTREAD_PREPROCESSING.out.versions )
    } else {
        ch_shortreads_preprocessed = ch_input.fastq
    }

    if ( params.perform_longread_qc ) {
        ch_longreads_preprocessed = LONGREAD_PREPROCESSING ( ch_input.nanopore ).reads
                                        .map { it -> [ it[0], [it[1]] ] }
        ch_versions = ch_versions.mix( LONGREAD_PREPROCESSING.out.versions )
    } else {
        ch_longreads_preprocessed = ch_input.nanopore
    }
    
    if ( params.perform_runmerging ) {

        ch_reads_for_cat_branch = ch_shortreads_preprocessed
            .mix( ch_longreads_preprocessed )
            .map {
                meta, reads ->
                    def meta_new = meta - meta.subMap('run_accession')
                    [ meta_new, reads ]
            }
            .groupTuple()
            .map {
                meta, reads ->
                    [ meta, reads.flatten() ]
            }
            .branch {
                meta, reads ->
                // we can't concatenate files if there is not a second run, we branch
                // here to separate them out, and mix back in after for efficiency
                cat: ( meta.single_end && reads.size() > 1 ) || ( !meta.single_end && reads.size() > 2 )
                skip: true
            }

        ch_reads_runmerged = MERGE_RUNS ( ch_reads_for_cat_branch.cat ).reads
            .mix( ch_reads_for_cat_branch.skip )
            .map {
                meta, reads ->
                [ meta, [ reads ].flatten() ]
            }
            .mix( ch_input.fasta_short, ch_input.fasta_long)

        ch_versions = ch_versions.mix(MERGE_RUNS.out.versions)

    } else {
        ch_reads_runmerged = ch_shortreads_preprocessed
            .mix( ch_longreads_preprocessed, ch_input.fasta_short, ch_input.fasta_long )
    }
    
    /*
        MODULE: Run fastq-scan
    */
    FASTQSCAN_TRIM (
        ch_reads_runmerged
    )
    ch_fastqscantrim_fastqscanparse = FASTQSCAN_RAW.out.json
    ch_fastqscantrim_readstats      = FASTQSCAN_RAW.out.json
    ch_versions                     = ch_versions.mix(FASTQSCAN_TRIM.out.versions.first())
    
    /*
        MODULE: Run fastqscanparse
    */
    FASTQSCANPARSE_TRIM (
        ch_fastqscantrim_fastqscanparse.collect{it[1]}.ifEmpty([])
    )
    ch_versions = ch_versions.mix(FASTQSCANPARSE_TRIM.out.versions.first())
    
    /*
        MODULE: PERFORM SUBSAMPLING
    */
     
    if ( params.perform_subsampling ) {
        ch_reads_subsampled = RASUSA( ch_reads_runmerged, params.genome_size, subsampling_depth_cutoff ).reads
        ch_versions = ch_versions.mix( RASUSA.out.versions )
    } else {
        ch_reads_subsampled = ch_reads_runmerged
    }
    
    /*
        MODULE: Run fastq-scan
    */
    FASTQSCAN_SUBSAMPLE (
        ch_reads_subsampled
    )
    ch_fastqscansubsample_fastqscanparse = FASTQSCAN_SUBSAMPLE.out.json
    ch_fastqscansubsample_readstats      = FASTQSCAN_SUBSAMPLE.out.json
    ch_versions                          = ch_versions.mix(FASTQSCAN_SUBSAMPLE.out.versions.first())
    
    /*
        MODULE: Run fastqscanparse
    */
    FASTQSCANPARSE_SUBSAMPLE (
        ch_fastqscansubsample_fastqscanparse.collect{it[1]}.ifEmpty([])
    )
    ch_versions = ch_versions.mix(FASTQSCANPARSE_SUBSAMPLE.out.versions.first())
    
    /*
        MODULE: Map reads
    */
       
    ch_shortreads_bam = MINIMAP2_ALIGN_SHORT ( ch_reads_subsampled.fastq, ch_reference_meta, true, "bai", false, false ).bam
    ch_versions = ch_versions.mix( MINIMAP2_ALIGN_SHORT.out.versions )

    ch_longreads_bam = MINIMAP2_ALIGN_LONG ( ch_reads_subsampled.nanopore, ch_reference_meta, true, "bai", false, false ).bam
                                        .map { it -> [ it[0], [it[1]] ] }
    ch_versions = ch_versions.mix( MINIMAP2_ALIGN_LONG.out.versions )

    ch_bam = ch_shortreads_bam
        .mix( ch_longreads_bam )

    /*
        SUBWORKFLOW: Mark duplicate reads
    */

    if ( params.perform_markduplicates ) {
        ch_bam_deduplicated = BAM_MARKDUPLICATES_PICARD (
            ch_bam,
            ch_reference_meta,
            BWAMEM2_INDEX.index
        ).bam
        ch_versions = ch_versions.mix( BAM_MARKDUPLICATES_PICARD.out.versions )
    } else {
        ch_bam_deduplicated = ch_bam
    }
    
    /*
        MODULE: Samtools index
    */
    SAMTOOLS_INDEX (
        ch_bam_deduplicated
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
      
    /*
        MODULE: Picard metrics
    */
    PICARD_COLLECTMULTIPLEMETRICS (
        ch_bam,
        ch_reference_meta,
        []
    )
    ch_versions = ch_versions.mix( PICARD_COLLECTMULTIPLEMETRICS.out.versions )

    /*
        MODULE: Calculate read stats
    */
    ch_fastqscantrim_readstats                          // tuple val(meta), path(json)
        .join( FASTQSCAN_SUBSAMPLE.out.json )           // tuple val(meta), path(json) 
        .join( BAM_MARKDUPLICATES_PICARD.out.depth )    // tuple val(meta), path(depth)
        .join( BAM_MARKDUPLICATES_PICARD.out.mapreads ) // tuple val(meta), path(mapreads)
        .set { ch_readstats }                           // tuple val(meta), path(json), path(json), path(depth), path(mapreads)

    READ_STATS (
        ch_readstats
    )
    ch_readstats_readstatsparse = READ_STATS.out.csv
    ch_versions                 = ch_versions.mix(READ_STATS.out.versions.first())

    /*
        MODULE: Summarise read stats outputs
    */
    READSTATS_PARSE (
        ch_readstats_readstatsparse.collect{it[1]}.ifEmpty([])
    )
    ch_versions = ch_versions.mix(READSTATS_PARSE.out.versions.first())

    /*
        MODULE: Index reference file with Samtools faidx
    */
    
    SAMTOOLS_FAIDX (
        ch_reference_meta,
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    /*
        SUBWORKFLOW: Call variants with FreeBayes
    */
    freebayes_input = PICARD_MARKDUPLICATES.out.bam     // channel: [ val(meta), path(bam) ]
        .join(SAMTOOLS_INDEX.out.bai)                   // channel: [ val(meta), path(bam), path(bam_index)]
            .multiMap{
                meta, bam, bai ->
                   reads: [ meta, bam, bai, [], [], [] ]
            }

    ch_reference_meta                                   // channel: [ val(meta), path(reference), path(fai)]
        .join(SAMTOOLS_FAIDX.fai)
        .set { ch_freebayes_ref }                       
       
    BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS (freebayes_input.reads,
                        ch_freebayes_ref,
                        [],
                        [],
                        []
    )
    ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.versions.first())

    /*
        MODULE: Filter variants with bcftools filter
    */
    BCFTOOLS_FILTER ( 
        BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.vcf  
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())
    
    /*
        MODULE: Create consensus sequence with bcftools consensus
    */
    consensus_mask = []
    
    bcftools_consensus_input = BCFTOOLS_FILTER.out.vcf // channel: [ val(meta), path(tbi) ]
        .join(BCFTOOLS_FILTER.out.tbi)                 // channel: [ val(meta), path(vcf), path(tbi)]
        .mix(ch_reference)                             // channel: [ val(meta), path(vcf), path(tbi), path(fasta)]
        .join(consensus_mask)                          // channel: [ val(meta), path(vcf), path(tbi), path(fasta), path(mask)]

    BCFTOOLS_CONSENSUS ( 
        bcftools_consensus_input  
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())

    /*
        MODULE: Calculate number of mapped positions in pseudogenome
    */
    SEQTK_COMP (
        BCFTOOLS_CONSENSUS.out.fasta
    )
    ch_seqtk_seqtkparse = SEQTK_COMP.out.tsv
    ch_versions = ch_versions.mix(SEQTK_COMP.out.versions.first())

    /*
        MODULE: Summarise seqtk outputs
    */
    SEQTK_PARSE (
        ch_seqtk_seqtkparse.collect{it[1]}.ifEmpty([])
    )
    ch_seqtk_metadata = SEQTK_PARSE.out.tsv
    ch_versions       = ch_versions.mix(SEQTK_PARSE.out.versions.first())
    
    /*
        MODULE: Make pseudogenome alignment
    */
    ALIGNPSEUDOGENOMES (
        BCFTOOLS_CONSENSUS.out.fasta.map { fasta -> fasta[1] }.collect(),
        ch_reference
    )
    ch_versions = ch_versions.mix(ALIGNPSEUDOGENOMES.out.versions.first())

    ALIGNPSEUDOGENOMES.out.aligned_pseudogenomes
        .branch {
            aligned_pseudogenomes ->
            ALIGNMENT_NUM_PASS: aligned_pseudogenomes[0].toInteger() >= 4
            ALIGNMENT_NUM_FAIL: aligned_pseudogenomes[0].toInteger() < 4
        }
        .set { aligned_pseudogenomes_branch }

    // Don't proceeed further if two few genonmes
    aligned_pseudogenomes_branch.ALIGNMENT_NUM_FAIL.view { "Insufficient (${it[0]}) genomes after filtering to continue. Check results/pseudogenomes/low_quality_pseudogenomes.tsv for details"}

    aligned_pseudogenomes_branch.ALIGNMENT_NUM_PASS
        .map{ it[1] }
        .set { aligned_pseudogenomes }
    
    /*
        MODULE: Extract SNPs from masked alignment
    */
    SNPSITES (
        aligned_pseudogenomes
    )
    ch_versions = ch_versions.mix(SNPSITES.out.versions.first())
    
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
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

    if ( !params.skip_preprocessing_qc ) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    }

    if (params.perform_shortread_qc) {
        ch_multiqc_files = ch_multiqc_files.mix( SHORTREAD_PREPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    if (params.perform_longread_qc) {
        ch_multiqc_files = ch_multiqc_files.mix( LONGREAD_PREPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    if (params.perform_markduplicates) {
        ch_multiqc_files = ch_multiqc_files.mix( BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]) )
        ch_multiqc_files = ch_multiqc_files.mix( BAM_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]) )
    }

    ch_multiqc_files = ch_multiqc_files.mix( PICARD_COLLECTMULTIPLEMETRICS.out.metrics.collect{it[1]}.ifEmpty([]) )
    
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
