/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
*/

//
// MODULE & SUBWORKFLOW: Installed directly from nf-core/modules & nf-core/subworkflows
//

include { VSEARCH_CLUSTER } from '../modules/nf-core/vsearch/cluster/main'

//
// MODULE: Installed directly from nf-core/modules
//

include { RENAME_RAW_DATA_FILES } from '../modules/local/rename_raw_data_files'
include { DADA2_ERR } from '../modules/local/dada2_err'
include { NOVASEQ_ERR } from '../modules/local/novaseq_err'
include { DADA2_DENOISING } from '../modules/local/dada2_denoising'
include { DADA2_RMCHIMERA } from '../modules/local/dada2_rmchimera'
include { DADA2_STATS } from '../modules/local/dada2_stats'
include { DADA2_MERGE } from '../modules/local/dada2_merge'
include { FILTER_LEN as FILTER_LEN_ASV } from '../modules/local/filter_len'
include { MERGE_STATS as MERGE_STATS_FILTERLENASV } from '../modules/local/merge_stats'
include { MERGE_STATS as MERGE_STATS_STD } from '../modules/local/merge_stats'
include { FILTER_CLUSTERS } from '../modules/local/filter_clusters'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { PARSE_INPUT } from '../subworkflows/local/parse_input'
include { DADA2_PREPROCESSING } from '../subworkflows/local/dada2_preprocessing'
include { CUTADAPT_WORKFLOW } from '../subworkflows/local/cutadapt_workflow'

//
// FUNCTIONS
//
include { samplesheetToList } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { makeComplement } from '../subworkflows/local/utils_nfcore_ampliseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
*/

workflow AMPLISEQ {

    main:
    //
    // INPUT AND VARIABLES
    //
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    single_end = params.single_end
    if (params.pacbio || params.iontorrent) {
        single_end = true
    }

    trunclenf = params.trunclenf ?: 0
    trunclenr = params.trunclenr ?: 0
    if ( !params.single_end && !params.illumina_pe_its && (params.trunclenf == null || params.trunclenr == null) && !params.input_fasta ) {
        find_truncation_values = true
        log.warn "No DADA2 read truncation cutoffs were specified (`--trunclenf` & `--trunclenr`), therefore reads will be truncated where median quality drops below ${params.trunc_qmin} (defined by `--trunc_qmin`) but at least a fraction of ${params.trunc_rmin} (defined by `--trunc_rmin`) of the reads will be retained.\n The chosen cutoffs do not account for required overlap for merging, therefore DADA2 might have poor merging efficiency or even fail.\nThe cutoffs are chosen before any quality score-based read truncation (using `--truncq`) is performed.\n"
    } else { find_truncation_values = false }

    //
    // Create input channels
    //
    ch_input_fasta = Channel.empty()
    ch_input_reads = Channel.empty()
    if ( params.input ) {
        // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
        ch_input_reads = Channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")) // meta: meta.sample, meta.run
            .map{ meta, readfw, readrv ->
                meta.single_end = single_end.toBoolean()
                def reads = single_end ? readfw : [readfw,readrv]
                if ( !meta.single_end && !readrv ) { error("Entry `reverseReads` is missing in $params.input for $meta.sample, either correct the samplesheet or use `--single_end`, `--pacbio`, or `--iontorrent" ) } // make sure that reverse reads are present when single_end isn't specified
                if ( !meta.single_end && ( readfw.getSimpleName() == meta.sample || readrv.getSimpleName() == meta.sample ) ) { error("Entry `sampleID` cannot be identical to simple name of `forwardReads` or `reverseReads`, please change `sampleID` in $params.input for sample $meta.sample") } // sample name and any file name without extensions aren't identical, because rename_raw_data_files.nf would forward 3 files (2 renamed +1 input) instead of 2 in that case
                if ( meta.single_end && ( readfw.getSimpleName() == meta.sample+"_1" || readfw.getSimpleName() == meta.sample+"_2" ) ) { error("Entry `sampleID`+ `_1` or `_2` cannot be identical to simple name of `forwardReads`, please change `sampleID` in $params.input for sample $meta.sample") } // sample name and file name without extensions aren't identical, because rename_raw_data_files.nf would forward 2 files (1 renamed +1 input) instead of 1 in that case
                return [meta, reads]
            }
    } else if ( params.input_folder ) {
        PARSE_INPUT ( params.input_folder, single_end, params.multiple_sequencing_runs, params.extension )
        ch_input_reads = PARSE_INPUT.out.reads
    } else {
        error("One of `--input` or `--input_folder` must be provided!")
    }

    //
    // Add primer info to sequencing files
    //
    ch_input_reads
        .map{ info, reads ->
            def meta = info +
                [region: null, region_length: null] +
                [fw_primer: params.FW_primer, rv_primer: params.RV_primer] +
                [id: info.sample] +
                [fw_primer_revcomp: params.FW_primer ? makeComplement(params.FW_primer.reverse()) : null] +
                [rv_primer_revcomp: params.RV_primer ? makeComplement(params.RV_primer.reverse()) : null]
            return [ meta, reads ]
        }
        .set { ch_input_reads }


    //Filter empty files
    ch_input_reads.dump(tag:'ch_input_reads')
        .branch { it ->
            failed: it[0].single_end ? it[1].countFastq() < params.min_read_counts : it[1][0].countFastq() < params.min_read_counts || it[1][1].countFastq() < params.min_read_counts
            passed: true
        }
        .set { ch_reads_result }
    ch_reads_result.passed.set { ch_reads }
    ch_reads_result.failed
        .map { meta, reads -> [ meta.id ] }
        .collect()
        .subscribe { it ->
            def samples = it.join("\n")
            if (params.ignore_empty_input_files) {
                log.warn "At least one input file for the following sample(s) had too few reads (<$params.min_read_counts):\n$samples\nThe threshold can be adjusted with `--min_read_counts`. Ignoring failed samples and continue!\n"
            } else {
                error("At least one input file for the following sample(s) had too few reads (<$params.min_read_counts):\n$samples\nEither remove those samples, adjust the threshold with `--min_read_counts`, or ignore that samples using `--ignore_empty_input_files`.")
            }
        }
    ch_reads.dump(tag: 'ch_reads')

    //
    // MODULE: Rename files
    //
    RENAME_RAW_DATA_FILES ( ch_reads )
    ch_versions = ch_versions.mix(RENAME_RAW_DATA_FILES.out.versions.first())

    //
    // MODULE: Cutadapt
    //
    CUTADAPT_WORKFLOW (
        RENAME_RAW_DATA_FILES.out.fastq,
        params.illumina_pe_its,
        params.double_primer
    ).reads.set { ch_trimmed_reads }
    ch_versions = ch_versions.mix(CUTADAPT_WORKFLOW.out.versions)

    //
    // SUBWORKFLOW: Read preprocessing & QC plotting with DADA2
    //
    DADA2_PREPROCESSING (
        ch_trimmed_reads,
        single_end,
        find_truncation_values,
        trunclenf,
        trunclenr
    ).reads.set { ch_filt_reads }
    ch_versions = ch_versions.mix(DADA2_PREPROCESSING.out.versions)

    //
    // MODULES: ASV generation with DADA2
    //

    //run error model
    if ( !params.illumina_novaseq ) {
        DADA2_ERR ( ch_filt_reads )
        ch_errormodel = DADA2_ERR.out.errormodel
        ch_versions = ch_versions.mix(DADA2_ERR.out.versions)
    } else {
        DADA2_ERR ( ch_filt_reads )
        NOVASEQ_ERR ( DADA2_ERR.out.errormodel )
        ch_errormodel = NOVASEQ_ERR.out.errormodel
        ch_versions = ch_versions.mix(DADA2_ERR.out.versions)
    }

    //group by meta
    ch_filt_reads
        .join( ch_errormodel )
        .set { ch_derep_errormodel }
    DADA2_DENOISING ( ch_derep_errormodel.dump(tag: 'into_denoising')  )
    ch_versions = ch_versions.mix(DADA2_DENOISING.out.versions)

    DADA2_RMCHIMERA ( DADA2_DENOISING.out.seqtab )
    ch_versions = ch_versions.mix(DADA2_RMCHIMERA.out.versions)

    //group by sequencing run & group by meta
    DADA2_PREPROCESSING.out.logs
        .join( DADA2_DENOISING.out.denoised )
        .join( DADA2_DENOISING.out.mergers )
        .join( DADA2_RMCHIMERA.out.rds )
        .set { ch_track_numbers }
    DADA2_STATS ( ch_track_numbers )
    ch_versions = ch_versions.mix(DADA2_STATS.out.versions)

    //merge if several runs, otherwise just publish
    DADA2_MERGE (
        DADA2_STATS.out.stats.map { meta, stats -> stats }.collect(),
        DADA2_RMCHIMERA.out.rds.map { meta, rds -> rds }.collect() )
    ch_versions = ch_versions.mix(DADA2_MERGE.out.versions)

    //merge cutadapt_summary and dada_stats files
    MERGE_STATS_STD (CUTADAPT_WORKFLOW.out.summary, DADA2_MERGE.out.dada2stats)
    ch_stats = MERGE_STATS_STD.out.tsv
    ch_versions = ch_versions.mix(MERGE_STATS_STD.out.versions)

    ch_dada2_fasta = DADA2_MERGE.out.fasta
    ch_dada2_asv = DADA2_MERGE.out.asv

    //
    // MODULE : ASV post-clustering with VSEARCH
    //
    if (params.vsearch_cluster) {
        ch_fasta_for_clustering = ch_dada2_fasta
            .map {
                fasta ->
                    def meta = [: ]
                    meta.id = "ASV_post_clustering"
                    [ meta, fasta ]
            }
        VSEARCH_CLUSTER ( ch_fasta_for_clustering )
        ch_versions = ch_versions.mix(VSEARCH_CLUSTER.out.versions)
        FILTER_CLUSTERS ( VSEARCH_CLUSTER.out.clusters, ch_dada2_asv )
        ch_versions = ch_versions.mix(FILTER_CLUSTERS.out.versions)
        ch_dada2_fasta = FILTER_CLUSTERS.out.fasta
        ch_dada2_asv = FILTER_CLUSTERS.out.asv
    }

    ch_unfiltered_fasta = ch_dada2_fasta

    ch_barrnapsummary = Channel.empty()
    ch_dada2_fasta = ch_unfiltered_fasta

    //
    // Modules : amplicon length filtering
    //
    if ( (params.min_len_asv || params.max_len_asv) ) {
        FILTER_LEN_ASV ( ch_dada2_fasta, ch_dada2_asv.ifEmpty( [] ) )
        ch_versions = ch_versions.mix(FILTER_LEN_ASV.out.versions)
        MERGE_STATS_FILTERLENASV ( ch_stats, FILTER_LEN_ASV.out.stats )
        ch_stats = MERGE_STATS_FILTERLENASV.out.tsv
        ch_dada2_fasta = FILTER_LEN_ASV.out.fasta
        ch_dada2_asv = FILTER_LEN_ASV.out.asv
        // Make sure that not all sequences were removed
        ch_dada2_fasta.subscribe { it -> if (it.countLines() == 0) error("ASV length filtering activated by '--min_len_asv' or '--max_len_asv' removed all ASVs, please adjust settings.") }
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    multiqc_report = []      // empty channel
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
*/
