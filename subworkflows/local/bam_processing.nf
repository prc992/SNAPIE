nextflow.enable.dsl=2

include {sort_bam} from '../../modules/local/sort_bam'
include {lib_complex_preseq} from '../../modules/local/lib_complex_preseq'
include {unique_sam} from '../../modules/local/unique_sam'
include {quality_filter} from '../../modules/local/quality_filter'
include {filter_properly_paired} from '../../modules/local/filter_properly_paired'
include {createStatsSamtoolsfiltered} from '../../modules/local/createStatsSamtoolsfiltered'
include {dedup} from '../../modules/local/dedup'
include {dac_exclusion} from '../../modules/local/dac_exclusion'
include {index_sam} from '../../modules/local/index_sam'
include {createSMaSHFingerPrint} from '../../modules/local/snp_smash_fingerprint'
include {createSMaSHFingerPrintPlot} from '../../modules/local/snp_smash_fingerprint'
include {multiqc} from '../../modules/local/multiqc'
include {moveSoftFiles} from '../../modules/local/moveSoftFiles'
include {alignment_stats_report} from '../../modules/local/alignment_stats_report'

workflow BAM_PROCESSING {

    take:
    chSampleInfo
    chGenome
    chGenomeIndex
    chChromSizes
    chDACFileRef
    chSNPS_ref
    chAlign
    chTrimmingCounts
    chSkipAlignment
    chSNPSMaSH
    chSNPSMaSHPyPlot
    chMultiQCConfig
    chFilesReportInitialization
    chInitReport
    chFilesReportAlignment
    chAlignmentReport

 
    main:
    chDACFilteredFiles = Channel.of("NO_DATA")
    chIndexFiles = Channel.of("NO_DATA")
    chSNPSMaSHPlot = Channel.of("NO_DATA")
    chLibComplexPreseq = Channel.of("NO_DATA")
    chFilesReportBamProcessing = Channel.of("NO_DATA")
    chBAMProcessReport = Channel.of("NO_DATA")

    if (params.deduped_bam) {

        chDedup = chAlign.map { sampleId,enrichment_mark,control,read_method,bam,alignYml -> 
            tuple(sampleId,enrichment_mark, control,read_method,bam, null, null)}

        chSortBam = Channel.of("NO_DATA")
        chFilterPP = Channel.of("NO_DATA")
        chLibComplexPreseq = Channel.of("NO_DATA")
        chUniqueSam = Channel.of("NO_DATA")
        chFilterQuality = Channel.of("NO_DATA")
        chFilteredFiles = Channel.of("NO_DATA")
    }
    else{
        chSortBam = sort_bam(chAlign)
        chFilterPP = filter_properly_paired(chSortBam)
        chLibComplexPreseq = lib_complex_preseq(chFilterPP)
        chUniqueSam = unique_sam(chFilterPP)
        chFilterQuality = quality_filter(chUniqueSam)
        chDedup = dedup(chFilterQuality)
    }

    // Filter the DAC files
    if (params.exclude_dac_regions) {
        chDACFilteredFiles = dac_exclusion(chDedup,chDACFileRef)
        //chDACFilteredFiles = chDedup
    } else {
        chDACFilteredFiles = chDedup
    }

    if (!params.deduped_bam && !chSkipAlignment) {
        chTrimmingCountFiles = chTrimmingCounts.map { sampleId, counts -> counts }.collect()
        chMappedBamFiles = chFilterQuality.map { sampleId, enrichment_mark, control, read_method, bam, version -> bam }.collect()
        chFinalBamFiles = chDACFilteredFiles.map { row -> row[4] }.collect()
        chAlignmentStatsReport = alignment_stats_report(chTrimmingCountFiles, chMappedBamFiles, chFinalBamFiles)
    } else {
        chAlignmentStatsReport = Channel.of("NO_DATA")
    }

    chCreateStatsSamtoolsfiltered = createStatsSamtoolsfiltered(chDACFilteredFiles)
    chIndexFiles = index_sam(chDACFilteredFiles)

    //SNP Fingerprint using SMaSH ************************************************
    chAllIndexFiles = chIndexFiles.collect()
    chAllBAMandBAIIndexFiles = chAllIndexFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.bam') || it.toString().endsWith('.bai') }}

    chSMaSHOutout = createSMaSHFingerPrint(chSNPSMaSH,chSNPS_ref,chAllBAMandBAIIndexFiles)
    chSNPSMaSHPlot = createSMaSHFingerPrintPlot(chSMaSHOutout,chSNPSMaSHPyPlot)

    // Collect all the files to generate the MultiQC report
    chSortBamAll = chSortBam.collect()
    chFilterPPAll = chFilterPP.collect()
    chLibComplexPreseqAll = chLibComplexPreseq.collect()
    chUniqueSamAll = chUniqueSam.collect()
    chFilterQualityAll = chFilterQuality.collect()
    chCreateStatsSamtoolsfilteredAll = chCreateStatsSamtoolsfiltered.collect()
    chDedupAll = chDedup.collect()
    chDACFilteredFilesAll = chDACFilteredFiles.collect()
    chSMaSHOutoutAll = chSMaSHOutout.collect()
    chSNPSMaSHPlotAll = chSNPSMaSHPlot.collect()
    chFilesReportInitializationAll = chFilesReportInitialization.collect()
    chFilesReportAlignmentAll = chFilesReportAlignment.collect()

    // Combine all the channels
    chAllChannels = chSortBamAll
        .combine(chSortBamAll)
        .combine(chFilterPPAll)
        .combine(chLibComplexPreseqAll)
        .combine(chUniqueSamAll)
        .combine(chCreateStatsSamtoolsfilteredAll)
        .combine(chFilterQualityAll)
        .combine(chCreateStatsSamtoolsfilteredAll)
        .combine(chDedupAll)
        .combine(chDACFilteredFilesAll)
        .combine(chAlignmentStatsReport)
        .combine(chSMaSHOutoutAll)
        .combine(chSNPSMaSHPlotAll)
        .combine(chFilesReportInitializationAll)
        .combine(chFilesReportAlignmentAll)
    
    // Filter only the files that will be used in the MultiQC report and remove duplicates
    chOnlyFiles = chAllChannels
    .flatten() // Make sure the files are in a single flow
    .collect() // Joins all files before processing them
    .map { files -> 
        def uniqueFiles = [:] as LinkedHashMap
        files.findAll { it instanceof Path } // Keeps only files (Path)
             .each { file -> uniqueFiles.putIfAbsent(file.getName(), file) } // Keeps only the first occurrence of the name
        return uniqueFiles.values()  // Returns only unique files
    } 
    .flatten()
    chFilesReportBamProcessing = chOnlyFiles.collect()

    // Create the MultiQC report and move the soft files only if this is the last process
    if (params.until == 'BAM_PROCESSING') {
        chBAMProcessReport = multiqc(chInitReport,chFilesReportBamProcessing,chMultiQCConfig)
        moveSoftFiles(chBAMProcessReport)
    } else {
        chBAMProcessReport = Channel.of("NO_DATA")
    }

    emit: bam_processed = chDACFilteredFiles
    emit: bam_processed_index = chIndexFiles
    emit: bam_and_bai_files = chAllBAMandBAIIndexFiles
    emit: report_SNP_SMaSH = chSNPSMaSHPlot
    emit: lib_complex = chLibComplexPreseq
    emit: files_report_bam_processing = chFilesReportBamProcessing
    emit: bam_process_report = chBAMProcessReport
}
