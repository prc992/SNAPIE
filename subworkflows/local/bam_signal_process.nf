nextflow.enable.dsl=2

include {call_peaks} from '../../modules/local/call_peaks'
include {peaks_annotations} from '../../modules/local/peaks_annotations.nf'
include {peaks_report} from '../../modules/local/peaks_report.nf'
include {bam_to_bedgraph} from '../../modules/local/bam_to_bedgraph'
include {bedgraph_to_bigwig} from '../../modules/local/bedgraph_to_bigwig'
include {igv_reports} from '../../modules/local/igv_reports'
include {igv_sample_reports} from '../../modules/local/igv_reports'
include {igv_consolidate_report} from '../../modules/local/igv_reports'
include {igv_session} from '../../modules/local/igv_reports'
include {enrichment} from '../../modules/local/enrichment'
include {enrichmentReport} from '../../modules/local/enrichmentReport'
//include {signalIntensityCalculation} from '../../modules/local/signal_intesity.nf'
include {chromatin_count_normalization_single} from '../../modules/local/chromatin_count_normalization.nf'
include {chromatin_count_normalization_batch} from '../../modules/local/chromatin_count_normalization.nf'
include {merge_enrichment_reports} from '../../modules/local/merge_enrichment_reports'
include {merge_signal_reports} from '../../modules/local/merge_signal_reports'
include {quality_report_lite} from '../../modules/local/quality_report_lite'
include {signal_report_lite} from '../../modules/local/signal_report_lite'
include {multiqc} from '../../modules/local/multiqc'
include {moveSoftFiles} from '../../modules/local/moveSoftFiles'

workflow BAM_SIGNAL_PROCESSING {

    take:
    chSampleInfo
    chGenome
    chGenomeIndex
    chChromSizes
    chSkipAlignment
    chBAMProcessedFiles
    chBAMProcessedIndexFiles
    chGenomesInfo
    chMultiQCHousekeepingHeader
    chMultiQCEnrichmentHeader
    chMultiQCPeaksHeader
    chIGVFilestoSessions
    chEnrichmentColors
    chEnrichmentFilesFolder
    chMultiQCConfig
    chEnrichmentScript
    chRGenomicAnnotation
    chPileUpBED
    chReportPeaks
    chReportEnrichment
    chMergeReportEnrichment
    chReportQualityLite
    chCTFragleFilesReport
    chFilesReportInitialization
    chFilesReportBamProcessing
    chFilesReportFragmentsProcess
    chBedFiles
    chDACFileRef
    chTSSPromoterPeaks_ref
    chRMEDIPSignalCalculation
    chRMARKSSignalCalculation
    chRegions_of_interest_MEDIP_signal
    chRegions_of_interest_MARKS_signal
    chHousekeeping_MEDIP_signal
    chHousekeeping_H3K4ME3_signal
    chHousekeeping_H3K27AC_signal
    chMultiQCSignalHeader
    chMergeReportSignal
    chFragsProcessReport
    

    main:

    chBedGraphFiles = bam_to_bedgraph(chBAMProcessedIndexFiles)
    chBigWig = bedgraph_to_bigwig(chBedGraphFiles,chChromSizes)

    //Pileups ****************************************************************
    chIGVReportsHtml = igv_sample_reports(chBedGraphFiles,chPileUpBED,chGenome,chGenomeIndex).collect()
    chIGVReportMerged = igv_consolidate_report(chIGVReportsHtml,chMultiQCHousekeepingHeader)

    chBigWigAllFiles = chBigWig.collect()
    chBigWigOnlyFiles = chBigWigAllFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.bw') }} // Filter the bw files
    chIGVSession = igv_session(chBigWigOnlyFiles,chIGVFilestoSessions,chEnrichmentColors,chGenomesInfo,chPileUpBED)

    // Match the samples with the controls
    def fake_control = file(params.dummy_control_file)
    chSamplesListCombine = chBAMProcessedFiles.combine(chBAMProcessedFiles)

    chSamplesListFilter = chSamplesListCombine.filter { row -> row[2] == row[6] }.map { row -> [row[0],row[1],row[3], row[4], row[10]] }
    chSamplesListNoControl = chBAMProcessedFiles.filter { row -> !row[2] }.map { row -> [row[0],row[1],row[3], row[4], fake_control] }
    chSamplesListMix = chSamplesListFilter.mix(chSamplesListNoControl)

    chPeakFiles = call_peaks(chSamplesListMix)

    //CHROMATIN COUNT NORMALIZATION *******************************
    chReferenceSitesCCN = params.chromatin_count_reference ? \
    Channel.fromPath(params.chromatin_count_reference, checkIfExists: true) : Channel.fromPath(params.dummy_control_file)
    chTargetSitesCCN = params.chromatin_count_target_sites ? \
    Channel.fromPath(params.chromatin_count_target_sites, checkIfExists: true) : chTSSPromoterPeaks_ref

    chromatin_count_mode = params.chromatin_count_mode.toLowerCase()
    if (chromatin_count_mode == "single") {
        chChromatinCountNormalization = chromatin_count_normalization_single(chBedFiles,chReferenceSitesCCN,chTargetSitesCCN)
    } else if (chromatin_count_mode == "batch") {

        chBedsAllFiles = chBedFiles.collect()
        chOnlyBedsFiles = chBedsAllFiles.map { collectedFiles ->
        collectedFiles.findAll { it.toString().endsWith('.bed') }} // Filter the bed files
        chChromatinCountNormalization = chromatin_count_normalization_batch(chOnlyBedsFiles,chReferenceSitesCCN,chTargetSitesCCN)
    }
    //***************************************************************

    chPeakAllFiles = chPeakFiles.collect()
    chNarrowPeakFiles = chPeakAllFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.narrowPeak') }} // Filter the narrowPeak files

    chPeaksFilesReport = peaks_report(chNarrowPeakFiles,chMultiQCPeaksHeader,chReportPeaks)

     if (params.report_peak_genomic_annotation == true){
        chPeaksAnnotationReport = peaks_annotations(chNarrowPeakFiles,chRGenomicAnnotation)
     } else {
        chPeaksAnnotationReport = Channel.of("NO_DATA")
     }

    //chPeaksAnnotationReport = peaks_annotations(chNarrowPeakFiles,chRGenomicAnnotation)
    
    //ENRICHMENT *********************************************************************
    chEnrichmentFilesCSV = enrichment(chBAMProcessedFiles,chEnrichmentFilesFolder,chEnrichmentScript).collect()
    chEnrichmentFilesReport = enrichmentReport(chSampleInfo,chEnrichmentFilesCSV,chReportEnrichment).collect()
    chMergedEnrichmentReport = merge_enrichment_reports(chEnrichmentFilesReport,chMultiQCEnrichmentHeader,chMergeReportEnrichment).collect()
    
    //Signal Intensity Calculation ********************************
    /*chSignalFilesReport = signalIntensityCalculation(chBedFiles,chDACFileRef,\
                                chRMEDIPSignalCalculation,chRMARKSSignalCalculation,
                                chRegions_of_interest_MEDIP_signal,chRegions_of_interest_MARKS_signal,\
                                chHousekeeping_MEDIP_signal,chHousekeeping_H3K4ME3_signal,chHousekeeping_H3K27AC_signal).collect()
    chMergedSignalReport = merge_signal_reports(chSignalFilesReport,chMultiQCSignalHeader,chMergeReportSignal,chSampleInfo).collect()*/
    //********************************
    //********************************

    quality_report_lite(chReportQualityLite,chEnrichmentFilesReport,chPeaksFilesReport,chFragsProcessReport,chCTFragleFilesReport)
    //chReportQualityLite.view()
    //signal_report_lite(chMergedSignalReport)
    
    // Collect all the files to generate the MultiQC report
    chBedGraphFilesAll = chBedGraphFiles.collect()
    chBigWigAllFiles = chBigWig.collect()
    chIGVReportsHtmlAll = chIGVReportsHtml.collect()
    chIGVReportMergedAll = chIGVReportMerged.collect()
    chBigWigAllFilesAll = chBigWigAllFiles.collect()
    chIGVSessionAll = chIGVSession.collect()
    chPeakFilesAll = chPeakFiles.collect()
    chPeaksFilesReportAll = chPeaksFilesReport.collect()
    chEnrichmentFilesCSVAll = chEnrichmentFilesCSV.collect()
    chEnrichmentFilesReportAll = chEnrichmentFilesReport.collect()
    chMergedEnrichmentReportAll = chMergedEnrichmentReport.collect()
    //chMergedSignalReportAll = chMergedSignalReport.collect()
    

    // Combine all the channels
    chAllChannelsProcessing = chBedGraphFilesAll
        .combine(chBigWigAllFiles)
        .combine(chIGVReportsHtmlAll)
        .combine(chIGVReportMergedAll)
        .combine(chBigWigAllFilesAll)
        .combine(chIGVSessionAll)
        .combine(chPeakFilesAll)
        .combine(chPeaksFilesReportAll)
        .combine(chPeaksAnnotationReport)
        .combine(chEnrichmentFilesCSVAll)
        .combine(chEnrichmentFilesReportAll)
        .combine(chMergedEnrichmentReportAll)
        //.combine(chMergedSignalReportAll)
        .combine(chFilesReportBamProcessing)
        .combine(chFilesReportInitialization)
        .combine(chFilesReportFragmentsProcess)
    
    chOnlyFilesProcessing = chAllChannelsProcessing
    .flatten() // Garante que os arquivos estejam em um único fluxo
    .collect() // Junta todos os arquivos antes de processá-los
    .map { files -> 
        def uniqueFiles = [:] as LinkedHashMap
        files.findAll { it instanceof Path } // 🔹 Mantém apenas arquivos (Path)
             .each { file -> uniqueFiles.putIfAbsent(file.getName(), file) } // Mantém apenas a primeira ocorrência do nome
        return uniqueFiles.values() // Retorna apenas os arquivos únicos
    }
    .flatten() // Garante que cada arquivo seja emitido separadamente no canal

    chFilesReportSignalProcess = chOnlyFilesProcessing.collect()
    chBAMSignalReport = multiqc(chFragsProcessReport,chFilesReportSignalProcess,chMultiQCConfig)
    moveSoftFiles(chBAMSignalReport)
}
