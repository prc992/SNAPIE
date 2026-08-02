nextflow.enable.dsl=2

// Import the required processes from the modules
include {trim} from '../../modules/local/trim'
include {trim_fastp} from '../../modules/local/trim'
include {align} from '../../modules/local/align'
include {multiqc} from '../../modules/local/multiqc'
include {moveSoftFiles} from '../../modules/local/moveSoftFiles'
include {trimming_stats} from '../../modules/local/alignment_stats_report'

workflow ALIGNMENT {


    take:
    chSampleInfo
    chGenome
    chGenomeIndex
    chFilesReportInitialization
    chInitReport
    chMultiQCConfig

    main:
    if (params.trim_method == 'FASTP') {
        chTrim = trim_fastp(chSampleInfo)
    } else {
        chTrim = trim(chSampleInfo)
    }

    chTrimmingReports = chTrim.map { sampleId, enrichment_mark, control, read_method, trimmedFiles, trimmingReport, trimmingHtml ->
        tuple(sampleId, trimmingReport)
    }
    chTrimmingCounts = trimming_stats(chTrimmingReports)

    chAlign = align(chTrim,chGenome,chGenomeIndex)

    // Collect all the files to generate the MultiQC report
    chTrimAll = chTrim.collect()
    chAlignAll = chAlign.collect()

    // Combine all the channels
    chAllChannels = chTrimAll
        .combine(chAlignAll)
        .combine(chFilesReportInitialization)
    
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
    chFilesReportAlignment = chOnlyFiles.collect()

    // Create the MultiQC report and move the soft files only if this is the last process
    if (params.until == 'ALIGNMENT') {
        chAlignmentReport = multiqc(chAlignAll,chFilesReportAlignment,chMultiQCConfig)
        moveSoftFiles(chAlignmentReport)
    } else {
        chAlignmentReport = Channel.of("NO_DATA")
    }

    emit: align = chAlign
    emit: trimming_counts = chTrimmingCounts
    emit: files_report_alignment = chFilesReportAlignment
    emit: aligment_report = chAlignmentReport

}
