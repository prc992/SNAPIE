 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

//local modules

//subworkflows
include { PREPROCESSING } from './subworkflows/local/preprocessing.nf'
include { DOWNLOAD_REFERENCES } from './subworkflows/local/download_references'
include { ALIGNMENT } from './subworkflows/local/alignment'
include { BAM_PROCESSING } from './subworkflows/local/bam_processing'
include { FRAGMENTS_PROCESSING } from './subworkflows/local/fragments_processing'
include { BAM_SIGNAL_PROCESSING } from './subworkflows/local/bam_signal_process'


workflow  {
    // Static information about the pipeline
    def githubPath = "https://github.com/prc992/SNAPIE"
    def releaseVersion = "v1.5.18"

    // ASCII art for SNAPIE
    def asciiArt = """
      ███████╗ ███╗   ██╗ █████╗ ██████╗  ██╗ ██████╗
      ██╔════╝ ████╗  ██║██╔══██╗██╔══██╗ ██║ ██╗
      ███████╗ ██╔██╗ ██║███████║██████╔╝ ██║ ██████╗
      ╚════██║ ██║╚██╗██║██╔══██║██╔═══╝  ██║ ██╗
      ███████╗ ██║ ╚████║██║  ██║██║      ██║ ██████╗
      ╚══════╝ ╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝      ╚═╝ ╚═════╝
    """

    // Print the introductory message
    println asciiArt
    println "SNAPIE pipeline running, created by BacaLab. https://bacalab.dana-farber.org/"
    println "SNAPIE: Streamlined Nextflow Analysis Pipeline for Immunoprecipitation-Based Epigenomic Profiling of Circulating Chromatin"
    println "GitHub repository: ${githubPath}"
    println "Release version: ${releaseVersion}"
    println "Running pipeline with profile: ${workflow.profile}"

    //Auxiliar code
    chEnrichmentScript= Channel.fromPath("$params.pathEnrichmentScript")
    chRPileups= Channel.fromPath("$params.pathRPileups")
    chSNPSMaSH = Channel.fromPath("$params.pathSNPSMaSH")
    chSNPSMaSHPyPlot = Channel.fromPath("$params.pathSNPSMaSHPlot")
    chReportFragHist = Channel.fromPath("$params.pathReportFragHist")
    chReportFrags = Channel.fromPath("$params.pathReportFrags")
    chReportPeaks = Channel.fromPath("$params.pathReportPeaks")
    chReportCT = Channel.fromPath("$params.pathReportCT")
    chReportEnrichment = Channel.fromPath("$params.pathReportEnrichment")
    chReportQualityLite = Channel.fromPath("$params.pathReportQualityLite")
    chMergeReportEnrichment = Channel.fromPath("$params.pathMergeReportEnrichment")
    chMergeReportSignal = Channel.fromPath("$params.pathMergeReportSignal")
    chIGVFilestoSessions = Channel.fromPath("$params.pathIGVFilestoSessions")
    chRGenomicAnnotation = Channel.fromPath("$params.pathRGenomicAnnotation")

    chRMEDIPSignalCalculation= Channel.fromPath("$params.pathRMEDIPSignalCalculation")
    chRMARKSSignalCalculation= Channel.fromPath("$params.pathRMARKSSignalCalculation")

    //Assets
    chPileUpBED = Channel.fromPath("$params.genes_pileup_report")

    chRegions_of_interest_MEDIP_signal = Channel.fromPath("$params.regions_of_interest_MEDIP_signal")
    chRegions_of_interest_MARKS_signal = Channel.fromPath("$params.regions_of_interest_MARKS_signal")
    chHousekeeping_MEDIP_signal = Channel.fromPath("$params.housekeeping_MEDIP_signal")
    chHousekeeping_H3K4ME3_signal = Channel.fromPath("$params.housekeeping_H3K4ME3_signal")
    chHousekeeping_H3K27AC_signal = Channel.fromPath("$params.housekeeping_H3K27AC_signal")

    chMultiQCHousekeepingHeader = Channel.fromPath("$params.multiqc_housekeeping_header")
    chMultiQCFragsHeader = Channel.fromPath("$params.multiqc_tot_frag_header")
    chMultiQCPeaksHeader = Channel.fromPath("$params.multiqc_tot_peaks_header")
    chMultiQCCTHeader = Channel.fromPath("$params.multiqc_tot_ct_header")
    chMultiQCEnrichmentHeader = Channel.fromPath("$params.multiqc_enrichment_header")
    chMultiQCSignalHeader = Channel.fromPath("$params.multiqc_signal_header")
    chEnrichmentColors = Channel.fromPath("$params.enrichment_states_colors")
    chEnrichmentFilesFolder = Channel.fromPath("$params.enrichment_states_ref")

    if (params.report_peak_genomic_annotation == true) {
        chMultiQCConfig = Channel.fromPath("$params.multiqc_config_with_peak_annotation")
    } else {
        chMultiQCConfig = Channel.fromPath("$params.multiqc_config_no_peak_annotation")
    }

    if (params.samplesheetBams || params.sample_dir_bam) {
        skip_alignment = true
    } else {
        skip_alignment = false
    }

    def steps = ['PREPROCESSING', 'DOWNLOAD_REFERENCES','ALIGNMENT', 'BAM_PROCESSING', 'FRAGMENTS_PROCESSING','BAM_SIGNAL_PROCESSING']
    def run_steps = steps.takeWhile { it != params.until } + params.until
    
    if ('PREPROCESSING' in run_steps) {
        PREPROCESSING(chMultiQCConfig,skip_alignment)
        chGenomesInfo = PREPROCESSING.out.genomes_info
        refDir = PREPROCESSING.out.ref_dir
        chSampleInfo = PREPROCESSING.out.sample_info
        chFastaQC = PREPROCESSING.out.fastqc_files
        chFilesReportInitialization = PREPROCESSING.out.files_report_initialization
        chInitReport = PREPROCESSING.out.init_report
        }

    if ('DOWNLOAD_REFERENCES' in run_steps) {
        DOWNLOAD_REFERENCES(chGenomesInfo,refDir)

        chGenome = DOWNLOAD_REFERENCES.out.genome
        chGenomeIndex = DOWNLOAD_REFERENCES.out.genome_index
        chChromSizes = DOWNLOAD_REFERENCES.out.chrom_sizes
        chDACFileRef = DOWNLOAD_REFERENCES.out.dac_file_ref
        chSNPS_ref = DOWNLOAD_REFERENCES.out.snp_ref
        chTSSPromoterPeaks_ref = DOWNLOAD_REFERENCES.out.tss_promoter_peaks_ref
        }

    
    if ('ALIGNMENT' in run_steps) {
        
        if (skip_alignment) {
            chAlign = chSampleInfo.map { sampleId, enrichment_mark, bam,control,read_method -> 
                        tuple(sampleId, enrichment_mark, control, read_method,file(bam),null) 
                        }
            chFilesReportAlignment = Channel.of("NO_DATA")
            chAlignmentReport = Channel.of("NO_DATA")
            chTrimmingCounts = Channel.of("NO_DATA")

        } else {
            ALIGNMENT (chSampleInfo,chGenome,chGenomeIndex,\
                    chFilesReportInitialization,chInitReport,\
                    chMultiQCConfig)
            chAlign = ALIGNMENT.out.align
            chTrimmingCounts = ALIGNMENT.out.trimming_counts
            chFilesReportAlignment = ALIGNMENT.out.files_report_alignment
            chAlignmentReport = ALIGNMENT.out.aligment_report
        }
        

        
        }

    if ('BAM_PROCESSING' in run_steps) {
        BAM_PROCESSING (chSampleInfo,chGenome, chGenomeIndex,chChromSizes,chDACFileRef,chSNPS_ref,
                        chAlign,
                        chTrimmingCounts,skip_alignment,
                        chSNPSMaSH,chSNPSMaSHPyPlot,chMultiQCConfig,
                        chFilesReportInitialization,chInitReport,
                        chFilesReportAlignment,chAlignmentReport)

        chBAMProcessedFiles = BAM_PROCESSING.out.bam_processed
        chBAMProcessedIndexFiles = BAM_PROCESSING.out.bam_processed_index
        chBAMBAIProcessedFiles = BAM_PROCESSING.out.bam_and_bai_files
        chSNPSMaSHPlot = BAM_PROCESSING.out.report_SNP_SMaSH
        chLibComplexPreseq = BAM_PROCESSING.out.lib_complex
        chFilesReportBamProcessing = BAM_PROCESSING.out.files_report_bam_processing
        chBAMProcessReport = BAM_PROCESSING.out.bam_process_report
        }
    
    if ('FRAGMENTS_PROCESSING' in run_steps) {

        FRAGMENTS_PROCESSING(chSampleInfo,chGenome,chGenomeIndex,
                            chBAMProcessedFiles,chBAMProcessedIndexFiles,chBAMBAIProcessedFiles,
                            chMultiQCFragsHeader,chMultiQCCTHeader,chReportFrags,
                            chMultiQCConfig,
                            chFilesReportInitialization,chFilesReportBamProcessing,
                            chBAMProcessReport,chReportCT)

        chFragmentsSizeFiles = FRAGMENTS_PROCESSING.out.frag_size_files
        chFragReport = FRAGMENTS_PROCESSING.out.frag_report
        chFragsProcessReport = FRAGMENTS_PROCESSING.out.frag_process_report
        chFilesReportFragmentsProcess = FRAGMENTS_PROCESSING.out.files_report_fragments_processing
        chCTFragleFilesReport = FRAGMENTS_PROCESSING.out.files_fragle_report
        chBedFiles = FRAGMENTS_PROCESSING.out.files_sample_bed
        }

    if ('BAM_SIGNAL_PROCESSING' in run_steps) {
        BAM_SIGNAL_PROCESSING (chSampleInfo, chGenome, chGenomeIndex,chChromSizes,skip_alignment,
                            chBAMProcessedFiles,chBAMProcessedIndexFiles,
                            chGenomesInfo,chMultiQCHousekeepingHeader,chMultiQCEnrichmentHeader,chMultiQCPeaksHeader,chIGVFilestoSessions,chEnrichmentColors,chEnrichmentFilesFolder,
                            chMultiQCConfig,chEnrichmentScript,chRGenomicAnnotation,chPileUpBED,chReportPeaks,chReportEnrichment,chMergeReportEnrichment,chReportQualityLite,chCTFragleFilesReport,
                            chFilesReportInitialization,chFilesReportBamProcessing,chFilesReportFragmentsProcess,
                            chBedFiles,chDACFileRef,chTSSPromoterPeaks_ref,
                            chRMEDIPSignalCalculation,chRMARKSSignalCalculation,
                            chRegions_of_interest_MEDIP_signal,chRegions_of_interest_MARKS_signal,
                            chHousekeeping_MEDIP_signal,chHousekeeping_H3K4ME3_signal,chHousekeeping_H3K27AC_signal,
                            chMultiQCSignalHeader,chMergeReportSignal,
                            chFragReport)



        }
    
}
