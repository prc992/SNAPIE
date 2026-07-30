process merge_enrichment_reports {
    label 'low_cpu_low_mem'
    container = params.containers.python
    tag "All samples"

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode: 'copy'

    input:
    path (chEnrichmentFilesReport)
    each path (chMultiQCEnrichmentHeader)
    each path (chMergeReportEnrichment)
    
    output:
    path ("*.csv")

    script:
    """
    python $chMergeReportEnrichment
    
    # If no file *_mqc.csv was created because there is no enrichment mark availabe
    # create an empty one to satisfy the output requirement but that is not going to be used by multiqc

    if ! ls *_mqc.csv 1> /dev/null 2>&1; then
        touch no_enrichment_data.csv
    fi
    """
}
