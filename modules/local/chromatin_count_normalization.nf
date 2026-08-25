process chromatin_count_normalization_batch {
  label 'med_cpu_med_mem'
  container params.containers.chromatin_count_normalization

  tag "All Samples"
  publishDir "${workflow.projectDir}/${params.outputFolder}/chromatin_count_normalization/batch", mode: 'copy'

  input:
  tuple path(bedFiles), path(referenceSitesFile), val(useReferenceSites)
  each path (targetSitesFile)

  output:
  path "output", type: 'dir'

  script:
  def ref_arg = useReferenceSites ? "--reference-sites ${referenceSitesFile}" : ""

  """
  # Guard: if target sites file is missing/empty, create empty outputs and exit cleanly
  if [[ ! -f "${targetSitesFile}" || ! -s "${targetSitesFile}" ]]; then
    echo "[WARNING] targetSitesFile missing or empty: ${targetSitesFile}"
    echo "[WARNING] Skipping chromatin_count_norm_v2.R and generating empty outputs."

    mkdir -p output
    : > output/EMPTY_TARGET_SITES
    echo "target_sites is missing/empty: ${targetSitesFile}" > output/EMPTY_REASON.txt
    exit 0
  fi

  echo "Creating sample_name list with header..."
  > sample_name
  echo "sample_name" >> sample_name
  for f in ${bedFiles}; do
    bn=\$(basename "\$f")
    name="\${bn%.bed}"
    printf "%s\\n" "\$name" >> sample_name
  done

  Rscript /workspace/chromatin_count_norm_v2.R \
    --samplesheet sample_name \
    --target-sites ${targetSitesFile} \
    --frags-dir . \
    ${ref_arg} \
    --verbose
  """
}

process chromatin_count_normalization_single {
  label 'med_cpu_med_mem'
  container params.containers.chromatin_count_normalization

  tag "Sample - ${sampleId}"
  publishDir "${workflow.projectDir}/${params.outputFolder}/chromatin_count_normalization/${sampleId}", mode: 'copy'

  input:
  tuple val(sampleId), val(_), val(_), val(_), path(bedFile), val(_), path(referenceSitesFile), val(useReferenceSites)
  each path (targetSitesFile)

  output:
  path "output", type: 'dir'

  script:
  def ref_arg = useReferenceSites ? "--reference-sites ${referenceSitesFile}" : ""

  """
  # Guard: if target sites file is missing/empty, create empty outputs and exit cleanly
  if [[ ! -f "${targetSitesFile}" || ! -s "${targetSitesFile}" ]]; then
    echo "[WARNING] targetSitesFile missing or empty: ${targetSitesFile}"
    echo "[WARNING] Skipping chromatin_count_norm_v2.R and generating empty outputs."

    mkdir -p output
    : > output/EMPTY_TARGET_SITES
    echo "target_sites is missing/empty: ${targetSitesFile}" > output/EMPTY_REASON.txt
    exit 0
  fi

  Rscript /workspace/chromatin_count_norm_v2.R \
    --sample-name ${sampleId} \
    --fragment-file ${bedFile} \
    --target-sites ${targetSitesFile} \
    ${ref_arg} \
    --verbose
  """
}
