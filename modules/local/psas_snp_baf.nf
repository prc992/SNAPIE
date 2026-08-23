process psas_snp_baf {
  label 'low_cpu_low_mem'

  container = params.containers.psas

  tag "Sample - ${sampleId}"

  input:
  tuple val(sampleId), path(vcf), path(vcfIndex)
  path snpBafScript

  output:
  tuple val(sampleId), path("${sampleId}.psas.baf.tsv")

  script:
  """
  bash ${snpBafScript} ${vcf} ${sampleId}.psas.baf.tsv
  """
}
