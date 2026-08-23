process psas_genotype {
  label 'process_medium'

  container = params.containers.psas

  tag "Sample - ${sampleId}"

  input:
  tuple val(sampleId), val(enrichment_mark), val(control), val(read_method), path(sampleBam), val(_), path(genomeFile), path(genomeFai), path(genomeDict)
  path genotypeScript

  output:
  tuple val(sampleId), path("${sampleId}.psas.vcf.gz"), path("${sampleId}.psas.vcf.gz.tbi")

  script:
  """
  bash ${genotypeScript} ${sampleBam} ${genomeFile} ${genomeFai} ${genomeDict} ${sampleId}.psas.vcf.gz ${task.cpus}
  """
}
