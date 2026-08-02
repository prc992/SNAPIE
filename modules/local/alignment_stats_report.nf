process trimming_stats {
  label 'low_cpu_low_mem'
  container = params.containers.samtools

  tag "Sample - $sampleId"

  input:
  tuple val(sampleId), path(trimmingReports)

  output:
  tuple val(sampleId), path("${sampleId}.trimming_counts.tsv")

  script:
  """
  if [[ "${params.trim_method}" == "FASTP" ]]; then
      total_reads=\$(awk '
        /"before_filtering"[[:space:]]*:/ { section=1; next }
        section && /"total_reads"[[:space:]]*:/ {
            gsub(/[^0-9]/, "", \$0); print; exit
        }
      ' $trimmingReports)

      filtered_reads=\$(awk '
        /"after_filtering"[[:space:]]*:/ { section=1; next }
        section && /"total_reads"[[:space:]]*:/ {
            gsub(/[^0-9]/, "", \$0); print; exit
        }
      ' $trimmingReports)

      if [[ -z "\$total_reads" || -z "\$filtered_reads" ]]; then
          echo "Could not read before/after filtering totals for $sampleId" >&2
          exit 1
      fi

      printf '%s\\t%s\\t%s\\n' '$sampleId' "\$total_reads" "\$filtered_reads" > ${sampleId}.trimming_counts.tsv
  else
      awk -F: '
        /^Total reads processed:/ {
            gsub(/[ ,]/, "", \$2); total += \$2
        }
        /^Reads written/ && /passing filters/ {
            gsub(/^[[:space:]]+/, "", \$2)
            split(\$2, value, /[[:space:]]+/)
            gsub(/,/, "", value[1]); filtered += value[1]
        }
        END {
            if (total == 0 && filtered == 0) exit 1
            printf "%s\\t%.0f\\t%.0f\\n", "$sampleId", total, filtered
        }
      ' $trimmingReports > ${sampleId}.trimming_counts.tsv
  fi
  """
}

process alignment_stats_report {
  label 'low_cpu_low_mem'
  container = params.containers.samtools

  tag "All Samples"

  input:
  path trimmingCountFiles
  path mappedBamFiles
  path finalBamFiles

  output:
  path "alignment_statistics_mqc.tsv"

  script:
  """
  printf 'Sample\\tTotal reads\\tFiltered reads\\tMapped Reads\\tDeduplicated reads\\n' > alignment_statistics_mqc.tsv

  for counts in *.trimming_counts.tsv; do
      sample=\${counts%.trimming_counts.tsv}
      mapped_bam="\${sample}.filtered.unique.sorted.bam"
      if [[ "${params.exclude_dac_regions}" == "true" ]]; then
          final_bam="\${sample}.dac_filtered.dedup.unique.sorted.bam"
      else
          final_bam="\${sample}.dedup.unique.sorted.bam"
      fi

      if [[ ! -f "\$mapped_bam" || ! -f "\$final_bam" ]]; then
          echo "Missing mapped or final processed BAM for sample \$sample" >&2
          exit 1
      fi

      total_reads=\$(cut -f2 "\$counts")
      filtered_reads=\$(cut -f3 "\$counts")

      mapped_reads=\$(samtools view -c "\$mapped_bam")
      deduplicated_reads=\$(samtools view -c "\$final_bam")

      printf '%s\\t%s\\t%s\\t%s\\t%s\\n' \
          "\$sample" "\$total_reads" "\$filtered_reads" "\$mapped_reads" "\$deduplicated_reads" \
          >> alignment_statistics_mqc.tsv
  done
  """
}
