process update_pangolin {

  tag { should_update.toString() }

  executor 'local'

  input:
  val(should_update)

  output:
  val(did_update)

  script:
  did_update = should_update
  should_update_string = should_update ? "true" : "false"
  """
  should_update=${should_update_string}
  if [ "\$should_update" = true ]
  then
    pangolin --update
  fi
  """
}

process download_ncov_tools {

  tag { version }

  executor 'local'
  
  input:
  val(version)
  
  output:
  path("ncov-tools", type: 'dir')

  script:
  """
  wget https://github.com/BCCDC-PHL/ncov-tools/archive/v${version}.tar.gz
  tar -xzf v${version}.tar.gz
  mv ncov-tools-${version} ncov-tools
  """
}

process download_artic_ncov2019 {

  tag { version }

  executor 'local'
  
  input:
  tuple val(version), val(primer_scheme)
  
  output:
  path("resources", type: 'dir')

  script:
  """
  wget ${params.artic-ncov2019}/v${version}.tar.gz
  tar -xzf v${version}.tar.gz
  mkdir resources
  cp artic-ncov2019-${version}/primer_schemes/nCoV-2019/${primer_scheme}/nCoV-2019.reference.fasta resources
  cp artic-ncov2019-${version}/primer_schemes/nCoV-2019/${primer_scheme}/nCoV-2019.bed resources
  """
}

process index_reference_genome {

  tag { "MN908947.3" }

  executor 'local'
  
  input:
  path(resources)
  
  output:
  path("resources", type: 'dir')

  script:
  """
  samtools faidx resources/nCoV-2019.reference.fasta
  """
}

process get_library_plate_ids {

  tag { params.run_name }

  executor 'local'
  
  input:
  path(artic_analysis_dir)
  
  output:
  stdout

  script:
  """
  tail -n+2 ${artic_analysis_dir}/*.qc.csv | \
  grep -v 'POS*' | grep -v 'NEG*' |
  cut -f 1 -d ',' | cut -f 2 -d '-' | sort | uniq
  """
}

process prepare_data_root {

  tag {params.run_name}

  executor 'local'
  
  input:
  tuple path(ncov2019_artic_nf_analysis_dir), path(primer_scheme_dir), path(metadata), val(library_plate_id)
  
  output:
  tuple val(library_plate_id), path("ncov-tools-input", type: 'dir')

  script:
  def metadata = metadata.name != 'NO_FILE' ? "cp ${metadata} ncov-tools-input" : ''
  def filename_glob = "*"
  def link_pre_downsampled_bams = params.pre_downsampled ? "ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_trimPrimerSequences/${filename_glob} ." : ''
  def link_freebayes_consensus = params.ivar_consensus ? '' : "ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_callConsensusFreebayes/${filename_glob}.fa ."
  def link_ivar_consensus = params.ivar_consensus ? "ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_makeConsensus/${filename_glob}.fa ." : ''
  def link_ivar_variants = params.ivar_variants ? "ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_callVariants/${filename_glob}.tsv ." : ''
  """
  mkdir ncov-tools-input
  ${metadata}
  pushd ncov-tools-input
  ${link_ivar_consensus}
  ${link_freebayes_consensus}
  ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_readMapping/${filename_glob} .
  ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_downsampleAmplicons/${filename_glob} .
  ${link_pre_downsampled_bams}
  ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_callConsensusFreebayes/${filename_glob}.vcf .
  ${link_ivar_variants}
  ln -sfn ../${primer_scheme_dir}/nCoV-2019.reference.fasta .
  ln -sfn ../${primer_scheme_dir}/nCoV-2019.bed .
  popd
  """
}

process create_sample_id_list {

  tag {params.run_name}

  executor 'local'

  input:
  tuple val(library_plate_id), path(data_root)

  output:
  tuple val(library_plate_id), path("sample_id_list.tsv")

  script:
  """
  find ${data_root}/ -name '*.variants.tsv' | xargs -n 1 basename | sed 's/\\.variants\\.tsv//' | sort > sample_id_list.tsv
  """
}

process find_negative_control {

  tag {params.run_name}

  executor 'local'
  
  input:
  tuple val(library_plate_id), path(data_root)
  
  output:
  tuple val(library_plate_id), path("neg_control_sample_id.txt")

  script:
  def filename_glob = "*"
  """
  find ${data_root}/ -name NEG${filename_glob}.consensus.fa -printf "%f" | cut -d '.' -f 1 > neg_control_sample_id.txt
  """
}

process create_config_yaml {

  tag {params.run_name}

  executor 'local'

  input:
  tuple path(negative_control_sample_id), val(metadata)
  
  output:
  tuple path("config.yaml")

  script:
  def metadata = metadata.name != 'NO_FILE' ? "metadata: \\\"{data_root}/metadata.tsv\\\"" : ''
  """
  echo "data_root: ncov-tools-input" >> config.yaml
  echo "run_name: ${params.run_name}" >> config.yaml
  if [[ \$( wc -l < ${negative_control_sample_id} ) -ge 1 ]]; then echo "negative_control_samples: [ \\"\$( cat ${negative_control_sample_id} )\\" ]" >> config.yaml; fi
  echo "${metadata}" >> config.yaml
  echo "reference_genome: \\"resources/nCoV-2019.reference.fasta\\"" >> config.yaml
  echo "primer_bed: \\"resources/nCoV-2019.bed\\"" >> config.yaml
  echo "bam_pattern: \\"{data_root}/{sample}.sorted.bam\\"" >> config.yaml
  echo "consensus_pattern: \\"{data_root}/{sample}.consensus.fa\\"" >> config.yaml
  echo "variants_pattern: \\"{data_root}/{sample}.variants.tsv\\"" >> config.yaml
  echo "platform: ${params.platform}" >> config.yaml
  echo "bed_type: unique_amplicons" >> config.yaml
  echo "offset: 0" >> config.yaml
  echo "completeness_threshold: ${params.completeness_threshold}" >> config.yaml
  echo "assign_lineages: true" >> config.yaml
  """
}


process ncov_tools {

  tag { params.run_name }

  publishDir "${params.outdir}", mode: 'copy', pattern: "config.yaml",                  
  publishDir "${params.outdir}", mode: 'copy', pattern: "bed",                          
  publishDir "${params.outdir}", mode: 'copy', pattern: "lineages/${params.run_name}*", 
  publishDir "${params.outdir}", mode: 'copy', pattern: "plots",                        
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_analysis",                  
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_reports/*.tsv",             
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_sequencing",                
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_annotation",                

  input:
  tuple val(library_plate_id), path(config_yaml), path(data_root), path(resources), path(ncov_tools), val(pangolin_updated)
  
  output:
  path("config.yaml")
  path("bed")
  path("lineages/${params.run_name}*_lineage_report.csv"), emit: lineage_report
  path("lineages/${params.run_name}*_pangolin_version.txt"), emit: pangolin_version
  path("plots")
  path("qc_analysis")
  path("qc_reports/*.tsv"), emit: qc_reports
  path("qc_sequencing")
  path("qc_annotation")

  script:
  """
  snakemake -s ./ncov-tools/workflow/Snakefile --cores ${task.cpus} all
  snakemake -s ./ncov-tools/workflow/Snakefile --cores ${task.cpus} all_qc_annotation
  """
}

process combine_all_qc_summaries_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}_summary_qc.tsv"

  input:
  path(summaries)

  output:
  path("${params.run_name}_summary_qc.tsv")

  script:
  """
  head -qn 1 *_summary_qc.tsv | uniq > header.tsv
  tail -qn+2 *_summary_qc.tsv | sort -k1,1 -k2,2 > data.tsv
  cat header.tsv data.tsv > "${params.run_name}_summary_qc_incorrect_run_name.tsv"
  tail -qn+2 "${params.run_name}_summary_qc_incorrect_run_name.tsv" | awk -F '\t' 'BEGIN {OFS=FS}; {split(\$2, run_id_split, "_"); \$2 = run_id_split[1]"_"run_id_split[2]"_"run_id_split[3]"_"run_id_split[4]; print}' > data.tsv
  cat header.tsv data.tsv > "${params.run_name}_summary_qc.tsv"
  """
}

process combine_all_mixture_reports_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}_mixture_report.tsv"

  input:
  path(mixture_reports)

  output:
  path("${params.run_name}_mixture_report.tsv")

  script:
  """
  head -qn 1 *_mixture_report.tsv | uniq > header.tsv
  tail -qn+2 *_mixture_report.tsv | sort -k1,1 -k2,2 > data.tsv
  cat header.tsv data.tsv > "${params.run_name}_mixture_report.tsv"
  """
}

process combine_ambiguous_position_reports_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}_ambiguous_position_report.tsv"

  input:
  path(ambiguous_position_reports)

  output:
  path("${params.run_name}_ambiguous_position_report.tsv")

  script:
  """
  head -qn 1 *_ambiguous_position_report.tsv | uniq > header.tsv
  tail -qn+2 *_ambiguous_position_report.tsv | sort -k1,1 -k2,2 > data.tsv
  cat header.tsv data.tsv > "${params.run_name}_ambiguous_position_report.tsv"
  """
}

process combine_all_lineage_reports_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/lineages", mode: 'copy', pattern: "${params.run_name}_lineage_report.csv"

  input:
  path(lineage_reports)

  output:
  path("${params.run_name}_lineage_report.csv")

  script:
  """
  head -qn 1 *_lineage_report.csv | uniq > header.csv
  tail -qn+2 *_lineage_report.csv | sort -k1,1 -k2,2 > data.csv
  cat header.csv data.csv > "${params.run_name}_lineage_report.csv"
  """
}

process get_pangolin_version_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/lineages", mode: 'copy', pattern: "${params.run_name}_pangolin_version.txt"

  input:
  path(lineage_reports)

  output:
  path("${params.run_name}_pangolin_version.txt")

  script:
  """
  cat ${params.run_name}*_pangolin_version.txt > ${params.run_name}_pangolin_version.txt
  """
}
