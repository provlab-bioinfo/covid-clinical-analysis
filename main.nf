#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { update_pangolin } from './modules/ncov-tools.nf'
include { download_ncov_tools } from './modules/ncov-tools.nf'
include { download_artic_ncov2019 } from './modules/ncov-tools.nf'
include { index_reference_genome } from './modules/ncov-tools.nf'
include { get_library_plate_ids } from './modules/ncov-tools.nf'
include { prepare_data_root } from './modules/ncov-tools.nf'
include { create_sample_id_list } from './modules/ncov-tools.nf'
include { find_negative_control } from './modules/ncov-tools.nf'
include { create_config_yaml } from './modules/ncov-tools.nf'
include { ncov_tools } from './modules/ncov-tools.nf'
include { combine_all_qc_summaries_for_run } from './modules/ncov-tools.nf'
include { combine_all_lineage_reports_for_run } from './modules/ncov-tools.nf'
include { combine_all_mixture_reports_for_run } from './modules/ncov-tools.nf'
include { combine_ambiguous_position_reports_for_run } from './modules/ncov-tools.nf'
include { get_pangolin_version_for_run } from './modules/ncov-tools.nf'

workflow {

  ch_primer_scheme_version = Channel.of(params.primer_scheme_version)
  ch_primer_scheme_name = Channel.of(params.primer_scheme_name)
  ch_ncov_tools_version = Channel.of(params.ncov_tools_version)
  ch_run_name = Channel.of(params.run_name)
  ch_samplesheet = Channel.fromPath(params.samplesheet, type: 'file')
  ch_metadata = Channel.fromPath(params.metadata, type: 'file')

  update_pangolin(Channel.value(params.update_pangolin))
  download_ncov_tools(ch_ncov_tools_version)
  download_artic_ncov2019(ch_primer_scheme_version.combine(ch_primer_scheme_name))
  index_reference_genome(download_artic_ncov2019.out)

  prepare_data_root(ch_artic_analysis_dir.combine(download_artic_ncov2019.out).combine(ch_metadata))
  
  create_sample_id_list(prepare_data_root.out)
  find_negative_control(prepare_data_root.out)
  create_config_yaml(find_negative_control.out.combine(ch_metadata))

  ncov_tools(create_config_yaml.out.join(prepare_data_root.out).combine(index_reference_genome.out).combine(download_ncov_tools.out).combine(update_pangolin.out))

}
