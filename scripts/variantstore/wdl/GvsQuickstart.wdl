version 1.0

import "GvsAssignIds.wdl" as GvsAssignIds
import "GvsImportGenomes.wdl" as GvsImportGenomes

workflow GvsQuickstart {
    input {
        # Begin GvsAssignIds
        String dataset_name
        String project_id

        Array[String] external_sample_names
        Boolean samples_are_controls = false

        File? gatk_override
        String? service_account_json_path
        # End GvsAssignIds

        # Begin GvsImportGenomes
        Array[File] input_vcfs
        Array[File] input_vcf_indexes

        File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

        # Is this specific version of the gatk still needed? We'll find out!
        # File? load_data_gatk_override = "gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/jars/ah_var_store_20220406/gatk-package-4.2.0.0-480-gb62026a-SNAPSHOT-local.jar"
        # End GvsImportGenomes
    }

    call GvsAssignIds.GvsAssignIds {
        input:
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = external_sample_names,
            samples_are_controls = samples_are_controls,
            assign_ids_gatk_override = gatk_override,
            service_account_json_path = service_account_json_path
    }

    call GvsImportGenomes.GvsImportGenomes {
        input:
            go = GvsAssignIds.gvs_ids_created,
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = external_sample_names,
            input_vcfs = input_vcfs,
            input_vcf_indexes = input_vcf_indexes,
            interval_list = interval_list,
            load_data_preemptible_override = load_data_preemptible_override,
            load_data_gatk_override = gatk_override,
            service_account_json_path = service_account_json_path
    }
}
