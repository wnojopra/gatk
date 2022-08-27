version 1.0

workflow GvsRepopulateFilterSet {
  input {
    String dataset_name
    String project_id

    String filter_set_name
    File indel_recal_file
    File indel_recal_file_index
    File snp_recal_file
    File snp_recal_file_index
    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index
    File snp_recal_tranches
    File indel_recal_tranches

    Array[String] indel_recalibration_annotation_values = ["AS_FS", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD", "AS_SOR"]
    Array[String] snp_recalibration_annotation_values = ["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"]

    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File gatk_override = "gs://gvs_quickstart_storage/jars/gatk-package-4.2.0.0-552-g0f9780a-SNAPSHOT-local.jar"
  }

  Array[String] snp_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ]

  # reference files
  File axiomPoly_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
  File axiomPoly_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"
  File dbsnp_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
  File dbsnp_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
  File hapmap_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz"
  File hapmap_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi"
  File mills_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  File mills_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
  File omni_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz"
  File omni_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi"
  File one_thousand_genomes_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  File one_thousand_genomes_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
  File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
  File reference_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
  File reference_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

  # fully-qualified table names
  String fq_sample_table = "~{project_id}.~{dataset_name}.sample_info"
  String fq_alt_allele_table = "~{project_id}.~{dataset_name}.alt_allele"
  String fq_info_destination_table = "~{project_id}.~{dataset_name}.filter_set_info"
  String fq_tranches_destination_table = "~{project_id}.~{dataset_name}.filter_set_tranches"
  String fq_filter_sites_destination_table = "~{project_id}.~{dataset_name}.filter_set_sites"

  call PopulateFilterSetInfo {
    input:
      gatk_override = gatk_override,
      filter_set_name = filter_set_name,
      snp_recal_file = snp_recal_file,
      snp_recal_file_index = snp_recal_file_index,
      indel_recal_file = indel_recal_file,
      indel_recal_file_index = indel_recal_file_index,
      fq_info_destination_table = fq_info_destination_table,
      query_project = project_id
  }

  call PopulateFilterSetSites {
    input:
      gatk_override = gatk_override,
      filter_set_name = filter_set_name,
      sites_only_variant_filtered_vcf = sites_only_variant_filtered_vcf,
      sites_only_variant_filtered_vcf_index = sites_only_variant_filtered_vcf_index,
      fq_filter_sites_destination_table = fq_filter_sites_destination_table,
      query_project = project_id
  }

  call PopulateFilterSetTranches {
    input:
      gatk_override = gatk_override,
      filter_set_name = filter_set_name,
      snp_recal_tranches = snp_recal_tranches,
      indel_recal_tranches = indel_recal_tranches,
      fq_tranches_destination_table = fq_tranches_destination_table,
      query_project = project_id
  }

  output {
    Boolean done = true
  }
}

################################################################################
task PopulateFilterSetInfo {
  input {
    String filter_set_name
    String fq_info_destination_table

    File snp_recal_file
    File snp_recal_file_index
    File indel_recal_file
    File indel_recal_file_index

    String query_project

    File? gatk_override
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    set -eo pipefail

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    echo "Creating SNPs reacalibration file"
    gatk --java-options "-Xmx1g" \
    CreateFilteringFiles \
    --ref-version 38 \
    --filter-set-name ~{filter_set_name} \
    -mode SNP \
    -V ~{snp_recal_file} \
    -O ~{filter_set_name}.snps.recal.tsv

    echo "Creating INDELs reacalibration file"
    gatk --java-options "-Xmx1g" \
    CreateFilteringFiles \
    --ref-version 38 \
    --filter-set-name ~{filter_set_name} \
    -mode INDEL \
    -V ~{indel_recal_file} \
    -O ~{filter_set_name}.indels.recal.tsv

    # merge into a single file
    echo "Merging SNP + INDELs"
    cat ~{filter_set_name}.snps.recal.tsv ~{filter_set_name}.indels.recal.tsv | grep -v filter_set_name | grep -v "#"  > ~{filter_set_name}.filter_set_load.tsv

    # BQ load likes a : instead of a . after the project
    bq_table=$(echo ~{fq_info_destination_table} | sed s/\\./:/)

    echo "Loading combined TSV into ~{fq_info_destination_table}"
    bq load --project_id=~{query_project} --skip_leading_rows 0 -F "tab" \
    --range_partitioning=location,0,26000000000000,6500000000 \
    --clustering_fields=location \
    --schema "filter_set_name:string,type:string,location:integer,ref:string,alt:string,vqslod:float,culprit:string,training_label:string,yng_status:string" \
    ${bq_table} \
    ~{filter_set_name}.filter_set_load.tsv > status_load_filter_set_info
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    memory: "3500 MB"
    disks: "local-disk 200 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }

  output {
    String status_load_filter_set_info = read_string("status_load_filter_set_info")
  }
}

task PopulateFilterSetSites {
  input {
    String filter_set_name
    String fq_filter_sites_destination_table

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    String query_project

    File? gatk_override
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    set -eo pipefail

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    echo "Generating filter set sites TSV"
    gatk --java-options "-Xmx1g" \
    CreateSiteFilteringFiles \
    --ref-version 38 \
    --filter-set-name ~{filter_set_name} \
    -V ~{sites_only_variant_filtered_vcf} \
    -O ~{filter_set_name}.filter_sites_load.tsv

    # BQ load likes a : instead of a . after the project
    bq_table=$(echo ~{fq_filter_sites_destination_table} | sed s/\\./:/)

    echo "Loading filter set sites TSV into ~{fq_filter_sites_destination_table}"
    bq load --project_id=~{query_project} --skip_leading_rows 1 -F "tab" \
    --range_partitioning=location,0,26000000000000,6500000000 \
    --clustering_fields=location \
    --schema "filter_set_name:string,location:integer,filters:string" \
    ${bq_table} \
    ~{filter_set_name}.filter_sites_load.tsv > status_load_filter_set_sites
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    memory: "3500 MB"
    disks: "local-disk 200 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }

  output {
    String status_load_filter_set_sites = read_string("status_load_filter_set_sites")

  }
}

task PopulateFilterSetTranches {
  input {
    File? gatk_override

    String filter_set_name
    String fq_tranches_destination_table

    File snp_recal_tranches
    File indel_recal_tranches

    String query_project
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    set -eo pipefail

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    cat ~{snp_recal_tranches} ~{indel_recal_tranches} | grep -v targetTruthSensitivity | grep -v "#" | awk -v CALLSET=~{filter_set_name} '{ print CALLSET "," $0 }' > ~{filter_set_name}.tranches_load.csv

    # BQ load likes a : instead of a . after the project
    bq_table=$(echo ~{fq_tranches_destination_table} | sed s/\\./:/)

    echo "Loading combined tranches CSV into ~{fq_tranches_destination_table}"
    bq load --project_id=~{query_project} --skip_leading_rows 0 -F "," \
    --schema "filter_set_name:string,target_truth_sensitivity:float,num_known:integer,num_novel:integer,known_ti_tv:float,novel_ti_tv:float,min_vqslod:float,filter_name:string,model:string,accessible_truth_sites:integer,calls_at_truth_sites:integer,truth_sensitivity:float" \
    ${bq_table} \
    ~{filter_set_name}.tranches_load.csv > status_load_filter_set_tranches
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    memory: "3500 MB"
    disks: "local-disk 200 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }

  output {
    String status_load_filter_set_tranches = read_string("status_load_filter_set_tranches")
  }
}
