import argparse
import json
import re
from collections import defaultdict


def generate_avro_dict(avro_prefix, gcs_listing):
    """
    Generate a dictionary of the Avro arguments for the `hl.import_gvs` invocation. The keys will match the formal
    parameters of the `hl.import_gvs` function and should include:

    * vets (list of lists, one outer list per GVS superpartition of 4000 samples max)
    * refs (list of lists, one outer list per GVS superpartition of 4000 samples max)
    * sample_mapping_data (list)
    * site_filtering_data (list)
    * vqsr_filtering_data (list)
    * vqsr_tranche_data (list)
    """

    avro_file_arguments = defaultdict(list)
    superpartitioned_keys = {'vets', 'refs'}
    slashed_avro_prefix = avro_prefix + "/" if not avro_prefix[-1] == '/' else avro_prefix

    for full_path in gcs_listing.readlines():
        full_path = full_path.strip()
        if not full_path.endswith(".avro"):
            continue
        relative_path = full_path[len(slashed_avro_prefix):]
        parts = relative_path.split('/')
        key = parts[0]

        if key in superpartitioned_keys:
            # Get the zero based index from the `vet_001` component
            index = int(parts[1].split('_')[-1]) - 1
            if len(avro_file_arguments[key]) == index:
                avro_file_arguments[key].append([])
            avro_file_arguments[key][index].append(full_path)
        else:
            avro_file_arguments[key].append(full_path)

    return dict(avro_file_arguments)


def generate_avro_args(avro_dict):
    """
    Stringifies an Avro arguments dictionary such that it can be interpolated into a `hl.gvs_import` invocation to
    encompass all required Avro parameters.
    """
    avro_args = []
    for k, v in avro_dict.items():
        s = f'{k}={json.dumps(v, indent=4)}'
        avro_args.append(s)
    return ',\n'.join(avro_args)


def generate_gvs_import_script(avro_prefix, gcs_listing, vds_output_path, temp_dir):
    """
    Generate a Python script that when executed will:

    * import GVS Avro files into a Hail VDS
    * export from the Hail VDS to VCF

    """
    avro_dict = generate_avro_dict(avro_prefix, gcs_listing)
    avro_args = generate_avro_args(avro_dict)
    indented_avro_args = re.sub('\n', '\n    ', avro_args)
    hail_script = f"""
# The following instructions can be used from the terminal of a Terra notebook to import GVS QuickStart Avro files
# and generate a VDS.
#
# Copy the appropriate Hail wheel locally first:
#
# gsutil -m cp 'gs://gvs-internal-scratch/hail-wheels/2022-08-18-01f7b77ebbcc/hail-0.2.97-py3-none-any.whl' .
#
# If running locally (non-Spark cluster) set this in the environment before launching Python:
# export PYSPARK_SUBMIT_ARGS='--driver-memory 16g --executor-memory 16g pyspark-shell'
#
# Hail wants Java 8, Java 11+ will not do. Make sure you have a Java 8 in your path with `java -version`.
#
# pip install --force-reinstall hail-0.2.97-py3-none-any.whl
# gcloud auth application-default login
# curl -sSL https://broad.io/install-gcs-connector | python3

import hail as hl

rg38 = hl.get_reference('GRCh38')
rg38.add_sequence('gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz',
                  'gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai')

hl.import_gvs(
    {indented_avro_args},
    final_path="{vds_output_path}",
    tmp_dir="{temp_dir}",
    reference_genome=rg38,
)

"""
    return hail_script



def generate_vat_inputs_script(vds_output_path, vcf_output_path, sites_only_vcf_output_path, vat_custom_annotations_tsv_path, ancestry_file_path):
    """
    Generate a Python script that when executed will:

    * import GVS Avro files into a Hail VDS
    * export from the Hail VDS to VCF

    """
    hail_script = f"""
    
## Get the original VDS    
vds = hl.vds.read_vds('{vds_output_path}')

from datetime import datetime
start = datetime.now()
current_time = start.strftime("%H:%M:%S")
print("Start Time =", current_time)

# Now parse the ancestry file to get it ready for the subpopulation work

import csv

sample_id_to_sub_population_map = {}
with open('{ancestry_file_path}', 'r') as file:
  reader = csv.reader(file, delimiter='\t')
  next(reader) # skip header
  for row in reader:
    key = row[0]
    value = row[4]
    sample_id_to_sub_population_map[key] = value

# print(sample_id_to_sub_population_map)
# need to make it looks like: sample_id_to_sub_population_map = {"ERS4367795":"eur","ERS4367796":"eas","ERS4367797":"eur","ERS4367798":"afr","ERS4367799":"oth","ERS4367800":"oth","ERS4367801":"oth","ERS4367802":"oth","ERS4367803":"oth","ERS4367804":"oth","ERS4367805":"oth"}

## * Hard filter out non-passing sites !!!!!TODO ok wait, DO WE want to do this? I guess it depends on what we are handing off and when.
# note: no AC/AN and AF for filtered out positions
vd = vds.variant_data
filtered_vd = vd.filter_rows(hl.len(vd.filters)==0)
filtered_vds = hl.vds.VariantDataset(vds.reference_data, filtered_vd) # now we apply it back to the vds 

## * Replace LGT with GT ( for easier calculations later )
filtered_vd = filtered_vds.variant_data
filtered_vd = filtered_vd.annotate_entries(GT=hl.vds.lgt_to_gt(filtered_vd.LGT, filtered_vd.LA) )
filtered_vds = hl.vds.VariantDataset(filtered_vds.reference_data, filtered_vd)

## * Respect the FT flag by setting all failing GTs to a no call
# Logic for assigning non passing GTs as no-calls to ensure that AC,AN and AF respect the FT flag
# filtered_vd.FT is True ⇒ GT keeps its current value
# filtered_vd.FT is False ⇒ GT assigned no-call
# filtered_vd.FT is missing ⇒ GT keeps its current value

filtered_vd = filtered_vd.annotate_entries(GT=hl.or_missing(hl.coalesce(filtered_vd.FT, True), filtered_vd.GT))
# TODO drop GT after it is used since it will not line up with the LGT


## * Turn the GQ0s into no calls so that ANs are correct
rd = filtered_vds.reference_data
rd = rd.filter_entries(rd.GQ > 0) ## would be better to drop these once its a dense mt? 
filtered_vds = hl.vds.VariantDataset(rd, filtered_vd)

## * Create a DENSE MATRIX TABLE to calculate AC, AN, AF and TODO: Sample Count
mt = hl.vds.to_dense_mt(filtered_vds)

mt = mt.annotate_cols(pop = hl.literal(sample_id_to_sub_population_map)[mt.s])
mt = mt.select_rows(
    ac_an_af =  hl.agg.call_stats(mt.GT, mt.alleles), 
    call_stats_by_pop = hl.agg.group_by(mt.pop, hl.agg.call_stats(mt.GT, mt.alleles))
)

qc_data = mt.rows() ## this is what we will use to create the TSV for Nirvana custom annotations

## Now we join this back to the VDS
filtered_vd = filtered_vds.variant_data
filtered_vd = filtered_vd.annotate_rows(ac_an_af = qc_data[filtered_vd.row_key])
final_vds = hl.vds.VariantDataset(rd, filtered_vd)

## save the VDS to GCP
(unclear that we even need this if we are only going to create a VAT with what is left)
# ???

## Create the VAT inputs:

# 1. Create a Sites only VCF
# TODO we will want to drop some cols because there's no reason to carry around some of this stuff
hl.export_vcf(final_vds.variant_data.rows(), '{sites_only_vcf_output_path}')

# 2. Create a VAT TSV file with subpopulation data
# TODO create the desired TSV-- We're gonna pull in Tim for this
# Do we need to track the dropped sites?
# filter out sites with too many alt alleles and trim extraneous INFO and FORMAT fields
# normalize, left align and split multi allelic sites to new lines, remove duplicate lines
# filter out spanning deletions and variants with an AC of 0
# drop GTs if we have not already since LGT doesn't match
# figure out how to calculate the Sample Count (is it just AC_ref + AC_var?)
# 
hl.export_table(qc_data, '{vat_custom_annotations_tsv_path}')



## The following can be used for extracting a VCF
# mt = hl.vds.to_dense_mt(vds)
# fail_case = 'FAIL'
# mt = mt.annotate_entries(FT=hl.if_else(mt.FT, 'PASS', fail_case))
# hl.export_vcf(mt, '{vcf_output_path}')
"""
    return hail_script


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Generate Hail gvs_import invocation script')

    parser.add_argument('--avro_prefix', type=str, help='Avro prefix', required=True)
    parser.add_argument('--avro_listing_file', type=str, help='File containing a recursive listing under `avro_prefix',
                        required=True)
    parser.add_argument('--vds_output_path', type=str, help='GCS location for VDS output', required=True)
    parser.add_argument('--vcf_output_path', type=str, help='GCS location for VCF output generated from VDS',
                        required=True)
    parser.add_argument('--sites_only_vcf_output_path', type=str, help='GCS location for Sites Only VCF output generated from VDS',
                        required=True)
    parser.add_argument('--gcs_temporary_path', type=str, help='GCS location under which to create temporary files',
                        required=True)

    args = parser.parse_args()

    with open(args.avro_listing_file, 'r') as listing:
        script = generate_gvs_import_script(args.avro_prefix,
                                            listing,
                                            args.vds_output_path,
                                            args.vcf_output_path,
                                            args.gcs_temporary_path)
        print(script)
        vat_script = generate_vat_inputs_script(args.vds_output_path,
                                                args.vcf_output_path,
                                                args.sites_only_vcf_output_path,
                                                args.vat_custom_annotations_tsv_path,
                                                args.ancestry_file_path)
        print(vat_script)

