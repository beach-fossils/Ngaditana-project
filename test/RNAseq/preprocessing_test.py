import os

from bioblend import galaxy

from src.preprocessing import *
import json
import pprint

def setUp(api_key):
    gi = galaxy.GalaxyInstance(url='https://usegalaxy.org/', key=api_key)
    return PreProcessing(galaxy_instance = gi)






if __name__ == "__main__":
    os.chdir("../../data/")
    preprocessing = setUp('b3a37dd4f27915d58f6bbcfe5b9158f2')
    # preprocessing.create_history('try1')
    my_history = preprocessing.get_histories(name='try1')
    preprocessing.set_current_history(my_history[0]['id'])
    # preprocessing.upload_data_history("inputs/toy_SRR5152513_1.fastq", file_name = 'fastq_1')
    # preprocessing.upload_data_history("inputs/toy_SRR5152513_2.fastq", file_name = 'fastq_2')
    # preprocessing.upload_data_history("inputs/genome_sequence.fasta", file_name='sequence')
    # preprocessing.upload_data_history("inputs/genome_annotation.gff3", file_name='annotation', file_type='gff3')
    fastq_dataset_1 = preprocessing.get_datasets(history_id=my_history[0]['id'], name="fastq_1")
    preprocessing.print(fastq_dataset_1)
    # fastq_dataset_2 = preprocessing.get_datasets(history_id=my_history[0]['id'], name="fastq_2")
    # preprocessing.print(fastq_dataset_2)
    # sequence = preprocessing.get_datasets(history_id=my_history[0]['id'], name="sequence")
    # annotation = preprocessing.get_datasets(history_id=my_history[0]['id'], name="annotation")
    # print()
    fastqc_1 = preprocessing.fastqc(fastq_dataset_1[0]['id'])
    # fastqc_1 = {'outputs': [{'id': 'bbd44e69cb8906b52a8c161a85ec8aae', 'hda_ldda': 'hda', 'uuid': '9a792d22-a46a-4b45-9ceb-49dbdfa6389a', 'hid': 7, 'file_ext': 'html', 'peek': None, 'model_class': 'HistoryDatasetAssociation', 'name': 'FastQC on data 2: Webpage', 'deleted': False, 'purged': False, 'visible': True, 'state': 'new', 'history_content_type': 'dataset', 'file_size': 0, 'create_time': '2022-05-25T09:46:17.606439', 'update_time': '2022-05-25T09:46:17.606449', 'data_type': 'galaxy.datatypes.text.Html', 'genome_build': '?', 'validated_state': 'unknown', 'validated_state_message': None, 'misc_info': None, 'misc_blurb': 'queued', 'tags': [], 'history_id': 'b47d6ac96387d407', 'metadata_dbkey': '?', 'metadata_data_lines': 0, 'output_name': 'html_file'}, {'id': 'bbd44e69cb8906b5bbac704f157a0eb6', 'hda_ldda': 'hda', 'uuid': 'ec563a08-3273-4dbf-9752-bbbc4ff24c9a', 'hid': 8, 'file_ext': 'txt', 'peek': None, 'model_class': 'HistoryDatasetAssociation', 'name': 'FastQC on data 2: RawData', 'deleted': False, 'purged': False, 'visible': True, 'state': 'new', 'history_content_type': 'dataset', 'file_size': 0, 'create_time': '2022-05-25T09:46:17.606457', 'update_time': '2022-05-25T09:46:17.606462', 'data_type': 'galaxy.datatypes.data.Text', 'genome_build': '?', 'validated_state': 'unknown', 'validated_state_message': None, 'misc_info': None, 'misc_blurb': 'queued', 'tags': [], 'history_id': 'b47d6ac96387d407', 'metadata_dbkey': '?', 'metadata_data_lines': 0, 'output_name': 'text_file'}], 'output_collections': [], 'jobs': [{'model_class': 'Job', 'id': 'bbd44e69cb8906b5eba058e49165309d', 'state': 'new', 'exit_code': None, 'update_time': '2022-05-25T09:46:17.658634', 'create_time': '2022-05-25T09:46:17.552613', 'galaxy_version': '22.01', 'tool_id': 'toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.73+galaxy0', 'history_id': 'b47d6ac96387d407'}], 'implicit_collections': [], 'produces_entry_points': False}
    # preprocessing.get_job_status(fastqc_1)
    # fastqc_2 = preprocessing.fastqc(fastq_dataset_2[0]['id'])
    # preprocessing.get_job_status(fastqc_2)
    # tool_id = preprocessing.galaxy_instance.tools.show_tool("rna_star")["id"]
    # preprocessing.tool_params(tool_id)
    # preprocessing.rna_star(tool_params={
    # }
    #                        )