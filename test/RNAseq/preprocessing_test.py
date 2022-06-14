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
    # fastq_dataset_1 = preprocessing.get_datasets(history_id=my_history[0]['id'], name="fastq_1")
    # preprocessing.print(fastq_dataset_1)
    # fastq_dataset_2 = preprocessing.get_datasets(history_id=my_history[0]['id'], name="fastq_2")
    # preprocessing.print(fastq_dataset_2)
    # sequence = preprocessing.get_datasets(history_id=my_history[0]['id'], name="sequence")
    annotation = preprocessing.get_datasets(history_id=my_history[0]['id'], name="annotation")
    # fastqc_1 = preprocessing.fastqc(fastq_dataset_1[0]['id'])
    # fastqc_2 = preprocessing.fastqc(fastq_dataset_2[0]['id'])
    # rna_star =  preprocessing.rna_star(tool_params={"SinglePaired": "paired",
    #                                     "input1": fastq_dataset_1[0]['id'],
    #                                     "input2": fastq_dataset_2[0]['id'],
    #                                     "FASTAfile": sequence[0]['id'],
    #                                     "genomeSAindexNbases": "11",
    #                                     "GTFfile": annotation[0]['id']})
    dataset =  preprocessing.get_dataset(history_id=my_history[0]['id'], dataset_id="bbd44e69cb8906b5235a8c7874a2590b")
    preprocessing.feature_counts(tool_params={"input1":"bbd44e69cb8906b5235a8c7874a2590b", "GTFfile": annotation[0]['id']})