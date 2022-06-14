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
    preprocessing = setUp('879271b2a412ff7df93ae34a5d163eed')
    # preprocessing.create_history('try1')
    #my_history = preprocessing.get_histories(name='try1')

    #preprocessing.set_current_history(my_history[0]['id'])

    #show files in history
    #preprocessing.show_history('77d77febf1e1394a')

    #upload files to history
    #preprocessing.upload_data_history('inputs/toy_SRR5152513_1.fastq', file_name='fastq_1_ng', file_type='fastq',
                                      #history_id='77d77febf1e1394a') #UPLOADED

    #preprocessing.upload_data_history('inputs/toy_SRR5152513_2.fastq', file_name='fastq_2_ng', file_type='fastq',
                                      #history_id='77d77febf1e1394a') #UPLOADED

    #preprocessing.upload_data_history('inputs/genome_sequence.fasta', file_name='sequence_ng', file_type='fasta',
                                     # history_id='77d77febf1e1394a') #UPLOADED

    #preprocessing.upload_data_history("inputs/annotation.gtf", file_name='annotation', file_type='gtf',
                                      #history_id='77d77febf1e1394a') #uploaded
    # show detaisl about existing datasets in history
    #preprocessing.get_datasets(history_id='77d77febf1e1394a')

    # run tools: fastqc
    #
    #fastq_dataset_1 = preprocessing.get_dataset(history_id='77d77febf1e1394a', name="fastq_1_ng")
    #fastq_dataset_2 = preprocessing.get_datasets(history_id='77d77febf1e1394a', name="fastq_2_ng")

    # sequence = preprocessing.get_datasets(history_id=my_history[0]['id'], name="sequence")
    #annotation = preprocessing.get_datasets(history_id=my_history[0]['id'], name="annotation")
    #print(fastq_dataset_1)
    # run fastqc
    #fastqc_1 = preprocessing.fastqc(fastq_dataset_1[0]['id'])
    #fastqc_2 = preprocessing.fastqc(fastq_dataset_2[0]['id'])

    # print about fastqc
    #preprocessing.print(fastqc_1)
    #preprocessing.print(fastqc_2)

    # rna_star =  preprocessing.rna_star(tool_params={"SinglePaired": "paired",
    #                                     "input1": fastq_dataset_1[0]['id'],
    #                                     "input2": fastq_dataset_2[0]['id'],
    #                                     "FASTAfile": sequence[0]['id'],
    #                                     "genomeSAindexNbases": "11",
    #                                     "GTFfile": annotation[0]['id']})
    #dataset =  preprocessing.get_dataset(history_id=my_history[0]['id'], dataset_id="bbd44e69cb8906b5235a8c7874a2590b")
    #preprocessing.feature_counts(tool_params={"input1":"bbd44e69cb8906b5235a8c7874a2590b", "GTFfile": annotation[0]['id']})