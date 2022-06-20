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
    # save the information with the ids of the datasets (useful to run the tools)
    #datasets_info = preprocessing.get_datasets(history_id='77d77febf1e1394a')
    #print(datasets_info)

    # the dataset method is using the dataset_id key!
    #preprocessing.get_dataset(history_id='77d77febf1e1394a', dataset_id='bbd44e69cb8906b51a4ae4ab9294c008') #~working
    preprocessing.get_dataset(history_id='77d77febf1e1394a', file_name='fastq_1_ng') #~working
    #preprocessing.get_dataset(history_id='77d77febf1e1394a', file_name='fastq_2_ng') #~working

    # run tools: fastqc
    # id fastq_1_ng = bbd44e69cb8906b5e555944504c5c0e2
    # id fastq_2_ng = bbd44e69cb8906b59b3cb69dc1cd2969
    # at this point you should save the generated files name or save the information into a variable
    # preprocessing.fastqc('bbd44e69cb8906b5e555944504c5c0e2', history_id='77d77febf1e1394a') #~working
    # FastQC on data 74
    # preprocessing.fastqc('bbd44e69cb8906b59b3cb69dc1cd2969', history_id='77d77febf1e1394a') #~working
    # FastQC on data 75
    #fastq_output_1 = preprocessing.fastqc('bbd44e69cb8906b5e555944504c5c0e2', history_id='77d77febf1e1394a')
    #fastq_output_2 = preprocessing.fastqc('bbd44e69cb8906b59b3cb69dc1cd2969', history_id='77d77febf1e1394a')#~working
    # grab the id of the generated files
    #print(fastq_output_1)
    #print(fastq_output_2)

    # print about fastqc ~grab ids of the generated files
    # preprocessing.print(fastq_output_1)
    # preprocessing.print(fastq_output_2)

    # sequence = preprocessing.get_datasets(history_id="77d77febf1e1394a", name="sequence_ng")
    # sequence = preprocessing.get_dataset_raw(history_id="77d77febf1e1394a", file_name="sequence_ng")
    # print(sequence['id'])

    # to get more details: preprocessing.get_dataset(history_id="77d77febf1e1394a", dataset_id='bbd44e69cb8906b509e1d2e79906330c')
    # preprocessing.get_datasets(history_id='77d77febf1e1394a')

    # annotation = preprocessing.get_dataset(history_id="77d77febf1e1394a", file_name="annotation") #~grab the dataset id
    # annotation = preprocessing.get_dataset_raw("77d77febf1e1394a", file_name="annotation")
    # print(annotation['id'])



    #rna_star = preprocessing.rna_star(tool_params={"SinglePaired": "paired",
    #                                        "input1": "bbd44e69cb8906b5933e8be1b873a24c", #fastq1
    #                                        "input2": "bbd44e69cb8906b5400fadedb4c5aeca", #fastq2
    #                                        "FASTAfile": "bbd44e69cb8906b565cde45a6ba4ef9a", #sequence
    #                                        "genomeSAindexNbases": "11",
    #                                        "GTFfile": "bbd44e69cb8906b5694de2059cdb746d"}, #annotation
    #                                  history_id="77d77febf1e1394a")

    #preprocessing.get_datasets("77d77febf1e1394a")
    #preprocessing.print(preprocessing.get_dataset_raw("77d77febf1e1394a", "bbd44e69cb8906b5400fadedb4c5aeca"))

    #rna_star =  preprocessing.rna_star(tool_params={"SinglePaired": "paired",
    #                                     "input1": fastq_dataset_1[0]['id'],
    #                                     "input2": fastq_dataset_2[0]['id'],
    #                                     "FASTAfile": sequence[0]['id'],
    #                                     "genomeSAindexNbases": "11",
    #                                     "GTFfile": annotation[0]['id']})

    #dataset =  preprocessing.get_dataset(history_id=my_history[0]['id'], dataset_id="bbd44e69cb8906b5235a8c7874a2590b")
    #preprocessing.feature_counts(tool_params={"input1":"bbd44e69cb8906b5235a8c7874a2590b", "GTFfile": annotation[0]['id']})