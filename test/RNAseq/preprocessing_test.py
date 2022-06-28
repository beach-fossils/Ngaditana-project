import os

import cobra.io
import numpy as np
from bioblend import galaxy
from mewpy.omics import ExpressionSet
from src.integration import *
from src.preprocessing import *
import json
import pprint
from mewpy.omics.integration.gimme import GIMME

def setUp(api_key):
    gi = galaxy.GalaxyInstance(url='https://usegalaxy.org/', key=api_key)
    return PreProcessing(galaxy_instance=gi)




if __name__ == "__main__":
    os.chdir("../../data/")
    preprocessing = setUp('b3a37dd4f27915d58f6bbcfe5b9158f2')
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
    # preprocessing.get_dataset(history_id='77d77febf1e1394a', file_name='fastq_1_ng') #~working
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
    # preprocessing.feature_counts(tool_params={"input1":"bbd44e69cb8906b5235a8c7874a2590b", "GTFfile": annotation[0]['id']})
    files = ["wt.tsv", "nd.tsv"]
    old_dir = os.getcwd()
    os.chdir("../data/inputs")
    for raw_counts in files:
        new_name =raw_counts.replace(".tsv", "")
        preprocessing.calculate_tpm(feature_counts_output=f"fcounts_{raw_counts}", save_path ="./", name_of_file=f"{new_name}_tpm.tsv", gene_length=f"fcounts_length_{raw_counts}")
    os.chdir(old_dir)
    # expr = pd.read_csv('inputs/tpm.tsv', sep='\t')
    # expr["Geneid"] = expr["Geneid"] + "_RA"
    # n_genes = expr.shape[0]
    # identifiers = expr['Geneid'].tolist()
    # conditions = ['tpm']
    # expression = expr['tpm'].to_numpy()[:, np.newaxis]
    # set_expression = ExpressionSet(identifiers, conditions, expression)
    # model = cobra.io.read_sbml_model("inputs/model_ngaditana.xml")
    # # gimme = GIMME(model, set_expression, "e_Biomass__cytop", condition="tpm", parsimonious =True)
    # # print(gimme)
    # eflux  = eFlux(model, set_expression, condition="tpm", parsimonious=True)
    # print(eflux)