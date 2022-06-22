import os
import pprint

import pandas as pd
from bioblend import galaxy
from bioblend.galaxy.datasets import DatasetClient
from bioblend.galaxy.histories import HistoryClient


# from bioinfokit.analys import norm, get_data

class PreProcessing:
    def __init__(self, galaxy_instance: galaxy.GalaxyInstance, data_directory: str = os.getcwd()):
        self.galaxy_instance = galaxy_instance
        self.data_directory = data_directory
        self.current_history = None

    @staticmethod
    def print(content):
        pprint.pprint(content)

    @staticmethod
    def calculate_tpm(feature_counts_output, save_path, name_of_file, gene_length):
        """
        From the results of the feature_counts tool, calculate the tpm of the genes. The file obtained in feature_counts
        must be downloaded from the history. So it can be used as input for this method.
        :param name_of_file: name of the file that will be saved after the calculation
        :param save_path: path where the tpm file will be saved (local)
        :param feature_counts_output: path to the file obtained from feature_counts
        """
        # df = get_data(feature_counts_output).data
        df = pd.read_csv(feature_counts_output, sep='\t', header=0, index_col=0)

        gene_length_df = pd.read_csv(gene_length, sep='\t', header=0, index_col=0)

        df = pd.merge(df, gene_length_df, left_index=True, right_index=True)

        # make gene column as index column #confirmar isto
        if 'gene' in df.columns:
            df = df.set_index('gene')

        # normalize raw counts using TPM method #confirmar isto
        nm = norm()
        nm.tpm(df=df, gl='Length')
        tpm_df = nm.tpm_norm

        # save tpm_df file locally in computer
        complete_name = os.path.join(save_path, name_of_file)
        tpm_df.columns = ["tpm"]
        tpm_df.to_csv(complete_name, sep='\t', index=True)
        return tpm_df.head(5)

    def create_history(self, name=None):
        """Create a new history in galaxy, optionally setting the name. ID of history will be given!

        Args:
            name (str, optional): name given to the history. Defaults to None.

        Returns:
            dict: dictionary containing information about newly created history
        """
        if self.get_histories(name='try1'):
            print(f"History with name {name} already exists! Please chose another name.")
        else:
            hi = self.galaxy_instance.histories.create_history(name)
            if hi:
                print(f"History with name {name} has been created! Here are the details:{hi}")

        # return self.galaxy_instance.histories.create_history(name)

    def view_histories(self):
        """method to visualize the current histories in account

        Returns:
            list: of dictionaries containing basic metadata, including the id and name of each history
        """
        histories = self.galaxy_instance.histories.get_histories()
        i = 1
        for history in histories:
            if history['deleted'] == False:
                print(f"history {i}. {history['name']}, id: {history['id']}")
                i += 1

    def show_history(self, history_id=None, contents=True):
        history_client = HistoryClient(self.galaxy_instance)

        history = history_client.show_history(history_id or self.current_history, contents=contents)
        print("These files are in history:")
        i = 1
        for item in history:
            print(f"file {i}. name is: {item['name']}, deleted? {item['deleted']}, dataset id: {item['id']}")
            i += 1

    def upload_data_library(self, library_id, file_local_path, folder_id=None, file_type='auto', dbkey='?'):
        """method to upload local files into data libraries in galaxy

        Args:
            library_id (str): id of the destination Data Library
            file_local_path (str) : path of local file to upload
            folder_id (str, optional): id of the folder where to place the uploaded file. If not provided, the root folder will be used. Defaults to None.
            file_type (str, optional): Galaxy file format name. Defaults to 'auto'.
            dbkey (str, optional): Dbkey. Defaults to '?'.
        """

        self.galaxy_instance.libraries.upload_file_from_local_path(library_id, file_local_path, folder_id=folder_id,
                                                                   file_type=file_type, dbkey=dbkey)
        # gi.libraries.upload_file_from_local_path('xxxxx', GitHub/bioinformatics-project/data/inputs/xxxx.fastq', file_type='fastq')
        # saber status?

    def upload_data_history(self, file_local_path, history_id=None, **kwargs):
        """ method to upload a local file to certain history

        Args:
            file_local_path (str): path where our file is located
            history_id (str): id of history where the file will be uploaded
        Returns:
            dic: dictionary containing information about newly updated file in history
        """
        counts = self.galaxy_instance.tools.upload_file(file_local_path, str(history_id or self.current_history),
                                                        **kwargs)  # ~working
        return counts

    def get_dataset_raw(self, history_id=None, dataset_id=None, file_name=None):
        """Method to obtain raw content of a dataset given its name or id to be used as input for another tool. This
        can be useful to automate the process of running different tools."""
        history_id = str(history_id or self.current_history)
        data = self.galaxy_instance.datasets.get_datasets(history_id=history_id)

        if dataset_id is not None:
            for dataset in data:
                if dataset['id'] == dataset_id:
                    return dataset
        elif file_name is not None:
            for dataset in data:
                if dataset['name'] == file_name:
                    return dataset

    def get_dataset(self, history_id=None, dataset_id=None, file_name=None):
        """"method to obtain information of a dataset given its name or id"""
        history_id = str(history_id or self.current_history)

        if dataset_id is not None:
            dataset_id = str(dataset_id)
            data = self.galaxy_instance.datasets
            found = False

            for dataset in data.get_datasets(history_id=history_id):
                if dataset['id'] == dataset_id:
                    print(f"information about the file with the id {dataset_id} was found!")
                    for features in dataset:
                        print(f"{features}: {dataset[features]}")
                    found = True
            if not found:
                print(f"No dataset with id {dataset_id} was found.")

        if file_name is not None:
            file_name = str(file_name)
            data = self.galaxy_instance.datasets
            found = False

            for dataset in data.get_datasets(history_id=history_id):
                if dataset['name'] == file_name:
                    print(f"information about the file with the name {file_name} was found!")
                    for features in dataset:
                        print(f"{features}: {dataset[features]}")
                    found = True
            if not found:
                print(f"No dataset with name {file_name} was found.")

    def get_datasets(self, history_id=None, **kwargs):
        """ method to obtain information about the datasets existing in a given history

        Args:
            history_id: id of history
        Returns:
            dic: dictionary containing information about newly updated file in history
        """
        history_id = str(history_id or self.current_history)
        data = self.galaxy_instance.datasets.get_datasets(history_id=history_id, **kwargs)
        i = 1
        for dataset in data:
            print(
                f"dataset {i}. name is: {dataset['name']}, deleted: {dataset['deleted']}, dataset id: {dataset['id']}, "
                f"created at: {dataset['create_time']}")
            i += 1

        # return self.galaxy_instance.datasets.get_datasets(history_id=history_id, **kwargs)

    #    def get_jobs(self, history_id=None):
    #        """ method to obtain information about the datasets existing in a given history
    #
    #        Args:
    #        history_id (str): id of history
    #        Returns:
    #        dic: dictionary containing information about newly updated file in history
    #        """
    #        return self.galaxy_instance.jobs.get_jobs(str(history_id or self.current_history))

    def set_current_history(self, history_id):
        self.current_history = history_id

    def download_data(self, dataset_id, file_path=None, use_default_filename=True, require_ok_state=True,
                      maxwait=12000):
        """
        Method to download a dataset at Galaxy.
        :param dataset_id: dataset id of the file to be downlowded.
        :param file_path: if this argument is given, the dataset will be downlowded to the chosen path. Otherwise,
        the file will be saved in memory.
        :param use_default_filename: if this parameter is True, the exported file is saved as "file_path/%s", which
        "%s" is the name of the dataset. If False, it is assumed "file_path" contains the complete path file, inlcuding
        the name of the file.
        :param require_ok_state: if False, the dataset is exported even if the "status" is not 'ok'.
        :param maxwait: maximum waiting time, in seconds, to export the file. If time is exceeded, warning
        "DatasetTimeoutException" is thrown.
        :return: if "file_path" is not given, returns a dictionary containing file content. Otherwise nothing will be
        returned.
        """
        data_client = DatasetClient(self.galaxy_instance)
        return data_client.download_dataset(dataset_id, file_path, use_default_filename, require_ok_state, maxwait)

    def generate_report(self, tool_inputs, history_id=None):
        """MultiQC is the tool used to perform the task. 
        MultiQC aggregates results from bioinformatics analyses across many samples into a single report

        Args:
            history_id (_type_): _description_
            tool_inputs (_type_): _description_
        """
        pass
        # tool_id = self.galaxy_instance.tools.show_tool("multiqc")["id"]
        # return self.galaxy_instance.tools.run_tool(history_id, tool_id, tool_inputs)

    def run_job(self):
        pass

    def workflows(self):
        """method to get information on the workflows currently in the account

        Returns:
            list: of metadata dictionaries with detail about the current workflows
        """
        return self.galaxy_instance.workflows.get_workflows()

    def create_folder(self, library_id, folder_name, description=None, base_folder_id=None):
        """Method to create a folder in a library

        Args:
            library_id (str): id of the library where the folder will be created
            folder_name (str): the name to be given to the folder
            description (str, optional): if necessary a description can be made. Defaults to None.
            base_folder_id (str, optional): is the id of the folder where to create the new folder. If not provided, the root folder will be used. Defaults to None.
            
            Returns:
                dict: containing basic information about created folder
        """
        return self.galaxy_instance.libraries.create_folder(library_id, folder_name, description=description,
                                                            base_folder_id=base_folder_id)

    def fastqc(self, id_dataset, history_id=None):
        """method to use the tool fastqc

        Args:
            history_id (str): id of history where to run the tool
            id_dataset (str): id of dataset with files to upload to run the tool

        Returns:
            dic: dictionary containing basic information about the used tool
        """
        tool_id = self.galaxy_instance.tools.show_tool("fastqc")["id"]
        inputs = {'adapters': None,
                  'contaminants': None,
                  'input_file': {'values': [{'id': id_dataset,
                                             'src': 'hda'}]},
                  'kmers': '7',
                  'limits': None,
                  'min_length': '',
                  'nogroup': 'false'}

        results = self.galaxy_instance.tools.run_tool(str(history_id or self.current_history), tool_id, inputs)
        print(f"{results['outputs'][0]['name']} was created in history {results['outputs'][0]['history_id']} has "
              f"id: {results['outputs'][0]['id']}")
        print(f"{results['outputs'][1]['name']} was created in history {results['outputs'][1]['history_id']} has "
              f"id: {results['outputs'][1]['id']}")
        print("More information below:")

        return results['outputs'][0], results['outputs'][1]

    def rna_star(self, tool_params=None, history_id=None):
        """method to use the tool rna_star

        Args:
            history_id (str): id of history where the tool will be used
            tool_params (dict, optional): parameters to be used in the tool. Defaults to None.  If not provided, the
            default parameters will be used.
        Returns:
            dic: dictionary containing basic information about the used tool
        """
        if tool_params is None:
            inputs = {}
        else:
            inputs = {"singlePaired|sPaired": tool_params['SinglePaired'],
                      "singlePaired|input1": {"values": [{"id": tool_params['input1'], "src": "hda"}]},
                      "singlePaired|input2": {"values": [{"id": tool_params['input2'], "src": "hda"}]},
                      "refGenomeSource|geneSource": "history",
                      "genomeFastaFiles": {"values": [{"id": tool_params['FASTAfile'], "src": "hda"}]},
                      "refGenomeSource|genomeSAindexNbases": tool_params['genomeSAindexNbases'],
                      "refGenomeSource|GTFconditional|GTFselect": "with-gtf",
                      "refGenomeSource|sjdbGTFfile": {"values": [{"id": tool_params["GTFfile"], "src": "hda"}]}
                      }
        tool_id = self.galaxy_instance.tools.show_tool("rna_star")["id"]
        return self.galaxy_instance.tools.run_tool(str(history_id or self.current_history), tool_id, inputs)

    def feature_counts(self, tool_params=None, history_id=None):
        """
        Method to use the tool feature_counts.
        :param tool_params: parameters to be used in the tool. If not provided, the default parameters will be used.
        :param history_id: id of history where the tool will be used
        :return: dictionary containing basic information about the used tool
        """
        if tool_params is None:
            tool_params = {}
        inputs = {"alignment": {"values": [{"id": tool_params['input1'], "src": "hda"}]},
                  "anno|anno_select": "history",
                  "anno|reference_gene_sets_builtin": tool_params["GTFfile"]
                  }
        tool_id = self.galaxy_instance.tools.show_tool(
            "toolshed.g2.bx.psu.edu/repos/iuc/featurecounts/featurecounts/2.0.1+galaxy2")["id"]
        return self.galaxy_instance.tools.run_tool(str(history_id or self.current_history), tool_id, inputs)

    def get_histories(self, **kwargs):
        """
        Returns a list of dictionaries containing information about all histories.
        :return:
        """
        return self.galaxy_instance.histories.get_histories(**kwargs)

    def get_job_status(self, job_query, history_id=None):
        """
        Method to get the status of a job (running, ok, error, etc).
        :param  job_query: job_id or job_name
        :param  history_id: history_id where the job is located
        :return: job_status  (str)
        """
        job_id = job_query['jobs'][0]['id']
        jobs = self.galaxy_instance.jobs.get_jobs(history_id=str(history_id or self.current_history))
        for job in jobs:
            if job['id'] == job_id:
                state = job['state']
                self.print(state)
                return state
        self.print('Job not found!')
        return None
