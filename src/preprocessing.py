import os
import pprint

from bioblend import galaxy
from bioblend.galaxy.histories import HistoryClient
from bioblend.galaxy.tools import ToolClient


class PreProcessing: 
    def __init__(self, galaxy_instance: galaxy.GalaxyInstance, data_directory: str = os.getcwd()):
        self.galaxy_instance = galaxy_instance
        self.data_directory = data_directory
        self.current_history = None

    def print(self, content):
        pprint.pprint(content)

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
            return hi

        # return self.galaxy_instance.histories.create_history(name)
        
    def view_histories(self):
        """method to visualize the current histories in account

        Returns:
            list: of dictionaries containing basic metadata, including the id and name of each history
        """
        return self.galaxy_instance.histories.get_histories()

    def show_history(self, history_id=None, contents=True):
        history_client = HistoryClient(self.galaxy_instance)
        return history_client.show_history(str(history_id or self.current_history), contents=contents)

    def upload_data_library(self, library_id, file_local_path, folder_id=None, file_type='auto', dbkey='?'):
        """method to upload local files into data libraries in galaxy

        Args:
            library_id (str): id of the destination Data Library
            file_local_path (str) : path of local file to upload
            folder_id (str, optional): id of the folder where to place the uploaded file. If not provided, the root folder will be used. Defaults to None.
            file_type (str, optional): Galaxy file format name. Defaults to 'auto'.
            dbkey (str, optional): Dbkey. Defaults to '?'.
        """

        self.galaxy_instance.libraries.upload_file_from_local_path(library_id, file_local_path, folder_id=folder_id, file_type=file_type, dbkey=dbkey)
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
        counts = self.galaxy_instance.tools.upload_file(file_local_path, str(history_id or self.current_history), **kwargs)  # ~working
        return counts

    def get_dataset(self, history_id=None, dataset_id=None):
        history_id = str(history_id or self.current_history)
        return self.galaxy_instance.datasets.get_dataset(history_id=history_id, ds_id=dataset_id)

    def get_datasets(self, history_id=None, **kwargs):
        """ method to obtain information about the datasets existing in a given history

        Args:
            history_id: id of history
        Returns:
            dic: dictionary containing information about newly updated file in history
        """
        history_id = str(history_id or self.current_history)
        return self.galaxy_instance.datasets.get_datasets(history_id=history_id, **kwargs)

    def get_jobs(self, history_id=None):
        """ method to obtain information about the datasets existing in a given history

        Args:
        history_id: id of history
        Returns:
        dic: dictionary containing information about newly updated file in history
        """
        return self.galaxy_instance.jobs.get_jobs(str(history_id or self.current_history))

    def set_current_history(self, history_id):
        self.current_history = history_id

    def download_data(self, dataset_id, file_path=None, use_default_filename=True, require_ok_state=True, maxwait=12000):
        # bioblend.galaxy.datasets.DatasetClient(self.galaxy_instance).download_dataset()
        # download_dataset(dataset_id, file_path=None, use_default_filename=True, require_ok_state=True, maxwait=12000
        #
        pass

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
        return self.galaxy_instance.libraries.create_folder(library_id, folder_name, description=description, base_folder_id=base_folder_id)

    def fastqc(self, id_dataset, history_id=None):
        """method to use the tool fastqc

        Args:
            history_id (str): id of history where the we want to run the tool
            id_dataset (str): id of dataset with files to upload to run the tool

        Returns:
            dic: dictionary containing basic information about the used tool
        """
        tool_id = self.galaxy_instance.tools.show_tool("fastqc")["id"]

        return self.galaxy_instance.tools.run_tool(str(history_id or self.current_history), tool_id, {"": id_dataset})

    def rna_star(self, tool_params=None, history_id=None):
        """method to use the tool rna_star

        Args:
            history_id (str): id of history where the tool will be used
            id_dataset (str): id of dataset with files to upload to run the tool

        Returns:
            dic: dictionary containing basic information about the used tool
        """
        if tool_params is None:
            tool_params = {}
        tool_id = self.galaxy_instance.tools.show_tool("rna_star")["id"]
        return self.galaxy_instance.tools.run_tool(str(history_id or self.current_history), tool_id, tool_params)

    def tool_params(self, tool_id):
        tool_show = self.galaxy_instance.tools.show_tool(tool_id=tool_id, io_details=True)
        input_params = tool_show['inputs']
        self.print(input_params)
        return input_params

    def get_histories(self, **kwargs):
        """

        :return:
        """
        return self.galaxy_instance.histories.get_histories(**kwargs)

    def get_job_status(self, job_query, history_id=None):
        """

        :param job:
        :return:
        """
        job_id = job_query['jobs'][0]['id']
        jobs = self.galaxy_instance.jobs.get_jobs(history_id=str(history_id or self.current_history))
        for job in jobs:
            if job['id'] == job_id:
                state = job['state']
                self.print(state)
                return state
        self.print('Job not found!')