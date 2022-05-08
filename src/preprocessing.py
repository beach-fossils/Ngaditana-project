import os
from bioblend import galaxy
from bioblend.galaxy.histories import HistoryClient
from bioblend.galaxy.tools import ToolClient
from bioblend.galaxy.objects import *


class PreProcessing: 
    def __init__(self, galaxy_instance: galaxy.GalaxyInstance, data_directory: str = os.getcwd()):
        self.galaxy_instance = galaxy_instance
        self.data_directory = data_directory

    def create_history(self, name = None):
        """Create a new history in galaxy, optionally setting the name. ID of history will be given!

        Args:
            name (str, optional): name given to the history. Defaults to None.

        Returns:
            dict: dictionary containing information about newly created history
        """
        historyClient = HistoryClient(self.galaxy_instance)
        hi = historyClient.create_history(name)
        self.save = hi
        return hi

        #return self.galaxy_instance.histories.create_history(name)
        

    def view_histories(self):
        """method to visualize the current histories in account

        Returns:
            list: of dictionaries containing basic metadata, including the id and name of each history
        """
        return self.galaxy_instance.histories.get_histories()

    def upload_data_library(self, library_id, file_local_path, folder_id=None, file_type='auto', dbkey='?'):
        """method to upload local files into data libraries in galaxy

        Args:
            library_id (str): id of the destination Data Library
            file_local_path (str) : path of local file to upload
            folder_id (str, optional): id of the folder where to place the uploaded file. If not provided, the root folder will be used. Defaults to None.
            file_type (str, optional): Galaxy file format name. Defaults to 'auto'.
            dbkey (str, optional): Dbkey. Defaults to '?'.
        """

        self.galaxy_instance.libraries.upload_file_from_local_path(library_id, file_local_path)
        #gi.libraries.upload_file_from_local_path('xxxxx', GitHub/bioinformatics-project/data/inputs/xxxx.fastq', file_type='fastq')
        #saber status?

    def upload_data_history(self, file_local_path, history_id):
        """ method to upload a local file to certain history

        Args:
            file_local_path (str): path where our file is located
            history_id (str): id of history where the file will be uploaded
        Returns:
            dic: dictionary containing information about newly updated file in history
        """
        
        toolClient = ToolClient(self.galaxy_instance)
        counts = toolClient.upload_file(file_local_path, history_id) #~working
        return counts
     

    def download_data(self):
        #bioblend.galaxy.datasets.DatasetClient(self.galaxy_instance).download_dataset()
        pass

    def generate_report(self, history_id, tool_inputs):
        """MultiQC is the tool used to perform the task. 
        MultiQC aggregates results from bioinformatics analyses across many samples into a single report

        Args:
            history_id (_type_): _description_
            tool_inputs (_type_): _description_
        """
        pass
        #tool_id = self.galaxy_instance.tools.show_tool("multiqc")["id"]
        #return self.galaxy_instance.tools.run_tool(history_id, tool_id, tool_inputs)


    def run_job(self):    
        pass

    def workflows(self):
        """method to get information on the workflows currently in the account

        Returns:
            list: of metadata dictionaries with detail about the current workflows
        """
        return self.galaxy_instance.GalaxyInstance.workflows.get_workflows()

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
        return self.galaxy_instance.libraries.create_folder(library_id, folder_name, description=None, base_folder_id=None)


    def fastqc(self, history_id, tool_inputs):
        """method to use the tool fastqc

        Args:
            history_id (str): id of history where the tool will be used
            tool_inputs (str): files to be used to run the desired tool

        Returns:
            dic: dictionary containing basic information about the used tool
        """
        tool_id = self.galaxy_instance.tools.show_tool("fastqc")["id"]
        return self.galaxy_instance.tools.run_tool(history_id, tool_id, tool_inputs)

    def star(self, history_id, tool_inputs):
        """method to use the tool star

        Args:
            history_id (str): id of history where the tool will be used
            tool_inputs (str): files to be used to run the desired tool

        Returns:
            dic: dictionary containing basic information about the used tool
        """
        tool_id = self.galaxy_instance.tools.show_tool("star")["id"]
        return self.galaxy_instance.tools.run_tool(history_id, tool_id, tool_inputs)

