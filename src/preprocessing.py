import os
from bioblend import galaxy

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
        return self.galaxy_instance.histories.HistoryClient(name)
        

    def view_histories(self):
        """method to visualize the current histories in account

        Returns:
            list: of dictionaries containing basic metadata, including the id and name of each history
        """
        return self.histories.get_histories()

    def upload_data_library(self, library_id, file_local_path, file_type='auto'):
        """method to upload local files into data libraries in galaxy

        Args:
            library_id (str): id of the destination Data Library
            file_local_path (str): path where our file is located
            file_type (str, optional): default value is 'auto', which Galaxy will attempt to guess the dataset type. Defaults to 'auto'.
        """
        
        self.galaxy.libraries.LibraryClient.upload_file_from_local_path(library_id, file_local_path)
        #gi.libraries.upload_file_from_local_path('xxxxx', GitHub/bioinformatics-project/data/inputs/xxxx.fastq', file_type='fastq')
        #saber status?

    def upload_data_history(self, file_local_path, history_id):
        """ method to upload a local file to certain history
        Args:
            file_local_path (str): path where our file is located
            history_id (str): id of history where the file will be uploaded
        """
        self.tools.upload_from_ftp(file_local_path, history_id)
        
        #gi.tools.upload_from_ftp('testing.fastq', 'xxxxxxxx')

    def download_data(self):
        #bioblend.galaxy.datasets.DatasetClient(self.galaxy_instance).download_dataset()
        pass

    def generate_report(self):
        """
        MultiQC
        :return:
        """
        pass

    def run_job(self):    
        pass

    def workflows(self):
        """method to get information on the workflows currently in the account

        Returns:
            list: of metadata dictionaries with detail about the current workflows
        """
        return self.GalaxyInstance.workflows.get_workflows()

