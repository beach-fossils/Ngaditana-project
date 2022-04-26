import os
from bioblend import galaxy

class PreProcessing:
    def __init__(self, galaxy_instance: galaxy.GalaxyInstance, data_directory: str = os.getcwd()):
        self.galaxy_instance = galaxy_instance
        self.data_directory = data_directory

    def upload_data(self, filename):
        """
        Method to upload data to galaxy
        :param filename: name of the file to be uploaded
        :return:
        """
        pass

    def download_data(self):
        pass

    def generate_report(self):
        pass

    def run_job(self):
        pass
