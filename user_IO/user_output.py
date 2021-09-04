from logger_factory import LoggerFactory
import os
from zipfile import ZipFile

# initialize the logger object
logger = LoggerFactory.create_logger("user_output")


class UserOutputModule(object):
    _ZIP_FILE_NAME = "results.zip"
    _FINAL_OPTIMIZED_SEQUENCE_FILE_NAME = "optimized_sequence.fasta"

    @staticmethod
    def get_name() -> str:
        return "User Output"

    @classmethod
    def run_module(cls, cds_sequence, zscore, zip_directory=None):
        logger.info('###########################')
        logger.info('# USER OUTPUT INFORMATION #')
        logger.info('###########################')

        # TODO - get zip_directory from the user
        cls._create_final_zip(zip_directory, cds_sequence)

        # TODO - fix the dict according to spec
        return {
            'final sequence: ': cds_sequence,  # str
            'Zscore': zscore,  # int
            'final promoter': None,  # str
            'promoter score': None,  # int
            'promoter fasta': None,  # fasta file
        }

    @classmethod
    def _create_final_zip(cls, zip_directory, cds_sequence):
        zip_directory = zip_directory if zip_directory else "."
        # Create a ZipFile Object
        zip_file_path = os.path.join(zip_directory, cls._ZIP_FILE_NAME)
        with ZipFile(zip_file_path, 'w') as zip_object:
            # Add multiple files to the zip

            # Add log files to Zip
            zip_object.write(os.path.join(LoggerFactory.LOG_DIRECTORY, "user_input.log"))
            zip_object.write(os.path.join(LoggerFactory.LOG_DIRECTORY, "RE.log"))
            zip_object.write(os.path.join(LoggerFactory.LOG_DIRECTORY, "user_output.log"))

            # Add Fasta files to zip
            zip_object.writestr(cls._FINAL_OPTIMIZED_SEQUENCE_FILE_NAME, cds_sequence)
