from modules.logger_factory import LoggerFactory
import os
from zipfile import ZipFile

# initialize the logger object
logger = LoggerFactory.create_logger("user_output")


class UserOutputModule(object):
    _ZIP_FILE_NAME = "results.zip"
    _LOGS_SUBDIRECTORY = "logs"
    _FINAL_OPTIMIZED_SEQUENCE_FILE_NAME = "optimized_sequence.fasta"

    @staticmethod
    def get_name() -> str:
        return "User Output"

    @classmethod
    def run_module(cls, cds_sequence, zscore, zip_directory=None):
        logger.info('###########################')
        logger.info('# USER OUTPUT INFORMATION #')
        logger.info('###########################')

        logger.info("Output zip file directory path: %s", zip_directory)

        cls._create_final_zip(zip_directory, cds_sequence)

        # TODO - fix the dict according to spec
        return {
            'final_sequence: ': cds_sequence,  # str
            'Zscore': zscore,  # int
            'final_promoter': None,  # str
            'promoter_score': None,  # int
            'promoter_fasta': None,  # fasta file
        }

    @classmethod
    def _create_final_zip(cls, zip_directory, cds_sequence):
        zip_directory = zip_directory if zip_directory else "."
        # Create a ZipFile Object
        zip_file_path = os.path.join(zip_directory, cls._ZIP_FILE_NAME)
        with ZipFile(zip_file_path, 'w') as zip_object:
            # Add multiple files to the zip

            # Add log files to Zip
            cls._write_log_file(zip_object=zip_object, log_file_name="user_input.log")
            cls._write_log_file(zip_object=zip_object, log_file_name="RE.log")
            cls._write_log_file(zip_object=zip_object, log_file_name="user_output.log")

            # Add Fasta files to zip
            zip_object.writestr(cls._FINAL_OPTIMIZED_SEQUENCE_FILE_NAME, cds_sequence)

    @classmethod
    def _write_log_file(cls, zip_object, log_file_name):
        zip_object.write(filename=os.path.join(LoggerFactory.LOG_DIRECTORY, log_file_name),
                         arcname=os.path.join(cls._LOGS_SUBDIRECTORY, log_file_name))
