import os
from pathlib import Path
import typing
from zipfile import ZipFile

from modules.logger_factory import LoggerFactory
from modules.shared_functions_and_vars import write_fasta

# initialize the logger object
logger = LoggerFactory.create_logger("user_output")


class UserOutputModule(object):
    _ZIP_FILE_NAME = "communique_results.zip"
    _LOGS_SUBDIRECTORY = "logs"
    _FINAL_OPTIMIZED_SEQUENCE_FILE_NAME = "orf_sequence"
    _FASTA_SUFFIX = ".fasta"

    @staticmethod
    def get_name() -> str:
        return "User Output"

    @classmethod
    def run_module(cls,
                   cds_sequence,
                   zscore,
                   weakest_score,
                   zip_directory=None):
        logger.info('###########################')
        logger.info('# USER OUTPUT INFORMATION #')
        logger.info('###########################')

        logger.info("Output zip file directory path: %s", zip_directory)

        #todo: no need to zip the results anymore
        zip_file_path = cls._create_final_zip(zip_directory=zip_directory,
                                              cds_sequence=cds_sequence)

        user_output_dict = {
            "final_sequence": cds_sequence,  # str
            "optimization_score": zscore,
            "score_for_weakest_pair": weakest_score,
            "zip_file_path": zip_file_path,
        }

        return user_output_dict

    @classmethod
    def _create_final_zip(cls,
                          zip_directory: typing.Optional[str],
                          cds_sequence: str,
                          ) -> str:
        zip_directory = zip_directory if zip_directory else "."
        # Create a ZipFile Object
        zip_file_path = os.path.join(zip_directory, cls._ZIP_FILE_NAME)
        # Add multiple files to the zi
        with ZipFile(zip_file_path, 'w') as zip_object:
            # Add Fasta files
            cls._add_sequence_fasta(zip_object=zip_object, cds_sequence=cds_sequence)

            # Add MEME files

            # Add log files
            cls._write_log_file(zip_object=zip_object, log_file_name="user_input.log")
            cls._write_log_file(zip_object=zip_object, log_file_name="ORF.log")
            cls._write_log_file(zip_object=zip_object, log_file_name="user_output.log")

        return zip_file_path

    @classmethod
    def _write_log_file(cls, zip_object, log_file_name):
        zip_object.write(filename=os.path.join(LoggerFactory.LOG_DIRECTORY, log_file_name),
                         arcname=os.path.join(cls._LOGS_SUBDIRECTORY, log_file_name))

    @classmethod
    def _add_sequence_fasta(cls, zip_object, cds_sequence: str) -> None:
        file_name = os.path.join(str(Path(LoggerFactory.LOG_DIRECTORY).parent.resolve()),
                                 cls._FINAL_OPTIMIZED_SEQUENCE_FILE_NAME)
        write_fasta(fid=file_name, list_seq=(cds_sequence,), list_name=("ORF_sequence", ))
        zip_object.write(filename=cls._add_fasta_suffix(file_name),
                         arcname=cls._add_fasta_suffix(cls._FINAL_OPTIMIZED_SEQUENCE_FILE_NAME))

    @classmethod
    def _add_fasta_suffix(cls, file_name: str) -> str:
        return file_name + cls._FASTA_SUFFIX
