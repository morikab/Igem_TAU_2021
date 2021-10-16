from modules.logger_factory import LoggerFactory
import os
from pathlib import Path
from zipfile import ZipFile

# initialize the logger object
logger = LoggerFactory.create_logger("user_output")


class UserOutputModule(object):
    _ZIP_FILE_NAME = "results.zip"
    _LOGS_SUBDIRECTORY = "logs"
    _FINAL_OPTIMIZED_SEQUENCE_FILE_NAME = "optimized_sequence.fasta"
    _MAST_RESULT_FILE_NAME = "mast_results.html"

    @staticmethod
    def get_name() -> str:
        return "User Output"

    @classmethod
    def run_module(cls,
                   cds_sequence,
                   zscore,
                   weakest_score,
                   p_name,
                   native_prom,
                   synth_promoter,
                   evalue,
                   zip_directory=None):
        logger.info('###########################')
        logger.info('# USER OUTPUT INFORMATION #')
        logger.info('###########################')

        logger.info("Output zip file directory path: %s", zip_directory)

        zip_file_path = cls._create_final_zip(zip_directory, cds_sequence)

        # TODO - fix the dict according to spec
        user_output_dict = {
            'final_sequence: ': cds_sequence,  # str
            'optimization_score': zscore,  # int
            'score_for_weakest_pair': weakest_score, #int
            'promoter_name_and_description': p_name, #str
            'final_promoter': native_prom,  # str
            'promoter_synthetic_version': synth_promoter,
            'promoter_score': evalue,  # int
            'promoter_fasta': None,  # fasta file, should include synthetic and regular options ranked (or just the mast file? might take less time)
        }

        return user_output_dict, zip_file_path

    @classmethod
    def _create_final_zip(cls, zip_directory, cds_sequence) -> str:
        zip_directory = zip_directory if zip_directory else "."
        # Create a ZipFile Object
        zip_file_path = os.path.join(zip_directory, cls._ZIP_FILE_NAME)
        with ZipFile(zip_file_path, 'w') as zip_object:
            # Add multiple files to the zip

            # Add log files to Zip
            cls._write_log_file(zip_object=zip_object, log_file_name="user_input.log")
            cls._write_log_file(zip_object=zip_object, log_file_name="RE.log")
            cls._write_log_file(zip_object=zip_object, log_file_name="user_output.log")
            cls._write_log_file(zip_object=zip_object, log_file_name="promoters.log")

            # Add Fasta files to zip
            zip_object.writestr(cls._FINAL_OPTIMIZED_SEQUENCE_FILE_NAME, cds_sequence)
            zip_object.write(filename=os.path.join(str(Path(LoggerFactory.LOG_DIRECTORY).parent.resolve()),
                                                   cls._MAST_RESULT_FILE_NAME),
                             arcname=cls._MAST_RESULT_FILE_NAME)
        return zip_file_path

    @classmethod
    def _write_log_file(cls, zip_object, log_file_name):
        zip_object.write(filename=os.path.join(LoggerFactory.LOG_DIRECTORY, log_file_name),
                         arcname=os.path.join(cls._LOGS_SUBDIRECTORY, log_file_name))
