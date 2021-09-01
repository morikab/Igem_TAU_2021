from RE.re_main import final_cds_nt, cds_nt
from logger_factory import LoggerFactory

# initialize the logger object
logger = LoggerFactory.create_logger("user_IO")

#todo: to finish and run everything from a mainscript- not using imports.
usr_output = {'final sequence: ': final_cds_nt,  # str
              'Zscore': None,  # int
              'final promoter': None,  # str
              'promoter score': None,  # int
              'promoter fasta': None,  # fasta file
              }

logger.info(final_cds_nt)
logger.info(cds_nt)
