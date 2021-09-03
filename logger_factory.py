import logging
import os
from pathlib import Path


class LoggerFactory(object):
    _LOG_DIRECTORY = "logs"
    _LOG_SUFFIX = "log"

    @classmethod
    def __create_logger(cls, log_file_name: str) -> logging.Logger:
        """
        A private method that interacts with the python
        logging module
        """
        logger = logging.getLogger(log_file_name)
        if logger.handlers:
            # If we already set the handlers for the logger, just return the initialized logger.
            return logger
        logger.setLevel(logging.INFO)
        logger.addHandler(logging.StreamHandler())
        # Make sure logs directory exists
        Path(cls._LOG_DIRECTORY).mkdir(parents=True, exist_ok=True)
        # Generate log file
        log_file_path = os.path.join(cls._LOG_DIRECTORY, log_file_name + "." + cls._LOG_SUFFIX)
        log_file_handler = logging.FileHandler(log_file_path)
        log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        log_file_handler.setFormatter(logging.Formatter(log_format))
        logger.addHandler(log_file_handler)
        return logger

    @staticmethod
    def create_logger(log_file_name: str) -> logging.Logger:
        """
        A static method called by other modules to initialize logger in
        their own module
        """
        logger = LoggerFactory.__create_logger(log_file_name)
        return logger
