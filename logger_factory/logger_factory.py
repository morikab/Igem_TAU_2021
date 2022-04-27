import logging
import os
import typing
from pathlib import Path
from queue import Queue


class QueueHandler(logging.Handler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.log_queue = Queue()

    def emit(self, record: str) -> None:
        self.log_queue.put(record)


class LoggerFactory(object):
    LOG_DIRECTORY = str(os.path.join(Path(__file__).parent.resolve(), "artifacts", "logs"))
    _LOG_SUFFIX = "log"
    _QUEUE_HANDLER_NAME = "queue_handler"

    @classmethod
    def _create_logger(
            cls,
            log_file_name: str,
    ) -> logging.Logger:
        # logger = logging.getLogger(log_file_name)
        logger = logging.getLogger()
        if logger.hasHandlers():
            # Handler is already initialized
            return logger
        logger.setLevel(logging.INFO)
        # -------------
        # Stream Handler
        # -------------
        # TODO - un-comment if running only modules.
        logger.addHandler(logging.StreamHandler())
        # -------------
        # Queue Handler
        # -------------
        queue_handler = QueueHandler()
        queue_handler.name = cls._QUEUE_HANDLER_NAME
        logger.addHandler(queue_handler)
        # -------------
        # File Handler
        # -------------
        # Make sure logs directory exists
        Path(cls.LOG_DIRECTORY).mkdir(parents=True, exist_ok=True)
        # Generate log file
        # TODO - Create a single log file instead of multiple file (use a single logger)
        log_file_path = os.path.join(cls.LOG_DIRECTORY, log_file_name + "." + cls._LOG_SUFFIX)
        log_file_handler = logging.FileHandler(log_file_path, mode="w")
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
        logger = LoggerFactory._create_logger(log_file_name)
        return logger

    @classmethod
    def get_queue_handler(cls, logger: logging.Logger) -> typing.Optional[QueueHandler]:
        for handler in logger.handlers:
            if handler.name == cls._QUEUE_HANDLER_NAME:
                return handler
        return None
