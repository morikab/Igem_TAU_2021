import logging
import os
import typing
from pathlib import Path

import tkinter as tk
from tkinter import scrolledtext

APP_LOGGER_NAME = "app"


class TextHandler(logging.Handler):
    def __init__(self, text: typing.Union[tk.Text, scrolledtext.ScrolledText]):
        logging.Handler.__init__(self)
        self.text = text

    def emit(self, record: str) -> None:
        msg = self.format(record)

        def append() -> None:
            self.text.configure(state="normal")
            self.text.insert(tk.END, msg + '\n')
            self.text.configure(state="disabled")
            # Auto-scroll to the bottom
            self.text.yview(tk.END)

        self.text.after(0, append)


class LoggerFactory(object):
    LOG_DIRECTORY = str(os.path.join(Path(__file__).parent.resolve(), "artifacts", "logs"))
    _LOG_SUFFIX = "log"

    @classmethod
    def _create_logger(
            cls,
            log_file_name: str,
            text: typing.Optional[typing.Union[tk.Text, scrolledtext.ScrolledText]] = None,
    ) -> logging.Logger:
        """
        A private method that interacts with the python
        logging module
        """
        # TODO - try to run without adding the file_name..? Do we always get the same logger that way?
        logger = logging.getLogger(log_file_name)
        if logger.handlers:
            # If we already set the handlers for the logger, just return the initialized logger.
            return logger
        logger.setLevel(logging.INFO)
        logger.addHandler(logging.StreamHandler())    # TODO - un-comment if running only modules.
        if text is not None:
            logger.addHandler(TextHandler(text))
        # Make sure logs directory exists
        Path(cls.LOG_DIRECTORY).mkdir(parents=True, exist_ok=True)
        # Generate log file
        log_file_path = os.path.join(cls.LOG_DIRECTORY, log_file_name + "." + cls._LOG_SUFFIX)
        log_file_handler = logging.FileHandler(log_file_path, mode="w")
        log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        log_file_handler.setFormatter(logging.Formatter(log_format))
        logger.addHandler(log_file_handler)
        return logger

    @staticmethod
    def create_logger(log_file_name: str,
                      text: typing.Optional[typing.Union[tk.Text, scrolledtext.ScrolledText]] = None) -> logging.Logger:
        """
        A static method called by other modules to initialize logger in
        their own module
        """
        logger = LoggerFactory._create_logger(log_file_name, text)

        # Add TextHandler
        app_logger = logging.getLogger(APP_LOGGER_NAME)
        for handle in app_logger.handlers:
            if isinstance(handle, TextHandler):
                logger.addHandler(handle)
                break

        return logger
