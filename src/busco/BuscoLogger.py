#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: BuscoLogger
   :synopsis: base logger customization for the analysis pipeline
.. versionadded:: 3.0.0
.. versionchanged:: 5.4.0

This is a logger for the pipeline that extends the default Python logger class

Copyright (c) 2016-2022, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""

import logging
import logging.handlers
import sys
import io
import os
from busco.Exceptions import BatchFatalError, BuscoError

from configparser import NoOptionError
from configparser import NoSectionError


class LogDecorator:

    _log_once_keywords = {}

    def __init__(
        self,
        msg,
        logger,
        on_func_exit=False,
        func_arg=None,
        attr_name=None,
        iswarn=False,
        debug=False,
        apply=None,
        log_once=False,
    ):
        self.msg = msg
        self.logger = logger
        self.on_func_exit = on_func_exit
        self.func_arg = func_arg
        self.attr_name = attr_name
        self.iswarn = iswarn
        self.debug = debug
        self.apply = apply
        self.log_once = log_once
        self.retval = None

    def __call__(self, func):
        def wrapped_func(*args, **kwargs):
            try:
                if "{" in self.msg and self.on_func_exit:
                    self.retval = func(*args, **kwargs)
                    self.format_string(*args)
                else:
                    self.format_string(*args)
                    self.retval = func(*args, **kwargs)
                return self.retval
            except (BuscoError, BatchFatalError):
                raise

        return wrapped_func

    def format_string(self, *args):
        if self.log_once:
            if self.attr_name in type(self)._log_once_keywords:
                return
            else:
                type(self)._log_once_keywords[self.attr_name] = 1

        if self.attr_name == "retvalue":
            string_arg = self.retval
            if self.apply == "join" and isinstance(string_arg, tuple):
                string_arg = " ".join(list(string_arg))
            elif self.apply == "basename" and isinstance(string_arg, str):
                string_arg = os.path.basename(string_arg)
            log_msg = self.msg.format(string_arg)

        elif self.attr_name is not None:
            try:
                try:
                    obj_inst = args[0]
                except IndexError:
                    self.logger.error("No arguments passed to function.")
                    return

                try:
                    string_arg = getattr(obj_inst, self.attr_name)
                    if self.apply == "join" and isinstance(string_arg, list):
                        string_arg = [
                            str(arg) for arg in string_arg
                        ]  # Ensure all parameters are joinable strings
                        string_arg = " ".join(string_arg)
                    elif self.apply == "basename" and isinstance(string_arg, str):
                        string_arg = os.path.basename(string_arg)

                    log_msg = self.msg.format(string_arg)
                except TypeError:  # if there are multiple attributes specified
                    string_args = (getattr(obj_inst, attr) for attr in self.attr_name)
                    log_msg = self.msg.format(*string_args)

            except AttributeError:
                self.logger.error("No such attribute {}".format(self.attr_name))

            except IndexError:
                self.logger.error(
                    "Index out of range for attribute {}".format(self.attr_name)
                )

        elif self.func_arg is not None:
            try:
                string_arg = repr(args[self.func_arg])
                log_msg = self.msg.format(string_arg)

            except IndexError:
                self.logger.error(
                    "Index out of range for function argument {}".format(self.func_arg)
                )

        else:
            log_msg = self.msg

        if self.iswarn:
            self.logger.warning(log_msg)
        elif self.debug:
            self.logger.debug(log_msg)
        else:
            self.logger.info(log_msg)
        return


class BuscoLogger(logging.getLoggerClass()):
    """
    This class customizes the _logger class
    """

    _level = logging.DEBUG
    _has_warning = False
    warn_output = io.StringIO()
    ppid = str(os.getppid())
    pid = str(os.getpid())
    quiet = False
    quiet_msg_logged = False

    def __init__(self, name):
        """
        :param name: the name of the BuscoLogger instance to be created
        :type name: str
        """
        super(BuscoLogger, self).__init__(name)
        self.setLevel(BuscoLogger._level)
        self._normal_formatter = logging.Formatter(
            "%(asctime)s %(levelname)s:\t%(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
        self._verbose_formatter = logging.Formatter(
            "%(asctime)s %(levelname)s:%(name)s\t%(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        self._external_formatter = logging.Formatter("%(message)s")

        self._out_hdlr = logging.StreamHandler(sys.stdout)
        self._out_hdlr.addFilter(LessThanFilter(logging.ERROR))
        self._out_hdlr.setLevel(logging.INFO)
        self._out_hdlr.setFormatter(self._normal_formatter)
        self.addHandler(self._out_hdlr)

        self._err_hdlr = logging.StreamHandler()
        self._err_hdlr.setLevel(logging.ERROR)
        self._err_hdlr.setFormatter(self._normal_formatter)
        self.addHandler(self._err_hdlr)

        try:
            # Identify main log file with process ID and have all spawned processes log to the correct file
            log_filename = "busco_{}.log".format(type(self).ppid)
            if not os.path.exists(log_filename):
                log_filename = "busco_{}.log".format(type(self).pid)

            # Process id used in filename to avoid complications for parallel BUSCO runs.
            self._file_hdlr = logging.FileHandler(log_filename, mode="a")
        except IOError as e:
            errStr = (
                "No permission to write in the current directory: {}".format(
                    os.getcwd()
                )
                if e.errno == 13
                else "IO error({0}): {1}".format(e.errno, e.strerror)
            )
            raise BatchFatalError(errStr)

        self._file_hdlr.setLevel(logging.DEBUG)
        self._file_hdlr.setFormatter(self._verbose_formatter)
        self.addHandler(self._file_hdlr)

        self._warn_hdlr = logging.StreamHandler(type(self).warn_output)
        self._warn_hdlr.setLevel(logging.WARNING)
        self._warn_hdlr.setFormatter(self._verbose_formatter)
        self.addHandler(self._warn_hdlr)

    def __call__(self):
        pass

    def _log(self, *args, **kwargs):
        if type(self).quiet:
            if not type(self).quiet_msg_logged:
                super()._log(
                    20,
                    "Quiet mode selected. All subsequent INFO log messages to stdout suppressed. "
                    "Detailed log still available in logs/busco.log",
                    [],
                )
                type(self).quiet_msg_logged = True
            self.removeHandler(self._out_hdlr)
        super()._log(*args, **kwargs)

    @classmethod
    def reset(cls):
        cls._level = logging.DEBUG
        cls._has_warning = False
        cls.warn_output = io.StringIO()
        return

    @staticmethod
    def get_logger(name, config=None):
        """
        :param name: the name of the logger to be returned
        :type name: str
        :param config: the parameters of the analysis
        :type config: PipeConfig
        :return: a BuscoLogger, new or existing, corresponding to the provided name
        :rtype: BuscoLogger
        """
        try:
            if config and config.getboolean("busco_run", "quiet"):
                BuscoLogger._level = logging.ERROR
        except NoOptionError:
            pass
        except NoSectionError:
            pass

        logging.setLoggerClass(BuscoLogger)
        return logging.getLogger(name)

    def warn(self, msg, *args, **kwargs):
        """
        This function redirects the obsolete logging class method "warn"
        :param msg: the message to log
        :type msg: str
        """
        self.warning(msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        """
        This function overrides the _logger class warning
        :param msg: the message to log
        :type msg: str
        """
        type(self)._has_warning = True
        super().warning(msg, *args, **kwargs)

    def has_warning(self):
        """
        :return: whether any _logger encountered any log warnings
        :rtype: boolean
        """
        return type(self)._has_warning


# Code from https://stackoverflow.com/a/31459386/4844311


class LessThanFilter(logging.Filter):
    def __init__(self, exclusive_maximum, name=""):
        super(LessThanFilter, self).__init__(name)
        self.max_level = exclusive_maximum

    def filter(self, record):
        # non-zero return means we log this message
        return 1 if record.levelno < self.max_level else 0
