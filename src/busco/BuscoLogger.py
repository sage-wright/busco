#!/usr/bin/env python
# coding: utf-8
"""
.. module:: BuscoLogger
   :synopsis: base logger customization for the analysis pipeline
.. versionadded:: 3.0.0
.. versionchanged:: 4.0.0

This is a logger for the pipeline that extends the default Python logger class

Copyright (c) 2016-2017, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""

import logging
import logging.handlers
import sys
import io
import os
import select
import time
import threading
import random

from configparser import NoOptionError
from configparser import NoSectionError


class LogDecorator:

    _log_once_keywords = {}

    def __init__(self, msg, logger,
                 on_func_exit=False, func_arg=None, attr_name=None, iswarn=False, debug=False, apply=None, log_once=False):
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
                if '{' in self.msg and self.on_func_exit:
                    self.retval = func(*args, **kwargs)
                    self.format_string(*args)
                else:
                    self.format_string(*args)
                    self.retval = func(*args, **kwargs)
                return self.retval
            except SystemExit:
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
            if self.apply == 'join' and isinstance(string_arg, tuple):
                string_arg = ' '.join(list(string_arg))
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
                    # string_arg = attr
                    if self.apply == 'join' and isinstance(string_arg, list):
                        string_arg = ' '.join(string_arg)
                    elif self.apply == "basename" and isinstance(string_arg, str):
                        string_arg = os.path.basename(string_arg)
                    log_msg = self.msg.format(string_arg)
                except TypeError:  # if there are multiple attributes specified
                    string_args = (getattr(obj_inst, attr) for attr in self.attr_name)
                    log_msg = self.msg.format(*string_args)

            except AttributeError:
                self.logger.error("No such attribute {}".format(self.attr_name))

            except IndexError:
                self.logger.error("Index out of range for attribute {}".format(self.attr_name))

        elif self.func_arg is not None:
            try:
                string_arg = repr(args[self.func_arg])
                log_msg = self.msg.format(string_arg)

            except IndexError:
                self.logger.error("Index out of range for function argument {}".format(self.func_arg))

        else:
            log_msg = self.msg

        if self.iswarn:
            self.logger.warning(log_msg)
        elif self.debug:
            self.logger.debug(log_msg)
        else:
            self.logger.info(log_msg)
        return


class ToolLogger(logging.getLoggerClass()):

    _level = logging.DEBUG

    def __init__(self, filename):
        super().__init__(filename)
        self.setLevel(type(self)._level)
        self._external_formatter = logging.Formatter('%(message)s')
        self._file_hdlr = logging.FileHandler(filename, mode='a', encoding="UTF-8")
        self._file_hdlr.setFormatter(self._external_formatter)
        self.addHandler(self._file_hdlr)

# The following code was created by combining code based on a SO answer here:
# https://stackoverflow.com/questions/4713932/decorate-delegate-a-file-object-to-add-functionality/4838875#4838875
# with code from the Python docs here:
# https://docs.python.org/3.7/howto/logging-cookbook.html#network-logging
# and is used to log multiple processes to the same file by way of a SocketHandler and separate the streams of
# external tools into stdout and stderr.


class StreamLogger(io.IOBase):
    def __init__(self, level, logger, augustus_out=False):
        self.logger = logger
        self.level = level
        self.augustus_out = augustus_out
        self._run = None
        self.pipe = os.pipe()
        if self.augustus_out:
            self.gene_found = False
            self.output_complete = False
            self.thread = threading.Thread(target=self._flusher_augustus_out)
        else:
            self.thread = threading.Thread(target=self._flusher)
        self.thread.start()

    def _flusher_augustus_out(self):
        self._run = True
        buf = b''
        timeout = 10
        read_only = select.POLLIN | select.POLLPRI | select.POLLHUP | select.POLLERR
        # Switched from select.select() to select.poll() using examples at https://pymotw.com/2/select/
        # This is necessary to handle greater than 1024 file descriptors, but makes BUSCO incompatible with Windows.
        poller = select.poll()
        server = self.pipe[0]
        poller.register(server, read_only)
        while self._run or (self.gene_found and not self.output_complete):
            events = poller.poll(timeout)
            for fd, flag in events:
                if flag & (select.POLLIN | select.POLLPRI):
                    buf += os.read(fd, 4096)
                    while b"\n" in buf:
                        if b"start gene" in buf:
                            self.gene_found = True
                        if b"command line" in buf:
                            self.output_complete = True
                        data, buf = buf.split(b'\n', 1)
                        self.write(data.decode())

        self._run = None

    def _flusher(self):
        self._run = True
        buf = b''
        timeout = 10
        read_only = select.POLLIN | select.POLLPRI | select.POLLHUP | select.POLLERR
        # Switched from select.select() to select.poll() using examples at https://pymotw.com/2/select/
        # This is necessary to handle greater than 1024 file descriptors, but makes BUSCO incompatible with Windows.
        poller = select.poll()
        server = self.pipe[0]
        poller.register(server, read_only)
        while self._run:
            events = poller.poll(timeout)
            for fd, flag in events:
                if flag & (select.POLLIN | select.POLLPRI):
                    buf += os.read(fd, 4096)
                    while b"\n" in buf:
                        data, buf = buf.split(b'\n', 1)
                        self.write(data.decode())

        self._run = None

    def write(self, data):
        return self.logger.log(self.level, data)

    def fileno(self):
        return self.pipe[1]

    def close(self):
        if self._run:
            self._run = False
            while self._run is not None:
                time.sleep(0.01)
            os.close(self.pipe[0])
            os.close(self.pipe[1])


class BuscoLogger(logging.getLoggerClass()):
    """
    This class customizes the _logger class
    """

    _level = logging.DEBUG
    _has_warning = False
    warn_output = io.StringIO()
    random_id = str(random.getrandbits(32))

    def __init__(self, name):
        """
        :param name: the name of the BuscoLogger instance to be created
        :type name: str
        """
        super(BuscoLogger, self).__init__(name)
        self.setLevel(BuscoLogger._level)
        self._normal_formatter = logging.Formatter('%(levelname)s:\t%(message)s')
        self._verbose_formatter = logging.Formatter('%(levelname)s:%(name)s\t%(message)s')
        self._external_formatter = logging.Formatter('%(message)s')

        self._out_hdlr = logging.StreamHandler(sys.stdout)
        self._out_hdlr.addFilter(LessThanFilter(logging.ERROR))
        self._out_hdlr.setLevel(logging.INFO)
        self._out_hdlr.setFormatter(self._normal_formatter)
        self.addHandler(self._out_hdlr)

        self._err_hdlr = logging.StreamHandler()
        self._err_hdlr.setLevel(logging.ERROR)
        self._err_hdlr.setFormatter(self._normal_formatter)
        self.addHandler(self._err_hdlr)

        # Random id used in filename to avoid complications for parallel BUSCO runs.
        self._file_hdlr = logging.FileHandler("busco_{}.log".format(type(self).random_id), mode="a")
        self._file_hdlr.setLevel(logging.DEBUG)
        self._file_hdlr.setFormatter(self._verbose_formatter)
        self.addHandler(self._file_hdlr)

        self._warn_hdlr = logging.StreamHandler(type(self).warn_output)
        self._warn_hdlr.setLevel(logging.WARNING)
        self._warn_hdlr.setFormatter(self._verbose_formatter)
        self.addHandler(self._warn_hdlr)

    def __call__(self):
        pass

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
        #non-zero return means we log this message
        return 1 if record.levelno < self.max_level else 0