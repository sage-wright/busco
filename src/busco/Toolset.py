#!/usr/bin/env python
# coding: utf-8
"""
.. module:: Toolset
   :synopsis: the interface to OS enables to run executables / scripts
   in external processes
.. versionadded:: 3.0.0
.. versionchanged:: 4.0.0

Copyright (c) 2016-2017, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
import os
import subprocess
import threading
import time
from shutil import which
from busco.BuscoLogger import BuscoLogger, ToolLogger
from busco.BuscoLogger import LogDecorator as log
from busco.BuscoLogger import StreamLogger
import logging

logger = BuscoLogger.get_logger(__name__)

class Job(threading.Thread):
    """
    Build and executes one work item in an external process
    """

    def __init__(self, tool_name, cmd, job_outlogger, job_errlogger, **kwargs):
        """
        :param name: a name of an executable / script ("a tool") to be run
        :type cmd: list
        :param thread_id: an int id for the thread
        :type thread_id: int
        """
        # initialize parent
        super().__init__()

        self.tool_name = tool_name
        self.cmd_line = [cmd]
        self.job_outlogger = job_outlogger
        self.job_errlogger = job_errlogger
        self.kwargs = kwargs

    def add_parameter(self, parameter):
        """
        Append parameter to the command line
        :parameter: a parameter
        :type parameter: str
        """
        self.cmd_line.append(parameter)

    @log('cmd call: {}', logger, attr_name='cmd_line', apply='join', debug=True)
    def run(self):
        """
        Start external process and block the current thread's execution
        till the process' run is over
        """
        with StreamLogger(logging.DEBUG, self.job_outlogger, **self.kwargs) as out:  # kwargs only provided to out to capture augustus stdout
            with StreamLogger(logging.ERROR, self.job_errlogger) as err:
                # Stick with Popen(), communicate() and wait() instead of just run() to ensure compatibility with
                # Python versions < 3.5.
                p = subprocess.Popen(self.cmd_line, shell=False, stdout=out, stderr=err)
                p.wait()
        self.job_outlogger._file_hdlr.close()
        self.job_outlogger.removeHandler(self.job_outlogger._file_hdlr)
        self.job_errlogger._file_hdlr.close()
        self.job_errlogger.removeHandler(self.job_errlogger._file_hdlr)

class ToolException(Exception):
    """
    Module-specific exception
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


class Tool:
    """
    Collection of utility methods used by all tools
    """

    def __init__(self, name, config, **kwargs):
        """
        Initialize job list for a tool
        :param name: the name of the tool to execute
        :type name: str
        :param config: initialized instance of ConfigParser
        :type config: configparser.ConfigParser
        """

        self.name = name
        self.config = config
        self.cmd = None
        if not self.check_tool_available():
            raise ToolException("{} tool cannot be found. Please check the 'path' and 'command' parameters "
                                "provided in the config file. Do not include the command in the path!".format(self.name))

        self.logfile_path_out = self._set_logfile_path()
        self.logfile_path_err = self.logfile_path_out.replace('_out.log', '_err.log')
        self.kwargs = kwargs
        self.jobs_to_run = []
        self.jobs_running = []
        self.nb_done = 0
        self.total = 0
        self.count_jobs_created = True
        self.logged_header = False

    def check_tool_available(self):
        """
        Check tool's availability.
        1. The section ['name'] is available in the config
        2. This section contains keys 'path' and 'command'
        3. The string resulted from contatination of values of these two keys
        represents the full path to the command
        :param name: the name of the tool to execute
        :type name: str
        :param config: initialized instance of ConfigParser
        :type config: configparser.ConfigParser
        :return: True if the tool can be run, False if it is not the case
        :rtype: bool
        """
        if not self.config.has_section(self.name):
            raise ToolException("Section for the tool [{}] is not present in the config file".format(self.name))

        if not self.config.has_option(self.name, 'path'):
            raise ToolException("Key \'path\' in the section [{}] is not present in the config file".format(self.name))

        if self.config.has_option(self.name, 'command'):
            executable = self.config.get(self.name, 'command')
        else:
            executable = self.name

        self.cmd = os.path.join(self.config.get(self.name, 'path'), executable)

        return which(self.cmd) is not None # True if tool available

    def create_job(self):
        """
        Create one work item
        """
        self.tool_outlogger = ToolLogger(self.logfile_path_out)
        self.tool_errlogger = ToolLogger(self.logfile_path_err)
        job = Job(self.name, self.cmd[:], self.tool_outlogger, self.tool_errlogger, **self.kwargs)
        self.jobs_to_run.append(job)
        if self.count_jobs_created:
            self.total += 1
        return job

    def remove_job(self, job):
        """
        Remove one work item
        :param job: the Job to remove
        :type job: Job
        """
        self.jobs_to_run.remove(job)

    def _set_logfile_path(self):
        return os.path.join(self.config.get("busco_run", "main_out"), "logs", "{}_out.log".format(self.name))

    @log("Running {} job(s) on {}", logger, attr_name=['total', 'name'])
    def log_jobs_to_run(self):
        self.logged_header = True

    def run_jobs(self, max_threads):
        """
        This method run all jobs created for the Tool and redirect
        the standard output and error to the current logger
        :param max_threads: the number or threads to run simultaneously
        :type max_threads: int
        :param log_it: whether to log the progress for the tasks. Default True
        :type log_it: boolean
        """
        if not self.logged_header:
            self.log_jobs_to_run()

        # Wait for all threads to finish and log progress
        already_logged = 0
        while len(self.jobs_to_run) > 0 or len(self.jobs_running) > 0:
            time.sleep(0.001)
            for j in self.jobs_to_run:
                if len(self.jobs_running) < max_threads:
                    self.jobs_running.append(j)
                    self.jobs_to_run.remove(j)
                    j.start()
            for j in self.jobs_running:
                # j.join()
                if not j.is_alive():
                    self.jobs_running.remove(j)
                    self.nb_done += 1

            if (self.nb_done == self.total or int(self.nb_done % float(self.total/10)) == 0) and self.nb_done != already_logged:
                already_logged = self._track_progress()
        # self.total = 0  # Reset for tools that are run twice (tblastn, augustus)


    @log('[{0}]\t{1} of {2} task(s) completed', logger, attr_name=['name', 'nb_done', 'total'], on_func_exit=True)
    def _track_progress(self):
        return self.nb_done