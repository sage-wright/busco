#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: Toolset
   :synopsis: the interface to OS enables to run executables / scripts
   in external processes
.. versionadded:: 3.0.0
.. versionchanged:: 5.4.0

Copyright (c) 2016-2022, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
import os
import subprocess
from subprocess import TimeoutExpired
from multiprocessing import Process, Pool, Value, set_start_method
import time
from abc import ABCMeta, abstractmethod
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.Exceptions import BatchFatalError

logger = BuscoLogger.get_logger(__name__)


class Job(Process):
    """
    Build and executes one work item in an external process
    """

    def __init__(
        self, tool_name, cmd, job_outlogger, job_errlogger, timeout, cwd, **kwargs
    ):
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
        self.timeout = timeout
        self.cwd = cwd
        self.kwargs = kwargs

    def add_parameter(self, parameter):
        """
        Append parameter to the command line
        :parameter: a parameter
        :type parameter: str
        """
        self.cmd_line.append(parameter)

    @log("cmd call: {}", logger, attr_name="cmd_line", apply="join", debug=True)
    def run(self):
        """
        Start external process and block the current thread's execution
        till the process' run is over
        """

        with open(self.job_outlogger, "wb") as f_out:
            with open(self.job_errlogger, "wb") as f_err:
                try:
                    process = subprocess.run(
                        self.cmd_line,
                        stdout=subprocess.PIPE,  # stdout and stderr streams are stored and written to file after job completion
                        stderr=subprocess.PIPE,
                        cwd=self.cwd,
                        shell=False,
                        timeout=self.timeout,
                    )
                    f_out.write(process.stdout)
                    f_err.write(process.stderr)

                except TimeoutExpired:
                    logger.warning(
                        "The following job was killed as it was taking too long (>1hr) to "
                        "complete.\n{}".format(" ".join(self.cmd_line))
                    )

        with cnt.get_lock():
            cnt.value += 1


class Tool(metaclass=ABCMeta):
    """
    Collection of utility methods used by all tools
    """

    set_start_method("spawn", force=True)

    def __init__(self):
        """
        Initialize job list for a tool
        """
        if self.name == "augustus":
            self.kwargs = {"augustus_out": True}
            self.timeout = 3600  # Possibly no longer necessary from 5.2.0 with the new logging system in place,
            # but no harm to leave it here
        else:
            self.kwargs = {}
            self.timeout = None
        self.jobs_to_run = []
        self.jobs_running = []
        self.nb_done = 0
        self.total = 0
        self.cpus = None
        self.chunksize = None
        self.cwd = os.getcwd()

    @abstractmethod
    def configure_job(self, *args):
        pass

    @abstractmethod
    def generate_job_args(self):
        pass

    @property
    @abstractmethod
    def name(self):
        raise NotImplementedError

    @abstractmethod
    def write_checkpoint_file(self):
        pass

    def create_job(self):
        """
        Create one work item
        """
        job = Job(
            self.name,
            self.cmd[:],
            self.logfile_path_out,
            self.logfile_path_err,
            self.timeout,
            self.cwd,
            **self.kwargs
        )
        self.jobs_to_run.append(job)
        return job

    def remove_job(self, job):
        """
        Remove one work item
        :param job: the Job to remove
        :type job: Job
        """
        self.jobs_to_run.remove(job)

    def log_jobs_to_run(self):
        logger.info(
            "Running {} job(s) on {}, starting at {}".format(
                self.total, self.name, time.strftime("%m/%d/%Y %H:%M:%S")
            )
        )
        return

    @log("No jobs to run on {}", logger, attr_name="name", iswarn=True)
    def log_no_jobs(self):
        return

    def run_jobs(self):
        if self.total > 0:
            self.log_jobs_to_run()
        else:
            self.log_no_jobs()
            return

        if self.cpus is None:
            raise BatchFatalError("Number of CPUs not specified.")
        elif self.cpus == 0:
            raise BatchFatalError("Number of CPUs must be greater than 0.")

        with Pool(
            self.cpus, initializer=type(self).init_globals, initargs=(Value("i", 0),)
        ) as job_pool:
            job_pool.map(
                self.run_job, self.generate_job_args(), chunksize=self.chunksize
            )
        self.write_checkpoint_file()

    def run_job(self, args):
        args = (
            (args,) if isinstance(args, str) else tuple(args or (args,))
        )  # Ensure args are tuples that can be unpacked. If no args, args=None, which is falsy,
        # and this evaluates to (None,)
        job = self.configure_job(*args)
        job.run()
        self.nb_done = cnt.value
        if (
            self.nb_done == self.total
            or int(self.nb_done % float(self.total / 10)) == 0
        ):
            self._track_progress()

    @log(
        "[{0}]\t{1} of {2} task(s) completed",
        logger,
        attr_name=["name", "nb_done", "total"],
        on_func_exit=True,
    )
    def _track_progress(self):
        return

    @classmethod
    def init_globals(cls, counter):
        """Counter code adapted from the answer here: https://stackoverflow.com/a/53621343/4844311"""
        global cnt
        cnt = counter
