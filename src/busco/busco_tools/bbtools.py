# coding: utf-8
"""
bbtools.py

Module for running BBTools.

Author(s): Matthew Berkeley

Copyright (c) 2015-2024, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

from busco.busco_tools.base import BaseRunner
from busco.BuscoLogger import BuscoLogger
import os
import subprocess

logger = BuscoLogger.get_logger(__name__)


class BBToolsRunner(BaseRunner):

    name = "bbtools"
    cmd = "stats.sh"
    metrics = {}

    def __init__(self):
        super().__init__()
        self.contig_break = self.config.get("busco_run", "contig_break")
        self.scaffold_composition = self.config.getboolean(
            "busco_run", "scaffold_composition"
        )
        self._output_folder = os.path.join(
            self.run_folder,
            "{}bbtools_output".format("." * int(not self.scaffold_composition)),
        )
        self.create_dirs(self._output_folder)

    def configure_runner(self):
        super().configure_runner()
        self.run_number += 1

    def configure_job(self, *args):
        bbtools_job = self.create_job()
        bbtools_job.add_parameter("format=2")
        bbtools_job.add_parameter("in={}".format(self.input_file))
        bbtools_job.add_parameter("n={}".format(self.contig_break))
        if self.scaffold_composition:
            bbtools_job.add_parameter(
                "gc={}".format(
                    os.path.join(self._output_folder, "scaffold_composition.txt")
                )
            )
        return bbtools_job

    def check_tool_dependencies(self):
        pass

    def get_version(self):
        try:
            bbtools_version = subprocess.check_output(
                [self.cmd, "--version"], stderr=subprocess.STDOUT, shell=False
            )
            lines = bbtools_version.decode("utf-8").split("\n")
        except subprocess.CalledProcessError as cpe:
            if cpe.output.decode("utf-8").startswith("BBMap version"):
                lines = cpe.output.decode("utf-8").split("\n")
            else:
                raise

        for line in lines:
            if line.startswith("BBMap version"):
                version = line.split("BBMap version ")[1]
                return version

    def generate_job_args(self):
        yield

    @property
    def output_folder(self):
        return self._output_folder

    @classmethod
    def reset(cls):
        super().reset()
        cls.metrics = {}

    def run(self):
        super().run()
        self.total = 1
        self.run_jobs()

    def parse_output(self):

        with open(os.path.join(self.log_folder, "bbtools_out.log"), "r") as f:
            lines = f.readlines()

        for line in lines:
            if line.startswith("Main genome scaffold total:"):
                self.metrics["Number of scaffolds"] = line.split(":")[-1].strip()
            elif line.startswith("Main genome contig total:"):
                self.metrics["Number of contigs"] = line.split(":")[-1].strip()
            elif line.startswith("Main genome scaffold sequence total:"):
                self.metrics["Total length"] = line.split(":")[-1].strip()
            elif line.startswith("Main genome contig sequence total:"):
                self.metrics["Percent gaps"] = (
                    line.split("\t")[-1].strip().strip(" gap")
                )
            elif line.startswith("Main genome scaffold N/L50:"):
                nl50 = line.split(":")[-1].strip().split("/")
                if int(nl50[0]) < int(
                    nl50[1]
                ):  # The N50/L50 values are inverted. Add a condition so if this is
                    # fixed in bbtools in future versions, it will still work.
                    self.metrics["Scaffold N50"] = nl50[1].strip()
                else:
                    self.metrics["Scaffold N50"] = nl50[0].strip()
            elif line.startswith("Main genome contig N/L50:"):
                nl50 = line.split(":")[-1].strip().split("/")
                if int(nl50[0]) < int(nl50[1]):
                    self.metrics["Contigs N50"] = nl50[1].strip()
                else:
                    self.metrics["Contigs N50"] = nl50[0].strip()
