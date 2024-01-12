# coding: utf-8
"""
base.py

A collection of classes and methods used by all tools.

Author(s): Matthew Berkeley, Mathieu Seppey, Mose Manni, Felipe Simao, Rob Waterhouse

Copyright (c) 2015-2024, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

import os
from busco.BuscoLogger import BuscoLogger
from busco.busco_tools.Toolset import Tool
from shutil import which
from abc import ABCMeta, abstractmethod
from busco.BuscoConfig import BuscoConfigAuto
from busco.Exceptions import BuscoError
import time
import gzip
import shutil
from collections import defaultdict
import numpy as np

logger = BuscoLogger.get_logger(__name__)


class ToolException(Exception):
    """
    Module-specific exception
    """

    def __init__(self, value=""):
        self.value = value

    def __str__(self):
        return self.value


class BaseRunner(Tool, metaclass=ABCMeta):

    config = None
    tool_versions = {}

    def __init__(self):
        super().__init__()
        self.run_number = 0
        self.input_file = self.config.get("busco_run", "in")
        # type(self).summary["versions"]["busco"] = busco.__version__
        self.main_out = self.config.get("busco_run", "main_out")
        self.working_dir = (
            os.path.join(self.main_out, "auto_lineage")
            if isinstance(self.config, BuscoConfigAuto)
            else self.main_out
        )
        self.lineage_results_dir = self.config.get("busco_run", "lineage_results_dir")
        self.run_folder = os.path.join(self.working_dir, self.lineage_results_dir)
        self.log_folder = os.path.join(self.main_out, "logs")
        self.cpus = self.config.getint("busco_run", "cpu")
        self.lineage_dataset = self.config.get("busco_run", "lineage_dataset")
        self.domain = self.config.get("busco_run", "domain")

        try:
            self.check_tool_available()
        except ToolException:
            raise ToolException(
                "{} tool cannot be found. Please check the 'path' and 'command' parameters "
                "provided in the config file or make sure the tool is available in your working environment.".format(
                    self.name
                )
            )
        self.version = self.get_version()
        type(self).tool_versions[self.name] = self.version
        self.check_tool_dependencies()

        self.checkpoint_file = None

        self.logfile_path_out = os.path.join(
            self.config.get("busco_run", "main_out"),
            "logs",
            "{}_out.log".format(self.name),
        )
        self.logfile_path_err = (
            self.logfile_path_out.rpartition("_out.log")[0] + "_err.log"
        )
        self.add_args = {}

    def init_checkpoint_file(self):
        self.checkpoint_file = os.path.join(self.output_folder, ".checkpoint")

    def write_checkpoint_file(self, additional_args=[]):
        with open(self.checkpoint_file, "a") as cpt_file:
            cpt_file.write("Tool: {}\n".format(self.name))
            cpt_file.write("Version: {}\n".format(self.version))
            cpt_file.write("Run: {}\n".format(self.run_number))
            for args in additional_args:
                cpt_file.write("{}: {}\n".format(args[0], args[1]))
            current_time = time.strftime("%d/%m/%Y %H:%M:%S")
            cpt_file.write("Time: {}\n".format(current_time))
            self.config.run_stats["{}_start_time".format(self.name)] = current_time
            cpt_file.write("Completed {} jobs\n\n".format(self.total))

    def check_previous_completed_run(self, additional_args=[]):
        if not os.path.exists(self.checkpoint_file):
            return False
        else:
            with open(self.checkpoint_file, "r") as cpt_file:
                block_size = 6 + len(additional_args)
                lines = cpt_file.readlines()
                tool_names = [s.strip().split(": ")[1] for s in lines[0::block_size]]
                tool_versions = [s.strip().split(": ")[1] for s in lines[1::block_size]]
                tool_run_numbers = [
                    s.strip().split(": ")[1] for s in lines[2::block_size]
                ]
                self.add_args = {}
                for a, arg in enumerate(additional_args):
                    self.add_args[arg] = [
                        s.strip().split(": ")[1]
                        for s in lines[a + 3 :: block_size]
                        if s.strip().split(": ")[0] == arg
                    ]
                try:
                    start_search = 0
                    while True:
                        tool_ind = tool_names.index(self.name, start_search)
                        if int(tool_run_numbers[tool_ind]) == int(self.run_number):
                            if str(self.version) != str(tool_versions[tool_ind]):
                                logger.warning(
                                    "A previous run used {} version {}. "
                                    "The restarted run is using {} version "
                                    "{}".format(
                                        self.name,
                                        tool_versions[tool_ind],
                                        self.name,
                                        self.version,
                                    )
                                )
                            return True
                        elif int(tool_run_numbers[tool_ind]) < int(self.run_number):
                            start_search = tool_ind + 1
                        else:
                            raise BuscoError(
                                "Something went wrong. Information for {} run {} missing but "
                                "information for run {} found.".format(
                                    self.name,
                                    self.run_number,
                                    tool_run_numbers[tool_ind],
                                )
                            )

                except ValueError:
                    return False

                except TypeError:
                    logger.warning(
                        "Unable to parse {} file. Restart mode not available.".format(
                            self.checkpoint_file
                        )
                    )

    @abstractmethod
    def check_tool_dependencies(self):
        pass

    @abstractmethod
    def configure_job(self, *args):
        pass

    @abstractmethod
    def configure_runner(self, *args):
        self.init_checkpoint_file()

    @abstractmethod
    def generate_job_args(self):
        pass

    @property
    @abstractmethod
    def output_folder(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def name(self):
        raise NotImplementedError

    @abstractmethod
    def run(self):
        if self.version is not None:
            logger.debug("Tool: {}".format(self.name))
            logger.debug("Version: {}".format(self.version))

    @staticmethod
    def create_dirs(dirnames, overwrite=False):
        """
        Create all required directories

        :param dirnames: list of paths already constructed
        :return:
        """
        if isinstance(dirnames, str):
            dirnames = [dirnames]
        elif isinstance(dirnames, list):
            pass
        else:
            raise TypeError("'dirnames' should be either a str or a list")

        for d in dirnames:
            if overwrite:
                if os.path.exists(d):
                    try:
                        shutil.rmtree(d)
                    except OSError:
                        logger.warning("Unable to remove {}".format(d))
                        pass
            os.makedirs(d, exist_ok=True)

    def check_tool_available(self):
        """
        Check tool's availability.


        :return: True if the tool can be run, False if it is not the case
        :rtype: bool
        """
        try:
            self.get_tool_from_config()
        except ToolException:
            try:
                self.get_tool_from_environment()
            except ToolException:
                raise

        return which(self.cmd) is not None  # True if tool available

    def get_tool_from_environment(self):
        which_tool = which(self.cmd)
        if not which_tool:
            raise ToolException()

    def get_tool_from_config(self):
        """
        1. The section ['name'] is available in the config
        2. This section contains keys 'path' and 'command'
        3. The string resulted from concatenation of values of these two keys
        represents the full path to the command
        :return:
        """
        if not self.config.has_section(self.name):
            raise ToolException()

        if not self.config.has_option(self.name, "path") or not self.config.get(
            self.name, "path"
        ):
            raise ToolException()

        if self.config.has_option(self.name, "command"):
            executable = self.config.get(self.name, "command")
        else:
            executable = self.name

        self.cmd = os.path.join(self.config.get(self.name, "path"), executable)

        return

    @abstractmethod
    def get_version(self):
        return

    @classmethod
    def reset(cls):
        BaseRunner.config = None
        BaseRunner.tool_versions = {}


class GenePredictor(BaseRunner, metaclass=ABCMeta):
    def __init__(self):
        super().__init__()
        self.gene_details = defaultdict(dict)
        self.sequences_aa = {}
        self.sequences_nt = {}

    @abstractmethod
    def record_gene_details(self):
        """
        Record all required match details from gene predictor output. The exact contents of the dictionary will vary
        depending on the gene predictor and pipeline.
        :return:
        """
        return

    @staticmethod
    def get_matches(results_seq):
        results_sorted = np.sort(
            results_seq, order=["low_coord"]
        )  # sort to facilitate a single-pass coordinate check
        for row1 in results_sorted:
            low_coord = row1["low_coord"]
            high_coord = row1["high_coord"]
            current_seqid = row1["gene_id"]
            matches = results_sorted[
                (results_sorted["low_coord"] >= low_coord)
                & (results_sorted["low_coord"] < high_coord)
            ]  # find entries with a start coordinate between the current exon start and end
            yield row1, current_seqid, matches

    def detect_overlap(self, record_arr, seq):
        results_seq = record_arr[record_arr["contig_id"] == seq]
        overlaps = []
        entries_to_remove = []
        handled_gene_ids = set()
        match_finder = self.get_matches(results_seq)
        for match in match_finder:
            row1, current_seqid, matches = match
            exon_id1 = "{}|{}-{}".format(
                row1["gene_id"], row1["low_coord"], row1["high_coord"]
            )
            if exon_id1 in entries_to_remove:
                continue
            for row2 in matches:
                exon_id2 = "{}|{}-{}".format(
                    row2["gene_id"], row2["low_coord"], row2["high_coord"]
                )
                if exon_id2 in handled_gene_ids:
                    continue
                if exon_id1 == exon_id2:  # don't consider overlaps with self
                    continue
                elif (
                    (
                        (
                            row1["gene_id"].split("|", 1)[-1]
                            != row2["gene_id"].split("|", 1)[-1]
                        )  # for efficiency skip overlaps that are
                        # actually the same gene matching multiple BUSCOs, as this will be dealt with in the
                        # filtering step later
                    )
                    and (
                        row1["strand"]
                        == row2["strand"]  # check overlaps are on the same strand
                    )
                    and (
                        row1["low_coord"] % 3
                        == row2["low_coord"]
                        % 3  # check overlaps are in the same reading frame
                    )
                ):
                    if row1["score"] > row2["score"]:
                        entries_to_remove.append(exon_id2)
                        overlaps.append((row2, row1))  # first entry has the lower score
                    else:
                        entries_to_remove.append(exon_id1)
                        overlaps.append((row1, row2))
            handled_gene_ids.add(exon_id1)
        return overlaps

    def find_overlaps(self, structured_arr):
        seq_ids = np.unique(structured_arr["contig_id"])
        overlaps = []
        for seq_id in seq_ids:
            overlaps += self.detect_overlap(structured_arr, seq_id)

        max_bitscores = []
        for overlap in overlaps:
            max_bitscores.append(max([x["score"] for x in overlap]))
        sort_order = np.argsort(max_bitscores)[::-1]
        overlaps = [overlaps[i] for i in sort_order]

        return overlaps

    @staticmethod
    def decompress_refseq_file(
        gzip_file,
    ):
        unzipped_filename = gzip_file.split(".gz")[0]
        if not os.path.exists(unzipped_filename):
            with gzip.open(gzip_file, "rb") as compressed_file:
                with open(unzipped_filename, "wb") as decompressed_file:
                    for line in compressed_file:
                        decompressed_file.write(line)
        if os.path.exists(gzip_file):
            try:
                os.remove(gzip_file)
            except OSError:
                logger.warning(
                    "Unable to remove compressed refseq file in dataset download"
                )
                pass
        return unzipped_filename


class NoGenesError(Exception):
    def __init__(self, gene_predictor):
        self.gene_predictor = gene_predictor


class NoRerunFile(Exception):
    def __init__(self):
        pass
