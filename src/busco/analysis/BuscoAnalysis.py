#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: BuscoAnalysis
   :synopsis: BuscoAnalysis implements general BUSCO analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 5.4.0

Copyright (c) 2016-2022, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.
"""

from abc import ABCMeta, abstractmethod
from busco.BuscoConfig import BaseConfig, BuscoConfigAuto
from busco.busco_tools.hmmer import HMMERRunner
import os
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.Exceptions import BatchFatalError, BuscoError

logger = BuscoLogger.get_logger(__name__)


class BuscoAnalysis(metaclass=ABCMeta):
    """
    This abstract base class (ABC) defines methods required for most of BUSCO analyses and has to be extended
    by each specific analysis class
    """

    config = None

    def __init__(self):
        """
        1) load parameters
        2) load and validate tools
        3) check data and dataset integrity
        4) Ready for analysis
        """
        super().__init__()

        # Get paths
        self._lineage_results_dir = self.config.get("busco_run", "lineage_results_dir")
        self.main_out = self.config.get("busco_run", "main_out")
        self._working_dir = (
            os.path.join(self.main_out, "auto_lineage")
            if isinstance(self.config, BuscoConfigAuto)
            else self.main_out
        )
        self.run_folder = os.path.join(self._working_dir, self._lineage_results_dir)
        self._log_folder = os.path.join(self.main_out, "logs")

        # Get other useful variables
        self.input_file = self.config.get("busco_run", "in")
        self.lineage_dataset = self.config.get("busco_run", "lineage_dataset")
        self.lineage_name = os.path.basename(self.lineage_dataset)
        self.domain = self.config.get("busco_run", "domain")
        self._has_variants_file = os.path.exists(
            os.path.join(self.lineage_dataset, "ancestral_variants")
        )
        self.has_dataset_config_file = os.path.exists(
            os.path.join(self.lineage_dataset, "dataset.cfg")
        )
        self._dataset_creation_date = self.config.get("busco_run", "creation_date")
        self.restart = self.config.getboolean("busco_run", "restart")

        self.gene_details = (
            {}
        )  # Dictionary containing coordinate information for predicted genes.
        self.headers = set()

        self._lineages_download_path = os.path.join(
            self.config.get("busco_run", "download_path"), "lineages"
        )

        self.hmmer_runner = None

    @abstractmethod
    def cleanup(self):
        # Delete any non-decompressed files in busco_downloads
        try:
            for dataset_name in os.listdir(self._lineages_download_path):
                if dataset_name.endswith((".gz", ".tar")):
                    os.remove(dataset_name)
        except OSError:
            pass

    @abstractmethod
    @log(
        "Running BUSCO using lineage dataset {0} ({1}, {2})",
        logger,
        attr_name=["lineage_name", "domain", "_dataset_creation_date"],
        on_func_exit=True,
    )
    def run_analysis(self):
        """
        Abstract method, override to call all needed steps for running the child analysis.
        """
        self._create_dirs()
        self.init_tools()
        self._check_data_integrity()

    @log("***** Run HMMER on gene sequences *****", logger)
    def run_hmmer(self, input_sequences, busco_ids=None):
        """
        This function runs hmmsearch.
        """
        if not busco_ids:
            files = sorted(os.listdir(os.path.join(self.lineage_dataset, "hmms")))
            busco_ids = [
                os.path.splitext(f)[0] for f in files
            ]  # Each Busco ID has a HMM file of the form "<busco_id>.hmm"
        self.hmmer_runner.configure_runner(
            input_sequences, busco_ids, self._mode, self.gene_details
        )
        if self.restart and self.hmmer_runner.check_previous_completed_run():
            logger.info("Skipping HMMER run as output already processed")
        elif self.restart and os.path.exists(
            "{}.tar.gz".format(self.hmmer_runner.output_folder)
        ):
            raise BuscoError(
                "Restart mode incompatible with a previously compressed (--tar) run. Please decompress the HMMER "
                "results folder and try again."
            )
        elif len(os.listdir(self.hmmer_runner.results_dir)) > 0:
            raise BuscoError(
                "HMMER results directory not empty. If you are running in restart mode, make sure you are "
                "using the same eukaryotic gene predictor (metaeuk/augustus) as before."
            )
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.hmmer_runner.run()
        self.hmmer_runner.process_output()
        self.validate_output()
        self.hmmer_runner.filter()
        self.hmmer_runner.consolidate_busco_lists()
        output = self.hmmer_runner.create_output_content()
        self.hmmer_runner.write_hmmer_results(output)
        self.hmmer_runner.record_results()
        self.hmmer_runner.produce_hmmer_summary()
        return

    def validate_output(
        self,
    ):  # Transparent method that can be overwritten by child classes
        return

    @log("Checking dataset for HMM profiles", logger, debug=True)
    def _check_dataset_integrity(self):
        """
        Check the input dataset for hmm profiles, both files and folder are available
        Note: score and length cutoffs are checked when read by hmmer_runner: see _load_scores and _load_lengths
        Note: dataset.cfg file is not mandatory for offline mode

        :raises BuscoError: if the dataset is missing files or folders
        """

        # Check hmm files exist
        files = os.listdir(os.path.join(self.lineage_dataset, "hmms"))
        if not files:
            raise BuscoError(
                "The dataset you provided lacks hmm profiles in {}".format(
                    os.path.join(self.lineage_dataset, "hmms")
                )
            )

        if self.domain == "eukaryota":
            # Check prfl folder exists and contains profiles
            for dirpath, dirnames, files in os.walk(
                os.path.join(self.lineage_dataset, "prfl")
            ):
                if not files:
                    raise BuscoError(
                        "The dataset you provided lacks elements in {}".format(
                            os.path.join(self.lineage_dataset, "prfl")
                        )
                    )

        if not self._has_variants_file:
            logger.warning(
                "The dataset you provided does not contain the file ancestral_variants, likely because it "
                'is an old version. All blast steps will use the file "ancestral" instead'
            )
        if not self.has_dataset_config_file:
            raise BuscoError(
                "The dataset you provided is missing the dataset.cfg file and is therefore corrupted."
            )

        return

    def _check_data_integrity(self):
        self._check_dataset_integrity()
        if not os.stat(self.input_file).st_size > 0:
            raise BuscoError("Input file is empty.")
        try:
            with open(self.input_file) as f:
                for line in f:
                    if line.startswith(">"):
                        self._check_fasta_header(line)
                        self._check_seq_uniqueness(line)
        except UnicodeDecodeError as ude:
            raise BuscoError(ude.__str__())
        return

    def _check_seq_uniqueness(self, line):
        seq_id = line.split(" ")[0]
        if seq_id in self.headers:
            raise BuscoError("Duplicate of sequence {} in input file".format(seq_id))
        self.headers.add(seq_id)
        return

    @staticmethod
    def _check_fasta_header(header):
        """
        This function checks problematic characters in fasta headers,
        and warns the user and stops the execution
        :param header: a fasta header to check
        :type header: str
        :raises BuscoError: if a problematic character is found
        """
        for char in BaseConfig.FORBIDDEN_HEADER_CHARS:
            if char in header:
                raise BuscoError(
                    'The character "%s" is present in the fasta header %s, '
                    "which will crash BUSCO. Please clean the header of your "
                    "input file." % (char, header.strip())
                )

        for char in BaseConfig.FORBIDDEN_HEADER_CHARS_BEFORE_SPLIT:
            if char in header.split()[0]:
                raise BuscoError(
                    'The character "%s" is present in the fasta header %s, '
                    "which will crash Reader. Please clean the header of your"
                    " input file." % (char, header.split()[0].strip())
                )

        if header.split()[0] == ">":
            raise BuscoError(
                "A space is present in the fasta header %s, directly after "
                '">" which will crash Reader. Please clean the header of '
                "your input file." % (header.strip())
            )

    def _create_dirs(self):
        """
        Create the run (main) directory, log directory and the temporary directories
        :return:
        """
        self._create_main_dir()
        self._create_log_dir()
        # self._create_tmp_dir()

    def _create_log_dir(self):
        """
        Create a subfolder of the main output folder that contains all log files from BUSCO and the external tools used.
        :return:
        """
        if not os.path.exists(self._log_folder):
            os.mkdir(self._log_folder)
        return

    def _create_main_dir(self):
        """
        This function creates the run (main) directory
        :raises BatchFatalError: if write permissions are not available to the specified location
        """
        try:
            os.makedirs(self.run_folder)
        except FileExistsError:
            if not self.restart:
                raise BatchFatalError(
                    "Something went wrong. BUSCO stopped before overwriting run folder "
                    "{}".format(self.run_folder)
                )
        except PermissionError:
            raise BatchFatalError(
                "Cannot write to the output directory, please make sure "
                "you have write permissions to {}".format(self.run_folder)
            )
        return

    @log("Check all required tools are accessible...", logger, debug=True)
    def init_tools(self):
        """
        Init the tools needed for the analysis. HMMER is needed for all BUSCO analysis types.
        """
        self.hmmer_runner = HMMERRunner()
        return

    def reset(self):
        if (
            self.hmmer_runner
        ):  # If final run has already been run, then the hmmer_runner object in the final runner object will still
            # be set to None
            self.hmmer_runner.reset()

    @property
    @abstractmethod
    def _mode(self):
        pass
