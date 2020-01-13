#!/usr/bin/env python
# coding: utf-8
"""
.. module:: BuscoAnalysis
   :synopsis: BuscoAnalysis implements general BUSCO analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 3.0.1

Copyright (c) 2016-2017, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""

from abc import ABCMeta, abstractmethod
import busco
from busco.BuscoConfig import BuscoConfig, BuscoConfigMain, BuscoConfigAuto
from busco.BuscoTools import HMMERRunner
import inspect
import os
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.Toolset import Tool
import subprocess
from Bio import SeqIO

logger = BuscoLogger.get_logger(__name__)


class BuscoAnalysis(metaclass=ABCMeta):
    """
    This abstract base class (ABC) defines methods required for most of BUSCO analyses and has to be extended
    by each specific analysis class
    """


    def __init__(self, config):
        """
        1) load parameters
        2) load and validate tools
        3) check data and dataset integrity
        4) Ready for analysis

        :param config: Values of all parameters to be used during the analysis
        :type config: BuscoConfig
        """
        self._config = config

        # Get paths
        self.lineage_results_dir = self._config.get("busco_run", "lineage_results_dir")
        self.main_out = self._config.get("busco_run", "main_out")  # todo: decide which are hidden attributes
        self.working_dir = (os.path.join(self.main_out, "auto_lineage")
                            if isinstance(self._config, BuscoConfigAuto)
                            else self.main_out)
        self.run_folder = os.path.join(self.working_dir, self.lineage_results_dir)
        self.log_folder = os.path.join(self.main_out, "logs")
        self._input_file = self._config.get("busco_run", "in")
        self._lineage_dataset = self._config.get("busco_run", "lineage_dataset")
        self._lineage_name = os.path.basename(self._lineage_dataset)
        self._datasets_version = self._config.get("busco_run", "datasets_version")
        super().__init__()

        # Get other useful variables
        self._cpus = self._config.getint("busco_run", "cpu")
        self._domain = self._config.get("busco_run", "domain")
        self._has_variants_file = os.path.exists(os.path.join(self._lineage_dataset, "ancestral_variants"))
        self._dataset_creation_date = self._config.get("busco_run", "creation_date")
        self._dataset_nb_species = self._config.get("busco_run", "number_of_species")
        self._dataset_nb_buscos = self._config.get("busco_run", "number_of_BUSCOs")

        # Get Busco downloader
        self.downloader = self._config.downloader

        # Create optimized command line call for the given input
        self.busco_type = "main" if isinstance(self._config, BuscoConfigMain) else "auto"
        # if self.busco_type == "main":
        #     self.set_rerun_busco_command(self._config.clargs)  # todo: rework rerun command

        # Variables storing BUSCO results
        self._missing_busco_list = []
        self._fragmented_busco_list = []
        self._gene_details = None  # Dictionary containing coordinate information for predicted genes.
        self.s_percent = None
        self.d_percent = None
        self.f_percent = None
        self._log_count = 0  # Dummy variable used to skip logging for intermediate eukaryote pipeline results.

    # TODO: catch unicode encoding exception and report invalid character line instead of doing content validation
    # todo: check config file exists before parsing

    @abstractmethod
    def _cleanup(self):
        # Delete any non-decompressed files in busco_downloads
        try:
            for dataset_name in os.listdir(os.path.join(self._config.get("busco_run", "download_path"), "lineages")):
                if dataset_name.endswith((".gz", ".tar")):
                    os.remove(dataset_name)
        except OSError:
            pass

    def _check_data_integrity(self):
        self._check_dataset_integrity()
        if not os.stat(self._input_file).st_size > 0:
            raise SystemExit("Input file is empty.")
        with open(self._input_file) as f:
            for line in f:
                if line.startswith(">"):
                    self._check_fasta_header(line)
        return

    def get_checkpoint(self): # TODO: rework checkpoint system
        """
        This function return the checkpoint if the checkpoint.tmp file exits or None if absent
        :return: the checkpoint name
        :rtype: int
        """
        checkpt_name = None
        checkpoint_file = os.path.join(self.run_folder, "checkpoint.tmp")
        if os.path.exists(checkpoint_file):
            with open(checkpoint_file, "r") as check_file:
                line = check_file.readline()
                self._random = line.split(".")[-1] # Reset random suffix
            checkpt_name = int(line.split(".")[0])
        return checkpt_name

    @abstractmethod
    @log("Running BUSCO using lineage dataset {0} ({1}, {2})", logger, attr_name=["_lineage_name", "_domain", "_dataset_creation_date"], on_func_exit=True)
    def run_analysis(self):
        """
        Abstract method, override to call all needed steps for running the child analysis.
        """
        self.create_dirs()
        self.init_tools()
        self.check_tool_dependencies()
        self._check_data_integrity()

    @log("Checking dataset for HMM profiles", logger, debug=True)
    def _check_dataset_integrity(self):
        """
        Check the input dataset for hmm profiles, both files and folder are available
        Note: score and length cutoffs are checked when read,
        See _load_scores and _load_lengths

        :raises SystemExit: if the dataset is missing files or folders
        """

        # Check hmm files exist
        files = os.listdir(os.path.join(self._lineage_dataset, "hmms"))
        if not files:
            raise SystemExit("The dataset you provided lacks hmm profiles in {}".format(
                os.path.join(self._lineage_dataset, "hmms")))

        if self._domain == "eukaryota":
            # Check prfl folder exists and contains profiles
            for dirpath, dirnames, files in os.walk(os.path.join(self._lineage_dataset, "prfl")):
                if not files:
                    raise SystemExit("The dataset you provided lacks elements in {}".format(
                        os.path.join(self._lineage_dataset, "prfl")))

        # note: score and length cutoffs are checked when read,
        # see _load_scores and _load_lengths
        # ancestral would cause blast to fail, and be detected, see _blast() # TODO: clarify comment
        # dataset.cfg is not mandatory

        if not self._has_variants_file:
            logger.warning("The dataset you provided does not contain the file ancestral_variants, likely because it "
                           "is an old version. All blast steps will use the file \"ancestral\" instead")

        return

    def _check_fasta_header(self, header):
        """
        This function checks problematic characters in fasta headers,
        and warns the user and stops the execution
        :param header: a fasta header to check
        :type header: str
        :raises SystemExit: if a problematic character is found
        """
        for char in BuscoConfig.FORBIDDEN_HEADER_CHARS:
            if char in header:
                raise SystemExit(
                    "The character \"%s\" is present in the fasta header %s, "
                    "which will crash BUSCO. Please clean the header of your "
                    "input file." % (char, header.strip()))


        for char in BuscoConfig.FORBIDDEN_HEADER_CHARS_BEFORE_SPLIT:
            if char in header.split()[0]:
                raise SystemExit(
                    "The character \"%s\" is present in the fasta header %s, "
                    "which will crash Reader. Please clean the header of your"
                    " input file." % (char, header.split()[0].strip()))

        if header.split()[0] == ">":
            raise SystemExit(
                "A space is present in the fasta header %s, directly after "
                "\">\" which will crash Reader. Please clean the header of "
                "your input file." % (header.strip()))

    def check_tool_dependencies(self):
        """
        check dependencies on tools
        :raises SystemExit: if a Tool version is not supported
        """
        # check hmm version
        if not self._get_hmmer_version() >= BuscoConfig.HMMER_VERSION:
            raise SystemExit(
                "HMMer version detected is not supported, please use HMMer v.{} +".format(BuscoConfig.HMMER_VERSION))
        return

    @abstractmethod
    def create_dirs(self):
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
        if not os.path.exists(self.log_folder):
            os.mkdir(self.log_folder)
        return

    def _create_main_dir(self):
        """
        This function creates the run (main) directory
        :raises SystemExit: if write permissions are not available to the specified location
        """
        try:
            os.makedirs(self.run_folder)
        except FileExistsError:
            raise SystemExit("Something went wrong. BUSCO stopped before overwriting run folder {}".format(self.run_folder))
        except PermissionError:
            raise SystemExit(
                "Cannot write to the output directory, please make sure "
                "you have write permissions to {}".format(self.run_folder))
        return

    # @log("Temp directory is {}", logger, attr_name="_tmp", on_func_exit=True)
    # def _create_tmp_dir(self):
    #     """
    #     This function creates the tmp directory
    #     :raises
    #     SystemExit: if the user cannot write in the tmp directory
    #     """
    #     try:
    #         if not os.path.exists(self._tmp):
    #             os.makedirs(self._tmp)
    #
    #     except OSError:
    #         raise SystemExit(
    #             "Cannot write to the temp directory, please make sure "
    #             "you have write permissions to {}".format(self._tmp))
    #     return



    def _get_busco_percentages(self):
        self.single_copy = len(self.hmmer_runner.single_copy_buscos)  # int
        self.multi_copy = len(self.hmmer_runner.multi_copy_buscos)  # int
        self.only_fragments = len(self.hmmer_runner.fragmented_buscos)  # int
        self.total_buscos = len(self.hmmer_runner.cutoff_dict)

        # Get percentage of each kind of BUSCO match
        self.s_percent = abs(round((self.single_copy / self.total_buscos) * 100, 1))
        self.d_percent = abs(round((self.multi_copy / self.total_buscos) * 100, 1))
        self.f_percent = abs(round((self.only_fragments / self.total_buscos) * 100, 1))

        return self.single_copy, self.multi_copy, self.only_fragments, self.total_buscos

    def _get_hmmer_version(self):
        """
        check the Tool has the correct version
        :raises SystemExit: if the version is not correct
        """
        hmmer_version = subprocess.check_output([self._hmmer_tool.cmd, "-h"], shell=False)
        hmmer_version = hmmer_version.decode("utf-8")
        try:
            hmmer_version = hmmer_version.split("\n")[1].split()[2]
            hmmer_version = float(hmmer_version[:3])
        except ValueError:
            # to avoid a crash with a super old version
            hmmer_version = hmmer_version.split("\n")[1].split()[1]
            hmmer_version = float(hmmer_version[:3])
        finally:
            return hmmer_version

    @log("Check all required tools are accessible...", logger, debug=True)
    def init_tools(self):
        """
        Init the tools needed for the analysis. HMMER is needed for all BUSCO analysis types.
        """
        try:
            assert(isinstance(self._hmmer_tool, Tool))
        except AttributeError:
            self._hmmer_tool = Tool("hmmsearch", self._config)
        except AssertionError:
            raise SystemExit("HMMer should be a tool")

        return

    @property
    @abstractmethod
    def _mode(self):
        pass

    # @log("This is not an incomplete run that can be restarted", logger, iswarn=True)
    # # Todo: decide if mini functions are necessary to facilitate decorator logging
    # def _not_incomplete_run(self):
    #     self._restart = False

    def _produce_hmmer_summary(self):
        single_copy, multi_copy, only_fragments, total_buscos = self._get_busco_percentages()

        self.hmmer_results_lines = []
        self.hmmer_results_lines.append("***** Results: *****\n\n")
        self.one_line_summary = "C:{}%[S:{}%,D:{}%],F:{}%,M:{}%,n:{}\t{}\n".format(
            round(self.s_percent + self.d_percent, 1), self.s_percent, self.d_percent, self.f_percent,
            round(100 - self.s_percent - self.d_percent - self.f_percent, 1), total_buscos, "   ")
        self.hmmer_results_lines.append(self.one_line_summary)
        self.hmmer_results_lines.append("{}\tComplete BUSCOs (C)\t\t\t{}\n".format(single_copy + multi_copy, "   "))
        self.hmmer_results_lines.append("{}\tComplete and single-copy BUSCOs (S)\t{}\n".format(single_copy, "   "))
        self.hmmer_results_lines.append("{}\tComplete and duplicated BUSCOs (D)\t{}\n".format(multi_copy, "   "))
        self.hmmer_results_lines.append("{}\tFragmented BUSCOs (F)\t\t\t{}\n".format(only_fragments, "   "))
        self.hmmer_results_lines.append("{}\tMissing BUSCOs (M)\t\t\t{}\n".format(
            total_buscos - single_copy - multi_copy - only_fragments, "   "))
        self.hmmer_results_lines.append("{}\tTotal BUSCO groups searched\t\t{}\n".format(total_buscos, "   "))

        with open(os.path.join(self.run_folder, "short_summary.txt"), "w") as summary_file:

            self._write_output_header(summary_file, no_table_header=True)
            summary_file.write("# Summarized benchmarking in BUSCO notation for file {}\n"
                               "# BUSCO was run in mode: {}\n\n".format(self._input_file, self._mode))

            for line in self.hmmer_results_lines:
                summary_file.write("\t{}".format(line))

        if isinstance(self._config, BuscoConfigAuto):  # todo: rework this if/else block
            self._one_line_hmmer_summary()
        elif self._domain == "eukaryota" and self._log_count == 0:
            self._log_count += 1
            self._produce_full_hmmer_summary_debug()
        else:
            self._one_line_hmmer_summary()
        return

    @log("{}", logger, attr_name="hmmer_results_lines", apply="join", on_func_exit=True)
    def _produce_full_hmmer_summary(self):
        return

    @log("{}", logger, attr_name="hmmer_results_lines", apply="join", on_func_exit=True, debug=True)
    def _produce_full_hmmer_summary_debug(self):
        return

    @log("{}", logger, attr_name="one_line_summary", on_func_exit=True)
    def _one_line_hmmer_summary(self):
        self.one_line_summary = "Results:\t{}".format(self.one_line_summary)
        return

    @log("***** Run HMMER on gene sequences *****", logger)
    def run_hmmer(self, input_sequences):
        """
        This function runs hmmsearch.
        """
        self._hmmer_tool.total = 0
        self._hmmer_tool.nb_done = 0
        hmmer_output_dir = os.path.join(self.run_folder, "hmmer_output")
        if not os.path.exists(hmmer_output_dir):
            os.makedirs(hmmer_output_dir)

        files = sorted(os.listdir(os.path.join(self._lineage_dataset, "hmms")))
        busco_ids = [os.path.splitext(f)[0] for f in files]  # Each Busco ID has a HMM file of the form "<busco_id>.hmm"

        self.hmmer_runner = HMMERRunner(self._hmmer_tool, input_sequences, busco_ids, hmmer_output_dir,
                            self._lineage_dataset, self._mode, self._cpus, self._gene_details, self._datasets_version)
        self.hmmer_runner.load_buscos()
        self.hmmer_runner.run()
        self.hmmer_runner.process_output()
        self._write_hmmer_results()
        self._produce_hmmer_summary()
        return

    def _write_buscos_to_file(self, sequences_aa, sequences_nt=None):
        """
        Write BUSCO matching sequences to output fasta files. Each sequence is printed in a separate file and both
        nucleotide and amino acid versions are created.
        :param busco_type: one of ["single_copy", "multi_copy", "fragmented"]
        :return:
        """
        for busco_type in ["single_copy", "multi_copy", "fragmented"]:
            if busco_type == "single_copy":
                output_dir = os.path.join(self.run_folder, "busco_sequences", "single_copy_busco_sequences")
                busco_matches = self.hmmer_runner.single_copy_buscos
            elif busco_type == "multi_copy":
                output_dir = os.path.join(self.run_folder, "busco_sequences", "multi_copy_busco_sequences")
                busco_matches = self.hmmer_runner.multi_copy_buscos
            elif busco_type == "fragmented":
                output_dir = os.path.join(self.run_folder, "busco_sequences", "fragmented_busco_sequences")
                busco_matches = self.hmmer_runner.fragmented_buscos

            if not os.path.exists(output_dir):  # todo: move all create_dir commands to one place
                os.makedirs(output_dir)

            for busco, gene_matches in busco_matches.items():
                try:
                    aa_seqs, nt_seqs = zip(*[(sequences_aa[gene_id], sequences_nt[gene_id]) for gene_id in gene_matches])
                    with open(os.path.join(output_dir, "{}.fna".format(busco)), "w") as f2:
                        SeqIO.write(nt_seqs, f2, "fasta")
                except TypeError:
                    aa_seqs = [sequences_aa[gene_id] for gene_id in gene_matches]
                with open(os.path.join(output_dir, "{}.faa".format(busco)), "w") as f1:
                    SeqIO.write(aa_seqs, f1, "fasta")

        return

    # def _run_tarzip_hmmer_output(self):  # todo: rewrite using tarfile
    #     """
    #     This function tarzips "hmmer_output" results folder
    #     """
    #     self._p_open(["tar", "-C", "%s" % self.run_folder, "-zcf", "%shmmer_output.tar.gz" % self.run_folder,
    #                   "hmmer_output", "--remove-files"], "bash", shell=False)



    @log("To reproduce this run: {}", logger, attr_name="_rerun_cmd", on_func_exit=True)
    def set_rerun_busco_command(self, clargs):  # todo: reconfigure
        """
        This function sets the command line to call to reproduce this run
        """

        # Find python script path
        entry_point = ""
        frame_ind = -1
        while "run_BUSCO.py" not in entry_point:
            entry_point = inspect.stack()[frame_ind].filename
            frame_ind -= 1

        # Add required parameters and other options
        self._rerun_cmd = "python %s -i %s -o %s -l %s -m %s -c %s" % (entry_point, self._input_file, os.path.basename(self.main_out),
                                                                       self._lineage_dataset, self._mode, self._cpus)

        try:
            if self._long:
                self._rerun_cmd += " --long"
            if self._region_limit != BuscoConfig.DEFAULT_ARGS_VALUES["limit"]:
                self._rerun_cmd += " --limit %s" % self._region_limit
            # if self._tmp != BuscoConfig.DEFAULT_ARGS_VALUES["tmp_path"]:
            #     self._rerun_cmd += " -t %s" % self._tmp
            if self._ev_cutoff != BuscoConfig.DEFAULT_ARGS_VALUES["evalue"]:
                self._rerun_cmd += " -e %s" % self._ev_cutoff
            # if self._tarzip:
            #     self._rerun_cmd += " -z"
        except AttributeError:
            pass

        # Include any command line arguments issued by the user
        # arg_aliases = {"-i": "--in", "-o": "--out", "-l": "--lineage_dataset", "-m": "--mode", "-c": "--cpu",
        #                "-e": "--evalue", "-f": "--force", "-sp": "--species", "-z": "--tarzip",
        #                "-r": "--restart", "-q": "--quiet", "-v": "--version", "-h": "--help"}
        arg_aliases.update(dict(zip(arg_aliases.values(), arg_aliases.keys())))
        for a, arg in enumerate(clargs):
            if arg.startswith("-") and not arg in self._rerun_cmd:
                if arg in arg_aliases:
                    if arg_aliases[arg] in self._rerun_cmd:
                        continue
                if a + 1 < len(clargs) and not clargs[a + 1].startswith("-"):
                    self._rerun_cmd += " %s %s" % (arg, clargs[a + 1])
                else:
                    self._rerun_cmd += " %s" % arg
        return

    def _write_hmmer_results(self):
        """
        Create two output files: one with information on all BUSCOs for the given dataset and the other with a list of
        all BUSCOs that were not found.
        :return:
        """

        with open(os.path.join(self.run_folder, "full_table.tsv"), "w") as f_out:
            self._write_output_header(f_out)

            output_lines = self.hmmer_runner._create_output_content()

            with open(os.path.join(self.run_folder, "missing_busco_list.tsv"), "w") as miss_out:

                self._write_output_header(miss_out, missing_list=True)

                missing_buscos_lines, missing_buscos = self.hmmer_runner._list_missing_buscos()
                output_lines += missing_buscos_lines

                for missing_busco in sorted(missing_buscos):
                    miss_out.write("{}\n".format(missing_busco))

            sorted_output_lines = self._sort_lines(output_lines)
            for busco in sorted_output_lines:
                f_out.write(busco)
        return

    @staticmethod
    def _sort_lines(lines):
        sorted_lines = sorted(lines, key=lambda x: int(x.split("\t")[0].split("at")[0]))
        return sorted_lines




    def _write_output_header(self, file_object, missing_list=False, no_table_header=False):
        """
        Write a standardized file header containing information on the BUSCO run.
        :param file_object: Opened file object ready for writing
        :type file_object: file
        :return:
        """
        file_object.write("# BUSCO version is: {} \n"
                          "# The lineage dataset is: {} (Creation date: {}, number of species: {}, number of BUSCOs: {}"
                          ")\n".format(busco.__version__, self._lineage_name, self._dataset_creation_date,
                                     self._dataset_nb_species, self._dataset_nb_buscos))
        # if isinstance(self._config, BuscoConfigMain):  # todo: wait until rerun command properly implemented again
        #     file_object.write("# To reproduce this run: {}\n#\n".format(self._rerun_cmd))

        if no_table_header:
            pass
        elif missing_list:
            file_object.write("# Busco id\n")
        elif self._mode == "proteins" or self._mode == "transcriptome":
            if self.hmmer_runner.extra_columns:
                file_object.write("# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
            else:
                file_object.write("# Busco id\tStatus\tSequence\tScore\tLength\n")
        elif self._mode == "genome":
            if self.hmmer_runner.extra_columns:
                file_object.write("# Busco id\tStatus\tSequence\tGene Start\tGene End\tScore\tLength\tOrthoDB url\tDescription\n")
            else:
                file_object.write("# Busco id\tStatus\tSequence\tGene Start\tGene End\tScore\tLength\n")

        return

