import os
import re
from collections import defaultdict
import busco
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.Toolset import Tool
from busco.BuscoConfig import BuscoConfig, BuscoConfigMain
from Bio import SeqIO
import shutil
import csv
import numpy as np
from shutil import which
from abc import ABCMeta, abstractmethod
from configparser import NoOptionError
import subprocess
from busco.BuscoConfig import BuscoConfigAuto
import time

# todo: docstrings
logger = BuscoLogger.get_logger(__name__)


class ToolException(Exception):
    """
    Module-specific exception
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


class BaseRunner(Tool, metaclass=ABCMeta):

    config = None

    def __init__(self):
        super().__init__()
        self.run_number = 0
        self.input_file = self.config.get("busco_run", "in")
        self.main_out = self.config.get("busco_run", "main_out")
        self.working_dir = (os.path.join(self.main_out, "auto_lineage")
                            if isinstance(self.config, BuscoConfigAuto)
                            else self.main_out)
        self.lineage_results_dir = self.config.get("busco_run", "lineage_results_dir")
        self.run_folder = os.path.join(self.working_dir, self.lineage_results_dir)
        self.log_folder = os.path.join(self.main_out, "logs")
        self.cpus = self.config.getint("busco_run", "cpu")
        self.lineage_dataset = self.config.get("busco_run", "lineage_dataset")
        self.domain = self.config.get("busco_run", "domain")

        if not self.check_tool_available():
            raise ToolException("{} tool cannot be found. Please check the 'path' and 'command' parameters "
                                "provided in the config file. Do not include the command in the "
                                "path!".format(self.name))
        self.version = self.get_version()
        self.check_tool_dependencies()

        self.checkpoint_file = None

        self.logfile_path_out = os.path.join(self.config.get("busco_run", "main_out"), "logs",
                                             "{}_out.log".format(self.name))
        self.logfile_path_err = self.logfile_path_out.replace('_out.log', '_err.log')

    def init_checkpoint_file(self):
        self.checkpoint_file = os.path.join(self.output_folder, ".checkpoint")

    def write_checkpoint_file(self):
        with open(self.checkpoint_file, "a") as cpt_file:
            cpt_file.write("Tool: {}\n".format(self.name))
            cpt_file.write("Version: {}\n".format(self.version))
            cpt_file.write("Run: {}\n".format(self.run_number))
            cpt_file.write("Time: {}\n".format(time.strftime('%m/%d/%Y %H:%M:%S')))
            cpt_file.write("Completed {} jobs\n\n".format(self.total))

    def check_previous_completed_run(self):
        if not os.path.exists(self.checkpoint_file):
            return False
        else:
            with open(self.checkpoint_file, "r") as cpt_file:
                lines = cpt_file.readlines()
                tool_names = [s.strip().split(": ")[1] for s in lines[0::6]]
                tool_versions = [s.strip().split(": ")[1] for s in lines[1::6]]
                tool_run_numbers = [s.strip().split(": ")[1] for s in lines[2::6]]
                try:
                    start_search = 0
                    while True:
                        tool_ind = tool_names.index(self.name, start_search)
                        if str(self.version) != str(tool_versions[tool_ind]):
                            logger.warning("A previous run used {} version {}. "
                                           "The restarted run is using {} version "
                                           "{}".format(self.name, tool_versions[tool_ind], self.name, self.version))
                        if int(tool_run_numbers[tool_ind]) == int(self.run_number):
                            return True
                        elif int(tool_run_numbers[tool_ind]) < int(self.run_number):
                            start_search = tool_ind + 1
                        else:
                            raise SystemExit("Something went wrong. Information for {} run {} missing but "
                                             "information for run {} found.".format(self.name, self.run_number,
                                                                                    tool_run_numbers[tool_ind]))

                except ValueError:
                    return False

                except TypeError:
                    logger.warning("Unable to parse {} file. Restart mode not available.".format(self.checkpoint_file))

    @abstractmethod
    def check_tool_dependencies(self):
        pass

    @abstractmethod
    def configure_job(self, *args):
        pass

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
    def create_dirs(dirnames):
        """
        Create all required directories

        :param dirnames: list of paths already constructed
        :return:
        """
        if isinstance(dirnames, str):
            os.makedirs(dirnames, exist_ok=True)
        elif isinstance(dirnames, list):
            for d in dirnames:
                os.makedirs(d, exist_ok=True)
        else:
            raise TypeError("'dirnames' should be either a str or a list")

    def check_tool_available(self):
        """
        Check tool's availability.
        1. The section ['name'] is available in the config
        2. This section contains keys 'path' and 'command'
        3. The string resulted from concatenation of values of these two keys
        represents the full path to the command

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

        return which(self.cmd) is not None  # True if tool available

    @abstractmethod
    def get_version(self):
        return


class ProdigalRunner(BaseRunner):

    name = "prodigal"

    _gc_run_results = defaultdict(dict)

    def __init__(self):
        super().__init__()
        self._output_folder = os.path.join(self.main_out, "prodigal_output")
        self._pred_genes_dir = os.path.join(self._output_folder, "predicted_genes")
        self._tmp_path = os.path.join(self._pred_genes_dir, "tmp")
        self.cpus = 1
        self.create_dirs([self._pred_genes_dir, self._tmp_path])

        # Get genetic_code from dataset.cfg file
        # bacteria/archaea=11; Entomoplasmatales,Mycoplasmatales=4
        try:
            self._genetic_code = self.config.get("prodigal", "prodigal_genetic_code").split(",")
        except NoOptionError:
            self._genetic_code = ["11"]

        # Set the ambiguous coding density range
        try:
            self._cd_upper = float(self.config.get("prodigal", "ambiguous_cd_range_upper")) \
                if len(self._genetic_code) > 1 else 0
        except NoOptionError:
            raise SystemExit("Dataset config file does not contain required information. Please upgrade datasets.")

        self.current_gc = None
        self._current_run_mode = None
        self._tmp_name = None

        self.output_faa = os.path.join(self._pred_genes_dir, "predicted.faa")
        self._output_fna = os.path.join(self._pred_genes_dir, "predicted.fna")
        self.sequences_aa = {}
        self.sequences_nt = {}
        self.gene_details = {}

        self._input_length = self._get_genome_length()
        self._run_mode = ["single", "meta"] if self._input_length > 100000 else ["meta"]

        self.init_checkpoint_file()
        self.run_number += 1

    @property
    def output_folder(self):
        return self._output_folder

    def check_tool_dependencies(self):
        pass

    def generate_job_args(self):
        yield

    @log("Genetic code {} selected as optimal", logger, attr_name="current_gc", on_func_exit=True)
    def run(self):
        """
        1) If genome length > 100000 run in "single" mode, then "meta" mode if there are no gene predictions. Otherwise
        just run in "meta" mode. This is based on the recommendations in the Prodigal docs.
        2) Run once using genetic code 11. This can be overridden if the user includes a spceific genetic code in the
        config file.
        3) Check the genome coding density. If the coding density is above the ambiguous range (typically 0.73-0.8)
        then continue with the current genetic code. The ambiguous range was determined based on analysis done by Mose
        Manni. Otherwise run again on the next genetic code specified.
        4) If the next run still has a genetic density within the ambiguous range, read the stdout log files (formerly
        the GFF files) and extract the scores assigned to each gene prediction. Whichever genetic code yields the
        greatest mean score is selected.
        :return:
        """
        super().run()
        tmp_files = []

        for ix, m in enumerate(self._run_mode):
            self._current_run_mode = m
            for g in self._genetic_code:
                self.current_gc = g

                file_id = os.path.join(self._tmp_path,
                                       "prodigal_mode_{0}_code_{1}".format(self._current_run_mode, self.current_gc))
                self._tmp_name = "{}.faa".format(file_id)
                self.logfile_path_out = "{}_out.log".format(file_id)
                self.logfile_path_err = "err".join(
                    self.logfile_path_out.rsplit("out", 1))  # Replace only the last occurence of "out" substring
                self._gc_run_results[self.current_gc].update({"tmp_name": self._tmp_name,
                                                             "log_file": self.logfile_path_out})

                if os.path.exists(self._tmp_name):  # Check to see if has already been run with these parameters
                    coding_density = self._gc_run_results[g]["cd"]
                else:
                    logger.info("Running Prodigal with genetic code {} in {} mode".format(self.current_gc,
                                                                                          self._current_run_mode))
                    self.total = 1
                    self.run_jobs()
                    coding_length = self._get_coding_length(self.logfile_path_out)
                    coding_density = coding_length / self._input_length
                    self._gc_run_results[self.current_gc].update({"cd": coding_density})

                logger.debug("Coding density is {}".format(coding_density))
                tmp_files.append(self._gc_run_results[self.current_gc]["tmp_name"])

                # if the coding density is above the ambiguous range, then continue with these parameters
                if coding_density >= self._cd_upper:
                    break

            # If output files from both runs in "single" mode are empty, run again in "meta" mode, else raise Exception.
            if not any([os.stat(tmp_file).st_size > 0 for tmp_file in tmp_files]):
                if ix + 1 == len(self._run_mode):
                    raise NoGenesError("Prodigal")
                else:
                    continue

            # if only one genetic code to consider, proceed with it
            # if there is more than one possible set of parameters, decide which to use
            self.current_gc = self._select_best_gc() if len(tmp_files) > 1 else self._genetic_code[0]

            selected_logfile = self._gc_run_results[self.current_gc]["log_file"]
            selected_tmpfile = self._gc_run_results[self.current_gc]["tmp_name"]

            self._organize_prodigal_files(selected_tmpfile, selected_logfile)
            self.get_gene_details()
            self._gc_run_results[self.current_gc].update({"seqs_aa": self.sequences_aa, "seqs_nt": self.sequences_nt,
                                                          "gene_details": self.gene_details})
            break

        return

    def configure_job(self, *args):
        tmp_name_nt = self._tmp_name.replace("faa", "fna")

        prodigal_job = self.create_job()
        prodigal_job.add_parameter("-p")
        prodigal_job.add_parameter("%s" % self._current_run_mode)
        prodigal_job.add_parameter("-f")
        prodigal_job.add_parameter("gff")
        prodigal_job.add_parameter("-g")
        prodigal_job.add_parameter("%s" % self.current_gc)
        prodigal_job.add_parameter("-a")
        prodigal_job.add_parameter("%s" % self._tmp_name)
        prodigal_job.add_parameter("-d")
        prodigal_job.add_parameter("%s" % tmp_name_nt)
        prodigal_job.add_parameter("-i")
        prodigal_job.add_parameter("%s" % self.input_file)

        return prodigal_job

    def get_gene_details(self):
        self.gene_details = defaultdict(list)

        with open(self._output_fna, "rU") as f:
            for record in SeqIO.parse(f, "fasta"):
                gene_name = record.id
                self.sequences_nt[gene_name] = record
                gene_start = int(record.description.split()[2])
                gene_end = int(record.description.split()[4])
                self.gene_details[gene_name].append({"gene_start": gene_start, "gene_end": gene_end})

        with open(self.output_faa, "rU") as f:
            for record in SeqIO.parse(f, "fasta"):
                self.sequences_aa[record.id] = record

        return

    @staticmethod
    def _get_coding_length(out_logfile):
        total_coding_length = 0
        with open(out_logfile, "r") as f:
            for line in f:
                if not line.startswith('#'):
                    try:
                        start = int(line.split('\t')[3])
                        stop = int(line.split('\t')[4])
                        total_coding_length += (stop-start)
                    except IndexError:
                        continue
                    except ValueError:
                        continue
        return total_coding_length

    def _get_genome_length(self):
        length_seqs = 0
        for line in open(self.input_file):
            if not line.startswith(">"):
                length_seqs += len(line)
        return length_seqs

    def _get_mean_score(self, gc):
        logfile = self._gc_run_results[gc]["log_file"]
        scores = []
        with open(logfile, "r") as f:
            for line in f:
                try:
                    score = re.search(";score=(.+?);", line).group(1)
                    scores.append(float(score))
                except AttributeError:
                    continue
        mean_score = (sum(scores) / len(scores))
        return mean_score

    def _organize_prodigal_files(self, tmp_file, tmp_logfile):

        shutil.copy(tmp_file, self.output_faa)
        shutil.copy(tmp_file.replace(".faa", ".fna"), self._output_fna)

        # copy selected log files from tmp/ to logs/
        new_logname = os.path.join(self.log_folder, "prodigal_out.log")
        shutil.copy(tmp_logfile, new_logname)
        shutil.copy(tmp_logfile.replace("_out.log", "_err.log"), new_logname.replace("_out.log", "_err.log"))
        return

    def _select_best_gc(self):
        gcs, cds = zip(*[[gc, self._gc_run_results[gc]["cd"]] for gc in self._genetic_code])
        sort_order = np.argsort(np.array(cds))[::-1]
        gcs_sorted = np.array(gcs)[sort_order]
        cds_sorted = np.array(cds)[sort_order]
        if abs(cds_sorted[0] - cds_sorted[1]) <= 0.05:
            mean_score1 = self._get_mean_score(gcs_sorted[0])
            mean_score2 = self._get_mean_score(gcs_sorted[1])
            gc = gcs_sorted[int(mean_score2 > mean_score1)]
        else:
            gc = gcs_sorted[0]
        return gc

    def get_version(self):
        prodigal_version = subprocess.check_output([self.cmd, "-v"], stderr=subprocess.STDOUT, shell=False)
        prodigal_version = prodigal_version.decode("utf-8")
        prodigal_version = prodigal_version.split("\n")[1].split(":")[0]
        prodigal_version = prodigal_version.replace("Prodigal V", "")
        return prodigal_version


class NoGenesError(Exception):

    def __init__(self, gene_predictor):
        self.gene_predictor = gene_predictor


class HMMERRunner(BaseRunner):

    name = "hmmsearch"

    def __init__(self):
        super().__init__()
        self._hmmer_output_folder = os.path.join(self.run_folder, "hmmer_output")
        self.datasets_version = self.config.get("busco_run", "datasets_version")
        self.dataset_creation_date = self.config.get("busco_run", "creation_date")
        self.dataset_nb_species = self.config.get("busco_run", "number_of_species")
        self.dataset_nb_buscos = self.config.get("busco_run", "number_of_BUSCOs")
        self.domain = self.config.get("busco_run", "domain")

        self.single_copy_sequences_folder = os.path.join(self.run_folder, "busco_sequences",
                                                         "single_copy_busco_sequences")
        self.multi_copy_sequences_folder = os.path.join(self.run_folder, "busco_sequences",
                                                        "multi_copy_busco_sequences")
        self.fragmented_sequences_folder = os.path.join(self.run_folder, "busco_sequences",
                                                        "fragmented_busco_sequences")
        self.short_summary_file = os.path.join(self.run_folder, "short_summary.txt")
        self.cutoff_dict = {}
        self.extra_columns = False
        self.log_count = 0  # Dummy variable used to skip logging for intermediate eukaryote pipeline results.
        self.one_line_summary = None

        # to be initialized before run time
        self.input_sequences = None
        self.busco_ids = None
        self.mode = None
        self.gene_details = None
        self.results_dir = None

        self.create_dirs([self._hmmer_output_folder, self.single_copy_sequences_folder,
                          self.multi_copy_sequences_folder, self.fragmented_sequences_folder])
        if self.domain == "eukaryota":
            self.initial_results_dir = os.path.join(self._hmmer_output_folder, "initial_run_results")
            self.rerun_results_dir = os.path.join(self._hmmer_output_folder, "rerun_results")
            self.create_dirs([self.initial_results_dir, self.rerun_results_dir])

        self.init_checkpoint_file()

    def configure_runner(self, input_sequences, busco_ids, mode, gene_details):
        self.run_number += 1
        self.input_sequences = input_sequences
        self.busco_ids = busco_ids
        self.mode = mode
        self.single_copy_buscos = {}
        self.multi_copy_buscos = {}
        self.fragmented_buscos = {}
        self.is_complete = {}
        self.is_fragment = {}
        self.is_very_large = {}
        self.matched_bitscores = {}
        self.matched_genes_complete = {}
        self.matched_genes_vlarge = {}
        self.matched_genes_fragment = {}
        self._already_used_genes = set()
        self.hmmer_results_lines = []
        self.missing_buscos = []
        self.gene_details = gene_details
        if len(self.cutoff_dict) == 0:
            self.load_buscos()

        if self.domain == "eukaryota":
            if self.run_number == 1:
                self.results_dir = self.initial_results_dir
            elif self.run_number == 2:
                self.results_dir = self.rerun_results_dir
            else:
                raise ValueError("HMMER should not be run more than twice in the same Run instance.")
        else:
            self.results_dir = self._hmmer_output_folder
        # gene_details can only be None for proteins mode. In the other modes the gene locations are written to a file
        # after the coordinates are loaded from this attribute

    def configure_job(self, busco_id, seq_filename, output_filename):

        hmmer_job = self.create_job()
        hmmer_job.add_parameter("--domtblout")
        hmmer_job.add_parameter(os.path.join(self.results_dir, output_filename))
        hmmer_job.add_parameter("--cpu")
        hmmer_job.add_parameter("1")
        hmmer_job.add_parameter(os.path.join(self.lineage_dataset, "hmms", "{}.hmm".format(busco_id)))
        hmmer_job.add_parameter(seq_filename)
        return hmmer_job

    def generate_job_args(self):
        for busco_id in self.busco_ids:
            if busco_id in self.cutoff_dict:
                if isinstance(self.input_sequences, str):
                    output_filename = "{}.out".format(busco_id)
                    yield busco_id, self.input_sequences, output_filename
                elif isinstance(self.input_sequences, list):
                    input_files = [f for f in self.input_sequences if os.path.basename(f).startswith(busco_id)]
                    for seq_filename in input_files:
                        output_filename = os.path.basename(seq_filename).replace("faa", "out")
                        yield busco_id, seq_filename, output_filename

    @property
    def output_folder(self):
        return self._hmmer_output_folder

    def load_buscos(self):
        """
        Load all BUSCOs for the lineage, along with their cutoff lengths and scores.
        :return:
        """
        self.cutoff_dict = defaultdict(dict)
        self._load_length()
        self._load_score()
        self.cutoff_dict = dict(self.cutoff_dict)
        return

    def run(self):
        """
        Create a HMMER job for each BUSCO. Each job searches the input sequence file for matches for the BUSCO gene.
        :return:
        """
        super().run()
        self.total = self._count_jobs()
        self.run_jobs()

    def _count_jobs(self):
        n = 0
        for busco_id in self.busco_ids:
            if busco_id in self.cutoff_dict:
                if isinstance(self.input_sequences, str):
                    n += 1
                elif isinstance(self.input_sequences, list):
                    input_files = [f for f in self.input_sequences if os.path.basename(f).startswith(busco_id)]
                    n += len(input_files)
        return n

    def get_version(self):
        """
        check the Tool has the correct version
        :raises SystemExit: if the version is not correct
        """
        hmmer_version = subprocess.check_output([self.cmd, "-h"], stderr=subprocess.STDOUT, shell=False)
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

    def check_tool_dependencies(self):
        """
        check dependencies on tools
        :raises SystemExit: if a Tool version is not supported
        """
        # check hmm version
        if not self.version >= BuscoConfig.HMMER_VERSION:
            raise SystemExit(
                "HMMer version detected is not supported, please use HMMer v.{} +".format(BuscoConfig.HMMER_VERSION))
        return

    def process_output(self):
        # Re-initialize dictionaries as defaultdicts - necessary because defaultdicts are not picklable and so they
        # cannot appear in the __init__ when using multiprocessing within the class
        self.matched_genes_complete = defaultdict(list)
        self.matched_genes_vlarge = defaultdict(list)
        self.matched_genes_fragment = defaultdict(list)
        self.matched_bitscores = defaultdict(lambda: defaultdict(list))
        self.is_complete = defaultdict(lambda: defaultdict(list))  # dict of a dict of lists of dicts
        self.is_fragment = defaultdict(lambda: defaultdict(list))
        self.is_very_large = defaultdict(lambda: defaultdict(list))

        self._load_matched_genes()
        self._filter()
        self._consolidate_busco_lists()

        self.matched_bitscores = dict(self.matched_bitscores)
        self.matched_genes_complete = dict(self.matched_genes_complete)
        self.matched_genes_vlarge = dict(self.matched_genes_vlarge)
        self.matched_genes_fragment = dict(self.matched_genes_fragment)
        self.is_complete = dict(self.is_complete)
        self.is_fragment = dict(self.is_fragment)
        self.is_very_large = dict(self.is_very_large)
        return

    @staticmethod
    def _get_matched_lengths(nested_dict):
        """
        For each entry in a nested dictionary, return a dict with the total lengths of all gene matches for each entry.
        :param nested_dict:
        :type nested_dict:
        :return:
        :rtype:
        """
        total_len = defaultdict(int)
        for entry in nested_dict:
            for hit in nested_dict[entry]:
                total_len[entry] += hit[1] - hit[0]
        return total_len

    def _parse_hmmer_output(self, filename, busco_query):
        """
        Read and parse HMMER output file.
        :param filename: Name of HMMER output file
        :param busco_query: Basename of file, used to identify BUSCO
        :type filename: str
        :type busco_query: str
        :return: Dictionary of (gene_id, total_matched_length) pairs
        :rtype: dict
        """
        matched_lengths = defaultdict(int)

        with open(filename, "r") as f:

            # Read HMMER output file
            for line in f:
                if line.startswith("#"):
                    continue
                else:
                    try:
                        line = line.strip().split()
                        gene_id = line[0]
                        bit_score = float(line[7])
                        hmm_start = int(line[15])
                        hmm_end = int(line[16])

                        # Store bitscore matches for each gene match. If match below cutoff, discard.
                        if bit_score >= float(self.cutoff_dict[busco_query]["score"]):
                            # todo: introduce upper bound - consult to see what a reasonable value would be
                            self.matched_bitscores[busco_query][gene_id].append(bit_score)
                        else:
                            continue

                        matched_lengths[gene_id] += (hmm_end - hmm_start)

                    except IndexError as e:
                        SystemExit(e, "Cannot parse HMMER output file {}".format(filename))
        return matched_lengths

    def _sort_matches(self, matched_lengths, busco_query):
        """
        The HMMER gene matches are sorted into "complete", "v_large" and "fragmented" matches based on a comparison
        with the cutoff value specified in the dataset cutoff_scores file
        :param matched_lengths: dict of (gene_id, total_matched_length) pairs
        :param busco_query: BUSCO identifier
        :type matched_lengths: dict
        :type busco_query: str
        :return: busco_complete, busco_vlarge, busco_fragment - three dictionaries of the form
        {gene_id: [{"bitscore": float, "length": int}, {...}, ...], ...}
        :rtype: dict
        """
        busco_complete = defaultdict(list)
        busco_vlarge = defaultdict(list)
        busco_fragment = defaultdict(list)

        # Determine whether matched gene represents a complete, very_large or fragment of a BUSCO
        for gene_id, size in matched_lengths.items():

            # Kind of like a z-score, but it is compared with a cutoff value, not a mean
            zeta = (self.cutoff_dict[busco_query]["length"] - size) \
                   / self.cutoff_dict[busco_query]["sigma"]

            # gene match can only be either complete, v_large or fragment
            if -2 <= zeta <= 2:
                busco_type = busco_complete
                match_type = self.matched_genes_complete
            elif zeta < -2:
                busco_type = busco_vlarge
                match_type = self.matched_genes_vlarge
            else:
                busco_type = busco_fragment
                match_type = self.matched_genes_fragment

            # Add information about match to dict
            busco_type[gene_id].append(dict({"bitscore": max(self.matched_bitscores[busco_query][gene_id]),
                                             "length": matched_lengths[gene_id]}))
            # Reference which busco_queries are associated with each gene match
            match_type[gene_id].append(busco_query)

        return busco_complete, busco_vlarge, busco_fragment

    def _load_matched_genes(self):
        """
        Load all gene matches from HMMER output and sort into dictionaries depending on match quality
        (complete, v_large, fragment).
        :return:
        """
        if self.run_number == 1:
            hmmer_results_files = sorted([os.path.join(self.results_dir, f) for f in os.listdir(self.results_dir)])
        elif self.run_number == 2:
            hmmer_initial_run_files = [os.path.join(self.initial_results_dir, f)
                                       for f in os.listdir(self.initial_results_dir)]
            hmmer_rerun_files = [os.path.join(self.rerun_results_dir, f)
                                 for f in os.listdir(self.rerun_results_dir)]
            hmmer_results_files = sorted(hmmer_initial_run_files + hmmer_rerun_files)
        else:
            raise ValueError("HMMER should not be run more than twice in the same Run instance.")

        for filename in hmmer_results_files:
            busco_query = str(os.path.basename(filename).split(".")[0])
            matched_lengths = self._parse_hmmer_output(filename, busco_query)
            busco_complete, busco_vlarge, busco_fragment = self._sort_matches(matched_lengths, busco_query)

            # Add all information for this busco_id to the full dictionary
            if len(busco_complete) > 0:
                self.is_complete[busco_query].update(busco_complete)
            if len(busco_vlarge) > 0:
                self.is_very_large[busco_query].update(busco_vlarge)
            if len(busco_fragment) > 0:
                self.is_fragment[busco_query].update(busco_fragment)

        return

    def _update_used_gene_set(self, busco_dict):
        """
        Update set of already used genes to prevent processing the same gene twice.
        :param busco_dict: One of [self.is_complete, self.is_very_large, self.is_fragment]
        :type busco_dict: dict
        :return:
        """
        for entries in busco_dict.values():
            for gene_id in entries:
                self._already_used_genes.add(gene_id)
        return

    def _remove_lower_ranked_duplicates(self, busco_dict):
        """
        Remove any genes and/or busco matches from input dictionary if they have previously been assigned to a better
        quality match.
        :param busco_dict: one of [self.is_very_large, self.is_fragment]
        :type busco_dict: dict
        :return:
        """
        # Determine which match ranks to worry about
        if busco_dict == self.is_very_large:
            higher_rank_buscos = self.is_complete.keys()
            matched_genes = self.matched_genes_vlarge
        elif busco_dict == self.is_fragment:
            higher_rank_buscos = list(self.is_complete.keys()) + list(self.is_very_large.keys())
            matched_genes = self.matched_genes_fragment
        else:
            raise SystemExit("Unrecognized dictionary of BUSCOs.")

        for busco_id in list(busco_dict.keys()):
            matches = busco_dict[busco_id]
            # Remove any buscos that appear in higher ranking dictionaries
            if busco_id in higher_rank_buscos:
                busco_dict.pop(busco_id)
                for gene_id in matches:
                    matched_genes[gene_id] = [x for x in matched_genes[gene_id] if x != busco_id]  # Remove all occurences of busco_id
                    if len(matched_genes[gene_id]) == 0:
                        matched_genes.pop(gene_id)
                continue

            # Remove any genes that have previously been processed under a different and higher ranking busco match
            for gene_id in list(matches.keys()):
                if gene_id in self._already_used_genes:
                    busco_dict[busco_id].pop(gene_id)
                    matched_genes[gene_id] = [x for x in matched_genes[gene_id] if x != busco_id]  # Remove all occurences of busco_id
                    if len(busco_dict[busco_id]) == 0:
                        busco_dict.pop(busco_id)
                    if len(matched_genes[gene_id]) == 0:
                        matched_genes.pop(gene_id)

        return

    def _remove_duplicates(self):
        """
        Remove duplicate gene matches of lesser importance, i.e. keep the complete ones, then the very large ones and
        finally the fragments.
        Also remove duplicate BUSCO ID matches of lower importance.
        Then search for any duplicate gene matches within the same rank for different BUSCOs and keep only the highest
        scoring gene match.
        :return:
        """
        self._update_used_gene_set(self.is_complete)
        self._remove_lower_ranked_duplicates(self.is_very_large)
        self._update_used_gene_set(self.is_very_large)
        self._remove_lower_ranked_duplicates(self.is_fragment)
        self._remove_remaining_duplicate_matches(self.is_complete)
        self._remove_remaining_duplicate_matches(self.is_very_large)
        self._remove_remaining_duplicate_matches(self.is_fragment)
        return

    def _remove_remaining_duplicate_matches(self, busco_dict):
        """
        For any genes matched under more than one BUSCO, keep only the highest scoring match in the input dictionary.
        :param busco_dict: one of [self.is_complete, self.is_very_large, self.is_fragment]
        :type busco_dict: dict
        :return:
        """
        # For a given input dictionary {busco_id: gene_ids}, make sure we are using the corresponding dictionary
        # {gene_id: busco_matches}
        if busco_dict == self.is_complete:
            matched_genes = self.matched_genes_complete
        elif busco_dict == self.is_very_large:
            matched_genes = self.matched_genes_vlarge
        elif busco_dict == self.is_fragment:
            matched_genes = self.matched_genes_fragment
        else:
            raise SystemExit("Unrecognized dictionary of BUSCOs.")

        # Keep the best scoring gene if gene is matched by more than one busco with the same match rank
        for gene_id, buscos in matched_genes.items():
            if len(buscos) > 1:
                busco_bitscores = []
                busco_matches = []
                for busco in buscos:
                    matches = busco_dict[busco][gene_id]
                    for match in matches:
                        bitscore = match["bitscore"]
                        busco_bitscores.append(bitscore)
                        busco_matches.append(busco)

                if len(set(buscos)) == 1:  # If only one busco is matched twice (initial run and rerun), don't remove it
                    continue
                best_match_ind = max(range(len(busco_bitscores)), key=busco_bitscores.__getitem__)
                buscos = [x for x in buscos if x != busco_matches[best_match_ind]]
                # Remove lower scoring duplicates from dictionary.
                # Note for future development: the matched_genes dictionary is not updated in this method when
                # duplicates are removed from busco_dict
                for duplicate in list(set(buscos)):
                    # Use set to account for any duplicate entries (matched in both initial run and rerun)
                    busco_dict[duplicate].pop(gene_id)
                    if len(busco_dict[duplicate]) == 0:
                        busco_dict.pop(duplicate)
        return

    def _remove_low_scoring_matches(self, busco_dict):
        """
        Go through input dictionary and remove any gene matches that score less than 85% of the top gene match score
        for each BUSCO.
        :param busco_dict: one of [self.is_complete, self.is_very_large, self.is_fragment]
        :type busco_dict: dict
        :return:
        """
        empty_buscos = []

        # For each busco, keep only matches within 85% of top bitscore match for that busco
        for busco_id, matches in busco_dict.items():
            if len(matches) > 1:
                _, max_bitscore = self._get_best_scoring_match(matches)
                # Go through all matches again, removing any below the threshold
                for gene_id in list(matches.keys()):
                    match_info = matches[gene_id]
                    matches_to_remove = []
                    for m, match in enumerate(match_info):
                        if match["bitscore"] < 0.85*max_bitscore:
                            matches_to_remove.append(m)

                    # Remove dict from list of dicts. Safe way to delete without risking list size changing during
                    # iteration
                    for ind in sorted(matches_to_remove, reverse=True):
                        del match_info[ind]

                    # Record dictionary address of empty gene records
                    if len(busco_dict[busco_id][gene_id]) == 0:
                        empty_buscos.append((busco_id, gene_id))

        # Safe way to delete empty records without risking dictionary size changing while iterating
        for item in empty_buscos:
            busco_id, gene_id = item
            busco_dict[busco_id].pop(gene_id)

        return

    @staticmethod
    def _get_best_scoring_match(gene_matches):
        """
        Find the highest bitscore in all gene matches.
        :param gene_matches: dictionary of the form
        {gene_id: [{"bitscore": float, "length": int}, {"bitscore": float, "length": int}, ...], ...}
        :type gene_matches: dict
        :return: best_match_gene, best_match_bitscore
        :rtype: str, float
        """
        match_scores = []
        match_genes = []
        for gene_id, matches in gene_matches.items():
            for match in matches:
                bitscore = match["bitscore"]
                match_scores.append(bitscore)
                match_genes.append(gene_id)
        best_match_ind = max(range(len(match_scores)), key=match_scores.__getitem__)
        best_match_gene = match_genes[best_match_ind]
        best_match_bitscore = match_scores[best_match_ind]
        return best_match_gene, best_match_bitscore

    def _filter(self):
        """
        Remove all duplicate matches and any matches below 85% of the top match for each BUSCO.
        :return:
        """
        self._remove_duplicates()
        self._remove_low_scoring_matches(self.is_complete)
        self._remove_low_scoring_matches(self.is_very_large)
        self._remove_low_scoring_matches(self.is_fragment)
        return

    def _consolidate_busco_lists(self):
        """
        Sort BUSCO matches into single-copy, multi-copy and fragments.
        Only the highest scoring fragment for each BUSCO is kept.
        :return:
        """
        for busco_dict in [self.is_complete, self.is_very_large]:
            for busco_id, gene_matches in busco_dict.items():
                if len(gene_matches) == 1:
                    self.single_copy_buscos[busco_id] = busco_dict[busco_id]
                else:
                    self.multi_copy_buscos[busco_id] = busco_dict[busco_id]

        for busco_id, gene_matches in self.is_fragment.items():
            if len(gene_matches) > 1:
                best_fragment, _ = self._get_best_scoring_match(gene_matches)
                self.fragmented_buscos[busco_id] = {best_fragment: self.is_fragment[busco_id][best_fragment]}
            else:
                self.fragmented_buscos[busco_id] = gene_matches
        return

    def load_links_info(self):
        links_info = defaultdict(dict)
        links_file = os.path.join(self.lineage_dataset, "links_to_{}.txt".format(self.datasets_version.upper()))
        if os.path.exists(links_file):
            with open(links_file, newline='') as f:
                contents = csv.reader(f, delimiter="\t")
                for row in contents:
                    busco_id, description, link = row
                    links_info[busco_id]["description"] = description
                    links_info[busco_id]["link"] = link
        return links_info

    def _format_output_lines(self, busco_dict, label):
        """
        Format BUSCO matches from input dictionary into output lines for writing to a file.
        :param busco_dict: one of [self.single_copy_buscos, self.multi_copy_buscos, self.fragmented_buscos]
        :type busco_dict: dict
        :return: output_lines
        :rtype: list
        """
        output_lines = []

        links_info = self.load_links_info()

        for busco, matches in busco_dict.items():
            for gene_id, match_info in matches.items():
                for m, match in enumerate(match_info):
                    bit_score = match["bitscore"]
                    match_length = match["length"]

                    if self.mode == "proteins" or self.mode == "transcriptome":
                        try:
                            desc = links_info[busco]["description"]
                            link = links_info[busco]["link"]
                            self.extra_columns = True
                            output_lines.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(busco, label, gene_id, bit_score,
                                                                                      match_length, link, desc))
                        except KeyError:
                            output_lines.append("{}\t{}\t{}\t{}\t{}\n".format(busco, label, gene_id, bit_score,
                                                                              match_length))
                    elif self.mode == "genome":
                        scaffold = self.gene_details[gene_id][m]
                        location_pattern = ":{}-{}".format(scaffold["gene_start"], scaffold["gene_end"])
                        if gene_id.endswith(location_pattern):
                            gene_id = gene_id.replace(location_pattern, "")
                        try:
                            desc = links_info[busco]["description"]
                            link = links_info[busco]["link"]
                            self.extra_columns = True
                            output_lines.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                busco, label, gene_id, scaffold["gene_start"], scaffold["gene_end"], bit_score,
                                match_length, link, desc))
                        except KeyError:
                            output_lines.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                busco, label, gene_id, scaffold["gene_start"], scaffold["gene_end"], bit_score,
                                match_length))
        return output_lines

    def _create_output_content(self):
        """
        Format output for all BUSCO matches.
        :return: output_lines
        :rtype: list
        """
        output_lines = []
        dict_labels = {"Complete": self.single_copy_buscos,
                       "Duplicated": self.multi_copy_buscos,
                       "Fragmented": self.fragmented_buscos}
        for label, busco_dict in dict_labels.items():
            output_lines += self._format_output_lines(busco_dict, label)

        return output_lines

    def _list_missing_buscos(self):
        """
        Create a list of all BUSCOs that are missing after processing the HMMER output.
        :return: output_lines, missing_buscos
        :rtype: list, list
        """
        output_lines = []
        for busco_group in self.cutoff_dict:
            if not any(busco_group in d for d in [self.is_complete, self.is_very_large, self.is_fragment]):
                output_lines.append("{}\tMissing\n".format(busco_group))
                self.missing_buscos.append(busco_group)

        if len(self.missing_buscos) == len(self.cutoff_dict):
            logger.warning("BUSCO did not find any match. Make sure to check the log files if this is unexpected.")

        return output_lines, self.missing_buscos

    def _load_length(self):
        """
        This function loads the length cutoffs file
        :raises SystemExit: if the lengths_cutoff file cannot be read
        """
        lengths_cutoff_file = os.path.join(self.lineage_dataset, "lengths_cutoff")
        try:
            with open(lengths_cutoff_file, "r") as f:
                for line in f:
                    line = line.strip().split()
                    try:
                        taxid = line[0]
                        sd = float(line[2])
                        length = float(line[3])

                        self.cutoff_dict[taxid]["sigma"] = sd
                        # there is an arthropod profile with sigma 0
                        # that causes a crash on divisions
                        if sd == 0.0:
                            self.cutoff_dict[taxid]["sigma"] = 1
                        self.cutoff_dict[taxid]["length"] = length
                    except IndexError as e:
                        raise SystemExit(e, "Error parsing the lengths_cutoff file.")
        except IOError:
            raise SystemExit("Impossible to read the lengths in {}".format(os.path.join(lengths_cutoff_file)))
        return

    def _load_score(self):
        """
        This function loads the score cutoffs file
        :raises SystemExit: if the scores_cutoff file cannot be read
        """
        scores_cutoff_file = os.path.join(self.lineage_dataset, "scores_cutoff")
        try:
            # open target scores file
            with open(scores_cutoff_file, "r") as f:
                for line in f:
                    line = line.strip().split()
                    try:
                        taxid = line[0]
                        score = float(line[1])
                        self.cutoff_dict[taxid]["score"] = score
                    except IndexError as e:
                        raise SystemExit(e, "Error parsing the scores_cutoff file.")
        except IOError:
            raise SystemExit("Impossible to read the scores in {}".format(scores_cutoff_file))
        return

    def write_buscos_to_file(self, sequences_aa, sequences_nt=None):
        """
        Write BUSCO matching sequences to output fasta files. Each sequence is printed in a separate file and both
        nucleotide and amino acid versions are created.
        :param sequences_aa: dict
        :param sequences_nt: dict
        :return:
        """
        for busco_type in ["single_copy", "multi_copy", "fragmented"]:
            if busco_type == "single_copy":
                output_dir = self.single_copy_sequences_folder
                busco_matches = self.single_copy_buscos
            elif busco_type == "multi_copy":
                output_dir = self.multi_copy_sequences_folder
                busco_matches = self.multi_copy_buscos
            elif busco_type == "fragmented":
                output_dir = self.fragmented_sequences_folder
                busco_matches = self.fragmented_buscos

            for busco, gene_matches in busco_matches.items():
                try:
                    aa_seqs, nt_seqs = zip(*[(sequences_aa[gene_id], sequences_nt[gene_id])
                                             for gene_id in gene_matches])
                    with open(os.path.join(output_dir, "{}.fna".format(busco)), "w") as f2:
                        SeqIO.write(nt_seqs, f2, "fasta")
                except TypeError:
                    aa_seqs = [sequences_aa[gene_id] for gene_id in gene_matches]
                with open(os.path.join(output_dir, "{}.faa".format(busco)), "w") as f1:
                    SeqIO.write(aa_seqs, f1, "fasta")

    def write_hmmer_results(self):
        """
        Create two output files: one with information on all BUSCOs for the given dataset and the other with a list of
        all BUSCOs that were not found.
        :return:
        """

        with open(os.path.join(self.run_folder, "full_table.tsv"), "w") as f_out:

            output_lines = self._create_output_content()
            self._write_output_header(f_out)

            with open(os.path.join(self.run_folder, "missing_busco_list.tsv"), "w") as miss_out:

                self._write_output_header(miss_out, missing_list=True)

                # todo: move to calculate busco percentages
                missing_buscos_lines, missing_buscos = self._list_missing_buscos()
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

    def produce_hmmer_summary(self):
        single_copy, multi_copy, only_fragments, total_buscos = self._get_busco_percentages()

        self.hmmer_results_lines.append("***** Results: *****\n\n")
        self.one_line_summary = "C:{}%[S:{}%,D:{}%],F:{}%,M:{}%,n:{}\t{}\n".format(
            round(self.s_percent + self.d_percent, 1), self.s_percent, self.d_percent,
            self.f_percent, abs(round(100 - self.s_percent - self.d_percent - self.f_percent, 1)), total_buscos, "   ")
        self.hmmer_results_lines.append(self.one_line_summary)
        self.hmmer_results_lines.append("{}\tComplete BUSCOs (C)\t\t\t{}\n".format(single_copy + multi_copy, "   "))
        self.hmmer_results_lines.append("{}\tComplete and single-copy BUSCOs (S)\t{}\n".format(single_copy, "   "))
        self.hmmer_results_lines.append("{}\tComplete and duplicated BUSCOs (D)\t{}\n".format(multi_copy, "   "))
        self.hmmer_results_lines.append("{}\tFragmented BUSCOs (F)\t\t\t{}\n".format(only_fragments, "   "))
        self.hmmer_results_lines.append("{}\tMissing BUSCOs (M)\t\t\t{}\n".format(
            total_buscos - single_copy - multi_copy - only_fragments, "   "))
        self.hmmer_results_lines.append("{}\tTotal BUSCO groups searched\t\t{}\n".format(total_buscos, "   "))

        if isinstance(self.config, BuscoConfigAuto):
            self._one_line_hmmer_summary()
        elif self.domain == "eukaryota" and self.log_count == 0:
            self.log_count += 1
            self._produce_full_hmmer_summary_debug()
        else:
            self._one_line_hmmer_summary()

        with open(self.short_summary_file, "w") as summary_file:

            self._write_output_header(summary_file, no_table_header=True)
            summary_file.write("# Summarized benchmarking in BUSCO notation for file {}\n"
                               "# BUSCO was run in mode: {}\n\n".format(self.input_file, self.mode))

            for line in self.hmmer_results_lines:
                summary_file.write("\t{}".format(line))

            if self.config.getboolean("busco_run", "auto-lineage") and isinstance(self.config, BuscoConfigMain) \
                    and hasattr(self.config, "placement_files"):
                summary_file.write("\nPlacement file versions:\n")
                for placement_file in self.config.placement_files:
                    summary_file.write("{}\n".format(placement_file))

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

    def _write_output_header(self, file_object, missing_list=False, no_table_header=False):
        """
        Write a standardized file header containing information on the BUSCO run.
        :param file_object: Opened file object ready for writing
        :type file_object: file
        :return:
        """
        file_object.write("# BUSCO version is: {} \n"
                          "# The lineage dataset is: {} (Creation date: {}, number of species: {}, number of BUSCOs: {}"
                          ")\n".format(busco.__version__, os.path.basename(self.lineage_dataset),
                                       self.dataset_creation_date, self.dataset_nb_species, self.dataset_nb_buscos))
        # if isinstance(self._config, BuscoConfigMain):  # todo: wait until rerun command properly implemented again
        #     file_object.write("# To reproduce this run: {}\n#\n".format(self._rerun_cmd))

        if no_table_header:
            pass
        elif missing_list:
            file_object.write("# Busco id\n")
        elif self.mode == "proteins" or self.mode == "transcriptome":
            if self.extra_columns:
                file_object.write("# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
            else:
                file_object.write("# Busco id\tStatus\tSequence\tScore\tLength\n")
        elif self.mode == "genome":
            if self.extra_columns:
                file_object.write(
                    "# Busco id\tStatus\tSequence\tGene Start\tGene End\tScore\tLength\tOrthoDB url\tDescription\n")
            else:
                file_object.write("# Busco id\tStatus\tSequence\tGene Start\tGene End\tScore\tLength\n")

        return

    def _get_busco_percentages(self):
        self.single_copy = len(self.single_copy_buscos)  # int
        self.multi_copy = len(self.multi_copy_buscos)  # int
        self.only_fragments = len(self.fragmented_buscos)  # int
        self.total_buscos = len(self.cutoff_dict)

        # Get percentage of each kind of BUSCO match
        self.s_percent = abs(round((self.single_copy / self.total_buscos) * 100, 1))
        self.d_percent = abs(round((self.multi_copy / self.total_buscos) * 100, 1))
        self.f_percent = abs(round((self.only_fragments / self.total_buscos) * 100, 1))

        return self.single_copy, self.multi_copy, self.only_fragments, self.total_buscos


class MKBLASTRunner(BaseRunner):

    name = "makeblastdb"

    def __init__(self):
        super().__init__()
        self.db_path = os.path.join(self.config.get("busco_run", "main_out"), "blast_db")
        self.output_db = os.path.join(self.db_path, os.path.basename(self.input_file))
        self.create_dirs(self.db_path)
        self.total = 1
        self.init_checkpoint_file()
        self.run_number += 1

    @log("Creating BLAST database with input file", logger)
    def configure_job(self, *args):
        mkblast_job = self.create_job()
        mkblast_job.add_parameter("-in")
        mkblast_job.add_parameter(self.input_file)
        mkblast_job.add_parameter("-dbtype")
        mkblast_job.add_parameter("nucl")
        mkblast_job.add_parameter("-out")
        mkblast_job.add_parameter(self.output_db)
        return mkblast_job

    def run(self):
        super().run()
        if os.path.exists(self.db_path) and len(os.listdir(self.db_path)) > 0:
            return

        self.run_jobs()

    def generate_job_args(self):
        yield

    def get_version(self):
        mkblastdb_version_call = subprocess.check_output([self.cmd, "-version"], stderr=subprocess.STDOUT, shell=False)
        mkblastdb_version = ".".join(mkblastdb_version_call.decode("utf-8").split("\n")[0].split()[1].rsplit("."))

        return mkblastdb_version

    def check_tool_dependencies(self):
        pass

    @property
    def output_folder(self):
        return self.db_path


class TBLASTNRunner(BaseRunner):

    name = "tblastn"

    MAX_FLANK = 20000

    def __init__(self):
        self.coords = {}
        super().__init__()
        self._output_folder = os.path.join(self.run_folder, "blast_output")
        self.output_seqs = os.path.join(self._output_folder, "sequences")
        self.create_dirs([self._output_folder, self.output_seqs])
        self.total = 1

        self.e_v_cutoff = self.config.getfloat("busco_run", "evalue")
        self.region_limit = self.config.getint("busco_run", "limit")
        self.flank = self._define_flank()

        self.init_checkpoint_file()

    def configure_runner(self, blast_db, missing_and_frag_only, ancestral_variants, incomplete_buscos):
        self.run_number += 1
        self.blast_db = blast_db
        self.missing_and_frag_only = missing_and_frag_only
        self.ancestral_variants = ancestral_variants
        self.incomplete_buscos = incomplete_buscos

        self.ancestral_sfx = "_variants" if self.ancestral_variants else ""
        self.ancestral_file = os.path.join(self.lineage_dataset, "ancestral{}".format(self.ancestral_sfx))
        self.query_file = os.path.join(self.lineage_dataset, "ancestral{}".format(self.ancestral_sfx))
        self.output_suffix = "_missing_and_frag_rerun" if self.missing_and_frag_only else ""
        self.rerun_query_file = os.path.join(self._output_folder,
                                             "ancestral{}{}".format(self.ancestral_sfx, self.output_suffix))
        if self.missing_and_frag_only and self.ancestral_variants:
            self._extract_incomplete_buscos_ancestral()

        self.blast_filename = os.path.join(self._output_folder, "tblastn{}.tsv".format(self.output_suffix))
        self.coords_filename = os.path.join(self._output_folder, "coordinates{}.tsv".format(self.output_suffix))

    def configure_job(self, *args):
        tblastn_job = self.create_job()
        tblastn_job.add_parameter("-evalue")
        tblastn_job.add_parameter(str(self.e_v_cutoff))
        tblastn_job.add_parameter("-num_threads")
        tblastn_job.add_parameter(str(self.cpus))
        tblastn_job.add_parameter("-query")
        tblastn_job.add_parameter(self.query_file)
        tblastn_job.add_parameter("-db")
        tblastn_job.add_parameter(self.blast_db)
        tblastn_job.add_parameter("-out")
        tblastn_job.add_parameter(self.blast_filename)
        tblastn_job.add_parameter("-outfmt")
        tblastn_job.add_parameter("7")
        return tblastn_job

    @property
    def output_folder(self):
        return self._output_folder

    def _define_flank(self):
        """
        TODO: Add docstring
        :return:
        """
        try:
            size = os.path.getsize(self.input_file) / 1000  # size in mb
            flank = int(size / 50)  # proportional flank size
            # Ensure value is between 5000 and MAX_FLANK
            flank = min(max(flank, 5000), type(self).MAX_FLANK)
        except IOError:  # Input data is only validated during run_analysis. This will catch any IO issues before that.
            raise SystemExit("Impossible to read the fasta file {}".format(self.input_file))

        return flank

    @log("Running a BLAST search for BUSCOs against created database", logger)
    def run(self):
        super().run()
        self.run_jobs()
        self._check_output()
        return

    def check_tool_dependencies(self):
        if ".".join(self.version.split(".")[:-1]) not in ["2.2", "2.3"] and self.version != "2.10.1+":
            # Known problems with multithreading on BLAST 2.4-2.10.0.
            logger.warning("You are using BLAST version {}. This is known to yield inconsistent results when "
                           "multithreading. BLAST will run on a single core as a result. For performance improvement, "
                           "please upgrade to BLAST 2.10.1+.".format(self.version))
            self.cpus = 1

    def get_version(self):
        tblastn_version_call = subprocess.check_output([self.cmd, "-version"], stderr=subprocess.STDOUT, shell=False)
        tblastn_version = ".".join(tblastn_version_call.decode("utf-8").split("\n")[0].split()[1].rsplit("."))

        return tblastn_version

    def generate_job_args(self):
        yield

    def _check_output(self):
        # check that blast worked
        if not os.path.exists(self.blast_filename):
            raise SystemExit("tblastn failed!")

        # check that the file is not truncated
        with open(self.blast_filename, "r") as f:
            try:
                if "processed" not in f.readlines()[-1]:
                    raise SystemExit("tblastn has ended prematurely (the result file lacks the expected final line), "
                                     "which will produce incomplete results in the next steps ! This problem likely "
                                     "appeared in blast+ 2.4 and seems not fully fixed in 2.6. It happens only when "
                                     "using multiple cores. You can use a single core (-c 1) or downgrade to "
                                     "blast+ 2.2.x, a safe choice regarding this issue. See blast+ documentation for "
                                     "more information.")

            except IndexError:
                # if the tblastn result file is empty, for example in phase 2
                # if 100% was found in phase 1
                pass
        return

    def _extract_incomplete_buscos_ancestral(self):

        logger.info("Extracting missing and fragmented buscos from the file {}...".format(
            os.path.basename(self.ancestral_file)))

        matched_seqs = []
        busco_ids_retrieved = set()
        with open(self.ancestral_file, "rU") as anc_file:

            for record in SeqIO.parse(anc_file, "fasta"):
                if any(record.id.startswith(b) for b in self.incomplete_buscos):
                    # Remove the ancestral variant identifier ("_1" etc) so it matches all other BUSCO IDs.
                    # The identifier is still present in the "name" and "description" Sequence Record attributes.
                    logger.debug("Found ancestral proteins for {}".format(record.id))
                    record.id = record.id.split("_")[0]
                    busco_ids_retrieved.add(record.id)
                    matched_seqs.append(record)

        unmatched_incomplete_buscos = list(set(self.incomplete_buscos) - set(busco_ids_retrieved))
        if len(unmatched_incomplete_buscos) > 0:
            logger.debug("The BUSCO ID(s) {} were not found in the file {}".format(
                unmatched_incomplete_buscos, os.path.basename(self.ancestral_file)))

        self.query_file = self.rerun_query_file
        with open(self.query_file, "w") as out_file:  # Create new query file for second tblastn run
            SeqIO.write(matched_seqs, out_file, "fasta")

        return

    def _get_all_boundaries(self, locations):
        sorted_locs = sorted(locations, key=lambda x: int(x[0]))
        all_boundaries = [sorted_locs[0]]
        for loc in sorted_locs[1:]:
            overlap, boundary = self._get_overlap(all_boundaries[-1], loc)
            if overlap > 0:
                all_boundaries[-1] = boundary
            else:
                all_boundaries.append(boundary)
        return all_boundaries

    def get_coordinates(self):
        self.coords = self._parse_blast_output()
        if self.ancestral_variants:
            self.coords = self._select_busco_variants()
        self._prune()
        return

    def _get_largest_regions(self, candidate_contigs, coords, busco_group):
        size_lists = []

        for contig in candidate_contigs:
            potential_locations = coords[busco_group][contig]["busco_coords"]
            final_regions = self._get_all_boundaries(potential_locations)

            # Get sum of all potential match sizes for a contig
            size_lists.append(self._sum_all_region_sizes(final_regions))

        return size_lists

    @staticmethod
    def _get_overlap(a, b):
        """
        This function checks whether two regions overlap and returns the length of the overlap region along with the
        boundaries of both regions combined as a [start, stop] list.

        :param a: first region, start and end
        :type a: list
        :param b: second region, start and end
        :type b: list
        :returns: overlap, boundary
        :rtype: int, list
        """
        a_start, a_end = a
        b_start, b_end = b
        overlap = min(a_end, b_end) - max(a_start, b_start)
        if overlap > 0:
            boundary = [min(a_start, b_start), max(a_end, b_end)]
        elif b_start > a_start:
            boundary = b
        else:
            boundary = a
        return max(0, overlap), boundary

    def _parse_blast_output(self):
        """
        Read the Blast output
        """
        coords = defaultdict(lambda: defaultdict(defaultdict))  # dict of busco_id -> contig_id -> {info}
        with open(self.blast_filename, "r") as blast_file:
            for line in blast_file:
                if line.startswith("#"):
                    continue
                else:
                    try:
                        line = line.strip().split()
                        busco_name = line[0]
                        contig_id = line[1]
                        busco_start = int(line[6])
                        busco_end = int(line[7])
                        contig_start = int(line[8])
                        contig_end = int(line[9])
                        blast_eval = float(line[10])
                    except (IndexError, ValueError):
                        continue

                    # for minus-strand genes, invert coordinates for convenience
                    if contig_end < contig_start:
                        contig_end, contig_start = contig_start, contig_end

                    # Add all matches to dictionary. The top matches are selected out later.
                    if contig_id not in coords[busco_name]:
                        coords[busco_name][contig_id] = {"contig_start": contig_start, "contig_end": contig_end,
                                                         "busco_coords": [[busco_start, busco_end]],
                                                         "blast_eval": blast_eval}

                    elif contig_id in coords[busco_name]:  # i.e. if the same gene matched the busco more than once.
                        # now update coordinates
                        coords = self._update_coordinates(coords, busco_name, contig_id, busco_start, busco_end,
                                                          contig_start, contig_end, blast_eval)

        return dict(coords)

    def _select_busco_variants(self):
        """
        Filter contig matches to prevent multiple BUSCO variants matching the same contig.
        The current behaviour combines all contig matches for all BUSCO variants, as long as the contig matches are
        different. There is an open question over whether or not we should only return the contig matches for a single
        BUSCO variant instead of all of them combined. This should only be an issue for the Transcriptome mode.
        :return:
        """
        selected_coords = defaultdict(lambda: defaultdict(defaultdict))
        for busco_name, contigs in self.coords.items():
            busco_basename = busco_name.split("_")[0]
            if busco_basename in selected_coords:
                for contig_id in contigs:
                    if contig_id in selected_coords[busco_basename]:
                        if contigs[contig_id]["blast_eval"] < selected_coords[busco_basename][contig_id]["blast_eval"]:
                            selected_coords[busco_basename][contig_id] = contigs[contig_id]
                    else:
                        selected_coords[busco_basename][contig_id] = contigs[contig_id]
            else:
                selected_coords[busco_basename] = contigs

        return selected_coords

    def _prune(self):
        for busco_name, contigs in self.coords.items():
            if len(contigs) > self.region_limit:
                # Sort by blast eval, then isolate smallest values leaving just "region_limit" number of contigs per
                # busco_name
                contigs_to_remove = sorted(
                    contigs, key=lambda contig: contigs[contig]["blast_eval"])[self.region_limit:]
                for c in contigs_to_remove:
                    self.coords[busco_name].pop(c)
        return

    @staticmethod
    def _sum_all_region_sizes(deck):
        """
        Sum all interval sizes in input list
        :param deck:
        :type deck: list
        :return:
        :rtype: int
        """
        total = 0
        for entry in deck:
            total += entry[1] - entry[0]
        return total

    @staticmethod
    def _update_coordinates(coords, busco_name, contig, busco_start, busco_end, contig_start, contig_end, blast_eval):
        """
        If a contig match starts or ends withing 50 kb of a previous match, extend the recorded start and end positions
        of the contig match, and record the start/end locations of the busco match.
        If the contig match is entirely within a previous match, just record the start/end locations of the busco match.
        If the match is outside 50 kb of a previous match, ignore it. The tblastn output file ranks matches in order of
        bitscore (inverse order of eval) so these subsequent matches at different locations are guaranteed not to be
        better than the ones already recorded for that contig.
        :param coords: # todo: fill in details
        :param busco_name:
        :param contig:
        :param busco_start:
        :param busco_end:
        :param contig_start:
        :param contig_end:
        :param blast_eval:
        :return:
        """
        append_busco_coords = False

        # Check if contig starts before and within 50kb of current position
        if 0 <= coords[busco_name][contig]["contig_start"] - contig_start <= 50000:
            coords[busco_name][contig]["contig_start"] = contig_start
            append_busco_coords = True

        # Check if contig ends after and within 50 kbs of current position
        if 0 <= contig_end - coords[busco_name][contig]["contig_end"] <= 50000:
            coords[busco_name][contig]["contig_end"] = contig_end
            append_busco_coords = True
        # Else, check if contig starts inside current coordinates
        elif coords[busco_name][contig]["contig_end"] >= contig_start >= coords[busco_name][contig]["contig_start"]:
            # If contig ends inside current coordinates, just add alignment positions to list
            if contig_end <= coords[busco_name][contig]["contig_end"]:
                append_busco_coords = True

            # If contig ends after current coordinates, extend contig end
            else:
                coords[busco_name][contig]["contig_end"] = contig_end
                append_busco_coords = True

        # moved to its own "if" statement to avoid multiple appends from the "if" statements above
        if append_busco_coords:
            coords[busco_name][contig]["busco_coords"].append([busco_start, busco_end])

            if blast_eval < coords[busco_name][contig]["blast_eval"]:
                coords[busco_name][contig]["blast_eval"] = blast_eval

        return coords

    def filter_best_matches(self):

        # Get a list of all start and stop positions of possible busco locations, merging overlapping regions
        for busco_group in self.coords:
            candidate_contigs = list(self.coords[busco_group].keys())
            size_lists = self._get_largest_regions(candidate_contigs, self.coords, busco_group)
            max_size = max(size_lists)  # Get largest match size for a busco group
            # Include all location matches for a busco as long as they are within 70% of the maximum size match
            size_cutoff = int(0.7 * max_size)
            for c, contig_name in enumerate(candidate_contigs):
                if size_lists[c] < size_cutoff:
                    self.coords[busco_group].pop(contig_name)
        return

    def write_coordinates_to_file(self):

        with open(self.coords_filename, "w") as out:
            for busco_group, contig_matches in self.coords.items():
                for contig_name in contig_matches:
                    self.coords[busco_group][contig_name]["contig_start"] = \
                        max(int(self.coords[busco_group][contig_name]["contig_start"]) - self.flank, 0)
                    contig_start = self.coords[busco_group][contig_name]["contig_start"]
                    self.coords[busco_group][contig_name]["contig_end"] += self.flank
                    contig_end = int(self.coords[busco_group][contig_name]["contig_end"])
                    out.write("{}\t{}\t{}\t{}\n".format(busco_group, contig_name, contig_start, contig_end))
        return

    def write_contigs(self):
        # Extract all contig identifiers
        contig_names = []
        for contig_info in self.coords.values():
            for contig in contig_info:
                contig_names.append(contig)

        # Write sequences that match contig ids
        with open(self.input_file, "rU") as f:
            for record in SeqIO.parse(f, "fasta"):
                if record.id in list(set(contig_names)):
                    with open(os.path.join(self.output_seqs, "{}.temp".format(record.id)), "w") as out:
                        SeqIO.write(record, out, "fasta")
        return


class MetaeukParsingError(Exception):

    def __init__(self):
        pass


class MetaeukRunner(BaseRunner):

    name = "metaeuk"

    ACCEPTED_PARAMETERS = ["comp-bias-corr", "add-self-matches", "seed-sub-mat", "s", "k", "k-score", "alph-size",
                           "max-seqs", "split", "split-mode", "split-memory-limit", "diag-score", "exact-kmer-matching",
                           "mask", "mask-lower-case", "min-ungapped-score", "spaced-kmer-mode", "spaced-kmer-pattern",
                           "local-tmp", "disk-space-limit", "a", "alignment-mode", "wrapped-scoring", "e", "min-seq-id",
                           "min_aln_len", "seq_id_mode", "alt-ali", "c", "cov-mode", "realign", "max-rejected",
                           "max-accept", "score-bias", "gap-open", "gap-extend", "zdrop", "pca", "pcb", "mask-profile",
                           "e-profile", "wg", "filter-msa", "max-seq-id", "qid", "qsc", "cov", "diff", "num-iterations",
                           "slice-search", "rescore-mode", "allow-deletion", "min-length", "max-length", "max-gaps",
                           "contig-start-mode", "contig-end-mode", "orf-start-mode", "forward-frames", "reverse-frames",
                           "translation-table", "translate", "use-all-table-starts", "id-offset", "add-orf-stop",
                           "search-type", "start-sens", "sens-steps", "metaeuk-eval", "metaeuk-tcov", "min_intron",
                           "max_overlap", "set-gap-open", "set-gap-extend", "overlap", "protein", "target-key",
                           "reverse-fragments", "sub-mat", "db-load-mode", "force-reuse", "remove-tmp-files",
                           "filter-hits", "sort-results", "omit-consensus", "create-lookup", "chain-alignments",
                           "merge-query", "strand", "compressed", "v"]

    def __init__(self):
        super().__init__()
        self._output_folder = os.path.join(self.run_folder, "metaeuk_output")
        self._tmp_folder = os.path.join(self._output_folder, "tmp")
        self._output_basename = os.path.join(self._output_folder, os.path.basename(self.input_file))
        self.ancestral_file = os.path.join(self.lineage_dataset, "ancestral")
        self.ancestral_variants_file = os.path.join(self.lineage_dataset, "ancestral_variants")
        self.refseq_db = os.path.join(self.lineage_dataset, "refseq_db.faa.gz")
        self.max_intron = self.config.get("busco_run", "max_intron")
        self.max_seq_len = self.config.get("busco_run", "max_seq_len")
        self.min_exon_aa = 10
        try:
            self.extra_params = self.config.get("busco_run", "metaeuk_parameters").replace(',', ' ')
        except NoOptionError:
            self.extra_params = ""
        self.param_keys = []
        self.param_values = []
        self.create_dirs(self._output_folder)
        self.gene_details = None

        self.headers_file = "{}.headersMap.tsv".format(self._output_basename)
        self.codon_file = "{}.codon.fas".format(self._output_basename)
        self.pred_protein_seqs = "{}.fas".format(self._output_basename)
        self.pred_protein_seqs_modified = self.pred_protein_seqs.replace(".fas", ".modified.fas")

        self.sequences_aa = {}

        self.init_checkpoint_file()
        self.run_number += 1

    def check_tool_dependencies(self):
        pass

    def configure_job(self, *args):

        metaeuk_job = self.create_job()
        metaeuk_job.add_parameter("easy-predict")
        metaeuk_job.add_parameter("--threads")
        metaeuk_job.add_parameter(str(self.cpus))
        metaeuk_job.add_parameter(self.input_file)
        metaeuk_job.add_parameter(self.refseq_db)
        metaeuk_job.add_parameter(self._output_basename)
        metaeuk_job.add_parameter(self._tmp_folder)
        metaeuk_job.add_parameter("--max-intron")
        metaeuk_job.add_parameter(str(self.max_intron))
        metaeuk_job.add_parameter("--max-seq-len")
        metaeuk_job.add_parameter(str(self.max_seq_len))
        metaeuk_job.add_parameter("--min-exon-aa")
        metaeuk_job.add_parameter(str(self.min_exon_aa))
        for k, key in enumerate(self.param_keys):
            metaeuk_job.add_parameter("--{}".format(key))
            metaeuk_job.add_parameter("{}".format(str(self.param_values[k])))

        return metaeuk_job

    def generate_job_args(self):
        yield

    @property
    def output_folder(self):
        return self._output_folder

    def run(self):
        super().run()
        if self.extra_params:
            logger.info("Additional parameters for Metaeuk are {}: ".format(self.extra_params))
            self.param_keys, self.param_values = self.parse_parameters()

        # self.cwd = self._output_folder
        self.total = 1
        self.run_jobs()

    def get_version(self):
        help_output = subprocess.check_output([self.cmd, "-h"], stderr=subprocess.STDOUT, shell=False)
        lines = help_output.decode("utf-8").split("\n")
        version = None
        for line in lines:
            if line.startswith("metaeuk Version:"):
                version = line.strip().split(" ")[-1]
        return version

    def parse_parameters(self):
        accepted_keys = []
        accepted_values = []
        if self.extra_params:
            self.extra_params = self.extra_params.strip("\" \'")
            try:
                if self.extra_params.startswith("--"):
                    key_val_pairs = self.extra_params.split(" --")
                    for kv in key_val_pairs:
                        key_vals = kv.strip("- ").split("=")
                        if len(key_vals) == 2:
                            key, val = key_vals
                            if key in type(self).ACCEPTED_PARAMETERS:
                                accepted_keys.append(key.strip())
                                accepted_values.append(val.strip())
                            else:
                                logger.warning("{} is not an accepted parameter for Metaeuk.".format(key))
                        else:
                            raise MetaeukParsingError
                else:
                    raise MetaeukParsingError
            except MetaeukParsingError:
                logger.warning(
                    "Metaeuk parameters are not correctly formatted. Please enter them as follows: "
                    "\"--param1=value1 --param2=value2\" etc. Proceeding without additional parameters.")
                return [], []
        return accepted_keys, accepted_values

    @staticmethod
    def parse_header(header):
        header_parts = header.split("|")
        T_acc = header_parts[0]
        C_acc = header_parts[1]
        strand = header_parts[2]
        bitscore = header_parts[3]
        eval = header_parts[4]
        num_exons = header_parts[5]
        low_coord = header_parts[6]
        high_coord = header_parts[7]
        exon_coords = header_parts[8:]

        all_low_exon_coords = []
        all_taken_low_exon_coords = []
        all_high_exon_coords = []
        all_taken_high_exon_coords = []
        all_exon_nucl_len = []
        all_taken_exon_nucl_len = []
        for exon in exon_coords:
            low_exon_coords, high_exon_coords, nucl_lens = exon.split(":")

            low_exon_coord, taken_low_exon_coord = low_exon_coords.split("[")
            all_low_exon_coords.append(low_exon_coord)
            all_taken_low_exon_coords.append(taken_low_exon_coord.strip("]"))

            high_exon_coord, taken_high_exon_coord = high_exon_coords.split("[")
            all_high_exon_coords.append(high_exon_coord)
            all_taken_high_exon_coords.append(taken_high_exon_coord.strip("]"))

            nucl_len, taken_nucl_len = nucl_lens.split("[")
            all_exon_nucl_len.append(nucl_len)
            all_taken_exon_nucl_len.append(taken_nucl_len.strip().rstrip("]"))

        gene_id = "{}:{}-{}".format(C_acc, low_coord, high_coord)

        details = {"T_acc": T_acc, "C_acc": C_acc, "S": strand, "bitscore": bitscore, "e-value": eval,
                   "num_exons": num_exons, "low_coord": low_coord, "high_coord": high_coord,
                   "all_low_exon_coords": all_low_exon_coords, "all_taken_low_exon_coords": all_taken_low_exon_coords,
                   "all_high_exon_coords": all_high_exon_coords, "all_taken_high_exon_coords": all_taken_high_exon_coords,
                   "all_exon_nucl_len": all_exon_nucl_len, "all_taken_exon_nucl_len": all_taken_exon_nucl_len,
                   "gene_id": gene_id}
        return details

    def edit_file_header(self):
        all_records = []
        with open(self.pred_protein_seqs, "rU") as f:
            for record in SeqIO.parse(f, "fasta"):
                header_details = self.parse_header(record.id)
                record.id = header_details["gene_id"]
                record.name = header_details["gene_id"]
                record.description = header_details["gene_id"]
                all_records.append(record)
                self.sequences_aa[record.id] = record

        with open(self.pred_protein_seqs_modified, "w") as f_mod:
            SeqIO.write(all_records, f_mod, "fasta")

    def get_gene_details(self):
        self.gene_details = defaultdict(list)
        with open(self.headers_file, "r") as f:
            lines = f.readlines()

        try:
            for line in lines:
                header = line.split("\t")[-1]
                header_details = self.parse_header(header)

                self.gene_details[header_details["gene_id"]].append({"gene_start": header_details["low_coord"],
                                                                     "gene_end": header_details["high_coord"]})


        except KeyError:
            raise SystemExit("*headersMap.tsv file could not be parsed.")


class SEPPRunner(BaseRunner):

    name = "sepp"

    def __init__(self):
        super().__init__()
        self._output_folder = os.path.join(self.main_out, "auto_lineage", self.lineage_results_dir)
        self.placement_folder = os.path.join(self._output_folder, "placement_files")
        self.datasets_version = self.config.get("busco_run", "datasets_version")

        self.init_checkpoint_file()

    def configure_runner(self, tree_nwk_file, tree_metadata_file, supermatrix_file, downloader):
        self.run_number += 1
        self.tree_nwk_file = tree_nwk_file
        self.tree_metadata_file = tree_metadata_file
        self.supermatrix_file = supermatrix_file
        self.downloader = downloader

    def generate_job_args(self):
        yield

    def run(self):
        super().run()
        self.total = 1
        self.run_jobs()

    def configure_job(self, *args):
        sepp_job = self.create_job()
        sepp_job.add_parameter("--cpu")
        sepp_job.add_parameter(str(self.cpus))
        sepp_job.add_parameter("--outdir")
        sepp_job.add_parameter(self.placement_folder)
        sepp_job.add_parameter("-t")
        sepp_job.add_parameter(self.tree_nwk_file)
        sepp_job.add_parameter("-r")
        sepp_job.add_parameter(self.tree_metadata_file)
        sepp_job.add_parameter("-a")
        sepp_job.add_parameter(self.supermatrix_file)
        sepp_job.add_parameter("-f")
        sepp_job.add_parameter(os.path.join(self.placement_folder, "marker_genes.fasta"))
        sepp_job.add_parameter("-F")
        sepp_job.add_parameter("15")
        sepp_job.add_parameter("-m")
        sepp_job.add_parameter("amino")
        return sepp_job

    def check_tool_dependencies(self):
        pass

    def get_version(self):
        sepp_version = subprocess.check_output([self.cmd, "-v"], stderr=subprocess.STDOUT, shell=False)
        sepp_version = sepp_version.decode("utf-8")
        sepp_version = sepp_version.strip().split(" ")[1]
        return sepp_version

    @property
    def output_folder(self):
        return self._output_folder
