from busco.busco_tools.base import BaseRunner, NoGenesError
import os
import re
from collections import defaultdict
from busco.BuscoLogger import BuscoLogger
from Bio import SeqIO
import shutil
import numpy as np
from configparser import NoOptionError
import subprocess
from busco.Exceptions import BuscoError

logger = BuscoLogger.get_logger(__name__)


class ProdigalRunner(BaseRunner):

    name = "prodigal"
    cmd = "prodigal"

    _gc_run_results = defaultdict(dict)

    def __init__(self):
        super().__init__()
        self._output_folder = os.path.join(self.main_out, "prodigal_output")
        self._pred_genes_dir = os.path.join(self._output_folder, "predicted_genes")
        self._tmp_path = os.path.join(self._pred_genes_dir, "tmp")
        self.cpus = 1
        self.create_dirs([self._pred_genes_dir, self._tmp_path])
        self.opt_gc_file = os.path.join(
            self._output_folder, ".optimal_gc"
        )  # hidden file useful for restart mode

        # Get genetic_code from dataset.cfg file
        # bacteria/archaea=11; Entomoplasmatales,Mycoplasmatales=4
        try:
            self._genetic_code = self.config.get(
                "busco_run", "prodigal_genetic_code"
            ).split(",")
        except NoOptionError:
            self._genetic_code = ["11"]

        # Set the ambiguous coding density range
        try:
            self._cd_upper = (
                float(self.config.get("busco_run", "ambiguous_cd_range_upper"))
                if len(self._genetic_code) > 1
                else 0
            )
        except NoOptionError:
            raise BuscoError(
                "Dataset config file does not contain required information. Please upgrade datasets."
            )

        self.current_gc = None
        self._current_run_mode = None
        self._tmp_name_faa = None
        self._tmp_name_fna = None

        self.output_faa = os.path.join(self._pred_genes_dir, "predicted.faa")
        self._output_fna = os.path.join(self._pred_genes_dir, "predicted.fna")
        self.sequences_aa = {}
        self.sequences_nt = {}
        self.gene_details = {}

        self._input_length = self._get_genome_length()
        self._run_mode = ["single", "meta"] if self._input_length > 100000 else ["meta"]

        self.run_number += 1
        self.good_run_found = None

    def configure_runner(self):
        super().configure_runner()

    @property
    def output_folder(self):
        return self._output_folder

    def check_tool_dependencies(self):
        pass

    def generate_job_args(self):
        yield

    def run(self, restart=False, run_check=False):
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

        number_of_runs = 0
        self.good_run_found = False

        for ix, m in enumerate(self._run_mode):
            self._current_run_mode = m
            for g in self._genetic_code:
                number_of_runs += 1

                self.current_gc = g

                (
                    self._tmp_name_fna,
                    self._tmp_name_faa,
                    self.logfile_path_out,
                    self.logfile_path_err,
                ) = self.get_output_filenames()

                self._gc_run_results[self.current_gc].update(
                    {"tmp_name": self._tmp_name_faa, "log_file": self.logfile_path_out}
                )

                if (
                    os.path.exists(self._tmp_name_faa) and not restart
                ):  # Check to see if has already been run with these parameters
                    coding_density = self._gc_run_results[g]["cd"]
                elif restart:
                    coding_length = self._get_coding_length(self.logfile_path_out)
                    coding_density = coding_length / self._input_length
                    self._gc_run_results[self.current_gc].update({"cd": coding_density})
                else:
                    logger.info(
                        "Running Prodigal with genetic code {} in {} mode".format(
                            self.current_gc, self._current_run_mode
                        )
                    )
                    self.total = 1
                    if not restart:
                        self.run_jobs()
                    coding_length = self._get_coding_length(self.logfile_path_out)
                    coding_density = coding_length / self._input_length
                    self._gc_run_results[self.current_gc].update({"cd": coding_density})

                logger.debug("Coding density is {}".format(coding_density))
                tmp_files.append(self._gc_run_results[self.current_gc]["tmp_name"])

                # if the coding density is above the ambiguous range, then continue with these parameters
                if coding_density >= self._cd_upper:
                    self.good_run_found = True
                    break

            if (
                self.good_run_found and run_check
            ):  # run check is for restart mode to check if an existing run is good enough to skip any outstanding GCs.
                break

            # If output files from both runs in "single" mode are empty, run again in "meta" mode, else raise Exception.
            if not any([os.stat(tmp_file).st_size > 0 for tmp_file in tmp_files]):
                if ix + 1 == len(self._run_mode):
                    raise NoGenesError("Prodigal")
                else:
                    continue

            # if only one genetic code to consider, proceed with it
            # if there is more than one possible set of parameters, decide which to use
            self.current_gc = (
                self._select_best_gc() if number_of_runs > 1 else self._genetic_code[0]
            )

            selected_logfile = self._gc_run_results[self.current_gc]["log_file"]
            selected_tmpfile = self._gc_run_results[self.current_gc]["tmp_name"]

            self._organize_prodigal_files(selected_tmpfile, selected_logfile)
            self.get_gene_details()
            self._gc_run_results[self.current_gc].update(
                {
                    "seqs_aa": self.sequences_aa,
                    "seqs_nt": self.sequences_nt,
                    "gene_details": self.gene_details,
                }
            )
            break
        if not run_check:
            logger.info("Genetic code {} selected as optimal".format(self.current_gc))
        return

    def get_output_filenames(self):
        file_id = os.path.join(
            self._tmp_path,
            "prodigal_mode_{0}_code_{1}".format(
                self._current_run_mode, self.current_gc
            ),
        )
        tmp_name_fna = "{}.fna".format(file_id)
        tmp_name_faa = "{}.faa".format(file_id)
        logfile_path_out = "{}_out.log".format(file_id)
        logfile_path_err = "err".join(
            logfile_path_out.rsplit("out", 1)
        )  # Replace only the last occurence of "out" substring
        return tmp_name_fna, tmp_name_faa, logfile_path_out, logfile_path_err

    def reset(self):
        super().reset()
        type(self)._gc_run_results = defaultdict(dict)

    def write_checkpoint_file(self):
        super().write_checkpoint_file(
            additional_args=[("GC", self.current_gc)]
        )  # The GC saved in the checkpoint file is used on restart to determine if all the GCs required have
        # already been run

    def check_previous_completed_run(self):
        """
        This checks to see if an existing run is sufficient to proceed. For example, an existing run for GC=11 might
        exist, but the current GC list is ["11", "4"]. Before breaking out of restart mode and running GC=4 it is
        important to check the coding density of GC=11 to see if there is any need for GC=4 to run.
        :return:
        """
        completed = super().check_previous_completed_run(additional_args=["GC"])
        if completed:
            higher_priority_runs = []
            for g in self._genetic_code:
                gcs_completed = self.add_args["GC"]
                if g not in gcs_completed:
                    gcs_completed = [
                        x for x in gcs_completed if x in higher_priority_runs
                    ]
                    if len(gcs_completed) > 0:
                        self.check_completed_runs(gcs_completed)
                        return self.good_run_found
                    else:
                        return False
                higher_priority_runs.append(g)
        return completed

    def check_completed_runs(self, gc_list):
        original_gc_list = self._genetic_code
        self._genetic_code = gc_list
        self.run(restart=True, run_check=True)
        self._genetic_code = original_gc_list
        return

    def configure_job(self, *args):

        prodigal_job = self.create_job()
        prodigal_job.add_parameter("-p")
        prodigal_job.add_parameter("%s" % self._current_run_mode)
        prodigal_job.add_parameter("-f")
        prodigal_job.add_parameter("gff")
        prodigal_job.add_parameter("-g")
        prodigal_job.add_parameter("%s" % self.current_gc)
        prodigal_job.add_parameter("-a")
        prodigal_job.add_parameter("%s" % self._tmp_name_faa)
        prodigal_job.add_parameter("-d")
        prodigal_job.add_parameter("%s" % self._tmp_name_fna)
        prodigal_job.add_parameter("-i")
        prodigal_job.add_parameter("%s" % self.input_file)

        return prodigal_job

    def get_gene_details(self):
        self.gene_details = defaultdict(list)
        (
            output_filename_fna,
            output_filename_faa,
            self.logfile_path_out,
            self.logfile_path_err,
        ) = self.get_output_filenames()

        with open(output_filename_fna, "rU") as f:
            for record in SeqIO.parse(f, "fasta"):
                gene_name = record.id
                self.sequences_nt[gene_name] = record
                description_parts = record.description.split()
                gene_start = int(description_parts[2])
                gene_end = int(description_parts[4])
                strand = "+" if int(description_parts[6]) == 1 else "-"
                self.gene_details[gene_name].append(
                    {"gene_start": gene_start, "gene_end": gene_end, "strand": strand}
                )

        with open(output_filename_faa, "rU") as f:
            for record in SeqIO.parse(f, "fasta"):
                self.sequences_aa[record.id] = record

        return

    @staticmethod
    def _get_coding_length(out_logfile):
        total_coding_length = 0
        with open(out_logfile, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    try:
                        start = int(line.split("\t")[3])
                        stop = int(line.split("\t")[4])
                        total_coding_length += stop - start
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
        mean_score = sum(scores) / len(scores)
        return mean_score

    def _organize_prodigal_files(self, tmp_file, tmp_logfile):

        shutil.copy(tmp_file, self.output_faa)
        shutil.copy(tmp_file.rpartition(".faa")[0] + ".fna", self._output_fna)

        # copy selected log files from tmp/ to logs/
        new_logname = os.path.join(self.log_folder, "prodigal_out.log")
        shutil.copy(tmp_logfile, new_logname)
        shutil.copy(
            tmp_logfile.rpartition("_out.log")[0] + "_err.log",
            new_logname.rpartition("_out.log")[0] + "_err.log",
        )
        return

    def _select_best_gc(self):
        gcs, cds = zip(
            *[[gc, self._gc_run_results[gc]["cd"]] for gc in self._genetic_code]
        )
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
        prodigal_version = subprocess.check_output(
            [self.cmd, "-v"], stderr=subprocess.STDOUT, shell=False
        )
        prodigal_version = prodigal_version.decode("utf-8")
        prodigal_version = prodigal_version.split("\n")[1].split(":")[0]
        prodigal_version = prodigal_version.replace("Prodigal V", "")
        return prodigal_version
