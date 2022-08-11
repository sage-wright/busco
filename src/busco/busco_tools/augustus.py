from busco.busco_tools.base import BaseRunner, NoGenesError
import os
import re
from collections import defaultdict
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
import numpy as np
from configparser import NoOptionError
import subprocess
from busco.Exceptions import BuscoError

logger = BuscoLogger.get_logger(__name__)


class AugustusParsingError(Exception):
    def __init__(self):
        pass


class AugustusRunner(BaseRunner):

    ACCEPTED_PARAMETERS = [
        "strand",
        "genemodel",
        "singlestrand",
        "hintsfile",
        "extrinsicCfgFile",
        "maxDNAPieceSize",
        "protein",
        "introns",
        "start",
        "stop",
        "cds",
        "AUGUSTUS_CONFIG_PATH",
        "alternatives-from-evidence",
        "alternatives-from-sampling",
        "sample",
        "minexonintronprob",
        "minmeanexonintronprob",
        "maxtracks",
        "gff3",
        "UTR",
        "outfile",
        "noInFrameStop",
        "noprediction",
        "contentmodels",
        "translation_table",
        "temperature",
        "proteinprofile",
        "progress",
        "predictionStart",
        "predictionEnd",
        "uniqueGeneId",
    ]

    name = "augustus"
    cmd = "augustus"

    def __init__(self):
        self.gene_details = None
        try:
            self._augustus_config_path = os.environ.get("AUGUSTUS_CONFIG_PATH")
            self.config.set(
                "busco_run", "augustus_config_path", self._augustus_config_path
            )
        except TypeError:
            raise BuscoError(
                "AUGUSTUS_CONFIG_PATH environment variable has not been set"
            )

        try:
            self._target_species = self.config.get("busco_run", "augustus_species")
        except KeyError:
            raise SystemExit(
                "Something went wrong. Eukaryota datasets should specify an augustus species."
            )

        super().__init__()
        self._output_folder = os.path.join(self.run_folder, "augustus_output")
        self.tmp_dir = os.path.join(self._output_folder, "tmp")
        self.pred_genes_dir = None
        self.pred_genes_dir_initial = os.path.join(
            self._output_folder, "predicted_genes_initial_run"
        )
        self.pred_genes_dir_rerun = os.path.join(
            self._output_folder, "predicted_genes_rerun"
        )
        self.extracted_prot_dir = os.path.join(
            self._output_folder, "extracted_proteins"
        )
        self.err_logfile = os.path.join(self.log_folder, "augustus_err.log")

        try:
            self.extra_params = self.config.get(
                "busco_run", "augustus_parameters"
            ).replace(",", " ")
        except NoOptionError:
            self.extra_params = ""
        self.chunksize = 10

        self.gff_dir = os.path.join(self._output_folder, "gff")
        self.err_logfiles = []
        self.any_gene_found = False
        self.param_keys = []
        self.param_values = []
        self.output_sequences = []
        self.seqs_path = None
        self.coords = None
        self.sequences_aa = None
        self.sequences_nt = None

        self.create_dirs([self._output_folder, self.extracted_prot_dir, self.gff_dir])

    def configure_runner(self, seqs_path, coords, sequences_aa, sequences_nt):
        super().configure_runner()
        self.run_number += 1

        # Placed here to allow reconfiguration for rerun
        self._target_species = self.config.get("busco_run", "augustus_species")

        self.check_tool_dependencies()
        self.gene_details = defaultdict(list)
        self.output_sequences = []

        self.seqs_path = seqs_path
        self.coords = coords

        self.sequences_aa = sequences_aa
        self.sequences_nt = sequences_nt

        self.pred_genes_dir = (
            self.pred_genes_dir_rerun
            if self.run_number == 2
            else self.pred_genes_dir_initial
        )

        # self.tmp_dir placed here to allow it to be recreated during reconfiguration for rerun
        self.create_dirs([self.pred_genes_dir, self.tmp_dir])

    @property
    def output_folder(self):
        return self._output_folder

    def check_tool_dependencies(self):
        """
        check dependencies on files and folders
        properly configured.
        :raises BuscoError: if Augustus config path is not writable or
        not set at all
        :raises BuscoError: if Augustus config path does not contain
        the needed species
        present
        """
        try:
            augustus_species_dir = os.path.join(self._augustus_config_path, "species")
            if not os.access(augustus_species_dir, os.W_OK):
                raise BuscoError(
                    "Cannot write to Augustus species folder, please make sure you have write "
                    "permissions to {}".format(augustus_species_dir)
                )

        except TypeError:
            raise BuscoError("The environment variable AUGUSTUS_CONFIG_PATH is not set")

        if not os.path.exists(os.path.join(augustus_species_dir, self._target_species)):
            # Exclude the case where this is a restarted run and the retraining parameters have already been moved.
            if (
                self.config.getboolean("busco_run", "restart")
                and self.run_number == 2
                and os.path.exists(
                    os.path.join(
                        self._output_folder,
                        "retraining_parameters",
                        self._target_species,
                    )
                )
            ):
                pass
            else:
                raise BuscoError(
                    'Impossible to locate the species "{0}" in Augustus species folder'
                    " ({1}), check that AUGUSTUS_CONFIG_PATH is properly set"
                    " and contains this species. \n\t\tSee the help if you want "
                    "to provide an alternative species".format(
                        self._target_species, augustus_species_dir
                    )
                )

    @log(
        "Running Augustus prediction using {} as species:",
        logger,
        attr_name="_target_species",
    )
    def run(self):
        super().run()
        if self.extra_params:
            logger.info(
                "Additional parameters for Augustus are {}: ".format(self.extra_params)
            )
            self.param_keys, self.param_values = self.parse_parameters()

        self.total = self._count_jobs()
        self.run_jobs()

    def reset(self):
        super().reset()

    def process_output(self):
        logger.info("Extracting predicted proteins...")
        files = [
            f
            for f in sorted(os.listdir(self.pred_genes_dir))
            if any(busco_id in f for busco_id in self.coords)
        ]
        for filename in files:
            self._extract_genes_from_augustus_output(filename)

        if not self.any_gene_found and self.run_number == 1:
            raise NoGenesError("Augustus")

        self.gene_details = dict(self.gene_details)

        self._merge_stderr_logs()
        self._remove_individual_err_logs()

        return

    def _count_jobs(self):
        n = 0
        for busco_group, contigs in self.coords.items():
            for _ in contigs:
                n += 1
        return n

    def sort_jobs(self):
        jobs_size_info = []
        for busco_group, contigs in self.coords.items():

            for contig_name, contig_info in contigs.items():
                contig_start = contig_info["contig_start"]
                contig_end = contig_info["contig_end"]
                pred_size = int(contig_end) - int(contig_start)
                jobs_size_info.append(
                    {
                        "busco_group": busco_group,
                        "contig_name": contig_name,
                        "contig_start": contig_start,
                        "contig_end": contig_end,
                        "pred_size": pred_size,
                    }
                )
        job_sizes = [item["pred_size"] for item in jobs_size_info]
        new_job_order = np.argsort(job_sizes)[::-1]
        ordered_jobs = [jobs_size_info[i] for i in new_job_order]
        return ordered_jobs

    def generate_job_args(self):
        contig_ordinal_inds = defaultdict(int)
        njobs = 0

        ordered_jobs = self.sort_jobs()

        for job_info in ordered_jobs:
            contig_name = job_info["contig_name"]
            busco_group = job_info["busco_group"]
            contig_start = job_info["contig_start"]
            contig_end = job_info["contig_end"]
            contig_tmp_file = "{}.temp".format(
                contig_name[:100]
            )  # Avoid very long filenames
            contig_ordinal_inds[busco_group] += 1
            output_index = contig_ordinal_inds[busco_group]
            out_filename = os.path.join(
                self.pred_genes_dir, "{}.out.{}".format(busco_group, output_index)
            )
            njobs += 1

            yield busco_group, contig_tmp_file, contig_start, contig_end, out_filename

    @log(
        "Additional parameters for Augustus are {}: ",
        logger,
        attr_name="_target_species",
    )
    def parse_parameters(self):
        accepted_keys = []
        accepted_values = []
        if self.extra_params:
            self.extra_params = self.extra_params.strip("\" '")
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
                                logger.warning(
                                    "{} is not an accepted parameter for Augustus.".format(
                                        key
                                    )
                                )
                        else:
                            raise AugustusParsingError
                else:
                    raise AugustusParsingError
            except AugustusParsingError:
                logger.warning(
                    "Augustus parameters are not correctly formatted. Please enter them as follows: "
                    '"--param1=value1 --param2=value2" etc. Proceeding without additional parameters.'
                )
                return [], []
        return accepted_keys, accepted_values

    def _merge_stderr_logs(self):
        with open(self.err_logfile, "a") as f:
            for err_logfile in self.err_logfiles:
                with open(err_logfile, "r") as g:
                    content = g.readlines()
                    f.writelines(content)
        return

    def _remove_individual_err_logs(self):
        shutil.rmtree(self.tmp_dir)
        return

    def get_version(self):
        augustus_help_output = subprocess.check_output(
            [self.cmd, "--version"], stderr=subprocess.STDOUT, shell=False
        )
        augustus_help_output = augustus_help_output.decode("utf-8")
        s = augustus_help_output.split("\n")[0]
        augustus_version = s[s.find("(") + 1 : s.find(")")]
        return augustus_version

    def configure_job(
        self, busco_group, contig_tmp_file, contig_start, contig_end, out_filename
    ):
        # Augustus does not provide an option to write to an output file, so have to change the pipe target from the
        # log file to the desired output file
        self.logfile_path_out = out_filename
        err_logfile = os.path.join(
            self.tmp_dir, os.path.basename(out_filename.rpartition(".out")[0] + ".err")
        )
        self.logfile_path_err = err_logfile
        self.err_logfiles.append(err_logfile)

        augustus_job = self.create_job()
        augustus_job.add_parameter("--codingseq=1")
        augustus_job.add_parameter(
            "--proteinprofile={}".format(
                os.path.join(
                    self.lineage_dataset, "prfl", "{}.prfl".format(busco_group)
                )
            )
        )
        augustus_job.add_parameter("--predictionStart={}".format(contig_start))
        augustus_job.add_parameter("--predictionEnd={}".format(contig_end))
        augustus_job.add_parameter("--species={}".format(self._target_species))
        for k, key in enumerate(self.param_keys):
            augustus_job.add_parameter("--{}={}".format(key, self.param_values[k]))
        augustus_job.add_parameter(os.path.join(self.seqs_path, contig_tmp_file))
        return augustus_job

    def _extract_genes_from_augustus_output(self, filename):
        # todo: consider parallelizing this and other parsing functions

        gene_id = None
        gene_info = []
        sequences_aa = []
        sequences_nt = []
        gene_found = False
        completed_record = False

        with open(
            os.path.join(self.pred_genes_dir, filename), "r", encoding="utf-8"
        ) as f:
            # utf-8 encoding needed to handle the umlaut in the third line of the file.
            gene_info_section = False
            nt_sequence_section = False
            aa_sequence_section = False
            nt_sequence_parts = []
            aa_sequence_parts = []

            for line in f:

                if aa_sequence_section:
                    if "]" in line:
                        line = line.strip().lstrip("# ").rstrip("]")
                        aa_sequence_parts.append(line)
                        aa_sequence_section = False
                        completed_record = True
                        if gene_id is not None:
                            aa_sequence = "".join(aa_sequence_parts)
                            nt_sequence = "".join(nt_sequence_parts)
                            seq_record_aa = SeqRecord(
                                Seq(aa_sequence.upper()), id=gene_id
                            )
                            seq_record_nt = SeqRecord(
                                Seq(nt_sequence.upper()), id=gene_id
                            )
                            sequences_aa.append(seq_record_aa)
                            sequences_nt.append(seq_record_nt)
                            aa_sequence_parts = []
                            nt_sequence_parts = []
                            gene_id = None
                        continue

                    else:
                        line = line.strip().lstrip("# ").rstrip("]")
                        aa_sequence_parts.append(line)
                        continue

                if line.startswith("# protein"):
                    nt_sequence_section = False
                    aa_sequence_section = True
                    if "]" in line:
                        line = line.strip().rstrip("]").split("[")
                        aa_sequence_parts.append(line[1])
                        aa_sequence_section = False
                        completed_record = True
                        if gene_id is not None:
                            aa_sequence = "".join(aa_sequence_parts)
                            nt_sequence = "".join(nt_sequence_parts)
                            seq_record_aa = SeqRecord(
                                Seq(aa_sequence.upper()), id=gene_id
                            )
                            seq_record_nt = SeqRecord(
                                Seq(nt_sequence.upper()), id=gene_id
                            )
                            sequences_aa.append(seq_record_aa)
                            sequences_nt.append(seq_record_nt)
                            aa_sequence_parts = []
                            nt_sequence_parts = []
                            gene_id = None
                    else:
                        line = line.strip().rstrip("]").split("[")
                        aa_sequence_parts.append(line[1])
                    continue

                if nt_sequence_section:
                    line = line.strip().lstrip("# ").rstrip("]")
                    nt_sequence_parts.append(line)
                    continue

                if line.startswith("# coding sequence"):
                    gene_info = []
                    gene_info_section = False
                    nt_sequence_section = True
                    line = (
                        line.strip().rstrip("]").split("[")
                    )  # Extract sequence part of line
                    nt_sequence_parts.append(line[1])
                    continue

                if gene_info_section:
                    line = line.strip().split()
                    seq_name = line[0].strip()
                    gene_start = line[3].strip()
                    gene_end = line[4].strip()
                    strand = line[6].strip()
                    if not gene_id:
                        gene_id = "{}:{}-{}".format(seq_name, gene_start, gene_end)
                        self.gene_details[gene_id].append(
                            {
                                "gene_start": gene_start,
                                "gene_end": gene_end,
                                "strand": strand,
                            }
                        )
                    gene_info.append("\t".join(line))
                    continue

                if line.startswith("# start gene"):
                    gene_found = True
                    self.any_gene_found = True
                    gene_info_section = True
                    completed_record = False
                    continue

            if gene_found and not completed_record:
                logger.warning("Augustus output file {} truncated".format(filename))

        self.sequences_aa.update({record.id: record for record in sequences_aa})
        self.sequences_nt.update({record.id: record for record in sequences_nt})
        if gene_found:
            self._write_sequences_to_file(filename, sequences_nt, sequences_aa)

        return

    def make_gff_files(self, single_copy_buscos):

        for b in single_copy_buscos:
            gene_info = []
            busco_files = [
                f for f in os.listdir(self.pred_genes_dir) if f.startswith(b)
            ]
            pred_genes_dir_current = self.pred_genes_dir
            if (
                len(busco_files) == 0
                and pred_genes_dir_current == self.pred_genes_dir_rerun
            ):
                pred_genes_dir_current = self.pred_genes_dir_initial
                busco_files = [
                    f for f in os.listdir(pred_genes_dir_current) if f.startswith(b)
                ]
            gff_filename = os.path.join(self.gff_dir, "{}.gff".format(b))
            single_copy_busco_gene = list(single_copy_buscos[b].keys())[0]
            gene_id_parts = single_copy_busco_gene.split(":")
            if (
                len(gene_id_parts) > 2
            ):  # if a ":" is present in the gene id, we don't want to break it up
                gene_id_parts = [":".join(gene_id_parts[:-1]), gene_id_parts[-1]]
            single_copy_busco_gene_id = gene_id_parts[0]
            (
                single_copy_busco_gene_start_coord,
                single_copy_busco_gene_end_coord,
            ) = gene_id_parts[1].split("-")
            gene_found = False
            for filename in busco_files:
                match_number = filename.split(".")[-1]
                with open(
                    os.path.join(pred_genes_dir_current, filename),
                    "r",
                    encoding="utf-8",
                ) as f:
                    gene_info_section = False
                    for line in f:
                        if gene_info_section and line.startswith("# coding sequence"):
                            with open(gff_filename, "w") as g:
                                g.write("\n".join(gene_info) + "\n")
                            gene_info = []
                            break

                        if line.startswith("# start gene"):
                            gene_info_section = True
                            continue

                        if gene_info_section:
                            line = line.strip().split()
                            seq_name = line[0]
                            gene_start = line[3]
                            gene_end = line[4]
                            if gene_found or (
                                seq_name == single_copy_busco_gene_id
                                and gene_start == single_copy_busco_gene_start_coord
                                and gene_end == single_copy_busco_gene_end_coord
                            ):
                                gene_found = True
                                gene_id_info = line[-1]
                                line[-1] = self.edit_gene_identifier(
                                    gene_id_info, match_number
                                )
                                if len(line) == 12:
                                    gene_id_info_2 = line[-3]
                                    line[-3] = self.edit_gene_identifier(
                                        gene_id_info_2, match_number
                                    )
                                gene_info.append("\t".join(line))
                            else:
                                gene_info_section = False
                            continue
                if gene_found:
                    break
            if not gene_found and self.run_number == 1:
                raise BuscoError(
                    "Unable to find single copy BUSCO gene in Augustus output."
                )

        return

    def edit_gene_identifier(self, orig_str, match_num):
        modified_str = re.sub(
            r"g([0-9])", r"r{}.m{}.g\1".format(self.run_number, match_num), orig_str
        )
        return modified_str

    def _write_sequences_to_file(self, filename, sequences_nt, sequences_aa):

        filename_parts = filename.rpartition(".out")
        output_fna = os.path.join(
            self.extracted_prot_dir, filename_parts[0] + ".fna" + filename_parts[-1]
        )
        output_faa = os.path.join(
            self.extracted_prot_dir, filename_parts[0] + ".faa" + filename_parts[-1]
        )
        self.output_sequences.append(output_faa)

        with open(output_fna, "w") as out_fna:
            SeqIO.write(sequences_nt, out_fna, "fasta")
        with open(output_faa, "w") as out_faa:
            SeqIO.write(sequences_aa, out_faa, "fasta")

        return

    def move_retraining_parameters(self):
        """
        This function moves retraining parameters from augustus species folder
        to the run folder
        """
        if self._target_species.startswith("BUSCO"):
            augustus_species_path = os.path.join(
                self._augustus_config_path, "species", self._target_species
            )
            if os.path.exists(augustus_species_path):
                new_path = os.path.join(
                    self._output_folder, "retraining_parameters", self._target_species
                )
                shutil.move(augustus_species_path, new_path)
            elif self.config.getboolean("busco_run", "restart") and os.path.exists(
                os.path.join(
                    self._output_folder, "retraining_parameters", self._target_species
                )
            ):
                pass
            else:
                logger.warning("Augustus did not produce a retrained species folder.")
        return


class GFF2GBRunner(BaseRunner):

    name = "gff2gbSmallDNA.pl"
    cmd = "gff2gbSmallDNA.pl"

    def __init__(self):
        super().__init__()
        self._output_folder = os.path.join(self.run_folder, "augustus_output")
        self.gff_folder = os.path.join(self._output_folder, "gff")
        self.gb_folder = os.path.join(self._output_folder, "gb")
        self.create_dirs([self.gff_folder, self.gb_folder])
        self.single_copy_buscos = None

    def configure_runner(self, single_copy_buscos):
        super().configure_runner()
        self.run_number += 1
        self.single_copy_buscos = single_copy_buscos

    def run(self):
        super().run()
        self.total = self._count_jobs()
        self.run_jobs()

    def reset(self):
        super().reset()

    def _count_jobs(self):
        n = len(self.single_copy_buscos)
        return n

    def generate_job_args(self):
        for busco_id in self.single_copy_buscos:
            yield busco_id

    def configure_job(self, busco_id):
        gff2_gb_small_dna_pl_job = self.create_job()
        gff2_gb_small_dna_pl_job.add_parameter(
            os.path.join(self.gff_folder, "{}.gff".format(busco_id))
        )
        gff2_gb_small_dna_pl_job.add_parameter(self.input_file)
        gff2_gb_small_dna_pl_job.add_parameter("1000")
        gff2_gb_small_dna_pl_job.add_parameter(
            os.path.join(self.gb_folder, "{}.raw.gb".format(busco_id))
        )
        return gff2_gb_small_dna_pl_job

    def check_tool_dependencies(self):
        pass

    def get_version(self):
        return

    @property
    def output_folder(self):
        return self._output_folder


class NewSpeciesRunner(BaseRunner):

    name = "new_species.pl"
    cmd = "new_species.pl"

    def __init__(self):
        super().__init__()
        self._output_folder = os.path.join(self.run_folder, "augustus_output")
        self.new_species_name = "BUSCO_{}".format(os.path.basename(self.main_out))
        self.run_number += 1

    def run(self):
        super().run()
        self.total = 1
        self.run_jobs()

    def configure_runner(self, *args):
        super().configure_runner(args)

    def reset(self):
        super().reset()

    def configure_job(self, *args):

        new_species_pl_job = self.create_job()
        # bacteria clade needs to be flagged as "prokaryotic"
        if self.domain == "prokaryota":
            new_species_pl_job.add_parameter("--prokaryotic")
        new_species_pl_job.add_parameter(
            "--species={}".format(os.path.basename(self.new_species_name))
        )
        return new_species_pl_job

    def check_tool_dependencies(self):
        pass

    def generate_job_args(self):
        yield

    def get_version(self):
        return

    @property
    def output_folder(self):
        return self._output_folder


class ETrainingRunner(BaseRunner):

    name = "etraining"
    cmd = "etraining"

    def __init__(self):
        super().__init__()
        self._output_folder = os.path.join(self.run_folder, "augustus_output")
        self._gb_folder = os.path.join(self._output_folder, "gb")
        self.augustus_config_path = self.config.get("busco_run", "augustus_config_path")
        self._training_file = os.path.join(self._output_folder, "training_set.db")
        self.new_species_name = None

    def configure_runner(self, new_species_name):
        super().configure_runner()
        self.run_number += 1
        self.new_species_name = new_species_name
        self._merge_gb_files()

    def run(self):
        super().run()
        self.total = 1
        self.run_jobs()
        self._validate_run()

    def reset(self):
        super().reset()

    def check_tool_dependencies(self):
        pass

    def generate_job_args(self):
        yield

    def _merge_gb_files(self):
        """Concatenate all GB files into one large file"""
        with open(self._training_file, "w") as outfile:
            for fname in os.listdir(self._gb_folder):
                with open(os.path.join(self._gb_folder, fname), "r") as infile:
                    outfile.writelines(infile.readlines())
        return

    def _validate_run(self):
        species_filepath = os.path.join(
            self.augustus_config_path, "species", self.new_species_name
        )
        if os.path.exists(species_filepath) and any(
            "exon_probs" in f for f in os.listdir(species_filepath)
        ):
            return
        else:
            raise BuscoError(
                "Retraining did not complete correctly. Check your Augustus config path environment variable."
            )

    def configure_job(self, *args):
        etraining_job = self.create_job()
        etraining_job.add_parameter("--species={}".format(self.new_species_name))
        etraining_job.add_parameter(
            os.path.join(self.run_folder, "augustus_output", "training_set.db")
        )
        return etraining_job

    def get_version(self):
        return

    @property
    def output_folder(self):
        return self._output_folder


class OptimizeAugustusRunner(BaseRunner):

    name = "optimize_augustus.pl"
    cmd = "optimize_augustus.pl"

    def __init__(self):
        super().__init__()
        self._output_folder = None
        self.training_set_db = None
        self.new_species_name = None

    def configure_runner(self, output_folder, new_species_name):
        self._output_folder = output_folder
        super().configure_runner()
        self.run_number += 1
        self.training_set_db = os.path.join(self._output_folder, "training_set.db")
        self.new_species_name = new_species_name

    def configure_job(self, *args):
        optimize_augustus_pl_job = self.create_job()
        optimize_augustus_pl_job.add_parameter("--cpus={}".format(self.cpus))
        optimize_augustus_pl_job.add_parameter(
            "--species={}".format(self.new_species_name)
        )
        optimize_augustus_pl_job.add_parameter(self.training_set_db)
        return optimize_augustus_pl_job

    def run(self):
        super().run()
        self.total = 1
        self.run_jobs()

    def reset(self):
        super().reset()

    def generate_job_args(self):
        yield

    def check_tool_dependencies(self):
        pass

    def get_version(self):
        return

    @property
    def output_folder(self):
        return self._output_folder
