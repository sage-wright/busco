# coding: utf-8
"""
metaeuk.py

Module for running Metaeuk.

Author(s): Matthew Berkeley

Copyright (c) 2015-2024, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

from busco.busco_tools.base import GenePredictor, NoRerunFile
import os
from busco.BuscoLogger import BuscoLogger
from Bio import SeqIO
import shutil
from configparser import NoOptionError
import subprocess
import numpy as np
import re
from busco.Exceptions import BuscoError
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = BuscoLogger.get_logger(__name__)


class MetaeukParsingError(Exception):
    def __init__(self):
        pass


class MetaeukRunner(GenePredictor):

    name = "metaeuk"
    cmd = "metaeuk"

    ACCEPTED_PARAMETERS = [
        "comp-bias-corr",
        "add-self-matches",
        "seed-sub-mat",
        "s",
        "k",
        "k-score",
        "alph-size",
        "max-seqs",
        "split",
        "split-mode",
        "split-memory-limit",
        "diag-score",
        "exact-kmer-matching",
        "mask",
        "mask-lower-case",
        "min-ungapped-score",
        "spaced-kmer-mode",
        "spaced-kmer-pattern",
        "local-tmp",
        "disk-space-limit",
        "a",
        "alignment-mode",
        "wrapped-scoring",
        "e",
        "min-seq-id",
        "min-aln-len",
        "seq-id-mode",
        "alt-ali",
        "c",
        "cov-mode",
        "realign",
        "max-rejected",
        "max-accept",
        "score-bias",
        "gap-open",
        "gap-extend",
        "zdrop",
        "pca",
        "pcb",
        "mask-profile",
        "e-profile",
        "wg",
        "filter-msa",
        "max-seq-id",
        "qid",
        "qsc",
        "cov",
        "diff",
        "num-iterations",
        "slice-search",
        "rescore-mode",
        "allow-deletion",
        "min-length",
        "max-length",
        "max-gaps",
        "contig-start-mode",
        "contig-end-mode",
        "orf-start-mode",
        "forward-frames",
        "reverse-frames",
        "translation-table",
        "translate",
        "use-all-table-starts",
        "id-offset",
        "add-orf-stop",
        "search-type",
        "start-sens",
        "sens-steps",
        "metaeuk-eval",
        "metaeuk-tcov",
        "min-intron",
        "min-exon-aa",
        "max-overlap",
        "set-gap-open",
        "set-gap-extend",
        "overlap",
        "protein",
        "target-key",
        "reverse-fragments",
        "sub-mat",
        "db-load-mode",
        "force-reuse",
        "remove-tmp-files",
        "filter-hits",
        "sort-results",
        "omit-consensus",
        "create-lookup",
        "chain-alignments",
        "merge-query",
        "strand",
        "compressed",
        "v",
        "max-intron",
        "max-seq-len",
    ]

    PRESET_PARAMETERS = [
        "--max-intron",
        "--max-seq-len",
        "--min-exon-aa",
        "--max-overlap",
        "--min-intron",
        "--overlap",
    ]

    def __init__(self):
        super().__init__()
        self._output_folder = os.path.join(self.run_folder, "metaeuk_output")
        self._initial_results_folder = os.path.join(
            self._output_folder, "initial_results"
        )
        self._rerun_results_folder = os.path.join(self._output_folder, "rerun_results")
        self._tmp_folder = os.path.join(self._output_folder, "tmp")

        self.ancestral_file = os.path.join(self.lineage_dataset, "ancestral")
        self.ancestral_variants_file = os.path.join(
            self.lineage_dataset, "ancestral_variants"
        )
        try:
            self.max_intron = self.config.get("busco_run", "max_intron")
            self.max_seq_len = self.config.get("busco_run", "max_seq_len")
        except NoOptionError:
            raise BuscoError(
                "{} is an old dataset version and is not compatible with Metaeuk. Please update by using "
                "the '--update-data' command line option".format(self.lineage_dataset)
            )
        self.overlap = 1
        self.s_set = False

        self.extra_params = None
        self.param_keys = []
        self.param_values = []
        self.create_dirs(
            [
                self._output_folder,
                self._initial_results_folder,
                self._rerun_results_folder,
            ]
        )

        self.headers_file = None
        self.codon_file = None
        self.pred_protein_seqs = None
        self.pred_protein_seqs_modified = None
        self.incomplete_buscos = None
        self.sequences_aa = {}
        self.sequences_nt = {}
        self.saved_gene_details = {}

        self.pred_protein_files = []
        self.pred_protein_mod_files = []
        self.codon_files = []
        self.codon_file_modified = None
        self.codon_mod_files = []
        self.headers_files = []
        self.gff_file = None
        self.gff_files = []
        self.combined_pred_protein_seqs = os.path.join(
            self._output_folder, "combined_pred_proteins.fas"
        )
        self.combined_nucleotides_seqs = os.path.join(
            self._output_folder, "combined_nucl_seqs.fas"
        )
        self.refseq_db_rerun = None
        self._output_basename = None
        self.refseq_db = None
        self.min_exon_aa = None
        self.max_overlap = None
        self.min_intron = None

    def configure_runner(self, incomplete_buscos=None):
        super().configure_runner()
        self.saved_gene_details = self.gene_details.copy()
        self.gene_details = defaultdict(dict)
        self.run_number += 1
        self.incomplete_buscos = incomplete_buscos

        if self.run_number > 1:
            self._output_basename = os.path.join(
                self._rerun_results_folder, os.path.basename(self.input_file)
            )
            self.refseq_db_rerun = os.path.join(
                self._output_folder, "refseq_db_rerun.faa"
            )
            self.min_exon_aa = 5
            self.max_overlap = 5
            self.min_intron = 1
        else:
            self._output_basename = os.path.join(
                self._initial_results_folder, os.path.basename(self.input_file)
            )
            gzip_refseq = os.path.join(self.lineage_dataset, "refseq_db.faa.gz")
            self.refseq_db = self.decompress_refseq_file(gzip_refseq)
            self.min_exon_aa = 15
            self.max_overlap = 15
            self.min_intron = 5

        try:
            if self.run_number == 1:
                self.extra_params = self.config.get(
                    "busco_run", "metaeuk_parameters"
                ).replace(",", " ")
            else:
                self.extra_params = self.config.get(
                    "busco_run", "metaeuk_rerun_parameters"
                ).replace(",", " ")
        except NoOptionError:
            self.extra_params = ""

        self.headers_file = "{}.headersMap.tsv".format(self._output_basename)
        self.headers_files.append(self.headers_file)
        self.gff_file = "{}.gff".format(self._output_basename)
        self.gff_files.append(self.gff_file)
        self.codon_file = "{}.codon.fas".format(self._output_basename)
        self.codon_files.append(self.codon_file)
        self.codon_file_modified = ".modified.codon.fas".join(
            self.codon_file.rsplit(".codon.fas", 1)
        )
        self.pred_protein_seqs = "{}.fas".format(self._output_basename)
        self.pred_protein_files.append(self.pred_protein_seqs)
        self.pred_protein_seqs_modified = ".modified.fas".join(
            self.pred_protein_seqs.rsplit(".fas", 1)
        )
        self.pred_protein_mod_files.append(self.pred_protein_seqs_modified)
        self.codon_mod_files.append(self.codon_file_modified)

    def combine_gene_details(self):
        for gene_id, details in self.saved_gene_details.items():
            if gene_id in self.gene_details:
                pass
            else:
                self.gene_details[gene_id] = details
        return

    def combine_run_results(self):
        with open(self.combined_pred_protein_seqs, "w") as combined_output:
            for run_result in self.pred_protein_mod_files:
                with open(run_result, "r") as f:
                    shutil.copyfileobj(f, combined_output)

        with open(self.combined_nucleotides_seqs, "w") as combined_output_nucl:
            for run_result_codons in self.codon_mod_files:
                with open(run_result_codons, "r") as f:
                    shutil.copyfileobj(f, combined_output_nucl)

        return

    @staticmethod
    def select_higher_bitscore_ind(matches):
        bitscores = []
        pos_pattern = "|+|"
        neg_pattern = "|-|"
        for match in matches:
            strand_pattern = pos_pattern if pos_pattern in match else neg_pattern
            score = re.search(
                r"{}(.*?)\|".format(re.escape(strand_pattern)), match
            ).group(1)
            bitscores.append(int(score))
        max_ind = bitscores.index(max(bitscores))
        return max_ind

    def extract_exon_coords(self, match):
        parts = str(match).split("\t")
        header = parts[5]
        details = self.parse_header(header)
        return (
            details["all_taken_low_exon_coords"],
            details["all_taken_high_exon_coords"],
        )

    def check_tool_dependencies(self):
        pass

    def configure_job(self, *args):
        self.logfile_path_out = self.logfile_path_out.replace(
            "metaeuk_out.log", "metaeuk_run{}_out.log".format(self.run_number)
        )
        self.logfile_path_err = self.logfile_path_err.replace(
            "metaeuk_err.log", "metaeuk_run{}_err.log".format(self.run_number)
        )

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
        metaeuk_job.add_parameter("--max-overlap")
        metaeuk_job.add_parameter(str(self.max_overlap))
        metaeuk_job.add_parameter("--min-intron")
        metaeuk_job.add_parameter(str(self.min_intron))
        metaeuk_job.add_parameter("--overlap")
        metaeuk_job.add_parameter(str(self.overlap))
        if self.run_number > 1 and not self.s_set:
            metaeuk_job.add_parameter("-s")
            metaeuk_job.add_parameter("6")
        elif self.run_number == 1 and not self.s_set:
            metaeuk_job.add_parameter("-s")
            metaeuk_job.add_parameter("4.5")
        for k, key in enumerate(self.param_keys):
            dashes = "-" if len(key) == 1 else "--"
            metaeuk_job.add_parameter("{}{}".format(dashes, key))
            metaeuk_job.add_parameter("{}".format(str(self.param_values[k])))
        return metaeuk_job

    def generate_job_args(self):
        yield

    @property
    def output_folder(self):
        return self._output_folder

    def remove_tmp_files(self):
        shutil.rmtree(self._tmp_folder)

    def run(self):
        super().run()
        if self.run_number > 1:
            self._extract_incomplete_buscos_ancestral()
            self.refseq_db = self.refseq_db_rerun
        if self.extra_params:
            logger.info(
                "Additional parameters for Metaeuk are {}: ".format(self.extra_params)
            )
            self.param_keys, self.param_values = self.parse_parameters()
        else:
            self.param_keys = self.param_values = []

        self.total = 1
        self.run_jobs()

    def _extract_incomplete_buscos_ancestral(self):

        logger.info(
            "Extracting missing and fragmented buscos from the file {}...".format(
                os.path.basename(self.refseq_db)
            )
        )

        matched_seqs = []
        busco_ids_retrieved = set()
        with open(self.refseq_db, "r") as refseq_file:

            for record in SeqIO.parse(refseq_file, "fasta"):
                if any(record.id.startswith(b) for b in self.incomplete_buscos):
                    # Remove the ancestral variant identifier ("_1" etc) so it matches all other BUSCO IDs.
                    # The identifier is kept in the Sequence Record ID.
                    busco_ids_retrieved.add(record.id.split("_")[0])
                    matched_seqs.append(record)

        unmatched_incomplete_buscos = list(
            set(self.incomplete_buscos) - set(busco_ids_retrieved)
        )
        if len(unmatched_incomplete_buscos) > 0:
            logger.debug(
                "The BUSCO ID(s) {} were not found in the file {}".format(
                    unmatched_incomplete_buscos, os.path.basename(self.refseq_db)
                )
            )

        with open(
            self.refseq_db_rerun, "w"
        ) as out_file:  # Create new query file for second tblastn run
            SeqIO.write(matched_seqs, out_file, "fasta")

        return

    def get_version(self):
        help_output = subprocess.check_output(
            [self.cmd, "-h"], stderr=subprocess.STDOUT, shell=False
        )
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
            self.extra_params = self.extra_params.strip("\" '")
            try:
                if self.extra_params.startswith("-"):
                    key_val_pairs = self.extra_params.split(" -")
                    for kv in key_val_pairs:
                        key_vals = kv.strip("- ").split("=")
                        if len(key_vals) == 2:
                            key, val = key_vals
                            if key in type(self).ACCEPTED_PARAMETERS:
                                if key == "min-exon-aa":
                                    self.min_exon_aa = val.strip()
                                    continue
                                elif key == "max-intron":
                                    self.max_intron = val.strip()
                                    continue
                                elif key == "max-seq-len":
                                    self.max_seq_len = val.strip()
                                    continue
                                elif key == "max-overlap":
                                    self.max_overlap = val.strip()
                                    continue
                                elif key == "min-intron":
                                    self.min_intron = val.strip()
                                    continue
                                elif key == "overlap":
                                    self.overlap = val.strip()
                                    continue
                                elif key == "s":
                                    self.s_set = True
                                accepted_keys.append(key.strip())
                                accepted_values.append(val.strip())
                            else:
                                logger.warning(
                                    "{} is not an accepted parameter for Metaeuk.".format(
                                        key
                                    )
                                )
                        else:
                            raise MetaeukParsingError
                else:
                    raise MetaeukParsingError
            except MetaeukParsingError:
                logger.warning(
                    "Metaeuk parameters are not correctly formatted. Please enter them as follows: "
                    '"--param1=value1 --param2=value2" etc. Proceeding without additional parameters.'
                )
                return [], []
        return accepted_keys, accepted_values

    @staticmethod
    def parse_header(header):
        header_parts = header.split("|")
        if not header_parts[2] in [
            "+",
            "-",
        ]:  # Deal with sequence IDs that contain the symbol "|"
            try:
                strand_ind = header_parts.index("+")
            except ValueError:
                strand_ind = header_parts.index("-")
            header_parts[1] = "|".join(header_parts[1:strand_ind])
            for i in range(2, strand_ind):
                header_parts.pop(i)

        T_acc = header_parts[0]
        C_acc = header_parts[1]
        strand = header_parts[2]
        bitscore = float(header_parts[3])
        evalue = float(header_parts[4])
        num_exons = int(header_parts[5])
        low_coord = int(header_parts[6])
        high_coord = int(header_parts[7])
        exon_details = header_parts[8:]

        gene_id = "{}|{}:{}-{}".format(T_acc, C_acc, low_coord, high_coord)

        all_low_exon_coords = []
        all_taken_low_exon_coords = []
        all_high_exon_coords = []
        all_taken_high_exon_coords = []
        all_exon_nucl_len = []
        all_taken_exon_nucl_len = []
        recorded_exon_coords = []
        for exon in exon_details:

            low_exon_coords, high_exon_coords, nucl_lens = exon.split(":")

            low_exon_coord, taken_low_exon_coord = low_exon_coords.split("[")
            all_low_exon_coords.append(int(low_exon_coord))
            taken_low_exon_coord = int(taken_low_exon_coord.strip("]"))

            high_exon_coord, taken_high_exon_coord = high_exon_coords.split("[")
            all_high_exon_coords.append(int(high_exon_coord))
            taken_high_exon_coord = int(taken_high_exon_coord.strip("]"))

            nucl_len, taken_nucl_len = nucl_lens.split("[")
            all_exon_nucl_len.append(int(nucl_len))
            taken_nucl_len = int(taken_nucl_len.strip().rstrip("]"))

            recorded_exon_coords.append(
                (taken_low_exon_coord, taken_high_exon_coord, strand)
            )

            # Need to fix the metaeuk coordinate problem
            if strand == "-":
                if int(taken_high_exon_coord) + int(taken_nucl_len) - 1 != int(
                    taken_low_exon_coord
                ):
                    taken_low_exon_coord = (
                        int(taken_high_exon_coord) + int(taken_nucl_len) - 1
                    )

            all_taken_low_exon_coords.append(taken_low_exon_coord)
            all_taken_high_exon_coords.append(taken_high_exon_coord)
            all_taken_exon_nucl_len.append(taken_nucl_len)

        details = {
            "T_acc": T_acc,
            "C_acc": C_acc,
            "S": strand,
            "bitscore": bitscore,
            "e-value": evalue,
            "num_exons": num_exons,
            "low_coord": low_coord,
            "high_coord": high_coord,
            "all_low_exon_coords": all_low_exon_coords,
            "all_taken_low_exon_coords": all_taken_low_exon_coords,
            "all_high_exon_coords": all_high_exon_coords,
            "all_taken_high_exon_coords": all_taken_high_exon_coords,
            "all_exon_nucl_len": all_exon_nucl_len,
            "all_taken_exon_nucl_len": all_taken_exon_nucl_len,
            "recorded_exon_coords": recorded_exon_coords,
            "gene_id": gene_id,
        }
        return details

    def parse_output(self):
        record_dict = {}
        structured_arr_contents = []
        try:
            with open(self.pred_protein_seqs, "r") as aa_file:
                for record in SeqIO.parse(aa_file, "fasta"):
                    header_details = self.parse_header(record.id)
                    if header_details["gene_id"] in record_dict:
                        raise BuscoError(
                            "{} already exists in dictionary".format(
                                header_details["gene_id"]
                            )
                        )
                    record_dict[header_details["gene_id"]] = record
                    record.id = header_details[
                        "gene_id"
                    ]  # simplify the header to just the match ID without all the scores and exon info
            with open(self.codon_file, "r") as nt_file:
                for record in SeqIO.parse(nt_file, "fasta"):
                    header_details = self.parse_header(record.id)
                    record.id = header_details["gene_id"]
                    record.name = header_details["gene_id"]
                    record.description = header_details["gene_id"]
                    seq_nt = record.seq
                    seq_aa = record_dict[record.id].seq
                    self.record_gene_details(
                        {
                            "target_id": header_details["T_acc"],
                            "contig_id": header_details["C_acc"],
                            "contig_start": header_details["low_coord"],
                            "contig_end": header_details["high_coord"],
                            "strand": header_details["S"],
                            "score": header_details["bitscore"],
                            "exon_coords": header_details["recorded_exon_coords"],
                            "aa_seq": seq_aa,
                            "nt_seq": seq_nt,
                            "run_number": self.run_number,
                        }
                    )
                    record.id = header_details["gene_id"]
                    structured_arr_contents.append(
                        (
                            record.id,
                            self.gene_details[record.id]["target_id"],
                            self.gene_details[record.id]["contig_id"],
                            self.gene_details[record.id]["contig_start"],
                            self.gene_details[record.id]["contig_end"],
                            self.gene_details[record.id]["strand"],
                            self.gene_details[record.id]["score"],
                            self.gene_details[record.id]["run_number"],
                        )
                    )
        except FileNotFoundError:
            raise NoRerunFile

        return structured_arr_contents

    def record_gene_details(self, details={}):
        """
        Record all required match details from gene predictor output. The exact contents of the dictionary will vary
        depending on the gene predictor and pipeline.
        :param details:
        :return:
        """

        gene_id = "{}|{}:{}-{}".format(
            details["target_id"],
            details["contig_id"],
            details["contig_start"],
            details["contig_end"],
        )

        aa_seq = SeqRecord(Seq(details["aa_seq"]), id=gene_id, description="")
        nt_seq = SeqRecord(Seq(details["nt_seq"]), id=gene_id, description="")

        self.gene_details[gene_id] = {
            "target_id": details["target_id"],
            "contig_id": details["contig_id"],
            "contig_start": details["contig_start"],
            "contig_end": details["contig_end"],
            "strand": details["strand"],
            "score": details["score"],
            "exon_coords": details["exon_coords"],
            "nt_seq": nt_seq,
            "aa_seq": aa_seq,
            "run_number": details["run_number"],
        }

        self.sequences_aa[gene_id] = aa_seq
        self.sequences_nt[gene_id] = nt_seq

    def create_filtered_sequence_files(self, structured_arr_contents):

        if len(self.gene_details) > 0:
            structured_arr = np.array(
                structured_arr_contents,
                dtype=[
                    ("gene_id", "U100"),
                    ("target_id", "U30"),
                    ("contig_id", "U70"),
                    ("low_coord", "i4"),
                    ("high_coord", "i4"),
                    ("strand", "U1"),
                    ("score", "f8"),
                    ("run_number", "i1"),
                ],
            )
            overlaps = self.find_overlaps(structured_arr)
            entries_to_remove = []
            for overlap in overlaps:
                match1 = overlap[0]
                match2 = overlap[1]
                gene_id1 = match1[0]
                gene_id2 = match2[0]
                pattern = r"(\d+at\d+)"  # match the BUSCO ID
                busco_id1 = re.search(pattern, match1[0]).group(1)
                busco_id2 = re.search(pattern, match2[0]).group(1)
                score1 = float(match1[6])
                score2 = float(match2[6])
                if gene_id1 in entries_to_remove or gene_id2 in entries_to_remove:
                    continue
                elif busco_id1 == busco_id2:
                    if score1 > score2:  # score
                        entry_to_remove = gene_id2
                    else:
                        entry_to_remove = gene_id1
                    entries_to_remove.append(entry_to_remove)

            filtered_records_faa = []
            filtered_records_fna = []

            for gene_id, details in self.gene_details.items():
                if gene_id in entries_to_remove:
                    continue
                filtered_records_faa.append(details["aa_seq"])
                filtered_records_fna.append(details["nt_seq"])

            with open(self.pred_protein_seqs_modified, "w") as f_mod:
                SeqIO.write(filtered_records_faa, f_mod, "fasta")
            with open(self.codon_file_modified, "w") as g_mod:
                SeqIO.write(filtered_records_fna, g_mod, "fasta")

        return

    def write_gff_files(self, sc_folder, mc_folder, frag_folder):
        try:
            with open(self.gff_file, "r") as gf:
                lines = gf.readlines()
        except FileNotFoundError:
            logger.warning(
                "Metaeuk did not create a GFF file. Please use Metaeuk version 5-34c21f2 or higher for GFF results."
            )
            return

        id_pattern = re.compile(r".*Target_ID=(.*?)[_;]")
        current_busco = ""
        f = None

        sc_busco_ids = set([s.split(".f")[0] for s in os.listdir(sc_folder)])
        mc_busco_ids = set([s.split(".f")[0] for s in os.listdir(mc_folder)])
        frag_busco_ids = set([s.split(".f")[0] for s in os.listdir(frag_folder)])
        unused_busco_ids = set()
        used_busco_ids = (
            set()
        )  # Needed to track which BUSCOs are used on the rerun vs initial run

        for line in lines:
            m = id_pattern.match(line)
            if m:
                busco_id = m.group(1)
                if busco_id in unused_busco_ids:
                    continue
                elif busco_id in sc_busco_ids:
                    output_folder = sc_folder
                elif busco_id in mc_busco_ids:
                    output_folder = mc_folder
                elif busco_id in frag_busco_ids:
                    output_folder = frag_folder
                else:
                    unused_busco_ids.add(busco_id)
                    continue

                if busco_id != current_busco:
                    if f:
                        if not f.closed:
                            f.close()
                    gff_filename = "{}.gff".format(
                        os.path.join(output_folder, busco_id)
                    )

                    if busco_id not in used_busco_ids and os.path.exists(
                        gff_filename
                    ):  # don't overwrite the result from the rerun with the inital run GFF file
                        unused_busco_ids.add(busco_id)
                        current_busco = busco_id
                        continue

                    else:
                        f = open(gff_filename, "a")
                        # I avoided using a context manager for the file IO to avoid multiple file open/close
                        # operations and keep a single pass through the large GFF file. Efficiency!
                        current_busco = busco_id
                        used_busco_ids.add(busco_id)
                f.write(line)
        if f:
            if not f.closed:
                f.close()

    def reset(self):
        super().reset()
