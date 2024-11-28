# coding: utf-8
"""
hmmer.py

Module for running HMMER.

Author(s): Matthew Berkeley, Mathieu Seppey, Mose Manni, Felipe Simao, Rob Waterhouse

Copyright (c) 2015-2024, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

from busco.busco_tools.base import BaseRunner
import os
from collections import defaultdict
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.BuscoConfig import BaseConfig
from Bio import SeqIO
import csv
import subprocess
from busco.BuscoConfig import BuscoConfigAuto
from busco.Exceptions import BatchFatalError, BuscoError
import busco
import numpy as np

logger = BuscoLogger.get_logger(__name__)


class HMMERRunner(BaseRunner):

    name = "hmmsearch"
    cmd = "hmmsearch"

    def __init__(self):
        super().__init__()
        self._hmmer_output_folder = os.path.join(self.run_folder, "hmmer_output")
        self.datasets_version = self.config.get("busco_run", "datasets_version")
        self.dataset_creation_date = self.config.get("busco_run", "creation_date")
        self.dataset_nb_species = self.config.get("busco_run", "number_of_species")
        self.dataset_nb_buscos = self.config.get("busco_run", "number_of_BUSCOs")
        self.domain = self.config.get("busco_run", "domain")

        self.single_copy_sequences_folder = os.path.join(
            self.run_folder, "busco_sequences", "single_copy_busco_sequences"
        )
        self.multi_copy_sequences_folder = os.path.join(
            self.run_folder, "busco_sequences", "multi_copy_busco_sequences"
        )
        self.fragmented_sequences_folder = os.path.join(
            self.run_folder, "busco_sequences", "fragmented_busco_sequences"
        )

        self.cutoff_dict = {}
        self.single_copy_buscos = {}
        self.multi_copy_buscos = {}
        self.fragmented_buscos = {}
        self.extra_columns = False
        self.log_count = 0  # variable used to skip logging for intermediate eukaryote pipeline results.
        self.one_line_summary = None
        self.one_line_summary_raw = None

        # to be initialized before run time
        self.input_sequences = None
        self.busco_ids = None
        self.mode = self.config.get("busco_run", "mode")
        self.gene_details = {}
        self.results_dir = None

        self.matched_genes_complete = {}
        self.matched_genes_vlarge = {}
        self.matched_genes_fragment = {}
        self.is_complete = {}
        self.is_fragment = {}
        self.is_very_large = {}

        self.create_dirs(self._hmmer_output_folder)
        self.create_dirs(
            [
                self.single_copy_sequences_folder,
                self.multi_copy_sequences_folder,
                self.fragmented_sequences_folder,
            ],
            overwrite=True,  # avoid a clash with a previous run result in restart mode
        )
        if (
            self.domain == "eukaryota" and self.mode != "euk_genome_min"
        ):  # miniprot pipeline only runs once
            self.initial_results_dir = os.path.join(
                self._hmmer_output_folder, "initial_run_results"
            )
            self.rerun_results_dir = os.path.join(
                self._hmmer_output_folder, "rerun_results"
            )
            self.create_dirs([self.initial_results_dir, self.rerun_results_dir])

        else:
            self.initial_results_dir = self._hmmer_output_folder
            self.rerun_results_dir = None

        self.single_copy = 0
        self.multi_copy = 0
        self.only_fragments = 0
        self.total_buscos = 0

        # Get percentage of each kind of BUSCO match
        self.s_percent = 0
        self.d_percent = 0
        self.f_percent = 0
        self.e_percent = 0
        self.complete_percent = 0
        self.missing_percent = 0

        self.hmmer_results_lines = None
        self._already_used_genes = None
        self.missing_buscos = None

        self.miniprot_pipeline = False

    def configure_runner(self, input_sequences, busco_ids, mode, gene_details):
        super().configure_runner()
        self.run_number += 1
        self.input_sequences = input_sequences
        self.busco_ids = busco_ids
        self.mode = mode

        self.single_copy_buscos = {}
        self.multi_copy_buscos = {}
        self.fragmented_buscos = {}

        self._already_used_genes = set()
        self.hmmer_results_lines = []
        self.missing_buscos = []
        self.gene_details = gene_details
        if len(self.cutoff_dict) == 0:
            self.load_buscos()

        if (
            self.domain == "eukaryota" and self.mode != "euk_genome_min"
        ):  # miniprot pipeline only runs once
            if self.run_number == 1:
                self.results_dir = self.initial_results_dir
            elif self.run_number == 2:
                self.results_dir = self.rerun_results_dir
            else:
                raise ValueError(
                    "HMMER should not be run more than twice in the same Run instance."
                )
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
        hmmer_job.add_parameter(
            os.path.join(self.lineage_dataset, "hmms", "{}.hmm".format(busco_id))
        )
        hmmer_job.add_parameter(seq_filename)
        return hmmer_job

    def generate_job_args(self):
        for busco_id in self.busco_ids:
            if busco_id in self.cutoff_dict:
                if isinstance(self.input_sequences, str):
                    output_filename = "{}.out".format(busco_id)
                    yield busco_id, self.input_sequences, output_filename
                elif isinstance(self.input_sequences, list):
                    input_files = [
                        f
                        for f in self.input_sequences
                        if os.path.basename(f).startswith(busco_id)
                    ]
                    for seq_filename in input_files:
                        filename_parts = os.path.basename(seq_filename).rpartition(
                            ".faa"
                        )
                        output_filename = (
                            filename_parts[0] + ".out" + filename_parts[-1]
                        )
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
        if self.config.get("busco_run", "datasets_version") == "odb10":
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
                    input_files = [
                        f
                        for f in self.input_sequences
                        if os.path.basename(f).startswith(busco_id)
                    ]
                    n += len(input_files)
        return n

    def get_version(self):
        """
        check the Tool has the correct version
        """
        hmmer_version = subprocess.check_output(
            [self.cmd, "-h"], stderr=subprocess.STDOUT, shell=False
        )
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
        :raises BatchFatalError: if a Tool version is not supported
        """
        # check hmm version
        if not self.version >= BaseConfig.HMMER_VERSION:
            raise BatchFatalError(
                "HMMer version detected is not supported, please use HMMer v.{} +".format(
                    BaseConfig.HMMER_VERSION
                )
            )
        return

    def merge_dicts(self):
        merged_dict = defaultdict(lambda: defaultdict(list))
        for hmmer_dict in [self.is_complete, self.is_very_large, self.is_fragment]:
            for busco_id, busco_matches in hmmer_dict.items():
                for gene_id, matches in busco_matches.items():
                    merged_dict[busco_id][gene_id].extend(matches)
        return dict(merged_dict)

    def sort_hmmer_results(self, arr):
        hmm_match_lengths = {}
        match_scores = {}
        matches = np.unique(arr["target_name"])
        for match in matches:
            all_domain_hits = arr[arr["target_name"] == match]
            scores = np.unique(
                all_domain_hits["score"]
            )  # in transcriptome mode the same gene id can have multiple reading frame matches. Take the highest scoring match.
            top_scoring_domain_hits = all_domain_hits[
                all_domain_hits["score"] == max(scores)
            ]
            hmm_coords = []
            for hit in top_scoring_domain_hits:
                hmm_coords.append((hit["hmm_coord_from"], hit["hmm_coord_to"]))
            match_len = self.sum_hmm_len(hmm_coords)
            hmm_match_lengths[match] = match_len
            match_scores[match] = max(
                scores
            )  # use score as a tiebreaker in case of identical match lengths

        target_hits_sorted = sorted(
            matches, key=lambda x: (-hmm_match_lengths[x], -match_scores[x])
        )  # sort descending by match length, then by score
        return target_hits_sorted, hmm_match_lengths

    def parse_hmmer_output(self, filename, busco_query):

        processed_genes = []
        matched_genes = []
        records = defaultdict(list)

        data = []
        with open(filename, "r") as f:
            lines = f.readlines()

        for line in lines:
            if line.startswith("#"):
                continue
            else:
                values = line.split()
                description = "".join(values[22:])
                row_data = tuple(values[:22]) + (description,)
                data.append(row_data)

        arr = np.array(
            data,
            dtype=[
                ("target_name", "U500"),
                ("target_accession", "U500"),
                ("tlen", "i4"),
                ("query_name", "U500"),
                ("query_accession", "U500"),
                ("qlen", "i4"),
                ("eval", "f8"),
                ("score", "f8"),
                ("bias", "f8"),
                ("dom_num", "i4"),
                ("ndom", "i4"),
                ("c-eval", "f8"),
                ("i-eval", "f8"),
                ("dom_score", "f8"),
                ("dom_bias", "f8"),
                ("hmm_coord_from", "i4"),
                ("hmm_coord_to", "i4"),
                ("ali_coord_from", "i4"),
                ("ali_coord_to", "i4"),
                ("env_coord_from", "i4"),
                ("env_coord_to", "i4"),
                ("acc", "f8"),
                ("description", "U500"),
            ],
        )

        target_hits_sorted, hmm_match_lengths = self.sort_hmmer_results(arr)

        for target in target_hits_sorted:
            all_domain_hits = arr[arr["target_name"] == target]
            gene_id = self.get_gene_id(target)
            score = all_domain_hits[0]["score"]

            if not self.mode == "proteins" and not target in processed_genes:
                if self._check_overlap(matched_genes, gene_id):
                    continue

            # Extract frame information (present in transcriptome mode)
            frame = (
                all_domain_hits[0]["description"]
                if "frame" in all_domain_hits[0]["description"]
                else None
            )
            # Store bitscore matches for each gene match. If match below cutoff, discard.
            if score < float(self.cutoff_dict[busco_query]["score"]):
                continue

            hmm_coords = [
                (hit["hmm_coord_from"], hit["hmm_coord_to"]) for hit in all_domain_hits
            ]
            env_coords = [
                (hit["env_coord_from"], hit["env_coord_to"]) for hit in all_domain_hits
            ]

            records[gene_id].append(
                {
                    "hmm_matched_length": hmm_match_lengths[target],
                    "hmm_profile_length": all_domain_hits[0]["qlen"],
                    "hmm_coords": hmm_coords,
                    "env_coords": env_coords,
                    "score": score,
                    "frame": frame,
                    "orig gene ID": target,
                    "ref gene ID": target,
                }
            )

            if gene_id not in matched_genes:
                matched_genes.append(gene_id)
                processed_genes.append(target)
        return records

    def get_gene_id(self, target):
        if self.mode == "proteins":
            return target
        else:
            parts = target.split("|")
            if len(parts) > 2:
                return "|".join(parts[1:-1])  # allow for contig IDs to contain pipes
            else:
                return parts[-2]

    @staticmethod
    def _check_overlap(matched_genes, gene2):
        overlaps = []
        contig2, coords2 = gene2.rsplit(":", 1)
        for gene1 in matched_genes:
            contig1, coords1 = gene1.rsplit(":", 1)
            if contig1 != contig2:
                continue
            start1, end1 = coords1.split(
                "-"
            )  # start should always be a lower number than end at this point
            start2, end2 = coords2.split("-")
            if int(end2) - int(start2) > int(end1) - int(start1):
                start1, end1, start2, end2 = start2, end2, start1, end1
            if int(start1) <= int(start2) <= int(end1) or int(start1) <= int(
                end2
            ) <= int(end1):
                overlaps.append(True)
            else:
                overlaps.append(False)
        return any(overlaps)

    def sum_hmm_len(self, hmm_coords):
        hmm_regions_sorted = sorted(hmm_coords, key=lambda x: x[0])
        hmm_regions_used = []
        for region in hmm_regions_sorted:
            for used in hmm_regions_used:
                if (
                    used[0] <= region[0] <= used[1]
                ):  # if any new coord contained in previous region, expand if necessary
                    if region[1] > used[1]:
                        used[1] = region[1]
                    break
            else:
                hmm_regions_used.append(list(region))

        return sum([region[1] - region[0] + 1 for region in hmm_regions_used])

    def _sort_matches(self, matched_record, busco_query):
        """
        The HMMER gene matches are sorted into "complete", "v_large" and "fragmented" matches based on a comparison
        with the cutoff value specified in the dataset cutoff_scores file
        :param matched_record: dict of (gene_id, total_matched_length) pairs
        :param busco_query: BUSCO identifier
        :type matched_record: dict
        :type busco_query: str
        :return: busco_complete, busco_vlarge, busco_fragment - three dictionaries of the form
        {gene_id: [{"bitscore": float, "length": int}, {...}, ...], ...}
        :rtype: dict
        """
        busco_complete = defaultdict(list)
        busco_vlarge = defaultdict(list)
        busco_fragment = defaultdict(list)
        matched_genes_complete = defaultdict(list)
        matched_genes_vlarge = defaultdict(list)
        matched_genes_fragment = defaultdict(list)

        # Determine whether matched gene represents a complete, very_large or fragment of a BUSCO
        for gene_id, record in matched_record.items():
            orig_gene_ids = set([r["orig gene ID"] for r in record])
            for orig_gene_id in orig_gene_ids:
                subrecord = [
                    r for r in record if r["orig gene ID"] == orig_gene_id
                ]  # needed for miniprot pipeline,
                # where identical gene coords can match different reference sequences for the same BUSCO

                size = subrecord[0]["hmm_matched_length"]

                if self.config.get("busco_run", "datasets_version") == "odb10":

                    # Kind of like a z-score, but it is compared with a cutoff value, not a mean
                    zeta = (
                        self.cutoff_dict[busco_query]["length"] - size
                    ) / self.cutoff_dict[busco_query]["sigma"]

                    # gene match can only be either complete, v_large or fragment
                    if -2 <= zeta <= 2:
                        busco_type = busco_complete
                        match_type = matched_genes_complete
                    elif zeta < -2:
                        busco_type = busco_vlarge
                        match_type = matched_genes_vlarge
                    else:
                        busco_type = busco_fragment
                        match_type = matched_genes_fragment

                else:  # logic is made more concise for odb12 datasets and later

                    profile_length = subrecord[0]["hmm_profile_length"]

                    if size >= 0.8 * profile_length:
                        busco_type = busco_complete
                        match_type = matched_genes_complete

                    else:
                        busco_type = busco_fragment
                        match_type = matched_genes_fragment

                # Add information about match to dict
                busco_type[gene_id].append(
                    dict({"bitscore": subrecord[0]["score"], "length": size})
                )
                for key, val in subrecord[
                    0
                ].items():  # allows miniprot-specific keyvals to be added simply without
                    # having to reference other modules
                    if not key in busco_type[gene_id][-1]:
                        busco_type[gene_id][-1][key] = val

                # Reference which busco_queries are associated with each gene match
                match_type[gene_id].append(busco_query)

        return (
            busco_complete,
            busco_vlarge,
            busco_fragment,
            matched_genes_complete,
            matched_genes_vlarge,
            matched_genes_fragment,
        )

    def process_output(self, gene_id_lookup=None):
        """
        Load all gene matches from HMMER output and sort into dictionaries depending on match quality
        (complete, v_large, fragment).
        :return:
        """
        if self.run_number == 1:
            hmmer_results_files = sorted(
                [
                    os.path.join(self.results_dir, f)
                    for f in os.listdir(self.results_dir)
                    if not f.startswith(".")
                ]
            )
        elif self.run_number == 2:
            hmmer_rerun_files = [
                os.path.join(self.rerun_results_dir, f)
                for f in os.listdir(self.rerun_results_dir)
                if not f.startswith(".")
            ]
            hmmer_results_files = sorted(hmmer_rerun_files)
        else:
            raise ValueError(
                "HMMER should not be run more than twice in the same Run instance."
            )
        self.gene_id_lookup = gene_id_lookup

        hmmer_records = []
        for (
            filename
        ) in hmmer_results_files:  # quicker to run in serial than in parallel
            results = self.load_results_from_file(filename)
            hmmer_records.append(results)

        self.unpack_hmmer_records(hmmer_records)

    def unpack_hmmer_records(self, hmmer_records):

        self.is_complete = defaultdict(
            lambda: defaultdict(list), self.is_complete
        )  # dict of a dict of lists of dicts
        self.is_fragment = defaultdict(lambda: defaultdict(list), self.is_fragment)
        self.is_very_large = defaultdict(lambda: defaultdict(list), self.is_very_large)
        self.matched_genes_complete = defaultdict(list, self.matched_genes_complete)
        self.matched_genes_vlarge = defaultdict(list, self.matched_genes_vlarge)
        self.matched_genes_fragment = defaultdict(list, self.matched_genes_fragment)

        for records in hmmer_records:
            (
                busco_query,
                busco_complete,
                busco_vlarge,
                busco_fragment,
                matched_genes_complete,
                matched_genes_vlarge,
                matched_genes_fragment,
            ) = records

            # Add all information for this busco_id to the full dictionary
            if len(busco_complete) > 0:
                self.is_complete[busco_query].update(busco_complete)
            if len(busco_vlarge) > 0:
                self.is_very_large[busco_query].update(busco_vlarge)
            if len(busco_fragment) > 0:
                self.is_fragment[busco_query].update(busco_fragment)

            for i in range(3):
                matched_genes_dict_small = [
                    matched_genes_complete,
                    matched_genes_vlarge,
                    matched_genes_fragment,
                ][i]
                matched_genes_dict_large = [
                    self.matched_genes_complete,
                    self.matched_genes_vlarge,
                    self.matched_genes_fragment,
                ][i]
                for gene_id in matched_genes_dict_small:
                    if gene_id in matched_genes_dict_large:
                        matched_genes_dict_large[gene_id].extend(
                            matched_genes_dict_small[gene_id]
                        )
                    else:
                        matched_genes_dict_large[gene_id] = matched_genes_dict_small[
                            gene_id
                        ]

        self.is_complete = dict(self.is_complete)
        self.is_fragment = dict(self.is_fragment)
        self.is_very_large = dict(self.is_very_large)
        self.matched_genes_complete = dict(self.matched_genes_complete)
        self.matched_genes_vlarge = dict(self.matched_genes_vlarge)
        self.matched_genes_fragment = dict(self.matched_genes_fragment)

        return

    def load_results_from_file(self, filename):
        busco_query = str(os.path.basename(filename).split(".")[0])
        matched_record = self.parse_hmmer_output(filename, busco_query)

        (
            busco_complete,
            busco_vlarge,
            busco_fragment,
            matched_genes_complete,
            matched_genes_vlarge,
            matched_genes_fragment,
        ) = self._sort_matches(matched_record, busco_query)
        return (
            busco_query,
            busco_complete,
            busco_vlarge,
            busco_fragment,
            matched_genes_complete,
            matched_genes_vlarge,
            matched_genes_fragment,
        )

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
            higher_rank_buscos = list(self.is_complete.keys()) + list(
                self.is_very_large.keys()
            )
            matched_genes = self.matched_genes_fragment
        else:
            raise BuscoError("Unrecognized dictionary of BUSCOs.")

        for busco_id in list(busco_dict.keys()):
            matches = busco_dict[busco_id]
            # Remove any buscos that appear in higher ranking dictionaries
            if busco_id in higher_rank_buscos:
                busco_dict.pop(busco_id)
                for gene_id in matches:
                    matched_genes[gene_id] = [
                        x for x in matched_genes[gene_id] if x != busco_id
                    ]  # Remove all occurences of busco_id
                    if len(matched_genes[gene_id]) == 0:
                        matched_genes.pop(gene_id)
                continue

            # Remove any genes that have previously been processed under a different and higher ranking busco match
            for gene_id in list(matches.keys()):
                if gene_id in self._already_used_genes:
                    busco_dict[busco_id].pop(gene_id)
                    matched_genes[gene_id] = [
                        x for x in matched_genes[gene_id] if x != busco_id
                    ]  # Remove all occurences of busco_id
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
        matched_genes = self.get_matched_genes_dict(busco_dict)

        busco_matches_to_remove = []
        # Keep the best scoring gene if gene is matched by more than one busco with the same match rank
        for gene_id, buscos in matched_genes.items():
            if len(buscos) > 1:
                busco_bitscores = []
                busco_matches = []
                for busco in buscos:
                    if not busco in busco_matches:
                        matches = busco_dict[busco][gene_id]
                        for match in matches:
                            bitscore = match["bitscore"]
                            busco_bitscores.append(bitscore)
                        busco_matches.append(busco)

                if (
                    len(set(buscos)) == 1
                ):  # If only one busco is matched twice (initial run and rerun), don't remove it
                    continue
                best_match_ind = max(
                    range(len(busco_bitscores)), key=busco_bitscores.__getitem__
                )
                buscos = [x for x in buscos if x != busco_matches[best_match_ind]]
                # Remove lower scoring duplicates from dictionary.

                for duplicate in list(set(buscos)):
                    # Use set to account for any duplicate entries (matched in both initial run and rerun)
                    busco_dict[duplicate].pop(gene_id)
                    if len(busco_dict[duplicate]) == 0:
                        busco_dict.pop(duplicate)
                    busco_matches_to_remove.append((gene_id, duplicate))

        for gene_busco_pair in busco_matches_to_remove:
            gene_id, busco_id = gene_busco_pair
            matched_genes[gene_id].remove(busco_id)
            if len(matched_genes[gene_id]) == 0:
                matched_genes.pop(gene_id)

        return

    def get_matched_genes_dict(self, busco_dict):
        if busco_dict == self.is_complete:
            matched_genes = self.matched_genes_complete
        elif busco_dict == self.is_very_large:
            matched_genes = self.matched_genes_vlarge
        elif busco_dict == self.is_fragment:
            matched_genes = self.matched_genes_fragment
        else:
            raise BuscoError("Unrecognized dictionary of BUSCOs.")
        return matched_genes

    def _remove_low_scoring_matches(self, busco_dict):
        """
        Go through input dictionary and remove any gene matches that score less than 85% of the top gene match score
        for each BUSCO.
        :param busco_dict: one of [self.is_complete, self.is_very_large, self.is_fragment]
        :type busco_dict: dict
        :return:
        """
        empty_buscos = []

        matched_genes = self.get_matched_genes_dict(busco_dict)

        # For each busco, keep only matches within 85% of top bitscore match for that busco
        for busco_id, matches in busco_dict.items():
            if len(matches) > 1:
                _, max_bitscore = self._get_best_scoring_match(matches)
                # Go through all matches again, removing any below the threshold
                for gene_id in list(matches.keys()):
                    match_info = matches[gene_id]
                    matches_to_remove = []
                    for m, match in enumerate(match_info):
                        if match["bitscore"] < 0.85 * max_bitscore:
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
            matched_genes[gene_id].remove(busco_id)
            if len(matched_genes[gene_id]) == 0:
                matched_genes.pop(gene_id)
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

    def filter(self):
        """
        Remove all duplicate matches and any matches below 85% of the top match for each BUSCO.
        :return:
        """
        self._remove_duplicates()
        self._remove_low_scoring_matches(self.is_complete)
        self._remove_low_scoring_matches(self.is_very_large)
        self._remove_low_scoring_matches(self.is_fragment)
        return

    def consolidate_busco_lists(self):
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
                self.fragmented_buscos[busco_id] = {
                    best_fragment: self.is_fragment[busco_id][best_fragment]
                }
            else:
                self.fragmented_buscos[busco_id] = gene_matches
        return

    def load_links_info(self):
        links_info = defaultdict(dict)
        links_file = os.path.join(
            self.lineage_dataset,
            "links_to_{}.txt".format(self.datasets_version.upper()),
        )
        if os.path.exists(links_file):
            with open(links_file, newline="") as f:
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
                            output_lines.append(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    busco,
                                    label,
                                    gene_id,
                                    bit_score,
                                    match_length,
                                    link,
                                    desc,
                                )
                            )
                        except KeyError:
                            output_lines.append(
                                "{}\t{}\t{}\t{}\t{}\n".format(
                                    busco, label, gene_id, bit_score, match_length
                                )
                            )
                    elif self.mode == "genome":
                        scaffold = self.gene_details[match["ref gene ID"]]
                        if "strand" in scaffold:
                            location_pattern = ":{}-{}".format(
                                scaffold["contig_start"], scaffold["contig_end"]
                            )
                            if gene_id.endswith(location_pattern):
                                gene_id = gene_id.replace(location_pattern, "")

                        start, end = sorted(
                            [scaffold["contig_start"], scaffold["contig_end"]],
                            reverse=bool(scaffold["strand"] == "-"),
                        )
                        try:
                            desc = links_info[busco]["description"]
                            link = links_info[busco]["link"]
                            self.extra_columns = True
                            output_lines.append(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    busco,
                                    label,
                                    gene_id,
                                    start,
                                    end,
                                    scaffold["strand"],
                                    bit_score,
                                    match_length,
                                    link,
                                    desc,
                                )
                            )
                        except KeyError:
                            output_lines.append(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    busco,
                                    label,
                                    gene_id,
                                    start,
                                    end,
                                    scaffold["strand"],
                                    bit_score,
                                    match_length,
                                )
                            )
        return output_lines

    def create_output_content(self):
        """
        Format output for all BUSCO matches.
        :return: output_lines
        :rtype: list
        """
        output_lines = []
        dict_labels = {
            "Complete": self.single_copy_buscos,
            "Duplicated": self.multi_copy_buscos,
            "Fragmented": self.fragmented_buscos,
        }
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
            if not any(
                busco_group in d
                for d in [
                    self.single_copy_buscos,
                    self.multi_copy_buscos,
                    self.fragmented_buscos,
                ]
            ):
                output_lines.append("{}\tMissing\n".format(busco_group))
                self.missing_buscos.append(busco_group)

        if len(self.missing_buscos) == len(self.cutoff_dict):
            logger.warning(
                "BUSCO did not find any match. Make sure to check the log files if this is unexpected."
            )

        return output_lines, self.missing_buscos

    def _load_length(self):
        """
        This function loads the length cutoffs file
        :raises BuscoError: if the lengths_cutoff file cannot be read
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
                        logger.error("Error parsing the lengths_cutoff file.")
                        raise BuscoError(e)
        except IOError:
            raise BuscoError(
                "Impossible to read the lengths in {}".format(
                    os.path.join(lengths_cutoff_file)
                )
            )
        return

    def _load_score(self):
        """
        This function loads the score cutoffs file
        :raises BuscoError: if the scores_cutoff file cannot be read
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
                        logger.error("Error parsing the scores_cutoff file.")
                        raise BuscoError(e)
        except IOError:
            raise BuscoError(
                "Impossible to read the scores in {}".format(scores_cutoff_file)
            )
        return

    def write_buscos_to_file(self):
        """
        Write BUSCO matching sequences to output fasta files. Each sequence is printed in a separate file and both
        nucleotide and amino acid versions are created.
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
                aa_seqs = []
                nt_seqs = []
                try:
                    for gene_id, matches in gene_matches.items():
                        for match in matches:
                            try:
                                nt_seqs.append(
                                    self.gene_details[match["ref gene ID"]]["nt_seq"]
                                )
                            except KeyError:
                                pass
                            aa_seqs.append(
                                self.gene_details[match["ref gene ID"]]["aa_seq"]
                            )
                    if len(nt_seqs) > 0:
                        with open(
                            os.path.join(output_dir, "{}.fna".format(busco)), "w"
                        ) as f2:
                            SeqIO.write(nt_seqs, f2, "fasta")
                except TypeError:
                    for gene_id, matches in gene_matches.items():
                        for match in matches:
                            aa_seqs.append(
                                self.gene_details[match["ref gene ID"]]["aa_seq"]
                            )
                with open(os.path.join(output_dir, "{}.faa".format(busco)), "w") as f1:
                    SeqIO.write(aa_seqs, f1, "fasta")
        return

    def write_hmmer_results(self, output_lines):
        """
        Create two output files: one with information on all BUSCOs for the given dataset and the other with a list of
        all BUSCOs that were not found.
        :return:
        """

        with open(os.path.join(self.run_folder, "full_table.tsv"), "w") as f_out:

            self.write_output_header(f_out)

            with open(
                os.path.join(self.run_folder, "missing_busco_list.tsv"), "w"
            ) as miss_out:

                self.write_output_header(miss_out, missing_list=True)

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
        internal_stop_codon_pattern = "(of which {} contain internal stop codons)"

        self.hmmer_results_lines.append("***** Results: *****\n\n")
        self.hmmer_results_lines.append(self.one_line_summary_raw)
        self.hmmer_results_lines.append(
            "{}\tComplete BUSCOs (C)\t{}\t\t{}\n".format(
                self.single_copy + self.multi_copy,
                internal_stop_codon_pattern.format(self.complete_stop_codon_count)
                if self.complete_stop_codon_count > 0
                else "",
                "   ",
            )
        )
        self.hmmer_results_lines.append(
            "{}\tComplete and single-copy BUSCOs (S)\t{}\n".format(
                self.single_copy, "   "
            )
        )
        self.hmmer_results_lines.append(
            "{}\tComplete and duplicated BUSCOs (D)\t{}\n".format(
                self.multi_copy, "   "
            )
        )
        self.hmmer_results_lines.append(
            "{}\tFragmented BUSCOs (F)\t\t\t{}\n".format(self.only_fragments, "   ")
        )
        self.hmmer_results_lines.append(
            "{}\tMissing BUSCOs (M)\t\t\t{}\n".format(
                self.total_buscos
                - self.single_copy
                - self.multi_copy
                - self.only_fragments,
                "   ",
            )
        )
        self.hmmer_results_lines.append(
            "{}\tTotal BUSCO groups searched\t\t{}\n".format(
                self.dataset_nb_buscos, "   "
            )
        )

        if isinstance(self.config, BuscoConfigAuto):
            self._log_one_line_hmmer_summary()
        elif self.domain == "eukaryota" and self.log_count == 0:
            self.log_count += 1
            self._produce_full_hmmer_summary_debug()
        else:
            self._log_one_line_hmmer_summary()

        return

    def record_results(self):
        self._get_busco_percentages()
        extra = ",E:{}%".format(self.e_percent) if self.e_percent else ""
        self.one_line_summary_raw = (
            "C:{}%[S:{}%,D:{}%],F:{}%,M:{}%,n:{}{}\t{}\n".format(
                self.complete_percent,
                self.s_percent,
                self.d_percent,
                self.f_percent,
                self.missing_percent,
                self.total_buscos,
                extra,
                "   ",
            )
        )
        self.one_line_summary = "Results:\t{}".format(self.one_line_summary_raw)

    @log("{}", logger, attr_name="hmmer_results_lines", apply="join", on_func_exit=True)
    def _produce_full_hmmer_summary(self):
        return

    @log(
        "{}",
        logger,
        attr_name="hmmer_results_lines",
        apply="join",
        on_func_exit=True,
        debug=True,
    )
    def _produce_full_hmmer_summary_debug(self):
        return

    @log("{}", logger, attr_name="one_line_summary")
    def _log_one_line_hmmer_summary(self):
        return

    def write_output_header(
        self, file_object, missing_list=False, no_table_header=False
    ):
        """
        Write a standardized file header containing information on the BUSCO run.
        :param file_object: Opened file object ready for writing
        :type file_object: file
        :param missing_list: Add list of missing BUSCOs
        :type missing_list: bool
        :param no_table_header: Include table header
        :type no_table_header: bool
        :return:
        """
        file_object.write(
            "# BUSCO version is: {} \n"
            "# The lineage dataset is: {} (Creation date: {}, number of genomes: {}, number of BUSCOs: {}"
            ")\n".format(
                busco.__version__,
                os.path.basename(self.lineage_dataset),
                self.dataset_creation_date,
                self.dataset_nb_species,
                self.dataset_nb_buscos,
            )
        )

        if no_table_header:
            pass
        elif missing_list:
            file_object.write("# Busco id\n")
        elif self.mode == "proteins" or self.mode == "transcriptome":
            if self.extra_columns:
                file_object.write(
                    "# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n"
                )
            else:
                file_object.write("# Busco id\tStatus\tSequence\tScore\tLength\n")
        elif self.mode == "genome":
            if self.extra_columns:
                file_object.write(
                    "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tOrthoDB url"
                    "\tDescription\n"
                )
            else:
                file_object.write(
                    "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"
                )

        return

    def _get_busco_percentages(self):
        self.single_copy = len(self.single_copy_buscos)  # int
        self.multi_copy = len(self.multi_copy_buscos)  # int
        self.only_fragments = len(self.fragmented_buscos)  # int
        self.total_buscos = len(self.cutoff_dict)
        self.num_missing = (
            self.total_buscos - self.single_copy - self.multi_copy - self.only_fragments
        )

        # Get percentage of each kind of BUSCO match
        self.s_percent = abs(round((self.single_copy / self.total_buscos) * 100, 1))
        self.d_percent = abs(round((self.multi_copy / self.total_buscos) * 100, 1))
        self.f_percent = abs(round((self.only_fragments / self.total_buscos) * 100, 1))
        self.complete_percent = abs(
            round(((self.single_copy + self.multi_copy) / self.total_buscos) * 100, 1)
        )
        self.missing_percent = abs(
            round((self.num_missing / self.total_buscos) * 100, 1)
        )

        if self.miniprot_pipeline:
            self._get_stop_codon_percent_and_avg_identity()
        else:
            self.e_percent = 0
            self.avg_identity = None
            self.complete_stop_codon_count = 0

        return

    def _get_stop_codon_percent_and_avg_identity(self):
        single_copy_stop_codon_count = 0
        multi_copy_stop_codon_count = 0
        identities = []
        for busco_id, matches in self.single_copy_buscos.items():
            for gene_id, match in matches.items():
                ref_gene_id = match[0]["ref gene ID"]
                identities.append(self.gene_details[ref_gene_id]["identity"])
                if self.gene_details[ref_gene_id]["stop_codon_count"] > 0:
                    single_copy_stop_codon_count += 1

        for busco_id, matches in self.multi_copy_buscos.items():
            for gene_id, match in matches.items():
                ref_gene_id = match[0]["ref gene ID"]
                identities.append(self.gene_details[ref_gene_id]["identity"])
                if self.gene_details[ref_gene_id]["stop_codon_count"] > 0:
                    pass
                else:
                    break
            else:
                multi_copy_stop_codon_count += 1

        for busco_id, matches in self.fragmented_buscos.items():
            for gene_id, match in matches.items():
                ref_gene_id = match[0]["ref gene ID"]
                identities.append(self.gene_details[ref_gene_id]["identity"])

        self.complete_stop_codon_count = (
            single_copy_stop_codon_count + multi_copy_stop_codon_count
        )

        if self.complete_stop_codon_count > 0:
            complete_count = len(self.single_copy_buscos) + len(self.multi_copy_buscos)
            self.e_percent = round(
                100 * self.complete_stop_codon_count / complete_count, 1
            )
            logger.warning(
                "{} of {} Complete matches ({}%) contain internal stop codons in Miniprot gene predictions".format(
                    self.complete_stop_codon_count, complete_count, self.e_percent
                )
            )
        if len(identities) > 0:
            self.avg_identity = round(sum(identities) / len(identities), 2)
            if self.avg_identity < 0.5:
                logger.warning(
                    "BUSCO gene predictions from Miniprot have low average identity ({}). You may want to repeat the analysis using the Metaeuk pipeline.".format(
                        self.avg_identity
                    )
                )
        else:
            self.avg_identity = None

        return

    def get_error_count(self):
        error_count = 0
        for busco_id, busco_matches in self.single_copy_buscos.items():
            for gene_id, matches in busco_matches.items():
                for match in matches:
                    stop_codon_count = self.gene_details[match["ref gene ID"]][
                        "stop_codon_count"
                    ]
                    error_count += bool(stop_codon_count)

        for busco_id, busco_matches in self.multi_copy_buscos.items():
            duplicate_contains_error = False
            for gene_id, matches in busco_matches.items():
                for match in matches:
                    stop_codon_count = self.gene_details[match["ref gene ID"]][
                        "stop_codon_count"
                    ]
                    if bool(stop_codon_count):
                        error_count += bool(stop_codon_count)
                        duplicate_contains_error = True
                if duplicate_contains_error:
                    break
