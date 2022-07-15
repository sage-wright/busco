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
import pandas as pd
import busco

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
        self.log_count = 0  # Dummy variable used to skip logging for intermediate eukaryote pipeline results.
        self.one_line_summary = None
        self.one_line_summary_raw = None

        # to be initialized before run time
        self.input_sequences = None
        self.busco_ids = None
        self.mode = None
        self.gene_details = None
        self.results_dir = None

        self.matched_genes_complete = {}
        self.matched_genes_vlarge = {}
        self.matched_genes_fragment = {}
        self.is_complete = {}
        self.is_fragment = {}
        self.is_very_large = {}

        self.create_dirs(
            [
                self._hmmer_output_folder,
                self.single_copy_sequences_folder,
                self.multi_copy_sequences_folder,
                self.fragmented_sequences_folder,
            ]
        )
        if self.domain == "eukaryota":
            self.initial_results_dir = os.path.join(
                self._hmmer_output_folder, "initial_run_results"
            )
            self.rerun_results_dir = os.path.join(
                self._hmmer_output_folder, "rerun_results"
            )
            self.create_dirs([self.initial_results_dir, self.rerun_results_dir])

        self.single_copy = 0
        self.multi_copy = 0
        self.only_fragments = 0
        self.total_buscos = 0

        # Get percentage of each kind of BUSCO match
        self.s_percent = 0
        self.d_percent = 0
        self.f_percent = 0
        self.complete_percent = 0
        self.missing_percent = 0

        self.hmmer_results_lines = None
        self._already_used_genes = None
        self.missing_buscos = None

    def configure_runner(self, input_sequences, busco_ids, mode, gene_details):
        super().configure_runner()
        self.run_number += 1
        self.input_sequences = input_sequences
        self.busco_ids = busco_ids
        self.mode = mode

        self.is_fragment = {}
        self.matched_genes_fragment = {}

        self.single_copy_buscos = {}
        self.multi_copy_buscos = {}
        self.fragmented_buscos = {}

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

    def process_output(self):
        self.is_complete = defaultdict(
            lambda: defaultdict(list), self.is_complete
        )  # dict of a dict of lists of dicts
        self.is_fragment = defaultdict(lambda: defaultdict(list), self.is_fragment)
        self.is_very_large = defaultdict(lambda: defaultdict(list), self.is_very_large)
        self.matched_genes_complete = defaultdict(list, self.matched_genes_complete)
        self.matched_genes_vlarge = defaultdict(list, self.matched_genes_vlarge)
        self.matched_genes_fragment = defaultdict(list, self.matched_genes_fragment)

        self._load_matched_genes()

        self.is_complete = dict(self.is_complete)
        self.is_fragment = dict(self.is_fragment)
        self.is_very_large = dict(self.is_very_large)
        self.matched_genes_complete = dict(self.matched_genes_complete)
        self.matched_genes_vlarge = dict(self.matched_genes_vlarge)
        self.matched_genes_fragment = dict(self.matched_genes_fragment)

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

    def merge_dicts(self):
        merged_dict = defaultdict(lambda: defaultdict(list))
        for hmmer_dict in [self.is_complete, self.is_very_large, self.is_fragment]:
            for busco_id, busco_matches in hmmer_dict.items():
                merged_dict[busco_id].update(busco_matches)
        return merged_dict

    def parse_hmmer_output(self, filename, busco_query):
        """
        Read and parse HMMER output file.
        :param filename: Name of HMMER output file
        :param busco_query: Basename of file, used to identify BUSCO
        :type filename: str
        :type busco_query: str
        :return: Dictionary of (gene_id, total_matched_length) pairs
        :rtype: dict
        """
        records = defaultdict(dict)

        with open(filename, "r") as f:

            # Read HMMER output file
            for line in f:
                if line.startswith("#"):
                    continue
                else:
                    try:
                        line = line.strip().split()
                        gene_id = line[0]
                        tlen = int(line[2])
                        bit_score = float(line[7])

                        # Extract frame information (present in transcriptome mode)
                        frame = str(line[-1]) if "frame" in str(line[-1]) else None

                        # Store bitscore matches for each gene match. If match below cutoff, discard.
                        if bit_score < float(self.cutoff_dict[busco_query]["score"]):
                            continue
                        if gene_id not in records:
                            records[gene_id] = {
                                "tlen": tlen,
                                "hmm_len": 0,
                                "env_coords": [],
                                "score": bit_score,
                                "frame": frame,
                            }
                        hmm_start = int(line[15])
                        hmm_end = int(line[16])
                        env_start = int(line[19])
                        env_end = int(line[20])
                        records[gene_id]["hmm_len"] += hmm_end - hmm_start
                        records[gene_id]["env_coords"].append((env_start, env_end))

                    except IndexError as e:
                        logger.error(
                            "Cannot parse HMMER output file {}".format(filename)
                        )
                        raise BuscoError(e)
        return records

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

        # Determine whether matched gene represents a complete, very_large or fragment of a BUSCO
        for gene_id, record in matched_record.items():
            size = record["hmm_len"]
            frame = record["frame"]

            # Kind of like a z-score, but it is compared with a cutoff value, not a mean
            zeta = (self.cutoff_dict[busco_query]["length"] - size) / self.cutoff_dict[
                busco_query
            ]["sigma"]

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
            busco_type[gene_id].append(
                dict({"bitscore": record["score"], "length": size, "frame": frame})
            )
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

        for filename in hmmer_results_files:
            busco_query = str(os.path.basename(filename).split(".")[0])
            matched_record = self.parse_hmmer_output(filename, busco_query)
            filtered_records = self.remove_overlaps(matched_record)
            busco_complete, busco_vlarge, busco_fragment = self._sort_matches(
                filtered_records, busco_query
            )

            # Add all information for this busco_id to the full dictionary
            if len(busco_complete) > 0:
                self.is_complete[busco_query].update(busco_complete)
            if len(busco_vlarge) > 0:
                self.is_very_large[busco_query].update(busco_vlarge)
            if len(busco_fragment) > 0:
                self.is_fragment[busco_query].update(busco_fragment)

        return

    @staticmethod
    def remove_overlaps(matched_records):
        seq_ids = []
        start_coords = []
        end_coords = []
        scores = []
        try:
            for record in matched_records:
                seq_id, coords = record.split(":")
                start_coord, end_coord = coords.split("-")
                seq_ids.append(seq_id)
                start_coords.append(start_coord)
                end_coords.append(end_coord)
                scores.append(matched_records[record]["score"])
        except ValueError:  # for protein sequences there is no ":<coords>" suffix, so skip the overlap filtering
            return matched_records

        records_df = pd.DataFrame(
            {
                "Sequence": seq_ids,
                "Start": start_coords,
                "Stop": end_coords,
                "Score": scores,
            }
        )
        results_grouped = records_df.groupby("Sequence")
        entries_to_remove = []
        seq_ids = results_grouped.groups.keys()
        for seq in seq_ids:
            g1 = results_grouped.get_group(seq)
            g1_sorted = g1.sort_values(
                "Start"
            )  # sort to facilitate a single-pass coordinate check
            for idx1, row1 in g1_sorted.iterrows():
                current_record = g1_sorted.loc[idx1]
                start_val = current_record["Start"]
                stop_val = current_record["Stop"]
                current_seqid = "{}:{}-{}".format(
                    current_record["Sequence"], start_val, stop_val
                )
                if (
                    current_seqid in entries_to_remove
                ):  # overlaps with removed entries don't count
                    continue

                matches = g1_sorted[g1_sorted["Start"] > start_val].loc[
                    g1_sorted["Start"] < stop_val
                ]  # find entries with a start coordinate between the current exon start and end
                for idx2 in matches.index.values:
                    if g1_sorted.loc[idx1]["Score"] >= g1_sorted.loc[idx2]["Score"]:
                        ind_to_remove = idx2
                    else:
                        ind_to_remove = idx1
                    record_to_remove = g1_sorted.loc[ind_to_remove]
                    entries_to_remove.append(
                        "{}:{}-{}".format(
                            record_to_remove["Sequence"],
                            record_to_remove["Start"],
                            record_to_remove["Stop"],
                        )
                    )

        filtered_records = {
            i: matched_records[i] for i in matched_records if i not in entries_to_remove
        }

        return filtered_records

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
        if busco_dict == self.is_complete:
            matched_genes = self.matched_genes_complete
        elif busco_dict == self.is_very_large:
            matched_genes = self.matched_genes_vlarge
        elif busco_dict == self.is_fragment:
            matched_genes = self.matched_genes_fragment
        else:
            raise BuscoError("Unrecognized dictionary of BUSCOs.")

        busco_matches_to_remove = []
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
                        scaffold = self.gene_details[gene_id][m]
                        if self.domain == "eukaryota":
                            location_pattern = ":{}-{}".format(
                                scaffold["gene_start"], scaffold["gene_end"]
                            )
                            if gene_id.endswith(location_pattern):
                                gene_id = gene_id.replace(location_pattern, "")
                        else:  # Remove suffix assigned by Prodigal
                            gene_id = gene_id.rsplit("_", 1)[0]
                        try:
                            desc = links_info[busco]["description"]
                            link = links_info[busco]["link"]
                            self.extra_columns = True
                            output_lines.append(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    busco,
                                    label,
                                    gene_id,
                                    scaffold["gene_start"],
                                    scaffold["gene_end"],
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
                                    scaffold["gene_start"],
                                    scaffold["gene_end"],
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
                for d in [self.is_complete, self.is_very_large, self.is_fragment]
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
                    aa_seqs, nt_seqs = zip(
                        *[
                            (sequences_aa[gene_id], sequences_nt[gene_id])
                            for gene_id in gene_matches
                        ]
                    )
                    with open(
                        os.path.join(output_dir, "{}.fna".format(busco)), "w"
                    ) as f2:
                        SeqIO.write(nt_seqs, f2, "fasta")
                except TypeError:
                    aa_seqs = [sequences_aa[gene_id] for gene_id in gene_matches]
                with open(os.path.join(output_dir, "{}.faa".format(busco)), "w") as f1:
                    SeqIO.write(aa_seqs, f1, "fasta")

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

        self.hmmer_results_lines.append("***** Results: *****\n\n")
        self.hmmer_results_lines.append(self.one_line_summary_raw)
        self.hmmer_results_lines.append(
            "{}\tComplete BUSCOs (C)\t\t\t{}\n".format(
                self.single_copy + self.multi_copy, "   "
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
        self.one_line_summary_raw = "C:{}%[S:{}%,D:{}%],F:{}%,M:{}%,n:{}\t{}\n".format(
            self.complete_percent,
            self.s_percent,
            self.d_percent,
            self.f_percent,
            self.missing_percent,
            self.total_buscos,
            "   ",
        )
        self.one_line_summary = "Results:\t{}".format(self.one_line_summary_raw)
        # type(self).summary["results"][
        #     "one_line_summary"
        # ] = self.one_line_summary_raw.strip()
        # type(self).summary["results"]["C"] = self.single_copy + self.multi_copy
        # type(self).summary["results"]["S"] = self.single_copy
        # type(self).summary["results"]["D"] = self.multi_copy
        # type(self).summary["results"]["F"] = self.only_fragments
        # type(self).summary["results"]["M"] = (
        #     self.total_buscos - self.single_copy - self.multi_copy - self.only_fragments
        # )
        # type(self).summary["parameters"]["mode"] = self.mode
        #
        # if self.mode == "genome":
        #     if self.config.get("busco_run", "domain") in ["prokaryota", "viruses"]:
        #         gene_predictor = "prodigal"
        #     elif self.config.getboolean("busco_run", "use_augustus"):
        #         gene_predictor = "augustus"
        #     else:
        #         gene_predictor = "metaeuk"
        #     type(self).summary["parameters"]["gene_predictor"] = gene_predictor
        # type(self).summary["lineage_dataset"]["name"] = self.lineage_dataset
        # type(self).summary["lineage_dataset"][
        #     "creation_date"
        # ] = self.dataset_creation_date
        # type(self).summary["lineage_dataset"][
        #     "number_species"
        # ] = self.dataset_nb_species
        # type(self).summary["lineage_dataset"]["total_buscos"] = self.dataset_nb_buscos

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

        # Get percentage of each kind of BUSCO match
        self.s_percent = abs(round((self.single_copy / self.total_buscos) * 100, 1))
        self.d_percent = abs(round((self.multi_copy / self.total_buscos) * 100, 1))
        self.f_percent = abs(round((self.only_fragments / self.total_buscos) * 100, 1))
        self.complete_percent = round(self.s_percent + self.d_percent, 1)
        self.missing_percent = abs(
            round(100 - self.s_percent - self.d_percent - self.f_percent, 1)
        )

        return
