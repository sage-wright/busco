import os
import re
from collections import defaultdict
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import shutil
import csv

# todo: docstrings
logger = BuscoLogger.get_logger(__name__)


class ProdigalRunner:

    def __init__(self, prodigal_tool, input_file, output_dir, genetic_code, ambiguous_coding_density_range, log_path):
        self.prodigal_tool = prodigal_tool
        self.input_file = input_file
        self.output_dir = os.path.join(output_dir, "predicted_genes")
        self.tmp_path = os.path.join(output_dir, "tmp")
        self.create_dirs()

        self.genetic_code = genetic_code  # list
        self.cd_lower, self.cd_upper = ambiguous_coding_density_range

        self.code_4_selected = False
        self.coding_density = None

        self.log_path = log_path
        self.output_faa = os.path.join(self.output_dir, "predicted.faa")
        self.output_fna = os.path.join(self.output_dir, "predicted.fna")
        self.sequences_aa = {}
        self.sequences_nt = {}
        self.gene_details = defaultdict(list)

        self.gc_run_results = defaultdict(dict)

        self.input_length = self._get_genome_length()
        if self.input_length > 100000:
            self.run_mode = ["single", "meta"]
        else:
            self.run_mode = ["meta"]

    def create_dirs(self):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.exists(self.tmp_path):
            os.makedirs(self.tmp_path)

    def run_prodigal(self, file_id, mode, genetic_code):
        logger.info("Running prodigal with genetic code {} in {} mode".format(genetic_code, mode))  # todo: replace with decorator
        tmp_name = "{}.faa".format(file_id)

        tmp_logfile_out = "{}_out.log".format(file_id)
        tmp_logfile_err = "err".join(tmp_logfile_out.rsplit("out", 1))  # Replace only the last occurence of "out" substring
        self.gc_run_results[genetic_code].update({"tmp_name": tmp_name, "log_file": tmp_logfile_out})
        self.prodigal_tool.logfile_path_out = tmp_logfile_out
        self.prodigal_tool.logfile_path_err = tmp_logfile_err
        self._configure_prodigal_job(mode, genetic_code, tmp_name)

        self.prodigal_tool.run_jobs(1)  # Prodigal only runs on 1 CPU and we're only running one job at a time
        coding_length = self._get_coding_length(tmp_logfile_out)
        coding_density = coding_length / self.input_length
        logger.debug("Coding density is {}".format(coding_density))
        return coding_density

    @log("Genetic code {} selected as optimal", logger, attr_name="gc", on_func_exit=True)
    def run(self):
        """
        1) If genome length > 100000 run in "single" mode, then "meta" mode if there are no gene predictions. Otherwise
        just run in "meta" mode. This is based on the recommendations in the Prodigal docs.
        2) Run once using genetic code 11. This can be overridden if the user includes a spceific genetic code in the
        config file.
        3) Check the genome coding density. If > 0.74 continue using genetic code 11. If <= 0.74 re-run using genetic
        code 4. The cut-off of 0.74 was chosen based on analysis done by Mose Manni.
        4) If the next run still has a genetic density < 0.74, read the stdout log files (formerly the GFF files) and
        extract the scores assigned to each gene prediction. Whichever genetic code yields the greatest mean score is
        selected.
        :return:
        """
        ambiguous_gcs = []
        tmp_files = []
        for ix, m in enumerate(self.run_mode):
            for g in self.genetic_code:
                file_id = os.path.join(self.tmp_path, "prodigal_mode_{0}_code_{1}".format(m, g))
                if os.path.exists("{}.faa".format(file_id)):
                    self.coding_density = self.gc_run_results[g]["cd"]
                else:
                    self.coding_density = self.run_prodigal(file_id, m, g)
                    self.gc_run_results[g].update({"cd": self.coding_density})

                if self.coding_density >= self.cd_upper:
                    self.prodigal_tool.total = 0
                    self.prodigal_tool.nb_done = 0
                    tmp_files.append(self.gc_run_results[g]["tmp_name"])
                    break
                elif self.coding_density > self.cd_lower:
                    ambiguous_gcs.append(g)

                tmp_files.append(self.gc_run_results[g]["tmp_name"])
                self.prodigal_tool.total = 0
                self.prodigal_tool.nb_done = 0

            # If output files from both runs in "single" mode are empty, run again in "meta" mode, else raise Exception.
            if not any([os.stat(tmp_file).st_size > 0 for tmp_file in tmp_files]):
                if ix + 1 == len(self.run_mode):
                    raise NoGenesError("Prodigal")
                else:
                    continue


            if len(tmp_files) == 1:  # todo: rework this if/else clause into a single file selection logic block
                self.gc = self.genetic_code[0]
                selected_tmpfile = tmp_files[0]
                selected_logfile = self.gc_run_results[g]["log_file"]
            else:
                self.gc, selected_tmpfile, selected_logfile = self._select_best_files(ambiguous_gcs=ambiguous_gcs)
            self._organize_prodigal_files(selected_tmpfile, selected_logfile)
            self._get_gene_details()
            self.gc_run_results[self.gc].update({"seqs_aa": self.sequences_aa, "seqs_nt": self.sequences_nt, "gene_details": self.gene_details})
            break
        return

    def select_optimal_results(self, genetic_code, ambiguous_range):
        """
        For selecting optimal results for final BUSCO run.
        :param genetic_code:
        :param ambiguous_range_lower:
        :param ambiguous_range_upper:
        :return:
        """
        ambiguous_range_lower, ambiguous_range_upper = ambiguous_range
        ambiguous_gcs = []
        for gc in genetic_code:
            if self.gc_run_results[gc]["cd"] > ambiguous_range_upper:
                self.gc = gc
                return
            else:
                ambiguous_gcs.append(gc)

        if len(ambiguous_gcs) > 1:
            ambiguous_files = [self.gc_run_results[gc]["log_file"] for gc in genetic_code if gc in ambiguous_gcs]
            self.gc = self._get_optimal_genetic_code(ambiguous_files)
        else:
            self.gc = ambiguous_gcs[0]
        return

    def _configure_prodigal_job(self, prodigal_mode, genetic_code, tmp_name):
        tmp_name_nt = tmp_name.replace("faa", "fna")

        prodigal_job = self.prodigal_tool.create_job()
        prodigal_job.add_parameter("-p")
        prodigal_job.add_parameter("%s" % prodigal_mode)
        prodigal_job.add_parameter("-f")
        prodigal_job.add_parameter("gff")
        prodigal_job.add_parameter("-g")
        prodigal_job.add_parameter("%s" % genetic_code)
        prodigal_job.add_parameter("-a")
        prodigal_job.add_parameter("%s" % tmp_name)
        prodigal_job.add_parameter("-d")
        prodigal_job.add_parameter("%s" % tmp_name_nt)
        prodigal_job.add_parameter("-i")
        prodigal_job.add_parameter("%s" % self.input_file)

        return

    def _get_gene_details(self):

        with open(self.output_fna, "rU") as f:
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

    def _get_optimal_genetic_code(self, files):
        mean_scores = []
        for filename in files:
            scores = []
            with open(filename, "r") as f:
                for line in f:
                    try:
                        score = re.search(";score=(.+?);", line).group(1)
                        scores.append(float(score))
                    except AttributeError:
                        continue
            mean_scores.append(sum(scores) / len(scores))
        optimal_gc = self.genetic_code[mean_scores.index(max(mean_scores))]
        return optimal_gc

    def _organize_prodigal_files(self, tmp_file, tmp_logfile):

        shutil.copy(tmp_file, self.output_faa)
        shutil.copy(tmp_file.replace(".faa", ".fna"), self.output_fna)

        # copy selected log files from tmp/ to logs/
        new_logname = os.path.join(self.log_path, "prodigal_out.log")
        shutil.copy(tmp_logfile, new_logname)
        shutil.copy(tmp_logfile.replace("_out.log", "_err.log"), new_logname.replace("_out.log", "_err.log"))

        return

    def _select_best_files(self, ambiguous_gcs):  # todo: clean up this code and check the logic is all necessary
        if len(ambiguous_gcs) > 1:
            ambiguous_files = [self.gc_run_results[gc]["log_file"] for gc in self.genetic_code if gc in ambiguous_gcs]
            gc = self._get_optimal_genetic_code(ambiguous_files)
        elif len(self.genetic_code) > 1:
            gc = self.genetic_code[-1]  # Select the last run that was unambiguously better with respect to coding density
        else:
            gc = self.genetic_code[0]  # If there is only one file, just take that.

        if gc == "4":
            self.code_4_selected = True
        selected_tmpfile = self.gc_run_results[gc]["tmp_name"]
        selected_logfile = self.gc_run_results[gc]["log_file"]
        return gc, selected_tmpfile, selected_logfile

class NoGenesError(Exception):

    def __init__(self, gene_predictor):
        self.gene_predictor = gene_predictor

class HMMERRunner:

    def __init__(self, hmmer_tool, input_sequences, busco_ids, output_folder, lineage_path, mode, cpus, gene_details, datasets_version):
        self.hmmer = hmmer_tool
        self.input_sequences = input_sequences
        self.busco_ids = busco_ids
        self.lineage_path = lineage_path
        self.output_folder = output_folder
        self.mode = mode
        self.cpus = cpus
        self.datasets_version = datasets_version

        # gene_details can only be None for proteins mode. In the other modes the gene locations are written to a file
        # after the coordinates are loaded from this attribute
        self.gene_details = gene_details

        self.cutoff_dict = defaultdict(dict)
        self.matched_bitscores = defaultdict(list)
        self.matched_genes_complete = defaultdict(list)
        self.matched_genes_vlarge = defaultdict(list)
        self.matched_genes_fragment = defaultdict(list)
        self.is_complete = defaultdict(lambda: defaultdict(list))  # dict of a dict of lists of dicts
        self.is_fragment = defaultdict(lambda: defaultdict(list))
        self.is_very_large = defaultdict(lambda: defaultdict(list))
        self._already_used_genes = set()
        self.single_copy_buscos = {}
        self.multi_copy_buscos = {}
        self.fragmented_buscos = {}
        self.missing_buscos = []
        self.extra_columns = False

    def load_buscos(self):
        """
        Load all BUSCOs for the lineage, along with their cutoff lengths and scores.
        :return:
        """
        self._load_length()
        self._load_score()
        return

    def run(self):
        """
        Create a HMMER job for each BUSCO. Each job searches the input sequence file for matches for the BUSCO gene.
        :return:
        """
        self.hmmer.total = self._count_jobs()
        self.hmmer.count_jobs_created = False
        job_controller = self.generate_jobs()
        jobs_to_run = True
        while jobs_to_run:
            jobs_to_run = next(job_controller)

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

    def generate_jobs(self):
        njobs = 0
        for busco_id in self.busco_ids:
            if busco_id in self.cutoff_dict:
                if isinstance(self.input_sequences, str):
                    output_filename = "{}.out".format(busco_id)
                    self._configure_hmmer_job(busco_id, self.input_sequences, output_filename)
                    njobs += 1
                    if njobs >= 350:
                        self.hmmer.run_jobs(self.cpus)
                        yield True
                        njobs = 0
                elif isinstance(self.input_sequences, list):
                    input_files = [f for f in self.input_sequences if os.path.basename(f).startswith(busco_id)]
                    for seq_filename in input_files:
                        output_filename = os.path.basename(seq_filename).replace("faa", "out")
                        self._configure_hmmer_job(busco_id, seq_filename, output_filename)
                        njobs += 1
                        if njobs >= 350:
                            self.hmmer.run_jobs(self.cpus)
                            yield True
                            njobs = 0
        if njobs > 0:
            self.hmmer.run_jobs(self.cpus)
        yield False

    def process_output(self):
        self._load_matched_genes()
        self._filter()
        self._consolidate_busco_lists()
        return

    def _configure_hmmer_job(self, busco_id, seq_filename, output_filename):
        hmmer_job = self.hmmer.create_job()
        hmmer_job.add_parameter("--domtblout")
        hmmer_job.add_parameter(os.path.join(self.output_folder, output_filename))
        hmmer_job.add_parameter("--cpu")
        hmmer_job.add_parameter("1")
        hmmer_job.add_parameter(os.path.join(self.lineage_path, "hmms", "{}.hmm".format(busco_id)))
        hmmer_job.add_parameter(seq_filename)

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

        with open(os.path.join(self.output_folder, filename), "r") as f:

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
                        if bit_score >= float(self.cutoff_dict[busco_query]["score"]):  # todo: introduce upper bound - consult to see what a reasonable value would be
                            self.matched_bitscores[gene_id].append(bit_score)
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
            busco_type[gene_id].append(dict({"bitscore": max(self.matched_bitscores[gene_id]),
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
        hmmer_results_files = sorted(os.listdir(self.output_folder))

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

        # empty_buscos = []

        for busco_id in list(busco_dict.keys()):
            matches = busco_dict[busco_id]
            # Remove any buscos that appear in higher ranking dictionaries
            if busco_id in higher_rank_buscos:
                busco_dict.pop(busco_id)
                for gene_id in matches:
                    matched_genes[gene_id].remove(busco_id)
                continue

            # Remove any genes that have previously been processed under a different and higher ranking busco match
            for gene_id in list(matches.keys()):
                if gene_id in self._already_used_genes:
                    busco_dict[busco_id].pop(gene_id)
                    matched_genes[gene_id].remove(busco_id)
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

                best_match_ind = max(range(len(busco_bitscores)), key=busco_bitscores.__getitem__)
                buscos.remove(busco_matches[best_match_ind])
                # Remove lower scoring duplicates from dictionary.
                # Note for future development: the matched_genes dictionary is not updated in this method when
                # duplicates are removed from busco_dict
                for duplicate in buscos:
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
        links_file = os.path.join(self.lineage_path, "links_to_{}.txt".format(self.datasets_version.upper()))
        if os.path.exists(links_file):
            with open(links_file, newline='') as f:
                contents = csv.reader(f, delimiter="\t")
                for row in contents:
                    busco_id, description, link = row
                    links_info[busco_id]["description"] = description
                    links_info[busco_id]["link"] = link
        return links_info


    def _format_output_lines(self, busco_dict):
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
                            output_lines.append("{}\tComplete\t{}\t{}\t{}\t{}\t{}\n".format(busco, gene_id, bit_score,
                                                                                            match_length, link, desc))
                        except KeyError:
                            output_lines.append("{}\tComplete\t{}\t{}\t{}\n".format(busco, gene_id,bit_score,
                                                                                    match_length))
                    elif self.mode == "genome":
                        scaffold = self.gene_details[gene_id][m]
                        try:
                            desc = links_info[busco]["description"]
                            link = links_info[busco]["link"]
                            output_lines.append("{}\tComplete\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                busco, gene_id, scaffold["gene_start"], scaffold["gene_end"], bit_score, match_length,
                                link, desc))
                        except KeyError:
                            output_lines.append("{}\tComplete\t{}\t{}\t{}\t{}\t{}\n".format(
                                busco, gene_id, scaffold["gene_start"], scaffold["gene_end"], bit_score, match_length))
        return output_lines

    def _create_output_content(self):
        """
        Format output for all BUSCO matches.
        :return: output_lines
        :rtype: list
        """
        output_lines = []
        for busco_dict in [self.single_copy_buscos, self.multi_copy_buscos, self.fragmented_buscos]:
            output_lines += self._format_output_lines(busco_dict)

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
        lengths_cutoff_file = os.path.join(self.lineage_path, "lengths_cutoff")
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
        scores_cutoff_file = os.path.join(self.lineage_path, "scores_cutoff")
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


class MKBLASTRunner:

    def __init__(self, mkblast_tool, input_file, output_folder, cpus):
        self.mkblast_tool = mkblast_tool
        self.input_file = input_file
        self.db_path = os.path.join(output_folder, "blast_db")
        self.create_dirs()
        self.cpus = cpus
        self.output_db = os.path.join(self.db_path, os.path.basename(self.input_file))

    def create_dirs(self):
        if not os.path.exists(self.db_path):
            os.makedirs(self.db_path)

    @log("Creating BLAST database with input file", logger)
    def configure_mkblast_job(self):
        blast_job = self.mkblast_tool.create_job()
        blast_job.add_parameter("-in")
        blast_job.add_parameter(self.input_file)
        blast_job.add_parameter("-dbtype")
        blast_job.add_parameter("nucl")
        blast_job.add_parameter("-out")
        blast_job.add_parameter(self.output_db)
        return

    def run(self):
        if os.path.exists(self.db_path) and len(os.listdir(self.db_path)) > 0:
            return
        self.configure_mkblast_job()
        self.mkblast_tool.run_jobs(self.cpus)


class TBLASTNRunner:

    def __init__(self, tblastn_tool, input_file, output_folder, lineage_dataset, blast_db, e_v_cutoff, cpus,
                 region_limit, flank, missing_and_frag_only, ancestral_variants, incomplete_buscos):
        self.tblastn_tool = tblastn_tool
        self.input_file = input_file
        self.output_folder = os.path.join(output_folder, "blast_output")
        self.output_seqs = os.path.join(self.output_folder, "sequences")
        self.create_dirs()
        self.lineage_dataset = lineage_dataset
        self.blast_db = blast_db
        if not len(os.listdir(os.path.split(self.blast_db)[0])) > 0:
            raise SystemExit("DB {} not found for tblastn job".format(self.blast_db))
        self.e_v_cutoff = e_v_cutoff
        self.cpus = cpus

        self.incomplete_buscos = incomplete_buscos

        self.ancestral_variants = ancestral_variants
        self.ancestral_sfx = "_variants" if self.ancestral_variants else ""
        self.ancestral_file = os.path.join(self.lineage_dataset, "ancestral{}".format(self.ancestral_sfx))
        self.missing_and_frag_only = missing_and_frag_only
        self.output_suffix = "_missing_and_frag_rerun" if self.missing_and_frag_only else ""
        self.query_file = os.path.join(self.lineage_dataset, "ancestral{}".format(self.ancestral_sfx))
        self.rerun_query_file = os.path.join(self.output_folder, "ancestral{}{}".format(self.ancestral_sfx, self.output_suffix))
        if self.missing_and_frag_only and self.ancestral_variants:
            self._extract_incomplete_buscos_ancestral()

        self.blast_filename = os.path.join(self.output_folder, "tblastn{}.tsv".format(self.output_suffix))
        self.coords_filename = os.path.join(self.output_folder, "coordinates{}.tsv".format(self.output_suffix))

        self.region_limit = region_limit
        self.flank = flank

    def create_dirs(self):
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        if not os.path.exists(self.output_seqs):
            os.makedirs(self.output_seqs)

    def run(self):
        self.tblastn_tool.total = 0
        self.tblastn_tool.nb_done = 0
        self._configure_tblastn_job()
        self.tblastn_tool.run_jobs(self.cpus)  # tblastn manages available cpus by itself
        self._check_output()
        return

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

    def _configure_tblastn_job(self):
        tblastn_job = self.tblastn_tool.create_job()
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
        return

    def _extract_incomplete_buscos_ancestral(self):

        logger.info("Extracting missing and fragmented buscos from the file {}...".format(
            os.path.basename(self.ancestral_file)))

        matched_seqs = []
        busco_ids_retrieved = set()
        with open(self.ancestral_file, "rU") as anc_file:

            for record in SeqIO.parse(anc_file, "fasta"):
                if any(b in record.id for b in self.incomplete_buscos):
                    # Remove the ancestral variant identifier ("_1" etc) so it matches all other BUSCO IDs.
                    # The identifier is still present in the "name" and "description" Sequence Record attributes.
                    record.id = record.id.split("_")[0]
                    logger.debug("Found contig {}".format(record.id))
                    busco_ids_retrieved.add(record.id)
                    matched_seqs.append(record)

        unmatched_incomplete_buscos = list(set(self.incomplete_buscos) - set(busco_ids_retrieved))
        if len(unmatched_incomplete_buscos) > 0:
            logger.warning("The BUSCO ID(s) {} were not found in the file {}".format(
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

    def _get_coordinates(self):
        coords = self._parse_blast_output()
        if self.ancestral_variants:
            coords = self._select_busco_variants(coords)
        coords = self._prune(coords)
        return coords

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

        return coords

    @staticmethod
    def _select_busco_variants(coords):
        """
        Filter contig matches to prevent multiple BUSCO variants matching the same contig.
        The current behaviour combines all contig matches for all BUSCO variants, as long as the contig matches are
        different. There is an open question over whether or not we should only return the contig matches for a single
        BUSCO variant instead of all of them combined. This should only be an issue for the Transcriptome mode.
        :param coords:
        :return:
        """
        selected_coords = defaultdict(lambda: defaultdict(defaultdict))
        for busco_name, contigs in coords.items():
            busco_basename = busco_name.split("_")[0]
            if busco_basename in selected_coords:
                for contig_id in contigs:
                    if contig_id in selected_coords[busco_basename]:
                        if contigs[contig_id]["blast_eval"] < selected_coords[busco_basename][contig_id]["blast_eval"]:
                            selected_coords[busco_basename][contig_id] = contigs[contig_id]
                    else:
                        selected_coords[busco_basename][contig_id] = contigs[contig_id]
                        # logger.warning("Two different genes were matched by two BUSCO ancestral variants")
            else:
                selected_coords[busco_basename] = contigs

        return selected_coords

    def _prune(self, coords):
        for busco_name, contigs in coords.items():
            if len(contigs) > self.region_limit:
                # Sort by blast eval, then isolate smallest values leaving just "region_limit" number of contigs per
                # busco_name
                contigs_to_remove = sorted(
                    contigs, key=lambda contig: contigs[contig]["blast_eval"])[self.region_limit:]
                for c in contigs_to_remove:
                    coords[busco_name].pop(c)
        return coords

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

    def _filter_best_matches(self, coords):

        # Get a list of all start and stop positions of possible busco locations, merging overlapping regions
        for busco_group in coords:
            candidate_contigs = list(coords[busco_group].keys())
            size_lists = self._get_largest_regions(candidate_contigs, coords, busco_group)
            max_size = max(size_lists)  # Get largest match size for a busco group
            # Include all location matches for a busco as long as they are within 70% of the maximum size match
            size_cutoff = int(0.7 * max_size)
            for c, contig_name in enumerate(candidate_contigs):
                if size_lists[c] < size_cutoff:
                    coords[busco_group].pop(contig_name)
        return coords

    def _write_coordinates_to_file(self, coords):

        with open(self.coords_filename, "w") as out:
            for busco_group, contig_matches in coords.items():
                for contig_name in contig_matches:
                    coords[busco_group][contig_name]["contig_start"] = \
                        max(int(coords[busco_group][contig_name]["contig_start"]) - self.flank, 0)
                    contig_start = coords[busco_group][contig_name]["contig_start"]
                    coords[busco_group][contig_name]["contig_end"] += self.flank
                    contig_end = int(coords[busco_group][contig_name]["contig_end"])
                    out.write("{}\t{}\t{}\t{}\n".format(busco_group, contig_name, contig_start, contig_end))
        return

    def _write_contigs(self, coords):
        # Extract all contig identifiers
        contig_names = []
        for contig_info in coords.values():
            for contig in contig_info:
                contig_names.append(contig)

        # Write sequences that match contig ids
        with open(self.input_file, "rU") as f:
            for record in SeqIO.parse(f, "fasta"):
                if record.id in list(set(contig_names)):
                    with open(os.path.join(self.output_seqs, "{}.temp".format(record.id)), "w") as out:
                        SeqIO.write(record, out, "fasta")
        return


class AugustusRunner:

    def __init__(self, augustus_tool, output_folder, seqs_path, target_species, lineage_dataset, params, coords, cpus,
                 log_path, sequences_aa, sequences_nt):
        self.augustus_tool = augustus_tool
        self.output_folder = output_folder
        self.seqs_path = seqs_path
        self.tmp_dir = os.path.join(self.output_folder, "tmp")
        self.target_species = target_species
        self.lineage_dataset = lineage_dataset
        self.params = params
        self.coords = coords
        self.cpus = cpus

        self.extracted_prot_dir = os.path.join(self.output_folder, "extracted_proteins")
        self.pred_genes_dir = os.path.join(self.output_folder, "predicted_genes")
        self.gff_dir = os.path.join(self.output_folder, "gff")
        self.create_dirs()

        self.output_sequences = []
        self.sequences_aa = sequences_aa
        self.sequences_nt = sequences_nt
        self.gene_details = defaultdict(list)
        self.err_logfiles = []
        self.log_path = log_path
        self.any_gene_found = False

    def create_dirs(self):
        if not os.path.exists(self.extracted_prot_dir):
            os.makedirs(self.extracted_prot_dir)
        if not os.path.exists(self.pred_genes_dir):
            os.makedirs(self.pred_genes_dir)
        if not os.path.exists(self.gff_dir):
            os.makedirs(self.gff_dir)
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        return

    def run(self):
        # Todo: refactor logger calls into decorator pattern
        logger.info("Running Augustus prediction using {} as species:".format(self.target_species))
        if self.params:
            logger.info("Additional parameters for Augustus are {}: ".format(self.params))

        self.augustus_tool.total = self._count_jobs()
        self.augustus_tool.count_jobs_created = False
        job_controller = self.generate_jobs()
        jobs_to_run = True
        while jobs_to_run:
            jobs_to_run = next(job_controller)

        self._merge_stderr_logs()

        logger.info("Extracting predicted proteins...")
        files = [f for f in sorted(os.listdir(self.pred_genes_dir)) if any(busco_id in f for busco_id in self.coords)]
        for filename in files:
            self._extract_genes_from_augustus_output(filename)

        if not self.any_gene_found:
            raise NoGenesError("Augustus")

        return

    def _count_jobs(self):
        n = 0
        for busco_group, contigs in self.coords.items():
            for _ in contigs:
                n += 1
        return n

    def generate_jobs(self):
        contig_ordinal_inds = defaultdict(int)
        njobs = 0
        for busco_group, contigs in self.coords.items():

            for contig_name, contig_info in contigs.items():
                contig_tmp_file = "{}.temp".format(contig_name[:100])  # Avoid very long filenames
                contig_start = contig_info["contig_start"]
                contig_end = contig_info["contig_end"]
                contig_ordinal_inds[busco_group] += 1
                output_index = contig_ordinal_inds[busco_group]
                out_filename = os.path.join(self.pred_genes_dir, "{}.out.{}".format(busco_group, output_index))
                self._configure_augustus_job(busco_group, contig_tmp_file, contig_start, contig_end, out_filename)
                njobs += 1
                if njobs >= 350:
                    self.augustus_tool.run_jobs(self.cpus)
                    yield True
                    njobs = 0
        if njobs > 0:
            self.augustus_tool.run_jobs(self.cpus)
        yield False


    def _merge_stderr_logs(self):
        with open(os.path.join(self.log_path, "augustus_err.log"), "a") as f:
            for err_logfile in self.err_logfiles:
                with open(err_logfile, "r") as g:
                    content = g.readlines()
                    f.writelines(content)
        return

    def _configure_augustus_job(self, busco_group, contig_tmp_file, contig_start, contig_end, out_filename):
        # Augustus does not provide an option to write to an output file, so have to change the pipe target from the
        # log file to the desired output file
        self.augustus_tool.logfile_path_out = out_filename
        err_logfile = os.path.join(self.tmp_dir, os.path.basename(out_filename.replace("out", "err")))
        self.augustus_tool.logfile_path_err = err_logfile
        self.err_logfiles.append(err_logfile)

        augustus_job = self.augustus_tool.create_job()
        augustus_job.add_parameter("--codingseq=1")
        augustus_job.add_parameter("--proteinprofile={}".format(os.path.join(self.lineage_dataset,
                                                                             "prfl",
                                                                             "{}.prfl".format(busco_group))))
        augustus_job.add_parameter("--predictionStart={}".format(contig_start))
        augustus_job.add_parameter("--predictionEnd={}".format(contig_end))
        augustus_job.add_parameter("--species={}".format(self.target_species))
        for p in self.params.split():
            if len(p) > 2:
                augustus_job.add_parameter(p)
        augustus_job.add_parameter(os.path.join(self.seqs_path, contig_tmp_file))
        return

    def _extract_genes_from_augustus_output(self, filename):  # todo: consider parallelizing this and other parsing functions

        gene_id = None
        gene_info = []
        gff_filename = os.path.join(self.gff_dir, filename.split(".")[0] + ".gff")
        sequences_aa = []
        sequences_nt = []
        gene_found = False
        completed_record = False

        with open(os.path.join(self.pred_genes_dir, filename), "r", encoding="utf-8") as f:
            # utf-8 encoding needed to handle the umlaut in the third line of the file.
            gene_info_section = False
            nt_sequence_section = False
            aa_sequence_section = False
            nt_sequence_parts = []
            aa_sequence_parts = []

            for line in f:

                if line.startswith("# end gene"):
                    aa_sequence_section = False
                    completed_record = True
                    if gene_id is not None:
                        aa_sequence = "".join(aa_sequence_parts)
                        nt_sequence = "".join(nt_sequence_parts)
                        seq_record_aa = SeqRecord(Seq(aa_sequence.upper(), IUPAC.protein), id=gene_id)
                        seq_record_nt = SeqRecord(Seq(nt_sequence.upper(), IUPAC.unambiguous_dna), id=gene_id)
                        sequences_aa.append(seq_record_aa)
                        sequences_nt.append(seq_record_nt)
                        aa_sequence_parts = []
                        nt_sequence_parts = []
                        # self.aa_seq_records.append(seq_record_aa)
                        # self.nt_seq_records.append(seq_record_nt)
                        gene_id = None
                    continue

                if aa_sequence_section and line.startswith("# sequence of block"):
                    aa_sequence_section = False
                    continue

                if aa_sequence_section:
                    line = line.strip().lstrip("# ").rstrip("]")
                    aa_sequence_parts.append(line)
                    continue

                if line.startswith("# protein"):
                    nt_sequence_section = False
                    aa_sequence_section = True
                    line = line.strip().rstrip("]").split("[")
                    aa_sequence_parts.append(line[1])
                    continue

                if nt_sequence_section:
                    line = line.strip().lstrip("# ").rstrip("]")
                    nt_sequence_parts.append(line)
                    continue

                if line.startswith("# coding sequence"):
                    with open(gff_filename, "a") as g:
                        g.write("\n".join(gene_info) + "\n")
                    gene_info_section = False
                    nt_sequence_section = True
                    line = line.strip().rstrip("]").split("[")  # Extract sequence part of line
                    nt_sequence_parts.append(line[1])
                    continue

                if gene_info_section:
                    line = line.strip().split()
                    seq_name = line[0]
                    gene_start = line[3]
                    gene_end = line[4]
                    if not gene_id:
                        gene_id = "{}:{}-{}".format(seq_name, gene_start, gene_end)
                        self.gene_details[gene_id].append({"gene_start": gene_start, "gene_end": gene_end})
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

    def _write_sequences_to_file(self, filename, sequences_nt, sequences_aa):

        output_fna = os.path.join(self.extracted_prot_dir, filename.replace("out", "fna"))
        output_faa = os.path.join(self.extracted_prot_dir, filename.replace("out", "faa"))
        self.output_sequences.append(output_faa)

        with open(output_fna, "w") as out_fna:
            SeqIO.write(sequences_nt, out_fna, "fasta")
        with open(output_faa, "w") as out_faa:
            SeqIO.write(sequences_aa, out_faa, "fasta")

        return


class GFF2GBRunner:

    def __init__(self, gff2gbSmallDNA_tool, output_folder, input_file, single_copy_buscos, cpus):
        self.gff2gbSmallDNA_tool = gff2gbSmallDNA_tool
        self.output_folder = output_folder
        self.input_file = input_file
        self.single_copy_buscos = single_copy_buscos
        self.cpus = cpus
        self.gff_folder = os.path.join(self.output_folder, "augustus_output", "gff")
        self.gb_folder = os.path.join(self.output_folder, "augustus_output", "gb")
        self.create_dirs()

    def create_dirs(self):
        if not os.path.exists(self.gb_folder):
            os.makedirs(self.gb_folder)
        if not os.path.exists(self.gff_folder):
            os.makedirs(self.gff_folder)

    def run(self):
        self.gff2gbSmallDNA_tool.total = self._count_jobs()
        self.gff2gbSmallDNA_tool.count_jobs_created = False
        job_controller = self.generate_jobs()
        jobs_to_run = True
        while jobs_to_run:
            jobs_to_run = next(job_controller)

    def _count_jobs(self):
        n = len(self.single_copy_buscos)
        return n

    def generate_jobs(self):
        njobs = 0
        for busco_id in self.single_copy_buscos:
            self._configure_gff2gb_job(busco_id)
            njobs += 1
            if njobs >= 350:
                self.gff2gbSmallDNA_tool.run_jobs(self.cpus)
                yield True
                njobs = 0
        if njobs > 0:
            self.gff2gbSmallDNA_tool.run_jobs(self.cpus)
        yield False

    def _configure_gff2gb_job(self, busco_id):
        gff2_gb_small_dna_pl_job = self.gff2gbSmallDNA_tool.create_job()
        gff2_gb_small_dna_pl_job.add_parameter(os.path.join(self.gff_folder, "{}.gff".format(busco_id)))
        gff2_gb_small_dna_pl_job.add_parameter(self.input_file)
        gff2_gb_small_dna_pl_job.add_parameter("1000")
        gff2_gb_small_dna_pl_job.add_parameter(os.path.join(self.gb_folder, "{}.raw.gb".format(busco_id)))


class NewSpeciesRunner:

    def __init__(self, new_species_tool, domain, new_species_name, cpus):
        self.new_species_tool = new_species_tool
        self.domain = domain
        self.new_species_name = new_species_name
        self.cpus = cpus

    def run(self):
        self._configure_new_species_job()
        self.new_species_tool.run_jobs(self.cpus)

    def _configure_new_species_job(self):

        new_species_pl_job = self.new_species_tool.create_job()
        # bacteria clade needs to be flagged as "prokaryotic"
        if self.domain == "prokaryota":
            new_species_pl_job.add_parameter("--prokaryotic")
        new_species_pl_job.add_parameter("--species={}".format(os.path.basename(self.new_species_name)))


class ETrainingRunner:

    def __init__(self, etraining_tool, output_folder, cpus):
        self.etraining_tool = etraining_tool
        self.output_folder = output_folder
        self.cpus = cpus

    def run(self):
        self._configure_etraining_job()
        self.etraining_tool.run_jobs(self.cpus)

    def _configure_etraining_job(self):
        etraining_job = self.etraining_tool.create_job()
        etraining_job.add_parameter("--species=BUSCO_{}".format(os.path.basename(self.output_folder)))
        etraining_job.add_parameter(os.path.join(self.output_folder, "augustus_output", "training_set.db"))


class OptimizeAugustusRunner:

    def __init__(self, optimize_augustus_tool, output_folder, new_species_name, cpus):
        self.optimize_augustus_tool = optimize_augustus_tool
        self.augustus_output_folder = os.path.join(output_folder, "training_set.db")
        self.new_species_name = new_species_name
        self.cpus = cpus

    def _configure_optimize_augustus_job(self):
        optimize_augustus_pl_job = self.optimize_augustus_tool.create_job()
        optimize_augustus_pl_job.add_parameter("--cpus={}".format(self.cpus))
        optimize_augustus_pl_job.add_parameter("--species={}".format(self.new_species_name))
        optimize_augustus_pl_job.add_parameter(self.augustus_output_folder)

    def run(self):
        self._configure_optimize_augustus_job()
        self.optimize_augustus_tool.run_jobs(self.cpus)


class SEPPRunner:

    def __init__(self, sepp_tool, output_folder, placement_folder, tree_nwk_file, tree_metadata_file, supermatrix_file,
                 downloader, datasets_version, cpus):
        self.sepp_tool = sepp_tool
        self.output_folder = output_folder
        self.placement_folder = placement_folder
        self.tree_nwk_file = tree_nwk_file
        self.tree_metadata_file = tree_metadata_file
        self.supermatrix_file = supermatrix_file
        self.downloader = downloader
        self.datasets_version = datasets_version
        self.cpus = cpus

    def run(self):
        self._configure_sepp_job()
        self.sepp_tool.run_jobs(1)

    def _configure_sepp_job(self):
        sepp_job = self.sepp_tool.create_job()
        sepp_job.add_parameter("--cpu")
        sepp_job.add_parameter(self.cpus)
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
