from busco.busco_tools.base import BaseRunner
import os
from collections import defaultdict
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.Exceptions import BuscoError
from Bio import SeqIO
import subprocess

logger = BuscoLogger.get_logger(__name__)


class MKBLASTRunner(BaseRunner):

    name = "makeblastdb"
    cmd = "makeblastdb"

    def __init__(self):
        super().__init__()
        self.db_path = os.path.join(
            self.config.get("busco_run", "main_out"), "blast_db"
        )
        self.output_db = os.path.join(self.db_path, os.path.basename(self.input_file))
        self.create_dirs(self.db_path)
        self.total = 1
        self.run_number += 1

    def configure_runner(self, *args):
        super().configure_runner(args)

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

    def reset(self):
        super().reset()

    def generate_job_args(self):
        yield

    def get_version(self):
        mkblastdb_version_call = subprocess.check_output(
            [self.cmd, "-version"], stderr=subprocess.STDOUT, shell=False
        )
        mkblastdb_version = ".".join(
            mkblastdb_version_call.decode("utf-8").split("\n")[0].split()[1].rsplit(".")
        )

        return mkblastdb_version

    def check_tool_dependencies(self):
        pass

    @property
    def output_folder(self):
        return self.db_path


class TBLASTNRunner(BaseRunner):

    name = "tblastn"
    cmd = "tblastn"

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

        self.blast_db = None
        self.missing_and_frag_only = None
        self.ancestral_variants = None
        self.incomplete_buscos = None

        self.ancestral_sfx = None
        self.ancestral_file = None
        self.query_file = None
        self.output_suffix = None
        self.rerun_query_file = None
        self.blast_filename = None
        self.coords_filename = None

    def configure_runner(
        self, blast_db, missing_and_frag_only, ancestral_variants, incomplete_buscos
    ):
        super().configure_runner()
        self.run_number += 1
        self.blast_db = blast_db
        self.missing_and_frag_only = missing_and_frag_only
        self.ancestral_variants = ancestral_variants
        self.incomplete_buscos = incomplete_buscos

        self.ancestral_sfx = "_variants" if self.ancestral_variants else ""
        self.ancestral_file = os.path.join(
            self.lineage_dataset, "ancestral{}".format(self.ancestral_sfx)
        )
        self.query_file = os.path.join(
            self.lineage_dataset, "ancestral{}".format(self.ancestral_sfx)
        )
        self.output_suffix = (
            "_missing_and_frag_rerun" if self.missing_and_frag_only else ""
        )
        self.rerun_query_file = os.path.join(
            self._output_folder,
            "ancestral{}{}".format(self.ancestral_sfx, self.output_suffix),
        )
        if self.missing_and_frag_only and self.ancestral_variants:
            self._extract_incomplete_buscos_ancestral()

        self.blast_filename = os.path.join(
            self._output_folder, "tblastn{}.tsv".format(self.output_suffix)
        )
        self.coords_filename = os.path.join(
            self._output_folder, "coordinates{}.tsv".format(self.output_suffix)
        )

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
        try:
            size = os.path.getsize(self.input_file) / 1000  # size in mb
            flank = int(size / 50)  # proportional flank size
            # Ensure value is between 5000 and MAX_FLANK
            flank = min(max(flank, 5000), type(self).MAX_FLANK)
        except IOError:  # Input data is only validated during run_analysis. This will catch any IO issues before that.
            raise BuscoError(
                "Impossible to read the fasta file {}".format(self.input_file)
            )

        return flank

    @log("Running a BLAST search for BUSCOs against created database", logger)
    def run(self):
        super().run()
        self.run_jobs()
        self._check_output()
        return

    def reset(self):
        super().reset()

    def check_tool_dependencies(self):
        if (
            ".".join(self.version.split(".")[:-1]) not in ["2.2", "2.3"]
            and self.version != "2.10.1+"
            and float(".".join(self.version.split(".")[:-1])) < 2.11
        ):
            # Known problems with multithreading on BLAST 2.4-2.10.0.
            logger.warning(
                "You are using BLAST version {}. This is known to yield inconsistent results when "
                "multithreading. BLAST will run on a single core as a result. For performance improvement, "
                "please upgrade to BLAST 2.10.1+.".format(self.version)
            )
            self.cpus = 1

    def get_version(self):
        tblastn_version_call = subprocess.check_output(
            [self.cmd, "-version"], stderr=subprocess.STDOUT, shell=False
        )
        tblastn_version = ".".join(
            tblastn_version_call.decode("utf-8").split("\n")[0].split()[1].rsplit(".")
        )

        return tblastn_version

    def generate_job_args(self):
        yield

    def _check_output(self):
        # check that blast worked
        if not os.path.exists(self.blast_filename):
            raise BuscoError("tblastn failed!")

        # check that the file is not truncated
        with open(self.blast_filename, "r") as f:
            try:
                if "processed" not in f.readlines()[-1]:
                    raise BuscoError(
                        "tblastn has ended prematurely (the result file lacks the expected final line), "
                        "which will produce incomplete results in the next steps ! This problem likely "
                        "appeared in blast+ 2.4 and seems not fully fixed in 2.6. It happens only when "
                        "using multiple cores. You can use a single core (-c 1) or downgrade to "
                        "blast+ 2.2.x, a safe choice regarding this issue. See blast+ documentation for "
                        "more information."
                    )

            except IndexError:
                # if the tblastn result file is empty, for example in phase 2
                # if 100% was found in phase 1
                pass
        return

    def _extract_incomplete_buscos_ancestral(self):

        logger.info(
            "Extracting missing and fragmented buscos from the file {}...".format(
                os.path.basename(self.ancestral_file)
            )
        )

        matched_seqs = []
        busco_ids_retrieved = set()
        with open(self.ancestral_file, "rU") as anc_file:

            for record in SeqIO.parse(anc_file, "fasta"):
                if any(record.id.startswith(b) for b in self.incomplete_buscos):
                    # Remove the ancestral variant identifier ("_1" etc) so it matches all other BUSCO IDs.
                    # The identifier is still present in the "name" and "description" Sequence Record attributes.
                    record.id = record.id.split("_")[0]
                    busco_ids_retrieved.add(record.id)
                    matched_seqs.append(record)

        unmatched_incomplete_buscos = list(
            set(self.incomplete_buscos) - set(busco_ids_retrieved)
        )
        if len(unmatched_incomplete_buscos) > 0:
            logger.debug(
                "The BUSCO ID(s) {} were not found in the file {}".format(
                    unmatched_incomplete_buscos, os.path.basename(self.ancestral_file)
                )
            )

        self.query_file = self.rerun_query_file
        with open(
            self.query_file, "w"
        ) as out_file:  # Create new query file for second tblastn run
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
        coords = defaultdict(
            lambda: defaultdict(defaultdict)
        )  # dict of busco_id -> contig_id -> {info}
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
                        coords[busco_name][contig_id] = {
                            "contig_start": contig_start,
                            "contig_end": contig_end,
                            "busco_coords": [[busco_start, busco_end]],
                            "blast_eval": blast_eval,
                        }

                    elif (
                        contig_id in coords[busco_name]
                    ):  # i.e. if the same gene matched the busco more than once.
                        # now update coordinates
                        coords = self._update_coordinates(
                            coords,
                            busco_name,
                            contig_id,
                            busco_start,
                            busco_end,
                            contig_start,
                            contig_end,
                            blast_eval,
                        )

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
                        if (
                            contigs[contig_id]["blast_eval"]
                            < selected_coords[busco_basename][contig_id]["blast_eval"]
                        ):
                            selected_coords[busco_basename][contig_id] = contigs[
                                contig_id
                            ]
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
                    contigs, key=lambda contig: contigs[contig]["blast_eval"]
                )[self.region_limit :]
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
    def _update_coordinates(
        coords,
        busco_name,
        contig,
        busco_start,
        busco_end,
        contig_start,
        contig_end,
        blast_eval,
    ):
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
        elif (
            coords[busco_name][contig]["contig_end"]
            >= contig_start
            >= coords[busco_name][contig]["contig_start"]
        ):
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
            size_lists = self._get_largest_regions(
                candidate_contigs, self.coords, busco_group
            )
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
                    self.coords[busco_group][contig_name]["contig_start"] = max(
                        int(self.coords[busco_group][contig_name]["contig_start"])
                        - self.flank,
                        0,
                    )
                    contig_start = self.coords[busco_group][contig_name]["contig_start"]
                    self.coords[busco_group][contig_name]["contig_end"] += self.flank
                    contig_end = int(
                        self.coords[busco_group][contig_name]["contig_end"]
                    )
                    out.write(
                        "{}\t{}\t{}\t{}\n".format(
                            busco_group, contig_name, contig_start, contig_end
                        )
                    )
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
                    with open(
                        os.path.join(self.output_seqs, "{}.temp".format(record.id)), "w"
                    ) as out:
                        SeqIO.write(record, out, "fasta")
        return
