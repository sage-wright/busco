# coding: utf-8
"""
GenomeAnalysis.py

Contains classes controlling all genome mode pipelines.

Author(s): Matthew Berkeley, Mathieu Seppey, Mose Manni, Felipe Simao, Rob Waterhouse

Copyright (c) 2015-2024, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

from busco.analysis.BuscoAnalysis import BuscoAnalysis
from busco.analysis.Analysis import NucleotideAnalysis, BLASTAnalysis
from busco.busco_tools.prodigal import ProdigalRunner
from busco.busco_tools.metaeuk import MetaeukRunner
from busco.busco_tools.miniprot import MiniprotIndexRunner, MiniprotAlignRunner
from busco.busco_tools.bbtools import BBToolsRunner
from busco.busco_tools.augustus import (
    AugustusRunner,
    GFF2GBRunner,
    NewSpeciesRunner,
    ETrainingRunner,
    OptimizeAugustusRunner,
)
from busco.busco_tools.base import NoRerunFile, NoGenesError
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from abc import ABCMeta, abstractmethod
from configparser import NoOptionError
import time
import os
import glob
from collections import defaultdict
from busco.Exceptions import BuscoError
import numpy as np

logger = BuscoLogger.get_logger(__name__)


class GenomeAnalysis(NucleotideAnalysis, BuscoAnalysis, metaclass=ABCMeta):
    _mode = "genome"
    _bbtools_already_run = False

    def __init__(self):
        super().__init__()
        self.bbtools_runner = None

    @abstractmethod
    def run_analysis(self):
        super().run_analysis()
        if "genome" in self.config.get("busco_run", "mode"):
            # This avoids the child class euk_tran mode running bbtools
            if not self.config.getboolean("busco_run", "skip_bbtools"):
                if GenomeAnalysis._bbtools_already_run:
                    logger.info("Skipping BBTools as already run")
                else:
                    self._run_bbtools()
            else:
                logger.info("Skipping BBTools run as requested")

    def init_tools(self):
        """
        Initialize tools needed for Genome Analysis.
        :return:
        """
        super().init_tools()
        if "genome" in self.config.get("busco_run", "mode"):
            # This avoids the child class euk_tran mode running bbtools
            self.bbtools_runner = BBToolsRunner()

    def _run_bbtools(self):
        self.bbtools_runner.configure_runner()
        if self.bbtools_runner.check_previous_completed_run():
            logger.info("Skipping BBTools run as already run")
        else:
            self.bbtools_runner.run()
        GenomeAnalysis._bbtools_already_run = (
            True  # needed outside the if/else block for restart mode on incomplete runs
        )
        if len(self.bbtools_runner.metrics) == 0:
            self.bbtools_runner.parse_output()

    def reset(self):
        super().reset()
        if self.bbtools_runner:
            self.bbtools_runner.reset()
            GenomeAnalysis._bbtools_already_run = False


class GenomeAnalysisProkaryotes(GenomeAnalysis):
    """
    This class runs a BUSCO analysis on a genome.
    """

    def __init__(self):
        """
        Initialize an instance.
        """
        super().__init__()
        self.prodigal_runner = None

    def cleanup(self):
        super().cleanup()

    def run_analysis(self):
        """
        This function calls all needed steps for running the analysis.
        """
        super().run_analysis()
        self._run_prodigal()
        self.gene_details = self.prodigal_runner.gene_details
        self.run_hmmer(self.prodigal_runner.output_faa)
        self.hmmer_runner.write_buscos_to_file()
        return

    def init_tools(self):
        """
        Init the tools needed for the analysis
        """
        super().init_tools()
        self.prodigal_runner = ProdigalRunner()

    def reset(self):
        super().reset()
        if (
            self.prodigal_runner
        ):  # If final run has already been run, then the prodigal_runner object in the final runner object will
            # still be set to None
            self.prodigal_runner.reset()

    @log("***** Run Prodigal on input to predict and extract genes *****", logger)
    def _run_prodigal(self):
        """
        Run Prodigal on input file to detect genes.
        :return:
        """
        self.prodigal_runner.configure_runner()
        if self.restart and self.prodigal_runner.check_previous_completed_run():
            logger.info("Skipping Prodigal run as it has already completed")
            self.prodigal_runner.run(restart=self.restart)
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.prodigal_runner.run()
        return


class GenomeAnalysisEukaryotes(GenomeAnalysis):
    """
    This class runs a BUSCO analysis on a eukaryote genome.
    """

    def __init__(self):
        super().__init__()
        self.gene_details = {}
        self.gene_update_mapping = defaultdict(dict)

    def cleanup(self):
        """
        This function cleans temporary files
        """
        super().cleanup()

    def init_tools(self):
        """
        Initialize all required tools for Genome Eukaryote Analysis:
        metaeuk
        :return:
        """
        super().init_tools()

        return

    def reset(self):
        super().reset()

    @abstractmethod
    def run_analysis(self):
        super().run_analysis()

    def validate_output(self):
        """Need to run this for both initial and rerun because it is possible that metaeuk matches overlap"""

        hmmer_results = self.hmmer_runner.merge_dicts()

        if len(hmmer_results) > 0:
            exon_records = self.get_exon_records(hmmer_results)
            overlaps = self.find_overlaps(exon_records)
            if len(overlaps) > 0:
                entries_to_remove = self.handle_overlaps(
                    overlaps, exon_records, hmmer_results
                )
                entries_to_remove = np.array(
                    entries_to_remove, dtype=exon_records.dtype
                )
                complete, matched_genes_complete = self.reconstruct_hmmer_results(
                    entries_to_remove, exon_records, self.hmmer_runner.is_complete
                )
                v_large, matched_genes_v_large = self.reconstruct_hmmer_results(
                    entries_to_remove, exon_records, self.hmmer_runner.is_very_large
                )
                fragmented, matched_genes_fragmented = self.reconstruct_hmmer_results(
                    entries_to_remove, exon_records, self.hmmer_runner.is_fragment
                )

                # Update hmmer runner with new dictionaries
                self.hmmer_runner.is_complete = complete
                self.hmmer_runner.is_very_large = v_large
                self.hmmer_runner.is_fragment = fragmented
                self.hmmer_runner.matched_genes_complete = matched_genes_complete
                self.hmmer_runner.matched_genes_vlarge = matched_genes_v_large
                self.hmmer_runner.matched_genes_fragment = matched_genes_fragmented
                self.hmmer_runner.gene_details = self.gene_details
        return

    def get_exon_records(
        self, busco_dict
    ):  # Placed in the GenomeAnalysis module because it draws on both hmmer_runner and metaeuk_runner methods

        matched_records = {}
        exon_records = []
        for busco_id, gene_match in busco_dict.items():
            for gene_id, details in gene_match.items():
                matched_records[gene_id] = self.gene_details[details[0]["ref gene ID"]]
                exon_coords = matched_records[gene_id]["exon_coords"]
                for entry in exon_coords:
                    taken_coord1 = int(entry[0])
                    taken_coord2 = int(entry[1])
                    low_coord = min([taken_coord1, taken_coord2])
                    high_coord = max([taken_coord1, taken_coord2])
                    exon_records.append(
                        (
                            details[0]["ref gene ID"],
                            busco_id,
                            matched_records[gene_id]["contig_id"],
                            low_coord,
                            high_coord,
                            matched_records[gene_id]["strand"],
                            details[0]["score"],
                            matched_records[gene_id]["run_number"],
                        )
                    )

        exon_records_arr = np.array(
            exon_records,
            dtype=[
                ("gene_id", "U100"),
                ("target_id", "U30"),
                ("contig_id", "U70"),
                ("low_coord", "i4"),
                ("high_coord", "i4"),
                ("strand", "U1"),
                ("score", "f4"),
                ("run_number", "i4"),
            ],
        )
        exon_records_arr_sorted = np.sort(exon_records_arr, order=["low_coord"])
        return exon_records_arr_sorted

    def reconstruct_hmmer_results(
        self, entries_to_remove, exon_records, hmmer_result_dict
    ):

        hmmer_result_dict_new = defaultdict(dict)
        matched_genes_new = defaultdict(list)
        for busco_id, matches in hmmer_result_dict.items():
            try:
                busco_group = exon_records[exon_records["target_id"] == busco_id]
            except KeyError:  # if busco was removed during overlap filtering
                continue
            for gene_id, details in matches.items():
                orig_gene_id = details[0]["ref gene ID"]
                busco_gene_group = busco_group[busco_group["gene_id"] == orig_gene_id]
                if gene_id not in matches:
                    if busco_id not in self.gene_update_mapping:
                        continue
                    elif orig_gene_id not in self.gene_update_mapping[busco_id]:
                        continue
                    else:
                        orig_gene_id = self.gene_update_mapping[busco_id][orig_gene_id][
                            "new_gene_coords"
                        ]
                new_gene_start, new_gene_stop = gene_id.rsplit(":", 1)[-1].split(
                    "-"
                )  # this line just initializes values that will be changed if an exon is removed
                new_gene_low = new_gene_start  # initialize values
                new_gene_high = new_gene_stop
                start_trim = 0
                end_trim = 0
                intersection = np.intersect1d(busco_gene_group, entries_to_remove)
                if len(intersection) > 0:
                    if np.array_equal(
                        np.sort(intersection, order=["low_coord"]),
                        np.sort(busco_gene_group, order=["low_coord"]),
                    ):  # if intersection == busco_gene_group
                        continue  # remove entire gene - don't add to new dict
                    ordered_exons = np.sort(busco_gene_group, order=["low_coord"])
                    strand = ordered_exons[0]["strand"]

                    for exon in ordered_exons:
                        if exon in intersection:
                            if strand == "+":
                                start_trim += exon["high_coord"] - exon["low_coord"] + 1
                            else:
                                end_trim += exon["high_coord"] - exon["low_coord"] + 1
                        else:
                            new_gene_low = exon["low_coord"]
                            break

                    for exon in ordered_exons[::-1]:
                        if exon in intersection:
                            if strand == "+":
                                end_trim += exon["high_coord"] - exon["low_coord"] + 1
                            else:
                                start_trim += exon["high_coord"] - exon["low_coord"] + 1
                        else:
                            new_gene_high = exon["high_coord"]
                            break

                    new_long_gene_id = "{}:{}-{}".format(
                        orig_gene_id.rsplit(":", 1)[0], new_gene_low, new_gene_high
                    )
                    new_short_gene_id = "{}:{}-{}".format(
                        gene_id.rsplit(":", 1)[0], new_gene_low, new_gene_high
                    )
                else:
                    new_long_gene_id = orig_gene_id
                    new_short_gene_id = gene_id

                details = matches[gene_id]
                details[0]["ref gene ID"] = new_long_gene_id
                details[0]["env_coords"] = [
                    (x[0] - int(start_trim / 3), x[1] - int(start_trim / 3))
                    for x in details[0]["env_coords"]
                ]
                hmmer_result_dict_new[busco_id].update({new_short_gene_id: details})
                matched_genes_new[new_short_gene_id].append(busco_id)

                if new_short_gene_id != gene_id:

                    self.gene_details[new_long_gene_id] = self.gene_details[
                        orig_gene_id
                    ].copy()
                    self.gene_details[new_long_gene_id]["contig_start"] = new_gene_low
                    self.gene_details[new_long_gene_id]["contig_end"] = new_gene_high
                    self.gene_details[new_long_gene_id]["exon_coords"] = []
                    for exon_coords in self.gene_details[orig_gene_id]["exon_coords"]:
                        old_low = min((int(exon_coords[0]), int(exon_coords[1])))
                        old_high = max((int(exon_coords[0]), int(exon_coords[1])))
                        if old_low >= int(new_gene_low) and old_high <= int(
                            new_gene_high
                        ):
                            self.gene_details[new_long_gene_id]["exon_coords"].append(
                                exon_coords
                            )

                    trimmed_sequence_aa, trimmed_sequence_nt = self.trim_sequence(
                        orig_gene_id, new_long_gene_id, start_trim, end_trim
                    )
                    self.gene_update_mapping[busco_id][gene_id] = {
                        "new_gene_coords": new_long_gene_id,
                        "start_trim": start_trim,
                        "end_trim": end_trim,
                    }
                else:
                    trimmed_sequence_aa = self.gene_details[orig_gene_id]["aa_seq"]
                    try:
                        trimmed_sequence_nt = self.gene_details[orig_gene_id]["nt_seq"]
                    except KeyError:
                        trimmed_sequence_nt = None
                self.gene_details[new_long_gene_id]["aa_seq"] = trimmed_sequence_aa
                if trimmed_sequence_nt:
                    self.gene_details[new_long_gene_id]["nt_seq"] = trimmed_sequence_nt
        return hmmer_result_dict_new, matched_genes_new

    def trim_sequence(self, old_gene_match, new_gene_match, start_trim, end_trim):
        old_sequence_aa = self.gene_details[old_gene_match]["aa_seq"]
        try:
            old_sequence_nt = self.gene_details[old_gene_match]["nt_seq"]

            new_sequence_nt = old_sequence_nt[
                start_trim : len(old_sequence_nt) - end_trim
            ]
            new_sequence_nt.id = (
                new_sequence_nt.name
            ) = new_sequence_nt.description = new_gene_match
        except KeyError:
            new_sequence_nt = None

        new_sequence_aa = old_sequence_aa[
            int(start_trim / 3) : len(old_sequence_aa) - int(end_trim / 3)
        ]
        new_sequence_aa.id = (
            new_sequence_aa.name
        ) = new_sequence_aa.description = new_gene_match
        return new_sequence_aa, new_sequence_nt

    @abstractmethod
    def run_overlap_finder(self, *args):
        return []

    def find_overlaps(self, exon_records_arr):
        overlaps = self.run_overlap_finder(exon_records_arr)
        logger.info("{} candidate overlapping regions found".format(len(overlaps)))
        logger.info("{} exons in total".format(len(exon_records_arr)))
        return overlaps

    def handle_overlaps(self, overlaps, exon_records, hmmer_results):
        entries_to_remove = []
        for overlap_entry in overlaps:
            if (
                overlap_entry[0] in entries_to_remove
                or overlap_entry[1] in entries_to_remove
            ):
                continue
            else:
                entries_to_remove.extend(
                    self.handle_diff_busco_overlap(
                        overlap_entry, exon_records, hmmer_results
                    )
                )
        return entries_to_remove

    def handle_diff_busco_overlap(self, overlap, exon_records, hmmer_results):

        match1, match2 = overlap
        exons1 = exon_records[exon_records["gene_id"] == match1["gene_id"]]
        exons2 = exon_records[exon_records["gene_id"] == match2["gene_id"]]
        busco_match1 = match1["target_id"]
        busco_match2 = match2["target_id"]
        gene_id1 = match1["gene_id"].split("|", maxsplit=1)[1]
        gene_id2 = match2["gene_id"].split("|", maxsplit=1)[1]
        strand1 = match1["strand"]
        strand2 = match2["strand"]

        hmmer_match_details1 = hmmer_results[busco_match1][gene_id1]
        hmmer_match_details2 = hmmer_results[busco_match2][gene_id2]

        if (
            busco_match1 in self.gene_update_mapping
            and gene_id1 in self.gene_update_mapping[busco_match1]
        ):
            start_trim1 = self.gene_update_mapping[busco_match1][gene_id1]["start_trim"]
            end_trim1 = self.gene_update_mapping[busco_match1][gene_id1]["end_trim"]
        else:
            start_trim1 = 0
            end_trim1 = 0

        if (
            busco_match2 in self.gene_update_mapping
            and gene_id2 in self.gene_update_mapping[busco_match2]
        ):
            start_trim2 = self.gene_update_mapping[busco_match2][gene_id2]["start_trim"]
            end_trim2 = self.gene_update_mapping[busco_match2][gene_id2]["end_trim"]
        else:
            start_trim2 = 0
            end_trim2 = 0

        if hmmer_match_details1[0]["score"] > hmmer_match_details2[0]["score"]:
            priority_match = hmmer_match_details1
            secondary_match = hmmer_match_details2
            priority_exons = np.sort(exons1, order=["low_coord"])
            if strand1 == "-":
                priority_exons = priority_exons[::-1]

            secondary_exons = np.sort(exons2, order=["low_coord"])
            if strand2 == "-":
                secondary_exons = secondary_exons[::-1]

            priority_gene_trim = (start_trim1, end_trim1)
            secondary_gene_trim = (start_trim2, end_trim2)
        else:
            priority_match = hmmer_match_details2
            secondary_match = hmmer_match_details1
            priority_exons = np.sort(exons2, order=["low_coord"])
            if strand2 == "-":
                priority_exons = priority_exons[::-1]

            secondary_exons = np.sort(exons1, order=["low_coord"])
            if strand1 == "-":
                secondary_exons = secondary_exons[::-1]

            priority_gene_trim = (start_trim2, end_trim2)
            secondary_gene_trim = (start_trim1, end_trim1)

        priority_env_coords = priority_match[0]["env_coords"]
        secondary_env_coords = secondary_match[0]["env_coords"]

        priority_used_exons, priority_unused_exons = self.find_unused_exons(
            priority_env_coords, priority_exons, priority_gene_trim
        )
        secondary_used_exons, secondary_unused_exons = self.find_unused_exons(
            secondary_env_coords, secondary_exons, secondary_gene_trim
        )

        priority_used_exons = np.array(priority_used_exons, dtype=exon_records.dtype)
        priority_unused_exons = np.array(
            priority_unused_exons, dtype=exon_records.dtype
        )
        secondary_used_exons = np.array(secondary_used_exons, dtype=exon_records.dtype)
        secondary_unused_exons = np.array(
            secondary_unused_exons, dtype=exon_records.dtype
        )

        entries_to_remove = []
        # Check if secondary match uses priority match exons
        used_exons = np.concatenate((priority_used_exons, secondary_used_exons))
        overlaps = self.run_overlap_finder(used_exons)

        if len(overlaps) > 0:
            # Remove secondary match
            entries_to_remove = secondary_exons
            return entries_to_remove

        # Check to see if unused priority exons are used by secondary match
        entries_to_remove = self.get_entries_to_remove(
            secondary_used_exons, priority_unused_exons, entries_to_remove
        )

        # Check to see if unused secondary exons are used by priority match
        entries_to_remove = self.get_entries_to_remove(
            priority_used_exons, secondary_unused_exons, entries_to_remove
        )

        # Check to see if unused secondary exons overlap with priority unused exons
        entries_to_remove = self.get_entries_to_remove(
            priority_unused_exons, secondary_unused_exons, entries_to_remove
        )

        return entries_to_remove

    def get_entries_to_remove(
        self, priority_exons, secondary_exons, entries_to_remove=[]
    ):
        exons_to_check = np.concatenate((priority_exons, secondary_exons))

        overlaps = self.run_overlap_finder(exons_to_check)
        if len(overlaps) > 0:
            for overlap in overlaps:
                if overlap[0] in entries_to_remove or overlap[1] in entries_to_remove:
                    continue

                if overlap[0] in secondary_exons or overlap[1] in secondary_exons:
                    entries_to_remove.extend(secondary_exons)

        return entries_to_remove

    def find_unused_exons(self, env_coords, exons, gene_trim):
        remaining_hmm_region = 0
        unused_exons = []
        used_exons = []
        idx = 0
        hmm_coords = env_coords[idx]
        start_trim, end_trim = gene_trim[0] / 3, gene_trim[1] / 3
        exon_cumul_len = 0
        for entry in exons:
            exon_matched = False
            exon_size_nt = int(entry["high_coord"]) - int(entry["low_coord"]) + 1
            exon_size_aa = round(exon_size_nt / 3)
            exon_cumul_len += exon_size_aa
            if remaining_hmm_region > exon_size_aa:
                remaining_hmm_region -= exon_size_aa
                exon_matched = True

            elif remaining_hmm_region:
                exon_matched = True

            elif hmm_coords:
                while hmm_coords[0] - start_trim < exon_cumul_len + 1:
                    # hmm starts within exon
                    exon_matched = True
                    if hmm_coords[1] - start_trim <= exon_cumul_len + 1:
                        # hmm ends within exon; cycle to the next hmm region
                        try:
                            idx += 1
                            hmm_coords = env_coords[idx]
                        except IndexError:
                            hmm_coords = None
                            break
                        continue
                    else:
                        remaining_hmm_region = (
                            hmm_coords[1] - start_trim - exon_size_aa + 1
                        )
                        break
            if exon_matched:
                used_exons.append(entry)
            else:
                unused_exons.append(entry)
        if len(used_exons) > 0:
            used_exons, unused_exons = self.adjust_exon_categories(
                used_exons, unused_exons
            )
        return used_exons, unused_exons

    @staticmethod
    def adjust_exon_categories(used_exons, unused_exons):
        """
        Ensure there are no unused exons sandwiched between used exons
        :param used_exons:
        :param unused_exons:
        :return:
        """

        used_exons_start = [x["low_coord"] for x in used_exons]
        used_exons_end = [x["high_coord"] for x in used_exons]
        start = min(used_exons_start)
        stop = max(used_exons_end)
        exons_relabeled = set()
        for exon_unused in unused_exons:
            exon_id = "{}|{}-{}".format(
                exon_unused["gene_id"],
                exon_unused["low_coord"],
                exon_unused["high_coord"],
            )
            for exon_relabeled in exons_relabeled:
                if np.array_equal(exon_id, exon_relabeled):
                    continue
            if (
                exon_unused["low_coord"] >= start and exon_unused["low_coord"] < stop
            ) or (
                exon_unused["high_coord"] > start and exon_unused["high_coord"] < stop
            ):
                # find exons that either start or stop within the "used" range
                used_exons.append(exon_unused)
                exons_relabeled.add(exon_id)

        # remove exons that are now used from the unused list
        new_unused_exons = []
        for exon_unused in unused_exons:
            exon_id = "{}|{}-{}".format(
                exon_unused["gene_id"],
                exon_unused["low_coord"],
                exon_unused["high_coord"],
            )
            if exon_id not in exons_relabeled:
                new_unused_exons.append(exon_unused)

        return used_exons, new_unused_exons


class GenomeAnalysisEukaryotesAugustus(BLASTAnalysis, GenomeAnalysisEukaryotes):
    def __init__(self):
        super().__init__()
        self._long = self.config.getboolean("busco_run", "long")
        try:
            self._target_species_initial = self.config.get(
                "busco_run", "augustus_species"
            )
            self._target_species = self._target_species_initial
        except KeyError:
            raise BuscoError(
                "Something went wrong. Eukaryota datasets should specify an augustus species."
            )
        try:
            self._augustus_parameters = self.config.get(
                "busco_run", "augustus_parameters"
            ).replace(",", " ")
        except NoOptionError:
            self._augustus_parameters = ""
        self.mkblast_runner = None
        self.tblastn_runner = None
        self.augustus_runner = None
        self.gff2gb_runner = None
        self.new_species_runner = None
        self.etraining_runner = None
        self.optimize_augustus_runner = None

    def init_tools(self):
        super().init_tools()
        self.augustus_runner = AugustusRunner()
        self.gff2gb_runner = GFF2GBRunner()
        self.new_species_runner = NewSpeciesRunner()
        self.etraining_runner = ETrainingRunner()

        if self._long:
            self.optimize_augustus_runner = OptimizeAugustusRunner()

    def reset(self):
        super().reset()
        self.augustus_runner.reset()
        self.gff2gb_runner.reset()
        self.new_species_runner.reset()
        self.etraining_runner.reset()
        if self._long:
            self.optimize_augustus_runner.reset()

    def cleanup(self):
        try:
            self.augustus_runner.move_retraining_parameters()
            self.config.set(
                "busco_run", "augustus_species", self._target_species_initial
            )  # Reset parameter for batch mode
        except (OSError, AttributeError):
            pass
        super().cleanup()

    def run_analysis(self):
        """This function calls all needed steps for running the analysis."""
        super().run_analysis()
        self._run_mkblast()
        self._run_tblastn()
        self._run_augustus(self.tblastn_runner.coords)
        self.gene_details = self.augustus_runner.gene_details
        self.run_hmmer(self.augustus_runner.output_sequences)
        self._rerun_analysis()

    def run_overlap_finder(self, *args):  # placeholder for now
        return []

    def _rerun_augustus(self, coords):
        missing_and_fragmented_buscos = self.hmmer_runner.missing_buscos + list(
            self.hmmer_runner.fragmented_buscos.keys()
        )
        logger.info(
            "Re-running Augustus with the new metaparameters, number of target BUSCOs: {}".format(
                len(missing_and_fragmented_buscos)
            )
        )
        missing_and_fragmented_coords = {
            busco: coords[busco]
            for busco in coords
            if busco in missing_and_fragmented_buscos
        }
        logger.debug("Trained species folder is {}".format(self._target_species))
        self._run_augustus(missing_and_fragmented_coords)
        return

    @log(
        "Starting second step of analysis. The gene predictor Augustus is retrained using the results from the "
        "initial run to yield more accurate results.",
        logger,
    )
    def _rerun_analysis(self):

        self.augustus_runner.make_gff_files(self.hmmer_runner.single_copy_buscos)
        self._run_tblastn(
            missing_and_frag_only=True, ancestral_variants=self._has_variants_file
        )
        self._run_gff2gb()
        self._run_new_species()
        self.config.set(
            "busco_run", "augustus_species", self.new_species_runner.new_species_name
        )
        self._target_species = self.new_species_runner.new_species_name
        self._run_etraining()

        if self._long:
            self._run_optimize_augustus(self.new_species_runner.new_species_name)
            self._run_etraining()

        try:
            self._rerun_augustus(self.tblastn_runner.coords)
            self.gene_details.update(self.augustus_runner.gene_details)
            self.run_hmmer(self.augustus_runner.output_sequences)
            self.augustus_runner.make_gff_files(self.hmmer_runner.single_copy_buscos)
            self.augustus_runner.make_gff_files(self.hmmer_runner.multi_copy_buscos)
            self.augustus_runner.make_gff_files(self.hmmer_runner.fragmented_buscos)
            self.hmmer_runner.write_buscos_to_file()
        except NoGenesError:
            logger.warning("No genes found on Augustus rerun.")

        return

    @log("Running Augustus gene predictor on BLAST search results.", logger)
    def _run_augustus(self, coords):
        self.augustus_runner.configure_runner(
            self.tblastn_runner.output_seqs,
            coords,
        )

        if self.restart and self.augustus_runner.check_previous_completed_run():
            run = "2nd" if self.augustus_runner.run_number == 2 else "1st"
            logger.info(
                "Skipping {} augustus run as output already processed".format(run)
            )
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.augustus_runner.run()
        self.augustus_runner.process_output()

    def _run_etraining(self):
        """Train on new training set (complete single copy buscos)"""
        self.etraining_runner.configure_runner(self.new_species_runner.new_species_name)
        if self.restart and self.etraining_runner.check_previous_completed_run():
            logger.info("Skipping etraining as it has already been done")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.etraining_runner.run()
        return

    @log("Converting predicted genes to short genbank files", logger)
    def _run_gff2gb(self):
        self.gff2gb_runner.configure_runner(self.hmmer_runner.single_copy_buscos)
        if self.restart and self.gff2gb_runner.check_previous_completed_run():
            logger.info("Skipping gff2gb conversion as it has already been done")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.gff2gb_runner.run()
        return

    @log(
        "All files converted to short genbank files, now training Augustus using Single-Copy Complete BUSCOs",
        logger,
    )
    def _run_new_species(self):
        """Create new species config file from template"""
        self.new_species_runner.configure_runner()
        if self.restart and self.new_species_runner.check_previous_completed_run():
            logger.info("Skipping new species creation as it has already been done")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.new_species_runner.run()
        return

    def _run_optimize_augustus(self, new_species_name):
        """ long mode (--long) option - runs all the Augustus optimization scripts (adds ~1 day of runtime)"""
        logger.warning(
            "Optimizing augustus metaparameters, this may take a very long time, started at {}".format(
                time.strftime("%m/%d/%Y %H:%M:%S")
            )
        )
        self.optimize_augustus_runner.configure_runner(
            self.augustus_runner.output_folder, new_species_name
        )
        self.optimize_augustus_runner.run()
        return


class GenomeAnalysisEukaryotesMetaeuk(GenomeAnalysisEukaryotes):
    def __init__(self):
        super().__init__()
        self.metaeuk_runner = None

    def init_tools(self):
        super().init_tools()

        self.metaeuk_runner = MetaeukRunner()

    def reset(self):
        super().reset()
        self.metaeuk_runner.reset()

    def run_analysis(self):
        """This function calls all needed steps for running the analysis."""
        super().run_analysis()
        incomplete_buscos = None
        for i in range(2):
            try:
                self._run_metaeuk(incomplete_buscos)
                self.gene_details.update(self.metaeuk_runner.gene_details)
                self.run_hmmer(
                    self.metaeuk_runner.pred_protein_seqs_modified,
                    busco_ids=incomplete_buscos,
                )
                incomplete_buscos = self.hmmer_runner.missing_buscos + list(
                    self.hmmer_runner.fragmented_buscos.keys()
                )
                if len(incomplete_buscos) == 0:
                    break
            except NoRerunFile:
                if i == 1:
                    logger.info("Metaeuk rerun did not find any genes")
                else:
                    raise NoGenesError("Metaeuk")

        try:
            self.metaeuk_runner.combine_run_results()
        except FileNotFoundError:
            # This exception should only happen if the rerun file does not exist. If the initial run file was
            # missing there would have been a BatchFatalError call above. The index 0 sets the "combined" file to the
            # output of the initial run.
            self.metaeuk_runner.combined_pred_protein_seqs = (
                self.metaeuk_runner.pred_protein_mod_files[0]
            )
            self.metaeuk_runner.combined_nucleotides_seqs = (
                self.metaeuk_runner.codon_mod_files[0]
            )
        self.hmmer_runner.write_buscos_to_file()
        self.write_gff_files()

    def run_overlap_finder(self, exon_records_arr):
        overlaps = self.metaeuk_runner.find_overlaps(exon_records_arr)
        return overlaps

    def write_gff_files(self):
        busco_seqs_folders = [
            self.hmmer_runner.single_copy_sequences_folder,
            self.hmmer_runner.multi_copy_sequences_folder,
            self.hmmer_runner.fragmented_sequences_folder,
        ]
        for folder in busco_seqs_folders:
            existing_gff_files = glob.glob(os.path.join(folder, "*.gff"))
            for f in existing_gff_files:
                os.remove(f)

        for gff_file in self.metaeuk_runner.gff_files[
            ::-1
        ]:  # Write GFF results to file, starting with the rerun results and then using the initial results
            self.metaeuk_runner.gff_file = gff_file
            self.metaeuk_runner.write_gff_files(*busco_seqs_folders)

    def _run_metaeuk(self, incomplete_buscos):
        self.metaeuk_runner.configure_runner(incomplete_buscos)
        if self.restart and self.metaeuk_runner.check_previous_completed_run():
            logger.info("Skipping Metaeuk run as already run")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.metaeuk_runner.run()

        structured_arr_contents = self.metaeuk_runner.parse_output()
        self.metaeuk_runner.create_filtered_sequence_files(structured_arr_contents)
        self.metaeuk_runner.combine_gene_details()

    def cleanup(self):
        try:
            self.metaeuk_runner.remove_tmp_files()
        except OSError:
            pass
        super().cleanup()


class GenomeAnalysisEukaryotesMiniprot(GenomeAnalysisEukaryotes):
    def __init__(self):
        super().__init__()
        self.miniprot_index_runner = None
        self.miniprot_align_runner = None
        self.gene_details = {}
        self.gene_update_mapping = defaultdict(dict)
        self.cpus = int(self.config.get("busco_run", "cpu"))
        self.filtered_records = defaultdict(list)
        self.filtered_busco_hits = []

    def init_tools(self):
        super().init_tools()
        self.miniprot_index_runner = MiniprotIndexRunner()
        self.miniprot_align_runner = MiniprotAlignRunner()

    def reset(self):
        super().reset()
        self.miniprot_index_runner.reset()
        self.miniprot_align_runner.reset()

    def run_analysis(self):
        """This function calls all needed steps for running the analysis."""
        super().run_analysis()
        incomplete_buscos = None
        try:
            self.run_miniprot(incomplete_buscos)
            self.hmmer_runner.miniprot_pipeline = True
            self.gene_details = self.miniprot_align_runner.gene_details
            self.run_hmmer(
                self.miniprot_align_runner.output_sequences, busco_ids=incomplete_buscos
            )
        except NoRerunFile:
            raise NoGenesError("Miniprot")

        self.hmmer_runner.write_buscos_to_file()

    def run_miniprot(self, incomplete_buscos):
        self.miniprot_index_runner.configure_runner()
        if self.restart and self.miniprot_index_runner.check_previous_completed_run():
            logger.info("Skipping Miniprot indexing run as already run")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.miniprot_index_runner.run()

        self.miniprot_align_runner.configure_runner(incomplete_buscos)
        if self.restart and self.miniprot_align_runner.check_previous_completed_run():
            logger.info("Skipping Miniprot aligning run as already run")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.miniprot_align_runner.run()
        self.miniprot_align_runner.parse_output()
        self.miniprot_align_runner.filter()
        self.miniprot_align_runner.record_gene_details()
        self.miniprot_align_runner.write_protein_sequences_per_busco()

    def run_overlap_finder(self, exon_records_arr):
        overlaps = self.miniprot_align_runner.find_overlaps(exon_records_arr)
        return overlaps

    def record_results(self):
        self.hmmer_runner.record_results()

    def cleanup(self):
        super().cleanup()
