#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: TranscriptomeAnalysis
   :synopsis:TranscriptomeAnalysis implements genome analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 5.4.0

Copyright (c) 2016-2022, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
import os
from busco.analysis.BuscoAnalysis import BuscoAnalysis
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from Bio.Seq import reverse_complement, translate
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from busco.analysis.Analysis import NucleotideAnalysis, BLASTAnalysis
from busco.analysis.GenomeAnalysis import GenomeAnalysisEukaryotesMetaeuk
from abc import ABCMeta
from collections import defaultdict


logger = BuscoLogger.get_logger(__name__)


class TranscriptomeAnalysis(NucleotideAnalysis, metaclass=ABCMeta):
    """
    Analysis on a transcriptome.
    """

    _mode = "transcriptome"

    def __init__(self):
        """
        Initialize an instance.
        """
        super().__init__()


class TranscriptomeAnalysisProkaryotes(
    TranscriptomeAnalysis, BLASTAnalysis, BuscoAnalysis
):
    """
    Analysis on a transcriptome.
    """

    def __init__(self):
        """
        Initialize an instance.
        """
        super().__init__()
        self.all_sequences = defaultdict(lambda: defaultdict(dict))

        self.sequences_aa = {}
        self.sequences_nt = {}
        self.complete_seqs_nt = {}
        self.complete_seqs_aa = {}

        self.single_copy_proteins_file = os.path.join(
            self.run_folder, "single_copy_proteins.faa"
        )

    def run_analysis(self):
        """
        This function calls all needed steps for running the analysis.
        """

        super().run_analysis()

        self._run_mkblast()
        self._run_tblastn(ancestral_variants=self._has_variants_file)

        protein_seq_files = self._translate_seqs(self.tblastn_runner.coords)

        self.run_hmmer(protein_seq_files)

        self.prepare_sequences()
        self.write_complete_seqs()

        self.hmmer_runner.write_buscos_to_file(self.sequences_aa, self.sequences_nt)

        # if self._tarzip:
        #     self._run_tarzip_hmmer_output()
        #     self._run_tarzip_translated_proteins()
        return

    def init_tools(self):
        super().init_tools()

    def reset(self):
        super().reset()

    def cleanup(self):
        """
        This function cleans temporary files.
        """
        super().cleanup()

    def match_hmmer_results(self, hmmer_results):
        seqs_aa = {}
        seqs_nt = {}

        for busco_id, gene_matches in hmmer_results.items():
            for gene_id in gene_matches:
                gene_info = gene_matches[gene_id][0]
                frame = gene_info["frame"]
                seq_aa = self.all_sequences[busco_id][gene_id]["translations"][frame]
                seq_nt = self.all_sequences[busco_id][gene_id]["original"]
                seqs_aa[gene_id] = seq_aa
                seqs_nt[gene_id] = seq_nt
        return seqs_aa, seqs_nt

    def prepare_sequences(self):
        hmmer_results_complete = self.hmmer_runner.single_copy_buscos

        sc_seqs_aa, sc_seqs_nt = self.match_hmmer_results(hmmer_results_complete)
        self.complete_seqs_aa = [val for key, val in sc_seqs_aa.items()]

        hmmer_results_remainder = {
            **self.hmmer_runner.multi_copy_buscos,
            **self.hmmer_runner.fragmented_buscos,
        }

        remainder_seqs_aa, remainder_seqs_nt = self.match_hmmer_results(
            hmmer_results_remainder
        )

        self.sequences_aa = {**sc_seqs_aa, **remainder_seqs_aa}
        self.sequences_nt = {**sc_seqs_nt, **remainder_seqs_nt}

        return

    def write_complete_seqs(self):

        with open(self.single_copy_proteins_file, "w") as out_faa:
            SeqIO.write(self.complete_seqs_aa, out_faa, "fasta")

        return

    @staticmethod
    def six_frame_translation(record):
        """
        Gets the sixframe translation for the provided sequence
        :param record: the sequence record to be translated
        :type record: SeqRecord
        :return: the six translated sequences
        :rtype: list
        """
        descriptions = {
            1: "orig_seq_frame_1",
            2: "orig_seq_frame_2",
            3: "orig_seq_frame_3",
            -1: "rev_comp_frame_1",
            -2: "rev_comp_frame_2",
            -3: "rev_comp_frame_3",
        }

        # Based on code excerpt from https://biopython.org/DIST/docs/api/Bio.SeqUtils-pysrc.html#six_frame_translations
        anti = reverse_complement(record.seq)
        translated_seqs = {}
        for i in range(3):
            fragment_length = 3 * ((len(record.seq) - i) // 3)
            translated_seqs[descriptions[i + 1]] = SeqRecord(
                translate(record.seq[i : i + fragment_length], stop_symbol="X"),
                id=record.id,
                name=record.id,
                description=descriptions[i + 1],
            )
            translated_seqs[descriptions[-(i + 1)]] = SeqRecord(
                translate(anti[i : i + fragment_length], stop_symbol="X"),
                id=record.id,
                name=record.id,
                description=descriptions[i + 1],
            )
        return translated_seqs

    @staticmethod
    def _reformats_seq_id(seq_id):
        """
        This function reformats the sequence id to its original values
        :param seq_id: the seq id to reformats
        :type seq_id: str
        :return: the reformatted seq_id
        :rtype: str
        """
        return "_".join(seq_id.split("_")[:-1])

    @log("Translating candidate transcripts", logger)
    def _translate_seqs(self, coords):

        translated_proteins_dir = os.path.join(self.main_out, "translated_proteins")
        if not os.path.exists(translated_proteins_dir):
            os.makedirs(translated_proteins_dir)

        contig_names = []
        for contig_info in coords.values():
            for contig in contig_info:
                contig_names.append(contig)

        protein_seq_files = []
        for busco_id, contig_info in coords.items():
            output_filename = os.path.join(
                translated_proteins_dir, "{}.faa".format(busco_id)
            )
            protein_seq_files.append(output_filename)
            translated_records = []
            for contig_name in contig_info:
                tmp_filename = os.path.join(
                    self.tblastn_runner.output_seqs, "{}.temp".format(contig_name[:100])
                )  # Avoid very long filenames
                for record in SeqIO.parse(
                    tmp_filename, "fasta"
                ):  # These files will only ever have one sequence,
                    # but BioPython examples always parse them in an iterator.
                    translated_seqs = self.six_frame_translation(record)
                    gene_id = "{}:{}-{}".format(
                        record.id,
                        contig_info[contig_name]["contig_start"],
                        contig_info[contig_name]["contig_end"],
                    )
                    self.all_sequences[busco_id][gene_id][
                        "translations"
                    ] = translated_seqs
                    self.all_sequences[busco_id][gene_id]["original"] = record
                    for (
                        desc_id
                    ) in translated_seqs:  # There are six possible translated sequences
                        prot_seq = translated_seqs[desc_id]
                        prot_seq.id = gene_id
                        translated_records.append(prot_seq)

            with open(output_filename, "w") as out_faa:
                SeqIO.write(translated_records, out_faa, "fasta")

        return protein_seq_files


class TranscriptomeAnalysisEukaryotes(
    TranscriptomeAnalysis, GenomeAnalysisEukaryotesMetaeuk
):
    def __init__(self):
        super().__init__()
