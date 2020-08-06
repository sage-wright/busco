#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: TranscriptomeAnalysis
   :synopsis:TranscriptomeAnalysis implements genome analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 5.0.beta

Copyright (c) 2016-2020, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
import os
from busco.BuscoAnalysis import BuscoAnalysis
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from Bio.Seq import reverse_complement, translate
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from busco.Analysis import NucleotideAnalysis
from busco.BuscoTools import TBLASTNRunner, MKBLASTRunner


logger = BuscoLogger.get_logger(__name__)

# todo: catch multiple buscos on one transcript


class TranscriptomeAnalysis(NucleotideAnalysis, BuscoAnalysis):
    """
    Analysis on a transcriptome.
    """

    _mode = "transcriptome"

    def __init__(self):
        """
        Initialize an instance.
        """
        super().__init__()

    def run_analysis(self):
        """
        This function calls all needed steps for running the analysis.
        """

        super().run_analysis()

        self._run_mkblast()
        self._run_tblastn(ancestral_variants=self._has_variants_file)

        protein_seq_files = self._translate_seqs(self.tblastn_runner.coords)

        self.run_hmmer(protein_seq_files)
        # Note BUSCO matches are not written to file, as we have not yet developed a suitable protocol for
        # Transcriptomes
        # if self._tarzip:
        #     self._run_tarzip_hmmer_output()
        #     self._run_tarzip_translated_proteins()
        return

    def init_tools(self):
        super().init_tools()
        self.mkblast_runner = MKBLASTRunner()
        self.tblastn_runner = TBLASTNRunner()

        if self.mkblast_runner.version != self.tblastn_runner.version:
            logger.warning("You are using version {} of makeblastdb and version {} of tblastn.".format(
                self.mkblast_runner.version, self.tblastn_runner.version))

    def cleanup(self):
        """
        This function cleans temporary files.
        """
        super().cleanup()

    def _run_mkblast(self):
        if self.restart and self.mkblast_runner.check_previous_completed_run():
            logger.info("Skipping makeblastdb as BLAST DB already exists at {}".format(self.mkblast_runner.output_db))
        else:
            self.restart = False  # Turn off restart mode if this is the entry point
            self.config.set("busco_run", "restart", str(self.restart))
            self.mkblast_runner.run()
        if len(os.listdir(os.path.split(self.mkblast_runner.output_db)[0])) == 0:
            raise SystemExit("makeblastdb failed to create a BLAST DB at {}".format(self.mkblast_runner.output_db))

    def _run_tblastn(self, missing_and_frag_only=False, ancestral_variants=False):

        incomplete_buscos = (self.hmmer_runner.missing_buscos + list(self.hmmer_runner.fragmented_buscos.keys())
                             if missing_and_frag_only else None)  # This parameter is only used on the re-run

        self.tblastn_runner.configure_runner(self.mkblast_runner.output_db, missing_and_frag_only,
                                             ancestral_variants, incomplete_buscos)
        if self.restart and self.tblastn_runner.check_previous_completed_run():
            logger.info("Skipping tblastn as results already exist at {}".format(self.tblastn_runner.blast_filename))
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.tblastn_runner.run()
        self.tblastn_runner.get_coordinates()
        self.tblastn_runner.filter_best_matches()
        self.tblastn_runner.write_coordinates_to_file()  # writes to "coordinates.tsv"
        self.tblastn_runner.write_contigs()
        return

    @staticmethod
    def six_frame_translation(seq):
        """
        Gets the sixframe translation for the provided sequence
        :param seq: the sequence to be translated
        :type seq: str
        :return: the six translated sequences
        :rtype: list
        """
        descriptions = {1: "orig_seq_frame_1",
                        2: "orig_seq_frame_2",
                        3: "orig_seq_frame_3",
                        -1: "rev_comp_frame_1",
                        -2: "rev_comp_frame_2",
                        -3: "rev_comp_frame_3"}

        # Based on code excerpt from https://biopython.org/DIST/docs/api/Bio.SeqUtils-pysrc.html#six_frame_translations
        anti = reverse_complement(seq)
        translated_seqs = {}
        for i in range(3):
            fragment_length = 3 * ((len(seq) - i) // 3)
            translated_seqs[descriptions[i+1]] = (translate(seq[i:i + fragment_length], stop_symbol="X"))
            translated_seqs[descriptions[-(i+1)]] = (translate(anti[i:i + fragment_length], stop_symbol="X"))
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
            output_filename = os.path.join(translated_proteins_dir, "{}.faa".format(busco_id))
            protein_seq_files.append(output_filename)
            translated_records = []
            for contig_name in contig_info:
                tmp_filename = os.path.join(self.tblastn_runner.output_seqs, "{}.temp".format(
                    contig_name[:100]))  # Avoid very long filenames
                for record in SeqIO.parse(tmp_filename, "fasta"):  # These files will only ever have one sequence,
                    # but BioPython examples always parse them in an iterator.
                    translated_seqs = self.six_frame_translation(record.seq)
                    for desc_id in translated_seqs:  # There are six possible translated sequences
                        prot_seq = translated_seqs[desc_id]
                        translated_records.append(SeqRecord(prot_seq, id=record.id, description=desc_id))

            with open(output_filename, "w") as out_faa:
                SeqIO.write(translated_records, out_faa, "fasta")

        return protein_seq_files

    # def _run_tarzip_translated_proteins(self):
    #     """
    #     This function tarzips results folder
    #     """
    #     # translated_proteins # Todo: rewrite with tarfile module
    #     self._p_open(["tar", "-C", "%s" % self.mainout, "-zcf",
    #                  "%stranslated_proteins.tar.gz" % self.mainout, "translated_proteins", "--remove-files"], "bash",
    #                  shell=False)
