#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: TranscriptomeAnalysis
   :synopsis:TranscriptomeAnalysis implements genome analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 3.0.0

Copyright (c) 2016-2017, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
import os
import time

from busco.BuscoAnalysis import BuscoAnalysis
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from Bio.Seq import reverse_complement, translate
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from busco.Analysis import NucleotideAnalysis
from busco.Toolset import Tool


logger = BuscoLogger.get_logger(__name__)

# todo: catch multiple buscos on one transcript

class TranscriptomeAnalysis(NucleotideAnalysis, BuscoAnalysis):
    """
    Analysis on a transcriptome.
    """

    _mode = "transcriptome"

    def __init__(self, config):
        """
        Initialize an instance.
        :param config: Values of all parameters that have to be defined
        :type config: BuscoConfig
        """
        super().__init__(config)


    def run_analysis(self):
        """
        This function calls all needed steps for running the analysis.
        """

        super().run_analysis()

        # if self._restart:  # todo: reimplement restart mode
        #     checkpoint = self.get_checkpoint(reset_random_suffix=True)
        #     logger.warning("Restarting an uncompleted run")
        # else:
        #     checkpoint = 0  # all steps will be done
        # if checkpoint < 1:

        self._run_mkblast()
        coords = self._run_tblastn(ancestral_variants=self._has_variants_file)

        protein_seq_files = self._translate_seqs(coords)

        self.run_hmmer(protein_seq_files)
        # Note BUSCO matches are not written to file, as we have not yet developed a suitable protocol for Transcriptomes
        self._cleanup()
        # if self._tarzip:
        #     self._run_tarzip_hmmer_output()
        #     self._run_tarzip_translated_proteins()
        return

    def create_dirs(self): # todo: remove this as abstract method, review all abstract methods
        super().create_dirs()

    def init_tools(self): # todo: This should be an abstract method

        super().init_tools()
        try:
            assert(isinstance(self._mkblast_tool, Tool))
        except AttributeError:
            self._mkblast_tool = Tool("makeblastdb", self._config)
        except AssertionError:
            raise SystemExit("mkblast should be a tool")

        try:
            assert(isinstance(self._tblastn_tool, Tool))
        except AttributeError:
            self._tblastn_tool = Tool("tblastn", self._config)
        except AssertionError:
            raise SystemExit("tblastn should be a tool")

    def check_tool_dependencies(self):
        blast_version = self._get_blast_version()
        if blast_version not in ["2.2", "2.3"]:  # Known problems with multithreading on BLAST 2.4-2.9.
            if blast_version == "2.9" and self._tblastn_tool.cmd.endswith(
                    "tblastn_June13"):  # NCBI sent a binary with this name that avoids the multithreading problems.
                pass
            else:
                logger.warning("You are using BLAST version {}. This is known to yield inconsistent results when "
                               "multithreading. BLAST will run on a single core as a result. For performance improvement, "
                               "please revert to BLAST 2.2 or 2.3.".format(blast_version))
                self.blast_cpus = 1

    def _cleanup(self):
        """
        This function cleans temporary files.
        """
        super()._cleanup()

    def six_frame_translation(self, seq):
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

    def _reformats_seq_id(self, seq_id):
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
                tmp_filename = os.path.join(self.tblastn_runner.output_seqs, "{}.temp".format(contig_name[:100]))  # Avoid very long filenames
                for record in SeqIO.parse(tmp_filename, "fasta"):  # These files will only ever have one sequence, but BioPython examples always parse them in an iterator.
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
