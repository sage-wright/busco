#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: ViralAnalysis
   :synopsis: ViralAnalysis implements genome analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 3.0.0

Copyright (c) 2016-2017, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
from busco.BuscoAnalysis import BuscoAnalysis
from busco.BuscoLogger import BuscoLogger

logger = BuscoLogger.get_logger(__name__)

class ViralAnalysis(BuscoAnalysis):
    """
    This class runs a BUSCO analysis on a gene set.
    """

    _mode = "proteins"


    def __init__(self, params):
        """
        Initialize an instance.
        :param params: Values of all parameters that have to be defined
        :type params: PipeConfig
        """
        super().__init__(params)
        # data integrity checks not done by the parent class
        if self.check_protein_file():
            ViralAnalysis._logger.error("Please provide a genome file as input or run BUSCO in protein mode (--mode proteins)")
            raise SystemExit

    def run_analysis(self):
        """
        This function calls all needed steps for running the analysis.
        """
        super().run_analysis()
        self._translate_virus()
        self._sequences = self.translated_proteins
        self._run_hmmer()
        # if self._tarzip:
        #     self._run_tarzip_hmmer_output()

    def _init_tools(self):
        """
        Init the tools needed for the analysis
        """
        super()._init_tools()

    def _run_hmmer(self):
        """
        This function runs hmmsearch.
        """
        super()._run_hmmer()

    def _sixpack(self, seq):
        """
        Gets the sixframe translation for the provided sequence
        :param seq: the sequence to be translated
        :type seq: str
        :return: the six translated sequences
        :rtype: list
        """
        s1 = seq
        s2 = seq[1:]
        s3 = seq[2:]
        rev = ""
        for letter in seq[::-1]:
            try:
                rev += BuscoAnalysis.COMP[letter]
            except KeyError:
                rev += BuscoAnalysis.COMP["N"]
        r1 = rev
        r2 = rev[1:]
        r3 = rev[2:]
        transc = []
        frames = [s1, s2, s3, r1, r3, r2]
        for sequence in frames:
            part = ""
            new = ""
            for letter in sequence:
                if len(part) == 3:
                    try:
                        new += BuscoAnalysis.CODONS[part]
                    except KeyError:
                        new += "X"
                    part = ""
                    part += letter
                else:
                    part += letter
            if len(part) == 3:
                try:
                    new += BuscoAnalysis.CODONS[part]
                except KeyError:
                    new += "X"
            transc.append(new)
        return transc
    
    def _translate_virus(self):
        """
        Prepares viral genomes for a BUSCO
        protein analysis.
        1) Translate any sequences in 6 frames
        2) Split the sequences on the stops
        3) Remove sequences shorter than 50aa
        :return: file name 
        :rtype: string
        """
        with open(self._sequences, "r") as f1:
            with open(self.mainout + "translated_proteins.faa", "w") as o1:
                seqs = {}
                # parse file, retrieve all sequences
                for line in f1:
                    if line.startswith(">"):
                        header = line.strip()
                        seqs[header] = ""
                    else:
                        seqs[header] += line.strip()

                # feed sequences to 6 frame translator
                # then split on STOP codons
                ctg_ct = 1
                for seqid in seqs:
                    seq_6f = self._sixpack(seqs[seqid])
                    nb_frame = 1
                    for frame in seq_6f:
                        valid_ts_ct = 1
                        # chop at the stop codons
                        chopped_seqs = frame.split("X")
                        for short_seq in chopped_seqs:
                            # must have at least 50 A.A. to be considered further
                            if len(short_seq) >= 50:
                                # ctg nb, frame nb, transcript nb
                                o1.write(">seq_n%s_f%s_t%s\n" % (ctg_ct, nb_frame, valid_ts_ct))
                                o1.write("%s\n" % short_seq)
                                valid_ts_ct += 1
                            else:
                                pass
                        nb_frame += 1
                    ctg_ct += 1
        self.translated_proteins = self.main_out + "translated_proteins.faa"
