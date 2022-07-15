#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: GeneSetAnalysis
   :synopsis: GeneSetAnalysis implements genome analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 5.4.0

Copyright (c) 2016-2022, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
from busco.analysis.BuscoAnalysis import BuscoAnalysis
from busco.BuscoLogger import BuscoLogger
from busco.analysis.Analysis import ProteinAnalysis
from Bio import SeqIO

logger = BuscoLogger.get_logger(__name__)


class GeneSetAnalysis(ProteinAnalysis, BuscoAnalysis):
    """
    This class runs a BUSCO analysis on a gene set.
    """

    _mode = "proteins"

    def __init__(self):
        """
        Initialize an instance.
        """
        super().__init__()
        self.sequences_aa = {
            record.id: record for record in list(SeqIO.parse(self.input_file, "fasta"))
        }

    def cleanup(self):
        super().cleanup()

    def run_analysis(self):
        """
        This function calls all needed steps for running the analysis.
        """
        super().run_analysis()
        self.run_hmmer(self.input_file)
        self.hmmer_runner.write_buscos_to_file(self.sequences_aa)
        # if self._tarzip:
        #     self._run_tarzip_hmmer_output()
        return

    def reset(self):
        super().reset()
