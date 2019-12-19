#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: GeneSetAnalysis
   :synopsis: GeneSetAnalysis implements genome analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 3.0.0

Copyright (c) 2016-2017, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
from busco.BuscoAnalysis import BuscoAnalysis
from busco.BuscoLogger import BuscoLogger
from busco.Analysis import ProteinAnalysis
from Bio import SeqIO

logger = BuscoLogger.get_logger(__name__)


class GeneSetAnalysis(ProteinAnalysis, BuscoAnalysis):
    """
    This class runs a BUSCO analysis on a gene set.
    """
    _mode = 'proteins'

    def __init__(self, config):
        """
        Initialize an instance.
        :param params: Values of all parameters that have to be defined
        :type params: PipeConfig
        """
        super().__init__(config)
        self.sequences_aa = {record.id: record for record in list(SeqIO.parse(self._input_file, "fasta"))}

    def _cleanup(self):
        super()._cleanup()

    def run_analysis(self):
        """
        This function calls all needed steps for running the analysis.
        """
        super().run_analysis()
        self.run_hmmer(self._input_file)
        self._write_buscos_to_file(self.sequences_aa)
        self._cleanup()
        # if self._tarzip:
        #     self._run_tarzip_hmmer_output()
        return

    def create_dirs(self):
        super().create_dirs()
