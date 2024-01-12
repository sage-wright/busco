# coding: utf-8
"""
GeneSetAnalysis.py

Module called for proteins mode.

Author(s): Matthew Berkeley, Mathieu Seppey, Mose Manni, Felipe Simao, Rob Waterhouse

Copyright (c) 2015-2024, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

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
        self.gene_details = {
            record.id: {"aa_seq": record}
            for record in list(SeqIO.parse(self.input_file, "fasta"))
        }

    def cleanup(self):
        super().cleanup()

    def run_analysis(self):
        """
        This function calls all needed steps for running the analysis.
        """
        super().run_analysis()
        self.run_hmmer(self.input_file)
        self.hmmer_runner.write_buscos_to_file()
        return

    def reset(self):
        super().reset()
