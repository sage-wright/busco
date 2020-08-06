#!/usr/bin/env python
# coding: utf-8
"""
.. module:: GenomeAnalysis
   :synopsis: GenomeAnalysis implements genome analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 3.0.0

Copyright (c) 2016-2020, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
from busco.BuscoAnalysis import BuscoAnalysis
from busco.Analysis import NucleotideAnalysis
from busco.BuscoTools import ProdigalRunner, MetaeukRunner
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from abc import ABCMeta, abstractmethod

logger = BuscoLogger.get_logger(__name__)


class GenomeAnalysis(NucleotideAnalysis, BuscoAnalysis, metaclass=ABCMeta):

    _mode = "genome"

    def __init__(self):
        super().__init__()

    @abstractmethod
    def run_analysis(self):
        super().run_analysis()

    def init_tools(self):
        """
        Initialize tools needed for Genome Analysis.
        :return:
        """
        super().init_tools()

    # def _run_tarzip_augustus_output(self): # Todo: rewrite using tarfile
    #     """
    #     This function tarzips results folder
    #     """
    #     # augustus_output/predicted_genes
    #
    #     self._p_open(["tar", "-C", "%saugustus_output" % self.main_out,
    #                   "-zcf", "%saugustus_output/predicted_genes.tar.gz" %
    #                   self.main_out, "predicted_genes", "--remove-files"],
    #                  "bash", shell=False)
    #     # augustus_output/extracted_proteins
    #     self._p_open(["tar", "-C", "%saugustus_output" % self.main_out,
    #                   "-zcf", "%saugustus_output/extracted_proteins.tar.gz" %
    #                   self.main_out, "extracted_proteins", "--remove-files"],
    #                  "bash", shell=False)
    #     # augustus_output/gb
    #     self._p_open(["tar", "-C", "%saugustus_output" % self.main_out,
    #                   "-zcf", "%saugustus_output/gb.tar.gz" % self.main_out, "gb", "--remove-files"],
    #                  "bash", shell=False)
    #     # augustus_output/gffs
    #     self._p_open(["tar", "-C", "%saugustus_output" % self.main_out,
    #                   "-zcf", "%saugustus_output/gffs.tar.gz" %
    #                   self.main_out, "gffs", "--remove-files"], "bash", shell=False)
    #     # single_copy_busco_sequences
    #     self._p_open(["tar", "-C", "%s" % self.main_out, "-zcf",
    #                   "%ssingle_copy_busco_sequences.tar.gz" % self.main_out,
    #                   "single_copy_busco_sequences", "--remove-files"], "bash", shell=False)

    # def set_rerun_busco_command(self, clargs):
    #     """
    #     This function sets the command line to call to reproduce this run
    #     """
    #     clargs.extend(["-sp", self._target_species])
    #     super().set_rerun_busco_command(clargs)


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
        self.run_hmmer(self.prodigal_runner.output_faa)
        self.hmmer_runner.write_buscos_to_file(self.sequences_aa, self.sequences_nt)
        return

    def init_tools(self):
        """
        Init the tools needed for the analysis
        """
        super().init_tools()
        self.prodigal_runner = ProdigalRunner()

    @log("***** Run Prodigal on input to predict and extract genes *****", logger)
    def _run_prodigal(self):
        """
        Run Prodigal on input file to detect genes.
        :return:
        """
        if self.restart and self.prodigal_runner.check_previous_completed_run():
            logger.info("Skipping Prodigal run as it has already completed")
            self.prodigal_runner.get_gene_details()
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.prodigal_runner.run()
        self.gene_details = self.prodigal_runner.gene_details
        self.sequences_nt = self.prodigal_runner.sequences_nt
        self.sequences_aa = self.prodigal_runner.sequences_aa

        return


class GenomeAnalysisEukaryotes(GenomeAnalysis):
    """
    This class runs a BUSCO analysis on a eukaryote genome.
    """
    def __init__(self):
        super().__init__()

        self.sequences_nt = {}
        self.sequences_aa = {}

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

        self.metaeuk_runner = MetaeukRunner()

        return

    def run_analysis(self):
        """This function calls all needed steps for running the analysis."""
        super().run_analysis()
        self._run_metaeuk()
        self.gene_details = self.metaeuk_runner.gene_details
        self.sequences_aa = self.metaeuk_runner.sequences_aa
        self.run_hmmer(self.metaeuk_runner.pred_protein_seqs_modified)
        self.hmmer_runner.write_buscos_to_file(self.sequences_aa)

    def _run_metaeuk(self):
        if self.restart and self.metaeuk_runner.check_previous_completed_run():
            logger.info("Skipping Metaeuk run as already run")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.metaeuk_runner.run()
        self.metaeuk_runner.get_gene_details()
        self.metaeuk_runner.edit_file_header()

    # def set_rerun_busco_command(self, clargs):
    #     """
    #     This function sets the command line to call to reproduce this run
    #     """
    #     clargs.extend(["-sp", self._target_species])
    #     if self._augustus_parameters:
    #         clargs.extend(["--augustus_parameters", "\"%s\"" % self._augustus_parameters])
    #     super().set_rerun_busco_command(clargs)
