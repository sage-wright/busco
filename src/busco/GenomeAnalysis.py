#!/usr/bin/env python
# coding: utf-8
"""
.. module:: GenomeAnalysis
   :synopsis: GenomeAnalysis implements genome analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 5.beta.1

Copyright (c) 2016-2020, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
from busco.BuscoAnalysis import BuscoAnalysis
from busco.Analysis import NucleotideAnalysis, BLASTAnalysis
from busco.BuscoTools import ProdigalRunner, MetaeukRunner, NoRerunFile, AugustusRunner, GFF2GBRunner, \
    NewSpeciesRunner, ETrainingRunner, OptimizeAugustusRunner, NoGenesError
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from abc import ABCMeta, abstractmethod
from configparser import NoOptionError
import time

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

        return

    @abstractmethod
    def run_analysis(self):
        super().run_analysis()



    # def set_rerun_busco_command(self, clargs):
    #     """
    #     This function sets the command line to call to reproduce this run
    #     """
    #     clargs.extend(["-sp", self._target_species])
    #     if self._augustus_parameters:
    #         clargs.extend(["--augustus_parameters", "\"%s\"" % self._augustus_parameters])
    #     super().set_rerun_busco_command(clargs)


class GenomeAnalysisEukaryotesAugustus(BLASTAnalysis, GenomeAnalysisEukaryotes):

    def __init__(self):
        super().__init__()
        self._long = self.config.getboolean("busco_run", "long")
        try:
            self._target_species = self.config.get("busco_run", "augustus_species")
        except KeyError:
            raise SystemExit("Something went wrong. Eukaryota datasets should specify an augustus species.")
        try:
            self._augustus_parameters = self.config.get("busco_run", "augustus_parameters").replace(',', ' ')
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

    def cleanup(self):
        try:
            if self._target_species.startswith("BUSCO"):
                self.augustus_runner.move_retraining_parameters()
        except OSError:
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

    def _rerun_augustus(self, coords):
        missing_and_fragmented_buscos = self.hmmer_runner.missing_buscos + list(
            self.hmmer_runner.fragmented_buscos.keys())
        logger.info("Re-running Augustus with the new metaparameters, number of target BUSCOs: {}".format(
            len(missing_and_fragmented_buscos)))
        missing_and_fragmented_coords = {busco: coords[busco] for busco in coords if busco in
                                         missing_and_fragmented_buscos}
        logger.debug('Trained species folder is {}'.format(self._target_species))
        self._run_augustus(missing_and_fragmented_coords, rerun=True)
        return

    @log("Starting second step of analysis. The gene predictor Augustus is retrained using the results from the "
         "initial run to yield more accurate results.", logger)
    def _rerun_analysis(self):

        self.augustus_runner.make_gff_files(self.hmmer_runner.single_copy_buscos)
        self._run_tblastn(missing_and_frag_only=True, ancestral_variants=self._has_variants_file)
        self._run_gff2gb()
        self._run_new_species()
        self.config.set("busco_run", "augustus_species", self.new_species_runner.new_species_name)
        self._target_species = self.new_species_runner.new_species_name
        self._run_etraining()

        if self._long:
            self._run_optimize_augustus(self.new_species_runner.new_species_name)
            self._run_etraining()

        try:
            self._rerun_augustus(self.tblastn_runner.coords)
            self.gene_details.update(self.augustus_runner.gene_details)
            self.run_hmmer(self.augustus_runner.output_sequences)
            self.hmmer_runner.write_buscos_to_file(self.sequences_aa, self.sequences_nt)
        except NoGenesError:
            logger.warning("No genes found on Augustus rerun.")

        # if self._tarzip:  # todo: zip folders with a lot of output
        #     self._run_tarzip_augustus_output()
        #     self._run_tarzip_hmmer_output()
        # remove the checkpoint, run is done
        # self._set_checkpoint()
        return

    @log("Running Augustus gene predictor on BLAST search results.", logger)
    def _run_augustus(self, coords, rerun=False):
        self.augustus_runner.configure_runner(self.tblastn_runner.output_seqs, coords, self.sequences_aa,
                                              self.sequences_nt, rerun)

        if self.restart and self.augustus_runner.check_previous_completed_run():
            run = "2nd" if rerun else "1st"
            logger.info("Skipping {} augustus run as output already processed".format(run))
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.augustus_runner.run()
        self.augustus_runner.process_output()
        self.sequences_nt = self.augustus_runner.sequences_nt
        self.sequences_aa = self.augustus_runner.sequences_aa

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

    @log("All files converted to short genbank files, now training Augustus using Single-Copy Complete BUSCOs", logger)
    def _run_new_species(self):
        """Create new species config file from template"""
        if self.restart and self.new_species_runner.check_previous_completed_run():
            logger.info("Skipping new species creation as it has already been done")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.new_species_runner.run()
        return

    def _run_optimize_augustus(self, new_species_name):
        """ long mode (--long) option - runs all the Augustus optimization scripts (adds ~1 day of runtime)"""
        logger.warning("Optimizing augustus metaparameters, this may take a very long time, started at {}".format(
            time.strftime("%m/%d/%Y %H:%M:%S")))
        self.optimize_augustus_runner.configure_runner(self.augustus_runner.output_folder, new_species_name)
        self.optimize_augustus_runner.run()
        return


class GenomeAnalysisEukaryotesMetaeuk(GenomeAnalysisEukaryotes):

    def __init__(self):
        super().__init__()


    def init_tools(self):
        super().init_tools()

        self.metaeuk_runner = MetaeukRunner()

    def run_analysis(self):
        """This function calls all needed steps for running the analysis."""
        super().run_analysis()
        incomplete_buscos = None
        for i in range(2):
            try:
                self._run_metaeuk(incomplete_buscos)
                self.gene_details.update(self.metaeuk_runner.gene_details)
                self.sequences_aa.update(self.metaeuk_runner.sequences_aa)
                self.run_hmmer(self.metaeuk_runner.pred_protein_seqs_modified, busco_ids=incomplete_buscos)
                incomplete_buscos = (self.hmmer_runner.missing_buscos + list(self.hmmer_runner.fragmented_buscos.keys()))
            except NoRerunFile:
                if i == 1:
                    logger.info("Metaeuk rerun did not find any genes")
                else:
                    raise SystemExit("Metaeuk did not find any genes in the input file.")

        try:
            self.metaeuk_runner.combine_run_results()
        except FileNotFoundError:
            # This exception should only happen if the rerun file does not exist. If the initial run file was
            # missing there would have been a SystemExit call above. The index 0 sets the "combined" file to the
            # output of the initial run.
            self.metaeuk_runner.combined_pred_protein_seqs = self.metaeuk_runner.pred_protein_files[0]
        self.hmmer_runner.write_buscos_to_file(self.sequences_aa)


    def _run_metaeuk(self, incomplete_buscos):
        self.metaeuk_runner.configure_runner(incomplete_buscos)
        if self.restart and self.metaeuk_runner.check_previous_completed_run():
            logger.info("Skipping Metaeuk run as already run")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.metaeuk_runner.run()

        self.metaeuk_runner.get_gene_details()
        self.metaeuk_runner.edit_file_header()