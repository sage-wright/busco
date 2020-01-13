#!/usr/bin/env python
# coding: utf-8
"""
.. module:: GenomeAnalysis
   :synopsis: GenomeAnalysis implements genome analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 3.0.0

Copyright (c) 2016-2017, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
from busco.BuscoAnalysis import BuscoAnalysis
from busco.Analysis import NucleotideAnalysis
from busco.BuscoTools import ProdigalRunner, AugustusRunner, GFF2GBRunner, NewSpeciesRunner, ETrainingRunner, OptimizeAugustusRunner
from busco.BuscoConfig import BuscoConfigAuto
import os
import shutil
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.Toolset import Tool
import time
from abc import ABCMeta, abstractmethod
from configparser import NoOptionError



logger = BuscoLogger.get_logger(__name__)


class GenomeAnalysis(NucleotideAnalysis, BuscoAnalysis, metaclass=ABCMeta):

    _mode = "genome"

    def __init__(self, config):
        self._target_species = config.get("busco_run", "augustus_species")
        super().__init__(config)

    @abstractmethod
    def run_analysis(self):
        super().run_analysis()


    @abstractmethod
    def create_dirs(self):
        super().create_dirs()

    def check_tool_dependencies(self):
        """
        check dependencies on tools
        :raises SystemExit: if a Tool is not available
        """
        super().check_tool_dependencies()

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

    def set_rerun_busco_command(self, clargs):
        """
        This function sets the command line to call to reproduce this run
        """
        clargs.extend(["-sp", self._target_species])
        super().set_rerun_busco_command(clargs)

    def _write_full_table_header(self, out):
        """
        This function adds a header line to the full table file
        :param out: a full table file
        :type out: file
        """
        out.write("# Busco id\tStatus\tContig\tStart\tEnd\tScore\tLength\n")


class GenomeAnalysisProkaryotes(GenomeAnalysis):
    """
    This class runs a BUSCO analysis on a genome.
    """

    def __init__(self, config):
        """
        Initialize an instance.
        :param config: Values of all parameters that have to be defined
        :type config: PipeConfig
        """
        super().__init__(config)
        self.load_persistent_tools()

        # Get genetic_code from dataset.cfg file
        # bacteria/archaea=11; Entomoplasmatales,Mycoplasmatales=4
        try:
            self._genetic_code = self._config.get("prodigal", "prodigal_genetic_code").split(",")
        except NoOptionError:
            self._genetic_code = ["11"]

        if len(self._genetic_code) > 1:
            try:
                self.ambiguous_cd_range = [float(self._config.get("prodigal", "ambiguous_cd_range_lower")),
                                           float(self._config.get("prodigal", "ambiguous_cd_range_upper"))]
            except NoOptionError:
                raise SystemExit("Dataset config file does not contain required information. Please upgrade datasets.")

        else:
            self.ambiguous_cd_range = [None, 0]

        self.code_4_selected = False
        self.prodigal_output_dir = os.path.join(self.main_out, "prodigal_output")

    def _cleanup(self):
        # tmp_path = os.path.join(self.prodigal_output_dir, "tmp")
        # if os.path.exists(tmp_path):
        #     shutil.rmtree(tmp_path)
        super()._cleanup()

    def run_analysis(self):
        """
        This function calls all needed steps for running the analysis.
        """
        # Initialize tools and check dependencies
        super().run_analysis()

        if not os.path.exists(self.prodigal_output_dir):  # If prodigal has already been run on the input, don't run it again
            os.makedirs(self.prodigal_output_dir)
            self._run_prodigal()
            self._config.persistent_tools.append(self.prodigal_runner)

        elif any(g not in self.prodigal_runner.genetic_code for g in self._genetic_code):
            self.prodigal_runner.genetic_code = self._genetic_code
            self.prodigal_runner.cd_lower, self.prodigal_runner.cd_upper = self.ambiguous_cd_range
            self._run_prodigal()

        else:
            # Prodigal has already been run on input. Don't run again, just load necessary params.
            # First determine which GC to use
            self.prodigal_runner.select_optimal_results(self._genetic_code, self.ambiguous_cd_range)
            tmp_file = self.prodigal_runner.gc_run_results[self.prodigal_runner.gc]["tmp_name"]
            log_file = self.prodigal_runner.gc_run_results[self.prodigal_runner.gc]["log_file"]
            self.prodigal_runner._organize_prodigal_files(tmp_file, log_file)

        self.code_4_selected = self.prodigal_runner.gc == "4"
        self.sequences_nt = self.prodigal_runner.gc_run_results[self.prodigal_runner.gc]["seqs_nt"]
        self.sequences_aa = self.prodigal_runner.gc_run_results[self.prodigal_runner.gc]["seqs_aa"]
        self._gene_details = self.prodigal_runner.gc_run_results[self.prodigal_runner.gc]["gene_details"]
        self.run_hmmer(self.prodigal_runner.output_faa)
        self._write_buscos_to_file(self.sequences_aa, self.sequences_nt)
        return

    def load_persistent_tools(self):
        """
        For multiple runs, load Prodigal Runner in the same state as the previous run, to avoid having to run Prodigal
        on the input again.
        :return:
        """
        for tool in self._config.persistent_tools:
            if isinstance(tool, ProdigalRunner):
                self.prodigal_runner = tool
            else:
                raise SystemExit("Unrecognized persistent tool.")

    def create_dirs(self):
        super().create_dirs()

    def check_tool_dependencies(self):
        """
        check dependencies on tools
        :raises SystemExit: if a Tool is not available
        """
        super().check_tool_dependencies()

    def init_tools(self):
        """
        Init the tools needed for the analysis
        """
        super().init_tools()
        try:
            assert(isinstance(self._prodigal_tool, Tool))
        except AttributeError:
            self._prodigal_tool = Tool("prodigal", self._config)
        except AssertionError:
            raise SystemExit("Prodigal should be a tool")

    @log("***** Run Prodigal on input to predict and extract genes *****", logger)
    def _run_prodigal(self):
        """
        Run Prodigal on input file to detect genes.
        :return:
        """
        if not hasattr(self, "prodigal_runner"):
            self.prodigal_runner = ProdigalRunner(self._prodigal_tool, self._input_file, self.prodigal_output_dir,
                                                  self._genetic_code, self.ambiguous_cd_range, self.log_folder)
        self.prodigal_runner.run()
        self.code_4_selected = self.prodigal_runner.code_4_selected
        return

    def _write_full_table_header(self, out):
        """
        This function adds a header line to the full table file
        :param out: a full table file
        :type out: file
        """
        out.write("# Busco id\tStatus\tContig\tStart\tEnd\tScore\tLength\n")


class GenomeAnalysisEukaryotes(GenomeAnalysis):
    """
    This class runs a BUSCO analysis on a euk_genome.
    Todo: reintroduce restart mode with checkpoints
    """
    def __init__(self, config):
        """
        Retrieve the augustus config path, mandatory for genome
        Cannot be specified through config because some augustus perl scripts use it as well
        BUSCO could export it if absent, but do not want to mess up with the user env,
        let's just tell the user to do it for now.

        :param config: Values of all parameters that have to be defined
        :type config: PipeConfig
        """
        self._augustus_config_path = os.environ.get("AUGUSTUS_CONFIG_PATH")
        try:
            self._augustus_parameters = config.get("busco_run", "augustus_parameters").replace(',', ' ')
        except NoOptionError:
            self._augustus_parameters = ""
        super().__init__(config)
        self._check_file_dependencies()
        self.mkblast_runner = None
        self.tblastn_runner = None
        self.augustus_runner = None
        self.sequences_nt = {}
        self.sequences_aa = {}

    def create_dirs(self):
        super().create_dirs()

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
        super().check_tool_dependencies()

    def init_tools(self):
        """
        Initialize all required tools for Genome Eukaryote Analysis:
        MKBlast, TBlastn, Augustus and Augustus scripts: GFF2GBSmallDNA, new_species, etraining
        :return:
        """
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
        try:
            assert(isinstance(self._augustus_tool, Tool))
        except AttributeError:
            self._augustus_tool = Tool("augustus", self._config, augustus_out=True)
            # For some reason Augustus appears to send a return code before it writes to stdout, so we have to
            # sleep briefly to allow the output to be written to the file. Otherwise we have a truncated output which
            # will cause an error.
            # self._augustus_tool.sleep = 0.4
        except AssertionError:
            raise SystemExit("Augustus should be a tool")

        try:
            assert(isinstance(self._gff2gbSmallDNA_tool, Tool))
        except AttributeError:
            self._gff2gbSmallDNA_tool = Tool("gff2gbSmallDNA.pl", self._config)
        except AssertionError:
            raise SystemExit("gff2gbSmallDNA.pl should be a tool")

        try:
            assert(isinstance(self._new_species_tool, Tool))
        except AttributeError:
            self._new_species_tool = Tool("new_species.pl", self._config)
        except AssertionError:
            raise SystemExit("new_species.pl should be a tool")

        try:
            assert(isinstance(self._etraining_tool, Tool))
        except AttributeError:
            self._etraining_tool = Tool("etraining", self._config)
        except AssertionError:
            raise SystemExit("etraining should be a tool")

        if self._long:
            try:
                assert (isinstance(self._optimize_augustus_tool, Tool))
            except AttributeError:
                self._optimize_augustus_tool = Tool("optimize_augustus.pl", self._config)
            except AssertionError:
                raise SystemExit("optimize_augustus should be a tool")

        return

    @log("Running Augustus gene predictor on BLAST search results.", logger)
    def _run_augustus(self, coords):
        output_dir = os.path.join(self.run_folder, "augustus_output")
        if not os.path.exists(output_dir):  # TODO: consider grouping all create_dir calls into one function for all tools
            os.mkdir(output_dir)
        self.augustus_runner = AugustusRunner(self._augustus_tool, output_dir, self.tblastn_runner.output_seqs, self._target_species,
                                              self._lineage_dataset, self._augustus_parameters, coords,
                                              self._cpus, self.log_folder, self.sequences_aa, self.sequences_nt)
        self.augustus_runner.run()
        self.sequences_nt = self.augustus_runner.sequences_nt
        self.sequences_aa = self.augustus_runner.sequences_aa

    def _rerun_augustus(self, coords):
        self._augustus_tool.total = 0  # Reset job count
        self._augustus_tool.nb_done = 0
        missing_and_fragmented_buscos = self.hmmer_runner.missing_buscos + list(
            self.hmmer_runner.fragmented_buscos.keys())
        logger.info("Re-running Augustus with the new metaparameters, number of target BUSCOs: {}".format(
            len(missing_and_fragmented_buscos)))
        missing_and_fragmented_coords = {busco: coords[busco] for busco in coords if busco in missing_and_fragmented_buscos}
        logger.debug('Trained species folder is {}'.format(self._target_species))
        self._run_augustus(missing_and_fragmented_coords)
        return

    def _set_checkpoint(self, id=None):
        """
        This function update the checkpoint file with the provided id or delete
        it if none is provided
        :param id: the id of the checkpoint
        :type id: int
        """
        checkpoint_filename = os.path.join(self.run_folder, "checkpoint.tmp")
        if id:
            with open(checkpoint_filename, "w") as checkpt_file:
                checkpt_file.write("{}.{}".format(id, self._mode))
        else:
            if os.path.exists(checkpoint_filename):
                os.remove(checkpoint_filename)
        return

    def _run_gff2gb(self):
        self.gff2gb = GFF2GBRunner(self._gff2gbSmallDNA_tool, self.run_folder, self._input_file,
                                   self.hmmer_runner.single_copy_buscos, self._cpus)
        self.gff2gb.run()
        return

    def _run_new_species(self):
        new_species_name = "BUSCO_{}".format(os.path.basename(self.main_out))
        self.new_species_runner = NewSpeciesRunner(self._new_species_tool, self._domain, new_species_name, self._cpus)
        # create new species config file from template
        self.new_species_runner.run()
        return new_species_name

    def _run_etraining(self):
        # train on new training set (complete single copy buscos)
        self.etraining_runner = ETrainingRunner(self._etraining_tool, self.main_out, self._cpus)
        self.etraining_runner.run()
        return

    def run_analysis(self):
        """
                This function calls all needed steps for running the analysis.
                Todo: reintroduce checkpoints and restart option.
        """

        super().run_analysis()
        self._run_mkblast()
        coords = self._run_tblastn()
        self._run_augustus(coords)
        self._gene_details = self.augustus_runner.gene_details
        self.run_hmmer(self.augustus_runner.output_sequences)
        if self.busco_type == "main":
            self.rerun_analysis()

    def rerun_analysis(self):

        # self._fix_restart_augustus_folder()  # todo: reintegrate this when checkpoints are restored
        coords = self._run_tblastn(missing_and_frag_only=True, ancestral_variants=self._has_variants_file)

        logger.info("Training Augustus using Single-Copy Complete BUSCOs:")
        logger.info("Converting predicted genes to short genbank files")

        self._run_gff2gb()

        logger.info("All files converted to short genbank files, now running the training scripts")
        new_species_name = self._run_new_species()
        self._target_species = new_species_name

        self._merge_gb_files()

        self._run_etraining()

        if self._long:
            self._run_optimize_augustus(new_species_name)
            self._run_etraining()

        self._rerun_augustus(coords)
        self._gene_details.update(self.augustus_runner.gene_details)
        self.run_hmmer(self.augustus_runner.output_sequences)
        self._write_buscos_to_file(self.sequences_aa, self.sequences_nt)

        self._move_retraining_parameters()  # todo: clean species folder on systemexit
        # if self._tarzip:
        #     self._run_tarzip_augustus_output()
        #     self._run_tarzip_hmmer_output()
        # remove the checkpoint, run is done
        self._set_checkpoint()
        return

    def _check_file_dependencies(self):  # todo: currently only implemented for GenomeAnalysisEukaryotes, checking Augustus dirs. Does it need to be rolled out for all analyses?
        """
        check dependencies on files and folders
        properly configured.
        :raises SystemExit: if Augustus config path is not writable or
        not set at all
        :raises SystemExit: if Augustus config path does not contain
        the needed species
        present
        """
        try:
            augustus_species_dir = os.path.join(self._augustus_config_path, "species")
            if not os.access(augustus_species_dir, os.W_OK):
                raise SystemExit("Cannot write to Augustus species folder, please make sure you have write "
                                 "permissions to {}".format(augustus_species_dir))

        except TypeError:
            raise SystemExit(
                "The environment variable AUGUSTUS_CONFIG_PATH is not set")


        if not os.path.exists(os.path.join(augustus_species_dir, self._target_species)):
            raise SystemExit(
                "Impossible to locate the species \"{0}\" in Augustus species folder"
                " ({1}), check that AUGUSTUS_CONFIG_PATH is properly set"
                " and contains this species. \n\t\tSee the help if you want "
                "to provide an alternative species".format(self._target_species, augustus_species_dir))

    def set_rerun_busco_command(self, clargs):
        """
        This function sets the command line to call to reproduce this run
        """
        clargs.extend(["-sp", self._target_species])
        if self._augustus_parameters:
            clargs.extend(["--augustus_parameters", "\"%s\"" % self._augustus_parameters])
        super().set_rerun_busco_command(clargs)

    def _cleanup(self):
        """
        This function cleans temporary files
        """
        try:
            augustus_tmp = self.augustus_runner.tmp_dir
            if os.path.exists(augustus_tmp):
                shutil.rmtree(augustus_tmp)
        except:
            pass
        super()._cleanup()


    def _fix_restart_augustus_folder(self):
        """
        This function resets and checks the augustus folder to make a restart
        possible in phase 2
        :raises SystemExit: if it is not possible to fix the folders
        # Todo: reintegrate this when restart option is added back
        """
        if os.path.exists(os.path.join(self.augustus_runner.output_folder, "predicted_genes_run1")) \
                and os.path.exists(os.path.join(self.main_out, "hmmer_output_run1")):
            os.remove(os.path.join(self.main_out, "augustus_output", "predicted_genes", "*"))
            os.rmdir(os.path.join(self.main_out, "augustus_output", "predicted_genes"))

            os.rename(os.path.join(self.main_out, "augustus_output", "predicted_genes_run1"),
                      os.path.join(self.main_out, "augustus_output", "predicted_genes"))

            os.remove(os.path.join(self.main_out, "hmmer_output", "*"))
            os.rmdir(os.path.join(self.main_out, "hmmer_output"))

            os.rename(os.path.join(self.main_out, "hmmer_output_run1"), os.path.join(self.main_out, "hmmer_output"))


        elif (os.path.exists(os.path.join(self.main_out, "augustus_output", "predicted_genes"))
              and os.path.exists(os.path.join(self.main_out, "hmmer_output"))):
            pass
        else:
            raise SystemExit("Impossible to restart the run, necessary folders are missing. Use the -f option instead of -r")
        return

    def _move_retraining_parameters(self):
        """
        This function moves retraining parameters from augustus species folder
        to the run folder
        """
        augustus_species_path = os.path.join(self._augustus_config_path, "species", self._target_species)
        if os.path.exists(augustus_species_path):
            new_path = os.path.join(self.augustus_runner.output_folder, "retraining_parameters", self._target_species)
            shutil.move(augustus_species_path, new_path)
        else:
            logger.warning("Augustus did not produce a retrained species folder.")
        return

    def _merge_gb_files(self):
        logger.debug("concat all gb files...")
        # Concatenate all GB files into one large file
        with open(os.path.join(self.augustus_runner.output_folder, "training_set.db"), "w") as outfile:
            gb_dir_path = os.path.join(self.augustus_runner.output_folder, "gb")
            for fname in os.listdir(gb_dir_path):
                with open(os.path.join(gb_dir_path, fname), "r") as infile:
                    outfile.writelines(infile.readlines())
        return

    def _run_optimize_augustus(self, new_species_name):
        # long mode (--long) option - runs all the Augustus optimization
        # scripts (adds ~1 day of runtime)
        logger.warning("Optimizing augustus metaparameters, this may take a very long time, started at {}".format(
            time.strftime("%m/%d/%Y %H:%M:%S")))
        self.optimize_augustus_runner = OptimizeAugustusRunner(self._optimize_augustus_tool, self.augustus_runner.output_folder, new_species_name, self._cpus)
        self.optimize_augustus_runner.run()
        return
