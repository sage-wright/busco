# coding: utf-8
"""
BuscoRunner.py

Intermediary module for delegating control of run features to individual classes

Author(s): Matthew Berkeley

Copyright (c) 2015-2024, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

from busco.Exceptions import BatchFatalError, BuscoError
from busco.analysis.BuscoAnalysis import BuscoAnalysis
from busco.analysis.GenomeAnalysis import (
    GenomeAnalysisEukaryotesAugustus,
    GenomeAnalysisEukaryotesMetaeuk,
    GenomeAnalysisEukaryotesMiniprot,
    GenomeAnalysisProkaryotes,
)
from busco.analysis.TranscriptomeAnalysis import (
    TranscriptomeAnalysisProkaryotes,
    TranscriptomeAnalysisEukaryotes,
)
from busco.analysis.GeneSetAnalysis import GeneSetAnalysis
from busco.BuscoLogger import BuscoLogger
from busco.BuscoConfig import BuscoConfigMain
from busco.busco_tools.base import NoGenesError, BaseRunner
from busco.BuscoLogger import LogDecorator as log
from configparser import NoOptionError
from busco.busco_tools.base import ToolException
import busco
import os
import sys
import shutil
import time
import traceback
from collections import defaultdict
import json


logger = BuscoLogger.get_logger(__name__)


class SingleRunner:

    all_runners = set()

    def __init__(self, config_manager):
        self.start_time = time.time()
        self.config_manager = config_manager
        self.config = self.config_manager.config_main
        self.input_file = self.config.get("busco_run", "in")
        self.lineage_dataset = None
        self.runner = None

    def get_lineage(self):
        if self.config.getboolean("busco_run", "auto-lineage"):
            (
                lineage_dataset_fullpath,
                runner,
                parent_domain,
            ) = self.auto_select_lineage()  # full path
            self.config.set("busco_run", "lineage_dataset", lineage_dataset_fullpath)
            self.config.set("busco_run", "domain_run_name", parent_domain)
            self.runner = runner
            self.runner.set_parent_dataset()

        self.lineage_dataset = self.config.get(
            "busco_run", "lineage_dataset"
        )  # full path
        self.config.load_dataset(self.lineage_dataset)

        return

    @log("No lineage specified. Running lineage auto selector.\n", logger)
    def auto_select_lineage(self):
        from busco.AutoLineage import (
            AutoSelectLineage,
        )  # Import statement inside method to avoid circular imports

        try:
            asl = AutoSelectLineage(self.config_manager)
            asl.run_auto_selector()
            asl.get_lineage_dataset()
            asl.set_best_match_lineage()
            lineage_dataset = asl.best_match_lineage_dataset
            runner = asl.selected_runner
            parent_domain = runner.config.get("busco_run", "domain_run_name")
        finally:
            try:
                type(self).all_runners.update(asl.runners)
                asl.reset()
            except UnboundLocalError:
                pass
        return lineage_dataset, runner, parent_domain

    @staticmethod
    def log_error(err):
        logger.error(err)
        logger.debug(err, exc_info=True)
        logger.error("BUSCO analysis failed!")
        logger.error(
            "Check the logs, read the user guide (https://busco.ezlab.org/busco_userguide.html), "
            "and check the BUSCO issue board on https://gitlab.com/ezlab/busco/issues\n"
        )

    def reset(self):
        for runner in type(self).all_runners:
            runner.reset()
            runner.config.reset()
            runner.analysis.reset()
        type(self).all_runners = set()

    @log("Input file is {}", logger, attr_name="input_file")
    def run(self):
        try:
            self.get_lineage()

            lineage_basename = os.path.basename(
                self.config.get("busco_run", "lineage_dataset")
            )
            main_out_folder = self.config.get("busco_run", "main_out")
            lineage_results_folder = os.path.join(
                main_out_folder,
                "auto_lineage",
                self.config.get("busco_run", "lineage_results_dir"),
            )

            if not (
                self.config.getboolean("busco_run", "auto-lineage")
                and (
                    lineage_basename.startswith(("bacteria", "archaea", "eukaryota"))
                    or (
                        lineage_basename.startswith(
                            ("mollicutes", "mycoplasmatales", "entomoplasmatales")
                        )
                        and os.path.exists(lineage_results_folder)
                    )
                    or (self.config.get("busco_run", "domain") == "viruses")
                )
            ):
                # It is possible that the  lineages ("mollicutes", "mycoplasmatales", "entomoplasmatales") were
                # arrived at either by the Prodigal genetic code shortcut or by BuscoPlacer. If the former, the run
                # will have already been completed. If the latter it still needs to be done.
                # If the auto-lineage selected dataset is from the virus pipeline, the run has also already been
                # completed.

                self.runner = AnalysisRunner(self.config)

            type(self).all_runners.add(self.runner)

            if os.path.exists(lineage_results_folder):
                new_dest = os.path.join(
                    main_out_folder, self.config.get("busco_run", "lineage_results_dir")
                )
                if not os.path.exists(new_dest):
                    new_src = os.path.relpath(
                        lineage_results_folder, os.path.dirname(new_dest)
                    )
                    os.symlink(new_src, new_dest)
            else:
                self.runner.run_analysis()
                AnalysisRunner.selected_dataset = lineage_basename

            if self.config.getboolean("busco_run", "tar"):
                self.compress_folders()

            wall_time_elapsed = time.time() - self.start_time
            self.config.run_stats["wall_time"] = wall_time_elapsed
            self.runner.finish(wall_time_elapsed)

        except BuscoError:
            if self.runner is not None:
                type(self).all_runners.add(self.runner)
            raise

        except ToolException as e:
            raise BatchFatalError(e.value)

        except KeyboardInterrupt:
            raise BatchFatalError(
                "A signal was sent to kill the process. \nBUSCO analysis failed !"
            )

        except BaseException:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            logger.critical(
                "Unhandled exception occurred:\n{}\n".format(
                    "".join(
                        traceback.format_exception(exc_type, exc_value, exc_traceback)
                    )
                )
            )
            raise BatchFatalError(str(exc_value))

    def compress_folders(self):
        for runner in type(self).all_runners:
            folders_to_compress = [
                runner.analysis.hmmer_runner.single_copy_sequences_folder,
                runner.analysis.hmmer_runner.multi_copy_sequences_folder,
                runner.analysis.hmmer_runner.fragmented_sequences_folder,
                runner.analysis.hmmer_runner.output_folder,
            ]
            if self.config.getboolean("busco_run", "use_augustus"):
                folders_to_compress.extend(
                    [
                        runner.analysis.augustus_runner.pred_genes_dir_initial,
                        runner.analysis.augustus_runner.pred_genes_dir_rerun,
                        runner.analysis.augustus_runner.gff_dir,
                        runner.analysis.gff2gb_runner.gb_folder,
                    ]
                )
            for folder in folders_to_compress:
                try:
                    shutil.make_archive(
                        folder,
                        "gztar",
                        os.path.dirname(folder),
                        os.path.basename(folder),
                    )
                    shutil.rmtree(folder)
                except OSError:
                    raise


class BatchRunner:

    batch_results = []

    @log(
        "Running in batch mode. {} input files found in {}",
        logger,
        attr_name=["num_inputs", "input_dir"],
        on_func_exit=True,
    )
    def __init__(self, config_manager):
        self.config_manager = config_manager
        self.config = self.config_manager.config_main
        self.input_dir = self.config.get("busco_run", "in")
        self.input_files = [
            os.path.join(self.input_dir, f)
            for f in os.listdir(self.input_dir)
            if not f.startswith(".")
        ]
        self.num_inputs = len(self.input_files)

    def run(self):
        for i, input_file in enumerate(self.input_files):
            try:
                logger.info("\n\n" + "*" * 80 + "\n" + "Running BUSCO on file: {}\n".format(input_file) + "*" * 80 + "\n\n")
                input_id = os.path.basename(input_file)
                self.config.set("busco_run", "in", input_file)
                self.config.set(
                    "busco_run",
                    "main_out",
                    os.path.join(
                        self.config.get("busco_run", "out_path"),
                        self.config.get("busco_run", "out"),
                        input_id,
                    ),
                )
                single_run = SingleRunner(self.config_manager)
                single_run.run()

                run_summary = single_run.runner.format_run_summary()
                type(self).batch_results.append(run_summary)

            except NoOptionError as noe:
                raise BatchFatalError(noe.message)

            except BatchFatalError:
                raise

            except BuscoError as be:
                if "did not recognize any genes" in be.value:
                    type(self).batch_results.append(
                        "{}\tNo genes found\t\t\t\t\t\t\t\t\t\t{}{}\n".format(
                            os.path.basename(input_file),
                            "\t\t\t"
                            * int(
                                single_run.config.getboolean(
                                    "busco_run", "auto-lineage"
                                )
                            ),
                            "\t" * int(single_run.config.getboolean("busco_run", "use_miniprot")),
                        )
                    )
                else:
                    type(self).batch_results.append(
                        "{}\tRun failed; check logs\t\t\t\t\t\t\t\t\t\t{}{}\n".format(
                            os.path.basename(input_file),
                            "\t\t\t"
                            * int(
                                single_run.config.getboolean(
                                    "busco_run", "auto-lineage"
                                )
                            ),
                            "\t" * int(single_run.config.getboolean("busco_run", "use_miniprot")),
                        )
                    )
                logger.error(be.value)
                continue
            except BaseException as e:
                raise BatchFatalError(e)
            finally:
                AnalysisRunner.reset()
                BuscoLogger.reset()
                single_run.reset()

        try:
            self.write_batch_summary()
        except NoOptionError as noe:
            raise BatchFatalError(noe.message)

    def write_batch_summary(self):
        summary_file = os.path.join(
            self.config.get("busco_run", "out_path"),
            self.config.get("busco_run", "out"),
            "batch_summary.txt",
        )

        with open(summary_file, "w") as f:
            if "genome" in self.config.get("busco_run", "mode"):
                metrics_header = (
                    "\tScaffold N50\tContigs N50\tPercent gaps\tNumber of scaffolds"
                )
                if self.config.getboolean("busco_run", "use_miniprot"):
                    extra = "\tInternal stop codon percent"
                else:
                    extra = ""
            else:
                metrics_header = ""
                extra = ""
            if self.config.getboolean("busco_run", "auto-lineage"):
                f.write(
                    "Input_file\tDataset\tComplete\tSingle\tDuplicated\tFragmented\tMissing\tn_markers{}\t"
                    "Scores_archaea_odb10\tScores_bacteria_odb10\tScores_eukaryota_odb10{}\n".format(
                        extra,
                        metrics_header,
                    )
                )
            else:
                f.write(
                    "Input_file\tDataset\tComplete\tSingle\tDuplicated\tFragmented\tMissing\tn_markers{}{}\n".format(
                        extra,
                        metrics_header,
                    )
                )
            for line in type(self).batch_results:
                f.write(line)


class AnalysisRunner:

    mode_dict = {
        "euk_genome_min": GenomeAnalysisEukaryotesMiniprot,
        "euk_genome_met": GenomeAnalysisEukaryotesMetaeuk,
        "euk_genome_aug": GenomeAnalysisEukaryotesAugustus,
        "prok_genome_prod": GenomeAnalysisProkaryotes,
        "prok_genome_min": GenomeAnalysisEukaryotesMiniprot,
        "euk_tran": TranscriptomeAnalysisEukaryotes,
        "prok_tran": TranscriptomeAnalysisProkaryotes,
        "proteins": GeneSetAnalysis,
    }

    final_results = {}
    results_datasets = []
    all_results = defaultdict(dict)
    selected_dataset = None

    def __init__(self, config):

        self.config = config
        setattr(BaseRunner, "config", config)
        setattr(BuscoAnalysis, "config", config)

        self.input_file = self.config.get("busco_run", "in")
        self.mode = self.config.update_mode()
        self.domain = self.config.get("busco_run", "domain")
        self.lineage_basename = os.path.basename(
            self.config.get("busco_run", "lineage_dataset")
        )

        analysis_type = type(self).mode_dict[self.mode]
        self.analysis = analysis_type()
        self.summary = {
            "parameters": {},
            "lineage_dataset": {},
            "versions": {},
            "results": {},
        }
        self.prok_fail_count = (
            0  # Needed to check if both bacteria and archaea return no genes.
        )
        self.short_summary_file = os.path.join(
            self.analysis.run_folder, "short_summary.txt"
        )
        self.short_summary_file_json = os.path.join(
            self.analysis.run_folder, "short_summary.json"
        )
        self.cleaned_up = False

    @classmethod
    def reset(cls):
        cls.final_results = {}
        cls.results_datasets = []
        cls.all_results = defaultdict(dict)

    def save_results(self, exception_raised=False):

        if not self.config.getboolean("busco_run", "auto-lineage"):
            type(self).selected_dataset = self.lineage_basename
        try:
            self.final_results[
                self.lineage_basename
            ] = self.analysis.hmmer_runner.hmmer_results_lines
            self.all_results[self.lineage_basename][
                "one_line_summary"
            ] = self.analysis.hmmer_runner.one_line_summary_raw.strip()
            self.all_results[self.lineage_basename][
                "Complete percentage"
            ] = self.analysis.hmmer_runner.complete_percent
            self.all_results[self.lineage_basename]["Complete BUSCOs"] = (
                self.analysis.hmmer_runner.single_copy
                + self.analysis.hmmer_runner.multi_copy
            )
            self.all_results[self.lineage_basename][
                "Single copy percentage"
            ] = self.analysis.hmmer_runner.s_percent
            self.all_results[self.lineage_basename][
                "Single copy BUSCOs"
            ] = self.analysis.hmmer_runner.single_copy
            self.all_results[self.lineage_basename][
                "Multi copy percentage"
            ] = self.analysis.hmmer_runner.d_percent
            self.all_results[self.lineage_basename][
                "Multi copy BUSCOs"
            ] = self.analysis.hmmer_runner.multi_copy
            self.all_results[self.lineage_basename][
                "Fragmented percentage"
            ] = self.analysis.hmmer_runner.f_percent
            self.all_results[self.lineage_basename][
                "Fragmented BUSCOs"
            ] = self.analysis.hmmer_runner.only_fragments
            self.all_results[self.lineage_basename][
                "Missing percentage"
            ] = self.analysis.hmmer_runner.missing_percent
            self.all_results[self.lineage_basename]["Missing BUSCOs"] = len(
                self.analysis.hmmer_runner.missing_buscos
            )
            self.all_results[self.lineage_basename][
                "n_markers"
            ] = self.analysis.hmmer_runner.total_buscos
            self.all_results[self.lineage_basename]["avg_identity"] = self.analysis.hmmer_runner.avg_identity
            self.all_results[self.lineage_basename]["domain"] = self.analysis.domain
            try:  # these values will only be present in miniprot runs, and if internal stop codons are detected
                self.all_results[self.lineage_basename][
                    "internal_stop_codon_count"
                ] = self.analysis.hmmer_runner.complete_stop_codon_count
                self.all_results[self.lineage_basename][
                    "internal_stop_codon_percent"
                ] = self.analysis.hmmer_runner.e_percent
            except AttributeError:
                pass
            if "genome" in self.mode and not self.config.getboolean(
                "busco_run", "skip_bbtools"
            ):
                self.all_results[self.lineage_basename].update(
                    self.analysis.bbtools_runner.metrics
                )
            self.save_parent_dataset(self.lineage_basename)
            self.results_datasets.append(os.path.basename(self.lineage_basename))
            self.compile_summary()
            self.write_summary_files()
        except (AttributeError, TypeError):
            if exception_raised:
                try:
                    self.compile_summary()
                except:
                    pass
            else:
                raise BuscoError("Something went wrong. Results not recorded.")

    def compile_summary(self):
        for option in self.config.items("busco_run"):
            if option[0] in [
                "creation_date",
                "name",
                "number_of_buscos",
                "number_of_species",
            ]:
                self.summary["lineage_dataset"][option[0]] = option[1]
            else:
                self.summary["parameters"][option[0]] = option[1]

        self.config.run_stats["versions"] = BaseRunner.tool_versions
        self.config.run_stats["versions"]["python"] = sys.version_info
        self.config.run_stats["versions"]["busco"] = busco.__version__
        self.summary["versions"] = self.config.run_stats["versions"]

        if "genome" in self.mode:
            if self.config.get("busco_run", "domain") in ["prokaryota", "viruses"]:
                gene_predictor = "prodigal"
            elif self.config.getboolean("busco_run", "use_augustus"):
                gene_predictor = "augustus"
            elif self.config.getboolean("busco_run", "use_metaeuk"):
                gene_predictor = "metaeuk"
            elif self.config.getboolean("busco_run", "use_miniprot"):
                gene_predictor = "miniprot"
            else:
                gene_predictor = "miniprot"
            self.summary["parameters"]["gene_predictor"] = gene_predictor
            try:
                self.summary["metrics"] = self.analysis.bbtools_runner.metrics
            except AttributeError:  # may not exist if bbtools is skipped
                pass

        self.summary["results"] = type(self).all_results[self.lineage_basename]
        self.config.run_stats["results"] = type(self).all_results[self.lineage_basename]
        if (
            self.config.getboolean("busco_run", "auto-lineage")
            and isinstance(self.config, BuscoConfigMain)
            and hasattr(self.config, "placement_files")
        ):
            self.summary["placement_file_versions"] = self.config.placement_files

    def save_parent_dataset(self, lineage_basename):
        try:
            parent_dataset = self.config.get("busco_run", "parent_dataset")
            if os.path.basename(parent_dataset) == os.path.basename(
                self.config.get("busco_run", "lineage_dataset")
            ):
                pass
            else:
                self.all_results[os.path.basename(lineage_basename)][
                    "parent_dataset"
                ] = os.path.basename(parent_dataset)
        except NoOptionError:
            pass

    def set_parent_dataset(self):
        lineage_dataset = self.config.get("busco_run", "lineage_dataset")
        domain_run_name = self.config.get("busco_run", "domain_run_name")
        if not lineage_dataset.endswith(
            domain_run_name
        ):  # Don't add a parent if Busco Placer failed to place
            self.config.set(
                "busco_run",
                "parent_dataset",
                os.path.join(
                    os.path.dirname(
                        lineage_dataset,
                    ),
                    domain_run_name,
                ),
            )
        self.save_parent_dataset(lineage_dataset)

    def run_analysis(self, callback=(lambda *args: None)):
        exception_raised = True
        try:
            self.analysis.run_analysis()
            s_buscos = self.analysis.hmmer_runner.single_copy
            d_buscos = self.analysis.hmmer_runner.multi_copy
            f_buscos = self.analysis.hmmer_runner.only_fragments
            s_percent = self.analysis.hmmer_runner.s_percent
            d_percent = self.analysis.hmmer_runner.d_percent
            f_percent = self.analysis.hmmer_runner.f_percent
            self.cleanup()

        except NoGenesError as nge:
            no_genes_msg = (
                "{0} did not recognize any genes matching the dataset {1} in the input file. "
                "If this is unexpected, check your input file and your "
                "installation of {0}\n".format(
                    nge.gene_predictor, self.analysis.lineage_name
                )
            )
            fatal = (
                isinstance(self.config, BuscoConfigMain)
                or (
                    self.config.getboolean("busco_run", "auto-lineage-euk")
                    and self.mode.startswith("euk_genome")
                )
                or (
                    self.config.getboolean("busco_run", "auto-lineage-prok")
                    and self.mode.startswith("prok_genome")
                )
                and self.prok_fail_count == 1
            )
            if fatal:
                raise BuscoError(no_genes_msg)
            else:
                logger.warning(no_genes_msg)
                s_buscos = d_buscos = f_buscos = s_percent = d_percent = f_percent = 0.0
                if self.mode.startswith("prok_genome"):
                    self.prok_fail_count += 1

        except (BuscoError, BatchFatalError):
            self.cleanup()
            raise
        else:
            exception_raised = False

        finally:
            self.save_results(exception_raised)
        return callback(s_buscos, d_buscos, f_buscos, s_percent, d_percent, f_percent)

    def cleanup(self):
        self.analysis.cleanup()
        self.cleaned_up = True

    def report_tool_versions(self):
        lines = ["\nDependencies and versions:"]
        for key, value in self.analysis.hmmer_runner.tool_versions.items():
            lines.append("\t{}: {}".format(key, value))
        return lines

    def write_summary_files(self):
        self.write_txt_summary()
        self.write_json_summary()

    def write_json_summary(self):

        with open(self.short_summary_file_json, "w") as summary_file_json:
            json.dump(self.summary, summary_file_json, indent=4)

        return

    def write_txt_summary(self):
        with open(self.short_summary_file, "w") as summary_file:

            self.analysis.hmmer_runner.write_output_header(
                summary_file, no_table_header=True
            )
            summary_file.write(
                "# Summarized benchmarking in BUSCO notation for file {}\n"
                "# BUSCO was run in mode: {}\n".format(self.input_file, self.mode)
            )
            if "genome" in self.mode:
                summary_file.write(
                    "# Gene predictor used: {}\n".format(
                        self.summary["parameters"]["gene_predictor"]
                    )
                )
            summary_file.write("\n")

            for line in self.analysis.hmmer_runner.hmmer_results_lines:
                summary_file.write("\t{}".format(line))

            if "genome" in self.mode:
                summary_file.write("\nAssembly Statistics:\n")
                for key, val in self.analysis.bbtools_runner.metrics.items():
                    if "N50" in key:
                        if int(val) > 1000000:
                            val = "{} MB".format(int(val) // 1000000)
                        elif int(val) > 1000:
                            val = "{} KB".format(int(val) // 1000)
                    summary_file.write("\t{}\t{}\n".format(val, key))
                summary_file.write("\n")

            tool_versions_lines = self.report_tool_versions()
            for line in tool_versions_lines:
                summary_file.write(line + "\n")

            if (
                self.config.getboolean("busco_run", "auto-lineage")
                and isinstance(self.config, BuscoConfigMain)
                and hasattr(self.config, "placement_files")
            ):
                summary_file.write("\nPlacement file versions:\n")
                for placement_file in self.config.placement_files:
                    summary_file.write("\t{}\n".format(placement_file))

    def format_results(self):
        framed_output = []
        final_dataset = type(self).selected_dataset
        domain = type(self).all_results[final_dataset]["domain"]
        try:
            parent_dataset = os.path.basename(
                type(self).all_results[final_dataset]["parent_dataset"]
            )
        except KeyError:
            parent_dataset = None
        if parent_dataset is None or domain == "viruses":
            header1 = "Results from dataset {}\n".format(final_dataset)
            final_output_results1 = "".join(
                self._check_parasitic(type(self).final_results[final_dataset][1:])
            )
        else:
            header1 = "Results from generic domain {}\n".format(parent_dataset)
            final_output_results1 = "".join(
                self._check_parasitic(type(self).final_results[parent_dataset][1:])
            )

        sb1 = SmartBox()
        framed_lines1 = sb1.create_results_box(header1, final_output_results1)
        framed_output.append(framed_lines1)

        if domain != "viruses" and parent_dataset and parent_dataset != final_dataset:
            header2 = "Results from dataset {}\n".format(final_dataset)
            final_output_results2 = "".join(
                self._check_parasitic(type(self).final_results[final_dataset][1:])
            )
            sb2 = SmartBox()
            framed_lines2 = sb2.create_results_box(header2, final_output_results2)
            framed_output.append(framed_lines2)

        return "".join(framed_output)

    def format_run_summary(self):
        input_file = os.path.basename(self.config.get("busco_run", "in"))
        dataset = type(self).selected_dataset
        complete = type(self).all_results[dataset]["Complete percentage"]
        single = type(self).all_results[dataset]["Single copy percentage"]
        duplicated = type(self).all_results[dataset]["Multi copy percentage"]
        fragmented = type(self).all_results[dataset]["Fragmented percentage"]
        missing = type(self).all_results[dataset]["Missing percentage"]
        n_markers = type(self).all_results[dataset]["n_markers"]

        if "genome" in self.mode:
            try:
                scaffold_n50 = self.analysis.bbtools_runner.metrics["Scaffold N50"]
                contigs_n50 = self.analysis.bbtools_runner.metrics["Contigs N50"]
                percent_gaps = self.analysis.bbtools_runner.metrics["Percent gaps"]
                num_scaffolds = self.analysis.bbtools_runner.metrics["Number of scaffolds"]
                metrics_scores = "\t{}\t{}\t{}\t{}".format(
                    scaffold_n50, contigs_n50, percent_gaps, num_scaffolds
                )
            except KeyError:
                metrics_scores = "\t\t\t\t"
        else:
            metrics_scores = "\t\t\t\t"

        if self.config.getboolean("busco_run", "use_miniprot"):
            internal_stop_codon_percent = self.analysis.hmmer_runner.e_percent
            extra = "\t{}".format(internal_stop_codon_percent)
        else:
            extra = ""

        if self.config.getboolean("busco_run", "auto-lineage"):
            try:
                scores_arch = (
                    type(self).all_results["archaea_odb10"]["one_line_summary"].strip()
                )
            except KeyError:
                scores_arch = "Not run"
            try:
                scores_bact = (
                    type(self).all_results["bacteria_odb10"]["one_line_summary"].strip()
                )
            except KeyError:
                scores_bact = "Not run"
            try:
                scores_euk = (
                    type(self)
                    .all_results["eukaryota_odb10"]["one_line_summary"]
                    .strip()
                )
            except KeyError:
                scores_euk = "Not run"

            parent_lineage_scores = "\t{}\t{}\t{}".format(
                scores_arch, scores_bact, scores_euk
            )
        else:
            parent_lineage_scores = ""

        summary_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}{}{}\n".format(
            input_file,
            dataset,
            complete,
            single,
            duplicated,
            fragmented,
            missing,
            n_markers,
            extra,
            parent_lineage_scores,
            metrics_scores,
        )
        return summary_line

    def _check_parasitic(self, final_output_results):
        try:
            with open(
                os.path.join(self.analysis.lineage_dataset, "missing_in_parasitic.txt")
            ) as parasitic_file:
                missing_in_parasitic_buscos = [
                    entry.strip() for entry in parasitic_file.readlines()
                ]
            if (
                len(self.analysis.hmmer_runner.missing_buscos)
                >= 0.8 * len(missing_in_parasitic_buscos)
                and len(missing_in_parasitic_buscos) > 0
            ):
                intersection = [
                    mb
                    for mb in self.analysis.hmmer_runner.missing_buscos
                    if mb in missing_in_parasitic_buscos
                ]
                percent_missing_in_parasites = round(
                    100
                    * len(intersection)
                    / len(self.analysis.hmmer_runner.missing_buscos),
                    1,
                )
                if percent_missing_in_parasites >= 80.0:
                    corrected_summary = self._recalculate_parasitic_scores(
                        len(missing_in_parasitic_buscos)
                    )
                    positive_parasitic_line = (
                        "\n!!! The missing BUSCOs match the pattern of a parasitic-reduced "
                        "genome. {}% of your missing BUSCOs are typically missing in these. "
                        "A corrected score would be: \n{}\n".format(
                            percent_missing_in_parasites, corrected_summary
                        )
                    )
                    final_output_results.append(positive_parasitic_line)
                    if not self.config.getboolean("busco_run", "auto-lineage"):
                        auto_lineage_line = "\nConsider using the auto-lineage mode to select a more specific lineage."
                        final_output_results.append(auto_lineage_line)
                    with open(self.short_summary_file, "a") as short_summary_file:
                        short_summary_file.write(positive_parasitic_line)

        except OSError:
            pass

        return final_output_results

    def _recalculate_parasitic_scores(self, num_missing_in_parasitic):
        total_buscos = (
            self.analysis.hmmer_runner.total_buscos - num_missing_in_parasitic
        )
        single_copy = self.analysis.hmmer_runner.single_copy
        multi_copy = self.analysis.hmmer_runner.multi_copy
        fragmented_copy = self.analysis.hmmer_runner.only_fragments
        s_percent = abs(round(100 * single_copy / total_buscos, 1))
        d_percent = abs(round(100 * multi_copy / total_buscos, 1))
        f_percent = abs(round(100 * fragmented_copy / total_buscos, 1))

        one_line_summary = "C:{}%[S:{}%,D:{}%],F:{}%,M:{}%,n:{}\t\n".format(
            round(s_percent + d_percent, 1),
            s_percent,
            d_percent,
            f_percent,
            round(100 - s_percent - d_percent - f_percent, 1),
            total_buscos,
        )
        return one_line_summary

    @staticmethod
    def copy_summary_file(from_dir, main_dir, domain_dir, tag, ext):
        shutil.copyfile(
            os.path.join(from_dir, "short_summary.{}".format(ext)),
            os.path.join(
                main_dir,
                "short_summary.{}.{}.{}.{}".format(
                    tag, domain_dir.replace("run_", ""), os.path.basename(main_dir), ext
                ),
            ),
        )

    def organize_final_output(self):
        main_out_folder = self.config.get("busco_run", "main_out")

        try:
            domain_results_folder = self.config.get("busco_run", "domain_run_name")
            root_domain_output_folder = os.path.join(
                main_out_folder, "auto_lineage", "run_{}".format(domain_results_folder)
            )
            root_domain_output_folder_final = os.path.join(
                main_out_folder, "run_{}".format(domain_results_folder)
            )
            new_src = os.path.relpath(
                root_domain_output_folder,
                os.path.dirname(root_domain_output_folder_final),
            )
            os.symlink(new_src, root_domain_output_folder_final)
            self.copy_summary_file(
                root_domain_output_folder_final,
                main_out_folder,
                domain_results_folder,
                "generic",
                "txt",
            )
            self.copy_summary_file(
                root_domain_output_folder_final,
                main_out_folder,
                domain_results_folder,
                "generic",
                "json",
            )

        except NoOptionError:
            pass

        except OSError:  # FileExistsError in restart mode
            pass

        finally:
            lineage_results_folder = self.config.get("busco_run", "lineage_results_dir")
            lineage_results_path = os.path.join(main_out_folder, lineage_results_folder)
            self.copy_summary_file(
                lineage_results_path,
                main_out_folder,
                lineage_results_folder,
                "specific",
                "txt",
            )
            self.copy_summary_file(
                lineage_results_path,
                main_out_folder,
                lineage_results_folder,
                "specific",
                "json",
            )
        return

    @staticmethod
    def move_log_file(config):
        # This is deliberately a staticmethod so it can be called from run_BUSCO() even if BuscoRunner has not yet
        # been initialized.
        try:
            log_folder = os.path.join(
                config.get("busco_run", "out_path"),
                config.get("busco_run", "out"),
                "logs",
            )
            if not os.path.exists(log_folder):
                os.makedirs(log_folder)
            shutil.move(
                "busco_{}.log".format(BuscoLogger.pid),
                os.path.join(log_folder, "busco.log"),
            )
        except OSError:
            logger.warning(
                "Unable to move 'busco_{}.log' to the 'logs' folder.".format(
                    BuscoLogger.pid
                )
            )
        return

    def finish(self, elapsed_time):

        final_output_results = self.format_results()
        logger.info("".join(final_output_results))

        self.organize_final_output()

        if not logger.has_warning():
            logger.info(
                "BUSCO analysis done. Total running time: {} seconds".format(
                    str(round(elapsed_time))
                )
            )
        else:
            logger.info(
                "BUSCO analysis done with WARNING(s). Total running time: {} seconds\n\n"
                "***** Summary of warnings: *****".format(str(round(elapsed_time)))
            )
            for item in type(logger).warn_output.getvalue().split("\n"):
                print(item)

        logger.info("Results written in {}".format(self.analysis.main_out))
        logger.info(
            "For assistance with interpreting the results, please consult the userguide: "
            "https://busco.ezlab.org/busco_userguide.html\n"
        )
        logger.info(
            "Visit this page https://gitlab.com/ezlab/busco#how-to-cite-busco to see how to cite BUSCO"
        )


class SmartBox:
    def __init__(self):
        self.width = 50
        self.tab = 4 * " "

    def define_width(self, header_text, body_text):
        lines = body_text.split("\n")
        lens = [len(x) for x in lines]
        max_ind = lens.index(max(lens))
        longest_line = lines[max_ind]
        scores_line = lines[0]
        if len(header_text) > len(longest_line):
            longest_line = header_text
        if len(longest_line) + 2 < 80:  # add 2 for vertical frame lines
            self.width = max(max(50, len(scores_line) + 2), len(longest_line.expandtabs()))
        else:
            self.width = max(50, len(scores_line) + 2)  # scores lines with "E:x%" values can be longer than 50

    def wrap_header(self, header_text):
        header_text = header_text.strip().replace("\t", self.tab)
        if len(header_text) > 80:
            header_text = self.wrap_long_line(header_text)
        return header_text

    def wrap_long_line(self, line):
        words = line.split(" ")
        word_num = 0
        word_start = 0
        length = 0
        output_lines = []
        while word_num < len(words):
            while length < self.width:
                word_num += 1
                line = " ".join(words[word_start:word_num])
                length = len(line.expandtabs())
                if length >= self.width:
                    word_num -= 1
                    line = " ".join(words[word_start:word_num])
                    break
                if word_num == len(words):
                    break
            output_lines.append(line)
            length = 0
            word_start = word_num
        return "\n".join(output_lines)

    def wrap_text(self, text):
        lines = text.split("\n")
        output_lines = []
        for line in lines:
            line = line.strip().replace("\t", self.tab)
            if (
                len(line.expandtabs()) + 2 < self.width
            ):  # add 2 for vertical frame lines
                output_lines.append(line)
            else:
                output_lines.append(self.wrap_long_line(line))
        return "\n".join(output_lines)

    def add_vertical(self, lines):
        if isinstance(lines, str):
            lines = lines.strip().split("\n")
        formatted_lines = []
        for line in lines:
            line = "|{}".format(
                line.strip()
            )  # left bar needs to be added before expanding tabs
            whitespace = " " * (self.width - len(line.expandtabs()))
            format_line = "{}{}|".format(line, whitespace)
            formatted_lines.append(format_line)
        return formatted_lines

    def add_horizontal(self):
        return "-" * self.width

    def create_results_box(self, header_text, body_text):
        self.define_width(header_text, body_text)  # Called first to define width
        header = self.wrap_header(header_text)
        box_lines = list(["\n"])
        box_lines.append("{}{}".format(self.tab, self.add_horizontal()))
        framed_header = self.add_vertical(header)
        for line in framed_header:
            box_lines.append("{}{}".format(self.tab, line))
        box_lines.append("{}{}".format(self.tab, self.add_horizontal()))
        body = self.wrap_text(body_text)
        framed_body = self.add_vertical(body)
        for line in framed_body:
            box_lines.append("{}{}".format(self.tab, line))
        box_lines.append("{}{}".format(self.tab, self.add_horizontal()))
        return "\n".join(box_lines)