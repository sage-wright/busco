#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: run_BUSCO
   :synopsis:
.. versionadded:: 3.0.0
.. versionchanged:: 5.4.0

BUSCO - Benchmarking Universal Single-Copy Orthologs.
This is the BUSCO main script.

To get help, ``busco -h``. See also the user guide.

And visit our website `<http://busco.ezlab.org/>`_

Copyright (c) 2016-2022, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""


import argparse
from argparse import RawTextHelpFormatter
import busco
from busco.BuscoRunner import AnalysisRunner, BatchRunner, SingleRunner
from busco.Exceptions import BatchFatalError, BuscoError
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.ConfigManager import BuscoConfigManager
from busco.Actions import (
    ListLineagesAction,
    CleanHelpAction,
    CleanVersionAction,
    DirectDownload,
)
from busco.ConfigManager import BuscoConfigMain

import sys
import time


logger = BuscoLogger.get_logger(__name__)


@log(
    "***** Start a BUSCO v{} analysis, current time: {} *****".format(
        busco.__version__, time.strftime("%m/%d/%Y %H:%M:%S")
    ),
    logger,
)
class BuscoMaster:
    def __init__(self, params):
        self.params = params
        self.config_manager = BuscoConfigManager(self.params)
        self.config = self.config_manager.config_main

    def load_config(self):
        """
        Load a busco config file that will figure out all the params from all sources
        i.e. provided config file, dataset cfg, and user args
        """
        self.config_manager.load_busco_config_main()
        self.config = self.config_manager.config_main

    def check_batch_mode(self):
        return self.config.getboolean("busco_run", "batch_mode")

    def run(self):
        try:
            self.load_config()
            runner = (
                BatchRunner(self.config_manager)
                if self.check_batch_mode()
                else SingleRunner(self.config_manager)
            )
            runner.run()

        except BuscoError as be:
            SingleRunner.log_error(be)
            raise SystemExit(1)

        except BatchFatalError as bfe:
            SingleRunner.log_error(bfe)
            raise SystemExit(1)

        finally:
            try:
                AnalysisRunner.move_log_file(self.config)
            except:
                pass


@log("Command line: {}".format(" ".join(sys.argv[:])), logger, debug=True)
def _parse_args():
    """
    This function parses the arguments provided by the user
    :return: a dictionary having a key for each arguments
    :rtype: dict
    """

    parser = argparse.ArgumentParser(
        description="Welcome to BUSCO {}: the Benchmarking Universal Single-Copy Ortholog assessment tool.\n"
        "For more detailed usage information, please review the README file provided with "
        "this distribution and the BUSCO user guide. "
        "Visit this page https://gitlab.com/ezlab/busco#how-to-cite-busco to see how to cite BUSCO".format(
            busco.__version__
        ),
        usage="busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]",
        formatter_class=RawTextHelpFormatter,
        add_help=False,
    )

    optional = parser.add_argument_group("optional arguments")

    optional.add_argument(
        "-i",
        "--in",
        dest="in",
        required=False,
        metavar="SEQUENCE_FILE",
        help="Input sequence file in FASTA format. "
        "Can be an assembled genome or transcriptome (DNA), or protein sequences from an annotated gene set. "
        "Also possible to use a path to a directory containing multiple input files.",
    )

    optional.add_argument(
        "-o",
        "--out",
        dest="out",
        required=False,
        metavar="OUTPUT",
        help="Give your analysis run a recognisable short name. "
        "Output folders and files will be labelled with this name. "
        "The path to the output folder is set with --out_path.",
    )

    optional.add_argument(
        "-m",
        "--mode",
        dest="mode",
        required=False,
        metavar="MODE",
        help="Specify which BUSCO analysis mode to run.\n"
        "There are three valid modes:\n- geno or genome, for genome assemblies (DNA)\n- tran or "
        "transcriptome, "
        "for transcriptome assemblies (DNA)\n- prot or proteins, for annotated gene sets (protein)",
    )

    optional.add_argument(
        "-l",
        "--lineage_dataset",
        dest="lineage_dataset",
        required=False,
        metavar="LINEAGE",
        help="Specify the name of the BUSCO lineage to be used.",
    )

    optional.add_argument(
        "--augustus",
        dest="use_augustus",
        action="store_true",
        required=False,
        help="Use augustus gene predictor for eukaryote runs",
    )

    optional.add_argument(
        "--augustus_parameters",
        dest="augustus_parameters",
        metavar="--PARAM1=VALUE1,--PARAM2=VALUE2",
        required=False,
        help="Pass additional arguments to Augustus. All arguments should be contained within a "
        "single string with no white space, with each argument separated by a comma.",
    )

    optional.add_argument(
        "--augustus_species",
        dest="augustus_species",
        required=False,
        help="Specify a species for Augustus training.",
    )

    optional.add_argument(
        "--auto-lineage",
        dest="auto-lineage",
        action="store_true",
        required=False,
        help="Run auto-lineage to find optimum lineage path",
    )

    optional.add_argument(
        "--auto-lineage-euk",
        dest="auto-lineage-euk",
        action="store_true",
        required=False,
        help="Run auto-placement just on eukaryote tree to find optimum lineage path",
    )

    optional.add_argument(
        "--auto-lineage-prok",
        dest="auto-lineage-prok",
        action="store_true",
        required=False,
        help="Run auto-lineage just on non-eukaryote trees to find optimum lineage path",
    )

    optional.add_argument(
        "-c",
        "--cpu",
        dest="cpu",
        type=int,
        required=False,
        metavar="N",
        help="Specify the number (N=integer) " "of threads/cores to use.",
    )

    optional.add_argument(
        "--config", dest="config_file", required=False, help="Provide a config file"
    )

    optional.add_argument(
        "--contig_break",
        dest="contig_break",
        type=int,
        required=False,
        metavar="n",
        help="Number of contiguous Ns to signify a break between contigs. Default is n=10.",
    )

    optional.add_argument(
        "--datasets_version",
        dest="datasets_version",
        required=False,
        help="Specify the version of BUSCO datasets, e.g. odb10",
    )

    optional.add_argument(
        "--download",
        dest="download",
        required=False,
        type=str,
        metavar="dataset",
        action=DirectDownload,
        help='Download dataset. Possible values are a specific dataset name, "all", "prokaryota", "eukaryota", '
        'or "virus". If used together with other command line arguments, make sure to place this last.',
    )

    optional.add_argument(
        "--download_base_url",
        dest="download_base_url",
        required=False,
        help="Set the url to the remote BUSCO dataset location",
    )

    optional.add_argument(
        "--download_path",
        dest="download_path",
        required=False,
        help="Specify local filepath for storing BUSCO dataset downloads",
    )

    optional.add_argument(
        "-e",
        "--evalue",
        dest="evalue",
        required=False,
        metavar="N",
        type=float,
        help="E-value cutoff for BLAST searches. "
        "Allowed formats, 0.001 or 1e-03 (Default: {:.0e})".format(
            BuscoConfigMain.BLAST_ARGS["evalue"]
        ),
    )

    optional.add_argument(
        "-f",
        "--force",
        action="store_true",
        required=False,
        dest="force",
        help="Force rewriting of existing files. "
        "Must be used when output files with the provided name already exist.",
    )

    optional.add_argument(
        "-h", "--help", action=CleanHelpAction, help="Show this help message and exit"
    )

    optional.add_argument(
        "--limit",
        dest="limit",
        metavar="N",
        required=False,
        type=int,
        help="How many candidate regions (contig or transcript) to consider per BUSCO (default: {})".format(
            str(BuscoConfigMain.BLAST_ARGS["limit"])
        ),
    )

    optional.add_argument(
        "--list-datasets",
        action=ListLineagesAction,
        help="Print the list of available BUSCO datasets",
    )

    optional.add_argument(
        "--long",
        action="store_true",
        required=False,
        dest="long",
        help="Optimization Augustus self-training mode (Default: Off); adds considerably to the run "
        "time, but can improve results for some non-model organisms",
    )

    optional.add_argument(
        "--metaeuk_parameters",
        dest="metaeuk_parameters",
        metavar='"--PARAM1=VALUE1,--PARAM2=VALUE2"',
        required=False,
        help="Pass additional arguments to Metaeuk for the first run. All arguments should be contained within a "
        "single string with no white space, with each argument separated by a comma.",
    )

    optional.add_argument(
        "--metaeuk_rerun_parameters",
        dest="metaeuk_rerun_parameters",
        metavar='"--PARAM1=VALUE1,--PARAM2=VALUE2"',
        required=False,
        help="Pass additional arguments to Metaeuk for the second run. All arguments should be contained within a "
        "single string with no white space, with each argument separated by a comma.",
    )

    optional.add_argument(
        "--offline",
        dest="offline",
        action="store_true",
        required=False,
        help="To indicate that BUSCO cannot attempt to download files",
    )

    optional.add_argument(
        "--out_path",
        dest="out_path",
        required=False,
        metavar="OUTPUT_PATH",
        help="Optional location for results folder, excluding results folder name. "
        "Default is current working directory.",
    )

    optional.add_argument(
        "-q",
        "--quiet",
        dest="quiet",
        required=False,
        help="Disable the info logs, displays only errors",
        action="store_true",
    )

    optional.add_argument(
        "-r",
        "--restart",
        action="store_true",
        required=False,
        dest="restart",
        help="Continue a run that had already partially completed.",
    )

    optional.add_argument(
        "--scaffold_composition",
        dest="scaffold_composition",
        action="store_true",
        required=False,
        help="Writes ACGTN content per scaffold to a file scaffold_composition.txt",
    )

    optional.add_argument(
        "--tar",
        dest="tar",
        action="store_true",
        required=False,
        help="Compress some subdirectories with many files to save space",
    )

    optional.add_argument(
        "--update-data",
        dest="update-data",
        action="store_true",
        required=False,
        help="Download and replace with last versions all lineages datasets and files necessary"
        " to their automated selection",
    )

    optional.add_argument(
        "-v",
        "--version",
        action=CleanVersionAction,
        help="Show this version and exit",
        version="BUSCO {}".format(busco.__version__),
    )

    return vars(parser.parse_args())


def main():
    """
    This function runs a BUSCO analysis according to the provided parameters.
    See the help for more details:
    ``busco -h``
    :raises SystemExit: if any errors occur
    """
    params = _parse_args()
    busco_run = BuscoMaster(params)
    busco_run.run()


# Entry point
if __name__ == "__main__":
    __spec__ = None
    main()
