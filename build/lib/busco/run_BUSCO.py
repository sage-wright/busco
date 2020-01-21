#!/usr/bin/env python
# coding: utf-8
"""
.. module:: run_BUSCO
   :synopsis: BUSCO - Benchmarking Universal Single-Copy Orthologs.
.. versionadded:: 3.0.0
.. versionchanged:: 4.0.beta1

This is the BUSCO main script.

To get help, ``busco -h``. See also the user guide.

And visit our website `<http://busco.ezlab.org/>`_

Copyright (c) 2016-2020, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""

import time
import traceback
import sys
import argparse
import os
from argparse import RawTextHelpFormatter
import busco
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.ConfigManager import BuscoConfigManager
from busco.BuscoConfig import BuscoConfigMain
from busco.Toolset import ToolException
from busco.BuscoRunner import BuscoRunner
from busco.Actions import ListLineagesAction, CleanHelpAction, CleanVersionAction

logger = BuscoLogger.get_logger(__name__)


def _parse_args():
    """
    This function parses the arguments provided by the user
    :return: a dictionary having a key for each arguments
    :rtype: dict
    """

    # todo: keyword arg order
    parser = argparse.ArgumentParser(
        description='Welcome to BUSCO %s: the Benchmarking Universal Single-Copy Ortholog assessment tool.\n'
                    'For more detailed usage information, please review the README file provided with '
                    'this distribution and the BUSCO user guide.' % busco.__version__,
        usage='busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]',
        formatter_class=RawTextHelpFormatter, add_help=False)

    optional = parser.add_argument_group('optional arguments')

    optional.add_argument(
        '-i', '--in', dest='in', required=False, metavar='FASTA FILE', help='Input sequence file in FASTA format. '
        'Can be an assembled genome or transcriptome (DNA), or protein sequences from an annotated gene set.')

    optional.add_argument(
        '-c', '--cpu', dest='cpu', required=False, metavar='N', help='Specify the number (N=integer) '
                                                                     'of threads/cores to use.')
    optional.add_argument(
        '-o', '--out', dest='out', required=False, metavar='OUTPUT',
        help='Give your analysis run a recognisable short name. '
             'Output folders and files will be labelled with this name. WARNING: do not provide a path')

    optional.add_argument(
        '-e', '--evalue', dest='evalue', required=False, metavar='N', type=float,
        help='E-value cutoff for BLAST searches. '
             'Allowed formats, 0.001 or 1e-03 (Default: %.0e)' % BuscoConfigMain.DEFAULT_ARGS_VALUES['evalue'])

    optional.add_argument(
        '-m', '--mode', dest='mode', required=False, metavar='MODE',
        help='Specify which BUSCO analysis mode to run.\n'
             'There are three valid modes:\n- geno or genome, for genome assemblies (DNA)\n- tran or '
             'transcriptome, '
             'for transcriptome assemblies (DNA)\n- prot or proteins, for annotated gene sets (protein)')

    optional.add_argument(
        '-l', '--lineage_dataset', dest='lineage_dataset', required=False, metavar='LINEAGE',
        help='Specify the name of the BUSCO lineage to be used.')

    optional.add_argument(
        '-f', '--force', action='store_true', required=False, dest='force',
        help='Force rewriting of existing files. '
             'Must be used when output files with the provided name already exist.')

    optional.add_argument(
        '--limit', dest='limit', metavar='REGION_LIMIT', required=False,
        type=int, help='How many candidate regions (contig or transcript) to consider per BUSCO (default: %s)'
                       % str(BuscoConfigMain.DEFAULT_ARGS_VALUES['limit']))

    optional.add_argument(
        '--long', action='store_true', required=False, dest='long',
        help='Optimization mode Augustus '
             'self-training (Default: Off) adds considerably to the run time, '
             'but can improve results for some non-model organisms')

    optional.add_argument(
        '-q', '--quiet', dest='quiet', required=False, help='Disable the info logs, displays only errors',
        action="store_true")

    optional.add_argument('--augustus_parameters', dest='augustus_parameters', required=False,
                          help="Pass additional arguments to Augustus. All arguments should be contained within a "
                               "single pair of quotation marks, separated by commas. E.g. \'--param1=1,--param2=2\'")

    optional.add_argument('--augustus_species', dest='augustus_species', required=False,
                          help="Specify a species for Augustus training.")

    # optional.add_argument(
    #     '-z', '--tarzip', dest='tarzip', required=False, help='Tarzip the output folders likely to '
    #                                                           'contain thousands of files',
    #     action="store_true")

    optional.add_argument(
        '--auto-lineage', dest='auto-lineage', action="store_true", required=False, help='Run auto-lineage to find optimum lineage path')

    optional.add_argument(
        '--auto-lineage-prok', dest='auto-lineage-prok', action="store_true", required=False,
        help='Run auto-lineage just on non-eukaryote trees to find optimum lineage path')

    optional.add_argument(
        '--auto-lineage-euk', dest='auto-lineage-euk', action="store_true", required=False,
        help='Run auto-placement just on eukaryote tree to find optimum lineage path')

    optional.add_argument(
        '--update-data', dest='update-data', action="store_true", required=False,
        help='Download and replace with last versions all lineages datasets and files necessary'
             ' to their automated selection')

    optional.add_argument(
        '--offline', dest='offline', action="store_true", required=False,
        help='To indicate that BUSCO cannot attempt to download files')

    optional.add_argument(
        '--config', dest='config_file', required=False, help='Provide a config file')

    optional.add_argument('-v', '--version', action=CleanVersionAction, help="Show this version and exit",
                          version='BUSCO %s' % busco.__version__)

    optional.add_argument('-h', '--help', action=CleanHelpAction, help="Show this help message and exit")

    optional.add_argument('--list-datasets', action=ListLineagesAction, help="Print the list of available BUSCO datasets")

    return vars(parser.parse_args())


def main():
    """
    This function runs a BUSCO analysis according to the provided parameters.
    See the help for more details:
    ``busco -h``
    :raises SystemExit: if any errors occur
    """

    params = _parse_args()
    run_BUSCO(params)

@log('***** Start a BUSCO v{} analysis, current time: {} *****'.format(busco.__version__, time.strftime('%m/%d/%Y %H:%M:%S')), logger)
def run_BUSCO(params):
    start_time = time.time()

    try:
        # Load a busco config file that will figure out all the params from all sources
        # i.e. provided config file, dataset cfg, and user args
        config_manager = BuscoConfigManager(params)
        config_manager.load_busco_config(sys.argv)
        config = config_manager.config

        lineage_basename = os.path.basename(config.get("busco_run", "lineage_dataset"))
        main_out_folder = config.get("busco_run", "main_out")
        lineage_results_folder = os.path.join(main_out_folder, "auto_lineage", config.get("busco_run", "lineage_results_dir"))

        if config.getboolean("busco_run", "auto-lineage"):
            if lineage_basename.startswith(("bacteria", "archaea", "eukaryota")):
                busco_run = config_manager.runner
                if lineage_basename.startswith("eukaryota") and busco_run.mode == "genome":
                    busco_run.complete_eukaryote_run()
            # It is possible that the following lineages were arrived at either by the Prodigal genetic code shortcut or by
            # BuscoPlacer. If the former, the run will have already been completed. If the latter it still needs to be done.
            elif lineage_basename.startswith(("mollicutes", "mycoplasmatales", "entomoplasmatales")) and \
                    os.path.exists(lineage_results_folder):
                busco_run = config_manager.runner
            else:
                busco_run = BuscoRunner(config)
        else:
            busco_run = BuscoRunner(config)

        if os.path.exists(lineage_results_folder):
            os.rename(lineage_results_folder, os.path.join(main_out_folder, config.get("busco_run", "lineage_results_dir")))
            busco_run.finish(time.time()-start_time, root_lineage=True)
        else:
            busco_run.run_analysis()
            BuscoRunner.final_results.append(busco_run.analysis.hmmer_results_lines)
            BuscoRunner.results_datasets.append(lineage_basename)
            busco_run.finish(time.time()-start_time)

    except ToolException as e:
        logger.error(e)
        raise SystemExit

    except SystemExit as se:
        logger.error(se)
        logger.debug(se, exc_info=True)
        logger.error('BUSCO analysis failed !')
        logger.error(
            "Check the logs, read the user guide, and check the BUSCO issue board on "
            "https://gitlab.com/ezlab/busco/issues")
        try:
            BuscoRunner.move_log_file(config)
        except NameError:
            try:
                BuscoRunner.move_log_file(config_manager.config)
            except:
                pass
        except:
            pass
        raise SystemExit

    except KeyboardInterrupt:
        logger.exception('A signal was sent to kill the process. \nBUSCO analysis failed !')
        raise SystemExit

    except BaseException:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        logger.critical("Unhandled exception occurred:\n{}\n".format("".join(traceback.format_exception(exc_type, exc_value, exc_traceback))))
        raise SystemExit


# Entry point
if __name__ == "__main__":
    main()
