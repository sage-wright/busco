import argparse
from busco.BuscoConfig import PseudoConfig
from busco.ConfigManager import BuscoConfigManager
from busco.BuscoLogger import BuscoLogger
import os
import sys

logger = BuscoLogger.get_logger(__name__)

class ListLineagesAction(argparse.Action):

    def __init__(self, option_strings, dest, nargs=0, default="==SUPPRESS==", **kwargs):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        try:
            self.config_manager = BuscoConfigManager({})
        except SystemExit as se:
            logger.error("The config file is necessary here as it contains remote and local path locations for "
                         "downloading dataset information")
            logger.error(se)
            raise SystemExit
        self.config = PseudoConfig(self.config_manager.config_file)
        try:
            self.config.load()
            self.print_lineages()
        except SystemExit as se:
            logger.error(se)
        finally:
            os.remove("busco_{}.log".format(BuscoLogger.random_id))
            parser.exit()

    def print_lineages(self):
        if self.config.update:
            lineages_list_file = self.download_lineages_list()
        else:
            lineages_list_file = self.config.existing_downloads[0]
        with open(lineages_list_file, "r") as f:
            print("".join(f.readlines()))

    def download_lineages_list(self):
        lineages_list_file = self.config.downloader.get("lineages_list.txt", "information")
        return lineages_list_file


class CleanHelpAction(argparse.Action):

    def __init__(self, option_strings, dest, nargs=0, default="==SUPPRESS==", **kwargs):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        parser.print_help()
        try:
            os.remove("busco_{}.log".format(BuscoLogger.random_id))
        except OSError:
            pass
        parser.exit()


class CleanVersionAction(argparse.Action):

    def __init__(self, option_strings, version=None, dest="==SUPPRESS==", nargs=0, default="==SUPPRESS==", **kwargs):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)
        self.version = version

    def __call__(self, parser, namespace, values, option_string=None):
        version = self.version
        if version is None:
            version = parser.version
        formatter = parser._get_formatter()
        formatter.add_text(version)
        parser._print_message(formatter.format_help(), sys.stdout)
        try:
            os.remove("busco_{}.log".format(BuscoLogger.random_id))
        except OSError:
            pass
        parser.exit()
