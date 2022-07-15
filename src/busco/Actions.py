import argparse
from busco.BuscoConfig import PseudoConfig
from busco.ConfigManager import BuscoConfigManager
from busco.BuscoLogger import BuscoLogger
from busco.BuscoDownloadManager import BuscoDownloadManager
import os
import sys


class ListLineagesAction(argparse.Action):

    logger = BuscoLogger.get_logger(__name__)

    def __init__(self, option_strings, dest, nargs=0, default="==SUPPRESS==", **kwargs):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        try:
            self.config_manager = BuscoConfigManager(namespace.__dict__)
        except SystemExit as se:
            type(self).logger.error(se)
            raise SystemExit(1)

        self.config = PseudoConfig(
            self.config_manager.config_file, self.config_manager.params
        )
        try:
            self.config.load()
            self.print_lineages()
        except SystemExit as se:
            type(self).logger.error(se)
            raise SystemExit(1)
        finally:
            try:
                os.remove("busco_{}.log".format(BuscoLogger.pid))
            except FileNotFoundError:
                pass
            parser.exit()

    def print_lineages(self):
        if self.config.update:
            lineages_list_file = self.download_lineages_list()
        else:
            lineages_list_file = self.config.existing_downloads[0]
        with open(lineages_list_file, "r") as f:
            print("".join(f.readlines()))

    def download_lineages_list(self):
        lineages_list_file = self.config.downloader.get(
            "lineages_list.txt", "information"
        )
        return lineages_list_file


class DirectDownload(argparse.Action):
    logger = BuscoLogger.get_logger(__name__)

    def __init__(
        self, option_strings, dest, nargs="*", default="==SUPPRESS==", **kwargs
    ):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        try:
            self.config_manager = BuscoConfigManager(namespace.__dict__)
        except SystemExit as se:
            type(self).logger.error(se)
            raise SystemExit(1)

        self.config = PseudoConfig(
            self.config_manager.config_file, self.config_manager.params
        )
        try:
            self.config.load()
            bdm = BuscoDownloadManager(self.config)
            bdm.update_data = True
            self.download_datasets(bdm, values)

        except SystemExit as se:
            type(self).logger.error(se)
            raise SystemExit(1)
        finally:
            try:
                os.remove("busco_{}.log".format(BuscoLogger.pid))
            except FileNotFoundError:
                pass
            parser.exit()

    def download_datasets(self, bdm, values):
        for item in values:
            self.get(bdm, item)

    def get(self, bdm, item):
        if item == "all":
            files_to_get = bdm.version_files.index
        elif item == "prokaryota":
            files_to_get = [
                item
                for item in bdm.version_files.index
                if bdm.version_files.loc[item]["domain"] == "Prokaryota"
            ]
        elif item == "eukaryota":
            files_to_get = [
                item
                for item in bdm.version_files.index
                if bdm.version_files.loc[item]["domain"] == "Eukaryota"
            ]
        elif item == "virus":
            files_to_get = [
                item
                for item in bdm.version_files.index
                if bdm.version_files.loc[item]["domain"] == "Virus"
            ]
        else:
            try:
                if isinstance(item, str) and item in bdm.version_files.index:
                    files_to_get = [item]
                else:
                    raise KeyError
            except KeyError:
                type(self).logger.error("{} is not a recognized option".format(item))
                files_to_get = []

        filetypes = [bdm.version_files.loc[f]["type"] for f in files_to_get]
        for f, filename in enumerate(files_to_get):
            bdm.get(filename, filetypes[f])


class CleanHelpAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=0, default="==SUPPRESS==", **kwargs):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        parser.print_help()
        try:
            os.remove("busco_{}.log".format(BuscoLogger.pid))
        except OSError:
            pass
        parser.exit()


class CleanVersionAction(argparse.Action):
    def __init__(
        self,
        option_strings,
        version=None,
        dest="==SUPPRESS==",
        nargs=0,
        default="==SUPPRESS==",
        **kwargs
    ):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)
        self.version = version

    def __call__(self, parser, namespace, values, option_string=None):
        version = self.version
        if version is None:
            version = "unknown"
        formatter = parser._get_formatter()
        formatter.add_text(version)
        parser._print_message(formatter.format_help(), sys.stdout)
        try:
            os.remove("busco_{}.log".format(BuscoLogger.pid))
        except OSError:
            pass
        parser.exit()
