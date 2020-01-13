from busco.AutoLineage import AutoSelectLineage
from busco.BuscoConfig import BuscoConfigMain
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
import os
import sys

logger = BuscoLogger.get_logger(__name__)

# todo: finalize config file

class BuscoConfigManager:

    def __init__(self, params):
        self.params = params
        self.config_file = None
        self.config = None
        self.get_config_file()
        self.runner = None

    @log("Getting config file", logger, debug=True)
    def get_config_file(self):
        """
            Check for BUSCO config file specified as a command line argument;
            if not present check if defined as an environment variable;
            if not present use default config file.
            :return config: A BuscoConfig object containing all the required configuration parameters
            """
        try:
            self.config_file = self.params["config_file"]
            if self.config_file is not None:
                return
        except KeyError:
            pass
        if os.environ.get("BUSCO_CONFIG_FILE") and os.access(os.environ.get("BUSCO_CONFIG_FILE"), os.R_OK):
            self.config_file = os.environ.get("BUSCO_CONFIG_FILE")
        else:
            raise SystemExit("Please specify a BUSCO config file using either "
                             "(i) an environment variable by entering 'export BUSCO_CONFIG_FILE=/path/to/config.ini' "
                             "or (ii) using the command line flag --config /path/to/config.ini")
        return self.config_file

    @log("Configuring BUSCO with {}", logger, attr_name="config_file")
    def load_busco_config(self, clargs):
        self.config = BuscoConfigMain(self.config_file, self.params, clargs)
        self.config.validate()
        if not self.config.check_lineage_present():
            if not self.config.getboolean("busco_run", "auto-lineage") and not self.config.getboolean("busco_run", "auto-lineage-prok"):# and not self.config.getboolean("busco_run", "auto-lineage-euk"):
                logger.warning("Running Auto Lineage Selector as no lineage dataset was specified. This will take a "
                               "little longer than normal. If you know what lineage dataset you want to use, please "
                               "specify this in the config file or using the -l (--lineage-dataset) flag in the "
                               "command line.")
            self.config.set("busco_run", "auto-lineage", "True")
            lineage_dataset_fullpath = self.auto_select_lineage()  # full path
            self.config.set("busco_run", "lineage_dataset", lineage_dataset_fullpath)
            lineage_dataset = os.path.basename(lineage_dataset_fullpath)  # base name
        else:
            if self.config.getboolean("busco_run", "auto-lineage") or self.config.getboolean("busco_run", "auto-lineage-prok"):# or self.config.getboolean("busco_run", "auto-lineage-euk"):
                logger.warning("You have selected auto-lineage but you have also provided a lineage dataset. "
                               "BUSCO will proceed with the specified dataset. "
                               "To run auto-lineage do not specify a dataset.")
            self.config.set("busco_run", "auto-lineage", "False")
            self.config.set("busco_run", "auto-lineage-prok", "False")
            self.config.set("busco_run", "auto-lineage-euk", "False")
            lineage_dataset = self.config.get("busco_run", "lineage_dataset")  # full path

        self.config.set_results_dirname(lineage_dataset)  # function always only uses basename
        self.config.download_lineage_file(lineage_dataset)  # full path will return, base name will attempt download
        # Todo: clean up error messages
        self.config.load_dataset_config()
        return

    @log("No lineage specified. Running lineage auto selector.\n", logger)
    def auto_select_lineage(self):
        asl = AutoSelectLineage(self.config)
        asl.run_auto_selector()
        asl.get_lineage_dataset()
        lineage_dataset = asl.best_match_lineage_dataset
        self.runner = asl.selected_runner
        return lineage_dataset
