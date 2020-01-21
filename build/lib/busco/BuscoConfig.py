from configparser import ConfigParser
from configparser import NoOptionError
from configparser import ParsingError
from configparser import DuplicateOptionError
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from abc import ABCMeta, abstractmethod
from busco.BuscoDownloadManager import BuscoDownloadManager
import os
import shutil
import glob
import pprint

logger = BuscoLogger.get_logger(__name__)

class BaseConfig(ConfigParser):

    DEFAULT_ARGS_VALUES = {"out_path": os.getcwd(), "cpu": 1, "force": False, "evalue": 1e-3,
                           "limit": 3, "long": False, "quiet": False,
                           "download_path": os.path.join(os.getcwd(), "busco_downloads"), "datasets_version": "odb10",
                           "offline": False, "download_base_url": "https://busco-data.ezlab.org/v4/data/",
                           "auto-lineage": False, "auto-lineage-prok": False, "auto-lineage-euk": False,
                           "update-data": False}

    def __init__(self):
        super().__init__()

    def _load_config_file(self):
        """
        Load config file using ConfigParser.
        :return:
        """
        try:
            with open(self.conf_file) as cfg_file:
                self.read_file(cfg_file)
        except IOError:
            raise SystemExit("Config file {} cannot be found".format(self.conf_file))
        except ParsingError:
            raise SystemExit("Unable to parse the contents of config file {}".format(self.conf_file))
        except DuplicateOptionError:
            raise SystemExit("Duplicated entry in config file {}. Unable to load configuration.".format(self.conf_file))
        return

    def _init_downloader(self):
        """
        Initialize BuscoDownloadManager to control any file downloads from the BUSCO server.
        :return:
        """
        self.downloader = BuscoDownloadManager(self)
        return


class PseudoConfig(BaseConfig):

    def __init__(self, conf_file):
        super().__init__()
        self.conf_file = conf_file

    def load(self):
        self._load_config_file()
        self._fill_default_values()
        self._init_downloader()

    def _fill_default_values(self):
        self.set("busco_run", "offline", "False")

        try:
            self.get("busco_run", "download_base_url")
        except NoOptionError:
            self.set("busco_run", "download_base_url", type(self).DEFAULT_ARGS_VALUES["download_base_url"])

        try:
            self.get("busco_run", "download_path")
        except NoOptionError:
            self.set("busco_run", "download_path", type(self).DEFAULT_ARGS_VALUES["download_path"])

        try:
            self.update = self.getboolean("busco_run", "update-data")
            if not self.update:
                self.existing_downloads = sorted(glob.glob(os.path.join(self.get("busco_run", "download_path"), "information", "lineages_list*.txt")))[::-1]
                if self.existing_downloads:
                    logger.warning("The datasets list shown may not be up-to-date. To get current information, make sure "
                                   "you have set 'update-data=True' in your config file.")
                else:
                    raise SystemExit("Unable to download list of datasets. Please make sure you have set "
                                     "update-data=True in your config file.")

        except NoOptionError:
            self.set("busco_run", "update-data", "True")



class BuscoConfig(ConfigParser, metaclass=ABCMeta):
    """
    This class extends busco.PipeConfig to read the config.ini file. Furthermore, it uses extra args that can be
    provided through command line and information available in the dataset.cfg file to produce a single instance
    containing all correct parameters to be injected to a busco.BuscoAnalysis instance.
    """

    FORBIDDEN_HEADER_CHARS = ["ç", "¬", "¢", "´", "ê", "î", "ô", "ŵ", "ẑ", "û", "â", "ŝ", "ĝ", "ĥ", "ĵ", "ŷ",
                              "ĉ", "é", "ï", "ẅ", "ë", "ẅ", "ë", "ẗ,", "ü", "í", "ö", "ḧ", "é", "ÿ", "ẍ", "è", "é",
                              "à", "ä", "¨", "€", "£", "á"]

    FORBIDDEN_HEADER_CHARS_BEFORE_SPLIT = ["/", "\""]

    HMMER_VERSION = 3.1

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

    def download_lineage_file(self, lineage):
        """
        Download lineage dataset if not present using BuscoDownloadManager.
        :param str lineage: Basename of the lineage dataset
        :return str lineage_filepath: Local path to downloaded file
        """
        local_lineage_filepath = self.downloader.get(lineage, "lineages")
        self.set("busco_run", "lineage_dataset", local_lineage_filepath)
        return

    def load_dataset_config(self):
        """
        Load BUSCO dataset config file.
        :return:
        """
        try:
            with open(os.path.join(self.get("busco_run", "lineage_dataset"), "dataset.cfg"), "r") as target_species_file:
                dataset_kwargs = dict(line.strip().split("=") for line in target_species_file)
                for key, value in dataset_kwargs.items():
                    if key == "species":
                        try:
                            config_species = self.get("busco_run", "augustus_species")
                            if config_species != value:
                                logger.warning("An augustus species was mentioned in the config file or on the command "
                                               "line, dataset default species ({}) will be ignored".format(value))
                        except NoOptionError:
                            self.set("busco_run", "augustus_species", value)

                    elif key in ["prodigal_genetic_code", "ambiguous_cd_range_upper", "ambiguous_cd_range_lower"]:
                        self.set("prodigal", key, value)
                    else:
                        self.set("busco_run", key, value)

        except IOError:
            logger.warning("The dataset you provided does not contain the file dataset.cfg and is not valid for "
                           "BUSCO v4.0 or higher")
        return

    @abstractmethod
    def _create_required_paths(self, main_out):
        """
        Mandatory child class method. Overridden by child classes.
        :return:
        """
        if not os.path.exists(main_out):
            os.makedirs(main_out)

    def set_results_dirname(self, lineage):
        self.set("busco_run", "lineage_results_dir", "run_{}".format(os.path.basename(lineage)))
        return


class BuscoConfigAuto(BuscoConfig):

    def __init__(self, config, lineage, **kwargs):
        super().__init__(**kwargs)
        self._propagate_config(config)
        self.set_results_dirname(lineage)

        self._create_required_paths()
        self.download_lineage_file(lineage)
        self.load_dataset_config()

    def _create_required_paths(self):
        """
        Create directory for auto-lineage runs.
        :return:
        """
        main_out = os.path.join(self.get("busco_run", "main_out"), "auto_lineage")
        super()._create_required_paths(main_out)
        return

    def _propagate_config(self, config):
        """
        Copy all information from BuscoConfigMain sections into this BuscoConfigAuto object.
        Also copy BuscoDownloadManager object instead of instantiating a second one.
        :param config:
        :return:
        """
        for key, value in config.items():
            self[key] = value

        self.downloader = config.downloader
        return


class BuscoConfigMain(BuscoConfig, BaseConfig):

    MANDATORY_USER_PROVIDED_PARAMS = ["in", "out", "mode"]

    CONFIG_STRUCTURE = {"busco_run": ["in", "out", "out_path", "mode", "auto-lineage", "auto-lineage-prok",
                                  "auto-lineage-euk", "cpu", "force", "download_path", "datasets_version", "evalue",
                                  "limit", "long", "quiet", "offline", "download_base_url", "lineage_dataset",
                                  "update-data", "augustus_parameters", "augustus_species", "main_out"],
                        "tblastn": ["path", "command"],
                        "makeblastdb": ["path", "command"],
                        "prodigal": ["path", "command"],
                        "sepp": ["path", "command"],
                        "augustus": ["path", "command"],
                        "etraining": ["path", "command"],
                        "gff2gbSmallDNA.pl": ["path", "command"],
                        "new_species.pl": ["path", "command"],
                        "optimize_augustus.pl": ["path", "command"],
                        "hmmsearch": ["path", "command"]}

    def __init__(self, conf_file, params, clargs, **kwargs):
        """
        :param conf_file: a path to a config.ini file
        :type conf_file: str
        :param args: key and values matching BUSCO parameters to override config.ini values
        :type args: dict
        """
        super().__init__(**kwargs)
        self.conf_file = conf_file
        self._load_config_file()
        self.clargs = clargs
        # Update the config with args provided by the user, else keep config
        self._update_config_with_args(params)
        self._fill_default_values()



    def validate(self):
        self._check_mandatory_keys_exist()
        self._check_no_previous_run()
        self._create_required_paths()

        self._check_allowed_keys()
        self._cleanup_config()
        self._check_required_input_exists()

        self._init_downloader()
        self.persistent_tools = []

        self.log_config()

    def log_config(self):
        logger.debug("State of BUSCO config before run:")
        logger.debug(PrettyLog(vars(self)))

    def check_lineage_present(self):
        """
        Check if "lineage_dataset" parameter has been provided by user.
        :return: True if present, False if not
        :rtype: bool
        """
        try:
            lineage_dataset = self.get("busco_run", "lineage_dataset")
            datasets_version = self.get("busco_run", "datasets_version")
            if "_odb" in lineage_dataset:
                dataset_version = lineage_dataset.rsplit("_")[-1]
                if datasets_version != dataset_version:
                    logger.warning("There is a conflict in your config. You specified a dataset from {0} while "
                                   "simultaneously requesting the datasets_version parameter be {1}. Proceeding with "
                                   "the lineage dataset as specified from {0}".format(dataset_version,
                                                                                      datasets_version))
                self.set("busco_run", "datasets_version", dataset_version)
            else: # Make sure the ODB version is in the dataset name
                lineage_dataset = "_".join([lineage_dataset, datasets_version])
                self.set("busco_run", "lineage_dataset", lineage_dataset)

            datasets_version = self.get("busco_run", "datasets_version")
            if datasets_version != "odb10":
                raise SystemExit("BUSCO v4 only works with datasets from OrthoDB v10 (with the suffix '_odb10'). "
                                 "For a full list of available datasets, enter 'busco --list-datasets'. "
                                 "You can also run BUSCO using auto-lineage, to allow BUSCO to automatically select "
                                 "the best dataset for your input data.")
            return True
        except NoOptionError:
            return False

    def _check_evalue(self):
        """
        Warn the user if the config contains a non-standard e-value cutoff.
        :return:
        """
        if self.getfloat("busco_run", "evalue") != type(self).DEFAULT_ARGS_VALUES["evalue"]:
            logger.warning("You are using a custom e-value cutoff")
        return

    def _check_limit_value(self):
        """
        Check the value of limit. Ensure it is between 1 and 20, otherwise raise SystemExit.
        :return:
        """
        limit_val = self.getint("busco_run", "limit")
        if limit_val <= 0 or limit_val > 20:
            raise SystemExit("Limit must be an integer between 1 and 20 (you have used: {}). Note that this parameter "
                             "is not needed by the protein mode.".format(limit_val))
        return

    @log("Mode is {0}", logger, attr_name="_mode", on_func_exit=True, log_once=True)
    def _check_mandatory_keys_exist(self):
        """
        Make sure all mandatory user-provided parameters are present in the config.
        :return:
        """
        for param in type(self).MANDATORY_USER_PROVIDED_PARAMS:
            try:
                value = self.get("busco_run", param)
                if param == "mode":
                    synonyms = {"genome": ["genome", "geno", "genomes", "Genome", "Genomes", "Geno"],
                                "transcriptome": ["transcriptome", "tran", "transcriptomes", "trans", "Transcriptome", "Transcriptomes", "Tran", "Trans"],
                                "proteins": ["proteins", "prot", "protein", "Proteins", "Protein", "Prot", "proteome", "proteomes", "Proteome", "Proteomes"]}
                    if value not in list(synonyms["genome"] + synonyms["transcriptome"] + synonyms["proteins"]):
                        raise SystemExit("Unknown mode {}.\n'Mode' parameter must be one of "
                                         "['genome', 'transcriptome', 'proteins']".format(value))
                    if value in synonyms["genome"]:
                        self.set("busco_run", "mode", "genome")
                    elif value in synonyms["transcriptome"]:
                        self.set("busco_run", "mode", "transcriptome")
                    elif value in synonyms["proteins"]:
                        self.set("busco_run", "mode", "proteins")

                    self._mode = self.get("busco_run", "mode")

                if param == "out":
                    if value == "":
                        raise SystemExit("Please specify an output name for the BUSCO run. "
                                         "This can be done using the -o flag or in the config file")
            except NoOptionError:
                raise SystemExit("The parameter \"{} (--{})\" was not provided. "
                                 "Please add it in the config file or provide it "
                                 "through the command line".format(param, param))
        return

    def _check_allowed_keys(self):
        for section_name, options in type(self).CONFIG_STRUCTURE.items():
            for option in self.options(section_name):
                if option not in options:
                    raise SystemExit("Unrecognized option '{}' in config file".format(option))
        return

    def _check_out_value(self):
        """
        Prevent the user form using "/" in out name
        :return:
        """
        if "/" in self.get("busco_run", "out"):
            raise SystemExit("Please do not provide a full path in --out parameter, no slash. "
                             "Use out_path in the config.ini file to specify the full path.")
        return

    @log("Input file is {}", logger, attr_name="_input_filepath", on_func_exit=True, log_once=True)
    def _check_required_input_exists(self):
        """
        Test for existence of input file.
        :return:
        """
        self._input_filepath = self.get("busco_run", "in")
        if not os.path.exists(self._input_filepath):
            raise SystemExit("Input file {} does not exist".format(self._input_filepath))
        return

    def _check_no_previous_run(self):
        self.main_out = os.path.join(self.get("busco_run", "out_path"), self.get("busco_run", "out"))
        if os.path.exists(self.main_out):
            if self.getboolean("busco_run", "force"):
                self._force_remove_existing_output_dir(self.main_out)
            else:
                raise SystemExit("A run with the name {} already exists...\n"
                                 "\tIf you are sure you wish to overwrite existing files, "
                                 "please use the -f (force) option".format(self.main_out))


        return

    def _cleanup_config(self):
        """
        Collection of housekeeping functions to ensure configuration is suitable.
        :return:
        """
        self._check_out_value()
        self._check_limit_value()
        self._check_evalue()
        self._expand_all_paths()

    @staticmethod
    @log("'Force' option selected; overwriting previous results directory", logger)  # todo: review log messages
    def _force_remove_existing_output_dir(dirpath):
        """
        Remove main output folder from a previous BUSCO run.
        :return:
        """
        shutil.rmtree(dirpath)
        return

    def _create_required_paths(self):
        """
        Create main output directory and tmp directory.
        :return:
        """
        super()._create_required_paths(self.main_out)
        self.set("busco_run", "main_out", self.main_out)
        return

    def _expand_all_paths(self):
        """
        Convert relative pathnames beginning with "~" or "." into absolute paths.
        :return:
        """
        for key in self.sections():
            for item in self.items(key):
                if item[0].endswith("_path") or item[0] == "path" or item[0] == "in":
                    if item[1].startswith("~"):
                        self.set(key, item[0], os.path.expanduser(item[1]))
                    elif item[1].startswith("."):
                        self.set(key, item[0], os.path.abspath(item[1]))
        return

    def _fill_default_values(self):
        """
        Load default values into config if not provided in config file or on the command line.
        :return:
        """
        for param in list(type(self).DEFAULT_ARGS_VALUES.keys()):
            try:
                self.get("busco_run", param)
            except NoOptionError:
                self.set("busco_run", param, str(type(self).DEFAULT_ARGS_VALUES[param]))

        # Set auto-lineage to True if either auto-lineage-prok or auto-lineage-euk is selected
        if self.getboolean("busco_run", "auto-lineage-prok") or self.getboolean("busco_run", "auto-lineage-euk"):
            self.set("busco_run", "auto-lineage", "True")

        if self.getboolean("busco_run", "auto-lineage-prok") and self.getboolean("busco_run", "auto-lineage-euk"):
            logger.warning("You have specified both --auto-lineage-prok and --auto-lineage-euk. This has the same behaviour as --auto-lineage.")
            self.set("busco_run", "auto-lineage-prok", "False")
            self.set("busco_run", "auto-lineage-euk", "False")

        return

    def _update_config_with_args(self, args):
        """
        Include command line arguments in config. Overwrite any values given in the config file.
        :param args: Dictionary of parsed command line arguments. To see full list, inspect run_BUSCO.py or
        type busco -h
        :type args: dict
        :return:
        """
        for key, val in args.items():
            if key in type(self).CONFIG_STRUCTURE["busco_run"]:
                if val is not None and type(val) is not bool:
                    self.set("busco_run", key, str(val))
                elif val:  # if True
                    self.set("busco_run", key, "True")
        return


# Code taken from https://dave.dkjones.org/posts/2013/pretty-print-log-python/
class PrettyLog():
    def __init__(self, obj):
        self.obj = obj
    def __repr__(self):
        return pprint.pformat(self.obj)