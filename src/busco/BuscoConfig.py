from configparser import ConfigParser
from configparser import NoOptionError, NoSectionError
from configparser import ParsingError
from configparser import DuplicateOptionError
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.Exceptions import BatchFatalError
from abc import ABCMeta
from busco.BuscoDownloadManager import BuscoDownloadManager
import os
import shutil
import glob
import pprint

logger = BuscoLogger.get_logger(__name__)


class BaseConfig(ConfigParser, metaclass=ABCMeta):

    DEFAULT_ARGS_VALUES = {
        "out_path": os.getcwd(),
        "cpu": 1,
        "force": False,
        "restart": False,
        "quiet": False,
        "download_path": os.path.join(os.getcwd(), "busco_downloads"),
        "datasets_version": "odb10",
        "offline": False,
        "download_base_url": "https://busco-data.ezlab.org/v5/data/",
        "auto-lineage": False,
        "auto-lineage-prok": False,
        "auto-lineage-euk": False,
        "update-data": False,
        "use_augustus": False,
        "batch_mode": False,
        "tar": False,
    }

    AUGUSTUS_ARGS = {
        "long": False,
        "augustus_parameters": "",
        "augustus_species": "",
    }

    BLAST_ARGS = {
        "evalue": 1e-3,
        "limit": 3,
    }

    BBTOOLS_ARGS = {
        "contig_break": 10,
        "scaffold_composition": False,
    }

    METAEUK_ARGS = {
        "metaeuk_parameters": "",
        "metaeuk_rerun_parameters": "",
    }

    DEPENDENCY_SECTIONS = {
        "tblastn",
        "makeblastdb",
        "prodigal",
        "sepp",
        "metaeuk",
        "augustus",
        "etraining",
        "gff2gbSmallDNA.pl",
        "new_species.pl",
        "optimize_augustus.pl",
        "hmmsearch",
        "bbtools",
    }

    BASIC_ARGS = {
        "in",
        "out",
        "out_path",
        "mode",
        "auto-lineage",
        "auto-lineage-prok",
        "auto-lineage-euk",
        "cpu",
        "force",
        "restart",
        "download_path",
        "datasets_version",
        "quiet",
        "offline",
        "lineage_dataset",
        "update-data",
        "batch_mode",
        "tar",
        "download_base_url",
        "use_augustus",
        "download_base_url",
    }

    PERMITTED_OPTIONS = {
        "in",
        "out",
        "out_path",
        "mode",
        "auto-lineage",
        "auto-lineage-prok",
        "auto-lineage-euk",
        "cpu",
        "force",
        "restart",
        "download_path",
        "datasets_version",
        "quiet",
        "offline",
        "long",
        "augustus_parameters",
        "augustus_species",
        "download_base_url",
        "lineage_dataset",
        "update-data",
        "metaeuk_parameters",
        "metaeuk_rerun_parameters",
        "evalue",
        "limit",
        "use_augustus",
        "batch_mode",
        "tar",
        "contig_break",
        "scaffold_composition",
    }

    FORBIDDEN_HEADER_CHARS = [
        "ç",
        "¬",
        "¢",
        "´",
        "ê",
        "î",
        "ô",
        "ŵ",
        "ẑ",
        "û",
        "â",
        "ŝ",
        "ĝ",
        "ĥ",
        "ĵ",
        "ŷ",
        "ĉ",
        "é",
        "ï",
        "ẅ",
        "ë",
        "ẅ",
        "ë",
        "ẗ,",
        "ü",
        "í",
        "ö",
        "ḧ",
        "é",
        "ÿ",
        "ẍ",
        "è",
        "é",
        "à",
        "ä",
        "¨",
        "€",
        "£",
        "á",
    ]

    FORBIDDEN_HEADER_CHARS_BEFORE_SPLIT = ["/", '"']

    HMMER_VERSION = 3.1

    def __init__(self):
        super().__init__()
        config_dict = {"busco_run": type(self).DEFAULT_ARGS_VALUES}
        config_dict.update(
            {
                tool: {"path": "", "command": ""}
                for tool in type(self).DEPENDENCY_SECTIONS
            }
        )
        self.read_dict(config_dict)

    def download_lineage_file(self, lineage):
        """
        Download lineage dataset if not present using BuscoDownloadManager.
        :param str lineage: Basename of the lineage dataset
        :return str lineage_filepath: Local path to downloaded file
        """
        local_lineage_filepath = self.downloader.get(lineage, "lineages")
        self.set("busco_run", "lineage_dataset", local_lineage_filepath)
        return

    def load_dataset(self, lineage):
        self.set_results_dirname(lineage)  # function always only uses basename
        self.download_lineage_file(
            lineage
        )  # full path will return, base name will attempt download
        self.load_dataset_config()
        self.add_mode_specific_parameters()

    def update_mode(self):
        mode = self.get("busco_run", "mode")
        domain = self.get("busco_run", "domain")
        if mode == "genome":
            if domain in ["prokaryota", "viruses"]:
                self.set("busco_run", "mode", "prok_genome")
            elif domain == "eukaryota":
                if self.getboolean("busco_run", "use_augustus"):
                    self.set("busco_run", "mode", "euk_genome_aug")
                else:
                    self.set("busco_run", "mode", "euk_genome_met")
            else:
                raise BatchFatalError(
                    "Unrecognized mode {}".format(self.get("busco_run", "mode"))
                )

        elif mode == "transcriptome":
            if domain == "prokaryota":
                self.set("busco_run", "mode", "prok_tran")
            elif domain == "eukaryota":
                self.set("busco_run", "mode", "euk_tran")
            elif domain == "viruses":
                self.set(
                    "busco_run", "mode", "prok_genome"
                )  # Suggested by Mose - Prodigal may perform better on viruses than BLAST + HMMER.
            else:
                raise BatchFatalError("Unrecognized mode {}".format(mode))

    def add_mode_specific_parameters(self):
        mode = self.get("busco_run", "mode")
        if mode == "euk_genome_aug":
            self.specific_params = type(self).AUGUSTUS_ARGS
            self.specific_params.update(type(self).BLAST_ARGS)
            self.specific_params.update(type(self).BBTOOLS_ARGS)

        elif mode == "euk_genome_met":
            self.specific_params = type(self).METAEUK_ARGS
            self.specific_params.update(type(self).BBTOOLS_ARGS)
        elif mode == "prok_tran":
            self.specific_params = type(self).BLAST_ARGS
        elif mode == "euk_tran":
            self.specific_params = type(self).METAEUK_ARGS
        elif mode == "prok_genome":
            self.specific_params = type(self).BBTOOLS_ARGS
        else:
            self.specific_params = {}

        for key, val in self.specific_params.items():
            if not self.has_option("busco_run", key):
                self.set("busco_run", key, str(val))

        for option, value in self.items("busco_run"):
            if (
                option not in type(self).BASIC_ARGS
                and option not in self.specific_params
                and option in type(self).PERMITTED_OPTIONS
            ):
                logger.warning(
                    "Option {} was provided but is not used in the selected run mode, {}".format(
                        option, self.get("busco_run", "mode")
                    )
                )

    def load_dataset_config(self):
        """
        Load BUSCO dataset config file.
        :return:
        """
        try:
            with open(
                os.path.join(self.get("busco_run", "lineage_dataset"), "dataset.cfg"),
                "r",
            ) as target_species_file:
                dataset_kwargs = dict(
                    line.strip().split("=") for line in target_species_file
                )
                domain = dataset_kwargs[
                    "domain"
                ]  # Necessary to set domain kw first to enable augustus and prodigal arguments to be handled properly
                self.set("busco_run", "domain", domain)
                self.update_mode()
                del dataset_kwargs["domain"]
                for key, value in dataset_kwargs.items():
                    if key == "species":
                        try:
                            config_species = self.get("busco_run", "augustus_species")
                            if config_species != value:
                                logger.warning(
                                    "An augustus species was mentioned in the config file or on the command "
                                    "line, dataset default species ({}) will be ignored".format(
                                        value
                                    )
                                )
                        except NoOptionError:
                            if self.get("busco_run", "mode") == "euk_genome_aug":
                                self.set("busco_run", "augustus_species", value)

                    elif key in [
                        "prodigal_genetic_code",
                        "ambiguous_cd_range_upper",
                        "ambiguous_cd_range_lower",
                    ]:
                        if self.get("busco_run", "mode") in [
                            "prok_genome",
                            "prok_tran",
                        ]:
                            self.set("busco_run", key, value)
                    elif key in [
                        "max_seq_len",
                        "max_intron",
                    ]:
                        if self.get("busco_run", "mode") in [
                            "euk_genome_met",
                            "euk_tran",
                        ]:
                            self.set("busco_run", key, value)
                    else:
                        self.set("busco_run", key, value)

        except IOError:
            logger.warning(
                "The dataset you provided does not contain the file dataset.cfg and is not valid for "
                "BUSCO v4.0 or higher"
            )
        return

    def _create_required_paths(self, main_out):
        """
        :return:
        """
        if not os.path.exists(main_out):
            os.makedirs(main_out)

    def set_results_dirname(self, lineage):
        self.set(
            "busco_run",
            "lineage_results_dir",
            "run_{}".format(os.path.basename(lineage.rstrip("/"))),
        )
        return

    def _load_config_file(self, conf_file):
        """
        Load config file using ConfigParser.
        :return:
        """
        try:
            with open(conf_file) as cfg_file:
                self.read_file(cfg_file)
        except IOError:
            raise BatchFatalError("Config file {} cannot be found".format(conf_file))
        except ParsingError:
            raise BatchFatalError(
                "Unable to parse the contents of config file {}".format(conf_file)
            )
        except DuplicateOptionError:
            raise BatchFatalError(
                "Duplicated entry in config file {}. Unable to load configuration.".format(
                    conf_file
                )
            )
        return

    def _init_downloader(self):
        """
        Initialize BuscoDownloadManager to control any file downloads from the BUSCO server.
        :return:
        """
        self.downloader = BuscoDownloadManager(self)
        return

    def _update_config_with_args(self, args):
        """
        Include command line arguments in config. Overwrite any values given in the config file.
        :param args: Dictionary of parsed command line arguments. To see full list, type busco -h
        :type args: dict
        :return:
        """
        for key, val in args.items():
            if key in type(self).PERMITTED_OPTIONS:
                if val is not None and type(val) is not bool:
                    self.set("busco_run", key, str(val))
                elif val:  # if True
                    self.set("busco_run", key, "True")
        return


class PseudoConfig(BaseConfig):
    def __init__(self, conf_file, params):
        super().__init__()
        self.conf_file = conf_file
        self.params = params

    def load(self):
        if self.conf_file != "local environment":
            self._load_config_file(self.conf_file)
        else:
            self.set("busco_run", "update-data", "True")
        self._update_config_with_args(self.params)
        self._fill_default_values()
        self._init_downloader()

    def _fill_default_values(self):
        self.set("busco_run", "offline", "False")

        try:
            self.get("busco_run", "download_base_url")
        except NoOptionError:
            self.set(
                "busco_run",
                "download_base_url",
                type(self).DEFAULT_ARGS_VALUES["download_base_url"],
            )

        try:
            self.get("busco_run", "download_path")
        except NoOptionError:
            self.set(
                "busco_run",
                "download_path",
                type(self).DEFAULT_ARGS_VALUES["download_path"],
            )

        try:
            self.update = self.getboolean("busco_run", "update-data")
            if not self.update:
                self.existing_downloads = sorted(
                    glob.glob(
                        os.path.join(
                            self.get("busco_run", "download_path"),
                            "information",
                            "lineages_list*.txt",
                        )
                    )
                )[::-1]
                if self.existing_downloads:
                    logger.warning(
                        "The datasets list shown may not be up-to-date. To get current information, make sure "
                        "you have set 'update-data=True' in your config file."
                    )
                else:
                    raise BatchFatalError(
                        "Unable to download list of datasets. Please use the --update-data option to use the latest "
                        "available datasets."
                    )

        except NoOptionError:
            self.set("busco_run", "update-data", "True")


class BuscoConfigAuto(BaseConfig):
    def __init__(self, config, lineage):
        super().__init__()
        self._propagate_config(config)
        self.load_dataset(lineage)
        self._create_required_paths()

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


class BuscoConfigMain(BaseConfig):

    MANDATORY_USER_PROVIDED_PARAMS = ["in", "out", "mode"]

    def __init__(self, conf_file, params):
        """
        :param conf_file: a path to a config.ini file
        :type conf_file: str
        :param params: key and values matching BUSCO parameters to override config.ini values
        :type params: dict
        """
        super().__init__()
        self.conf_file = conf_file
        self.params = params
        self.main_out = None
        self._input_filepath = None

    def configure(self):
        if self.conf_file != "local environment":
            self._load_config_file(self.conf_file)
        # Update the config with args provided by the user, else keep config
        self._update_config_with_args(self.params)
        self._update_logger()
        self.harmonize_auto_lineage_settings()

    def _update_logger(self):
        if self.getboolean("busco_run", "quiet"):
            type(logger).quiet = True

    @classmethod
    def merge_two_dicts(cls, x, y):
        """Given two dictionaries, merge them into a new dict as a shallow copy."""
        z = x.copy()
        z.update(y)
        return z

    def validate(self):
        self._check_mandatory_keys_exist()
        self._cleanup_config()
        self._check_no_previous_run()
        self._check_allowed_keys()
        self._create_required_paths(self.main_out)
        self._check_required_input_exists()
        self._check_batch_mode()

        self._init_downloader()

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
                dataset_version = lineage_dataset.rsplit("_")[-1].rstrip("/")
                if datasets_version != dataset_version:
                    logger.warning(
                        "There is a conflict in your config. You specified a dataset from {0} while "
                        "simultaneously requesting the datasets_version parameter be {1}. Proceeding with "
                        "the lineage dataset as specified from {0}".format(
                            dataset_version, datasets_version
                        )
                    )
                self.set("busco_run", "datasets_version", dataset_version)
            else:  # Make sure the ODB version is in the dataset name
                lineage_dataset = "_".join([lineage_dataset, datasets_version])
                self.set("busco_run", "lineage_dataset", lineage_dataset)

            datasets_version = self.get("busco_run", "datasets_version")
            if datasets_version != "odb10":
                raise BatchFatalError(
                    "BUSCO v5 only works with datasets from OrthoDB v10 (with the suffix '_odb10'). "
                    "For a full list of available datasets, enter 'busco --list-datasets'. "
                    "You can also run BUSCO using --auto-lineage, to allow BUSCO to automatically select "
                    "the best dataset for your input data."
                )
            return True
        except NoOptionError:
            return False

    def _check_batch_mode(self):
        if os.path.isdir(self._input_filepath):
            self.set("busco_run", "batch_mode", "True")
        elif not os.path.isfile(self._input_filepath):
            raise BatchFatalError(
                "Unrecognized input type. Please use either a single file or a directory name (for batch mode)"
            )
        return

    def _check_evalue(self):
        """
        Warn the user if the config contains a non-standard e-value cutoff.
        :return:
        """
        if self.getfloat("busco_run", "evalue") != type(self).BLAST_ARGS["evalue"]:
            logger.warning("You are using a custom e-value cutoff")
        return

    def _check_limit_value(self):
        """
        Check the value of limit. Ensure it is between 1 and 20, otherwise raise BatchFatalError.
        :return:
        """
        limit_val = self.getint("busco_run", "limit")
        if limit_val <= 0 or limit_val > 20:
            raise BatchFatalError(
                "Limit must be an integer between 1 and 20 (you have used: {}). Note that this parameter "
                "is not needed by the protein mode.".format(limit_val)
            )
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
                    synonyms = {
                        "genome": [
                            "genome",
                            "geno",
                            "genomes",
                            "Genome",
                            "Genomes",
                            "Geno",
                        ],
                        "transcriptome": [
                            "transcriptome",
                            "tran",
                            "transcriptomes",
                            "trans",
                            "Transcriptome",
                            "Transcriptomes",
                            "Tran",
                            "Trans",
                        ],
                        "proteins": [
                            "proteins",
                            "prot",
                            "protein",
                            "Proteins",
                            "Protein",
                            "Prot",
                            "proteome",
                            "proteomes",
                            "Proteome",
                            "Proteomes",
                        ],
                    }
                    if value not in list(
                        synonyms["genome"]
                        + synonyms["transcriptome"]
                        + synonyms["proteins"]
                    ):
                        raise BatchFatalError(
                            "Unknown mode {}.\n'Mode' parameter must be one of "
                            "['genome', 'transcriptome', 'proteins']".format(value)
                        )
                    if value in synonyms["genome"]:
                        self.set("busco_run", "mode", "genome")
                    elif value in synonyms["transcriptome"]:
                        self.set("busco_run", "mode", "transcriptome")
                    elif value in synonyms["proteins"]:
                        self.set("busco_run", "mode", "proteins")

                    self._mode = self.get("busco_run", "mode")

                if param == "out":
                    if value == "":
                        raise BatchFatalError(
                            "Please specify an output name for the BUSCO run. "
                            "This can be done using the -o flag or in the config file"
                        )
            except NoOptionError:
                raise BatchFatalError(
                    'The parameter "{} (--{})" was not provided. '
                    "Please add it in the config file or provide it "
                    "through the command line".format(param, param)
                )
        return

    def _check_allowed_keys(self):
        full_dict = {"busco_run": type(self).PERMITTED_OPTIONS}
        full_dict.update(
            {
                dependency: ["path", "command"]
                for dependency in type(self).DEPENDENCY_SECTIONS
            }
        )

        for section_name, options in full_dict.items():
            try:
                for option in self.options(section_name):
                    if option not in options:
                        raise BatchFatalError(
                            "Unrecognized option '{}' in config file".format(option)
                        )
            except NoSectionError:
                logger.warning("Section {} not found".format(section_name))
        return

    def _check_out_value(self):
        """
        Previously prevented the user form using "/" in out name. Now just strips any trailing "/".
        :return:
        """
        out = self.get("busco_run", "out")
        if "/" in out:
            self.set("busco_run", "out", out.strip("/"))
        return

    def _check_required_input_exists(self):
        """
        Test for existence of input file.
        :return:
        """
        self._input_filepath = self.get("busco_run", "in")
        if not os.path.exists(self._input_filepath):
            raise BatchFatalError(
                "Input {} does not exist".format(self._input_filepath)
            )
        return

    def _check_no_previous_run(self):
        self.main_out = os.path.join(
            self.get("busco_run", "out_path"), self.get("busco_run", "out")
        )
        if os.path.exists(self.main_out):
            if self.getboolean("busco_run", "force"):
                self._force_remove_existing_output_dir(self.main_out)
            elif self.getboolean("busco_run", "restart"):
                logger.info(
                    "Attempting to restart the run using the following directory: {}".format(
                        self.main_out
                    )
                )
            else:
                raise BatchFatalError(
                    "A run with the name {} already exists...\n"
                    "\tIf you are sure you wish to overwrite existing files, "
                    "please use the -f (force) option".format(self.main_out)
                )
        elif self.getboolean("busco_run", "restart"):
            logger.warning(
                "Restart mode not available as directory {} does not exist.".format(
                    self.main_out
                )
            )
            self.set("busco_run", "restart", "False")

        return

    def _cleanup_config(self):
        """
        Collection of housekeeping functions to ensure configuration is suitable.
        :return:
        """
        self._check_out_value()
        self._expand_all_paths()
        self._harmonize_augustus_options()

    def _harmonize_augustus_options(self):
        augustus_selected = self.getboolean("busco_run", "use_augustus")
        if not augustus_selected:
            try:
                augustus_species = self.get("busco_run", "augustus_species")
            except NoOptionError:
                augustus_species = None

            try:
                augustus_params = self.get("busco_run", "augustus_parameters")
            except NoOptionError:
                augustus_params = None

            if (augustus_species or augustus_params) and not augustus_selected:
                self.set("busco_run", "use_augustus", "True")
        return

    def _update_augustus_options(self):

        if self.getboolean("busco_run", "use_augustus"):
            for key, val in type(self).AUGUSTUS_ARGS.items():
                try:
                    self.get("busco_run", key)
                except NoOptionError:
                    self.set("busco_run", key, val)
            self._check_limit_value()
            self._check_evalue()

    @staticmethod
    @log("'Force' option selected; overwriting previous results directory", logger)
    def _force_remove_existing_output_dir(dirpath):
        """
        Remove main output folder from a previous BUSCO run.
        :return:
        """
        shutil.rmtree(dirpath)
        return

    def _create_required_paths(self, main_out):
        """
        Create main output directory and tmp directory.
        :return:
        """
        super()._create_required_paths(main_out)
        self.set("busco_run", "main_out", main_out)
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
                    elif item[1]:
                        self.set(key, item[0], os.path.abspath(item[1]))

        return

    def harmonize_auto_lineage_settings(self):

        if not self.check_lineage_present():
            if (
                not self.getboolean("busco_run", "auto-lineage")
                and not self.getboolean("busco_run", "auto-lineage-prok")
                and not self.getboolean("busco_run", "auto-lineage-euk")
            ):
                logger.warning(
                    "Running Auto Lineage Selector as no lineage dataset was specified. This will take a "
                    "little longer than normal. If you know what lineage dataset you want to use, please "
                    "specify this in the config file or using the -l (--lineage-dataset) flag in the "
                    "command line."
                )

            elif self.getboolean("busco_run", "auto-lineage-prok") and self.getboolean(
                "busco_run", "auto-lineage-euk"
            ):
                logger.warning(
                    "You have specified both --auto-lineage-prok and --auto-lineage-euk. This has the same behaviour "
                    "as --auto-lineage."
                )
                self.set("busco_run", "auto-lineage-prok", "False")
                self.set("busco_run", "auto-lineage-euk", "False")

            self.set("busco_run", "auto-lineage", "True")

        else:
            if (
                self.getboolean("busco_run", "auto-lineage")
                or self.getboolean("busco_run", "auto-lineage-prok")
                or self.getboolean("busco_run", "auto-lineage-euk")
            ):
                logger.warning(
                    "You have selected auto-lineage but you have also provided a lineage dataset. "
                    "BUSCO will proceed with the specified dataset. "
                    "To run auto-lineage do not specify a dataset."
                )
            self.set("busco_run", "auto-lineage", "False")
            self.set("busco_run", "auto-lineage-prok", "False")
            self.set("busco_run", "auto-lineage-euk", "False")
        return


# Code taken from https://dave.dkjones.org/posts/2013/pretty-print-log-python/
class PrettyLog:
    def __init__(self, obj):
        self.obj = obj

    def __repr__(self):
        return pprint.pformat(self.obj)
