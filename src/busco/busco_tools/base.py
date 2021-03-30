import os
from busco.BuscoLogger import BuscoLogger
from busco.busco_tools.Toolset import Tool
from shutil import which
from abc import ABCMeta, abstractmethod
from busco.BuscoConfig import BuscoConfigAuto
from busco.Exceptions import BuscoError
import time

logger = BuscoLogger.get_logger(__name__)


class ToolException(Exception):
    """
    Module-specific exception
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


class BaseRunner(Tool, metaclass=ABCMeta):

    config = None
    tool_versions = {}

    def __init__(self):
        super().__init__()
        self.run_number = 0
        self.input_file = self.config.get("busco_run", "in")
        self.main_out = self.config.get("busco_run", "main_out")
        self.working_dir = (
            os.path.join(self.main_out, "auto_lineage")
            if isinstance(self.config, BuscoConfigAuto)
            else self.main_out
        )
        self.lineage_results_dir = self.config.get("busco_run", "lineage_results_dir")
        self.run_folder = os.path.join(self.working_dir, self.lineage_results_dir)
        self.log_folder = os.path.join(self.main_out, "logs")
        self.cpus = self.config.getint("busco_run", "cpu")
        self.lineage_dataset = self.config.get("busco_run", "lineage_dataset")
        self.domain = self.config.get("busco_run", "domain")

        if not self.check_tool_available():
            raise ToolException(
                "{} tool cannot be found. Please check the 'path' and 'command' parameters "
                "provided in the config file. Do not include the command in the "
                "path!".format(self.name)
            )
        self.version = self.get_version()
        type(self).tool_versions[self.name] = self.version
        self.check_tool_dependencies()

        self.checkpoint_file = None

        self.logfile_path_out = os.path.join(
            self.config.get("busco_run", "main_out"),
            "logs",
            "{}_out.log".format(self.name),
        )
        self.logfile_path_err = (
            self.logfile_path_out.rpartition("_out.log")[0] + "_err.log"
        )

    def init_checkpoint_file(self):
        self.checkpoint_file = os.path.join(self.output_folder, ".checkpoint")

    def write_checkpoint_file(self):
        with open(self.checkpoint_file, "a") as cpt_file:
            cpt_file.write("Tool: {}\n".format(self.name))
            cpt_file.write("Version: {}\n".format(self.version))
            cpt_file.write("Run: {}\n".format(self.run_number))
            cpt_file.write("Time: {}\n".format(time.strftime("%m/%d/%Y %H:%M:%S")))
            cpt_file.write("Completed {} jobs\n\n".format(self.total))

    def check_previous_completed_run(self):
        if not os.path.exists(self.checkpoint_file):
            return False
        else:
            with open(self.checkpoint_file, "r") as cpt_file:
                lines = cpt_file.readlines()
                tool_names = [s.strip().split(": ")[1] for s in lines[0::6]]
                tool_versions = [s.strip().split(": ")[1] for s in lines[1::6]]
                tool_run_numbers = [s.strip().split(": ")[1] for s in lines[2::6]]
                try:
                    start_search = 0
                    while True:
                        tool_ind = tool_names.index(self.name, start_search)
                        if str(self.version) != str(tool_versions[tool_ind]):
                            logger.warning(
                                "A previous run used {} version {}. "
                                "The restarted run is using {} version "
                                "{}".format(
                                    self.name,
                                    tool_versions[tool_ind],
                                    self.name,
                                    self.version,
                                )
                            )
                        if int(tool_run_numbers[tool_ind]) == int(self.run_number):
                            return True
                        elif int(tool_run_numbers[tool_ind]) < int(self.run_number):
                            start_search = tool_ind + 1
                        else:
                            raise BuscoError(
                                "Something went wrong. Information for {} run {} missing but "
                                "information for run {} found.".format(
                                    self.name,
                                    self.run_number,
                                    tool_run_numbers[tool_ind],
                                )
                            )

                except ValueError:
                    return False

                except TypeError:
                    logger.warning(
                        "Unable to parse {} file. Restart mode not available.".format(
                            self.checkpoint_file
                        )
                    )

    @abstractmethod
    def check_tool_dependencies(self):
        pass

    @abstractmethod
    def configure_job(self, *args):
        pass

    @abstractmethod
    def generate_job_args(self):
        pass

    @property
    @abstractmethod
    def output_folder(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def name(self):
        raise NotImplementedError

    @abstractmethod
    def run(self):
        if self.version is not None:
            logger.debug("Tool: {}".format(self.name))
            logger.debug("Version: {}".format(self.version))

    @staticmethod
    def create_dirs(dirnames):
        """
        Create all required directories

        :param dirnames: list of paths already constructed
        :return:
        """
        if isinstance(dirnames, str):
            os.makedirs(dirnames, exist_ok=True)
        elif isinstance(dirnames, list):
            for d in dirnames:
                os.makedirs(d, exist_ok=True)
        else:
            raise TypeError("'dirnames' should be either a str or a list")

    def check_tool_available(self):
        """
        Check tool's availability.


        :return: True if the tool can be run, False if it is not the case
        :rtype: bool
        """
        try:
            self.get_tool_from_config()
        except ToolException:
            self.get_tool_from_environment()

        return which(self.cmd) is not None  # True if tool available

    def get_tool_from_environment(self):
        which_tool = which(self.cmd)
        if not which_tool:
            raise ToolException("Tool {} not found".format(self.name))

    def get_tool_from_config(self):
        """
        1. The section ['name'] is available in the config
        2. This section contains keys 'path' and 'command'
        3. The string resulted from concatenation of values of these two keys
        represents the full path to the command
        :return:
        """
        if not self.config.has_section(self.name):
            raise ToolException(
                "Section for the tool [{}] is not present in the config file".format(
                    self.name
                )
            )

        if not self.config.has_option(self.name, "path") or not self.config.get(
            self.name, "path"
        ):
            raise ToolException(
                "Key 'path' in the section [{}] is not present in the config file".format(
                    self.name
                )
            )

        if self.config.has_option(self.name, "command"):
            executable = self.config.get(self.name, "command")
        else:
            executable = self.name

        self.cmd = os.path.join(self.config.get(self.name, "path"), executable)

        return

    @abstractmethod
    def get_version(self):
        return


class NoGenesError(Exception):
    def __init__(self, gene_predictor):
        self.gene_predictor = gene_predictor


class NoRerunFile(Exception):
    def __init__(self):
        pass
