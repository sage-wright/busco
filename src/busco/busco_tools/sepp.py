from busco.busco_tools.base import BaseRunner
import os
from busco.BuscoLogger import BuscoLogger
import shutil
import subprocess

logger = BuscoLogger.get_logger(__name__)


class SEPPRunner(BaseRunner):

    name = "sepp"
    cmd = "run_sepp.py"

    def __init__(self):
        super().__init__()
        self._output_folder = os.path.join(
            self.main_out, "auto_lineage", self.lineage_results_dir
        )
        self.placement_folder = os.path.join(self._output_folder, "placement_files")
        self._tmp_folder = os.path.join(self._output_folder, "sepp_tmp_files")
        self.datasets_version = self.config.get("busco_run", "datasets_version")
        self.create_dirs([self._output_folder, self.placement_folder, self._tmp_folder])

        self.tree_nwk_file = None
        self.tree_metadata_file = None
        self.supermatrix_file = None
        self.downloader = None

    def configure_runner(
        self, tree_nwk_file, tree_metadata_file, supermatrix_file, downloader
    ):
        super().configure_runner()
        self.run_number += 1
        self.tree_nwk_file = tree_nwk_file
        self.tree_metadata_file = tree_metadata_file
        self.supermatrix_file = supermatrix_file
        self.downloader = downloader

    def cleanup(self):
        shutil.rmtree(self._tmp_folder)

    def generate_job_args(self):
        yield

    def run(self):
        super().run()
        self.total = 1
        self.run_jobs()

    def configure_job(self, *args):
        sepp_job = self.create_job()
        sepp_job.add_parameter("--cpu")
        sepp_job.add_parameter(str(self.cpus))
        sepp_job.add_parameter("--outdir")
        sepp_job.add_parameter(self.placement_folder)
        sepp_job.add_parameter("-t")
        sepp_job.add_parameter(self.tree_nwk_file)
        sepp_job.add_parameter("-r")
        sepp_job.add_parameter(self.tree_metadata_file)
        sepp_job.add_parameter("-a")
        sepp_job.add_parameter(self.supermatrix_file)
        sepp_job.add_parameter("-f")
        sepp_job.add_parameter(
            os.path.join(self.placement_folder, "marker_genes.fasta")
        )
        sepp_job.add_parameter("-F")
        sepp_job.add_parameter("15")
        sepp_job.add_parameter("-m")
        sepp_job.add_parameter("amino")
        sepp_job.add_parameter("-p")
        sepp_job.add_parameter(self._tmp_folder)
        return sepp_job

    def check_tool_dependencies(self):
        pass

    def get_version(self):
        sepp_version = subprocess.check_output(
            [self.cmd, "-v"], stderr=subprocess.STDOUT, shell=False
        )
        sepp_version = sepp_version.decode("utf-8")
        sepp_version = sepp_version.strip().split(" ")[1]
        return sepp_version

    @property
    def output_folder(self):
        return self._output_folder
