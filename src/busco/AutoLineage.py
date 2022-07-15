import os
from busco.BuscoPlacer import BuscoPlacer
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.BuscoRunner import AnalysisRunner
from busco.Exceptions import BuscoError
import numpy as np

logger = BuscoLogger.get_logger(__name__)


class EmptyResultsError(BuscoError):

    def __init__(self):
        super().__init__("No genes were recognized by BUSCO. Please check the content of your input file.")


class AutoSelectLineage:
    """
    Class for selecting the best lineage dataset for the input data.
    Auto Selector works by running BUSCO using all available datasets and identifying the dataset that returns the
    highest BUSCO score.
    """

    runners = []

    @log("***** Starting Auto Select Lineage *****\n\t"
         "This process runs BUSCO on the generic lineage datasets for the domains archaea, bacteria and eukaryota. "
         "Once the optimal domain is selected, BUSCO automatically attempts to find the most appropriate BUSCO dataset "
         "to use based on phylogenetic placement.\n\t"
         "--auto-lineage-euk and --auto-lineage-prok are also available if you know your input assembly is, or is not, "
         "an eukaryote. See the user guide for more information.\n\tA reminder: Busco evaluations are valid when an "
         "appropriate dataset is used, i.e., the dataset belongs to the lineage of the species to test. "
         "Because of overlapping markers/spurious matches among domains, busco matches in another domain do not "
         "necessarily mean that your genome/proteome contains sequences from this domain. "
         "However, a high busco score in multiple domains might help you identify possible contaminations.", logger)
    def __init__(self, config_manager):
        self.config_manager = config_manager
        self.config = self.config_manager.config_main
        if self.config.getboolean("busco_run", "auto-lineage-prok"):
            self.all_lineages = ["archaea", "bacteria"]
        elif self.config.getboolean("busco_run", "auto-lineage-euk"):
            self.all_lineages = ["eukaryota"]
        else:
            self.all_lineages = ["archaea", "bacteria", "eukaryota"]
        self.dataset_version = self.config.get("busco_run", "datasets_version")
        self.callback = self.record_results
        self.s_buscos = []
        self.d_buscos = []
        self.f_buscos = []
        self.s_percents = []
        self.d_percents = []
        self.f_percents = []
        self.best_match_lineage_dataset = None
        self.current_lineage = None
        self.virus_pipeline = False
        self.selected_runner = None

    def record_results(self, s_buscos, d_buscos, f_buscos, s_percent, d_percent, f_percent):
        """
        Record results of BUSCO run.

        :param float s_buscos: Number of Single copy BUSCOs present in data
        :param float d_buscos: Number of Double copy BUSCOs present in data
        :param float f_buscos: Number of Fragmented BUSCOs present in data
        :param float s_percent: Percentage of Single copy BUSCOs present in data
        :param float d_percent: Percentage of Double copy BUSCOs present in data
        :param float f_percent: Percentage of Fragmented BUSCOs present in data
        :return: None
        """
        self.s_buscos.append(s_buscos)
        self.d_buscos.append(d_buscos)
        self.f_buscos.append(f_buscos)
        self.s_percents.append(s_percent)
        self.d_percents.append(d_percent)
        self.f_percents.append(f_percent)
        return

    @classmethod
    def reset(cls):
        cls.runners = []

    @log("Running auto selector", logger, debug=True)
    def run_auto_selector(self):
        """
        Run BUSCO on each lineage listed in all_lineages.
        :return:
        """
        root_runners = self.run_lineages_list(self.all_lineages)
        self.get_best_match_lineage(root_runners)
        if self.virus_check() and not self.config.getboolean("busco_run", "auto-lineage-euk"):
            self.virus_pipeline = True
            self.run_virus_datasets()
            self.get_best_match_lineage(type(self).runners)

        if (len(self.selected_runner.analysis.hmmer_runner.single_copy_buscos) ==
                len(self.selected_runner.analysis.hmmer_runner.multi_copy_buscos) ==
                len(self.selected_runner.analysis.hmmer_runner.fragmented_buscos) == 0):
            raise EmptyResultsError

        logger.info("{} selected\n".format(os.path.basename(self.best_match_lineage_dataset)))
        return

    def virus_check(self):
        return (self.selected_runner.analysis.hmmer_runner.s_percent < 3.0) & \
               (os.stat(self.selected_runner.analysis.input_file).st_size < 500000)

    @log("Running virus detection pipeline", logger)
    def run_virus_datasets(self):
        lineages_to_check = []
        virus_datasets = self.selected_runner.config.downloader.get("virus_datasets.txt", "information")
        with open(virus_datasets, "r") as vir_sets:
            for line in vir_sets:
                lineages_to_check.append(line.strip().split("_odb")[0])
        self.run_lineages_list(lineages_to_check)
        return

    def run_lineages_list(self, lineages_list):
        root_runners = []
        for l in lineages_list:
            self.current_lineage = "{}_{}".format(l, self.dataset_version)
            autoconfig = self.config_manager.load_busco_config_auto(self.current_lineage)
            busco_run = AnalysisRunner(autoconfig)
            busco_run.run_analysis(callback=self.callback)
            root_runners.append(busco_run)
            type(self).runners.append(busco_run)  # Save all root runs so they can be recalled if chosen
        return root_runners

    @staticmethod
    def get_max_ind(arr):
        """
        Return maximum ind(s) of array. If max value appears twice, two indices are returned.
        :param arr:
        :return:
        """
        inds = np.arange(len(arr))
        max_mask = arr == np.amax(arr)
        max_ind = inds[max_mask]
        return max_ind

    def evaluate(self, runners, use_percent=False):
        """
        Evaluate output scores from all BUSCO runs. Lineage with the highest number of complete (single + multiple)
        copy BUSCOs is assigned as the best_match_lineage.
        In case of a tie, the number of fragmented BUSCOs are used as a tiebreak.
        I this is still a tie and the number of matched BUSCOs is zero, then an error is raised.
        If there is a further nonzero tie, the tiebreak is the highest percentage of single copy BUSCOs.
        If still a tie, use the first match.
        :return
        """
        self.collate_results(runners)

        max_ind = self.get_max_ind(np.array(self.s_percents) + np.array(self.d_percents)) if use_percent \
            else self.get_max_ind(np.array(self.s_buscos) + np.array(self.d_buscos))
        if len(max_ind) > 1:
            max_ind2 = self.get_max_ind(np.array(self.f_percents)[max_ind]) if use_percent \
                else self.get_max_ind(np.array(self.f_buscos)[max_ind])
            max_ind = max_ind[max_ind2]
            if len(max_ind) > 1:
                if ((self.s_buscos[max_ind[0]] == 0.0)
                        and (self.d_buscos[max_ind[0]] == 0.0)
                        and (self.f_buscos[max_ind[0]] == 0.0)):
                    return int(0)
                else:
                    max_ind3 = self.get_max_ind(np.array(self.s_percents)[max_ind])
                    max_ind = max_ind[max_ind3]
                    if len(max_ind) > 1:
                        logger.warning("Two lineage runs scored exactly the same. Proceeding with the first.")
                        # I don't expect this error message will ever be used.
                        max_ind = max_ind[0]

        return int(max_ind)

    def collate_results(self, runners):
        self.s_buscos = [runner.analysis.hmmer_runner.single_copy for runner in runners]
        self.d_buscos = [runner.analysis.hmmer_runner.multi_copy for runner in runners]
        self.f_buscos = [runner.analysis.hmmer_runner.only_fragments for runner in runners]
        self.s_percents = [runner.analysis.hmmer_runner.s_percent for runner in runners]
        self.d_percents = [runner.analysis.hmmer_runner.d_percent for runner in runners]
        self.f_percents = [runner.analysis.hmmer_runner.f_percent for runner in runners]
        return

    def get_best_match_lineage(self, runners, use_percent=False, mollicutes_pathway=False):
        max_ind = self.evaluate(runners, use_percent)
        self.selected_runner = runners[int(max_ind)]
        self.best_match_lineage_dataset = self.selected_runner.config.get("busco_run", "lineage_dataset")
        if mollicutes_pathway:
            self.selected_runner.config.set("busco_run", "domain_run_name", os.path.basename(
                "bacteria_{}".format(self.dataset_version)))
        else:
            self.selected_runner.config.set("busco_run", "domain_run_name", os.path.basename(
                self.best_match_lineage_dataset))
        self.selected_runner.config.set("busco_run", "lineage_dataset", self.best_match_lineage_dataset)
        self.selected_runner.set_parent_dataset()
        runners.pop(int(max_ind))
        self.cleanup_disused_runs(runners)
        return

    @staticmethod
    def cleanup_disused_runs(disused_runners):
        for runner in disused_runners:
            if not runner.cleaned_up:
                runner.cleanup()

    def get_lineage_dataset(self):
        """
        Run the output of the auto selection through BuscoPlacer to obtain a more precise lineage dataset.
        :return str lineage_dataset: Local path to the optimal lineage dataset.
        """
        if self.selected_runner.domain == "eukaryota":
            self.run_busco_placer()
        elif (self.selected_runner.mode in ["proteins", "prok_tran"] and
              os.path.basename(self.selected_runner.config.get("busco_run", "lineage_dataset")).startswith("bacteria")):
            logger.info(
                "Certain mollicute clades use a different genetic code to the rest of bacteria. They are not part "
                "of the BUSCO placement tree and need to be tested separately. For more information, see the user "
                "guide.")
            use_percent = self.selected_runner.mode == "proteins"
            self.check_mollicutes(use_percent)
            if os.path.basename(self.selected_runner.config.get("busco_run", "lineage_dataset")).startswith("bacteria"):
                logger.info("Bacteria domain is a better match than the mollicutes subclade. "
                            "Continuing to tree placement.")
                self.run_busco_placer()
            else:
                logger.info("Mollicutes dataset is a better match for your data. Testing subclades...")
                self._run_3_datasets(self.selected_runner)

        elif ("geno" in self.selected_runner.mode
              and self.selected_runner.analysis.prodigal_runner.current_gc == "4"
              and os.path.basename(
                    self.selected_runner.config.get("busco_run", "lineage_dataset")).startswith("bacteria")):
            logger.info("The results from the Prodigal gene predictor indicate that your data belongs to the "
                        "mollicutes clade. Testing subclades...")
            self._run_3_datasets()
        elif self.selected_runner.domain == "viruses":
            pass
        else:
            self.run_busco_placer()
        return

    def set_best_match_lineage(self):
        AnalysisRunner.selected_dataset = os.path.basename(self.best_match_lineage_dataset)

    def check_mollicutes(self, use_percent=False):
        runners = self.run_lineages_list(["mollicutes"])
        runners.append(self.selected_runner)
        self.get_best_match_lineage(runners, use_percent=use_percent, mollicutes_pathway=True)
        return

    def run_busco_placer(self):
        if "genome" in self.selected_runner.mode:
            if self.selected_runner.domain == "prokaryota":
                protein_seqs = self.selected_runner.analysis.prodigal_runner.output_faa
            elif self.selected_runner.domain == "eukaryota":
                if self.config.getboolean("busco_run", "use_augustus"):
                    protein_seqs_dir = self.selected_runner.analysis.augustus_runner.extracted_prot_dir
                    protein_seqs = [os.path.join(protein_seqs_dir, f) for f in os.listdir(protein_seqs_dir)
                                    if f.split(".")[-2] == "faa"]
                else:
                    protein_seqs = self.selected_runner.analysis.metaeuk_runner.combined_pred_protein_seqs
        elif "tran" in self.selected_runner.mode:
            if self.selected_runner.mode == "euk_tran":
                protein_seqs = self.selected_runner.analysis.metaeuk_runner.combined_pred_protein_seqs
            elif self.selected_runner.mode == "prok_tran":
                protein_seqs = self.selected_runner.analysis.single_copy_proteins_file

        else:
            protein_seqs = self.selected_runner.config.get("busco_run", "in")
        out_path = self.config.get("busco_run", "main_out")
        run_folder = os.path.join(out_path, "auto_lineage",
                                  self.selected_runner.config.get("busco_run", "lineage_results_dir"))

        bp = BuscoPlacer(self.selected_runner.config, run_folder, protein_seqs,
                         self.selected_runner.analysis.hmmer_runner.single_copy_buscos)
        dataset_details, placement_file_versions = bp.define_dataset()
        self.config.placement_files = placement_file_versions
        # Necessary to pass these filenames to the final run to be recorded.
        lineage, supporting_markers, placed_markers = dataset_details
        self.best_match_lineage_dataset = lineage  # basename
        self.selected_runner.config.set("busco_run", "domain_run_name", os.path.basename(
            self.selected_runner.config.get("busco_run", "lineage_dataset")))
        self.selected_runner.config.set("busco_run", "lineage_dataset", self.best_match_lineage_dataset)
        self.selected_runner.set_parent_dataset()
        return

    def _run_3_datasets(self, mollicutes_runner=None):
        if mollicutes_runner:
            datasets = ["mycoplasmatales", "entomoplasmatales"]
            dataset_runners = [mollicutes_runner]
        else:
            datasets = ["mollicutes", "mycoplasmatales", "entomoplasmatales"]
            dataset_runners = []
        dataset_runners += self.run_lineages_list(datasets)
        self.get_best_match_lineage(dataset_runners, mollicutes_pathway=True)

        return
