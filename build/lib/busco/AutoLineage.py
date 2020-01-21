import os
from busco.BuscoConfig import BuscoConfigAuto
from busco.BuscoPlacer import BuscoPlacer
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.BuscoRunner import BuscoRunner
import numpy as np
import numpy.ma as ma

logger = BuscoLogger.get_logger(__name__)

class AutoSelectLineage:
    """
    Class for selecting the best lineage dataset for the input data.
    Auto Selector works by running BUSCO using all available datasets and identifying the dataset that returns the
    highest BUSCO score.
    """

    @log("***** Starting Auto Select Lineage *****\n\t"
         "This process runs BUSCO on the generic lineage datasets for the domains archaea, bacteria and eukaryota. "
         "Once the optimal domain is selected, BUSCO automatically attempts to find the most appropriate BUSCO dataset "
         "to use based on phylogenetic placement.\n\t"
         "--auto-lineage-euk and --auto-lineage-prok are also available if you know your input assembly is, or is not, an"
         " eukaryote. See the user guide for more information.\n\tA reminder: Busco evaluations are valid when an "
         "appropriate dataset is used, i.e., the dataset belongs to the lineage of the species to test. "
         "Because of overlapping markers/spurious matches among domains, busco matches in another domain do not "
         "necessarily mean that your genome/proteome contains sequences from this domain. "
         "However, a high busco score in multiple domains might help you identify possible contaminations.", logger)
    def __init__(self, config):
        self.config = config
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

    def record_results(self, s_buscos, d_buscos, f_buscos, s_percent, d_percent, f_percent):
        """
        Record results of BUSCO run.

        :param float s_buscos: Number of Single copy BUSCOs present in data
        :param float d_buscos: Number of Double copy BUSCOs present in data
        :param float f_buscos: Number of Fragmented BUSCOs present in data
        :return: None
        """
        self.s_buscos.append(s_buscos)
        self.d_buscos.append(d_buscos)
        self.f_buscos.append(f_buscos)
        self.s_percents.append(s_percent)
        self.d_percents.append(d_percent)
        self.f_percents.append(f_percent)
        return

    @log("Running auto selector", logger, debug=True)
    def run_auto_selector(self):
        """
        Run BUSCO on each lineage listed in all_lineages.
        :return:
        """
        root_runners = self.run_lineages_list(self.all_lineages)
        self.get_best_match_lineage(root_runners)
        self.config.set("busco_run", "domain_run_name", os.path.basename(self.best_match_lineage_dataset))
        BuscoRunner.final_results.append(self.selected_runner.analysis.hmmer_results_lines)
        BuscoRunner.results_datasets.append(os.path.basename(self.best_match_lineage_dataset))
        return

    def run_lineages_list(self, lineages_list):
        root_runners = []
        for l in lineages_list:
            self.current_lineage = "{}_{}".format(l, self.dataset_version)
            autoconfig = BuscoConfigAuto(self.config, self.current_lineage)
            # The following line creates a direct reference, so whenever one analysis run adds a tool to this list it
            # is automatically updated here too.
            autoconfig.persistent_tools = self.config.persistent_tools
            busco_run = BuscoRunner(autoconfig)
            busco_run.run_analysis(callback=self.callback)
            root_runners.append(busco_run)  # Save all root runs so they can be recalled if chosen
        return root_runners


    def evaluate(self):
        """
        Evaluate output scores from all BUSCO runs. Lineage with the highest number of complete (single + multiple)
        copy BUSCOs is assigned as the best_match_lineage.
        In case of a tie, the number of fragmented BUSCOs are used as a tiebreak.
        I this is still a tie and the number of matched BUSCOs is zero, then an error is raised.
        If there is a further nonzero tie, the tiebreak is the highest percentage of single copy BUSCOs.
        If still a tie, use the first match.
        :return
        """
        total_complete = np.array(self.s_buscos) + np.array(self.d_buscos)
        inds = np.arange(len(total_complete))

        max_mask = total_complete == np.amax(total_complete)
        max_ind = inds[max_mask]
        if len(max_ind) > 1:
            self.f_buscos = np.array(self.f_buscos)
            tie_break = ma.array(self.f_buscos, mask=~max_mask)
            max_mask &= tie_break == ma.max(tie_break)
            max_ind = inds[max_mask]
            if len(max_ind) > 1:
                if ((self.s_buscos[max_ind[0]] == 0.0)
                        and (self.d_buscos[max_ind[0]] == 0.0)
                        and (self.f_buscos[max_ind[0]] == 0.0)):
                    raise SystemExit("No genes were recognized by BUSCO. Please check the content of your input file.")
                else:
                    self.s_percents = np.array(self.s_percents)
                    tie_break = ma.array(self.s_percents, mask=~max_mask)
                    max_mask &= tie_break == ma.max(tie_break)
                    max_ind = inds[max_mask]
                    if len(max_ind) > 1:
                        logger.warning("Two lineage runs scored exactly the same. Proceeding with the first.")  # I don't expect this error message will ever be used.
                        max_ind = max_ind[0]

        return max_ind

    @log("{} selected\n", logger, attr_name="best_match_lineage_dataset", apply="basename", on_func_exit=True)
    def get_best_match_lineage(self, runners):
        max_ind = self.evaluate()
        self.selected_runner = runners[int(max_ind)]
        self.best_match_lineage_dataset = self.selected_runner.config.get("busco_run", "lineage_dataset")
        runners.pop(int(max_ind))
        self.cleanup_disused_runs(runners)
        return

    def cleanup_disused_runs(self, disused_runners):
        for runner in disused_runners:
            runner.analysis._cleanup()


    def get_lineage_dataset(self):  # todo: rethink structure after BuscoPlacer is finalized and protein mode with mollicutes is fixed.
        """
        Run the output of the auto selection through BuscoPlacer to obtain a more precise lineage dataset.
        :return str lineage_dataset: Local path to the optimal lineage dataset.
        """
        if self.selected_runner.domain == "eukaryota":
            self.run_busco_placer()
        elif (self.selected_runner.mode in ["proteins", "prot", "transcriptome", "tran"] and
              os.path.basename(self.selected_runner.config.get("busco_run", "lineage_dataset")).startswith("bacteria")):
            logger.info(
                "Certain mollicute clades use a different genetic code to the rest of bacteria. They are not part "
                "of the BUSCO placement tree and need to be tested separately. For more information, see the user "
                "guide.")
            self.check_mollicutes()
            if os.path.basename(self.selected_runner.config.get("busco_run", "lineage_dataset")).startswith("bacteria"):
                logger.info("Bacteria domain is a better match than the mollicutes subclade. Continuing to tree placement.")
                self.run_busco_placer()
            else:
                logger.info("Mollicutes dataset is a better match for your data. Testing subclades...")
                self._run_3_datasets(self.selected_runner)
                BuscoRunner.final_results.append(self.selected_runner.analysis.hmmer_results_lines)
                BuscoRunner.results_datasets.append(os.path.basename(self.best_match_lineage_dataset))
        elif ("geno" in self.selected_runner.mode and self.selected_runner.analysis.code_4_selected and
              os.path.basename(self.selected_runner.config.get("busco_run", "lineage_dataset")).startswith("bacteria")):
            logger.info("The results from the Prodigal gene predictor indicate that your data belongs to the "
                        "mollicutes clade. Testing subclades...")
            self._run_3_datasets()
            BuscoRunner.final_results.append(self.selected_runner.analysis.hmmer_results_lines)
            BuscoRunner.results_datasets.append(os.path.basename(self.best_match_lineage_dataset))
        else:
            self.run_busco_placer()
        return

    def check_mollicutes(self):
        self.s_buscos = []
        self.d_buscos = []
        self.f_buscos = []
        self.s_percents = []
        self.d_percents = []
        self.f_percents = []
        runners = self.run_lineages_list(["mollicutes"])
        runners.append(self.selected_runner)
        self.s_buscos.append(self.selected_runner.analysis.single_copy)
        self.d_buscos.append(self.selected_runner.analysis.multi_copy)
        self.f_buscos.append(self.selected_runner.analysis.only_fragments)
        self.s_percents.append(self.selected_runner.analysis.s_percent)
        self.d_percents.append(self.selected_runner.analysis.d_percent)
        self.f_percents.append(self.selected_runner.analysis.f_percent)
        self.get_best_match_lineage(runners)
        return

    def run_busco_placer(self):  # todo: revisit structure of this method after cleaning BuscoPlacer
        if "genome" in self.selected_runner.mode:
            if self.selected_runner.domain == "prokaryota":
                protein_seqs = self.selected_runner.analysis.prodigal_runner.output_faa
            elif self.selected_runner.domain == "eukaryota":
                protein_seqs = self.selected_runner.analysis.augustus_runner.output_sequences
        else:
            protein_seqs = self.selected_runner.config.get("busco_run", "in")
        out_path = self.config.get("busco_run", "main_out")
        run_folder = os.path.join(out_path, "auto_lineage", self.selected_runner.config.get("busco_run", "lineage_results_dir"))
        bp = BuscoPlacer(self.selected_runner.config, run_folder, protein_seqs, self.selected_runner.analysis.hmmer_runner.single_copy_buscos)
        dataset_details, placement_file_versions = bp.define_dataset()
        lineage, supporting_markers, placed_markers = dataset_details
        lineage = "{}_{}".format(lineage, self.config.get("busco_run", "datasets_version"))  # todo: this should probably be done in buscoplacer
        self.best_match_lineage_dataset = os.path.join(self.config.get("busco_run", "download_path"),
                                                       "lineages",
                                                       os.path.basename(lineage))
        self.record_placement_file_versions(run_folder, placement_file_versions)
        return

    def record_placement_file_versions(self, run_folder, placement_file_versions):
        try:
            with open(os.path.join(run_folder, "short_summary.txt"), "a") as summary_file:
                summary_file.write("\nPlacement file versions:\n")
                for placement_file in placement_file_versions:
                    summary_file.write("{}\n".format(placement_file))
        except OSError:
            pass
        return



    def _run_3_datasets(self, mollicutes_runner=None):
        if mollicutes_runner:
            datasets = ["mycoplasmatales", "entomoplasmatales"]
            self.s_buscos = [mollicutes_runner.analysis.single_copy]
            self.d_buscos = [mollicutes_runner.analysis.multi_copy]
            self.f_buscos = [mollicutes_runner.analysis.only_fragments]
            self.s_percents = [mollicutes_runner.analysis.s_percent]
            self.d_percents = [mollicutes_runner.analysis.d_percent]
            self.f_percents = [mollicutes_runner.analysis.f_percent]
            dataset_runners = [mollicutes_runner]
        else:
            datasets = ["mollicutes", "mycoplasmatales", "entomoplasmatales"]
            self.s_buscos = []
            self.d_buscos = []
            self.f_buscos = []
            self.s_percents = []
            self.d_percents = []
            self.f_percents = []
            dataset_runners = []
        dataset_runners += self.run_lineages_list(datasets)
        self.get_best_match_lineage(dataset_runners)
        return