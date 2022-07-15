import unittest
from unittest.mock import patch, call, Mock
from busco import AutoLineage, BuscoRunner
from busco.ConfigManager import BuscoConfigManager


class TestAutoLineage(unittest.TestCase):
    def setUp(self):
        pass

    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_init_autolineage(self, mock_config_manager):
        with self.assertLogs(AutoLineage.logger, 20):
            AutoLineage.AutoSelectLineage(mock_config_manager)

    @patch("busco.AutoLineage.AutoSelectLineage.virus_check", return_value=False)
    @patch("busco.AutoLineage.logger.info")
    @patch("busco.AutoLineage.os")
    @patch("busco.AutoLineage.AnalysisRunner")
    @patch("busco.AutoLineage.AutoSelectLineage.get_best_match_lineage")
    @patch("busco.AutoLineage.AutoSelectLineage.run_lineages_list")
    def test_run_auto_selector_lineage_lists_no_virus(
        self, mock_run_lineages_list, *args
    ):
        config_manager = BuscoConfigManager({})
        config_manager.config_main = Mock()
        config_manager.config_main.getboolean = Mock(
            side_effect=[True, False, True, False, False]
        )
        for _ in range(3):
            asl = AutoLineage.AutoSelectLineage(config_manager)
            asl.selected_runner = Mock()
            asl.selected_runner.analysis.hmmer_runner.single_copy_buscos = [
                0
            ]  # avoid SystemExit with empty HMMER results
            asl.selected_runner.analysis.hmmer_runner.multi_copy_buscos = [0]
            asl.selected_runner.analysis.hmmer_runner.fragmented_buscos = [0]
            asl.run_auto_selector()

        calls = [
            call(["archaea", "bacteria"]),
            call(["eukaryota"]),
            call(["archaea", "bacteria", "eukaryota"]),
        ]
        mock_run_lineages_list.assert_has_calls(calls, any_order=True)

    @patch("busco.AutoLineage.logger.info")
    @patch("__main__.AutoLineage_unittests.AutoLineage.AnalysisRunner")
    @patch("busco.ConfigManager.BuscoConfigAuto", autospec=True)
    def test_run_lineages_initializes_BuscoConfigAuto(self, mock_config_auto, *args):
        config_manager = BuscoConfigManager({})
        config_manager.config_main = Mock()
        asl = AutoLineage.AutoSelectLineage(config_manager)
        test_lineages = ["a", "b", "c"]
        test_dataset_version = "<dataset_version>"
        asl.dataset_version = test_dataset_version
        asl.run_lineages_list(test_lineages)
        calls = [
            call(
                config_manager.config_main,
                "{}_{}".format(test_lineages[0], test_dataset_version),
            ),
            call(
                config_manager.config_main,
                "{}_{}".format(test_lineages[1], test_dataset_version),
            ),
            call(
                config_manager.config_main,
                "{}_{}".format(test_lineages[2], test_dataset_version),
            ),
        ]
        mock_config_auto.assert_has_calls(calls, any_order=True)

    @patch("busco.AutoLineage.logger.info")
    @patch("__main__.AutoLineage_unittests.BuscoConfigManager.load_busco_config_auto")
    @patch("__main__.AutoLineage_unittests.AutoLineage.AnalysisRunner")
    def test_run_lineages_initializes_BuscoRunner(
        self, mock_runner, mock_config_auto, *args
    ):
        config_manager = BuscoConfigManager({})
        config_manager.config_main = Mock()
        asl = AutoLineage.AutoSelectLineage(config_manager)
        test_lineages = ["a", "b", "c"]
        test_dataset_version = "<dataset_version>"
        asl.dataset_version = test_dataset_version
        asl.run_lineages_list(test_lineages)
        mock_runner.assert_called_with(mock_config_auto.return_value)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoConfig.BuscoConfigAuto", autospec=True)
    @patch("busco.ConfigManager.BuscoConfigManager")
    @patch("__main__.AutoLineage_unittests.AutoLineage.AnalysisRunner")
    def test_run_lineages_runs_analysis(self, mock_runner, mock_config_manager, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_manager.config_main)
        test_lineages = ["a", "b", "c"]
        test_dataset_version = "<dataset_version>"
        asl.dataset_version = test_dataset_version
        asl.run_lineages_list(test_lineages)
        mock_runner.return_value.run_analysis.assert_called_with(
            callback=asl.record_results
        )

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoConfig.BuscoConfigAuto", autospec=True)
    @patch("busco.ConfigManager.BuscoConfigManager")
    @patch("__main__.AutoLineage_unittests.AutoLineage.AnalysisRunner")
    def test_run_lineages_returns_runners(
        self, mock_runner, mock_config_manager, *args
    ):
        asl = AutoLineage.AutoSelectLineage(mock_config_manager.config_main)
        test_lineages = ["a", "b", "c"]
        test_dataset_version = "<dataset_version>"
        asl.dataset_version = test_dataset_version
        runners = asl.run_lineages_list(test_lineages)
        self.assertEqual(
            runners,
            [
                mock_runner.return_value,
                mock_runner.return_value,
                mock_runner.return_value,
            ],
        )

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoConfig.BuscoConfigAuto", autospec=True)
    @patch("__main__.AutoLineage_unittests.AutoLineage.AnalysisRunner")
    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_record_results_first_run(self, mock_config_manager, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_manager.config_main)
        asl.record_results(0, 1, 2, 0, 1, 2)
        self.assertGreater(len(asl.s_buscos), 0)
        self.assertGreater(len(asl.d_buscos), 0)
        self.assertGreater(len(asl.f_buscos), 0)
        self.assertGreater(len(asl.s_percents), 0)
        self.assertGreater(len(asl.d_percents), 0)
        self.assertGreater(len(asl.f_percents), 0)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoConfig.BuscoConfigAuto", autospec=True)
    @patch("__main__.AutoLineage_unittests.AutoLineage.AnalysisRunner")
    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_record_results_multiple_runs(self, mock_config_manager, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_manager.config_main)
        asl.record_results(0, 1, 2, 0, 1, 2)
        asl.record_results(0, 1, 2, 0, 1, 2)
        self.assertGreater(len(asl.s_buscos), 1)
        self.assertGreater(len(asl.d_buscos), 1)
        self.assertGreater(len(asl.f_buscos), 1)
        self.assertGreater(len(asl.s_percents), 1)
        self.assertGreater(len(asl.d_percents), 1)
        self.assertGreater(len(asl.f_percents), 1)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_evaluate_single_runner(self, mock_config_manager, *args):
        runner1 = Mock()
        runner1.analysis.hmmer_runner.single_copy = 1
        runner1.analysis.hmmer_runner.multi_copy = 1
        runner1.analysis.hmmer_runner.only_fragments = 1
        asl = AutoLineage.AutoSelectLineage(mock_config_manager.config_main)
        max_ind = asl.evaluate([runner1])
        self.assertEqual(max_ind, 0)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_evaluate_multiple_runners(self, mock_config_manager, *args):
        runner1 = Mock()
        runner2 = Mock()
        runner3 = Mock()
        runner1.analysis.hmmer_runner.single_copy = 10
        runner1.analysis.hmmer_runner.multi_copy = 5
        runner2.analysis.hmmer_runner.single_copy = 15
        runner2.analysis.hmmer_runner.multi_copy = 5
        runner3.analysis.hmmer_runner.single_copy = 12
        runner3.analysis.hmmer_runner.multi_copy = 5
        asl = AutoLineage.AutoSelectLineage(mock_config_manager.config_main)
        max_ind = asl.evaluate([runner1, runner2, runner3])
        self.assertEqual(max_ind, 1)

        runner2.analysis.hmmer_runner.single_copy = 10
        runner2.analysis.hmmer_runner.multi_copy = 6
        runner3.analysis.hmmer_runner.single_copy = 10
        runner3.analysis.hmmer_runner.multi_copy = 7
        max_ind = asl.evaluate([runner1, runner2, runner3])
        self.assertEqual(max_ind, 2)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_evaluate_first_order_tiebreak(self, mock_config_manager, *args):
        runner1 = Mock()
        runner2 = Mock()
        runner3 = Mock()
        runner1.analysis.hmmer_runner.single_copy = 10
        runner1.analysis.hmmer_runner.multi_copy = 5
        runner1.analysis.hmmer_runner.only_fragments = 1
        runner2.analysis.hmmer_runner.single_copy = 10
        runner2.analysis.hmmer_runner.multi_copy = 5
        runner2.analysis.hmmer_runner.only_fragments = 2
        runner3.analysis.hmmer_runner.single_copy = 12
        runner3.analysis.hmmer_runner.multi_copy = 0
        runner3.analysis.hmmer_runner.only_fragments = 3
        asl = AutoLineage.AutoSelectLineage(mock_config_manager.config_main)
        max_ind = asl.evaluate([runner1, runner2, runner3])
        self.assertEqual(max_ind, 1)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_evaluate_second_order_tiebreak(self, mock_config_manager, *args):
        runner1 = Mock()
        runner2 = Mock()
        runner3 = Mock()
        runner4 = Mock()
        runner1.analysis.hmmer_runner.single_copy = 10
        runner1.analysis.hmmer_runner.multi_copy = 5
        runner1.analysis.hmmer_runner.only_fragments = 1
        runner1.analysis.hmmer_runner.s_percent = 20
        runner2.analysis.hmmer_runner.single_copy = 10
        runner2.analysis.hmmer_runner.multi_copy = 5
        runner2.analysis.hmmer_runner.only_fragments = 2
        runner2.analysis.hmmer_runner.s_percent = 40
        runner3.analysis.hmmer_runner.single_copy = 12
        runner3.analysis.hmmer_runner.multi_copy = 0
        runner3.analysis.hmmer_runner.only_fragments = 3
        runner3.analysis.hmmer_runner.s_percent = 60
        runner4.analysis.hmmer_runner.single_copy = 14
        runner4.analysis.hmmer_runner.multi_copy = 1
        runner4.analysis.hmmer_runner.only_fragments = 2
        runner4.analysis.hmmer_runner.s_percent = 80
        asl = AutoLineage.AutoSelectLineage(mock_config_manager.config_main)
        max_ind = asl.evaluate([runner1, runner2, runner3, runner4])
        self.assertEqual(max_ind, 3)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_evaluate_third_order_tiebreak(self, mock_config_manager, *args):
        runner1 = Mock()
        runner2 = Mock()
        runner3 = Mock()
        runner4 = Mock()
        runner1.analysis.hmmer_runner.single_copy = 10
        runner1.analysis.hmmer_runner.multi_copy = 5
        runner1.analysis.hmmer_runner.only_fragments = 1
        runner1.analysis.hmmer_runner.s_percent = 20
        runner2.analysis.hmmer_runner.single_copy = 10
        runner2.analysis.hmmer_runner.multi_copy = 5
        runner2.analysis.hmmer_runner.only_fragments = 2
        runner2.analysis.hmmer_runner.s_percent = 80
        runner3.analysis.hmmer_runner.single_copy = 12
        runner3.analysis.hmmer_runner.multi_copy = 0
        runner3.analysis.hmmer_runner.only_fragments = 3
        runner3.analysis.hmmer_runner.s_percent = 60
        runner4.analysis.hmmer_runner.single_copy = 14
        runner4.analysis.hmmer_runner.multi_copy = 1
        runner4.analysis.hmmer_runner.only_fragments = 2
        runner4.analysis.hmmer_runner.s_percent = 80
        asl = AutoLineage.AutoSelectLineage(mock_config_manager.config_main)
        with self.assertLogs(AutoLineage.logger, "WARNING"):
            max_ind = asl.evaluate([runner1, runner2, runner3, runner4])
        self.assertEqual(max_ind, 1)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoRunner.AnalysisRunner.set_parent_dataset")
    @patch("busco.AutoLineage.AutoSelectLineage.cleanup_disused_runs")
    @patch("__main__.AutoLineage_unittests.BuscoRunner.AnalysisRunner.mode_dict")
    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_get_best_match_lineage(
        self, mock_config_manager, fake_modedict, mock_cleanup, *args
    ):
        mock_config_manager.config_main.get.side_effect = [None]

        mock_config1 = Mock()
        mock_config2 = Mock()
        mock_config3 = Mock()
        mock_config1.get.side_effect = [
            "test_input1",
            "prok_tran",
            "prokaryota",
            "test_lineage1",
            "test_lineage1",
        ]
        mock_config2.get.side_effect = [
            "test_input2",
            "prok_tran",
            "prokaryota",
            "test_lineage2",
            "test_lineage2",
        ]
        mock_config3.get.side_effect = [
            "test_input3",
            "euk_tran",
            "eukaryota",
            "test_lineage3",
            "test_lineage3",
        ]

        mock_analysis1 = Mock()
        mock_analysis2 = Mock()
        mock_analysis3 = Mock()
        mock_analysis1.return_value.hmmer_runner.single_copy = 75
        mock_analysis2.return_value.hmmer_runner.single_copy = 85
        mock_analysis3.return_value.hmmer_runner.single_copy = 80
        mock_analysis1.return_value.hmmer_runner.multi_copy = 5
        mock_analysis2.return_value.hmmer_runner.multi_copy = 6
        mock_analysis3.return_value.hmmer_runner.multi_copy = 7

        fake_modedict.__getitem__.side_effect = [
            mock_analysis1,
            mock_analysis2,
            mock_analysis3,
        ]
        mock_analysis1.return_value.run_folder = ""
        mock_analysis2.return_value.run_folder = ""
        mock_analysis3.return_value.run_folder = ""

        runner1 = BuscoRunner.AnalysisRunner(mock_config1)
        runner2 = BuscoRunner.AnalysisRunner(mock_config2)
        runner3 = BuscoRunner.AnalysisRunner(mock_config3)

        asl = AutoLineage.AutoSelectLineage(mock_config_manager.config_main)
        asl.get_best_match_lineage([runner1, runner2, runner3])
        self.assertEqual(asl.best_match_lineage_dataset, "test_lineage2")
        self.assertEqual(asl.selected_runner, runner2)
        mock_cleanup.assert_called_with([runner1, runner3])

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_cleanup_disused_runs(self, mock_config_manager, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_manager)
        mock_runner1 = Mock()
        mock_runner2 = Mock()
        mock_runner1.cleaned_up = False
        mock_runner2.cleaned_up = True
        asl.cleanup_disused_runs([mock_runner1, mock_runner2])
        mock_runner1.cleanup.assert_called()
        mock_runner2.cleanup.assert_not_called()

    # Todo: add tests for get_lineage_datasets, 3_dataset check and busco_placer step

    def tearDown(self):
        pass
