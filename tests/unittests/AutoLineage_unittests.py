import unittest
from unittest.mock import patch, call, Mock
from busco import AutoLineage, BuscoRunner

# class Options:
#
#     @staticmethod
#     def auto_lineage_options(key, value):
#         if key == "busco_run" and value == "auto-lineage-prok":
#             return False
#         if key == "busco_run" and value == "auto-lineage-prok":
#             return False

class TestAutoLineage(unittest.TestCase):

    def setUp(self):
        pass

    @patch("busco.BuscoConfig.BuscoConfigMain", autospec=True)
    def test_init_autolineage(self, mock_config_main):
        with self.assertLogs(AutoLineage.logger, 20):
            AutoLineage.AutoSelectLineage(mock_config_main)

    @patch("busco.BuscoConfig.BuscoConfigMain.getboolean", side_effect=[True, False, True, False, False])
    @patch("busco.AutoLineage.logger.info")
    @patch("busco.AutoLineage.os")
    @patch("busco.AutoLineage.BuscoRunner")
    @patch("busco.AutoLineage.AutoSelectLineage.get_best_match_lineage")
    @patch("busco.AutoLineage.AutoSelectLineage.run_lineages_list")
    @patch("busco.BuscoConfig.BuscoConfigMain", autospec=True)
    def test_run_auto_selector_lineage_lists(self, mock_config_main, mock_run_lineages_list, *args):
        for _ in range(3):
            asl = AutoLineage.AutoSelectLineage(mock_config_main)
            asl.selected_runner = Mock()
            asl.run_auto_selector()

        calls = [call(["archaea", "bacteria"]),
                 call(["eukaryota"]),
                 call(["archaea", "bacteria", "eukaryota"])]
        mock_run_lineages_list.assert_has_calls(calls, any_order=True)

    @patch("busco.AutoLineage.logger.info")
    @patch("__main__.AutoLineage_unittests.AutoLineage.BuscoRunner")
    @patch("busco.AutoLineage.BuscoConfigAuto", autospec=True)
    @patch("busco.BuscoConfig.BuscoConfigMain")
    def test_run_lineages_initializes_BuscoConfigAuto(self, mock_config_main, mock_config_auto, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        test_lineages = ["a", "b", "c"]
        test_dataset_version = "<dataset_version>"
        asl.dataset_version = test_dataset_version
        asl.run_lineages_list(test_lineages)
        calls = [call(mock_config_main, "{}_{}".format(test_lineages[0], test_dataset_version)),
                 call(mock_config_main, "{}_{}".format(test_lineages[1], test_dataset_version)),
                 call(mock_config_main, "{}_{}".format(test_lineages[2], test_dataset_version))]
        mock_config_auto.assert_has_calls(calls, any_order=True)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.AutoLineage.BuscoConfigAuto", autospec=True)
    @patch("busco.BuscoConfig.BuscoConfigMain")
    @patch("__main__.AutoLineage_unittests.AutoLineage.BuscoRunner")
    def test_run_lineages_initializes_BuscoRunner(self, mock_runner, mock_config_main, mock_config_auto, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        test_lineages = ["a", "b", "c"]
        test_dataset_version = "<dataset_version>"
        asl.dataset_version = test_dataset_version
        asl.run_lineages_list(test_lineages)
        mock_runner.assert_called_with(mock_config_auto.return_value)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.AutoLineage.BuscoConfigAuto", autospec=True)
    @patch("busco.BuscoConfig.BuscoConfigMain")
    @patch("__main__.AutoLineage_unittests.AutoLineage.BuscoRunner")
    def test_run_lineages_runs_analysis(self, mock_runner, mock_config_main, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        test_lineages = ["a", "b", "c"]
        test_dataset_version = "<dataset_version>"
        asl.dataset_version = test_dataset_version
        asl.run_lineages_list(test_lineages)
        mock_runner.return_value.run_analysis.assert_called_with(callback=asl.record_results)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.AutoLineage.BuscoConfigAuto", autospec=True)
    @patch("busco.BuscoConfig.BuscoConfigMain")
    @patch("__main__.AutoLineage_unittests.AutoLineage.BuscoRunner")
    def test_run_lineages_returns_runners(self, mock_runner, mock_config_main, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        test_lineages = ["a", "b", "c"]
        test_dataset_version = "<dataset_version>"
        asl.dataset_version = test_dataset_version
        runners = asl.run_lineages_list(test_lineages)
        self.assertEqual(runners, [mock_runner.return_value, mock_runner.return_value, mock_runner.return_value])

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.AutoLineage.BuscoConfigAuto", autospec=True)
    @patch("__main__.AutoLineage_unittests.AutoLineage.BuscoRunner")
    @patch("busco.BuscoConfig.BuscoConfigMain")
    def test_record_results_first_run(self, mock_config_main, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        asl.record_results(0,1,2,0,1,2)
        self.assertGreater(len(asl.s_buscos), 0)
        self.assertGreater(len(asl.d_buscos), 0)
        self.assertGreater(len(asl.f_buscos), 0)
        self.assertGreater(len(asl.s_percents), 0)
        self.assertGreater(len(asl.d_percents), 0)
        self.assertGreater(len(asl.f_percents), 0)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.AutoLineage.BuscoConfigAuto", autospec=True)
    @patch("__main__.AutoLineage_unittests.AutoLineage.BuscoRunner")
    @patch("busco.BuscoConfig.BuscoConfigMain")
    def test_record_results_multiple_runs(self, mock_config_main, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        asl.record_results(0, 1, 2, 0, 1, 2)
        asl.record_results(0, 1, 2, 0, 1, 2)
        self.assertGreater(len(asl.s_buscos), 1)
        self.assertGreater(len(asl.d_buscos), 1)
        self.assertGreater(len(asl.f_buscos), 1)
        self.assertGreater(len(asl.s_percents), 1)
        self.assertGreater(len(asl.d_percents), 1)
        self.assertGreater(len(asl.f_percents), 1)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoConfig.BuscoConfigMain")
    def test_evaluate_single_runner(self, mock_config_main, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        asl.s_buscos = [1]
        asl.f_buscos = [1]
        asl.d_buscos = [1]
        max_ind = asl.evaluate()
        self.assertEqual(max_ind, 0)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoConfig.BuscoConfigMain")
    def test_evaluate_multiple_runners(self, mock_config_main, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        asl.s_buscos = [10, 15, 12]
        asl.d_buscos = [5, 5, 5]
        max_ind = asl.evaluate()
        self.assertEqual(max_ind, 1)
        asl.s_buscos = [10, 10, 10]
        asl.d_buscos = [5, 6, 7]
        max_ind = asl.evaluate()
        self.assertEqual(max_ind, 2)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoConfig.BuscoConfigMain")
    def test_evaluate_first_order_tiebreak(self, mock_config_main, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        asl.s_buscos = [10, 10, 12]
        asl.d_buscos = [5, 5, 0]
        asl.f_buscos = [1, 2, 3]
        max_ind = asl.evaluate()
        self.assertEqual(max_ind, 1)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoConfig.BuscoConfigMain")
    def test_evaluate_second_order_tiebreak(self, mock_config_main, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        asl.s_buscos = [10, 10, 12, 14]
        asl.d_buscos = [5, 5, 0, 1]
        asl.f_buscos = [1, 2, 3, 2]
        asl.s_percents = [20, 40, 60, 80]
        max_ind = asl.evaluate()
        self.assertEqual(max_ind, 3)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoConfig.BuscoConfigMain")
    def test_evaluate_third_order_tiebreak(self, mock_config_main, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        asl.s_buscos = [10, 10, 12, 14]
        asl.d_buscos = [5, 5, 0, 1]
        asl.f_buscos = [1, 2, 3, 2]
        asl.s_percents = [20, 80, 60, 80]
        with self.assertLogs(AutoLineage.logger, "WARNING"):
            max_ind = asl.evaluate()
        self.assertEqual(max_ind, 1)

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoConfig.BuscoConfigMain")
    def test_evaluate_zero_results(self, mock_config_main, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        asl.s_buscos = [0, 0, 0, 0]
        asl.d_buscos = [0, 0, 0, 0]
        asl.f_buscos = [0, 0, 0, 0]
        with self.assertRaises(SystemExit):
            asl.evaluate()

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.AutoLineage.AutoSelectLineage.cleanup_disused_runs")
    @patch("__main__.AutoLineage_unittests.BuscoRunner.BuscoRunner.mode_dict")
    @patch("busco.BuscoConfig.BuscoConfigMain", autospec=True)
    def test_get_best_match_lineage(self, mock_config_main, fake_modedict, mock_cleanup, *args):
        mock_config_main.get.side_effect = [None]

        mock_config1 = Mock()
        mock_config2 = Mock()
        mock_config3 = Mock()
        mock_config1.get.side_effect = ["tran", None, "test1"]
        mock_config2.get.side_effect = ["tran", None, "test2"]
        mock_config3.get.side_effect = ["tran", None, "test3"]


        mock_analysis1 = Mock()
        mock_analysis2 = Mock()
        mock_analysis3 = Mock()
        mock_analysis1.return_value.hmmer_runner.single_copy = 75
        mock_analysis2.return_value.hmmer_runner.single_copy = 85
        mock_analysis3.return_value.hmmer_runner.single_copy = 80
        mock_analysis1.return_value.hmmer_runner.multi_copy = 5
        mock_analysis2.return_value.hmmer_runner.multi_copy = 6
        mock_analysis3.return_value.hmmer_runner.multi_copy = 7

        fake_modedict.__getitem__.side_effect = [mock_analysis1, mock_analysis2, mock_analysis3]

        runner1 = BuscoRunner.BuscoRunner(mock_config1)
        runner2 = BuscoRunner.BuscoRunner(mock_config2)
        runner3 = BuscoRunner.BuscoRunner(mock_config3)

        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        asl.get_best_match_lineage([runner1, runner2, runner3])
        self.assertEqual(asl.best_match_lineage_dataset, "test2")
        self.assertEqual(asl.selected_runner, runner2)
        mock_cleanup.assert_called()

    @patch("busco.AutoLineage.logger.info")
    @patch("busco.BuscoConfig.BuscoConfigMain", autospec=True)
    def test_cleanup_disused_runs(self, mock_config_main, *args):
        asl = AutoLineage.AutoSelectLineage(mock_config_main)
        mock_runner1 = Mock()
        mock_runner2 = Mock()
        asl.cleanup_disused_runs([mock_runner1, mock_runner2])
        mock_runner1.analysis.cleanup.assert_called()
        mock_runner2.analysis.cleanup.assert_called()


    def tearDown(self):
        pass