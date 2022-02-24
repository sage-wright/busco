import unittest
from unittest.mock import patch, call, Mock
from busco import BuscoRunner, ConfigManager


class TestAutoLineage(unittest.TestCase):
    def setUp(self):
        pass

    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_config_updated_if_lineage_missing(self, mock_config_manager, *args):
        mock_config_manager.config_main.check_lineage_present = lambda: False
        runner = BuscoRunner.SingleRunner(mock_config_manager)
        analysis_runner = Mock()
        test_dataset_path = "/path/to/lineage_dataset"
        test_parent = "parent"
        with patch.object(
            runner,
            "auto_select_lineage",
            lambda: (test_dataset_path, analysis_runner, test_parent),
        ):
            runner.get_lineage()
        calls = [
            call("busco_run", "lineage_dataset", test_dataset_path),
        ]
        mock_config_manager.config_main.set.assert_has_calls(calls, any_order=True)

    @patch("busco.BuscoRunner.logger.info")
    @patch(
        "busco.AutoLineage.AutoSelectLineage",
        autospec=True,
    )
    def test_auto_select_lineage_call_function_initializes_asl(self, mock_asl, *args):
        config = Mock()
        mock_asl.return_value.best_match_lineage_dataset = Mock()
        mock_asl.return_value.selected_runner = Mock()
        runner = BuscoRunner.SingleRunner(config)
        runner.auto_select_lineage()
        mock_asl.assert_called()

    @patch("busco.BuscoRunner.logger.info")
    @patch(
        "busco.AutoLineage.AutoSelectLineage",
        autospec=True,
    )
    def test_auto_select_lineage_call_function_runs_asl(self, mock_asl, *args):
        config = Mock()
        mock_asl.return_value.best_match_lineage_dataset = Mock()
        mock_asl.return_value.selected_runner = Mock()
        runner = BuscoRunner.SingleRunner(config)
        runner.auto_select_lineage()
        mock_asl.return_value.run_auto_selector.assert_called()

    @patch("busco.BuscoRunner.logger.info")
    @patch(
        "busco.AutoLineage.AutoSelectLineage",
        autospec=True,
    )
    def test_auto_select_lineage_call_function_gets_lineage_dataset(
        self, mock_asl, *args
    ):
        config = Mock()
        mock_asl.return_value.best_match_lineage_dataset = Mock()
        mock_asl.return_value.selected_runner = Mock()
        runner = BuscoRunner.SingleRunner(config)
        runner.auto_select_lineage()
        mock_asl.return_value.get_lineage_dataset.assert_called()

    @patch("busco.BuscoRunner.logger.info")
    @patch(
        "busco.AutoLineage.AutoSelectLineage",
        autospec=True,
    )
    def test_auto_select_lineage_call_function_returns_lineage_dataset(
        self, mock_asl, *args
    ):
        config = Mock()
        lineage_dataset = "best_match_dataset"
        mock_asl.return_value.best_match_lineage_dataset = lineage_dataset
        mock_asl.return_value.selected_runner = Mock()
        runner = BuscoRunner.SingleRunner(config)
        retval_dataset, retval_runner, retval_parent = runner.auto_select_lineage()
        self.assertEqual(retval_dataset, lineage_dataset)

    @patch("busco.BuscoRunner.logger.info")
    @patch(
        "busco.AutoLineage.AutoSelectLineage",
        autospec=True,
    )
    def test_auto_select_lineage_call_function_selects_runner(self, mock_asl, *args):
        config = Mock()
        lineage_dataset = "best_match_dataset"
        mock_asl.return_value.best_match_lineage_dataset = lineage_dataset
        mock_asl.return_value.selected_runner = Mock()
        runner = BuscoRunner.SingleRunner(config)
        retval_dataset, retval_runner, retval_parent = runner.auto_select_lineage()
        self.assertEqual(retval_runner, mock_asl.return_value.selected_runner)

    @patch("busco.BuscoRunner.logger.info")
    @patch(
        "busco.AutoLineage.AutoSelectLineage",
        autospec=True,
    )
    def test_auto_select_lineage_call_function_returns_parent(self, mock_asl, *args):
        config = Mock()
        lineage_dataset = "best_match_dataset"
        mock_asl.return_value.best_match_lineage_dataset = lineage_dataset
        mock_asl.return_value.selected_runner = Mock()
        mock_asl.return_value.selected_runner.config.get.side_effect = ["parent"]
        runner = BuscoRunner.SingleRunner(config)
        retval_dataset, retval_runner, retval_parent = runner.auto_select_lineage()
        self.assertEqual(retval_parent, "parent")

    @patch("busco.ConfigManager.BuscoConfigManager")
    def test_config_set_parent_dataset_if_not_virus(self, mock_config_manager, *args):
        test_dataset_path = "/path/to/lineage_dataset"
        test_parent_dataset = "/path/to/parent_dataset"
        runner = BuscoRunner.SingleRunner(mock_config_manager)
        analysis_runner = Mock()
        analysis_runner.config.get.side_effect = ["bacteria", test_parent_dataset]
        with patch.object(
            runner,
            "auto_select_lineage",
            lambda: (test_dataset_path, analysis_runner, test_parent_dataset),
        ):
            runner.get_lineage()
        calls = [
            call("busco_run", "domain_run_name", test_parent_dataset),
        ]
        mock_config_manager.config_main.set.assert_has_calls(calls, any_order=True)

    def tearDown(self):
        pass
