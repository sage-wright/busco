import unittest
from unittest.mock import patch, call, Mock
from busco import ConfigManager
import os
import importlib


class TestConfigManager(unittest.TestCase):
    def setUp(self):
        self.base_config = "config/config.ini"
        self.params = {"config_file": "path/to/config.ini"}

        # In order to suppress logs, we need to patch the LogDecorator object that handles the log messages and the
        # loggers from each module. The following code was adapted from a StackOverflow answer here:
        # https://stackoverflow.com/a/37890916/4844311

        # Do cleanup first so it is ready if an exception is raised
        def kill_patches():  # Create a cleanup callback that undoes our patches
            patch.stopall()  # Stops all patches started with start()
            importlib.reload(
                ConfigManager
            )  # Reload our UUT module which restores the original decorator

        self.addCleanup(
            kill_patches
        )  # We want to make sure this is run so we do this in addCleanup instead of tearDown

        # Now patch the decorator where the decorator is being imported from
        patch(
            "busco.BuscoLogger.LogDecorator", lambda *args, **kwargs: lambda f: f
        ).start()  # The lambda makes our decorator into a pass-thru. Also, don't forget to call start()
        # HINT: if you're patching a decor with params use something like:
        # lambda *x, **y: lambda f: f
        importlib.reload(
            ConfigManager
        )  # Reloads the module which applies our patched decorator

    def test_get_config_from_params(self):

        config_manager = ConfigManager.BuscoConfigManager(self.params)
        self.assertEqual(config_manager.config_file, self.params["config_file"])

    def test_get_config_from_env(self):
        os.environ["BUSCO_CONFIG_FILE"] = "path/to/config.ini"
        with patch("os.access") as mockaccess:
            mockaccess.return_value = True
            config_manager = ConfigManager.BuscoConfigManager({})
            self.assertEqual(
                config_manager.config_file, os.environ.get("BUSCO_CONFIG_FILE")
            )

    @patch("__main__.ConfigManager_unittests.ConfigManager.logger", autospec=True)
    def test_config_validated(self, *args):
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        with patch(
            "__main__.ConfigManager_unittests.ConfigManager.BuscoConfigMain",
            autospec=True,
        ):
            config_manager.load_busco_config()
            config_manager.config.validate.assert_called()

    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.BuscoConfigMain", autospec=True
    )
    def test_log_warning_if_neither_lineage_nor_autolineage_specified(
        self, mock_config_main, *args
    ):
        mock_config_main.return_value.check_lineage_present = lambda: False
        mock_config_main.return_value.getboolean.side_effect = [
            False,
            False,
            True,
            False,
            True,
        ]  # These values comprise the three distinct logical paths into this part
        # of the code: run1: False, False; run2: True, (shortcircuit); run3: False, True.
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        test_dataset_path = "/path/to/lineage_dataset"
        with patch.object(
            config_manager, "auto_select_lineage", lambda: test_dataset_path
        ):
            with self.assertLogs(ConfigManager.logger, "WARNING"):
                config_manager.load_busco_config()
            with patch("busco.ConfigManager.logger.warning") as mock_logger:
                for _ in range(2):
                    config_manager.load_busco_config()
                    self.assertFalse(mock_logger.called)

    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.BuscoConfigMain", autospec=True
    )
    def test_log_warning_if_both_lineage_and_autolineage_specified(
        self, mock_config_main, *args
    ):
        mock_config_main.return_value.check_lineage_present = lambda: True
        mock_config_main.return_value.getboolean.side_effect = [
            True,
            False,
            True,
            False,
            False,
        ]
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        test_dataset_path = "/path/to/lineage_dataset"
        with patch.object(
            config_manager, "auto_select_lineage", lambda: test_dataset_path
        ):
            for _ in range(2):
                with self.assertLogs(ConfigManager.logger, "WARNING"):
                    config_manager.load_busco_config()
            with patch("busco.ConfigManager.logger.warning") as mock_logger:
                config_manager.load_busco_config()
                self.assertFalse(mock_logger.called)

    patch("busco.ConfigManager.logger.warning")

    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.BuscoConfigMain", autospec=True
    )
    def test_config_updated_if_lineage_missing(self, mock_config_main, *args):
        mock_config_main.return_value.check_lineage_present = lambda: False
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        test_dataset_path = "/path/to/lineage_dataset"
        with patch.object(
            config_manager, "auto_select_lineage", lambda: test_dataset_path
        ):
            config_manager.load_busco_config()
        calls = [
            call("busco_run", "auto-lineage", "True"),
            call("busco_run", "lineage_dataset", test_dataset_path),
        ]
        mock_config_main.return_value.set.assert_has_calls(calls, any_order=True)

    @patch("busco.ConfigManager.logger.warning")
    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.BuscoConfigMain", autospec=True
    )
    def test_config_updated_if_lineage_present(self, mock_config_main, *args):
        mock_config_main.return_value.check_lineage_present = lambda: True
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        test_dataset_path = "/path/to/lineage_dataset"
        with patch.object(
            config_manager, "auto_select_lineage", lambda: test_dataset_path
        ):
            config_manager.load_busco_config()
        calls = [
            call("busco_run", "auto-lineage", "False"),
            call("busco_run", "auto-lineage-prok", "False"),
            call("busco_run", "auto-lineage-euk", "False"),
        ]
        mock_config_main.return_value.set.assert_has_calls(calls, any_order=True)

    @patch("busco.ConfigManager.logger.warning")
    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.BuscoConfigMain", autospec=True
    )
    def test_update_dirname(self, mock_config_main, *args):
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        test_dataset_path = "/path/to/lineage_dataset"
        mock_config_main.return_value.get = lambda *args: test_dataset_path
        with patch.object(
            config_manager, "auto_select_lineage", lambda: test_dataset_path
        ):
            config_manager.load_busco_config()
        mock_config_main.return_value.set_results_dirname.assert_called_with(
            test_dataset_path
        )

    @patch("busco.ConfigManager.logger.warning")
    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.BuscoConfigMain", autospec=True
    )
    def test_lineage_downloaded(self, mock_config_main, *args):
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        test_dataset_path = "/path/to/lineage_dataset"
        mock_config_main.return_value.get = lambda *args: test_dataset_path
        with patch.object(
            config_manager, "auto_select_lineage", lambda: test_dataset_path
        ):
            config_manager.load_busco_config()
        mock_config_main.return_value.download_lineage_file.assert_called_with(
            test_dataset_path
        )

    @patch("busco.ConfigManager.logger.warning")
    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.BuscoConfigMain", autospec=True
    )
    def test_lineage_dataset_config_loaded(self, mock_config_main, *args):
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        test_dataset_path = "/path/to/lineage_dataset"
        with patch.object(
            config_manager, "auto_select_lineage", lambda: test_dataset_path
        ):
            config_manager.load_busco_config()
        mock_config_main.return_value.load_dataset_config.assert_called()

    @patch("busco.ConfigManager.logger.warning")
    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.BuscoConfigMain", autospec=True
    )
    def test_run_auto_select_if_no_lineage(self, mock_config_main, *args):
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        mock_config_main.return_value.check_lineage_present = lambda: False
        test_dataset_path = "/path/to/lineage_dataset"
        with patch.object(
            config_manager, "auto_select_lineage", return_value=test_dataset_path
        ):
            config_manager.load_busco_config()
            config_manager.auto_select_lineage.assert_called()

    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.AutoSelectLineage",
        autospec=True,
    )
    def test_auto_select_lineage_call_function_initializes_asl(self, mock_asl):
        mock_asl.return_value.best_match_lineage_dataset = Mock()
        mock_asl.return_value.selected_runner = Mock()
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        config_manager.auto_select_lineage()
        mock_asl.assert_called()

    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.AutoSelectLineage",
        autospec=True,
    )
    def test_auto_select_lineage_call_function_runs_asl(self, mock_asl):
        mock_asl.return_value.best_match_lineage_dataset = Mock()
        mock_asl.return_value.selected_runner = Mock()
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        config_manager.auto_select_lineage()
        mock_asl.return_value.run_auto_selector.assert_called()

    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.AutoSelectLineage",
        autospec=True,
    )
    def test_auto_select_lineage_call_function_gets_lineage_dataset(self, mock_asl):
        mock_asl.return_value.best_match_lineage_dataset = Mock()
        mock_asl.return_value.selected_runner = Mock()
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        config_manager.auto_select_lineage()
        mock_asl.return_value.get_lineage_dataset.assert_called()

    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.AutoSelectLineage",
        autospec=True,
    )
    def test_auto_select_lineage_call_function_returns_lineage_dataset(self, mock_asl):
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        lineage_dataset = "best_match_dataset"
        mock_asl.return_value.best_match_lineage_dataset = lineage_dataset
        mock_asl.return_value.selected_runner = Mock()
        retval = config_manager.auto_select_lineage()
        self.assertEqual(retval, lineage_dataset)

    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.AutoSelectLineage",
        autospec=True,
    )
    def test_auto_select_lineage_call_function_selects_runner(self, mock_asl):
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        mock_asl.return_value.best_match_lineage_dataset = Mock()
        mock_asl.return_value.selected_runner = "test"
        config_manager.auto_select_lineage()
        self.assertEqual("test", config_manager.runner)

    def tearDown(self):
        pass
