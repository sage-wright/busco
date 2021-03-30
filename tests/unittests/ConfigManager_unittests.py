import unittest
from unittest.mock import patch
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
    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.BuscoConfigMain", autospec=True
    )
    def test_config_main_configured(self, mock_config, *args):
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        config_manager.load_busco_config_main()
        mock_config.return_value.configure.assert_called()

    @patch("__main__.ConfigManager_unittests.ConfigManager.logger", autospec=True)
    @patch(
        "__main__.ConfigManager_unittests.ConfigManager.BuscoConfigMain", autospec=True
    )
    def test_config_main_validated(self, mock_config, *args):
        config_manager = ConfigManager.BuscoConfigManager(self.params)
        config_manager.load_busco_config_main()
        mock_config.return_value.validate.assert_called()

    def tearDown(self):
        pass
