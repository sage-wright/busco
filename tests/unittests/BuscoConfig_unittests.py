import unittest
from busco import BuscoConfig
import shutil
import os
from unittest.mock import Mock
from unittest.mock import patch, call


class TestBuscoConfig(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.maxDiff = None
        self.base_config = "config/config.ini"

        self.params = {
            "auto-lineage": False,
            "auto-lineage-euk": False,
            "auto-lineage-prok": False,
            "config_file": None,
            "cpu": None,
            "evalue": None,
            "force": False,
            "help": "==SUPPRESS==",
            "in": None,
            "limit": None,
            "lineage_dataset": None,
            "list_datasets": "==SUPPRESS==",
            "mode": None,
            "offline": False,
            "out": None,
            "out_path": None,
            "quiet": False,
            "restart": False,
            "metaeuk_parameters": None,
            "metaeuk_rerun_parameters": None,
            "use_augustus": False,
            "augustus_parameters": None,
            "augustus_species": None,
            "long": False,
            "datasets_version": None,
            "download_base_url": None,
            "download_path": None,
            "update-data": False,
            "version": "==SUPPRESS==",
            "tar": False,
        }

        self.test_params = {
            "in": "input_test",
            "out": "output_test",
            "mode": "mode_test",
        }

        self.config_structure = {
            "augustus": ["path", "command"],
            "busco_run": [
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
            ],
            "etraining": ["path", "command"],
            "gff2gbSmallDNA.pl": ["path", "command"],
            "hmmsearch": ["path", "command"],
            "makeblastdb": ["path", "command"],
            "metaeuk": ["path", "command"],
            "new_species.pl": ["path", "command"],
            "optimize_augustus.pl": ["path", "command"],
            "prodigal": ["path", "command"],
            "sepp": ["path", "command"],
            "tblastn": ["path", "command"],
        }

    def test_read_config_file(self):
        config = BuscoConfig.BaseConfig()
        config.conf_file = self.base_config
        config._load_config_file(config.conf_file)
        self.assertIn("busco_run", config.sections())

    def test_read_config_file_ioerror(self):
        with self.assertRaises(BuscoConfig.BatchFatalError):
            config = BuscoConfig.BaseConfig()
            config.conf_file = "/path/not/found"
            config._load_config_file(config.conf_file)

    def test_read_config_file_parseerror(self):
        config_path = "tests/config_parseerror_test.ini"
        test_config_contents = "in=input_file\n"
        with open(config_path, "w") as f:
            f.write(test_config_contents)

        with self.assertRaises(BuscoConfig.BatchFatalError):
            config = BuscoConfig.BaseConfig()
            config.conf_file = config_path
            config._load_config_file(config.conf_file)
        os.remove(config_path)

    def test_read_config_file_duplicateerror(self):
        config_path = "tests/config_duplicate_test.ini"
        test_config_contents = "[busco_run]\n" "in=input_file\n" "in=input_file\n"
        with open(config_path, "w") as f:
            f.write(test_config_contents)

        with self.assertRaises(BuscoConfig.BatchFatalError):
            config = BuscoConfig.BaseConfig()
            config.conf_file = config_path
            config._load_config_file(config.conf_file)
        os.remove(config_path)

    def test_config_update_args_bool(self):
        update_params = {
            "force": True,
            "offline": True,
            "quiet": True,
            "restart": True,
        }
        config = BuscoConfig.BuscoConfigMain(self.base_config, update_params)
        config.configure()
        self.assertEqual(
            update_params,
            {key: config.getboolean("busco_run", key) for key in update_params.keys()},
        )

    def test_config_update_args_nonbool(self):
        update_params = {
            "cpu": "10",
            "evalue": "0.01",
            "in": "input_file",
            "limit": "1",
            "lineage_dataset": "test_odb10",
            "mode": "test",
            "out": "test",
            "out_path": "test",
        }
        config = BuscoConfig.BuscoConfigMain(self.base_config, update_params)
        config.configure()
        self.assertEqual(
            update_params,
            {key: config.get("busco_run", key) for key in update_params.keys()},
        )

    def test_config_default_params(self):
        correct_default_params = {
            "auto-lineage": False,
            "auto-lineage-euk": False,
            "auto-lineage-prok": False,
            "cpu": "1",
            "datasets_version": "odb10",
            "download_base_url": "https://busco-data.ezlab.org/v5/data/",
            "download_path": os.path.join(os.getcwd(), "busco_downloads"),
            "force": False,
            "offline": False,
            "out_path": os.getcwd(),
            "quiet": False,
            "restart": False,
            "update-data": False,
            "use_augustus": False,
        }
        config = BuscoConfig.BuscoConfigMain(
            self.base_config, {"lineage_dataset": "test"}
        )
        config.configure()
        config_default_filled = {
            key: config.get("busco_run", key) for key in correct_default_params
        }

        self.assertEqual(
            {key: str(val) for key, val in correct_default_params.items()},
            config_default_filled,
        )

    @patch(
        "busco.BuscoConfig.BuscoConfigMain.check_lineage_present",
        side_effect=[False, False],
    )
    @patch(
        "busco.BuscoConfig.BuscoConfigMain.getboolean",
        side_effect=[False, True, True, False, False, False, True, False],
    )
    def test_config_auto_lineage_settings(self, *args):
        for _ in range(2):
            config = BuscoConfig.BuscoConfigMain(self.base_config, {})
            config.configure()
            self.assertEqual(config.get("busco_run", "auto-lineage"), "True")

    @patch(
        "busco.BuscoConfig.BuscoConfigMain.check_lineage_present",
        side_effect=[False],
    )
    @patch(
        "busco.BuscoConfig.BuscoConfigMain.getboolean",
        side_effect=[False, False, True, True, True],
    )
    def test_config_auto_lineage_both_selected_warning(self, *args):
        with self.assertLogs(BuscoConfig.logger, "WARNING"):
            config = BuscoConfig.BuscoConfigMain(self.base_config, {})
            config.configure()
        self.assertEqual(config.get("busco_run", "auto-lineage-euk"), "False")
        self.assertEqual(config.get("busco_run", "auto-lineage-prok"), "False")
        self.assertEqual(config.get("busco_run", "auto-lineage"), "True")

    @patch(
        "busco.BuscoConfig.BuscoConfigMain.check_lineage_present",
        side_effect=[False],
    )
    @patch(
        "busco.BuscoConfig.BuscoConfigMain.getboolean",
        side_effect=[False, False, False, False],
    )
    def test_config_auto_lineage_none_selected_no_lineage_warning(self, *args):
        with self.assertLogs(BuscoConfig.logger, "WARNING"):
            config = BuscoConfig.BuscoConfigMain(self.base_config, {})
            config.configure()
        self.assertEqual(config.get("busco_run", "auto-lineage-euk"), "False")
        self.assertEqual(config.get("busco_run", "auto-lineage-prok"), "False")
        self.assertEqual(config.get("busco_run", "auto-lineage"), "True")

    @patch(
        "busco.BuscoConfig.BuscoConfigMain.check_lineage_present",
        side_effect=[True, True, True],
    )
    @patch(
        "busco.BuscoConfig.BuscoConfigMain.getboolean",
        side_effect=[False, True, False, True, False, False, True],
    )
    def test_config_auto_lineage_dataset_specified_warning(self, *args):
        for _ in range(3):
            with self.assertLogs(BuscoConfig.logger, "WARNING"):
                config = BuscoConfig.BuscoConfigMain(self.base_config, {})
                config.configure()
            self.assertEqual(config.get("busco_run", "auto-lineage-euk"), "False")
            self.assertEqual(config.get("busco_run", "auto-lineage-prok"), "False")
            self.assertEqual(config.get("busco_run", "auto-lineage"), "False")

    def test_mandatory_keys_check_log(self):
        with self.assertLogs(BuscoConfig.logger, 20):
            params_test = {"in": "input_file", "out": "output_name", "mode": "genome"}
            config = BuscoConfig.BuscoConfigMain(self.base_config, params_test)
            config.configure()
            config._check_mandatory_keys_exist()

    def test_mandatory_keys_check_missing_param_in(self):
        with self.assertRaises(BuscoConfig.BatchFatalError):
            params_test = {"out": "output_name", "mode": "genome"}
            config = BuscoConfig.BuscoConfigMain(self.base_config, params_test)
            config.configure()
            config._check_mandatory_keys_exist()

    def test_mandatory_keys_check_missing_param_mode(self):
        with self.assertRaises(BuscoConfig.BatchFatalError):
            params_test = {"in": "input_file", "out": "output_name"}
            config = BuscoConfig.BuscoConfigMain(self.base_config, params_test)
            config.configure()
            config._check_mandatory_keys_exist()

    def test_mandatory_keys_check_missing_param_out(self):
        with self.assertRaises(BuscoConfig.BatchFatalError):
            params_test = {"in": "input_file", "mode": "genome"}
            config = BuscoConfig.BuscoConfigMain(self.base_config, params_test)
            config.configure()
            config._check_mandatory_keys_exist()

    def test_previous_run_check_without_existing_run(self):
        output_dir = os.path.join(os.getcwd(), self.test_params["out"])
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        self.assertIsNone(config._check_no_previous_run())

    def test_previous_run_check_with_existing_run_no_force(self):
        previous_run_name = "test_busco_run_dir"
        os.makedirs(previous_run_name, exist_ok=True)
        self.test_params["out"] = previous_run_name
        self.test_params["force"] = "False"
        with self.assertRaises(BuscoConfig.BatchFatalError):
            config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
            config.configure()
            config._check_no_previous_run()
        shutil.rmtree(previous_run_name)

    def test_previous_run_check_with_existing_run_with_force_and_log(self):
        previous_run_name = "test_busco_run_dir"
        os.makedirs(previous_run_name, exist_ok=True)
        self.test_params["out"] = previous_run_name
        self.test_params["force"] = "True"
        with self.assertLogs(BuscoConfig.logger, 20):
            config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
            config.configure()
            config._check_no_previous_run()
            self.assertFalse(os.path.exists(previous_run_name))

        try:  # In case of test failure, remove tmp folder anyway
            shutil.rmtree(previous_run_name)
        except FileNotFoundError:
            pass

    def test_previous_run_check_without_existing_run_and_restart(self):
        self.test_params["restart"] = "True"
        with self.assertLogs(BuscoConfig.logger, "WARNING"):
            config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
            config.configure()
            config._check_no_previous_run()
            self.assertEqual(config.getboolean("busco_run", "restart"), False)

    def test_previous_run_check_with_existing_run_and_restart(self):
        previous_run_name = "test_busco_run_dir"
        os.makedirs(previous_run_name, exist_ok=True)
        self.test_params.update({"out": previous_run_name, "restart": True})
        with self.assertLogs(BuscoConfig.logger, "INFO"):
            config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
            config.configure()
            config._check_no_previous_run()
            self.assertEqual(config.getboolean("busco_run", "restart"), True)
        shutil.rmtree(previous_run_name)

    def test_create_required_paths(self):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        config.main_out = os.path.join(
            config.get("busco_run", "out_path"), config.get("busco_run", "out")
        )
        config._create_required_paths(config.main_out)
        output_dir = os.path.join(os.getcwd(), self.test_params["out"])
        self.assertTrue(os.path.exists(output_dir))
        shutil.rmtree(output_dir)

    def test_config_structure(self):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        self.assertEqual(
            set(config.PERMITTED_OPTIONS), set(self.config_structure["busco_run"])
        )

    def test_catch_disallowed_keys(self):
        for section_name in self.config_structure:
            with self.assertRaises(BuscoConfig.BatchFatalError):
                config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
                config.configure()
                config.set(section_name, "forbidden_option", "forbidden_value")
                config._check_allowed_keys()

    def test_out_value_check_invalid(self):
        self.test_params["out"] = "output/"
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        config._check_out_value()
        self.assertTrue(config.get("busco_run", "out") == "output")

    def test_out_value_check_valid(self):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        self.assertIsNone(config._check_out_value())

    def test_limit_value_out_of_range(self):
        for lim_val in [-1, 0, 25]:
            self.test_params["limit"] = lim_val
            with self.assertRaises(BuscoConfig.BatchFatalError):
                config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
                config.configure()
                config._check_limit_value()

    def test_limit_value_within_range(self):
        for lim_val in [1, 20]:
            self.test_params["limit"] = lim_val
            config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
            config.configure()
            self.assertIsNone(config._check_limit_value())

    def test_evalue_nondefault(self):
        self.test_params["evalue"] = 1
        with self.assertLogs(BuscoConfig.logger, level="WARNING"):
            config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
            config.configure()
            config._check_evalue()

    @patch("__main__.BuscoConfig_unittests.BuscoConfig.logger.warning")
    def test_evalue_default(self, mock_logger):
        self.test_params["evalue"] = 0.001
        self.test_params["lineage_dataset"] = "test"
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        config._check_evalue()
        mock_logger.assert_not_called()

    def test_expand_all_paths_tilde(self):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        config.set("busco_run", "download_path", "~/test_download_path")
        config._expand_all_paths()
        self.assertEqual(
            config.get("busco_run", "download_path"),
            os.path.expanduser("~/test_download_path"),
        )

    def test_expand_all_paths_relative_path_current_dir(self):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        config.set("busco_run", "out_path", "./test_out_path")
        config._expand_all_paths()
        self.assertEqual(
            config.get("busco_run", "out_path"), os.path.abspath("./test_out_path")
        )

    def test_expand_all_paths_relative_path_parent_dir(self):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        config.set("busco_run", "in", "../test_input_file")
        config._expand_all_paths()
        self.assertEqual(
            config.get("busco_run", "in"), os.path.abspath("../test_input_file")
        )

    def test_expand_all_paths_hmmsearch(self):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        config.set("hmmsearch", "path", "~/test_hmmsearch_path")
        config._expand_all_paths()
        self.assertEqual(
            config.get("hmmsearch", "path"), os.path.expanduser("~/test_hmmsearch_path")
        )

    @patch(
        "__main__.BuscoConfig_unittests.BuscoConfig.os.path.isdir", return_value=True
    )
    def test_batch_mode_true(self, *args):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.set = Mock()
        config._check_batch_mode()
        calls = [call("busco_run", "batch_mode", "True")]
        config.set.assert_has_calls(calls)

    @patch(
        "__main__.BuscoConfig_unittests.BuscoConfig.os.path.isdir", return_value=False
    )
    @patch(
        "__main__.BuscoConfig_unittests.BuscoConfig.os.path.isfile", return_value=True
    )
    def test_batch_mode_false_with_file(self, *args):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.set = Mock()
        config._check_batch_mode()

    @patch(
        "__main__.BuscoConfig_unittests.BuscoConfig.os.path.isdir", return_value=False
    )
    @patch(
        "__main__.BuscoConfig_unittests.BuscoConfig.os.path.isfile", return_value=False
    )
    def test_batch_mode_false_with_error(self, *args):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.set = Mock()
        with self.assertRaises(BuscoConfig.BatchFatalError):
            config._check_batch_mode()

    def test_required_input_exists_false(self):
        input_filename = "test_input_file"
        if os.path.exists(input_filename):
            os.remove(input_filename)
        self.test_params["in"] = input_filename
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        with self.assertRaises(BuscoConfig.BatchFatalError):
            config._check_required_input_exists()

    @patch("__main__.BuscoConfig_unittests.BuscoConfig.BuscoDownloadManager")
    def test_downloader_initialized(self, mock_downloader):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        config._init_downloader()
        mock_downloader.assert_called()

    @patch("__main__.BuscoConfig_unittests.BuscoConfig.PrettyLog")
    def test_log_config(self, mock_pretty_log):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        with self.assertLogs(BuscoConfig.logger, level="DEBUG"):
            config.log_config()
        mock_pretty_log.assert_called()

    @patch.object(BuscoConfig.BuscoConfigMain, "log_config")
    @patch.object(BuscoConfig.BuscoConfigMain, "_init_downloader")
    @patch.object(BuscoConfig.BuscoConfigMain, "_check_batch_mode")
    @patch.object(BuscoConfig.BuscoConfigMain, "_check_required_input_exists")
    @patch.object(BuscoConfig.BuscoConfigMain, "_expand_all_paths")
    @patch.object(BuscoConfig.BuscoConfigMain, "_check_out_value")
    @patch.object(BuscoConfig.BuscoConfigMain, "_check_allowed_keys")
    @patch.object(BuscoConfig.BuscoConfigMain, "_create_required_paths")
    @patch.object(BuscoConfig.BuscoConfigMain, "_check_no_previous_run")
    @patch.object(BuscoConfig.BuscoConfigMain, "_check_mandatory_keys_exist")
    def test_validation(
        self,
        mock_check_mandatory_keys,
        mock_check_no_previous_run,
        mock_create_required_paths,
        mock_check_allowed_keys,
        mock_check_out_value,
        mock_expand_all_paths,
        mock_check_input,
        mock_check_batch,
        mock_init_downloader,
        mock_log_config,
    ):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        config.validate()
        mock_check_mandatory_keys.assert_called()
        mock_check_no_previous_run.assert_called()
        mock_create_required_paths.assert_called()
        mock_check_allowed_keys.assert_called()
        mock_check_out_value.assert_called()
        mock_expand_all_paths.assert_called()
        mock_check_input.assert_called()
        mock_check_batch.assert_called()
        mock_init_downloader.assert_called()
        mock_log_config.assert_called()

    def test_check_lineage_present_false(self):
        try:
            del self.test_params["lineage_dataset"]  # just in case, probably redundant
        except KeyError:
            pass
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        self.assertFalse(config.check_lineage_present())

    def test_check_lineage_present_true_with_dataset_version_correct(self):
        self.test_params["lineage_dataset"] = "test_dataset_odb10"
        self.test_params["datasets_version"] = "odb10"
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        self.assertEqual(
            config.get("busco_run", "datasets_version"),
            self.test_params["datasets_version"],
        )

    def test_check_lineage_present_true_with_dataset_version_mismatch(self):
        self.test_params["lineage_dataset"] = "test_dataset_odb10"
        self.test_params["datasets_version"] = "odb11"
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        with self.assertLogs(BuscoConfig.logger, level="WARNING"):
            config.configure()
        self.assertEqual(
            config.get("busco_run", "datasets_version"),
            self.test_params["lineage_dataset"].split("_")[-1],
        )

    def test_check_lineage_present_true_with_odb_missing(self):
        self.test_params["lineage_dataset"] = "test_dataset"
        self.test_params["datasets_version"] = "odb10"
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        self.assertEqual(
            config.get("busco_run", "lineage_dataset"),
            "{}_{}".format(
                self.test_params["lineage_dataset"],
                self.test_params["datasets_version"],
            ),
        )

    def test_check_lineage_present_true_with_invalid_dataset_version(self):
        self.test_params["lineage_dataset"] = "test_dataset"
        self.test_params["datasets_version"] = "odb11"
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        with self.assertRaises(BuscoConfig.BatchFatalError):
            config.configure()

    def test_set_results_dirname(self):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        test_dataset_path = "/path/to/lineage_dataset"
        with patch("busco.BuscoConfig.BaseConfig.set"):
            config.set_results_dirname(test_dataset_path)
            config.set.assert_called_with(
                "busco_run", "lineage_results_dir", "run_lineage_dataset"
            )

    @patch("busco.BuscoConfig.BuscoConfigAuto.add_mode_specific_parameters")
    @patch("busco.BuscoConfig.BuscoConfigAuto.load_dataset_config")
    @patch("busco.BuscoConfig.BuscoConfigAuto.download_lineage_file")
    @patch("busco.BuscoConfig.BuscoConfigAuto._create_required_paths")
    @patch("busco.BuscoConfig.BuscoConfigAuto.set_results_dirname")
    @patch("busco.BuscoConfig.BaseConfig")
    @patch("busco.BuscoConfig.BuscoConfigAuto.get", return_value="test")
    @patch("busco.BuscoConfig.BuscoConfigAuto._propagate_config")
    def test_autoconfig_init_propagates_mainconfig(self, mock_propagate, *args):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        BuscoConfig.BuscoConfigAuto(config, None)
        mock_propagate.assert_called_with(config)

    @patch("busco.BuscoConfig.BuscoConfigAuto.add_mode_specific_parameters")
    @patch("busco.BuscoConfig.BuscoConfigAuto.load_dataset_config")
    @patch("busco.BuscoConfig.BuscoConfigAuto.download_lineage_file")
    @patch("busco.BuscoConfig.BuscoConfigAuto._create_required_paths")
    @patch("busco.BuscoConfig.BaseConfig")
    @patch("busco.BuscoConfig.BuscoConfigAuto._propagate_config")
    @patch("busco.BuscoConfig.BuscoConfigAuto.get", return_value="test")
    @patch("busco.BuscoConfig.BuscoConfigAuto.set_results_dirname")
    def test_autoconfig_init_sets_results_dirname(self, mock_set_dirname, *args):
        BuscoConfig.BuscoConfigAuto(None, "lineage")
        mock_set_dirname.assert_called_with("lineage")

    @patch("busco.BuscoConfig.BuscoConfigAuto.load_dataset")
    @patch("busco.BuscoConfig.BaseConfig")
    @patch("busco.BuscoConfig.BuscoConfigAuto.add_mode_specific_parameters")
    @patch("busco.BuscoConfig.BuscoConfigAuto._propagate_config")
    @patch("busco.BuscoConfig.BuscoConfigAuto.get", return_value="test")
    @patch("busco.BuscoConfig.BuscoConfigAuto._create_required_paths")
    def test_autoconfig_init_creates_paths(self, mock_create_paths, *args):
        BuscoConfig.BuscoConfigAuto(None, None)
        mock_create_paths.assert_called()

    @patch("busco.BuscoConfig.BuscoConfigAuto.load_dataset_config")
    @patch("busco.BuscoConfig.BaseConfig")
    @patch("busco.BuscoConfig.BuscoConfigAuto.add_mode_specific_parameters")
    @patch("busco.BuscoConfig.BuscoConfigAuto._propagate_config")
    @patch("busco.BuscoConfig.BuscoConfigAuto.set_results_dirname")
    @patch("busco.BuscoConfig.BuscoConfigAuto._create_required_paths")
    @patch("busco.BuscoConfig.BuscoConfigAuto.get", return_value="test")
    @patch("busco.BuscoConfig.BuscoConfigAuto.download_lineage_file")
    def test_autoconfig_init_downloads_lineage(self, mock_download_lineage, *args):
        BuscoConfig.BuscoConfigAuto(None, "lineage")
        mock_download_lineage.assert_called_with("lineage")

    @patch("busco.BuscoConfig.BuscoConfigAuto.load_dataset_config")
    @patch("busco.BuscoConfig.BuscoConfigAuto.download_lineage_file")
    @patch("busco.BuscoConfig.BaseConfig")
    @patch("busco.BuscoConfig.BuscoConfigAuto._propagate_config")
    @patch("busco.BuscoConfig.BuscoConfigAuto.set_results_dirname")
    @patch("busco.BuscoConfig.BuscoConfigAuto._create_required_paths")
    @patch("busco.BuscoConfig.BuscoConfigAuto.get", return_value="test")
    @patch("busco.BuscoConfig.BuscoConfigAuto.add_mode_specific_parameters")
    def test_autoconfig_init_adds_specific_parameters(self, mock_add_parameters, *args):
        BuscoConfig.BuscoConfigAuto(None, "lineage")
        mock_add_parameters.assert_called()

    @patch("busco.BuscoConfig.BaseConfig")
    @patch("busco.BuscoConfig.BuscoConfigAuto.add_mode_specific_parameters")
    @patch("busco.BuscoConfig.BuscoConfigAuto._propagate_config")
    @patch("busco.BuscoConfig.BuscoConfigAuto.set_results_dirname")
    @patch("busco.BuscoConfig.BuscoConfigAuto._create_required_paths")
    @patch("busco.BuscoConfig.BuscoConfigAuto.download_lineage_file")
    @patch("busco.BuscoConfig.BuscoConfigAuto.get", return_value="test")
    @patch("busco.BuscoConfig.BuscoConfigAuto.load_dataset_config")
    def test_autoconfig_init_loads_lineage_config(self, mock_load_dataset, *args):
        BuscoConfig.BuscoConfigAuto(None, None)
        mock_load_dataset.assert_called()

    @patch("busco.BuscoConfig.BuscoConfigAuto._propagate_config")
    @patch("busco.BuscoConfig.BuscoConfigAuto._create_required_paths")
    @patch("busco.BuscoConfig.BuscoConfigAuto.load_dataset")
    @patch("busco.BuscoConfig.BuscoConfigAuto.get", return_value="test")
    @patch("busco.BuscoConfig.BaseConfig.__init__")
    def test_autoconfig_init_calls_super(self, mock_config_parent, *args):
        BuscoConfig.BuscoConfigAuto(None, None)
        mock_config_parent.assert_called()

    @patch("busco.BuscoConfig.BuscoConfigAuto.add_mode_specific_parameters")
    @patch("busco.BuscoConfig.BuscoConfigAuto._create_required_paths")
    @patch("busco.BuscoConfig.BuscoConfigAuto.download_lineage_file")
    @patch("busco.BuscoConfig.BuscoConfigAuto.load_dataset_config")
    def test_propagate_config(self, *args):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.params)
        config.configure()
        config.downloader = Mock()
        autoconfig = BuscoConfig.BuscoConfigAuto(config, "test")
        autoconfig._propagate_config(config)
        self.assertEqual(autoconfig, config)

    @patch("busco.BuscoConfig.BuscoConfigAuto.load_dataset_config")
    @patch("busco.BuscoConfig.BuscoConfigAuto.download_lineage_file")
    @patch("busco.BuscoConfig.BuscoConfigAuto._propagate_config")
    @patch("busco.BuscoConfig.BuscoConfigAuto.set_results_dirname")
    @patch("busco.BuscoConfig.BuscoConfigAuto.get", return_value="test")
    @patch("busco.BuscoConfig.BaseConfig._create_required_paths")
    def test_autolineage_create_path_method_calls_parent(
        self, mock_create_paths, *args
    ):
        config = BuscoConfig.BuscoConfigMain(self.base_config, self.test_params)
        config.configure()
        BuscoConfig.BuscoConfigAuto(config, None)
        mock_create_paths.assert_called_with("test/auto_lineage")

    def tearDown(self):
        self.test_params = {}
