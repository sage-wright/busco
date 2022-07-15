import unittest
from busco import run_BUSCO
import sys
import io
from unittest.mock import Mock, patch, call


class TestParams(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        pass

    def test_help_short(self):
        args = ["-h"]
        sys.argv[1:] = args
        with self.assertRaises(SystemExit) as cm:
            captured_output = io.StringIO()
            sys.stdout = captured_output
            try:
                run_BUSCO._parse_args()
            finally:
                sys.stdout = sys.__stdout__
        self.assertEqual(cm.exception.code, 0)

    def test_help_long(self):
        args = ["--help"]
        sys.argv[1:] = args
        with self.assertRaises(SystemExit) as cm:
            captured_output = io.StringIO()
            sys.stdout = captured_output
            try:
                run_BUSCO._parse_args()
            finally:
                sys.stdout = sys.__stdout__
        self.assertEqual(cm.exception.code, 0)

    def test_version_short(self):
        args = ["-v"]
        sys.argv[1:] = args
        with self.assertRaises(SystemExit) as cm:
            captured_output = io.StringIO()
            sys.stdout = captured_output
            try:
                run_BUSCO._parse_args()
            finally:
                sys.stdout = sys.__stdout__
        self.assertEqual(cm.exception.code, 0)

    def test_version_long(self):
        args = ["--version"]
        sys.argv[1:] = args
        with self.assertRaises(SystemExit) as cm:
            captured_output = io.StringIO()
            sys.stdout = captured_output
            try:
                run_BUSCO._parse_args()
            finally:
                sys.stdout = sys.__stdout__
        self.assertEqual(cm.exception.code, 0)

    def test_list_lineages(self):
        args = ["--list-datasets"]
        sys.argv[1:] = args
        with self.assertRaises(SystemExit) as cm:
            captured_output = io.StringIO()
            sys.stdout = captured_output
            try:
                run_BUSCO._parse_args()
            finally:
                sys.stdout = sys.__stdout__
        self.assertEqual(cm.exception.code, 0)

    @patch("busco.BuscoDownloadManager.logger.info")
    def test_direct_download(self, *args):
        args = ["--download", "archaea_odb10"]
        sys.argv[1:] = args
        with self.assertRaises(SystemExit) as cm:
            captured_output = io.StringIO()
            sys.stdout = captured_output
            try:
                run_BUSCO._parse_args()
            finally:
                sys.stdout = sys.__stdout__
        self.assertEqual(cm.exception.code, 0)

    def test_cmdline_options_short_minimum(self):
        params = run_BUSCO._parse_args()
        correct_parse = {
            "auto-lineage": False,
            "auto-lineage-euk": False,
            "auto-lineage-prok": False,
            "config_file": None,
            "contig_break": None,
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
            "download": "==SUPPRESS==",
            "download_base_url": None,
            "download_path": None,
            "update-data": False,
            "version": "==SUPPRESS==",
            "tar": False,
            "scaffold_composition": False,
        }
        self.assertDictEqual(params, correct_parse)

    def test_cmdline_options_all_short(self):
        input_file = "input_file"
        output_file = "output_file"
        mode = "mode"
        lineage_dataset = "lineage_dataset"
        cpus = 10
        evalue = 0.1

        arg_values = {
            "-i": input_file,
            "-o": output_file,
            "-m": mode,
            "-l": lineage_dataset,
            "-c": cpus,
            "-e": evalue,
        }
        flag_options = ["-f", "-q", "-r"]
        command_str = " ".join(
            [" ".join([key, str(value)]) for key, value in arg_values.items()]
            + flag_options
        )
        sys.argv[1:] = command_str.split(" ")
        params = run_BUSCO._parse_args()
        correct_parse = {
            "auto-lineage": False,
            "auto-lineage-euk": False,
            "auto-lineage-prok": False,
            "config_file": None,
            "contig_break": None,
            "cpu": cpus,
            "evalue": evalue,
            "force": True,
            "help": "==SUPPRESS==",
            "in": input_file,
            "limit": None,
            "lineage_dataset": lineage_dataset,
            "list_datasets": "==SUPPRESS==",
            "mode": mode,
            "offline": False,
            "out": output_file,
            "out_path": None,
            "quiet": True,
            "restart": True,
            "scaffold_composition": False,
            "metaeuk_parameters": None,
            "metaeuk_rerun_parameters": None,
            "use_augustus": False,
            "augustus_parameters": None,
            "augustus_species": None,
            "long": False,
            "datasets_version": None,
            "download": "==SUPPRESS==",
            "download_base_url": None,
            "download_path": None,
            "update-data": False,
            "version": "==SUPPRESS==",
            "tar": False,
        }
        self.assertDictEqual(params, correct_parse)

    def test_cmdline_options_all_long(self):
        input_file = "input_file"
        output_file = "output_file"
        mode = "mode"
        lineage_dataset = "lineage_dataset"
        cpus = 10
        evalue = 0.1
        limit = 1
        contig_break = 5
        augustus_parameters = "augustus_parameters"
        augustus_species = "augustus_species"
        config = "config"
        out_path = "out_path"
        download_path = "download_path"
        datasets_version = "datasets_version"
        download_base_url = "download_base_url"
        metaeuk_parameters = "metaeuk_parameters"
        metaeuk_rerun_parameters = "metaeuk_rerun_parameters"

        arg_values = {
            "--in": input_file,
            "--cpu": cpus,
            "--out": output_file,
            "--evalue": evalue,
            "--mode": mode,
            "--lineage_dataset": lineage_dataset,
            "--limit": limit,
            "--augustus_parameters": augustus_parameters,
            "--augustus_species": augustus_species,
            "--config": config,
            "--contig_break": contig_break,
            "--out_path": out_path,
            "--download_path": download_path,
            "--datasets_version": datasets_version,
            "--download_base_url": download_base_url,
            "--metaeuk_parameters": metaeuk_parameters,
            "--metaeuk_rerun_parameters": metaeuk_rerun_parameters,
        }
        flag_options = [
            "--force",
            "--restart",
            "--quiet",
            "--long",
            "--auto-lineage",
            "--auto-lineage-prok",
            "--auto-lineage-euk",
            "--augustus",
            "--update-data",
            "--offline",
            "--tar",
            "--scaffold_composition",
        ]
        command_str = " ".join(
            [" ".join([key, str(value)]) for key, value in arg_values.items()]
            + flag_options
        )
        sys.argv[1:] = command_str.split(" ")
        params = run_BUSCO._parse_args()
        correct_parse = {
            "augustus_parameters": augustus_parameters,
            "augustus_species": augustus_species,
            "metaeuk_parameters": metaeuk_parameters,
            "metaeuk_rerun_parameters": metaeuk_rerun_parameters,
            "datasets_version": datasets_version,
            "download": "==SUPPRESS==",
            "download_base_url": download_base_url,
            "download_path": download_path,
            "auto-lineage": True,
            "auto-lineage-euk": True,
            "auto-lineage-prok": True,
            "config_file": config,
            "contig_break": contig_break,
            "cpu": cpus,
            "evalue": evalue,
            "force": True,
            "restart": True,
            "use_augustus": True,
            "help": "==SUPPRESS==",
            "in": input_file,
            "limit": limit,
            "lineage_dataset": lineage_dataset,
            "list_datasets": "==SUPPRESS==",
            "long": True,
            "mode": mode,
            "offline": True,
            "out": output_file,
            "out_path": out_path,
            "quiet": True,
            "update-data": True,
            "version": "==SUPPRESS==",
            "tar": True,
            "scaffold_composition": True,
        }
        self.assertDictEqual(params, correct_parse)

    def test_command_line_cpu_type(self):
        bad_args_cpu = ["cpus", "1.5", None]

        for arg in bad_args_cpu:
            sys.argv[1:] = ["--cpu", arg]
            with self.assertRaises(SystemExit) as cm:
                captured_output = io.StringIO()
                sys.stderr = captured_output
                try:
                    run_BUSCO._parse_args()
                finally:
                    sys.stderr = sys.__stderr__
            self.assertEqual(cm.exception.code, 2)

    def test_command_line_evalue_type(self):
        bad_args_evalue = ["evalue", None]

        for arg in bad_args_evalue:
            sys.argv[1:] = ["--evalue", arg]
            with self.assertRaises(SystemExit) as cm:
                captured_output = io.StringIO()
                sys.stderr = captured_output
                try:
                    run_BUSCO._parse_args()
                finally:
                    sys.stderr = sys.__stderr__
            self.assertEqual(cm.exception.code, 2)

    def test_command_line_limit_type(self):
        bad_args_limit = ["limit", "1.5", None]
        for arg in bad_args_limit:
            sys.argv[1:] = ["--limit", arg]
            with self.assertRaises(SystemExit) as cm:
                captured_output = io.StringIO()
                sys.stderr = captured_output
                try:
                    run_BUSCO._parse_args()
                finally:
                    sys.stderr = sys.__stderr__
            self.assertEqual(cm.exception.code, 2)

    def tearDown(self):
        sys.argv = [sys.argv[0]]


class TestMaster(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.params = {}
        pass

    def tearDown(self):
        pass
