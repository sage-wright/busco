import unittest
from busco.analysis import GenomeAnalysis
from unittest.mock import patch, Mock


class TestConfigManager(unittest.TestCase):
    def setUp(self) -> None:
        pass

    # @patch('busco.analysis.GenomeAnalysis.BuscoAnalysis.config.get', return_value="test")
    # @patch('busco.analysis.GenomeAnalysis.BuscoAnalysis.config.getboolean', return_value=True)
    @patch("busco.analysis.GenomeAnalysis.BuscoAnalysis.config")
    @patch("busco.analysis.BuscoAnalysis.os.path")
    @patch("busco.analysis.GenomeAnalysis.NucleotideAnalysis.check_nucleotide_file")
    def test_init_eukaryota_augustus_checks_filetype(self, mock_check_nucl_file, *args):
        GenomeAnalysis.GenomeAnalysisEukaryotesAugustus()
        mock_check_nucl_file.assert_called()

    @patch("busco.analysis.GenomeAnalysis.NucleotideAnalysis.__init__")
    @patch(
        "busco.analysis.GenomeAnalysis.BuscoAnalysis.config.get",
        return_value="euk_genome_aug",
    )
    # @patch('busco.analysis.GenomeAnalysis.BuscoAnalysis.config.getboolean')
    @patch("busco.analysis.GenomeAnalysis.BuscoAnalysis.config")
    @patch("busco.analysis.Analysis.logger.warning")
    @patch("busco.analysis.GenomeAnalysis.BBToolsRunner")
    @patch("busco.analysis.GenomeAnalysis.OptimizeAugustusRunner")
    @patch("busco.analysis.GenomeAnalysis.ETrainingRunner")
    @patch("busco.analysis.GenomeAnalysis.NewSpeciesRunner")
    @patch("busco.analysis.GenomeAnalysis.GFF2GBRunner")
    @patch("busco.analysis.GenomeAnalysis.AugustusRunner")
    @patch("busco.analysis.Analysis.TBLASTNRunner")
    @patch("busco.analysis.Analysis.MKBLASTRunner")
    @patch("busco.analysis.BuscoAnalysis.HMMERRunner")
    def test_init_tools_eukaryota_augustus(
        self,
        mock_hmmer,
        mock_mkblast,
        mock_tblastn,
        mock_augustus,
        mock_gff2gb,
        mock_new_species,
        mock_etraining,
        mock_optimize_augustus,
        mock_bbtools,
        *args
    ):
        analysis = GenomeAnalysis.GenomeAnalysisEukaryotesAugustus()
        analysis.init_tools()
        mock_hmmer.assert_called()
        mock_mkblast.assert_called()
        mock_tblastn.assert_called()
        mock_augustus.assert_called()
        mock_gff2gb.assert_called()
        mock_new_species.assert_called()
        mock_etraining.assert_called()
        mock_optimize_augustus.assert_called()
        mock_bbtools.assert_called()

    @patch("busco.analysis.GenomeAnalysis.NucleotideAnalysis.__init__")
    @patch(
        "busco.analysis.GenomeAnalysis.BuscoAnalysis.config.get",
        return_value="euk_genome_met",
    )
    @patch("busco.analysis.GenomeAnalysis.BuscoAnalysis.config", autospec=True)
    @patch("busco.analysis.GenomeAnalysis.BBToolsRunner")
    @patch("busco.analysis.GenomeAnalysis.MetaeukRunner")
    @patch("busco.analysis.BuscoAnalysis.HMMERRunner")
    def test_init_tools_eukaryota_metaeuk(
        self, mock_hmmer, mock_metaeuk, mock_bbtools, *args
    ):
        analysis = GenomeAnalysis.GenomeAnalysisEukaryotesMetaeuk()
        analysis.init_tools()
        mock_hmmer.assert_called()
        mock_metaeuk.assert_called()
        mock_bbtools.assert_called()

    @patch("busco.analysis.GenomeAnalysis.NucleotideAnalysis.__init__")
    @patch(
        "busco.analysis.GenomeAnalysis.BuscoAnalysis.config.get",
        return_value="prok_genome",
    )
    @patch("busco.analysis.GenomeAnalysis.BuscoAnalysis.config", autospec=True)
    @patch("busco.analysis.GenomeAnalysis.BBToolsRunner")
    @patch("busco.analysis.GenomeAnalysis.ProdigalRunner")
    @patch("busco.analysis.BuscoAnalysis.HMMERRunner")
    def test_init_tools_prokaryota(
        self, mock_hmmer, mock_prodigal, mock_bbtools, *args
    ):
        analysis = GenomeAnalysis.GenomeAnalysisProkaryotes()
        analysis.init_tools()
        mock_hmmer.assert_called()
        mock_prodigal.assert_called()
        mock_bbtools.assert_called()

    @patch("busco.analysis.GenomeAnalysis.NucleotideAnalysis.__init__")
    @patch(
        "busco.analysis.GenomeAnalysis.BuscoAnalysis.config.get",
        return_value="euk_genome_met",
    )
    @patch("busco.analysis.GenomeAnalysis.BuscoAnalysis.config", autospec=True)
    @patch("busco.analysis.GenomeAnalysis.BuscoAnalysis.run_analysis")
    @patch("busco.analysis.GenomeAnalysis.GenomeAnalysis._run_bbtools")
    @patch("busco.analysis.GenomeAnalysis.BuscoAnalysis.run_hmmer")
    @patch("busco.analysis.GenomeAnalysis.GenomeAnalysisEukaryotesMetaeuk._run_metaeuk")
    def test_run_analysis_metaeuk(
        self, mock_run_metaeuk, mock_run_hmmer, mock_run_bbtools, *args
    ):
        analysis = GenomeAnalysis.GenomeAnalysisEukaryotesMetaeuk()
        analysis.bbtools_runner = Mock()
        analysis.metaeuk_runner = Mock()
        analysis.hmmer_runner = Mock(missing_buscos=[])
        analysis.hmmer_runner.fragmented_buscos.keys = Mock(return_value=[])
        analysis.gene_details = Mock(autospec=True)
        analysis.sequences_aa = Mock(autospec=True)
        analysis.sequences_nt = Mock(autospec=True)
        analysis.write_gff_files = Mock()
        analysis.run_analysis()
        mock_run_metaeuk.assert_called()
        mock_run_hmmer.assert_called()
        mock_run_bbtools.assert_called()
        analysis.write_gff_files.assert_called()

        # @patch('busco.GenomeAnalysis.GenomeAnalysisEukaryotesAugustus._rerun_analysis')
        # @patch('busco.GenomeAnalysis.GenomeAnalysisEukaryotesAugustus.run_hmmer')
        # @patch('busco.GenomeAnalysis.GenomeAnalysisEukaryotesAugustus._run_augustus')
        # @patch('busco.GenomeAnalysis.BLASTAnalysis._run_tblastn')
        # @patch('busco.GenomeAnalysis.BLASTAnalysis._run_mkblast')
        # mock_mkblast.assert_called()
        # mock_tblastn.assert_called()
        # mock_augustus.assert_called()
        # mock_hmmer.assert_called()
        # mock_rerun.assert_called()

    def tearDown(self) -> None:
        pass
