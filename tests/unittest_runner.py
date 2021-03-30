import unittest
from tests.unittests import run_BUSCO_unittests
from tests.unittests import ConfigManager_unittests
from tests.unittests import BuscoConfig_unittests
from tests.unittests import AutoLineage_unittests
from tests.unittests import GenomeAnalysis_unittests
from tests.unittests import BuscoRunner_unittests
import sys

loader = unittest.TestLoader()
suite = unittest.TestSuite()

suite.addTests(loader.loadTestsFromModule(run_BUSCO_unittests))
suite.addTests(loader.loadTestsFromModule(ConfigManager_unittests))
suite.addTests(loader.loadTestsFromModule(BuscoConfig_unittests))
suite.addTests(loader.loadTestsFromModule(AutoLineage_unittests))
suite.addTests(loader.loadTestsFromModule(GenomeAnalysis_unittests))
suite.addTests(loader.loadTestsFromModule(BuscoRunner_unittests))

runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)

ret = not result.wasSuccessful()
sys.exit(ret)
