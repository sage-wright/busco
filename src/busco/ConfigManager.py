#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: ConfigManager
   :synopsis: manage setup for BUSCO run configuration
.. versionadded:: 3.0.0
.. versionchanged:: 5.4.0

Copyright (c) 2016-2022, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""

from busco.BuscoConfig import BuscoConfigMain, BuscoConfigAuto
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
import os

logger = BuscoLogger.get_logger(__name__)


class BuscoConfigManager:
    def __init__(self, params):
        self.params = params
        self.config_file = None
        self.config_main = None
        self.get_config_file()
        self.runner = None

    @log("Getting config file", logger, debug=True)
    def get_config_file(self):
        """
        Check for BUSCO config file specified as a command line argument;
        if not present check if defined as an environment variable;
        if not present use default config file.
        :return config: A BuscoConfig object containing all the required configuration parameters
        """
        try:
            self.config_file = self.params["config_file"]
            if self.config_file is not None:
                return
        except KeyError:
            pass
        if os.environ.get("BUSCO_CONFIG_FILE") and os.access(
            os.environ.get("BUSCO_CONFIG_FILE"), os.R_OK
        ):
            self.config_file = os.environ.get("BUSCO_CONFIG_FILE")
        else:
            self.config_file = "local environment"
        return self.config_file

    @log("Configuring BUSCO with {}", logger, attr_name="config_file")
    def load_busco_config_main(self):
        self.config_main = BuscoConfigMain(self.config_file, self.params)
        self.config_main.configure()
        self.config_main.validate()
        return

    def load_busco_config_auto(self, lineage):
        autoconfig = BuscoConfigAuto(self.config_main, lineage)
        return autoconfig
