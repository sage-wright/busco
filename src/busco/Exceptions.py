# coding: utf-8
"""
Exceptions.py

Custom exceptions for BUSCO.

Author(s): Matthew Berkeley

Copyright (c) 2015-2024, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""


class BatchFatalError(Exception):
    """
    Error that prevents batch mode from running.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


class BuscoError(Exception):
    """
    Module-specific exception
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value
