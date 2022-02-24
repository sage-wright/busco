#!/usr/bin/env python3
# coding: utf-8
"""
.. versionadded:: 3.0.0
.. versionchanged:: 5.0.0

Copyright (c) 2016-2022, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

This script proceeds to the BUSCO packages installation

"""

from distutils.core import setup

version = {}
with open("src/busco/_version.py") as version_file:
    exec(version_file.read(), version)

setup(
    name="BUSCO",
    version=version["__version__"],
    author="ezlab",
    license="Licensed under the MIT license. See LICENSE.md file.",
    author_email="ez@ezlab.org",
    long_description="Assessing genome assembly and annotation completeness "
    "with Benchmarking Universal Single-Copy Orthologs ",
    url="https://busco.ezlab.org/",
    platforms="Unix like",
    packages=["busco", "busco.analysis", "busco.busco_tools"],
    package_dir={"busco": "src/busco"},
    scripts=["bin/busco"],
)
