#!/usr/bin/env python3
# coding: utf-8
"""
setup.py

Script for package installation with setuptools.

Copyright (c) 2015-2024, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

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
