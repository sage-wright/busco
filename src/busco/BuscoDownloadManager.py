#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: BuscoDownloadManager
   :synopsis: BuscoDownloadManager manage the version and download the most recent file
.. versionadded:: 4.0.0
.. versionchanged:: 5.4.0

Copyright (c) 2016-2022, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""

import os
import time
import glob
import tarfile
import hashlib
import urllib.request
from urllib.error import URLError
import gzip
import pandas as pd

from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.Exceptions import BatchFatalError, BuscoError

logger = BuscoLogger.get_logger(__name__)


class BuscoDownloadManager:
    """
    This class obtains and manages the version of data files required to run a BUSCO analysis.
    When the config parameter `offline` is set, no attempt to download is made.
    When the config parameter `auto_update_files` is set, new versions replace old versions.s
    Else, a warning is produced.
    """

    version_files = None

    def __init__(self, config):
        """
        :param config: Values of all parameters to be used during the analysis
        :type config: BuscoConfig
        """
        self.offline = config.getboolean("busco_run", "offline")
        self.update_data = config.getboolean("busco_run", "update-data")
        self.download_base_url = config.get("busco_run", "download_base_url")
        self.local_download_path = config.get("busco_run", "download_path")
        self._create_main_download_dir()
        if not type(self).version_files is not None and not self.offline:
            try:
                self._obtain_versions_file()
            except BatchFatalError as bfe:
                self.download_base_url = "https://busco-data2.ezlab.org/v5/data/"
                logger.info(
                    "{}. Retrying with backup url {}".format(
                        bfe.value, self.download_base_url
                    )
                )
                self._obtain_versions_file()

    def _create_main_download_dir(self):
        if not os.path.exists(self.local_download_path):
            # exist_ok=True to allow for multiple parallel BUSCO runs each trying to create this folder simultaneously
            os.makedirs(self.local_download_path, exist_ok=True)

    @log("Downloading information on latest versions of BUSCO data...", logger)
    def _obtain_versions_file(self):
        remote_filepath = os.path.join(self.download_base_url, "file_versions.tsv")
        local_filepath = os.path.join(self.local_download_path, "file_versions.tsv")

        for tsleep in [10, 100, None]:
            try:
                urllib.request.urlretrieve(remote_filepath, local_filepath)
                type(self).version_files = pd.read_csv(
                    local_filepath,
                    sep="\t",
                    names=["dataset", "date", "hash", "domain", "type"],
                    index_col="dataset",
                )
            except ValueError:
                raise BatchFatalError(
                    "Invalid URL. Cannot reach {}. Please provide a valid URL for --base_download_url".format(
                        remote_filepath
                    )
                )
            except URLError:
                if self.offline:
                    logger.warning(
                        "Unable to verify BUSCO datasets because of offline mode"
                    )
                    break
                elif tsleep:
                    logger.info(
                        "Download connection problem. Retrying in {} seconds".format(
                            tsleep
                        )
                    )
                    time.sleep(tsleep)
                else:
                    raise BatchFatalError("Cannot reach {}".format(remote_filepath))

        return

    def _create_category_dir(self, category):
        # if the category folder does not exist, create it
        category_folder = os.path.join(self.local_download_path, category)
        if not os.path.exists(category_folder):
            # exist_ok=True to allow for multiple parallel BUSCO runs, each trying to create this folder
            os.makedirs(category_folder, exist_ok=True)
        return

    @staticmethod
    def _extract_creation_date(dataset_config_file):
        dataset_date = None
        with open(dataset_config_file, "r") as data_config:
            for line in data_config:
                line = line.strip().split("=")
                if line[0] == "creation_date":
                    dataset_date = line[1]
                    break
        if not dataset_date:
            raise BuscoError(
                "Creation date could not be extracted from dataset.cfg file."
            )
        return dataset_date

    def _check_existing_version(self, local_filepath, category, data_basename):
        try:
            latest_update = type(self).version_files.loc[data_basename]["date"]
        except KeyError:
            raise BuscoError(
                "{} is not a valid option for '{}'".format(data_basename, category)
            )
        path_basename, extension = os.path.splitext(data_basename)

        if category == "lineages":
            latest_version = ".".join([path_basename, latest_update])
            try:
                dataset_date = self._extract_creation_date(
                    os.path.join(local_filepath, "dataset.cfg")
                )
                up_to_date = dataset_date == latest_update
                present = True
            except FileNotFoundError:
                up_to_date = False
                present = False

        else:
            latest_version = ".".join(
                [path_basename, latest_update, extension.lstrip(".")]
            )
            local_filepath = local_filepath.replace(data_basename, latest_version)
            up_to_date = os.path.exists(local_filepath)
            path_to_check, extension_to_check = os.path.splitext(local_filepath)
            present = (
                len(
                    glob.glob(
                        "{}.*.{}".format(path_to_check[0:-11], extension_to_check[1:])
                    )
                )
                > 0
            )

        hash = type(self).version_files.loc[data_basename]["hash"]

        return present, up_to_date, latest_version, local_filepath, hash

    def get(self, data_name, category):
        if os.path.exists(data_name) and "/" in data_name:
            logger.info("Using local {} directory {}".format(category, data_name))
            return data_name
        elif "/" in data_name:
            raise BuscoError("{} does not exist".format(data_name))
        if self.offline:
            if category == "lineages":
                local_dataset = os.path.join(
                    self.local_download_path, category, data_name
                )
                if os.path.exists(local_dataset):
                    return local_dataset
                else:
                    raise BuscoError(
                        "Unable to run BUSCO in offline mode. Dataset {} does not "
                        "exist.".format(local_dataset)
                    )
            else:
                basename, extension = os.path.splitext(data_name)
                placement_files = sorted(
                    glob.glob(
                        os.path.join(
                            self.local_download_path,
                            category,
                            "{}.*{}".format(basename, extension),
                        )
                    )
                )
                if len(placement_files) > 0:
                    return placement_files[-1]
                else:
                    raise BuscoError(
                        "Unable to run BUSCO placer in offline mode. Cannot find necessary placement "
                        "files in {}".format(self.local_download_path)
                    )
        data_basename = os.path.basename(data_name)
        local_filepath = os.path.join(self.local_download_path, category, data_basename)
        (
            present,
            up_to_date,
            latest_version,
            local_filepath,
            hash,
        ) = self._check_existing_version(local_filepath, category, data_basename)

        if (not up_to_date and self.update_data) or not present:
            # download
            self._create_category_dir(category)
            compression_extension = ".tar.gz"
            remote_filepath = os.path.join(
                self.download_base_url, category, latest_version + compression_extension
            )
            if present and category == "lineages":
                self._rename_old_version(local_filepath)

            for tsleep in [10, 100, None]:
                download_success = self._download_file(
                    remote_filepath, local_filepath + compression_extension, hash
                )
                if download_success:
                    local_filepath = self._decompress_file(
                        local_filepath + compression_extension
                    )
                    if present:
                        logger.warning(
                            "The file or folder {} was updated automatically.".format(
                                data_basename
                            )
                        )
                    break
                else:
                    if tsleep:
                        logger.info(
                            "Download connection problem. Retrying in {} seconds".format(
                                tsleep
                            )
                        )
                        time.sleep(tsleep)
            else:
                raise BuscoError(
                    "File download unsuccessful. Check the download url or consider running in offline mode"
                )
        elif not up_to_date:
            logger.warning(
                "The file or folder {} is not the last available version. "
                "To update all data files to the last version, add the parameter "
                "--update-data in your next run.".format(local_filepath)
            )

        return local_filepath

    @staticmethod
    def _rename_old_version(local_filepath):
        if os.path.exists(local_filepath):
            try:
                os.rename(local_filepath, "{}.old".format(local_filepath))
                logger.info(
                    "Renaming {} into {}.old".format(local_filepath, local_filepath)
                )
            except OSError:
                try:
                    timestamp = time.time()
                    os.rename(
                        local_filepath, "{}.old.{}".format(local_filepath, timestamp)
                    )
                    logger.info(
                        "Renaming {} into {}.old.{}".format(
                            local_filepath, local_filepath, timestamp
                        )
                    )
                except OSError:
                    pass
        return

    @log("Downloading file {}", logger, func_arg=1)
    def _download_file(self, remote_filepath, local_filepath, expected_hash):
        try:
            urllib.request.urlretrieve(remote_filepath, local_filepath)
            observed_hash = type(self)._md5(local_filepath)
            if observed_hash != expected_hash:
                logger.error(
                    "md5 hash is incorrect: {} while {} expected".format(
                        str(observed_hash), str(expected_hash)
                    )
                )
                logger.info("deleting corrupted file {}".format(local_filepath))
                os.remove(local_filepath)
                raise BuscoError(
                    "BUSCO was unable to download or update all necessary files"
                )
            else:
                logger.debug("md5 hash is {}".format(observed_hash))
        except URLError:
            logger.error("Cannot reach {}".format(remote_filepath))
            return False
        return True

    @staticmethod
    def _md5(fname):
        hash_md5 = hashlib.md5()
        with open(fname, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    @log("Decompressing file {}", logger, func_arg=1)
    def _decompress_file(self, local_filepath):
        unzipped_filename = local_filepath.replace(".gz", "")

        if os.path.splitext(local_filepath)[1] == ".gz":
            with gzip.open(os.path.join(local_filepath), "rb") as compressed_file:
                with open(unzipped_filename, "wb") as decompressed_file:
                    for line in compressed_file:
                        decompressed_file.write(line)
            try:
                os.remove(local_filepath)
            except FileNotFoundError:
                logger.warning(
                    "Downloaded gzip file was removed before this BUSCO run could remove it."
                )
            local_filepath = unzipped_filename

        if os.path.splitext(local_filepath)[1] == ".tar":
            untarred_filename = local_filepath.replace(".tar", "")
            with tarfile.open(local_filepath) as tar_file:
                tar_file.extractall(os.path.dirname(local_filepath))
            try:
                os.remove(local_filepath)
            except FileNotFoundError:
                logger.warning(
                    "Downloaded tarball was removed before this BUSCO run could remove it."
                )
            local_filepath = untarred_filename

        return local_filepath
