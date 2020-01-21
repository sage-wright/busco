#!/usr/bin/env python
# coding: utf-8

"""
.. module:: BuscoPlacer
   :synopsis: BuscoPlacer implements methods required for automatically selecting the appropriate dataset
   to be used during BUSCO analysis
.. versionadded:: 4.0.0
.. versionchanged:: 4.0.0

Copyright (c) 2016-2020, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
import glob
import json
import os
from busco.Toolset import Tool
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
from busco.BuscoDownloadManager import BuscoDownloadManager
from Bio import SeqIO
from busco.BuscoTools import SEPPRunner

logger = BuscoLogger.get_logger(__name__)


class BuscoPlacer:

    _logger = BuscoLogger.get_logger(__name__)

    @log("***** Searching tree for chosen lineage to find best taxonomic match *****\n", logger)
    def __init__(self, config, run_folder, protein_seqs, single_copy_buscos):
        self._config = config
        self._params = config
        self.mode = self._config.get("busco_run", "mode")
        self.cpus = self._config.get("busco_run", "cpu")
        self.run_folder = run_folder
        self.placement_folder = os.path.join(run_folder, "placement_files")
        os.mkdir(self.placement_folder)
        self.downloader = self._config.downloader
        self.datasets_version = self._config.get("busco_run", "datasets_version")
        self.protein_seqs = protein_seqs
        self.single_copy_buscos = single_copy_buscos  # dict
        self._init_tools()

    def _download_placement_files(self):
        self.ref_markers_file = self.downloader.get("list_of_reference_markers.{0}_{1}.txt".format(
            os.path.basename(self.run_folder).split("_")[-2],
            self.datasets_version),
            "placement_files")
        self.tree_nwk_file = self.downloader.get("tree.{0}_{1}.nwk".format(
            os.path.basename(self.run_folder).split("_")[-2],
            self.datasets_version),
            "placement_files")
        self.tree_metadata_file = self.downloader.get("tree_metadata.{0}_{1}.txt".format(
            os.path.basename(self.run_folder).split("_")[-2],
            self.datasets_version),
            "placement_files")
        self.supermatrix_file = self.downloader.get("supermatrix.aln.{0}_{1}.faa".format(
            os.path.basename(self.run_folder).split("_")[-2],
            self.datasets_version),
            "placement_files")
        self.taxid_busco_file = self.downloader.get("mapping_taxids-busco_dataset_name.{0}_{1}.txt".format(
            os.path.basename(self.run_folder).split("_")[-2],
            self.datasets_version),
            "placement_files")
        self.taxid_lineage_file = self.downloader.get("mapping_taxid-lineage.{0}_{1}.txt".format(
            os.path.basename(self.run_folder).split("_")[-2],
            self.datasets_version),
            "placement_files")
        return

    def _get_placement_file_versions(self):
        placement_file_versions = [
            os.path.basename(filepath)
            for filepath in [self.ref_markers_file, self.tree_nwk_file, self.tree_metadata_file,
                             self.supermatrix_file, self.taxid_busco_file, self.taxid_lineage_file]]
        return placement_file_versions


    @log("Extract markers...", logger)
    def define_dataset(self):
        # If mode is genome, substitute input with prodigal/augustus output
        self._download_placement_files()
        placement_file_versions = self._get_placement_file_versions()
        self._extract_marker_sequences()
        self._run_sepp()

        dataset = self._pick_dataset()

        return dataset, placement_file_versions

    def _init_tools(self):
        try:
            assert isinstance(self._sepp, Tool)
        except AttributeError:
            self._sepp = Tool("sepp", self._config)
        except AssertionError:
            raise SystemExit("SEPP should be a tool")

    def _pick_dataset(self):

        run_folder = self.run_folder

        # load busco dataset name by id in a dict {taxid:name}
        datasets_mapping = {}

        with open(self.taxid_busco_file) as f:
            for line in f:
                datasets_mapping.update(
                    {line.strip().split("\t")[0]: line.strip().split("\t")[1].split(",")[0]}
                )

        # load the lineage for each taxid in a dict {taxid:reversed_lineage}
        # lineage is 1:2:3:4:5:6 => {6:[6,5,4,3,2,1]}
        lineages = set()
        parents = {}
        taxid_dataset = {}
        for t in datasets_mapping:
            taxid_dataset.update({t: t})
        with open(self.taxid_lineage_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                lineages.add(line.strip().split("\t")[4])

                i = 0
                # for each line, e.g. 6\t1:2:3:4:5:6, create/update the lineage for each level
                # 6:[1,2,3,4,5,6], 5:[1,2,3,4,5], 4:[1,2,3,4], etc.
                for t in line.strip().split("\t")[4].split(","):
                    i += 1
                    parents.update({t: line.strip().split("\t")[4].split(",")[0:i]})

        for t in parents:
            for p in parents[t][::-1]: # reverse the order to get the deepest parent, not the root one
                if p in datasets_mapping:
                    taxid_dataset.update({t: p})
                    break
        # load json
        # load "tree" in a string
        # load placements
        # obtain a dict of taxid num of markers
        # figure out which taxid to use by using the highest number of markers and some extra rules

        try:
            with open(os.path.join(self.placement_folder, "output_placement.json")) as json_file:
                data = json.load(json_file)
            tree = data["tree"]
            placements = data["placements"]
        except FileNotFoundError:
            raise SystemExit("Placements failed. Try to rerun increasing the memory or select a lineage manually.")

        node_weight = {}
        n_p = 0
        for placement in placements:
            n_p += 1
            for individual_placement in placement["p"]:
                # find the taxid in tree
                node = individual_placement[0]
                import re

                match = re.findall(  # deal with weird character in the json file, see the output yourself.
                    # if this pattern is inconsistant with pplacer version, it may break buscoplacer.
                    "[^0-9][0-9]*:[0-9]*[^0-9]{0,1}[0-9]*[^0-9]{0,2}[0-9]*\[%s\]"
                    % node,
                    tree,
                )
                # extract taxid:
                try:
                    if re.match("^[A-Za-z]", match[0]):
                        taxid = match[0][7:].split(":")[0]
                    else:
                        taxid = match[0][1:].split(":")[0]
                except IndexError as e:
                    raise e
                if taxid_dataset[taxid] in node_weight:
                    node_weight[taxid_dataset[taxid]] += 1
                else:
                    node_weight[taxid_dataset[taxid]] = 1
                break  # Break here to keep only the best match. In my experience, keeping all does not change much.
        type(self)._logger.debug("Placements counts by node are: %s" % node_weight)

        # from here, define which placement can be trusted
        max_markers = 0
        choice = []

        for key in node_weight:
            type(self)._logger.debug('%s markers assigned to the taxid %s' % (node_weight[key],key))

        # taxid for which no threshold or minimal amount of placement should be considered.
        # If it is the best, go for it.
        no_rules = ["204428"]

        ratio = 2.5
        if run_folder.split("/")[-1].split("_")[-2] == "archaea":
            ratio = 1.2
        min_markers = 12

        node_with_max_markers = None
        for n in node_weight:
            if node_weight[n] > max_markers:
                max_markers = node_weight[n]
                node_with_max_markers = n
        if node_with_max_markers in no_rules:
            choice = [node_with_max_markers]
        else:
            for n in node_weight:
                # if the ration between the best and the current one is not enough, keep both
                if node_weight[n] * ratio >= max_markers:
                    choice.append(n)
        if len(choice) > 1:
            # more than one taxid should be considered, pick the common ancestor
            choice = self._get_common_ancestor(choice, parents)
            #print('last common')
            #print(choice)
        elif len(choice) == 0:
            if run_folder.split("/")[-1].split("_")[-2] == "bacteria":
                choice.append("2")
            elif run_folder.split("/")[-1].split("_")[-2] == "archaea":
                choice.append("2157")
            elif run_folder.split("/")[-1].split("_")[-2] == "eukaryota":
                choice.append("2759")
        if max_markers < min_markers and not (choice[0] in no_rules):
            if run_folder.split("/")[-1].split("_")[-2] == "bacteria":
                key_taxid = "2"
            elif run_folder.split("/")[-1].split("_")[-2] == "archaea":
                key_taxid = "2157"
            elif run_folder.split("/")[-1].split("_")[-2] == "eukaryota":
                key_taxid = "2759"
            else:
                key_taxid = None  # unexpected. Should throw an exception or use assert.
            type(self)._logger.info(
                "Not enough markers were placed on the tree (%s). Root lineage %s is kept" % (max_markers, datasets_mapping[taxid_dataset[key_taxid]])
            )
            return [
                datasets_mapping[taxid_dataset[key_taxid]],
                max_markers,
                sum(node_weight.values()),
            ]

        type(self)._logger.info('Lineage %s is selected, supported by %s markers out of %s' % (datasets_mapping[taxid_dataset[choice[0]]],max_markers,sum(node_weight.values())))

        return [
            datasets_mapping[taxid_dataset[choice[0]]],
            max_markers,
            sum(node_weight.values()),
        ]

    @staticmethod
    def _get_common_ancestor(choice, parents):
        # starts with the parents of the first choice
        all_ancestors = set(parents[choice[0]])
        # order will be lost with sets, so keep in a list the lineage of one entry to later pick the deepest ancestor
        ordered_lineage = []
        for c in choice:
            #print('c is %s' % c)
            if len(parents[c]) > len(
                ordered_lineage
            ):
            #    print('len parent c is %s' % len(parents[c]))
            #    print('len ordered lineage is %s' % ordered_lineage)
                # probably useless. Init with parents[choice[0] should work
                ordered_lineage = parents[c]
            #    print('ordered_lineage us %s' % ordered_lineage)
            # keep in set only entries that are in the currently explored lineage
            all_ancestors = all_ancestors.intersection(parents[c])

        # go through the ordered list of the deepest linage until you found a common ancestor.
        for parent in ordered_lineage[::-1]:
            if parent in all_ancestors:
                return [parent]

    @log("Place the markers on the reference tree...", logger)
    def _run_sepp(self):
        self.sepp_runner = SEPPRunner(self._sepp, self.run_folder, self.placement_folder, self.tree_nwk_file,
                                      self.tree_metadata_file, self.supermatrix_file, self.downloader,
                                      self.datasets_version, self.cpus)
        self.sepp_runner.run()

    def _extract_marker_sequences(self):
        """
        This function extracts all single copy BUSCO genes from a protein run folder
        :param run_folder: a BUSCO protein run folder
        :type: str
        :param protein_sequences: protein fasta used as input for the BUSCO run corresponding to the folder
        :type: str
        """



        with open(self.ref_markers_file, "r") as f:
            marker_list = [line.strip() for line in f]


        marker_genes_names = []
        for busco, gene_matches in self.single_copy_buscos.items():
            if busco in marker_list:
                marker_genes_names.append(list(gene_matches.keys())[0])  # The list should only have one entry because they are single copy buscos

        marker_genes_records = []
        if isinstance(self.protein_seqs, (str,)):
            list_protein_seqs = [self.protein_seqs]
        else:
            list_protein_seqs = self.protein_seqs

        for protein_seqs in list_protein_seqs:
            with open(protein_seqs, "r") as prot_seqs:
                for record in SeqIO.parse(prot_seqs, "fasta"):
                    if record.id in marker_genes_names:
                        record.seq = record.seq.rstrip("*")
                        record.description = ""
                        marker_genes_records.append(record)

        marker_genes_file = os.path.join(self.placement_folder, "marker_genes.fasta")
        with open(marker_genes_file, "w") as output:
            SeqIO.write(marker_genes_records, output, "fasta")

