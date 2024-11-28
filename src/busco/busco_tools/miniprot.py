# coding: utf-8
"""
miniprot.py

Module for running Miniprot.

Author(s): Matthew Berkeley

Copyright (c) 2015-2024, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

from busco.busco_tools.base import GenePredictor, BaseRunner
from busco.BuscoLogger import BuscoLogger
import subprocess
import os
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import numpy as np
import re
from multiprocessing import Pool

logger = BuscoLogger.get_logger(__name__)


class MiniprotRunner(BaseRunner):
    """
    Class to run Miniprot
    """

    name = "miniprot"
    cmd = "miniprot"

    def __init__(self):
        super().__init__()
        self._output_folder = os.path.join(self.run_folder, "miniprot_output")
        self.translated_proteins_folder = os.path.join(
            self._output_folder, "translated_proteins"
        )
        self.create_dirs(
            [
                self._output_folder,
                self.translated_proteins_folder,
            ]
        )
        self.index_file = os.path.join(self._output_folder, "ref.mpi")
        self.refseq_db = None
        self.incomplete_buscos = None
        self._output_basename = None

    def configure_runner(self, *args):
        """
        Configure Miniprot runner
        """
        super().configure_runner(*args)
        self.run_number += 1

    def check_tool_dependencies(self):
        pass

    def configure_job(self, *args):
        """
        Overridden by child classes
        """
        return

    def generate_job_args(self):
        yield

    def get_version(self):
        help_output = subprocess.check_output(
            [self.cmd, "--version"], stderr=subprocess.STDOUT, shell=False
        )
        version = help_output.decode("utf-8").strip()
        return version

    @property
    def output_folder(self):
        return self._output_folder

    def reset(self):
        super().reset()

    def run(self):
        super().run()
        self.total = 1
        self.run_jobs()


class MiniprotIndexRunner(MiniprotRunner):
    name = "miniprot_index"

    def generate_job_args(self):
        yield "index"

    def configure_job(self, *args):
        """
        Configure Miniprot job
        """

        miniprot_job = self.create_job()
        miniprot_job.add_parameter("-t")
        miniprot_job.add_parameter(str(self.cpus))
        miniprot_job.add_parameter("-d")
        miniprot_job.add_parameter(self.index_file)
        miniprot_job.add_parameter(self.input_file)

        return miniprot_job


class MiniprotAlignRunner(MiniprotRunner, GenePredictor):
    name = "miniprot_align"

    def __init__(self):
        super().__init__()
        self.output_gff = None

        self.sequences_aa = {}
        self.busco_matches = defaultdict(set)
        self.gene_matches = defaultdict(list)
        self.combined_pred_protein_seqs = os.path.join(
            self._output_folder, "combined_pred_proteins.fas"
        )
        self.output_sequences = []
        self.gene_nominal = 0
        self.gene_lookup = {}
        self.cigar_lookup = {}
        self.nominal_lookup = defaultdict(list)

        self.gene_details = defaultdict(dict)
        self.sequences_aa = {}
        self.sequences_nt = {}

    def generate_job_args(self):
        yield "align"

    def configure_job(self, *args):

        miniprot_job = self.create_job()
        miniprot_job.add_parameter("--trans")
        miniprot_job.add_parameter("-u")
        miniprot_job.add_parameter("-I")
        miniprot_job.add_parameter("--outs")
        miniprot_job.add_parameter("0.95")
        miniprot_job.add_parameter("-t")
        miniprot_job.add_parameter(str(self.cpus))
        miniprot_job.add_parameter("--gff")

        miniprot_job.add_parameter(self.index_file)
        miniprot_job.add_parameter(self.refseq_db)

        return miniprot_job

    def configure_runner(self, incomplete_buscos=None):
        super().configure_runner([])
        self.logfile_path_out = os.path.join(
            self.config.get("busco_run", "main_out"),
            "logs",
            "{}_{}_out.log".format(self.name, os.path.basename(self.lineage_dataset)),
        )
        self.logfile_path_err = (
            self.logfile_path_out.rpartition("_out.log")[0] + "_err.log"
        )

        self.incomplete_buscos = incomplete_buscos
        self._output_basename = os.path.join(
            self._output_folder, os.path.basename(self.input_file)
        )
        gzip_refseq = os.path.join(self.lineage_dataset, "refseq_db.faa.gz")
        self.refseq_db = self.decompress_refseq_file(gzip_refseq)
        self.output_gff = Path(self._output_folder).joinpath(
            "{}_{}{}".format(
                Path(self.input_file).stem,
                os.path.basename(self.lineage_dataset),
                ".gff",
            )
        )

    def create_symlink(self):
        if not self.output_gff.exists():
            Path(self.output_gff).symlink_to(
                os.path.relpath(self.logfile_path_out, self.output_gff.parent)
            )
        return

    def save_record(
        self,
        gene_id,
        target_id,
        contig_id,
        contig_start,
        contig_end,
        strand,
        score,
        exon_coords,
        ata_seq,
        protein_start,
        protein_end,
        protein_length,
        cigar_seq,
        frameshifts,
        frameshift_events,
        frameshift_lengths,
        stop_codon_count,
        identity,
        gff_start,
        gff_end,
    ):
        self.all_gff_records.append(
            np.array(
                (
                    gene_id,
                    target_id.split("_")[0],
                    target_id,
                    contig_id,
                    contig_start,
                    contig_end,
                    strand,
                    score,
                    exon_coords,
                    ata_seq,
                    protein_start,
                    protein_end,
                    protein_length,
                    cigar_seq,
                    frameshifts,
                    frameshift_events,
                    frameshift_lengths,
                    stop_codon_count,
                    identity,
                    gff_start,
                    gff_end,
                ),
                dtype=[
                    ("gene_id", "U500"),
                    ("busco_id", "U100"),
                    ("target_id", "U500"),
                    ("contig_id", "U500"),
                    ("contig_start", "i4"),
                    ("contig_end", "i4"),
                    ("strand", "U1"),
                    ("score", "i4"),
                    ("exon_coords", "O"),
                    ("aa_seq", "U10000"),
                    ("protein_start", "i4"),
                    ("protein_end", "i4"),
                    ("protein_length", "i4"),
                    ("cigar_seq", "U10000"),
                    ("frameshifts", "O"),
                    ("frameshift_events", "i4"),
                    ("frameshift_lengths", "i4"),
                    ("stop_codon_count", "i4"),
                    ("identity", "f4"),
                    ("gff_start", "i4"),
                    ("gff_end", "i4"),
                ],
            )
        )
        return

    def parse_output(self):
        self.create_symlink()
        exon_coords = []
        self.all_gff_records = []

        paf_block_started = False
        gene_id = ""
        identity = None
        line_no = -1
        gff_block = False
        gff_start = None

        with open(self.output_gff, "r") as gff:

            for line in gff:
                line_no += 1
                if line.startswith("##PAF"):
                    gff_block = False
                    gff_end = line_no
                    if identity and len(exon_coords) > 0 and gff_start:
                        self.save_record(
                            gene_id,
                            target_id,
                            contig_id,
                            contig_start,
                            contig_end,
                            strand,
                            score,
                            exon_coords,
                            ata_seq,
                            protein_start,
                            protein_end,
                            protein_length,
                            cigar_seq,
                            frameshifts,
                            frameshift_events,
                            frameshift_lengths,
                            stop_codon_count,
                            identity,
                            gff_start,
                            gff_end,
                        )
                        identity = None
                        exon_coords = []
                        paf_block_started = False

                    fields = line.strip().split("\t")[1:]
                    if fields[5] == "*":
                        ## Unmapped protein
                        continue
                    paf_block_started = True
                    exon_coords = []
                    target_id = fields[0]
                    protein_length = int(fields[1])
                    protein_start = int(fields[2])
                    protein_end = int(fields[3])
                    strand = fields[4]
                    contig_id = fields[5]
                    contig_start = int(fields[7])
                    contig_end = int(fields[8])
                    gene_id = "{}|{}:{}-{}|{}".format(
                        target_id, contig_id, contig_start, contig_end, strand
                    )
                    score = int(fields[13].strip().split(":")[2])

                    cigar_seq = str(fields[17].strip().split(":")[2])
                    (
                        frameshifts,
                        frameshift_events,
                        frameshift_lengths,
                    ) = self.decode_cigar(cigar_seq)
                    sta_line = gff.readline()
                    line_no += 1
                    sta_seq = sta_line.strip().split("\t")[1]
                    ata_seq = sta_seq.upper()
                    stop_codon_count = ata_seq.strip("*").count("*")
                    ata_seq = ata_seq.replace("*", "")

                elif paf_block_started:
                    if not gff_block:
                        gff_start = line_no
                    gff_block = True
                    fields = line.strip().split("\t")
                    if fields[2] == "CDS":
                        start, stop, score, strand = (
                            int(fields[3]),
                            int(fields[4]),
                            float(fields[5]),
                            fields[6],
                        )
                        exon_coords.append((start, stop, score, strand))
                    if fields[2] == "mRNA":
                        info_dict = dict(
                            v.split("=") for v in fields[8].split()[0].split(";")
                        )
                        identity = float(info_dict["Identity"])
            else:
                gff_end = line_no
                self.save_record(
                    gene_id,
                    target_id,
                    contig_id,
                    contig_start,
                    contig_end,
                    strand,
                    score,
                    exon_coords,
                    ata_seq,
                    protein_start,
                    protein_end,
                    protein_length,
                    cigar_seq,
                    frameshifts,
                    frameshift_events,
                    frameshift_lengths,
                    stop_codon_count,
                    identity,
                    gff_start,
                    gff_end,
                )
        self.gff_arr = np.array(self.all_gff_records)
        return

    def filter(self):
        contigs = np.unique(self.gff_arr["contig_id"])
        self.filtered_matches = np.array([], dtype=self.gff_arr.dtype)
        if (
            len(contigs) > 1700
        ):  # 1700 seems to be around the number where parallel post-processing is faster
            with Pool(self.cpus) as pool:
                results = pool.map(
                    self.filter_contig,
                    contigs,
                    chunksize=max(1, len(contigs) // self.cpus),
                )
            stacked_results = np.vstack([np.vstack(res) for res in results])
            self.filtered_matches = np.squeeze(stacked_results, axis=1)
        else:
            for contig in contigs:
                contig_filtered = self.filter_contig(contig)
                self.filtered_matches = np.concatenate(
                    (self.filtered_matches, contig_filtered)
                )

        return

    def filter_contig(self, contig):

        contig_matches = self.gff_arr[self.gff_arr["contig_id"] == contig]
        for i in range(len(contig_matches)):
            if contig_matches["contig_start"][i] < 0:
                continue  # Already removed

            overlap_mask = (
                (
                    (
                        (
                            contig_matches["contig_start"]
                            >= contig_matches["contig_start"][i]
                        )
                        & (
                            contig_matches["contig_start"]
                            < contig_matches["contig_end"][i]
                        )
                    )
                    | (
                        (
                            contig_matches["contig_end"]
                            > contig_matches["contig_start"][i]
                        )
                        & (
                            contig_matches["contig_end"]
                            <= contig_matches["contig_end"][i]
                        )
                    )
                )
                & (contig_matches["contig_start"] > 0)
                & (contig_matches["busco_id"] == contig_matches["busco_id"][i])
            )
            overlap_inds = np.arange(len(contig_matches))[overlap_mask]
            for j in overlap_inds:
                if i == j:
                    continue

                length1 = (
                    contig_matches["contig_end"][i] - contig_matches["contig_start"][i]
                )
                length2 = (
                    contig_matches["contig_end"][j] - contig_matches["contig_start"][j]
                )
                max_length = max(length1, length2)

                overlap_start = max(
                    contig_matches["contig_start"][i],
                    contig_matches["contig_start"][j],
                )
                overlap_end = min(
                    contig_matches["contig_end"][i], contig_matches["contig_end"][j]
                )
                overlap_length = overlap_end - overlap_start

                if overlap_length > 0.80 * max_length:
                    if contig_matches["score"][j] > contig_matches["score"][i]:
                        contig_matches["contig_start"][i] = -1  # Remove match
                        break
                    else:
                        contig_matches["contig_start"][j] = -1  # Remove match
        filtered_contig_matches = contig_matches[contig_matches["contig_start"] != -1]
        return filtered_contig_matches

    def record_gene_details(self):
        for match in self.filtered_matches:
            gene_id = match["gene_id"].item()
            target_id = match["target_id"].item()
            contig_id = match["contig_id"].item()
            contig_start = match["contig_start"].item()
            contig_end = match["contig_end"].item()
            strand = match["strand"].item()
            score = match["score"].item()
            exon_coords = match["exon_coords"]
            ata_seq = match["aa_seq"].item()
            protein_start = match["protein_start"].item()
            protein_end = match["protein_end"].item()
            protein_length = match["protein_length"].item()
            cigar_seq = match["cigar_seq"].item()
            frameshifts = match["frameshifts"]
            frameshift_events = match["frameshift_events"].item()
            frameshift_lengths = match["frameshift_lengths"].item()
            stop_codon_count = match["stop_codon_count"].item()
            identity = match["identity"].item()
            gff_start = match["gff_start"].item()
            gff_end = match["gff_end"].item()

            self.gene_details[gene_id].update(
                {
                    "gene_id": gene_id,
                    "target_id": target_id,
                    "contig_id": contig_id,
                    "contig_start": contig_start,
                    "contig_end": contig_end,
                    "strand": strand,
                    "score": score,
                    "exon_coords": exon_coords,
                    "aa_seq": SeqRecord(Seq(ata_seq), id=gene_id),
                    "protein_start": protein_start,
                    "protein_end": protein_end,
                    "protein_length": protein_length,
                    "cigar_seq": cigar_seq,
                    "frameshifts": frameshifts,
                    "frameshift_events": frameshift_events,
                    "frameshift_lengths": frameshift_lengths,
                    "stop_codon_count": stop_codon_count,
                    "identity": identity,
                    "run_number": self.run_number,
                    "gff_start": gff_start,
                    "gff_end": gff_end,
                }
            )
            self.sequences_aa[gene_id] = SeqRecord(Seq(ata_seq), id=gene_id)
            self.busco_matches[target_id.split("_")[0]].add(gene_id)
        return

    def check_overlap(self, contig_id, contig_start, contig_end, score, gene_id):
        keeper = gene_id
        matches_to_remove = []
        matches = np.array(self.genomic_regions[contig_id])
        overlap_matches = matches[
            (
                (matches["contig_start"] >= contig_start)
                & (matches["contig_start"] < contig_end)
            )
            | (
                (matches["contig_end"] > contig_start)
                & (matches["contig_end"] <= contig_end)
            )
        ]
        if len(overlap_matches) > 0:
            for match in overlap_matches:
                overlap_start = max(contig_start, match["contig_start"])
                overlap_end = min(contig_end, match["contig_end"])
                overlap_length = overlap_end - overlap_start
                if (overlap_length > 0.8 * (contig_end - contig_start)) or (
                    overlap_length > 0.8 * (match["contig_end"] - match["contig_start"])
                ):
                    if score > match["score"]:
                        keeper = gene_id
                        matches_to_remove.append(match)
                    else:
                        keeper = match["gene_id"]
        if len(matches_to_remove) > 0:
            self.genomic_regions[contig_id] = [
                m for m in matches if m not in matches_to_remove
            ]
        return keeper

    @staticmethod
    def decode_cigar(cigar):
        frameshifts = []
        frameshift_events = 0
        frameshift_lengths = 0
        pattern = r"[0-9]+[MIDFGNUV]"
        parts = list(re.finditer(pattern, cigar))
        for p, part in enumerate(parts):
            n, type = int(part.group(0)[:-1]), part.group(0)[-1]
            if type in ["M", "D", "U", "V", "N"]:
                continue
            elif type in ["F", "G"]:
                # left search
                q = p - 1
                left_match_cnt = 0
                while q >= 0:
                    part2 = parts[q]
                    n2, type2 = int(part2.group(0)[:-1]), part2.group(0)[-1]
                    if type2 == "M":
                        left_match_cnt += n2
                    elif type2 in ["N", "U", "V"]:
                        break
                    q -= 1
                # right search
                q = p + 1
                right_match_cnt = 0
                while q < len(parts):
                    part2 = parts[q]
                    n2, type2 = int(part2.group(0)[:-1]), part2.group(0)[-1]
                    if type2 == "M":
                        right_match_cnt += n2
                    elif type2 in ["N", "U", "V"]:
                        break
                    q += 1
                if left_match_cnt >= 20 and right_match_cnt >= 20:
                    frameshifts.append(str(n) + type)
                    frameshift_events += 1
                    frameshift_lengths += int(n)

        return frameshifts, frameshift_events, frameshift_lengths

    def write_protein_sequences_per_busco(self):
        for busco_id in self.busco_matches:
            seqs_to_write = []
            output_filename = os.path.join(
                self.translated_proteins_folder, "{}.faa".format(busco_id)
            )
            self.output_sequences.append(output_filename)
            with open(output_filename, "w") as f:
                for g in self.busco_matches[busco_id]:
                    if g in self.sequences_aa:
                        seqs_to_write.append(self.sequences_aa[g])
                SeqIO.write(seqs_to_write, f, "fasta")
