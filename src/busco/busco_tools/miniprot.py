from busco.busco_tools.base import BaseRunner
from busco.BuscoLogger import BuscoLogger
import subprocess
import os
from pathlib import Path
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import shutil
import numpy as np
import re

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


class MiniprotAlignRunner(MiniprotRunner):

    name = "miniprot_align"

    def __init__(self):
        super().__init__()
        self.output_gff = None

        self.gene_details = defaultdict(dict)
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
            "{}_{}{}".format(Path(self.input_file).stem, os.path.basename(self.lineage_dataset), ".gff"))

    def create_symlink(self):
        if not self.output_gff.exists():
            Path(self.output_gff).symlink_to(self.logfile_path_out)
        return

    def parse_output(self):
        self.create_symlink()
        self.ata_seq = ""
        self.target_id = ""
        self.contig_id = ""
        self.contig_start = 0
        self.contig_end = 0
        self.strand = ""
        self.score = 0
        self.exon_coords = defaultdict(list)
        self.cigar_seq = ""
        paf_block_started = False
        gene_id = ""

        with open(self.output_gff, "r") as gff:
            for line in gff:
                if line.startswith("##PAF"):
                    paf_block_started = True
                    fields = line.strip().split("\t")[1:]
                    if fields[5] == "*":
                        ## Unmapped protein
                        continue
                    self.target_id = fields[0]
                    busco_id = self.target_id.split("_")[0]
                    self.protein_length = int(fields[1])
                    self.protein_start = int(fields[2])
                    self.protein_end = int(fields[3])
                    self.strand = fields[4]
                    self.contig_id = fields[5]
                    self.contig_start = int(fields[7])
                    self.contig_end = int(fields[8])
                    gene_id = "{}|{}:{}-{}".format(self.target_id, self.contig_id, self.contig_start, self.contig_end)
                    self.score = int(fields[13].strip().split(":")[2])
                    self.cigar_seq = str(fields[17].strip().split(":")[2])
                    part_lengths, exon_lengths, match_lengths, group_types, ngroups, nexons, frameshifts, \
                        frameshift_events, frameshift_lengths = self.decode_cigar(self.cigar_seq)
                    sta_line = gff.readline()
                    sta_seq = sta_line.strip().split("\t")[1]
                    self.ata_seq = re.sub("\*", "", sta_seq.upper())

                    self.busco_matches[busco_id].add(gene_id)

                    self.gene_details[gene_id] = {"gene_start": self.contig_start, "gene_end": self.contig_end,
                                                       "strand": self.strand, "score": self.score,
                                                       "cigar": self.cigar_seq, "nexons": nexons,
                                                       "frameshift_events": frameshift_events,
                                                       "protein_start": self.protein_start,
                                                       "protein_end": self.protein_end,
                                                       "protein_length": self.protein_length}

                    self.sequences_aa[gene_id] = SeqRecord(Seq(self.ata_seq), id=gene_id, description=gene_id)

                elif paf_block_started:
                    fields = line.strip().split("\t")
                    if fields[2] == "CDS":
                        start, stop, score, strand = fields[3], fields[4], fields[5], fields[6]
                        self.exon_coords[gene_id].append((start, stop, score, strand))
                    if fields[2] == "mRNA":
                        info_dict = dict(v.split("=") for v in fields[8].split()[0].split(";"))
                        identity = float(info_dict["Identity"])
                        self.gene_details[gene_id].update({"identity": identity})
        for item in self.exon_coords:
            self.exon_coords[item] = np.array(self.exon_coords[item], dtype=[("start", "i4"), ("stop", "i4"), ("score", "f4"), ("strand", "U1")])
        return

    @staticmethod
    def decode_cigar(cigar):
        frameshifts = []
        frameshift_events = 0
        frameshift_lengths = 0
        pattern = r"[0-9]+[MIDFGNUV]"
        parts = list(re.finditer(pattern, cigar))
        part_lengths = []
        exon_lengths = []
        exon_length = 0
        match_lengths = {"M": 0, "I": 0, "D": 0, "F": 0, "G": 0, "N": 0, "U": 0, "V": 0}
        group_types = {"M": 0, "I": 0, "D": 0, "F": 0, "G": 0, "N": 0, "U": 0, "V": 0}
        ngroups = 0
        nexons = 0
        for p, part in enumerate(parts):
            ngroups += 1
            n, type = int(part.group(0)[:-1]), part.group(0)[-1]
            match_lengths[type] += n
            group_types[type] += 1
            if type in ["M", "D"]:
                exon_length += n
            elif type in ["U", "V"]:
                part_lengths.append(exon_length)
                exon_lengths.append(exon_length)
                nexons += 1
                part_lengths.append(1)
                exon_length = 0
            elif type == "N":
                part_lengths.append(exon_length)
                exon_lengths.append(exon_length)
                nexons += 1
                exon_length = 0
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

            # elif type == "G":
            #     exon_length += 1
        part_lengths.append(exon_length)
        exon_lengths.append(exon_length)
        nexons += 1
        return part_lengths, exon_lengths, match_lengths, group_types, ngroups, nexons, frameshifts, frameshift_events, frameshift_lengths

    def write_protein_sequences_per_busco(self):
        for busco_id in self.busco_matches:
            seqs_to_write = []
            output_filename = os.path.join(self.translated_proteins_folder, "{}.faa".format(busco_id))
            self.output_sequences.append(output_filename)
            with open(output_filename, "w") as f:
                for g in self.busco_matches[busco_id]:
                    if g in self.sequences_aa:
                        seqs_to_write.append(self.sequences_aa[g])
                SeqIO.write(seqs_to_write, f, "fasta")
