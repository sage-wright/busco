from Bio import SeqIO
from busco.BuscoTools import TBLASTNRunner, MKBLASTRunner
from busco.Toolset import Tool
from busco.BuscoLogger import BuscoLogger
from busco.BuscoLogger import LogDecorator as log
import subprocess
import os
from abc import ABCMeta, abstractmethod

logger = BuscoLogger.get_logger(__name__)


class NucleotideAnalysis(metaclass=ABCMeta):

    LETTERS = ["A", "C", "T", "G", "N"]

    # explanation of ambiguous codes found here: https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
    AMBIGUOUS_CODES = ["Y", "R", "W", "S", "K", "M", "D", "V", "H", "B"]

    MAX_FLANK = 20000

    def __init__(self, config):
        # Variables inherited from BuscoAnalysis
        self._config = None
        self._cpus = None
        self._input_file = None

        super().__init__(config)  # Initialize BuscoAnalysis
        self._long = self._config.getboolean("busco_run", "long")
        self._flank = self._define_flank()
        self._ev_cutoff = self._config.getfloat("busco_run", "evalue")
        self._region_limit = self._config.getint("busco_run", "limit")
        self.blast_cpus = self._cpus

        if not self.check_nucleotide_file(self._input_file):
            raise SystemExit("Please provide a nucleotide file as input")

    def check_nucleotide_file(self, filename):

        i = 0
        for record in SeqIO.parse(filename, "fasta"):
            for letter in record.seq.upper():
                if i > 5000:
                    break
                i += 1
                if letter not in type(self).LETTERS and letter not in type(self).AMBIGUOUS_CODES:
                    return False
            else:
                continue  # only continue to next record of 5000 has not been hit
            break  # If for loop exits with "break", the else clause is skipped and the outer loop also breaks.

        return True

    def _define_flank(self):
        """
        TODO: Add docstring
        :return:
        """
        try:
            size = os.path.getsize(self._input_file) / 1000  # size in mb
            flank = int(size / 50)  # proportional flank size
            # Ensure value is between 5000 and MAX_FLANK
            flank = min(max(flank, 5000), type(self).MAX_FLANK)
        except IOError:  # Input data is only validated during run_analysis. This will catch any IO issues before that.
            raise SystemExit("Impossible to read the fasta file {}".format(self._input_file))

        return flank

    @abstractmethod
    def init_tools(self):  # todo: This should be an abstract method
        """
        Initialize all required tools for Genome Eukaryote Analysis:
        MKBlast, TBlastn, Augustus and Augustus scripts: GFF2GBSmallDNA, new_species, etraining
        :return:
        """
        super().init_tools()


    def check_tool_dependencies(self):
        super().check_tool_dependencies()

    def _get_blast_version(self):
        blast_version_call = subprocess.check_output([self._tblastn_tool.cmd, "-version"], shell=False)
        blast_version = ".".join(blast_version_call.decode("utf-8").split("\n")[0].split()[1].rsplit(".")[:-1])
        return blast_version

    def _run_mkblast(self):
        self.mkblast_runner = MKBLASTRunner(self._mkblast_tool, self._input_file, self.main_out, self._cpus)
        self.mkblast_runner.run()

    @log("Running a BLAST search for BUSCOs against created database", logger)
    def _run_tblastn(self, missing_and_frag_only=False, ancestral_variants=False):

        incomplete_buscos = (self.hmmer_runner.missing_buscos + list(self.hmmer_runner.fragmented_buscos.keys())
                             if missing_and_frag_only else None)  # This parameter is only used on the re-run

        self.tblastn_runner = TBLASTNRunner(self._tblastn_tool, self._input_file, self.run_folder, self._lineage_dataset,
                                            self.mkblast_runner.output_db, self._ev_cutoff, self.blast_cpus,
                                            self._region_limit, self._flank, missing_and_frag_only, ancestral_variants,
                                            incomplete_buscos)

        self.tblastn_runner.run()
        coords = self.tblastn_runner._get_coordinates()
        coords = self.tblastn_runner._filter_best_matches(coords)  # Todo: remove underscores from non-hidden methods
        self.tblastn_runner._write_coordinates_to_file(coords)  # writes to "coordinates.tsv"
        self.tblastn_runner._write_contigs(coords)
        return coords


class ProteinAnalysis:

    LETTERS = ["F", "L", "I", "M", "V", "S", "P", "T", "A", "Y", "X", "H", "Q", "N", "K", "D", "E", "C", "W", "R", "G"]
    NUCL_LETTERS = ["A", "C", "T", "G", "N"]

    def __init__(self, config):
        super().__init__(config)
        if not self.check_protein_file(self._input_file):
            raise SystemExit('Please provide a protein file as input')

    def check_protein_file(self, filename):

        for i, record in enumerate(SeqIO.parse(filename, "fasta")):
            if i > 10:
                break
            for letter in record.seq:
                if letter.upper() not in type(self).NUCL_LETTERS and letter.upper() in type(self).LETTERS:
                    return True
                elif letter.upper() not in type(self).LETTERS:
                    return False
                else:
                    continue
        return False  # if file only contains "A", "T", "C", "G", "N", it is unlikely to be a protein file
