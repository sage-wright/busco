from Bio import SeqIO
from busco.BuscoLogger import BuscoLogger
from abc import ABCMeta
import os
from busco.busco_tools.blast import TBLASTNRunner, MKBLASTRunner
from busco.Exceptions import BuscoError

logger = BuscoLogger.get_logger(__name__)


class NucleotideAnalysis(metaclass=ABCMeta):

    LETTERS = ["A", "C", "T", "G", "N"]

    # explanation of ambiguous codes found here: https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
    AMBIGUOUS_CODES = ["Y", "R", "W", "S", "K", "M", "D", "V", "H", "B"]

    def __init__(self):

        super().__init__()  # Initialize BuscoAnalysis
        if not self.check_nucleotide_file(self.input_file):
            raise BuscoError("The input file does not contain nucleotide sequences.")

    def check_nucleotide_file(self, filename):
        i = 0
        num_records = 0
        for record in SeqIO.parse(filename, "fasta"):
            num_records += 1
            for letter in record.seq.upper():
                if i > 5000:
                    break
                i += 1
                if (
                    letter not in type(self).LETTERS
                    and letter not in type(self).AMBIGUOUS_CODES
                ):
                    return False
            else:
                continue  # only continue to next record if 5000 has not been hit
            break  # If for loop exits with "break", the else clause is skipped and the outer loop also breaks.

        if num_records == 0:
            return False

        return True

    def init_tools(self):
        super().init_tools()

    def reset(self):
        super().reset()


class ProteinAnalysis:

    LETTERS = [
        "F",
        "L",
        "I",
        "M",
        "V",
        "S",
        "P",
        "T",
        "A",
        "Y",
        "X",
        "H",
        "Q",
        "N",
        "K",
        "D",
        "E",
        "C",
        "W",
        "R",
        "G",
        "*",
    ]
    NUCL_LETTERS = ["A", "C", "T", "G", "N"]

    def __init__(self):
        super().__init__()
        if not self.check_protein_file(self.input_file):
            raise BuscoError("Please provide a protein file as input")

    def check_protein_file(self, filename):

        for i, record in enumerate(SeqIO.parse(filename, "fasta")):
            if i > 10:
                break
            for letter in record.seq:
                if (
                    letter.upper() not in type(self).NUCL_LETTERS
                    and letter.upper() in type(self).LETTERS
                ):
                    return True
                elif letter.upper() not in type(self).LETTERS:
                    return False
                else:
                    continue
        return False  # if file only contains "A", "T", "C", "G", "N", it is unlikely to be a protein file


class BLASTAnalysis(metaclass=ABCMeta):
    def __init__(self):
        super().__init__()

    def init_tools(self):
        super().init_tools()
        self.mkblast_runner = MKBLASTRunner()
        self.tblastn_runner = TBLASTNRunner()

        if self.mkblast_runner.version != self.tblastn_runner.version:
            logger.warning(
                "You are using version {} of makeblastdb and version {} of tblastn.".format(
                    self.mkblast_runner.version, self.tblastn_runner.version
                )
            )

    def reset(self):
        super().reset()
        self.mkblast_runner.reset()
        self.tblastn_runner.reset()

    def _run_mkblast(self):
        self.mkblast_runner.configure_runner()
        if self.restart and self.mkblast_runner.check_previous_completed_run():
            logger.info(
                "Skipping makeblastdb as BLAST DB already exists at {}".format(
                    self.mkblast_runner.output_db
                )
            )
        else:
            self.restart = False  # Turn off restart mode if this is the entry point
            self.config.set("busco_run", "restart", str(self.restart))
            self.mkblast_runner.run()
        if len(os.listdir(os.path.split(self.mkblast_runner.output_db)[0])) == 0:
            raise BuscoError(
                "makeblastdb failed to create a BLAST DB at {}".format(
                    self.mkblast_runner.output_db
                )
            )

    def _run_tblastn(self, missing_and_frag_only=False, ancestral_variants=False):

        incomplete_buscos = (
            self.hmmer_runner.missing_buscos
            + list(self.hmmer_runner.fragmented_buscos.keys())
            if missing_and_frag_only
            else None
        )  # This parameter is only used on the re-run

        self.tblastn_runner.configure_runner(
            self.mkblast_runner.output_db,
            missing_and_frag_only,
            ancestral_variants,
            incomplete_buscos,
        )
        if self.restart and self.tblastn_runner.check_previous_completed_run():
            logger.info(
                "Skipping tblastn as results already exist at {}".format(
                    self.tblastn_runner.blast_filename
                )
            )
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.tblastn_runner.run()
        self.tblastn_runner.get_coordinates()
        self.tblastn_runner.filter_best_matches()
        self.tblastn_runner.write_coordinates_to_file()  # writes to "coordinates.tsv"
        self.tblastn_runner.write_contigs()
        return
