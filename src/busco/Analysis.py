from Bio import SeqIO
from busco.BuscoLogger import BuscoLogger
from abc import ABCMeta

logger = BuscoLogger.get_logger(__name__)


class NucleotideAnalysis(metaclass=ABCMeta):

    LETTERS = ["A", "C", "T", "G", "N"]

    # explanation of ambiguous codes found here: https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
    AMBIGUOUS_CODES = ["Y", "R", "W", "S", "K", "M", "D", "V", "H", "B"]

    def __init__(self):

        super().__init__()  # Initialize BuscoAnalysis
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

    def init_tools(self):
        super().init_tools()


class ProteinAnalysis:

    LETTERS = ["F", "L", "I", "M", "V", "S", "P", "T", "A", "Y", "X", "H", "Q", "N", "K", "D", "E", "C", "W", "R", "G"]
    NUCL_LETTERS = ["A", "C", "T", "G", "N"]

    def __init__(self):
        super().__init__()
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
