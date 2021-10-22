from Bio import SeqIO
from SpacerSequence import SpacerSequence
from Sequence import Sequence

class MetaGenome:

    """Initialization method that takes in the local .FASTA filename of a metaGenome (and, optionally, a name for that 
    metagenome) and instantiates."""
    def __init__(self, metaGenome, name="Generic MetaGenome"):
        self.__Name = name
        self.__Sequences = []
        if not self.isValidFasta(metaGenome):
            print("The given metaGenome file is not a valid .fasta string file name:" + metaGenome)
        else:
            self.addSequences(metaGenome)

    """Getter method that returns the name of the MetaGenome."""
    def getName(self):
        return self.__Name

    """Takes in a string name and sets the Sequence name accordingly."""

    def setName(self, name):
        if not isinstance(name, str) or len(name) <= 0:
            print("Given name is not a valid string name: " + name)
        else:
            self.__Name = name

    """Method that returns the size of the Sequence array."""
    def size(self):
        return len(self.__Sequences)

    """Returns a list containing either the string sequence at the given key index or the string sequences with the 
    given string key sequence in __Sequences, ensuring the validity of the key input as either an index or a string 
    sequence."""
    def getSequence(self, key):
        if isinstance(key, str):
            if not SpacerSequence.isValidDNA(None, key) and not SpacerSequence.isValidRNA(None, key):
                print("Subsequence entered is neither valid RNA nor valid DNA.")
                return [None]
            else:
                matchingSequences = []
                for sequence in self.__Sequences:
                    if not sequence.hasSubsequence(key):
                        matchingSequences.append(sequence.getSequence())
                    if len(matchingSequences) == 0:
                        return [None]
                    else:
                        return matchingSequences
        elif isinstance(key, int):
            if key < 0 or key >= len(self.__Sequences):
                print("Invalid given index of " + key + "for __Sequences size of " + len(self.__Sequences))
                return [None]
            else:
                return [self.__Sequences[key].getSequence()]
        else:
            print("Invalid sequence key: must be int (index) or str (subsequence).")

    """Method that takes an input String spacerSequence and uses the findOffTargets method of each Sequence in 
    __Sequences, compiles them, and returns."""
    def findOffTargets(self, spacerSequence):
        offTargets = []
        if not len(spacerSequence) == 20:
            print("Warning: given spacerSequence is not 20 nucleotides long: " + spacerSequence)
        for sequence in self.__Sequences:
            if not offTargets:
                offTargets = [sequence.findOffTargets(spacerSequence)]
            else:
                offTargets.append(sequence.findOffTargets(spacerSequence))
        return offTargets

    """Method that takes in the local filename of a sequence file and adds sequences directly from it to __Sequences."""
    def addSequences(self, sequencesFile):
        count = 1
        for record in SeqIO.parse(sequencesFile, "fasta"):
            if not isinstance(record, SeqIO.SeqRecord):
                print("Error: record parsed from .fasta is not a SeqRecord.")
            else:
                self.__Sequences.append(Sequence(record, "Sequence " + str(count)))
            count = count + 1

    """Takes in a filename and returns whether the given variable is a valid fasta filename"""
    def isValidFasta(self, sequencesFile):
        return isinstance(sequencesFile,str) and ".fasta" in str.lower(sequencesFile) and \
               str.lower(sequencesFile[-6:]) == ".fasta" and len(sequencesFile) > 6


