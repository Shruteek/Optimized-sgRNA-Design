from Bio import SeqIO
from GenomeTools import *


class Sequence:

    def __init__(self, sequenceRecord, name="Generic Sequence"):
        """Initialization method that takes in a SeqRecord (and, optionally, a name for the sequence) and
        instantiates."""
        self.__MismatchStrictness = 3
        self.__SeedMismatchStrictness = 1
        self.__Name = name
        if not isinstance(sequenceRecord, SeqIO.SeqRecord):
            # print("The given sequenceRecord is not a valid SeqRecord object.")
            self.__Record = None
        else:
            self.__Record = sequenceRecord

    def setMismatchStrictness(self, mismatches):
        """Takes an integer representing the maximum number of mismatches and sets the __MismatchStrictness
        accordingly."""
        if isinstance(mismatches, int) and mismatches > 0:
            self.__MismatchStrictness = mismatches

    def setSeedMismatchStrictness(self, mismatches):
        """Takes an integer representing the maximum number of seed sequence mismatches and sets the
            __SeedMismatchStrictness accordingly."""
        if isinstance(mismatches, int) and mismatches > 0:
            self.__SeedMismatchStrictness = mismatches

    def setName(self, name):
        """Takes in a string name and sets the Sequence name accordingly."""
        if isinstance(name, str) and len(name) > 0:
            self.__Name = name

    def getName(self):
        """Getter method that returns the name of the stored Sequence."""
        return self.__Name

    def getSequence(self):
        """Getter method that returns the string representation of the stored Sequence."""
        return str(self.__Record.seq)

    def getMismatchStrictness(self):
        """Getter method that returns the mismatch strictness of the Sequence."""
        return self.__MismatchStrictness

    def getSeedMismatchStrictness(self):
        """Getter method that returns the seed mismatch strictness of the Sequence."""
        return self.__SeedMismatchStrictness

    def hasSubsequence(self, subsequence):
        """Method that returns whether or not the given string subsequence occurs in the stored Sequence."""
        return not self.__Record.seq.find(subsequence) == -1

    def findTargetsFromSpacer(self, spacerSequence):
        """Takes in a 20 bp RNA string spacerSequence and returns all off-target sequences found using the GenomeTools
         class' findMatchingSequences method."""
        return findTargetsFromSpacer(spacerSequence, str(self.__Record.seq), self.__MismatchStrictness)

