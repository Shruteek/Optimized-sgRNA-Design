from Bio import SeqIO
from GenomeTools import *


class Sequence:

    def __init__(self, sequenceRecord, name="Generic Sequence"):
        """Initialization method that takes in a SeqRecord (and, optionally, a name for the sequence) and instantiates."""
        self.__Record = None
        self.__MismatchStrictness = 3
        self.__SeedMismatchStrictness = 1
        self.__Name = name
        if not isinstance(sequenceRecord, SeqIO.SeqRecord):
            print("The given sequenceRecord is not a valid SeqRecord object.")
        else:
            self.__Record = sequenceRecord

    def setMismatchStrictness(self, mismatches):
        """Takes an integer representing the maximum number of mismatches and sets the __MismatchStrictness accordingly."""
        if not isinstance(mismatches, int) or mismatches < 0:
            print("Given maximum mismatch strictness is not valid: " + mismatches)
        else:
            self.__MismatchStrictness = mismatches

    def setSeedMismatchStrictness(self, mismatches):
        """Takes an integer representing the maximum number of seed sequence mismatches and sets the
            __SeedMismatchStrictness accordingly."""
        if not isinstance(mismatches, int) or mismatches < 0:
            print("Given maximum mismatch strictness is not valid: " + mismatches)
        else:
            self.__SeedMismatchStrictness = mismatches

    def setName(self, name):
        """Takes in a string name and sets the Sequence name accordingly."""
        if not isinstance(name, str) or len(name) <= 0:
            print("Given name is not a valid string name: " + name)
        else:
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
        return self.__Record.seq.find(subsequence) == -1

    def findOffTargets(self, spacerSequence):
        """Takes in a DNA string spacerSequence and cross-correlates its seed sequence against the stored genome,
        identifying all genome sequences within the specified SeedMismatchStrictness, then correlates the whole
        spacerSequence against the initially matching genome sequences to identify and return all genome sequences
        within the specified MismatchStrictness."""
        if not isValidRNA(spacerSequence):
            print("Given spacerSequence is not valid RNA: " + spacerSequence)
            return None
        else:
            targetSequence = complementaryDNA(spacerSequence)
            offTargets = []
            substrateSequence = str(self.__Record.seq)
            crosscorrelation = crossCorrelateSequences(targetSequence, substrateSequence)
            for shift in range(len(crosscorrelation)):
                if crosscorrelation[shift] >= len(targetSequence) - self.__MismatchStrictness:
                    print("Off-target sequence: " + substrateSequence[shift + 6:(shift + 6 + 20 + 3)]
                          + " against " + targetSequence)
                    print("with a correlation of " + str(crosscorrelation[shift]))
                    offTargets.append(substrateSequence[shift:(shift + 6 + 6 + 20 + 3)])
            return offTargets
