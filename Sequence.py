from Bio import SeqIO
from SpacerSequence import SpacerSequence


class Sequence:

    """Initialization method that takes in a SeqRecord (and, optionally, a name for the sequence) and instantiates."""
    def __init__(self, sequenceRecord, name="Generic Sequence"):
        self.__Record = None
        self.__MismatchStrictness = 3
        self.__SeedMismatchStrictness = 1
        self.__Name = name
        if not isinstance(sequenceRecord, SeqIO.SeqRecord):
            print("The given sequenceRecord is not a valid SeqRecord object.")
        else:
            self.__Record = sequenceRecord

    """Takes an integer representing the maximum number of mismatches and sets the __MismatchStrictness accordingly."""
    def setMismatchStrictness(self, mismatches):
        if not isinstance(mismatches, int) or mismatches < 0:
            print("Given maximum mismatch strictness is not valid: " + mismatches)
        else:
            self.__MismatchStrictness = mismatches

    """Takes an integer representing the maximum number of seed sequence mismatches and sets the 
    __SeedMismatchStrictness accordingly."""
    def setSeedMismatchStrictness(self, mismatches):
        if not isinstance(mismatches, int) or mismatches < 0:
            print("Given maximum mismatch strictness is not valid: " + mismatches)
        else:
            self.__SeedMismatchStrictness = mismatches

    """Takes in a string name and sets the Sequence name accordingly."""
    def setName(self, name):
        if not isinstance(name, str) or len(name) <= 0:
            print("Given name is not a valid string name: " + name)
        else:
            self.__Name = name

    """Getter method that returns the name of the stored Sequence."""
    def getName(self):
        return self.__Name

    """Getter method that returns the string representation of the stored Sequence."""
    def getSequence(self):
        return str(self.__Record.seq)

    """Getter method that returns the mismatch strictness of the Sequence."""
    def getMismatchStrictness(self):
        return self.__MismatchStrictness

    """Getter method that returns the seed mismatch strictness of the Sequence."""

    def getSeedMismatchStrictness(self):
        return self.__SeedMismatchStrictness

    """Method that returns whether or not the given string subsequence occurs in the stored Sequence."""
    def hasSubsequence(self, subsequence):
        return self.__Record.seq.find(subsequence) == -1

    """Takes in a string keySequence and a string substrateSequence, cross-correlates keySequence along
    substrateSequence, and returns the numerical results as a list of integers."""
    def correlateSequences(self, keySequence, substrateSequence):
        if not SpacerSequence.isValidDNA(None, keySequence) or not SpacerSequence.isValidDNA(None, substrateSequence):
            print("Given sequences are not valid DNA: " + keySequence + " and " + substrateSequence)
            return []
        else:
            crosscorrelation = []
            for shift in range(len(substrateSequence) - len(keySequence) + 1):
                correlation = 0
                for nucleotideIndex in range(len(keySequence)):
                    if keySequence[nucleotideIndex] == substrateSequence[nucleotideIndex + shift]:
                        correlation = correlation + 1
                    elif substrateSequence[nucleotideIndex + shift] == "n":
                        correlation = correlation + 1
                crosscorrelation.append(correlation)
            return crosscorrelation

    """Tales in a string spacerSequence and cross-correlates it against the stored genome, returning all genome
    sequences within the specified MismatchStrictness."""
    def findOffTargets(self, spacerSequence):
        if not SpacerSequence.isValidRNA(None, spacerSequence):
            print("Given spacerSequence is not valid RNA: " + spacerSequence)
            return None
        else:
            sequenceToMatch = SpacerSequence.complementaryDNA(None, spacerSequence)
            offTargets = []
            crosscorrelation = self.correlateSequences(sequenceToMatch, str(self.__Record.seq))
            for shift in range(len(crosscorrelation)):
                if crosscorrelation[shift] >= len(sequenceToMatch) - self.__MismatchStrictness:
                    offTargets.append(str(self.__Record.seq)[shift:(shift+len(sequenceToMatch))])
            return offTargets

    """Tales in a string spacerSequence and cross-correlates its seed sequence against the stored genome, identifying 
    all genome sequences within the specified SeedMismatchStrictness, then correlates the whole spacerSequence against 
    the initially matching genome sequences to identify and return all genome sequences within the specified 
    MismatchStrictness."""

    def findOffTargetsEfficiently(self, spacerSequence):
        if not SpacerSequence.isValidRNA(None, spacerSequence):
            print("Given spacerSequence is not valid RNA: " + spacerSequence)
            return None
        else:
            sequenceToMatch = SpacerSequence.complementaryDNA(None, spacerSequence)
            offTargets = []
            crosscorrelation = self.correlateSequences(sequenceToMatch, str(self.__Record.seq))
            for shift in range(len(crosscorrelation)):
                if crosscorrelation[shift] >= len(sequenceToMatch) - self.__MismatchStrictness:
                    offTargets.append(str(self.__Record.seq)[shift:(shift + len(sequenceToMatch))])
            return offTargets
