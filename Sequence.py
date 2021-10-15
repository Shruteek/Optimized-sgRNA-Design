from Bio import SeqIO
import SpacerSequence


class Sequence:
    __Name = ""
    __Record = None
    __MismatchStrictness = 3

    def __init__(self, sequenceRecord, name="Generic Sequence"):
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

    """Takes in a string name and sets the Sequence name accordingly."""
    def setName(self, name):
        if not isinstance(name, str) or len(name) <= 0:
            print("Given name is not a valid string name: " + name)
        else:
            self.__Name = name

    """Takes in a string keySequence and a string substrateSequence, cross-correlates keySequence along
    substrateSequence, and returns the numerical results as a list of integers."""
    def correlateSequences(self, keySequence, substrateSequence):
        if not SpacerSequence.isValidDNA(keySequence) or SpacerSequence.isValidDNA(substrateSequence):
            print("Given sequences are not valid DNA: " + keySequence + " and " + substrateSequence)
            return None
        else:
            crosscorrelation = []
            for shift in range(len(substrateSequence) - len(keySequence) + 1):
                correlation = 0
                for nucleotideIndex in keySequence:
                    if keySequence[nucleotideIndex] == substrateSequence[nucleotideIndex + shift]:
                        correlation = correlation + 1
                crosscorrelation.append(correlation)
            return crosscorrelation

    """"""
    def findOffTargets(self, spacerSequence):
        if not SpacerSequence.isValidRNA(spacerSequence):
            print("Given spacerSequence is not valid RNA: " + spacerSequence)
            return None
        elif not len(spacerSequence) == 20:
            print("Given spacerSequence is not 20 nucleotides long: " + spacerSequence)
            return None
        else:
            offTargets = []
            crosscorrelation = self.correlateSequences(spacerSequence, str(self.__Record.seq))
            for shift in range(len(crosscorrelation)):
                if crosscorrelation[shift] <= self.__MismatchStrictness:
                    offTargets.append(str(self.__Record.seq)[shift:(shift+len(spacerSequence))])
            return offTargets
