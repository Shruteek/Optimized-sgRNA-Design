from Bio import SeqIO
from GenomeTools import *
from Sequence import Sequence


class MetaGenome:

    def __init__(self, metaGenomePath, name="Generic MetaGenome"):
        """Initialization method that takes in the local .FASTA filename of a metaGenome (and, optionally, a name for that
            metagenome) and instantiates."""
        self.__Name = name
        if not isValidFasta(metaGenomePath):
            # print("The given metaGenome file is not a valid .fasta string file name:" + metaGenomePath)
            self.__Sequences = []
        else:
            self.__Sequences = []
            self.addSequences(metaGenomePath)

    def getName(self):
        """Getter method that returns the name of the MetaGenome."""
        return self.__Name

    def setName(self, name):
        """Takes in a string name and sets the Sequence name accordingly."""
        if isinstance(name, str) and len(name) > 0:
            self.__Name = name

    def size(self):
        """Method that returns the size of the Sequence array."""
        return len(self.__Sequences)

    def getSequence(self, key):
        """Returns a list containing either the string sequence at the given key index or the string sequences with the
            given string key sequence in __Sequences, ensuring the validity of the key input as either an index or a string
            sequence."""
        if isinstance(key, str):
            if not isValidDNA(key) and not isValidRNA(key):
                # print("Subsequence entered is neither valid RNA nor valid DNA.")
                return []
            else:
                matchingSequences = []
                for sequence in self.__Sequences:
                    if not sequence.hasSubsequence(key):
                        matchingSequences.append(sequence.getSequence())
                    if len(matchingSequences) == 0:
                        return matchingSequences
                    else:
                        return matchingSequences
        elif isinstance(key, int):
            if key < 0 or key >= len(self.__Sequences):
                # print("Invalid given index of " + key + "for __Sequences size of " + len(self.__Sequences))
                return []
            else:
                return [self.__Sequences[key].getSequence()]
        else:
            # print("Invalid sequence key: must be int (index) or str (subsequence).")
            return []

    def findOffTargets(self, spacerSequence):
        """Method that takes an input String RNA spacerSequence and uses the findOffTargets method of each Sequence in
            __Sequences, compiles them, and returns."""
        nestedOffTargets = []
        if not len(spacerSequence) == 20:
            # print("Warning: given spacerSequence is not 20 nucleotides long: " + spacerSequence)
            nestedOffTargets = []
        for sequence in self.__Sequences:
            if not nestedOffTargets:
                nestedOffTargets = [sequence.findOffTargets(spacerSequence)]
            else:
                nestedOffTargets.append(sequence.findOffTargets(spacerSequence))
        for offTargetsInSequence in nestedOffTargets:
            for offTargetIndex in len(offTargetsInSequence):
                if correlateSequences(complementaryDNA(spacerSequence), offTargetsInSequence[offTargetIndex][6:26]) == 20:
                    offTargetsInSequence.pop(offTargetIndex)
        return nestedOffTargets

    def addSequences(self, sequencesFile):
        """Method that takes in the filename of a sequence file and adds sequences directly from it to __Sequences."""
        count = 1
        for record in SeqIO.parse(sequencesFile, "fasta"):
            if isinstance(record, SeqIO.SeqRecord):
                self.__Sequences.append(Sequence(record, "Sequence " + str(count)))
            count = count + 1
