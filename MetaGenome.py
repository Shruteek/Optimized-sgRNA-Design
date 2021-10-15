from Bio import SeqIO
import SpacerSequence


class MetaGenome:
    __Name = ""
    __Sequences = []

    """Initialization function that """
    def __init__(self, metaGenome, name="MetaGenome"):
        self.__Name = name
        if not isinstance(metaGenome, str):
            print("The passed in metaGenome file name is not a string:" + metaGenome)
        elif not ".fasta" in metaGenome or not str.lower(metaGenome[-6:]) == ".fasta":
            print("The passed in metaGenome file name is not a .fasta file.")
        else:
            self.addSequences(metaGenome)

    """Getter method that returns the name of the MetaGenome."""
    def getName(self):
        return self.__Name

    """Method that uses the findOffTargets method of each Sequence in __Sequences, compiles them, and returns."""
    def findOffTargets(self, spacerSequence):
        offTargets = []
        for sequence in self.__Sequences:
            if not offTargets:
                offTargets = sequence.findOffTargets(spacerSequence)
            else:
                offTargets = [offTargets, sequence.findOffTargets(spacerSequence)]
        return offTargets

    """Adds sequences directly from a .FASTA file."""
    def addSequences(self, sequencesFile):
        for record in SeqIO.parse(sequencesFile, "fasta"):
            if isinstance(record, SeqIO.SeqRecord):
                print("Successful SeqRecord opened!")
            self.__Sequences.append(record)

    """Returns the sequence at the given index in __Sequences, ensuring the validity of the index."""
    def getSequenceAtIndex(self, index):
        if index < 0 or index >= len(self.__Sequences):
            print("Invalid index of " + index + "compared to __Sequences size of " + len(self.__Sequences))
            return None
        else:
            return self.__Sequences[index]

    """Returns the first sequence containing the given string in __Sequences, ensuring the validity of the sequence."""
    def getSequenceWithSubsequence(self, subsequence):
        if not SpacerSequence.isValidDNA(subsequence) and not SpacerSequence.isValidRNA(subsequence):
            print("Sequence entered is neither valid RNA nor valid DNA.")
            return None
        else:
            for sequence in self.__Sequences:
                if not sequence.seq.find(subsequence) == -1:
                    return sequence
            print("Sequence not find in __Sequences.")
            return None


